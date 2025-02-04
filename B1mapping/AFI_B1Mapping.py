"""
Script :
    AFI_B1Mapping.py
Description :
    Generate B1 map using the AFI method, registered on the FLASH FA = 90°.
Needs :
    2022-12-20_gkg singularity container
    ANTS
    FSL
Usage :
    Previously initiate FSL
    python AFI_B1Mapping.py ./Working/Directory/ AFI_filename.nii.gz FLASH_FA90_filename.nii.gz TR1_in_ms TR2_in_ms
Authors :
    Lucas Arcamone
"""

import math
import optparse
import os
import shutil
import subprocess

import nibabel
import nibabel.processing
import numpy as np
from scipy.ndimage import binary_erosion


def safeAcos(arr) -> np.ndarray:
    mask_neg = arr < -1
    mask_pos = arr > 1
    ret = np.arccos(arr, where=~(mask_neg | mask_pos))
    ret[mask_neg] = math.pi
    ret[mask_pos] = 0
    return ret


def B1(WIP_directory, Afi_filename, TR1, TR2) -> None:
    output = os.path.join(WIP_directory, "b1.nii.gz")
    meta = nibabel.load(Afi_filename)
    arr = meta.get_fdata()
    with np.errstate(divide="ignore", invalid="ignore"):
        r = arr[:, :, :, 1] / arr[:, :, :, 0]
    n = TR2 / TR1
    res = safeAcos((r * n - 1) / (n - r)) * 180 / (60 * np.pi)  # Yarnykh, V. L. (2007).
    nibabel.save(nibabel.Nifti1Image(res, meta.affine), output)


def GetMask(Afi_filename, maskfilename, maskEroded) -> None:
    subprocess.run(
        [
            "singularity",
            "exec",
            "--bind",
            "/neurospin:/neurospin:rw",
            "/neurospin/phcp/code/gkg/2022-12-20_gkg/2022-12-20_gkg.sif",
            "GkgExecuteCommand",
            "GetMask",
            "-i",
            Afi_filename,
            "-o",
            maskfilename,
            "-a",
            "2",
            "-verbose",
        ]
    )
    meta_mask = nibabel.load(maskfilename)
    arr_mask = meta_mask.get_fdata()
    kernel = np.asarray(
        [
            [[0, 1, 0], [1, 1, 1], [0, 1, 0]],
            [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
            [[0, 1, 0], [1, 1, 1], [0, 1, 0]],
        ]
    )
    erosion = binary_erosion(arr_mask, kernel, 2)
    ni_im = nibabel.Nifti1Image(np.where(erosion, 1.0, 0.0), 0.2 * meta_mask.affine)
    nibabel.save(ni_im, maskEroded)


def cropB1(maskErodeFilename, b1Filename, output) -> None:
    meta_maskcorr = nibabel.load(maskErodeFilename)
    meta_B1 = nibabel.load(b1Filename)
    arr_maskcorr = meta_maskcorr.get_fdata()
    arr_maskcorr = np.flip(arr_maskcorr, 2)
    arr_res = meta_B1.get_fdata() * arr_maskcorr
    ni_im = nibabel.Nifti1Image(arr_res, meta_B1.affine)
    nibabel.save(ni_im, output)


def DilM(b1CropFilename, b1CropDilMFilename) -> None:
    subprocess.run(["fslmaths", b1CropFilename] + ["-dilM"] * 9 + [b1CropDilMFilename])


def NeutralValue(b1CropDilMFilename, b1CropDilMNeutralFilename) -> None:
    meta = nibabel.load(b1CropDilMFilename)
    arr = (
        meta.get_fdata() * 900
    )  # Neutral value is 900 in Gkg toolbox (Siemens standards)
    arr_filtered = np.where(arr < 5.0, 900, arr)
    ni_im = nibabel.Nifti1Image(
        np.where(np.isnan(arr_filtered), 900, arr_filtered), meta.affine
    )
    nibabel.save(ni_im, b1CropDilMNeutralFilename)


def Cropped_and_smooth(b1CropDilMNeutralFilename, b1Smoothed) -> None:
    meta = nibabel.load(b1CropDilMNeutralFilename)
    fwhm = nibabel.processing.sigma2fwhm(2)
    B1AfterRegistration_smoothed = nibabel.processing.smooth_image(meta, fwhm)
    nibabel.save(B1AfterRegistration_smoothed, b1Smoothed)


def Rigid_registration(
    Afi_filename,
    Afi_extracted,
    FA90_filename,
    transformation,
    afi_registred,
    afi_registred_inv,
) -> None:
    meta = nibabel.load(Afi_filename)
    ni_im = nibabel.Nifti1Image(meta.get_fdata()[:, :, :, 0], meta.affine)
    nibabel.save(ni_im, Afi_extracted)
    subprocess.run(
        [
            "antsRegistration",
            "-d",
            "3",
            "-r",
            "[" + FA90_filename + "," + Afi_extracted + ",1]",
            "-t",
            "Rigid[0.1]",
            "-m",
            "MI[" + FA90_filename + "," + Afi_extracted + ",1,100,Regular,0.5]",
            "-o",
            "[" + transformation + "," + afi_registred + "," + afi_registred_inv + "]",
            "-c",
            "[100x100x100, 1e-6, 5]",
            "-s",
            "4x2x1",
            "-f",
            "4x2x1",
            "-v",
            "1",
        ]
    )


def Apply_registration(
    FA90_filename, transformation, b1Smoothed, b1SmoothedAfterRegistration
) -> None:
    subprocess.run(
        [
            "antsApplyTransforms",
            "-d",
            "3",
            "-i",
            b1Smoothed,
            "-r",
            FA90_filename,
            "-o",
            b1SmoothedAfterRegistration,
            "-n",
            "BSpline[3]",
            "-t",
            transformation,
            "-v",
        ]
    )


##  MAIN part


def main(MaterialDirectory, Afi_filename, FA90_filename, TR1, TR2):
    # creation of the material directory
    WIP_directory = os.path.join(MaterialDirectory, "WIP")
    os.mkdir(WIP_directory)

    # Creation of b1
    B1(WIP_directory, Afi_filename, int(TR1), int(TR2))

    # Creation of the eroded mask
    b1filename = os.path.join(WIP_directory, "b1.nii.gz")
    maskb1filename = os.path.join(WIP_directory, "mask_b1.nii.gz")
    maskErodedfilename = os.path.join(WIP_directory, "mask_b1_eroded.nii.gz")
    GetMask(Afi_filename, maskb1filename, maskErodedfilename)

    # b1 masking
    b1cropfilename = os.path.join(WIP_directory, "b1_crop.nii.gz")
    cropB1(maskErodedfilename, b1filename, b1cropfilename)

    # Median dilation
    b1cropdilmfilename = os.path.join(WIP_directory, "b1_crop_dilm.nii.gz")
    DilM(b1cropfilename, b1cropdilmfilename)

    # Neutrality
    b1cropdilmneutralfilename = os.path.join(
        WIP_directory, "b1_crop_dilm_neutral.nii.gz"
    )
    NeutralValue(b1cropdilmfilename, b1cropdilmneutralfilename)

    # Smoothing
    b1smoothedfilename = os.path.join(WIP_directory, "b1_smoothed.nii.gz")
    Cropped_and_smooth(b1cropdilmneutralfilename, b1smoothedfilename)

    # Registration
    Afi_extracted = os.path.join(WIP_directory, "afi_extracted.nii.gz")
    transformation = os.path.join(WIP_directory, "AfiExtractedToFA90")
    afi_registred = os.path.join(WIP_directory, "afi_registred.nii.gz")
    afi_registred_inv = os.path.join(WIP_directory, "afi_registred_inv.nii.gz")
    Rigid_registration(
        Afi_filename,
        Afi_extracted,
        FA90_filename,
        transformation,
        afi_registred,
        afi_registred_inv,
    )

    transformationfinal = os.path.join(
        WIP_directory, "AfiExtractedToFA900GenericAffine.mat"
    )
    b1registredfilename = os.path.join(MaterialDirectory, "b1_registred.nii.gz")
    Apply_registration(
        FA90_filename, transformationfinal, b1smoothedfilename, b1registredfilename
    )

    # Remove WIP directory
    shutil.rmtree(WIP_directory)


parser = optparse.OptionParser()
parser.add_option(
    "-m", "--material", dest="materialDirectory", help="Material directory"
)
parser.add_option("-a", "--afi", dest="AFI", help="AFIfilename")
parser.add_option("-f", "--FA90", dest="FA90", help="Fa90 filename")
parser.add_option("-t", "--TR1", dest="TR1", help="TR1")
parser.add_option("-u", "--TR2", dest="TR2", help="TR2")

(options, args) = parser.parse_args()


main(options.materialDirectory, options.AFI, options.FA90, options.TR1, options.TR2)
