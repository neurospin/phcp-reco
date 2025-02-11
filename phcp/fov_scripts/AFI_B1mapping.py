"""
Script :
    AFI_B1Mapping.py
Description :
    Generate B1 map using the AFI method
Needs :
    2022-12-20_gkg singularity container
    ANTS
    FSL
Usage :
    Previously initiate FSL
    python AFI_B1Mapping.py img_TB1AFI.nii.gz img_TB1map.nii.gz
Authors :
    Lucas Arcamone
"""

import logging
import math
import os
import shutil
import sys
import tempfile

import nibabel
import nibabel.processing
import numpy as np
from scipy.ndimage import binary_erosion

from phcp.fsl import run_fsl_command
from phcp.gkg import run_gkg_GetMask


logger = logging.getLogger(__name__)

DEFAULT_TR1 = 15.0
DEFAULT_TR2 = 45.0


def safeAcos(arr) -> np.ndarray:
    mask_neg = arr < -1
    mask_pos = arr > 1
    ret = np.arccos(arr, where=~(mask_neg | mask_pos))
    ret[mask_neg] = math.pi
    ret[mask_pos] = 0
    return ret


def compute_raw_B1map(
    Afi_filename: str, B1map_filename: str, TR1: float, TR2: float
) -> None:
    """Compute a B1 map from the AFI using the formula of Yarnykh et al. 2007."""
    meta = nibabel.load(Afi_filename)
    arr = meta.get_fdata()
    n = TR2 / TR1
    with np.errstate(divide="ignore", invalid="ignore"):
        r = arr[:, :, :, 1] / arr[:, :, :, 0]
        res = safeAcos((r * n - 1) / (n - r)) * (
            180 / (60 * np.pi)
        )  # Yarnykh, V. L. (2007).
    nibabel.save(nibabel.Nifti1Image(res, meta.affine), B1map_filename)


def compute_mask_from_AFI(
    Afi_filename: str, maskfilename: str, maskEroded: str
) -> None:
    """Compute a mask and eroded mask from an AFI image."""
    run_gkg_GetMask(
        [
            "-i",
            Afi_filename,
            "-o",
            maskfilename,
            "-a",
            "2",
            "-verbose",
        ],
        input_dirs=[os.path.dirname(Afi_filename)],
        output_dirs=[os.path.dirname(maskfilename)],
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
    # Gkg GetMask sets the sform incorrectly to the identity matrix, so we have
    # to use the qform explicitly instead of the "best affine".
    affine, qform_code = meta_mask.header.get_qform(coded=True)
    assert qform_code == 1  # ensure qform points to scanner-based coordinates
    ni_im = nibabel.Nifti1Image(np.where(erosion, 1.0, 0.0), affine)
    nibabel.save(ni_im, maskEroded)


def nibabel_orient_like(
    img: nibabel.spatialimages.SpatialImage, ref: nibabel.spatialimages.SpatialImage
) -> nibabel.spatialimages.SpatialImage:
    img_ornt = nibabel.orientations.io_orientation(img.affine)
    ref_ornt = nibabel.orientations.io_orientation(ref.affine)
    return img.as_reoriented(nibabel.orientations.ornt_transform(img_ornt, ref_ornt))


def apply_mask_to_B1(maskErodeFilename: str, b1Filename: str, output: str) -> None:
    meta_B1 = nibabel.load(b1Filename)
    meta_maskcorr = nibabel_orient_like(nibabel.load(maskErodeFilename), meta_B1)
    arr_maskcorr = meta_maskcorr.get_fdata()
    arr_res = meta_B1.get_fdata() * arr_maskcorr
    ni_im = nibabel.Nifti1Image(arr_res, meta_B1.affine)
    nibabel.save(ni_im, output)


def dilate_in_background(b1CropFilename: str, b1CropDilMFilename: str) -> None:
    run_fsl_command(
        ["fslmaths", b1CropFilename] + ["-dilM"] * 9 + [b1CropDilMFilename], check=True
    )


def scale_to_900_and_filter(
    b1CropDilMFilename: str, b1CropDilMNeutralFilename: str
) -> None:
    """Scale the B1map from relative to decidegrees.

    Values very close to zero and NaNs are replaced by the neutral value, 900.
    """
    meta = nibabel.load(b1CropDilMFilename)
    arr = (
        meta.get_fdata() * 900
    )  # Neutral value is 900 in Gkg toolbox (Siemens standards)
    arr_filtered = np.where(arr < 5.0, 900, arr)
    ni_im = nibabel.Nifti1Image(
        np.where(np.isnan(arr_filtered), 900, arr_filtered), meta.affine
    )
    nibabel.save(ni_im, b1CropDilMNeutralFilename)


def smooth_B1map(b1CropDilMNeutralFilename: str, b1Smoothed: str) -> None:
    """Apply a smoothing (2 millimetres Gaussian) to the B1 map."""
    meta = nibabel.load(b1CropDilMNeutralFilename)
    fwhm = nibabel.processing.sigma2fwhm(2)
    B1_smoothed = nibabel.processing.smooth_image(meta, fwhm)
    nibabel.save(B1_smoothed, b1Smoothed)


def AFI_B1Mapping(
    Afi_filename: str,
    B1map_filename: str,
    TR1: float = DEFAULT_TR1,
    TR2: float = DEFAULT_TR2,
    work_dir: str = None,
) -> None:
    # creation of the working directory
    try:
        if work_dir is None:
            delete_work_dir = True
            work_dir = tempfile.mkdtemp(prefix="AFI_B1Mapping_")
            logger.info("Creating temporary working directory %s", work_dir)
        else:
            delete_work_dir = False

        # Creation of b1
        b1filename = os.path.join(work_dir, "b1_raw.nii.gz")
        logger.info("Computing raw B1 map: %s", b1filename)
        compute_raw_B1map(Afi_filename, b1filename, TR1, TR2)

        # Creation of the eroded mask
        maskb1filename = os.path.join(work_dir, "mask_b1.nii.gz")
        maskErodedfilename = os.path.join(work_dir, "mask_b1_eroded.nii.gz")
        logger.info("Computing mask: %s", maskErodedfilename)
        compute_mask_from_AFI(Afi_filename, maskb1filename, maskErodedfilename)

        # b1 masking
        b1cropfilename = os.path.join(work_dir, "b1_crop.nii.gz")
        logger.info("Applying mask: %s", b1cropfilename)
        apply_mask_to_B1(maskErodedfilename, b1filename, b1cropfilename)

        # Median dilation
        b1cropdilmfilename = os.path.join(work_dir, "b1_crop_dilm.nii.gz")
        logger.info("Dilating B1 map: %s", b1cropdilmfilename)
        dilate_in_background(b1cropfilename, b1cropdilmfilename)

        # Neutrality
        b1cropdilmneutralfilename = os.path.join(
            work_dir, "b1_crop_dilm_neutral.nii.gz"
        )
        logger.info("Scaling and filtering B1 map: %s", b1cropdilmneutralfilename)
        scale_to_900_and_filter(b1cropdilmfilename, b1cropdilmneutralfilename)

        # Smoothing
        logger.info("Smoothing B1 map: %s", B1map_filename)
        smooth_B1map(b1cropdilmneutralfilename, B1map_filename)
    finally:
        if delete_work_dir and work_dir is not None:
            logger.info("Deleting temporary working directory %s", work_dir)
            shutil.rmtree(work_dir)


def parse_command_line(argv):
    """Parse the script's command line."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Calculate a map of transmit field inhomogeneity (B1+) "
        "based on an Actual Flip angle Imaging (AFI) acquisition",
    )
    parser.add_argument(
        "AFI", help="Input raw AFI scan (normally a BIDS file named *_TB1AFI.nii.gz"
    )
    parser.add_argument("TB1map", help="Output map of the Transmit B1 field")
    parser.add_argument(
        "-t", "--TR1", type=float, default=DEFAULT_TR1, help="TR1 in milliseconds"
    )
    parser.add_argument(
        "-u", "--TR2", type=float, default=DEFAULT_TR2, help="TR2 in milliseconds"
    )
    parser.add_argument("--work-dir", type=str, default=None)

    args = parser.parse_args()
    return args


def main(argv=sys.argv):
    """The script's entry point."""
    logging.basicConfig(level=logging.INFO)
    args = parse_command_line(argv)
    return (
        AFI_B1Mapping(
            args.AFI,
            args.TB1map,
            TR1=args.TR1,
            TR2=args.TR2,
            work_dir=args.work_dir,
        )
        or 0
    )


if __name__ == "__main__":
    sys.exit(main())
