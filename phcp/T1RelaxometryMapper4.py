#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 09:44:07 2023

@author: la272118
"""

import os
import json
import optparse
import math
import numpy as np
import nibabel as ni
import glob
import nibabel.processing as proc
import ants
from sklearn.mixture import GaussianMixture


""" Miscellaneous algorithms """


def print_message(message):
    print("==================================")
    print(message)
    print("==================================")


def ConversionGisToNifti(InputGisFilename, OutputNiftiFilename):
    os.system(
        "singularity exec --bind /neurospin:/neurospin:rw "
        + " /neurospin/phcp/code/gkg/2022-12-20_gkg/2022-12-20_gkg.sif "
        + "GkgExecuteCommand Gis2NiftiConverter -i "
        + InputGisFilename
        + " -o "
        + OutputNiftiFilename
        + " -verbose"
    )


def create_ArrayFromGisFilename(GisFilename):
    DimFilename = GisFilename[:-3] + "dim"
    TxtInDimFile = open(DimFilename, "r").read()
    TxtFirstLine = TxtInDimFile.split("\n")[0]
    shape = tuple(np.int16(TxtFirstLine.split(" ")))  # [:-1]
    arr = np.fromfile(GisFilename, np.float32())
    arr = np.reshape(arr, shape, order="F")
    return arr


def arrange_ArrayToFlipHeaderFormat(arr):
    newarr = np.swapaxes(arr, 1, 2)
    newarr = np.flip(newarr, 1)
    newarr = np.flip(newarr, 0)
    return newarr


def dearrange_ArrayToFlipHeaderFormat(arr):
    newarr = np.flip(arr, 0)
    newarr = np.flip(newarr, 1)
    newarr = np.swapaxes(newarr, 1, 2)
    return newarr


def VFAregister(img_filename, output_filename, verbose):
    if verbose:
        print("========= Loading data =========")

    img_4D = ni.load(img_filename)
    arr_4D = img_4D.get_fdata()
    volume_number = arr_4D.shape[-1]
    ANTS_fixed_image = ants.from_numpy(arr_4D[:, :, :, 0])
    registered_volumes = []
    registered_volumes.append(ANTS_fixed_image)

    for i in range(1, volume_number):
        print("========= Registration volume n°" + str(i) + " =========")
        ANTS_moving_image = ants.from_numpy(arr_4D[:, :, :, i])
        registration = ants.registration(
            ANTS_fixed_image,
            ANTS_moving_image,
            type_of_transform="TRSAA",
            aff_iterations=(200, 200, 200),
            aff_shrink_factors=(8, 4, 2),
            aff_smoothing_sigmas=(3, 2, 1),
            verbose=False,
        )
        warped_image = registration["warpedmovout"]
        registered_volumes.append(warped_image)

    registered_data = np.stack([img.numpy() for img in registered_volumes], axis=-1)
    ni_im = ni.Nifti1Image(registered_data, img_4D.affine)
    ni.save(ni_im, output_filename)


""" Single Compartment T1 Relaxometry Mapper """


def createCommand_SubVolume_gkgMethod(fileNameVFACat, t1FileName):
    command = (
        "GkgExecuteCommand SubVolume"
        + " -i "
        + fileNameVFACat
        + " -o "
        + t1FileName
        + " -tIndices "
        + str(0)
        + " -verbose true"
    )
    return command


def createCommand_Mask_gkgMethod(t1FileName, fileNameMask):
    command = (
        "GkgExecuteCommand GetMask"
        + " -i "
        + t1FileName
        + " -o "
        + fileNameMask
        + " -a 2 "
        + " -verbose "
    )
    return command


def createCommand_T1SingleCompartmentRelaxometryMapper(
    fileNameVFACat,
    fileNameMask,
    fileNameTRValues,
    fileNameFAValues,
    fileNameB1,
    fileNameProtonDensity,
    fileNameOutputT1,
    fileNameFittedVFA,
    verbose,
):
    # StringParameter =  ( str( fileNameTRValues ), str( fileNameFAValues ), str( fileNameB1 ) )
    command = (
        "GkgExecuteCommand SingleCompartmentRelaxometryMapper"
        + " -i %s" % fileNameVFACat
        + " -m %s" % fileNameMask
        + " -t t1-mapping-vfa-spgr"
        + " -optimizerParameters %s " % 50000
        + " -optimizerParameters %s " % 0.001
        + " -optimizerParameters %s" % 1
        + " -optimizerParameters %s" % 200
        + " -optimizerParameters %s" % 100
        + " -optimizerParameters %s" % 10
        + " -optimizerParameters %s" % 1000
        + " -scalarParameters %s" % 5
        + " -scalarParameters %s" % 50
        + " -scalarParameters %s" % 0
        + " -scalarParameters %s" % 5
        + " -scalarParameters %s" % 50
        + " -scalarParameters %s" % 500
        + " -scalarParameters %s" % 1
        + " -scalarParameters %s" % 20
        + " -scalarParameters %s" % 0.002
        + " -stringParameters "
        + str(fileNameTRValues)
        + " -stringParameters "
        + str(fileNameFAValues)
        + " -stringParameters "
        + str(fileNameB1)
        + " -op %s" % fileNameProtonDensity
        + " -ot %s" % fileNameOutputT1
        + " -f %s" % fileNameFittedVFA
        + " -ascii False"
        + " -format gis"
        + " -verbose %s" % verbose
    )
    return command


def run_command(command):
    os.system(
        "singularity exec --bind /neurospin:/neurospin:rw "
        + " /neurospin/phcp/code/gkg/2022-12-20_gkg/2022-12-20_gkg.sif "
        + command
    )
    return None


def write_TR_FA(fileNameFAValues, fileNameTRValues, FAFilenames, fileNameB1):
    with open(fileNameFAValues, "w") as file:
        with open(fileNameTRValues, "w") as file2:
            for metadata_filename in sorted(glob.iglob(FAFilenames)):
                with open(metadata_filename) as metadata_file:
                    metadata = json.load(metadata_file)
                if os.path.exists(fileNameB1):
                    file.write(str(metadata["FlipAngle"]) + "\n")
                else:
                    file.write(
                        str(math.radians(metadata["FlipAngle"])) + "\n"
                    )  # Si pas de carte  de B1 --> convertir les angles en radians
                file2.write(str(metadata["RepetitionTime"] * 1000) + "\n")

    return None


def create_ConfidenceMap(
    fileNameVFACat, fileNameFittedVFA, fileNameMask, fileNameT1nifti
):
    arr_fittedvfa = create_ArrayFromGisFilename(fileNameFittedVFA)
    arr_fittedvfa = arrange_ArrayToFlipHeaderFormat(arr_fittedvfa)

    meta_vfa = ni.load(fileNameVFACat)
    arr_vfa = meta_vfa.get_fdata()

    meta_mask = ni.load(fileNameMask)
    arr_mask = meta_mask.get_fdata()
    arr_mask = arrange_ArrayToFlipHeaderFormat(arr_mask)

    diff = arr_fittedvfa - arr_vfa
    std_diff = np.std(diff, axis=-1) * arr_mask

    metaT1nifti = ni.load(fileNameT1nifti)

    ni_im = ni.Nifti1Image(
        dearrange_ArrayToFlipHeaderFormat(std_diff),
        metaT1nifti.affine,
        metaT1nifti.header,
    )
    return ni_im


""" Correction Part """


def extract_BiasField(ants_T1img, ants_T1img_mask, verbose):
    return ants.n4_bias_field_correction(
        ants_T1img,
        ants_T1img_mask,
        spline_param=20,
        rescale_intensities=True,
        return_bias_field=True,
        verbose=verbose,
    )


def create_NewBiasField(biasField):
    flat_arr = biasField[biasField > 0].reshape(-1, 1)
    gmm = GaussianMixture(n_components=4)
    gmm.fit(flat_arr)

    means = gmm.means_.flatten().tolist()
    print(means)
    moy_inconsistency, index_inconsistency = min(means), means.index(min(means))
    means.remove(moy_inconsistency)
    moy_normal = min(means)  # index_normal = means.index(min(means))

    stds = np.sqrt(gmm.covariances_.flatten()).tolist()
    std_inconsistency = stds[index_inconsistency]

    biasField_numpy = biasField.numpy()
    return np.where(
        biasField_numpy < (moy_inconsistency + 2 * std_inconsistency),
        moy_normal / biasField_numpy,
        1,
    )


def correct_qt1(T1GisFilename, verbose=True):
    if verbose:
        print_message("Loading data")

    T1NiftiFilename = T1GisFilename.split(".")[0] + ".nii.gz"

    if not (os.path.exists(T1NiftiFilename)):
        ConversionGisToNifti(T1GisFilename, T1NiftiFilename)

    ants_T1img = ants.image_read(T1NiftiFilename)
    ants_T1img_mask = ants_T1img.get_mask()

    if verbose:
        print_message("Bias Field Processing")

    biasField = extract_BiasField(ants_T1img, ants_T1img_mask, verbose)

    ants.image_write(biasField, T1GisFilename.split(".")[0] + "BiasField1.nii.gz")

    if verbose:
        print_message("New Bias Field Generation")

    newbiasField = create_NewBiasField(biasField)

    meta = ni.load(T1NiftiFilename)
    fwhm = proc.sigma2fwhm(6)
    ni_im = ni.Nifti1Image(np.float32(newbiasField), meta.affine, meta.header)

    if verbose:
        print_message("Smoothing part")

    newbiasField_smooth = proc.smooth_image(ni_im, fwhm)
    newbiasField_smooth_Filename = T1GisFilename.split(".")[0] + "_BiasField.nii.gz"
    ni.save(newbiasField_smooth, newbiasField_smooth_Filename)

    newbiasField_smooth = ants.from_nibabel(newbiasField_smooth)

    unbiased_T1 = ants_T1img * newbiasField_smooth

    if verbose:
        print_message("Saving Part")

    T1Nifti_unbiased_Filename = T1GisFilename.split(".")[0] + "_rec-unbiased.nii.gz"

    ants.image_write(unbiased_T1, T1Nifti_unbiased_Filename)
    return None


def runT1RelaxometryMapper(
    subjectDirectoryGisConversion, FAFilenames, outputDirectory, verbose
):
    fileNameVFACat0 = os.path.join(subjectDirectoryGisConversion, "t1map.nii.gz")
    fileNameVFACat = os.path.join(subjectDirectoryGisConversion, "t1map-TRSAA.nii.gz")

    if not (os.path.exists(fileNameVFACat)):
        if verbose:
            print_message("T1 RELAXOMETRY : VFA REGISTRATION")
        VFAregister(fileNameVFACat0, fileNameVFACat, verbose)

    if verbose:
        print_message("SINGLE COMPARTMENT T1 RELAXOMETRY MAPPER")

    fileNameMask = os.path.join(subjectDirectoryGisConversion, "mask.nii.gz")
    t1FileName = os.path.join(subjectDirectoryGisConversion, "t1_extracted.nii.gz")
    fileNameB1 = os.path.join(subjectDirectoryGisConversion, "b1_registred.nii.gz")
    fileNameTRValues = os.path.join(outputDirectory, "tr.txt")
    fileNameFAValues = os.path.join(outputDirectory, "fa.txt")

    if not (os.path.exists(t1FileName)):
        command = createCommand_SubVolume_gkgMethod(fileNameVFACat, t1FileName)
        run_command(command)

    if not (os.path.exists(fileNameMask)):
        command = createCommand_Mask_gkgMethod(t1FileName, fileNameMask)
        run_command(command)

    fileNameT1nifti = os.path.join(outputDirectory, "T1.nii.gz")
    fileNameProtonDensity = os.path.join(outputDirectory, "proton-density.ima")
    fileNameFittedVFA = os.path.join(outputDirectory, "fitted-vfa.ima")
    fileNameOutputT1 = os.path.join(outputDirectory, "T1.ima")

    if not (os.path.exists(fileNameT1nifti)):
        write_TR_FA(fileNameFAValues, fileNameTRValues, FAFilenames, fileNameB1)

        command = createCommand_T1SingleCompartmentRelaxometryMapper(
            fileNameVFACat,
            fileNameMask,
            fileNameTRValues,
            fileNameFAValues,
            fileNameB1,
            fileNameProtonDensity,
            fileNameOutputT1,
            fileNameFittedVFA,
            verbose,
        )
        run_command(command)

        ConversionGisToNifti(fileNameOutputT1, fileNameT1nifti)

    fileNameConfidenceMap = os.path.join(outputDirectory, "T1ConfidenceMap.nii.gz")

    if not (os.path.exists(fileNameConfidenceMap)):
        ni_im = create_ConfidenceMap(
            fileNameVFACat, fileNameFittedVFA, fileNameMask, fileNameT1nifti
        )
        ni.save(ni_im, fileNameConfidenceMap)

    if verbose:
        print_message("Correction Part")

    correct_qt1(fileNameOutputT1, verbose)
    return None


################################################################################
# parser to get option(s)
################################################################################

parser = optparse.OptionParser()
parser.add_option("-i", "--input", dest="vfadirectory", help="VFA Directory")
parser.add_option(
    "-m",
    "--vfa",
    dest="VFAFilenames",
    help="String for glob research of MGE volumes. Ex : /phcp/rawdata/sub/ses/anat/sub_ses_flip*_VFA.json",
)
parser.add_option(
    "-o", "--outputDirectory", dest="outputDirectory", help="Output directory"
)
parser.add_option("-v", "--verbose", dest="verbose", default=True, help="verbose")

(options, args) = parser.parse_args()


################################################################################
# 1) Value Extractor
################################################################################

runT1RelaxometryMapper(
    options.vfadirectory, options.VFAFilenames, options.outputDirectory, options.verbose
)
