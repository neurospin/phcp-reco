#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 10:03:08 2023

@author: la272118
"""

import os
import json
import optparse
import numpy as np
import nibabel as ni
import glob


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


""" Single Compartment T2Star Relaxometry Mapper """


def createCommand_SubVolume_gkgMethod(fileNameMGE, t2starFileName):
    command = (
        "GkgExecuteCommand SubVolume"
        + " -i "
        + fileNameMGE
        + " -o "
        + t2starFileName
        + " -tIndices "
        + str(0)
        + " -verbose true"
    )
    return command


def createCommand_Mask_gkgMethod(fileNameMGE, fileNameMask):
    command = (
        "GkgExecuteCommand GetMask"
        + " -i "
        + fileNameMGE
        + " -o "
        + fileNameMask
        + " -a 2 "
        + " -verbose "
    )
    return command


def createCommand_T2starSingleCompartmentRelaxometryMapper(
    fileNameMGE,
    fileNameMask,
    fileNameTEValues,
    fileNameProtonDensity,
    fileNameT2star,
    fileNameFittedMGE,
    verbose,
):
    # optimizerParameters :: <P1> : NLP maximum iteration count  <P2> : NLP test size <P3> : 0->do not apply MCMC 1->apply MCMC <P4> : MCMC burnin count <P5> : MCMC sample count <P6> : MCMC interval count <P7> : MCMC maximum iteration coun
    # ' -optimizerParameters 50000 0.001 0 2000 300 50 10000' +\  ' -scalarParameters 0.75 50 0 0 8 180 1 20 5'  +\
    command = (
        "GkgExecuteCommand SingleCompartmentRelaxometryMapper"
        + " -i %s" % fileNameMGE
        + " -m %s" % fileNameMask
        + " -t t2-star-mapping-mgre"
        + " -optimizerParameters %s " % 50000
        + " -optimizerParameters %s " % 0.001
        + " -optimizerParameters %s" % 0
        + " -optimizerParameters %s" % 2000
        + " -optimizerParameters %s" % 300
        + " -optimizerParameters %s" % 50
        + " -optimizerParameters %s" % 10000
        + " -scalarParameters %s" % 0.3
        + " -scalarParameters %s" % 30
        + " -scalarParameters %s" % 0
        + " -scalarParameters %s" % 0
        + " -scalarParameters %s" % 8
        + " -scalarParameters %s" % 100
        + " -scalarParameters %s" % 1
        + " -scalarParameters %s" % 10
        + " -scalarParameters %s" % 2
        + " -stringParameters %s" % fileNameTEValues
        + " -op %s" % fileNameProtonDensity
        + " -ot %s" % fileNameT2star
        + " -f %s" % fileNameFittedMGE
        + " -ascii False"
        + " -format nifti"
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


def write_EchotTime(fileNameTEValues, MSMEFilenames):
    with open(fileNameTEValues, "w") as file:
        for metadata_filename in sorted(glob.iglob(MSMEFilenames)):
            with open(metadata_filename) as metadata_file:
                metadata = json.load(metadata_file)
            file.write(str(metadata["EchoTime"] * 1000) + "\n")
    return None


def create_ConfidenceMap(
    fileNameMGE, fileNameFittedMGE, fileNameMask, fileNameT2starnifti
):
    meta_fittedmge = ni.load(fileNameFittedMGE)
    arr_fittedmge = meta_fittedmge.get_fdata()
    arr_fittedmge = arrange_ArrayToFlipHeaderFormat(arr_fittedmge)

    meta_mge = ni.load(fileNameMGE)
    arr_mge = meta_mge.get_fdata()

    meta_mask = ni.load(fileNameMask)
    arr_mask = meta_mask.get_fdata()
    arr_mask = arrange_ArrayToFlipHeaderFormat(arr_mask)

    diff = arr_fittedmge - arr_mge
    std_diff = np.std(diff[:, :, :, :6], axis=-1) * arr_mask

    metaT2starnifti = ni.load(fileNameT2starnifti)

    ni_im = ni.Nifti1Image(
        dearrange_ArrayToFlipHeaderFormat(std_diff),
        metaT2starnifti.affine,
        metaT2starnifti.header,
    )
    return ni_im


def runT2StarRelaxometryMapper(
    subjectDirectoryGisConversion, MGEFilenames, outputDirectory, verbose
):
    if verbose:
        print_message("SINGLE COMPARTMENT T2Star RELAXOMETRY MAPPER")

    fileNameMGE = os.path.join(subjectDirectoryGisConversion, "t2star-mge.nii.gz")
    fileNameMask = os.path.join(subjectDirectoryGisConversion, "mask.nii.gz")
    fileNameTEValues = os.path.join(outputDirectory, "te.txt")
    t2starFileName = os.path.join(
        subjectDirectoryGisConversion, "t2star_extracted.nii.gz"
    )

    if not (os.path.exists(t2starFileName)):
        command = createCommand_SubVolume_gkgMethod(fileNameMGE, t2starFileName)
        run_command(command)

    if not (os.path.exists(fileNameMask)):
        command = createCommand_Mask_gkgMethod(fileNameMGE, fileNameMask)
        run_command(command)

    write_EchotTime(fileNameTEValues, MGEFilenames)

    fileNameProtonDensity = os.path.join(outputDirectory, "proton-density.nii.gz")
    fileNameFittedMGE = os.path.join(outputDirectory, "fitted-mge.nii.gz")
    fileNameT2star = os.path.join(outputDirectory, "T2star.nii.gz")

    command = createCommand_T2starSingleCompartmentRelaxometryMapper(
        fileNameMGE,
        fileNameMask,
        fileNameTEValues,
        fileNameProtonDensity,
        fileNameT2star,
        fileNameFittedMGE,
        verbose,
    )
    run_command(command)

    ni_im = create_ConfidenceMap(
        fileNameMGE, fileNameFittedMGE, fileNameMask, fileNameT2star
    )
    ni.save(ni_im, os.path.join(outputDirectory, "T2starConfidenceMap.nii.gz"))
    return None


################################################################################
# parser to get option(s)
################################################################################

parser = optparse.OptionParser()
parser.add_option("-i", "--input", dest="mgedirectory", help="MEGRE Directory")
parser.add_option(
    "-m",
    "--mge",
    dest="MGEFilenames",
    help="String for glob research of MGE volumes. Ex : /phcp/rawdata/sub/ses/anat/sub_ses_echo*_MEGRE.json",
)
parser.add_option(
    "-o", "--outputDirectory", dest="outputDirectory", help="Output directory"
)
parser.add_option("-v", "--verbose", dest="verbose", default=True, help="verbose")

(options, args) = parser.parse_args()


################################################################################
# 1) Value Extractor
################################################################################

runT2StarRelaxometryMapper(
    options.mgedirectory, options.MGEFilenames, options.outputDirectory, options.verbose
)
