#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 2025

@author: la272118
"""

import os
import json
import optparse
import subprocess
import numpy as np
import nibabel as ni
import glob


""" Miscellaneous algorithms """


def print_message(message):
    print("==================================")
    print(message)
    print("==================================")


def ConversionGisToNifti(InputGisFilename, OutputNiftiFilename):
    subprocess.run(
        [
            "singularity",
            "exec",
            "--bind",
            "/neurospin:/neurospin:rw",
            "/neurospin/phcp/code/gkg/2022-12-20_gkg/2022-12-20_gkg.sif",
            "GkgExecuteCommand",
            "Gis2NiftiConverter",
            "-i",
            InputGisFilename,
            "-o",
            OutputNiftiFilename,
            "-verbose",
        ]
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


""" Single Compartment T2 Relaxometry Mapper """


def createCommand_SubVolume_gkgMethod(fileNameMSME, t2FileName):
    command = (
        "GkgExecuteCommand SubVolume"
        + " -i "
        + fileNameMSME
        + " -o "
        + t2FileName
        + " -tIndices "
        + str(0)
        + " -verbose true"
    )
    return command


def createCommand_Mask_gkgMethod(fileNameMSME, fileNameMask):
    command = (
        "GkgExecuteCommand GetMask"
        + " -i "
        + fileNameMSME
        + " -o "
        + fileNameMask
        + " -a 2 "
        + " -verbose "
    )
    return command


def createCommand_SingleCompartmentRelaxometryMapper(
    fileNameMSME,
    fileNameMask,
    fileNameTEValues,
    fileNameProtonDensity,
    fileNameT2,
    fileNameFittedMSME,
    verbose,
):
    # optimizerParameters :: <P1> : NLP maximum iteration count  <P2> : NLP test size <P3> : 0->do not apply MCMC 1->apply MCMC <P4> : MCMC burnin count <P5> : MCMC sample count <P6> : MCMC interval count <P7> : MCMC maximum iteration coun
    # ' -optimizerParameters 50000 0.001 0 2000 300 50 10000' +\  ' -scalarParameters 0.75 50 0 0 8 180 1 20 5'  +\
    command = (
        "GkgExecuteCommand SingleCompartmentRelaxometryMapper"
        + " -i %s" % fileNameMSME
        + " -m %s" % fileNameMask
        + " -t t2-mapping-msme"
        + " -optimizerParameters %s " % 50000
        + " -optimizerParameters %s " % 0.001
        + " -optimizerParameters %s" % 0
        + " -optimizerParameters %s" % 2000
        + " -optimizerParameters %s" % 300
        + " -optimizerParameters %s" % 50
        + " -optimizerParameters %s" % 10000
        + " -scalarParameters %s" % 0.75
        + " -scalarParameters %s" % 50
        + " -scalarParameters %s" % 0
        + " -scalarParameters %s" % 0
        + " -scalarParameters %s" % 8
        + " -scalarParameters %s" % 180
        + " -scalarParameters %s" % 1
        + " -scalarParameters %s" % 20
        + " -scalarParameters %s" % 5
        + " -stringParameters %s" % fileNameTEValues
        + " -op %s" % fileNameProtonDensity
        + " -ot %s" % fileNameT2
        + " -f %s" % fileNameFittedMSME
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


def write_EchotTime(fileNameTEValues, MSMEFilenames):
    with open(fileNameTEValues, "w") as file:
        for metadata_filename in sorted(glob.iglob(MSMEFilenames)):
            with open(metadata_filename) as metadata_file:
                metadata = json.load(metadata_file)
            file.write(str(metadata["EchoTime"] * 1000) + "\n")
    return None


def create_ConfidenceMap(
    fileNameMSME, fileNameFittedMSME, fileNameMask, fileNameT2nifti
):
    arr_fittedmsme = create_ArrayFromGisFilename(fileNameFittedMSME)
    arr_fittedmsme = arrange_ArrayToFlipHeaderFormat(arr_fittedmsme)

    meta_t2msme = ni.load(fileNameMSME)
    arr_t2msme = meta_t2msme.get_fdata()

    meta_mask_t2 = ni.load(fileNameMask)
    arr_maskt2 = meta_mask_t2.get_fdata()
    arr_maskt2 = arrange_ArrayToFlipHeaderFormat(arr_maskt2)

    diff = arr_fittedmsme - arr_t2msme
    std_diff = np.std(diff[:, :, :, :10], axis=-1) * arr_maskt2

    meta_t2nifti = ni.load(fileNameT2nifti)
    ni_im = ni.Nifti1Image(
        dearrange_ArrayToFlipHeaderFormat(std_diff),
        meta_t2nifti.affine,
        meta_t2nifti.header,
    )
    return ni_im


def runT2RelaxometryMapper(
    subjectDirectoryGisConversion, MSMEFilenames, outputDirectory, verbose
):
    if verbose:
        print_message("SINGLE COMPARTMENT T2 RELAXOMETRY MAPPER")

    fileNameMSME = os.path.join(subjectDirectoryGisConversion, "t2-msme.nii.gz")
    fileNameMask = os.path.join(subjectDirectoryGisConversion, "mask.nii.gz")
    fileNameTEValues = os.path.join(outputDirectory, "te.txt")
    t2FileName = os.path.join(subjectDirectoryGisConversion, "t2_extracted.nii.gz")

    if not (os.path.exists(t2FileName)):
        command = createCommand_SubVolume_gkgMethod(fileNameMSME, t2FileName)
        run_command(command)

    if not (os.path.exists(fileNameMask)):
        command = createCommand_Mask_gkgMethod(fileNameMSME, fileNameMask)
        run_command(command)

    write_EchotTime(fileNameTEValues, MSMEFilenames)

    fileNameProtonDensity = os.path.join(outputDirectory, "proton-density.ima")
    fileNameFittedMSME = os.path.join(outputDirectory, "fitted-msme.ima")
    fileNameT2 = os.path.join(outputDirectory, "T2.ima")

    command = createCommand_SingleCompartmentRelaxometryMapper(
        fileNameMSME,
        fileNameMask,
        fileNameTEValues,
        fileNameProtonDensity,
        fileNameT2,
        fileNameFittedMSME,
        verbose,
    )
    run_command(command)

    fileNameT2nifti = os.path.join(outputDirectory, "T2.nii.gz")
    ConversionGisToNifti(fileNameT2, fileNameT2nifti)

    ni_im = create_ConfidenceMap(
        fileNameMSME, fileNameFittedMSME, fileNameMask, fileNameT2nifti
    )
    ni.save(ni_im, os.path.join(outputDirectory, "T2ConfidenceMap.nii.gz"))
    return None


################################################################################
# parser to get option(s)
################################################################################

parser = optparse.OptionParser()
parser.add_option("-i", "--input", dest="t2msmedirectory", help="T2-MSME Directory")
parser.add_option(
    "-m",
    "--msme",
    dest="MSMEFilenames",
    help="String for glob research of MSME volumes. Ex : /phcp/rawdata/sub/ses/anat/sub_ses_echo*_MESE.json",
)
parser.add_option(
    "-o", "--outputDirectory", dest="outputDirectory", help="Output directory"
)
parser.add_option("-v", "--verbose", dest="verbose", default=True, help="verbose")

(options, args) = parser.parse_args()


################################################################################
# 1) Value Extractor
################################################################################

runT2RelaxometryMapper(
    options.t2msmedirectory,
    options.MSMEFilenames,
    options.outputDirectory,
    options.verbose,
)
