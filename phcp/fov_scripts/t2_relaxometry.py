"""Single Compartment T2 Relaxometry Mapper"""

import glob
import json
import logging
import os
import sys

import numpy as np
import nibabel as ni

from phcp.gkg import (
    gkg_convert_gis_to_nifti,
    run_gkg_command,
    run_gkg_GetMask,
    run_gkg_SubVolume,
)


logger = logging.getLogger(__name__)


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


def createCommand_SingleCompartmentRelaxometryMapper(
    fileNameMSME,
    fileNameMask,
    fileNameTEValues,
    fileNameProtonDensity,
    fileNameT2,
    fileNameFittedMSME,
    verbose,
):
    command = [
        "GkgExecuteCommand",
        "SingleCompartmentRelaxometryMapper",
        "-i",
        fileNameMSME,
        "-m",
        fileNameMask,
        "-t",
        "t2-mapping-msme",
        "-optimizerParameters",
        "50000",  # <P1> : NLP maximum iteration count
        "-optimizerParameters",
        "0.001",  # <P2> : NLP test size
        "-optimizerParameters",
        "0",  # <P3> : 0->do not apply MCMC 1->apply MCMC
        "-optimizerParameters",
        "2000",  # <P4> : MCMC burnin count
        "-optimizerParameters",
        "300",  # <P5> : MCMC sample count
        "-optimizerParameters",
        "50",  # <P6> : MCMC interval count
        "-optimizerParameters",
        "10000",  # <P7> : MCMC maximum iteration count
        "-scalarParameters",
        "0.75",
        "-scalarParameters",
        "50",
        "-scalarParameters",
        "0",
        "-scalarParameters",
        "0",
        "-scalarParameters",
        "8",
        "-scalarParameters",
        "180",
        "-scalarParameters",
        "1",
        "-scalarParameters",
        "20",
        "-scalarParameters",
        "5",
        "-stringParameters",
        fileNameTEValues,
        "-op",
        fileNameProtonDensity,
        "-ot",
        fileNameT2,
        "-f",
        fileNameFittedMSME,
        "-ascii",
        "False",
        "-format",
        "gis",
        "-verbose",
        str(verbose),
    ]
    return command


def write_EchoTimes(fileNameTEValues, MSMEFilenames):
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
    logger.info("SINGLE COMPARTMENT T2 RELAXOMETRY MAPPER")

    fileNameMSME = os.path.join(subjectDirectoryGisConversion, "t2-msme.nii.gz")
    fileNameMask = os.path.join(subjectDirectoryGisConversion, "mask.nii.gz")
    fileNameTEValues = os.path.join(outputDirectory, "te.txt")
    t2FileName = os.path.join(subjectDirectoryGisConversion, "t2_extracted.nii.gz")

    if not (os.path.exists(t2FileName)):
        run_gkg_SubVolume(
            [
                "-i",
                fileNameMSME,
                "-o",
                t2FileName,
                "-tIndices",
                "0",
                "-verbose",
                "true",
            ],
            output_dirs=[subjectDirectoryGisConversion],
        )

    if not (os.path.exists(fileNameMask)):
        run_gkg_GetMask(
            ["-i", fileNameMSME, "-o", fileNameMask, "-a", "2", "-verbose"],
            output_dirs=[subjectDirectoryGisConversion],
        )

    write_EchoTimes(fileNameTEValues, MSMEFilenames)

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
    run_gkg_command(
        command,
        gkg_container_version="2022-12-20",
        input_dirs=[subjectDirectoryGisConversion],
        output_dirs=[outputDirectory],
    )

    fileNameT2nifti = os.path.join(outputDirectory, "T2.nii.gz")
    gkg_convert_gis_to_nifti(fileNameT2, fileNameT2nifti, verbose=True)

    ni_im = create_ConfidenceMap(
        fileNameMSME, fileNameFittedMSME, fileNameMask, fileNameT2nifti
    )
    ni.save(ni_im, os.path.join(outputDirectory, "T2ConfidenceMap.nii.gz"))
    return None


def parse_command_line(argv):
    """Parse the script's command line."""
    import argparse

    parser = argparse.ArgumentParser(
        description="reconstruct T2 relaxometry based on the Multi Echo Spin Echo method",
    )
    parser.add_argument(
        "-i", "--input", dest="t2msmedirectory", help="T2-MSME Directory"
    )
    parser.add_argument(
        "-m",
        "--msme",
        dest="MSMEFilenames",
        help="String for glob search of MSME volumes. Ex: /phcp/rawdata/sub/ses/anat/sub_ses_echo*_MESE.json",
    )
    parser.add_argument(
        "-o", "--outputDirectory", dest="outputDirectory", help="Output directory"
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=True,
        help="print detailed information during the fit",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        dest="verbose",
        action="store_false",
        help="do not print detailed information during the fit",
    )

    args = parser.parse_args()
    return args


def main(argv=sys.argv):
    """The script's entry point."""
    logging.basicConfig(level=logging.INFO)
    args = parse_command_line(argv)
    return (
        runT2RelaxometryMapper(
            args.t2msmedirectory,
            args.MSMEFilenames,
            args.outputDirectory,
            args.verbose,
        )
        or 0
    )


if __name__ == "__main__":
    sys.exit(main())
