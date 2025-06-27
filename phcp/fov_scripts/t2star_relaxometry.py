"""Single Compartment T2Star Relaxometry Mapper"""

import glob
import json
import logging
import os
import sys

import nibabel as ni
import numpy as np

from phcp.gkg import (
    run_gkg_command,
    run_gkg_GetMask,
    run_gkg_SubVolume,
)
from phcp.image import nibabel_orient_like

logger = logging.getLogger(__name__)


def createCommand_T2starSingleCompartmentRelaxometryMapper(
    fileNameMGE,
    fileNameMask,
    fileNameTEValues,
    fileNameProtonDensity,
    fileNameT2star,
    fileNameFittedMGE,
    verbose,
):
    command = [
        "GkgExecuteCommand",
        "SingleCompartmentRelaxometryMapper",
        "-i",
        fileNameMGE,
        "-m",
        fileNameMask,
        "-t",
        "t2-star-mapping-mgre",
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
        "0.3",
        "-scalarParameters",
        "30",
        "-scalarParameters",
        "0",
        "-scalarParameters",
        "0",
        "-scalarParameters",
        "8",
        "-scalarParameters",
        "100",
        "-scalarParameters",
        "1",
        "-scalarParameters",
        "10",
        "-scalarParameters",
        "2",
        "-stringParameters",
        fileNameTEValues,
        "-op",
        fileNameProtonDensity,
        "-ot",
        fileNameT2star,
        "-f",
        fileNameFittedMGE,
        "-ascii",
        "False",
        "-format",
        "nifti",
        "-verbose",
        str(verbose),
    ]
    return command


def write_EchoTimes(fileNameTEValues, MGEFilenames):
    with open(fileNameTEValues, "w") as file:
        for metadata_filename in sorted(glob.iglob(MGEFilenames)):
            with open(metadata_filename) as metadata_file:
                metadata = json.load(metadata_file)
            file.write(str(metadata["EchoTime"] * 1000) + "\n")
    return None


def create_ConfidenceMap(
    fileNameMGE, fileNameFittedMGE, fileNameMask, fileNameT2starnifti
):
    meta_fittedmge = ni.load(fileNameFittedMGE)
    arr_fittedmge = meta_fittedmge.get_fdata()

    meta_mge = nibabel_orient_like(ni.load(fileNameMGE), meta_fittedmge)
    arr_mge = meta_mge.get_fdata()

    meta_mask = nibabel_orient_like(ni.load(fileNameMask), meta_fittedmge)
    arr_mask = meta_mask.get_fdata()

    diff = arr_fittedmge - arr_mge
    std_diff = (
        np.std(diff[:, :, :, :6], axis=-1) * arr_mask
    )  # FIXME: magic parameter should be documented

    ni_im = ni.Nifti1Image(
        std_diff,
        meta_fittedmge.affine,
        meta_fittedmge.header,
    )
    return ni_im


def runT2StarRelaxometryMapper(
    subjectDirectoryGisConversion, MGEFilenames, outputDirectory, verbose
):
    logger.info("SINGLE COMPARTMENT T2Star RELAXOMETRY MAPPER")

    # TODO: start from the separate BIDS MGE files (incorporate the concatenation in this script)

    fileNameMGE = os.path.join(subjectDirectoryGisConversion, "t2star-mge.nii.gz")
    fileNameMask = os.path.join(subjectDirectoryGisConversion, "mask.nii.gz")
    fileNameTEValues = os.path.join(outputDirectory, "te.txt")
    t2starFileName = os.path.join(
        subjectDirectoryGisConversion, "t2star_extracted.nii.gz"
    )

    if not (os.path.exists(t2starFileName)):
        run_gkg_SubVolume(
            [
                "-i",
                fileNameMGE,
                "-o",
                t2starFileName,
                "-tIndices",
                "0",
                "-verbose",
                "true",
            ],
            output_dirs=[subjectDirectoryGisConversion],
        )

    if not (os.path.exists(fileNameMask)):
        run_gkg_GetMask(
            ["-i", fileNameMGE, "-o", fileNameMask, "-a", "2", "-verbose"],
            output_dirs=[subjectDirectoryGisConversion],
        )

    write_EchoTimes(fileNameTEValues, MGEFilenames)

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
    run_gkg_command(
        command,
        gkg_container_version="2022-12-20",
        input_dirs=[subjectDirectoryGisConversion],
        output_dirs=[outputDirectory],
    )
    # FIXME: outputs have incorrect qform (no translation despite seemingly
    # correct rotation) and grossly incorrect sform (no scale, no translation)

    ni_im = create_ConfidenceMap(
        fileNameMGE, fileNameFittedMGE, fileNameMask, fileNameT2star
    )
    ni.save(ni_im, os.path.join(outputDirectory, "T2starConfidenceMap.nii.gz"))
    return None


def parse_command_line(argv):
    """Parse the script's command line."""
    import argparse

    parser = argparse.ArgumentParser(
        description="reconstruct T2* relaxometry based on the Multi Echo GRadient Echo method",
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Directory with the input t2star-mge.nii.gz, "
        "typically fov/derivatives/T2starmapping/sub-${sub}/ses-${ses}/01-Materials",
    )
    parser.add_argument(
        "-m",
        "--mge",
        dest="MGEFilenames",
        required=True,
        help="String for glob search of MGE volume metadata, typically fov/rawdata/sub-${sub}/ses-${ses}/anat/'*_MEGRE.json'",
    )
    parser.add_argument(
        "-o",
        "--outputDirectory",
        required=True,
        help="Output directory, typically fov/derivatives/T2starmapping/sub-${sub}/ses-${ses}/02-Results",
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
        runT2StarRelaxometryMapper(
            args.input, args.MGEFilenames, args.outputDirectory, args.verbose
        )
        or 0
    )


if __name__ == "__main__":
    sys.exit(main())
