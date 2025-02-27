import glob
import logging
import os
import shutil
import sys

import nibabel as ni
import numpy as np


logger = logging.getLogger(__name__)


def flip_headers(sesDirectory, outputDirectory):
    SessionOutputDir = os.path.join(outputDirectory, sesDirectory.split("/")[-2])
    if not os.path.exists(SessionOutputDir):
        os.mkdir(SessionOutputDir)

    modalitiesDirectory = glob.iglob(os.path.join(sesDirectory, "*"))
    session = sesDirectory.split("/")[-2]
    logger.info("Processing session %s", session)
    for modality in modalitiesDirectory:
        ModOutputDir = os.path.join(outputDirectory, session, modality.split("/")[-1])
        if not os.path.exists(ModOutputDir):
            os.mkdir(ModOutputDir)

        filenamesDirectory = glob.iglob(os.path.join(modality, "*.nii.gz"))
        logger.info("Processing modality %s", modality.split("/")[-1])
        for volume in filenamesDirectory:
            logger.info("Converting %s", volume.split("/")[-1])
            outputFilename = os.path.join(
                outputDirectory,
                session,
                modality.split("/")[-1],
                volume.split("/")[-1],
            )
            if not os.path.exists(outputFilename):
                meta = ni.load(volume)
                matRotation = np.array(
                    [[1, 0, 0, 0], [0, 0, -1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]
                )
                ni_im = ni.Nifti1Image(
                    meta.get_fdata(), np.dot(matRotation, meta.affine), meta.header
                )
                ni.save(ni_im, outputFilename)

            JsonFilename = os.path.join(
                sesDirectory,
                modality.split("/")[-1],
                (volume.split("/")[-1]).split(".")[0] + ".json",
            )
            outputJsonFilename = os.path.join(
                outputDirectory,
                session,
                modality.split("/")[-1],
                (volume.split("/")[-1]).split(".")[0] + ".json",
            )
            if os.path.exists(JsonFilename) and not os.path.exists(outputJsonFilename):
                shutil.copy(JsonFilename, outputJsonFilename)


def parse_command_line(argv):
    """Parse the script's command line."""
    import argparse

    parser = argparse.ArgumentParser(
        description='Rewrite the Nifti header transforms for a whole session, so that the axes of the new "scanner-based" '
        "referential correspond to sensible anatomical orientations. This is done so that visualizing the images the usual "
        "viewers will display the brain in a sensible orientation. We cannot do this at acquisition time because the "
        "console software does not allow us to input the real orientation of the tissue block, so we default to choosing "
        '"head first supine" (HFS).',
    )
    parser.add_argument(
        "-i",
        "--sessionDirectory",
        help="input session directory, e.g. fov/rawdata/sub-${sub}/ses-${ses}",
    )
    parser.add_argument(
        "-o",
        "--outputDirectory",
        help="top-level output BIDS directory, e.g. fov/headerfliprawdata",
    )

    args = parser.parse_args()
    return args


def main(argv=sys.argv):
    """The script's entry point."""
    logging.basicConfig(level=logging.INFO)
    args = parse_command_line(argv)
    return flip_headers(args.sessionDirectory, args.outputDirectory) or 0


if __name__ == "__main__":
    sys.exit(main())
