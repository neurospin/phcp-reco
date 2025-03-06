import logging
import os
import pathlib
import shutil
import sys

import nibabel
import numpy

logger = logging.getLogger(__name__)


REORIENT_MATRIX = numpy.array([[1, 0, 0, 0], [0, 0, -1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])


def flip_headers(
    sesDirectory: str | os.PathLike, outputDirectory: str | os.PathLike
) -> None:
    input_session_dir = pathlib.Path(sesDirectory).resolve()
    output_bids_dir = pathlib.Path(outputDirectory)

    if (
        not input_session_dir.is_dir()
        or not input_session_dir.name.startswith("ses-")
        or not input_session_dir.parent.name.startswith("sub-")
    ):
        logger.fatal("Expecting a BIDS Session directory as first argument, aborting.")
        return 1

    # re-create sub-*/ses-*
    output_session_dir = (
        output_bids_dir / input_session_dir.parent.name / input_session_dir.name
    )
    output_session_dir.parent.mkdir(exist_ok=True)
    output_session_dir.mkdir(exist_ok=True)

    logger.info("Processing session %s", input_session_dir.name)
    for datatype_input_dir in input_session_dir.iterdir():
        if not datatype_input_dir.is_dir():
            logger.info("Skipping non-directory entry: %s", datatype_input_dir)
            continue
        datatype = datatype_input_dir.name
        output_datatype_dir = output_session_dir / datatype
        output_datatype_dir.mkdir(exist_ok=True)

        logger.info("Processing datatype %s", datatype)
        for input_volume_path in datatype_input_dir.glob("*.nii*"):
            input_volume_basename = input_volume_path.name.split(".")[0]
            output_volume_path = output_datatype_dir / input_volume_path.name
            input_json_path = datatype_input_dir / (input_volume_basename + ".json")
            output_json_path = output_datatype_dir / input_json_path.name

            if output_volume_path.exists() and output_json_path.exists():
                logger.info("Skipping %s which already exists", input_volume_path.name)
            else:
                logger.info("Processing %s...", input_volume_path.name)
                input_img = nibabel.load(input_volume_path)
                output_header = input_img.header.copy()
                output_header.set_dim_info(None, None, None)
                output_img = nibabel.Nifti1Image(
                    input_img.dataobj, REORIENT_MATRIX @ input_img.affine, output_header
                )
                nibabel.save(output_img, output_volume_path)
                shutil.copy(input_json_path, output_json_path)


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
        required=True,
        help="input session directory, e.g. fov/rawdata/sub-${sub}/ses-${ses}",
    )
    parser.add_argument(
        "-o",
        "--outputDirectory",
        required=True,
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
