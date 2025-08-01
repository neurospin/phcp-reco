#!/usr/bin/env python3
import itertools
import logging
import os
import subprocess
import sys
from pathlib import Path

import ants

import phcp.gkg as gkg

logger = logging.getLogger(__name__)


def convert_all_gis_files_to_nifti(
    moving_directory, fixed_directory, verbose: bool
) -> None:
    """Convert GIS files to NIFTI format."""
    GIS_file_to_convert_iterator = itertools.chain(
        moving_directory.glob("*.ima"), fixed_directory.glob("*.ima")
    )

    for gis_path in GIS_file_to_convert_iterator:
        nii_path = gis_path.with_suffix(".nii")
        niigz_path = gis_path.with_suffix(".nii.gz")

        if not nii_path.exists() and not niigz_path.exists():
            verbose and logger.info(f"Converting {gis_path.name} to NIFTI format...")

            gkg.gkg_convert_gis_to_nifti(
                gis_path.as_posix(),
                gis_path.with_suffix(".nii.gz").as_posix(),
                verbose=verbose,
            )

            verbose and logger.info(f"Conversion of {gis_path.name} completed.")
        else:
            verbose and logger.info(f"{gis_path.name} already exists in NIFTI format.")


def apply_n4_bias_correction_to_all_nifti_files(
    moving_directory, fixed_directory, verbose: bool
) -> None:
    """Apply N4 bias field correction to all NIFTI files in moving and fixed directories."""
    all_nifti_files_in_moving_and_fixed = itertools.chain(
        moving_directory.glob("*.nii"),
        moving_directory.glob("*.nii.gz"),
        fixed_directory.glob("*.nii"),
        fixed_directory.glob("*.nii.gz"),
    )

    for nifti_file in all_nifti_files_in_moving_and_fixed:
        if (
            not nifti_file.name.endswith("_N4.nii")
            and not nifti_file.name.endswith("_N4.nii.gz")
            and not Path(nifti_file.as_posix().split(".")[0] + "_N4.nii.gz").exists()
        ):
            verbose and logger.info(
                f"Applying N4 bias field correction to {nifti_file.name}..."
            )
            ants_volume = ants.image_read(nifti_file.as_posix())
            ants_bias_field = ants.n4_bias_field_correction(
                ants_volume,
                spline_param=20,
                rescale_intensities=True,
                return_bias_field=True,
            )
            ants.image_write(
                ants_volume / ants_bias_field,
                (
                    nifti_file.parent
                    / (
                        nifti_file.name.removesuffix(".nii.gz")
                        if nifti_file.name.endswith(".nii.gz")
                        else nifti_file.stem
                    )
                ).as_posix()
                + "_N4.nii.gz",
            )
            verbose and logger.info("N4 bias field correction applied.")
        else:
            verbose and logger.info(
                f"Skipping {nifti_file.name} as it already has N4 applied."
            )


def log_registration_info(moving_filename: str, fixed_filename: str) -> None:
    logger.info(
        "\n"
        "======================================================\n"
        "    Rigid / Affine / SyN personalized REGISTRATION :\n"
        f"        --> Moving file : {moving_filename}\n"
        f"        --> Fixed file : {fixed_filename}\n"
        "======================================================"
    )


def run_ants_registration(
    ref_filename: str, moving_filename: str, out_prefix: str
) -> None:
    """Run ANTs registration command."""
    ants_command = [
        "antsRegistration",
        "-d",
        "3",
        "-r",
        f"[{ref_filename},{moving_filename},1]",
        "-m",
        f"MI[{ref_filename},{moving_filename},1,8,Regular,0.5]",
        "-t",
        "Rigid[2]",
        "-c",
        "[500x500x500,1e-6,8]",
        "-s",
        "4x2x1vox",
        "-f",
        "4x2x2",
        "-m",
        f"MI[{ref_filename},{moving_filename},1,16,Regular,0.8]",
        "-t",
        "Affine[1]",
        "-c",
        "[500x500x500,1e-6,8]",
        "-s",
        "4x2x1vox",
        "-f",
        "4x2x2",
        "-m",
        f"mattes[{ref_filename},{moving_filename},1,32,Regular,0.8]",
        "-t",
        "SyN[0.5,3,1]",
        "-c",
        "[1200x1200x1000,1e-6,8]",
        "-s",
        "4x2x1vox",
        "-f",
        "4x2x2",
        "-o",
        f"[{out_prefix},{out_prefix}Regis.nii.gz,{out_prefix}INV_Regis.nii.gz]",
        "-v",
        "1",
    ]
    subprocess.run(ants_command, check=True)


def compute_and_save_jacobian(
    fixed_image: ants.ANTsImage, warp_filename: str, output_filename: str
) -> None:
    """Compute and save the Jacobian determinant image."""
    jacobian_image = ants.create_jacobian_determinant_image(
        fixed_image, warp_filename, do_log=False
    )
    ants.image_write(jacobian_image, output_filename)


def run_registrations_intrafov(working_directory: str | Path, verbose: bool) -> None:
    moving_folder = Path(working_directory) / "moving"
    fixed_folder = Path(working_directory) / "fixed"

    logger.info("Starting registration part...")

    # FIXED VOLUME INITIALIZATION
    t2w_recN4_filename = list(fixed_folder.glob("*N4.nii.gz"))[0].as_posix()
    all_moving_filename = moving_folder.glob("*N4.nii.gz")
    transform_files_folder = Path(working_directory) / "transform_files"
    transform_files_folder.mkdir(exist_ok=True)
    # REGISTRATION PART
    for file in all_moving_filename:
        name_moving_image = file.name.split(".")[0].split("_")[0]

        warp_filename = transform_files_folder / f"{name_moving_image}_To_T2w_1Warp.mat"
        affine_filename = (
            transform_files_folder / f"{name_moving_image}_To_T2w_0GenericAffine.mat"
        )

        if not (os.path.exists(warp_filename)) and not (
            os.path.exists(affine_filename)
        ):
            verbose and log_registration_info(
                file.name, t2w_recN4_filename.split("/")[-1]
            )

            out = transform_files_folder / f"{name_moving_image}_To_T2w_"

            run_ants_registration(t2w_recN4_filename, file.as_posix(), out.as_posix())
            verbose and logger.info("ANTS registration command executed.")

            if warp_filename.exists():
                verbose and logger.info("Jacobian determinant computation started...")
                ref_init_space = ants.image_read(t2w_recN4_filename)
                compute_and_save_jacobian(
                    ref_init_space,
                    warp_filename.as_posix(),
                    (
                        transform_files_folder
                        / f"{name_moving_image}_To_T2w_Jacobian.nii.gz"
                    ).as_posix(),
                )
                verbose and logger.info("Jacobian determinant computation completed.")


def run_pipeline(working_directory: str | Path, verbose: bool) -> None:
    moving_folder = Path(working_directory) / "moving"
    fixed_folder = Path(working_directory) / "fixed"

    logger.info("Starting data formatting and GIS to NIFTI conversion...")
    convert_all_gis_files_to_nifti(moving_folder, fixed_folder, verbose)
    logger.info("GIS to NIfTI conversion : done.")

    logger.info("Applying N4 bias field correction to all NIFTI files...")
    apply_n4_bias_correction_to_all_nifti_files(moving_folder, fixed_folder, verbose)
    logger.info("N4 bias field correction : done.")

    run_registrations_intrafov(working_directory, verbose)
    logger.info("Intrafov registration pipeline completed.")


def parse_command_line(argv):
    """Parse the script's command line."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Intrafov registration script for PHCP data."
    )
    parser.add_argument(
        "-i",
        "--workingDirectory",
        required=True,
        help="Working directory such as DirectoryName/01-Materials   /02-..   /03-..",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        required=False,
        default=False,
        help="do not print detailed",
    )

    args = parser.parse_args()
    return args


def main(argv=sys.argv):
    """The script's entry point."""
    logging.basicConfig(level=logging.INFO)
    args = parse_command_line(argv)
    return run_pipeline(args.workingDirectory, args.verbose) or 0


if __name__ == "__main__":
    sys.exit(main())
