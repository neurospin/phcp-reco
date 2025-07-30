#!/usr/bin/env python3
import glob
import itertools
import logging
import os
import sys
from pathlib import Path

import ants
import numpy as np

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


def Data_formating(working_directory: str | Path, verbose: bool) -> None:
    moving_folder = Path(working_directory) / "moving"
    fixed_folder = Path(working_directory) / "fixed"

    print(f"Moving folder: {moving_folder}")
    verbose and logger.info("Starting data formatting and GIS to NIFTI conversion...")
    convert_all_gis_files_to_nifti(moving_folder, fixed_folder, verbose)
    verbose and logger.info("GIS to NIfTI conversion : done.")

    verbose and logger.info("Applying N4 bias field correction to all NIFTI files...")
    apply_n4_bias_correction_to_all_nifti_files(moving_folder, fixed_folder, verbose)
    verbose and logger.info("N4 bias field correction : done.")


def Registration_Type1(MainDirectory, verbose):
    MaterialsDirectory = os.path.join(MainDirectory, "01-Materials")
    RegistrationDirectory = os.path.join(MainDirectory, "02-Registration")
    verbose = bool(verbose)

    with open(RegistrationDirectory + "/Jacobian_resume.txt", "a") as f:
        Jacobian_textfile = f.readlines()

    if verbose:
        print("======================================================")
        print("REGISTRATION Type 1: ")
        print("======================================================")

    # FIXED VOLUME INITIALIZATION
    Path_FixedFilename = (
        glob.glob(os.path.join(MaterialsDirectory, "Fixed") + "/*N4.nii.gz")
        + glob.glob(os.path.join(MaterialsDirectory, "Fixed") + "/*N4.nii")
    )[0]

    Path_MovingFile = os.path.join(MaterialsDirectory, "Moving")

    FileToRegis = glob.glob(os.path.join(Path_MovingFile, "*N4.nii")) + glob.glob(
        os.path.join(Path_MovingFile, "*N4.nii.gz")
    )

    os.makedirs(os.path.join(RegistrationDirectory, "TransformFiles"), exist_ok=True)

    Filename_Fixed = Path_FixedFilename.split("/")[-1].split(".")[0].split("_")[-2]

    Path_TransformFiles = os.path.join(RegistrationDirectory, "TransformFiles")

    Fixed = ants.image_read(Path_FixedFilename)

    # REGISTRATION PART
    for file in FileToRegis:
        Filename_Moving = file.split("/")[-1].split(".")[0].split("_")[-2]

        Path_WarpFilename = os.path.join(
            Path_TransformFiles,
            Filename_Moving + "To" + Filename_Fixed + "3Warp.nii.gz",
        )
        Path_AffineFilename = os.path.join(
            Path_TransformFiles, Filename_Moving + "To" + Filename_Fixed + "2Affine.mat"
        )

        Path_WarpFilename_2 = os.path.join(
            Path_TransformFiles,
            Filename_Moving + "To" + Filename_Fixed + "1Warp.nii.gz",
        )
        Path_AffineFilename_2 = os.path.join(
            Path_TransformFiles,
            Filename_Moving + "To" + Filename_Fixed + "0GenericAffine.mat",
        )

        if (
            not (os.path.exists(Path_WarpFilename))
            and not (os.path.exists(Path_AffineFilename))
        ) and (
            not (os.path.exists(Path_WarpFilename_2))
            and not (os.path.exists(Path_AffineFilename_2))
        ):
            if verbose:
                print("======================================================")
                print("    Rigid / Affine / SyN personalized REGISTRATION :")
                print("        --> Moving file : " + file.split("/")[-1])
                print("        --> Fixed file : " + Path_FixedFilename.split("/")[-1])
                print("======================================================")

            # grad = 0.25
            out = os.path.join(
                Path_TransformFiles, Filename_Moving + "To" + Filename_Fixed
            )
            os.system(
                "antsRegistration "
                f"-d {3} "
                f"-r [{Path_FixedFilename}, {file}, 1] "
                f"-m MI[{Path_FixedFilename},{file},{1},{8}, Regular, {0.5}] "
                f"-t Rigid[{2}] "
                f"-c [{500}x{500}x{500},{1e-6},{8}] "
                f"-s {4}x{2}x{1}vox "
                f"-f {4}x{2}x{2} "
                f"-m MI[{Path_FixedFilename},{file},{1},{16}, Regular, {0.8}] "
                f"-t Affine[{1}] "
                f"-c [{500}x{500}x{500},{1e-6},{8}] "
                f"-s {4}x{2}x{1}vox "
                f"-f {4}x{2}x{2} "
                f"-m mattes[{Path_FixedFilename},{file},{1},{32}, Regular, {0.8}] "
                f"-t SyN[{0.5},{3},{1}] "
                f"-c [{1200}x{1200}x{1000},{1e-6},{8}] "
                f"-s {4}x{2}x{1}vox "
                f"-f {4}x{2}x{2} "
                f"-o [{out},{out}_Regis.nii.gz,{out}INV_Regis.nii.gz] "
                " -x [NA,NA] "
                "-v 1"
            )

            if verbose:
                print("REGISTRATION DONE")

            ##TRANSFORM FILES SAVE
            if verbose:
                print("======================================================")
                print("    SAVING PART ")
                print("======================================================")

            # WarpFile = ants.image_read(Path_WarpFilename)
            # ants.image_write(WarpFile, Path_WarpFilename)

            if not (os.path.exists(Path_WarpFilename)) and not (
                os.path.exists(Path_AffineFilename)
            ):
                Path_WarpFilename = Path_WarpFilename_2
                Path_AffineFilename = Path_AffineFilename_2

            if os.path.exists(Path_WarpFilename) or os.path.exists(Path_WarpFilename_2):
                jac = ants.create_jacobian_determinant_image(
                    Fixed, Path_WarpFilename, do_log=False
                )
                Path_JacFilename = os.path.join(
                    Path_TransformFiles,
                    Filename_Moving + "To" + Filename_Fixed + "_Jacobian.nii.gz",
                )
                ants.image_write(jac, Path_JacFilename)

                Affine = ants.read_transform(Path_AffineFilename)
                # ants.write_transform(Affine,  Path_AffineFilename)

                Affine_matrix = (Affine.parameters)[0:9].reshape((3, 3))

                Jacobian_textfile.write(
                    "REGISTRATION type 1\n "
                    "\n          moving --> "
                    + str(file.split("/")[-1])
                    + "\n          fixed --> "
                    + Filename_Fixed
                    + "\n    Mean Jacobian Warp File = "
                    + str(jac.mean())
                    + "\n    Jacobian Affine matrix = "
                    + str(np.linalg.det(Affine_matrix))
                    + "\n \n"
                )
            else:
                Affine = ants.read_transform(Path_AffineFilename)
                # ants.write_transform(Affine,  Path_AffineFilename)

                Affine_matrix = (Affine.parameters)[0:9].reshape((3, 3))
                Jacobian_textfile.write(
                    "REGISTRATION type 1\n "
                    "\n          moving --> "
                    + str(file.split("/")[-1])
                    + "\n          fixed --> "
                    + Filename_Fixed
                    + "\n    Jacobian Affine matrix = "
                    + str(np.linalg.det(Affine_matrix))
                    + "\n \n"
                )

            if verbose:
                print("SAVING PART DONE")
        else:
            print("======================================================")
            print("    TRANSFORM FILES ALREADY EXIST")
            print("======================================================")


def run_pipeline(MainDirectory, verbose):
    Data_formating(MainDirectory, verbose)
    # Registration_Type1(MainDirectory, verbose)


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
