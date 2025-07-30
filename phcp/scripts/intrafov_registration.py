#!/usr/bin/env python3
"""
Created on Fri Aug  4 15:53:57 2023

@author: la272118

Objective: Registration of all modalities through sessions of acquisitions.
"""

import glob
import optparse
import os

import ants
import numpy as np

# from multiprocessing import Process


def ConversionGisToNifti(InputGisFilename, OutputNiftiFilename):
    os.system(
        "singularity exec --bind /neurospin:/neurospin:rw "
        " /directory/2022-12-20_gkg.sif "
        "GkgExecuteCommand Gis2NiftiConverter -i "
        + InputGisFilename
        + " -o "
        + OutputNiftiFilename
        + " -verbose"
    )


def Data_formating(MainDirectory, verbose):
    MaterialsDirectory = os.path.join(MainDirectory, "01-Materials")

    Path_MovingFile = os.path.join(MaterialsDirectory, "Moving")
    Path_FixedFile = os.path.join(MaterialsDirectory, "Fixed")

    FileToConvert = glob.glob(os.path.join(Path_MovingFile, "*.ima")) + glob.glob(
        os.path.join(Path_FixedFile, "*.ima")
    )

    if verbose:
        print("======================================================")
        print("DATA FORMATING, GIS TO NIFTI CONVERSION PART ")
        print("======================================================")

    # GIS CONVERSION
    for gisfile in FileToConvert:
        if not os.path.exists(gisfile.split(".")[0] + ".nii") and not os.path.exists(
            gisfile.split(".")[0] + ".nii.gz"
        ):
            if verbose:
                print("======================================================")
                print("    GIS TO NIFTI CONVERSION OF :" + gisfile.split("/")[-1])
                print("======================================================")
            ConversionGisToNifti(gisfile, gisfile.split(".")[0] + ".nii")
            if verbose:
                print("DONE")

    if verbose:
        print("======================================================")
        print("GIS TO NIFTY CONVERSION PART: DONE")
        print("======================================================")

    # N4 BIAS FIELD CORRECTION
    FileToN4 = (
        glob.glob(os.path.join(Path_MovingFile, "*.nii"))
        + glob.glob(os.path.join(Path_MovingFile, "*.nii.gz"))
        + glob.glob(os.path.join(Path_FixedFile, "*.nii"))
        + glob.glob(os.path.join(Path_FixedFile, "*.nii.gz"))
    )

    # print(FileToN4)
    for Niftifile in FileToN4:
        extension = ".nii.gz" if len(Niftifile.split(".")) > 2 else ".nii"

        if (
            not os.path.exists(Niftifile.split(".")[0] + "_N4.nii")
            and not os.path.exists(Niftifile.split(".")[0] + "_N4.nii.gz")
            and len(Niftifile.split("N4")) < 2
        ):
            if verbose:
                print("======================================================")
                print("    N4 BIAS CORRECTION FIELD :" + Niftifile.split("/")[-1])
                print("======================================================")
            ANTS_volume = ants.image_read(Niftifile)
            ANTS_bias_field = ants.n4_bias_field_correction(
                ANTS_volume,
                spline_param=20,
                rescale_intensities=True,
                return_bias_field=True,
            )
            ants.image_write(
                ANTS_volume / ANTS_bias_field,
                Niftifile.split(".")[0] + "_N4" + extension,
            )
            if verbose:
                print("DONE")
    if verbose:
        print("======================================================")
        print("DATA FORMATING: DONE")
        print("======================================================")


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


def ApplyRegistrationOnRawdata(MainDirectory, verbose):
    MaterialsDirectory = os.path.join(MainDirectory, "01-Materials")
    RegistrationDirectory = os.path.join(MainDirectory, "02-Registration")
    os.makedirs(RegistrationDirectory, exist_ok=True)
    output = os.path.join(MainDirectory, "03-Result")
    os.makedirs(output, exist_ok=True)

    # PROPEDEUTIC TO APPLY REGISTRATION ALGORITHM

    Path_fixedfile = (
        glob.glob(os.path.join(MaterialsDirectory, "Fixed") + "/*N4.nii.gz")
        + glob.glob(os.path.join(MaterialsDirectory, "Fixed") + "/*N4.nii")
    )[0]

    fixed = ants.image_read(Path_fixedfile)

    verbose = bool(verbose)

    if verbose:
        print("======================================================")
        print("APPLY REGISTRATION: ")
        print("======================================================")

    Path_movingfile = os.path.join(MaterialsDirectory, "Moving")
    Path_fixedSessionfile = os.path.join(MaterialsDirectory, "Fixed")
    Path_TransformFiles = os.path.join(RegistrationDirectory, "TransformFiles")

    Filename_fixed = (
        glob.glob(os.path.join(Path_fixedSessionfile, "*N4.nii"))
        + glob.glob(os.path.join(Path_fixedSessionfile, "*N4.nii.gz"))
    )[0]

    ### ATTENTION MODIFICATION POUR APPLIQUER SUR IMAGES BRUTES
    # Filenames_moving =  glob.glob(os.path.join(Path_movingfile, '*N4.nii')) +\
    #     glob.glob(os.path.join(Path_movingfile, '*N4.nii.gz'))
    Filenames_moving = glob.glob(os.path.join(Path_movingfile, "*.nii")) + glob.glob(
        os.path.join(Path_movingfile, "*.nii.gz")
    )
    Filenames_moving_N4 = glob.glob(
        os.path.join(Path_movingfile, "*N4.nii")
    ) + glob.glob(os.path.join(Path_movingfile, "*N4.nii.gz"))

    for name in Filenames_moving_N4:
        Filenames_moving.remove(name)

    Dict_Transforms = {}
    Filename_Fixed = Filename_fixed.split("/")[-1].split(".")[0].split("_")[-2]

    for file in Filenames_moving:
        Filename_Moving = file.split("/")[-1].split(".")[0].split("_")[-1]

        Path_WarpFilename = os.path.join(
            Path_TransformFiles,
            Filename_Moving + "To" + Filename_Fixed + "3Warp.nii.gz",
        )
        Path_AffineFilename = os.path.join(
            Path_TransformFiles, Filename_Moving + "To" + Filename_Fixed + "2Affine.mat"
        )
        Path_RigidFilename = os.path.join(
            Path_TransformFiles, Filename_Moving + "To" + Filename_Fixed + "1Rigid.mat"
        )
        Path_InitialFilename = os.path.join(
            Path_TransformFiles,
            Filename_Moving
            + "To"
            + Filename_Fixed
            + "0DerivedInitialMovingTranslation.mat",
        )

        Path_WarpFilename_2 = os.path.join(
            Path_TransformFiles,
            Filename_Moving + "To" + Filename_Fixed + "1Warp.nii.gz",
        )
        Path_AffineFilename_2 = os.path.join(
            Path_TransformFiles,
            Filename_Moving + "To" + Filename_Fixed + "0GenericAffine.mat",
        )
        # ATTENTION deux varialbles presque pareil Path_TransformFiles et Path_transformfile

        Dict_Transforms[file] = [
            Path_WarpFilename,
            Path_AffineFilename,
            Path_RigidFilename,
            Path_InitialFilename,
        ]

        if not (os.path.exists(Path_WarpFilename)) and not (
            os.path.exists(Path_AffineFilename)
        ):
            Dict_Transforms[file] = [Path_WarpFilename_2, Path_AffineFilename_2]

        if not os.path.exists(
            os.path.join(
                output, file.split("/")[-1].split(".")[0] + "_Registration.nii.gz"
            )
        ):
            if verbose:
                print("======================================================")
                print("    APPLY REGISTRATION, File: " + file.split("/")[-1])
                print("======================================================")
            Moving = ants.image_read(file)
            apply = ants.apply_transforms(
                fixed,
                Moving,
                Dict_Transforms[file],
                "lanczosWindowedSinc",
                imagetype=0,
                verbose=True,
            )

            if verbose:
                print("======================================================")
                print("    APPLY REGISTRATION DONE")
                print("======================================================")
            output_filename = os.path.join(
                output, file.split("/")[-1].split(".")[0] + "_Registration.nii.gz"
            )
            ants.image_write(apply, output_filename)
        else:
            if verbose:
                print("======================================================")
                print(
                    "    APPLY REGISTRATION ALREADY EXIST, File: " + file.split("/")[-1]
                )
                print("======================================================")

    if verbose:
        print("======================================================")
        print("APPLY REGISTRATION: DONE")
        print("======================================================")


def runmain(MainDirectory, verbose):
    Data_formating(MainDirectory, verbose)
    Registration_Type1(MainDirectory, verbose)
    # ApplyRegistrationOnRawdata(MainDirectory, verbose)


parser = optparse.OptionParser()
parser.add_option("-i", "--input", dest="MainDirectory", help="01-/02-/... directory")
parser.add_option(
    "-v",
    "--verbose",
    type="int",
    dest="verbose",
    default=True,
    help="1 --> True, 0 --> False. 1 by Default.",
)

(options, args) = parser.parse_args()


################################################################################
# 1) Value Extractor
################################################################################

runmain(options.MainDirectory, options.verbose)


# os.system("antsRegistration "+\
#             f"-d {3} "+\
#             f"-r [{Path_FixedFilename}, {file}, 1] "+\
#             "-n LanczosWindowedSinc "+\
#                 #
#             f"-m MI[{Path_FixedFilename},{file},{1},{8}, Regular, {0.5}] " +\
#             f"-t Rigid[{2}] "+\
#             f"-c [{500}x{500}x{500},{1e-6},{8}] "+\
#             f"-s {4}x{2}x{1}vox "+\
#             f"-f {4}x{2}x{2} "+\
#                 #
#             f"-m MI[{Path_FixedFilename},{file},{1},{16}, Regular, {0.8}] " +\
#             f"-t Affine[{1}] "+\
#             f"-c [{500}x{500}x{500},{1e-6},{8}] "+\
#             f"-s {4}x{2}x{1}vox "+\
#             f"-f {4}x{2}x{2} "+\
#                 #
#             f"-m mattes[{Path_FixedFilename},{file},{1},{32}, Regular, {0.8}] " +\
#             f"-t SyN[{0.5}] "+\
#             f"-c [{1200}x{1200}x{1000},{1e-6},{8}] "+\
#             f"-s {4}x{2}x{1}vox "+\
#             f"-f {4}x{2}x{2} "+\
#                 #
#             "-u 0 "+\
#             "-z 0 "+\
#             f"-o [{out},{out}_Regis.nii.gz,{out}INV_Regis.nii.gz] "+\
#             " -x [NA,NA] "+\
#             "--float 1 "+\
#             "--write-composite-transform 0 "+\
#             "-v 1")
