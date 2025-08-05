#!/usr/bin/env python3
import json
import logging
import os
import subprocess
import sys

import ants
import nibabel as ni
import numpy as np

from phcp.registration import run_ants_apply_registration

logger = logging.getLogger(__name__)


def make_header_a_deformation_field_header(arr, affine_matrix) -> ni.Nifti1Image:
    ni_im = ni.Nifti1Image(arr, affine_matrix)
    ni_im.header.set_intent(1007)
    ni_im.header["regular"] = b"r"
    ni_im.header.set_xyzt_units(2)
    sform = ni_im.header.get_sform()
    ni_im.header.set_qform(sform, code=1)
    ni_im.header.set_sform(sform, code=1)
    return ni_im


def Unwarp(
    InputFilename: str | os.PathLike, Warpfilename, headerFlip=True, verbose=True
) -> None:
    logger.info("========= Loading GNLC data =========")
    meta_warp = ni.load(Warpfilename)
    arr_warp = meta_warp.get_fdata()

    basename = os.path.join(
        os.path.dirname(InputFilename), os.path.basename(InputFilename).split(".")[0]
    )
    meta_input = ni.load(InputFilename)

    logger.info("========= Adjust GNLC Deformation Field Header =========")
    size_input = meta_input.header["dim"][1:4]
    size_warp = np.asarray(arr_warp.shape)[0:3]
    resolution_input = np.abs(meta_input.header["pixdim"][1:4])  # [x,y,z]
    resolution_warp = np.abs(meta_warp.header["pixdim"][1:4])  # [x,y,z]

    if headerFlip:
        translation_input = -np.abs(
            np.asarray(
                [
                    float(meta_input.header["qoffset_x"]),
                    float(meta_input.header["qoffset_z"]),
                    float(meta_input.header["qoffset_y"]),
                ]
            )
        )  # HeaderFlip swaps z and y axes
        Input_affine = meta_input.affine
        with np.errstate(divide="ignore", invalid="ignore"):
            Input_affine[0][0:3] = (
                Input_affine[0][0:3] / Input_affine[0][0:3] * resolution_warp[0]
            )
            Input_affine[1][0:3] = (
                -Input_affine[1][0:3] / Input_affine[1][0:3] * resolution_warp[2]
            )
            Input_affine[2][0:3] = (
                Input_affine[2][0:3] / Input_affine[2][0:3] * resolution_warp[1]
            )
            Input_affine[np.isnan(Input_affine)] = 0

    else:
        translation_input = -np.abs(
            np.asarray(
                [
                    float(meta_input.header["qoffset_x"]),
                    float(meta_input.header["qoffset_y"]),
                    float(meta_input.header["qoffset_z"]),
                ]
            )
        )
        Input_affine = meta_input.affine
        with np.errstate(divide="ignore", invalid="ignore"):
            Input_affine[0][0:3] = (
                Input_affine[0][0:3] / Input_affine[0][0:3] * resolution_warp[0]
            )
            Input_affine[1][0:3] = (
                Input_affine[1][0:3] / Input_affine[1][0:3] * resolution_warp[1]
            )
            Input_affine[2][0:3] = (
                Input_affine[2][0:3] / Input_affine[2][0:3] * resolution_warp[2]
            )
            Input_affine[np.isnan(Input_affine)] = 0

    ecart_input = np.int16(-size_input * resolution_input / 2 - translation_input)
    vx_of_interest = np.int16(resolution_input * size_input / resolution_warp)
    lower_bundaries = np.int16(
        (size_warp - vx_of_interest) / 2 - ecart_input / resolution_warp
    )
    higher_bundaries = lower_bundaries + vx_of_interest
    newarr = arr_warp[
        lower_bundaries[0] : higher_bundaries[0],
        lower_bundaries[1] : higher_bundaries[1],
        lower_bundaries[2] : higher_bundaries[2],
        :,
        :,
    ]

    ni_im = make_header_a_deformation_field_header(newarr, Input_affine)

    WarpInputSpaceFilename = basename + "_Warp.nii.gz"

    ni.save(ni_im, WarpInputSpaceFilename)

    logger.info("========= Unwarp Input File =========")

    OutputFilename = basename + "_rec-unwarp.nii.gz"

    run_ants_apply_registration(
        ref_space=InputFilename,
        input_img=InputFilename,
        output_filename=OutputFilename,
        transforms=[WarpInputSpaceFilename],
        interpolation="Linear",
        dimensionality=3,
        invert_flags=None,
    )


def Overlap_Mask_Extraction(
    File_directory: str | os.PathLike,
    fixedfilename: str,
    movingfilename: str,
    ty1: float,
    ty2: float,
) -> None:
    Preparation_directory = os.path.join(File_directory, "02-Preparation")

    meta_fixed = ni.load(fixedfilename)
    meta_moving = ni.load(movingfilename)

    ty = ty2 - ty1

    size_fixed = tuple(meta_fixed.header["dim"][1:4])  # arr_fixed.shape

    size_moving = tuple(meta_moving.header["dim"][1:4])  # arr_moving.shape

    resolution_fixed = np.abs(meta_fixed.header["pixdim"][1:4])  # [x,y,z]
    resolution_moving = np.abs(meta_moving.header["pixdim"][1:4])  # [x,y,z]

    # Safety band of 4 mm
    largeur = 4
    # number of slice in the fixed img
    bd_security_fixed = np.int16(largeur / resolution_fixed[2])
    # number of slice in the moving img
    bd_security_moving = np.int16(largeur / resolution_moving[2])

    seg_fixed, seg_moving = np.zeros(size_fixed), np.zeros(size_moving)

    # calculus of the overlap in mm
    overlap = np.int16(
        0.5
        * (size_fixed[2] * resolution_fixed[2] + size_moving[2] * resolution_moving[2])
        - ty
    )

    # now in vx
    nbrslice_overlap_moving = np.int16(overlap / resolution_moving[2])
    nbrslice_overlap_fixed = np.int16(overlap / resolution_fixed[2])

    json_file = dict()
    json_file["nbrslice_overlap_fixed"] = str(nbrslice_overlap_fixed)
    json_file["nbrslice_overlap_moving"] = str(nbrslice_overlap_moving)
    json_file["bd_security_fixed"] = str(bd_security_fixed)
    json_file["bd_security_moving"] = str(bd_security_moving)
    json_file["ty"] = ty

    json_filename = os.path.join(
        Preparation_directory,
        fixedfilename.split("/")[-1].split(".")[0]
        + "_"
        + movingfilename.split("/")[-1].split(".")[0]
        + ".json",
    )

    with open(json_filename, "w") as fp:
        json.dump(json_file, fp)

    seg_fixed[:, :, bd_security_fixed:nbrslice_overlap_fixed] = np.ones(
        size_fixed[:2] + (nbrslice_overlap_fixed - bd_security_fixed,)
    )
    seg_moving[
        :, :, size_moving[2] + bd_security_moving - nbrslice_overlap_moving :
    ] = np.ones(size_fixed[:2] + (nbrslice_overlap_moving - bd_security_moving,))

    seg_fixed_filename = os.path.join(
        Preparation_directory,
        fixedfilename.split("/")[-1].split(".")[0] + "_maskref.nii.gz",
    )
    seg_moving_filename = os.path.join(
        Preparation_directory,
        movingfilename.split("/")[-1].split(".")[0] + "_mask.nii.gz",
    )

    meta_fixed.affine[1][3] += ty1
    ni_im = ni.Nifti1Image(seg_fixed, meta_fixed.affine)
    ni.save(ni_im, seg_fixed_filename)

    meta_moving.affine[1][3] += ty2
    ni_im = ni.Nifti1Image(seg_moving, meta_moving.affine)
    ni.save(ni_im, seg_moving_filename)


def RASyN_overlap(
    File_directory: str | os.PathLike,
    fixedfilename: str,
    movingfilename: str,
    seg_fixed_filename: str,
    seg_moving_filename: str,
    dimensionality=3,
    SamplingStrategy="Random",
    Radius="NA",
    number_of_bins=16,
    interpolation="Linear",
) -> None:
    TransformFilesDirectory = os.path.join(File_directory, "03-TransformFiles")
    os.makedirs(TransformFilesDirectory, exist_ok=True)

    Prefix = os.path.join(
        TransformFilesDirectory,
        movingfilename.split("/")[-1].split("_")[1]
        + "To"
        + fixedfilename.split("/")[-1].split("_")[1],
    )
    # f"MeanSquares[{fixedfilename},{movingfilename},{1},{Radius},{SamplingStrategy},{0.5}]"
    command = [
        "antsRegistration",
        "--dimensionality",
        str(dimensionality),
        "--output",
        Prefix,
        "--interpolation",
        interpolation,
        "--metric",
        f"MI[{fixedfilename},{movingfilename},{1},{number_of_bins},{SamplingStrategy},{0.5}]",
        "--transform",
        f"Rigid[{1}]",
        "--convergence",
        f"[{500}x{200}x{100},{1e-6},{8}]",
        "--smoothing-sigmas",
        f"{4}x{2}x{1}vox",
        "--shrink-factors",
        f"{4}x{2}x{2}",
        "--masks",
        f"[{seg_fixed_filename},{seg_moving_filename}]",
        "--metric",
        f"MI[{fixedfilename},{movingfilename},{1},{number_of_bins},{SamplingStrategy},{0.5}]",
        "--transform",
        f"Affine[{1}]",
        "--convergence",
        f"[{500}x{200}x{100},{1e-6},{8}]",
        "--smoothing-sigmas",
        f"{4}x{2}x{1}vox",
        "--shrink-factors",
        f"{4}x{2}x{2}",
        "--masks",
        f"[{seg_fixed_filename},{seg_moving_filename}]",
        "--metric",
        f"MI[{fixedfilename},{movingfilename},{1},{2 * number_of_bins},{SamplingStrategy},{0.8}]",
        "--transform",
        f"SyN[{1}]",
        "--convergence",
        f"[{500}x{200}x{100},{1e-6},{8}]",
        "--smoothing-sigmas",
        f"{4}x{2}x{1}vox",
        "--shrink-factors",
        f"{4}x{2}x{2}",
        "--masks",
        f"[{seg_fixed_filename},{seg_moving_filename}]",
        "--verbose",
    ]
    logger.info("========= Run ANTs Registration =========")

    try:
        subprocess.run(command, check=True)
        logger.info(" ANTs Registration : done ")
    except subprocess.CalledProcessError as e:
        logger.info(" ERROR : ANTs Registration : %s", e)


def apply_NLMF_and_n4(
    img_filename: str, output_filename: str, verbose: int | bool
) -> None:
    ANTSimage = ants.image_read(img_filename)
    n4_antsimg = ants.n4_bias_field_correction(
        ANTSimage, spline_param=20, return_bias_field=True
    )
    ANTSimage = ANTSimage / n4_antsimg
    ANTSimage.denoise_image(r=8, noise_model="Rician", v=int(verbose))
    # Save correction
    ants.image_write(ANTSimage, output_filename)


def correct_header(
    DeformationFieldFilename: str, UnwarpedFilename: str, ty: float
) -> None:
    matrix = np.eye(4, 4)
    matrix[1][-1] = ty

    meta = ni.load(UnwarpedFilename)
    meta_warp = ni.load(DeformationFieldFilename)

    meta.affine[1][3] += ty
    meta_warp.affine[1][3] += ty

    ni.save(meta, UnwarpedFilename)
    ni.save(meta_warp, DeformationFieldFilename)


def init_intermediate_space(
    number_of_fovs: int,
    img_ref_filename: str | os.PathLike,
    intermediate_space_filename: str | os.PathLike,
) -> None:
    meta = ni.load(img_ref_filename)
    newarr = np.ones((460, 480, number_of_fovs * 450))
    resolution_img = meta.header["pixdim"]
    matrix = np.asarray(
        [
            [resolution_img[1], 0, 0, 0],
            [0, resolution_img[2], 0, 0],
            [0, 0, resolution_img[3], 0],
            [0, 0, 0, 1],
        ]
    )
    ni_im = ni.Nifti1Image(newarr, matrix, meta.header)
    ni.save(ni_im, intermediate_space_filename)


def load_translation_from_json(json_object):
    """Decode the translation from the description JSON file.

    Old encoding: ["10.0"]
    New encoding: 10 or 10.0
    """
    if isinstance(json_object, list):
        json_object = json_object[-1]
    return float(json_object)


def Preparation(
    File_directory: str | os.PathLike, JSON_filename: str, verbose=0
) -> None:
    verbose = bool(verbose)
    Materials_directory = os.path.join(File_directory, "01-Materials")
    Preparation_directory = os.path.join(File_directory, "02-Preparation")
    os.makedirs(Preparation_directory, exist_ok=True)
    TransformFilesDirectory = os.path.join(File_directory, "03-TransformFiles")
    os.makedirs(TransformFilesDirectory, exist_ok=True)

    with open(JSON_filename) as f:
        description = json.load(f)

    refsuperposition, maskref_filename, refprepasuperposition = "", "", ""
    refty = 0

    for key in description:
        ty = load_translation_from_json(description[key])
        key_directory = os.path.join(Materials_directory, key + ".nii.gz")
        keycorrected_directory = os.path.join(
            Preparation_directory, key + "_NLMF_N4.nii.gz"
        )

        mask_filename = os.path.join(Preparation_directory, key + "_mask.nii.gz")

        if not (os.path.exists(keycorrected_directory)):
            logger.info("NLMF & N4 correction: %s", key)
            apply_NLMF_and_n4(key_directory, keycorrected_directory, verbose)
            logger.info("NLMF & N4 correction: done.")

        if ty != 0.0:  # ty=0 correspond to the FOV of reference (posterior part)
            seg_fixed_filename = os.path.join(
                Preparation_directory,
                refprepasuperposition.split("/")[-1].split(".")[0] + "_maskref.nii.gz",
            )
            if not (os.path.exists(seg_fixed_filename)):
                logger.info("Creating masks for overlapping zones")
                Overlap_Mask_Extraction(
                    File_directory, refprepasuperposition, key_directory, refty, ty
                )
                logger.info("Masks : done.")

        Warpfilename = os.path.join(Materials_directory, "WarpHeaderFile_FWHM1.nii.gz")
        unwarpfilename = os.path.join(
            Preparation_directory, key + "_NLMF_N4_rec-unwarp.nii.gz"
        )

        if not (os.path.exists(unwarpfilename)):
            logger.info("========= Unwarp.py part =========")
            Unwarp(keycorrected_directory, Warpfilename)
            logger.info("Unwarp.py part : done.")

            if ty != 0.0:
                logger.info("Header correction: %s", unwarpfilename.split("/")[-1])
                WarpFile_directory = os.path.join(
                    Preparation_directory, key + "_NLMF_N4_Warp.nii.gz"
                )
                correct_header(WarpFile_directory, unwarpfilename, ty)
                logger.info("Header correction: done.")

        if ty != 0.0:
            Prefix = os.path.join(
                TransformFilesDirectory,
                unwarpfilename.split("/")[-1].split("_")[1]
                + "To"
                + refsuperposition.split("/")[-1].split("_")[1],
            )
            if not (os.path.exists(Prefix + "0GenericAffine.mat")):
                RASyN_overlap(
                    File_directory,
                    refsuperposition,
                    unwarpfilename,
                    maskref_filename,
                    mask_filename,
                )

        maskref_filename = os.path.join(Preparation_directory, key + "_maskref.nii.gz")
        refprepasuperposition = key_directory
        refsuperposition = unwarpfilename
        refty = ty

    intermediate_space = os.path.join(File_directory, "IntermediateSpace.nii.gz")
    if not (os.path.exists(intermediate_space)):
        logger.info("========= Init RefSpace =========")
        img_ref_name = next(
            (k for k, v in description.items() if load_translation_from_json(v) == 0.0),
            None,
        )
        pos_filename = os.path.join(
            Preparation_directory, img_ref_name + "_NLMF_N4_rec-unwarp.nii.gz"
        )

        init_intermediate_space(len(description), pos_filename, intermediate_space)
        logger.info("Init RefSpace : done.")
    logger.info(
        "Preparation step : done. WARNING : Before selecting the Fusion step, "
        " be sure to create the transformation matrix (.txt or .mat) "
        " so that the reference FOV is correctly positioned in the intermediate space"
    )


def parse_command_line(argv):
    """Parse the script's command line."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Prepare inter-FOV registration. "
        "Be sure to create the transformation matrix after running this preparation script "
        "so that the reference FOV is correctly positioned in the intermediate space",
    )
    parser.add_argument(
        "-i",
        "--workingDirectory",
        required=True,
        help="Working directory such as DirectoryName/01-Materials   /02-..   /03-..",
    )
    parser.add_argument(
        "-j",
        "--jsonFilename",
        required=True,
        help="JSON filename typically sub-${sub}_ses-{ses} : ['float']",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="print detailed status information",
    )

    args = parser.parse_args()
    return args


def main(argv=sys.argv):
    """The script's entry point."""
    logging.basicConfig(level=logging.INFO)
    args = parse_command_line(argv)
    return Preparation(args.workingDirectory, args.jsonFilename, args.verbose) or 0


if __name__ == "__main__":
    sys.exit(main())
