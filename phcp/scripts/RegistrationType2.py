#!/usr/bin/env python3
"""
Created on Fri Jan  5 10:06:57 2024

@author: la272118
"""

import json
import logging
import os
import subprocess
import sys

import ants
import nibabel as ni
import numpy as np
from Unwarp import Unwarp

logger = logging.getLogger(__name__)


def run_ant_apply_registration(
    ref_space,
    input_img,
    output_filename,
    transforms,
    interpolation="Linear",
    dimensionality=3,
    invert_flags=None,
) -> None:
    """
    Apply ANTs tranforms to a moving image.

    Parameters:
        ref_space (str) : Reference space.
        input_img (str) : Image to transform.
        output_filename (str) : Output filename.
        transforms (list) : List of transforms filenames.
        interpolation (str) : Interpolation methods.
        dimensionality (int) : dimensionality.
        invert_flags (list) : Optional list of booleans indicating wether to invert each transform.
    """

    command = [
        "antsApplyTransforms",
        "--dimensionality",
        str(dimensionality),
        "--input",
        input_img,
        "--reference-image",
        ref_space,
        "--output",
        output_filename,
        "--interpolation",
        interpolation,
    ]

    if invert_flags is None:
        invert_flags = [False] * len(transforms)

    for tfm, invert in zip(transforms, invert_flags, strict=False):
        if invert:
            command.extend(["--trannform", f"[{tfm},1]"])
        else:
            command.extend(["--trannform", tfm])
    logger.info("========= Run ANTs apply transforms =========")

    try:
        subprocess.run(command, check=True)
        logger.info(" ANTs apply transforms : done ")
    except subprocess.CalledProcessEror as e:
        logger.info(" ERROR : ANTs apply transforms : %s", e)


def Overlap_Mask_Extraction(
    File_directory, fixedfilename, movingfilename, ty1, ty2
) -> None:
    Preparation_directory = os.path.join(File_directory, "02-Preparation")

    meta_fixed = ni.load(fixedfilename)
    meta_moving = ni.load(movingfilename)

    ty = ty2 - ty1

    # arr_fixed = meta_fixed.get_fdata()
    size_fixed = tuple(meta_fixed.header["dim"][1:4])  # arr_fixed.shape

    # arr_moving = meta_moving.get_fdata()
    size_moving = tuple(meta_moving.header["dim"][1:4])  # arr_moving.shape

    resolution_fixed = np.abs(meta_fixed.header["pixdim"][1:4])  # [x,y,z]
    resolution_moving = np.abs(meta_moving.header["pixdim"][1:4])  # [x,y,z]

    # Bande de sécurité de 4 mm
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

    # seg_fixed[:,:,bd_security_fixed:nbrslice_overlap_fixed-bd_security_fixed] = \
    #     np.ones((size_fixed[:2]+(nbrslice_overlap_fixed-2*bd_security_fixed,)))
    # seg_moving[:,:,size_moving[2]+bd_security_moving-nbrslice_overlap_moving:-bd_security_moving] = \
    #     np.ones((size_fixed[:2]+(nbrslice_overlap_moving-2*bd_security_moving,)))
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


def apply_NLMF_and_n4(img_filename, output_filename, verbose) -> None:
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


def Preparation(File_directory, JSON_filename, verbose=0) -> None:
    verbose = bool(verbose)
    Materials_directory = os.path.join(File_directory, "01-Materials")
    Preparation_directory = os.path.join(File_directory, "02-Preparation")
    os.makedirs(Preparation_directory, exist_ok=True)

    with open(JSON_filename) as f:
        description = json.load(f)

    refsuperposition, maskref_filename, refprepasuperposition = "", "", ""
    refty = 0

    for key in description:
        ty = float(description[key][-1])  # Load translation from json file
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

        # $$$$$$$$$$$$$$$$ Application de la superposition ici $$$$$$$$$$$$$$$$$
        if ty != 0.0:
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


def parse_command_line(argv):
    """Parse the script's command line."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Algorithms for inter-FOVs registration. ",
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
    return Preparation(args.workingDirectory, args.jsonFilename, args.verbose) or 0


if __name__ == "__main__":
    sys.exit(main())
