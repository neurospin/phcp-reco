#!/usr/bin/env python3
"""
Created on Fri Mar 22 10:26:51 2024

@author: la272118
"""

import logging
import os
import subprocess

import nibabel as ni
import numpy as np

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
            command.extend(["--transform", f"[{tfm},1]"])
        else:
            command.extend(["--transform", tfm])
    logger.info("========= Run ANTs apply transforms =========")

    try:
        subprocess.run(command, check=True)
        logger.info(" ANTs apply transforms : done ")
    except subprocess.CalledProcessError as e:
        logger.info(" ERROR : ANTs apply transforms : %s", e)


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

    basename = InputFilename.split(".")[0]
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

    run_ant_apply_registration(
        ref_space=InputFilename,
        input_img=InputFilename,
        output_filename=OutputFilename,
        transforms=[WarpInputSpaceFilename],
        interpolation="Linear",
        dimensionality=3,
        invert_flags=None,
    )
