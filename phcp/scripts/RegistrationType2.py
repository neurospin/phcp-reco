#!usr/bin/env python3
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
import scipy.ndimage as nd
from Unwarp import Unwarp

logger = logging.getLogger(__name__)


def run_ants_apply_registration(
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


def Preparation(File_directory, JSON_filename, verbose=0) -> None:
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
        img_ref_name = next((k for k, v in description.items() if v == ["0"]), None)
        pos_filename = os.path.join(
            Preparation_directory, img_ref_name + "_NLMF_N4_rec-unwarp.nii.gz"
        )

        init_intermediate_space(len(description), pos_filename, intermediate_space)
        logger.info("Init RefSpace : done.")


def insert_transformation_file(
    transformationfile_directory,
    number_iter,
    act_session,
    prec_session,
    transformations_list,
) -> None:
    if number_iter > 0:
        prefix = "".join([act_session, "To", prec_session])
        warp_transformation_filename = os.path.join(
            transformationfile_directory, "".join([prefix, "1Warp.nii.gz"])
        )
        mat_transformation_filename = os.path.join(
            transformationfile_directory, "".join([prefix, "0GenericAffine.mat"])
        )
        if number_iter == 1:
            transformations_list.append(warp_transformation_filename)
            transformations_list.append(mat_transformation_filename)
        else:
            transformations_list.insert(
                -2 * (number_iter - 1), warp_transformation_filename
            )
            transformations_list.insert(
                -2 * (number_iter - 1), mat_transformation_filename
            )


def merge_FOVs(File_directory, JSON_filename) -> None:
    preparation_directory = os.path.join(File_directory, "02-Preparation")
    transformfile_directory = os.path.join(File_directory, "03-TransformFiles")
    distribution_directory = os.path.join(File_directory, "04-FovsDistribution")

    intermediate_space_filename = os.path.join(
        File_directory, "IntermediateSpace.nii.gz"
    )
    initspace_to_refspace_transform_filename = os.path.join(
        File_directory, "refFOV_to_intermediateSpace.txt"
    )

    with open(JSON_filename) as f:
        description = json.load(f)

    transforms_list = [initspace_to_refspace_transform_filename]
    precedent_session = ""

    for i, subject in enumerate(description):
        filename = subject + "_NLMF_N4_rec-unwarp.nii.gz"
        logger.info(f"========= Send {filename} to the intermediate space =========")
        input_filename = os.path.join(preparation_directory, filename)
        output_filename = os.path.join(
            distribution_directory, filename.replace(".nii.gz", "_IS.nii.gz")
        )
        actual_session = subject.split("_")[1]

        insert_transformation_file(
            transformfile_directory,
            i,
            actual_session,
            precedent_session,
            transforms_list,
        )

        run_ants_apply_registration(
            intermediate_space_filename,
            input_filename,
            output_filename,
            transforms_list,
        )  # Apply transforms_list first tranformation = last applied transformation

        # logger.info("Done") no need because run_ants_apply_registration already return something when ending

        logger.info("========= Creation of FOV mapping =========")
        meta = ni.load(input_filename)
        FOV_filename = input_filename.replace(
            "_NLMF_N4_rec-unwarp.nii.gz", "_FOV.nii.gz"
        )
        output_FOV_filename = os.path.join(
            distribution_directory,
            FOV_filename.split("/")[-1].replace(".nii.gz", "_IS.nii.gz"),
        )
        ni.save(
            ni.Nifti1Image(np.ones(meta.shape), meta.affine, meta.header), FOV_filename
        )
        run_ants_apply_registration(
            intermediate_space_filename,
            FOV_filename,
            output_FOV_filename,
            transforms_list,
            interpolation="NearestNeighbor",
        )
        # logger.info("Done")
        logger.info("========= Creation of FOV weight =========")
        calcul_weight_FOVs(output_FOV_filename)
        logger.info("Done")

        if not (i > 0):
            reconstructed_block_nifti_im = ni.load(output_filename)
        else:
            logger.info("========= Concat consecutive FOVs =========")
            reconstructed_block_nifti_im = merge_consecutive_FOVs(
                reconstructed_block_nifti_im,
                output_filename,
                precedent_session,
                sub=subject.split("_")[0],
                dist_directory=distribution_directory,
            )
            logger.info("Done")

        precedent_session = actual_session
    logger.info("========= Save the reconstructed block =========")
    ni.save(
        reconstructed_block_nifti_im,
        os.path.join(File_directory, "Reconstructed_block.nii.gz"),
    )


def merge_consecutive_FOVs(
    reconstruction_block_ni_im,
    actual_fov_in_IS_filename,
    prec_session,
    sub,
    dist_directory,
):
    meta_act_fov_is = ni.load(actual_fov_in_IS_filename)
    arr_act_fov_is = meta_act_fov_is.get_fdata()

    reconstructed_block_arr = reconstruction_block_ni_im.get_fdata()

    # init materials
    precedent_subject_name = "".join(
        [os.path.join(dist_directory, sub), "_", prec_session, "_T2w"]
    )
    prec_fov_mask_filename = "".join([precedent_subject_name, "_FOV_IS.nii.gz"])
    # prec_fov_weight_filename = "".join([precedent_subject_name, "_FOV_IS_weight.nii.gz"])

    actual_subject_filename = os.path.join(
        dist_directory,
        actual_fov_in_IS_filename.split("/")[-1].replace(
            "_NLMF_N4_rec-unwarp_IS.nii.gz", ""
        ),
    )
    act_fov_mask_filename = "".join([actual_subject_filename, "_FOV_IS.nii.gz"])
    act_fov_weight_filename = "".join(
        [actual_subject_filename, "_FOV_IS_weight.nii.gz"]
    )

    meta_mask_prec = ni.load(prec_fov_mask_filename)
    meta_mask_act = ni.load(act_fov_mask_filename)
    overlap = meta_mask_prec.get_fdata() * meta_mask_act.get_fdata()

    meta_act_weight = ni.load(act_fov_weight_filename)

    act_weight_in_overlap = overlap * meta_act_weight.get_fdata()
    prec_weight_in_overlap = 1 - act_weight_in_overlap

    values_in_overlap = (
        reconstructed_block_arr * prec_weight_in_overlap
        + arr_act_fov_is * act_weight_in_overlap
    )
    values_act_outside_overlap = arr_act_fov_is * (meta_mask_act.get_fdata() - overlap)
    return ni.Nifti1Image(
        reconstructed_block_arr * np.logical_not(overlap)
        + values_in_overlap
        + values_act_outside_overlap,
        reconstruction_block_ni_im.affine,
        reconstruction_block_ni_im.header,
    )


def sigmoid(x, k, x0):
    return 1 / (1 + np.exp(-k * (x - x0)))


def calcul_weight_FOVs(fov_mapping_filename) -> None:
    meta = ni.load(fov_mapping_filename)
    arr = meta.get_fdata()
    distance_mapping = nd.distance_transform_edt(arr, meta.header["pixdim"][1])
    weight_mapping = sigmoid(distance_mapping, 1.5, 9)
    ni.save(
        ni.Nifti1Image(weight_mapping, meta.affine, meta.header),
        fov_mapping_filename.replace(".nii.gz", "_weight.nii.gz"),
    )


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
        "-r",
        "--runSendFOVs",
        action="store_true",
        help="Run the second part of the Inter-FOVs registration pipeline if this flag is present",
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
    if args.runSendFOVs:
        return merge_FOVs(args.workingDirectory, args.jsonFilename) or 0
    else:
        return Preparation(args.workingDirectory, args.jsonFilename, args.verbose) or 0


if __name__ == "__main__":
    sys.exit(main())
