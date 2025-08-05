#!/usr/bin/env python3
import json
import logging
import os
import sys

import nibabel as ni
import numpy as np
import scipy.ndimage as nd
from sklearn.mixture import GaussianMixture

from phcp.registration import run_ants_apply_registration

logger = logging.getLogger(__name__)


def insert_transformation_file(
    transformationfile_directory: str | os.PathLike,
    number_iter: int,
    act_session: str,
    prec_session: str,
    transformations_list: list[str],
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


def merge_FOVs(File_directory: str | os.PathLike, JSON_filename: str) -> None:
    preparation_directory = os.path.join(File_directory, "02-Preparation")
    transformfile_directory = os.path.join(File_directory, "03-TransformFiles")
    distribution_directory = os.path.join(File_directory, "04-FovsDistribution")
    os.makedirs(distribution_directory, exist_ok=True)
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

        logger.info(f"Send {filename} to the intermediate space : Done")

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
        logger.info("Creation of FOV mapping : Done")
        logger.info("========= Creation of FOV weight =========")
        calcul_weight_FOVs(output_FOV_filename)
        logger.info("Creation of FOV weight : Done")

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
            logger.info("Concat consecutive FOVs : Done")

        precedent_session = actual_session
    logger.info("========= Save the reconstructed block =========")
    ni.save(
        reconstructed_block_nifti_im,
        os.path.join(File_directory, "Reconstructed_block.nii.gz"),
    )


def merge_consecutive_FOVs(
    reconstruction_block_ni_im: ni.Nifti1Image,
    actual_fov_in_IS_filename: str | os.PathLike,
    prec_session: str,
    sub: str,
    dist_directory: str | os.PathLike,
) -> ni.Nifti1Image:
    meta_act_fov_is = ni.load(actual_fov_in_IS_filename)
    arr_act_fov_is = meta_act_fov_is.get_fdata()

    reconstructed_block_arr = reconstruction_block_ni_im.get_fdata()

    # init materials
    precedent_subject_name = "".join(
        [os.path.join(dist_directory, sub), "_", prec_session, "_T2w"]
    )
    prec_fov_mask_filename = "".join([precedent_subject_name, "_FOV_IS.nii.gz"])

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

    coefficient = find_coefficient_to_correct_intensity(
        arr_act_fov_is, reconstructed_block_arr, overlap
    )

    arr_act_fov_is = arr_act_fov_is * coefficient

    values_in_overlap = (
        reconstructed_block_arr * prec_weight_in_overlap
        + arr_act_fov_is * act_weight_in_overlap
    ) * overlap
    values_act_outside_overlap = arr_act_fov_is * (meta_mask_act.get_fdata() - overlap)
    return ni.Nifti1Image(
        reconstructed_block_arr * np.logical_not(overlap)
        + values_in_overlap
        + values_act_outside_overlap,
        reconstruction_block_ni_im.affine,
        reconstruction_block_ni_im.header,
    )


def sigmoid(x: np.ndarray, k: float, x0: float) -> np.ndarray:
    return 1 / (1 + np.exp(-k * (x - x0)))


def find_coefficient_to_correct_intensity(
    actual_arr: np.ndarray, precedent_arr: np.ndarray, overlap: np.ndarray
) -> float:
    actual_arr_filtered = actual_arr[np.bool_(overlap)]
    precedent_arr_filtered = precedent_arr[np.bool_(overlap)]
    actual_arr_filtered = actual_arr_filtered[actual_arr_filtered > 0]
    precedent_arr_filtered = precedent_arr_filtered[precedent_arr_filtered > 0]
    gmm = GaussianMixture(2)
    gmm2 = GaussianMixture(2)
    gmm.fit(actual_arr_filtered.reshape((-1, 1)))
    gmm2.fit(precedent_arr_filtered.reshape((-1, 1)))
    return gmm2.means_.flatten().max() / gmm.means_.flatten().max()


def calcul_weight_FOVs(fov_mapping_filename: str | os.PathLike) -> None:
    meta = ni.load(fov_mapping_filename)
    arr = meta.get_fdata()
    distance_mapping = nd.distance_transform_edt(arr, meta.header["pixdim"][1])
    weight_mapping = sigmoid(distance_mapping, 2, 4)
    ni.save(
        ni.Nifti1Image(weight_mapping, meta.affine, meta.header),
        fov_mapping_filename.replace(".nii.gz", "_weight.nii.gz"),
    )


def parse_command_line(argv):
    """Parse the script's command line."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Perform inter-FOV registration. "
        "Be sure to run phcp-interfov-registration and create the transformation matrix "
        "before running this script.",
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
    return merge_FOVs(args.workingDirectory, args.jsonFilename) or 0


if __name__ == "__main__":
    sys.exit(main())
