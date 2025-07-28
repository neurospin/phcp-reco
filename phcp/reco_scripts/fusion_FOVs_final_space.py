import json
import logging
import subprocess
import sys
from pathlib import Path

import nibabel as ni
import nibabel.processing as proc
import numpy as np
import scipy.ndimage as nd

from phcp.gkg import gkg_convert_gis_to_nifti

logger = logging.getLogger(__name__)


QMAP_VALUES = [
    "QT1",
    "QT2",
    "QT2star",
    "ADC1500",
    "ADC4500",
    "ADC8000",
    "FA1500",
    "FA4500",
    "FA8000",
    "TransDiff1500",
    "TransDiff4500",
    "TransDiff8000",
    "ParaDiff1500",
    "ParaDiff4500",
    "ParaDiff8000",
    "GFA1500",
    "GFA4500",
    "GFA8000",
    "T2w",
]
CM_QMAP = ["CM_T1", "CM_T2", "CM_T2star"]
DMAP_VALUES = [
    "ADC1500",
    "ADC4500",
    "ADC8000",
    "FA1500",
    "FA4500",
    "FA8000",
    "TransDiff1500",
    "TransDiff4500",
    "TransDiff8000",
    "ParaDiff1500",
    "ParaDiff4500",
    "ParaDiff8000",
    "GFA1500",
    "GFA4500",
    "GFA8000",
]

""" Miscellaneous algorithms """


def load_jsonfile(json_filename: str | Path) -> dict:
    with open(str(json_filename)) as f:
        dictionary = json.load(f)
    return dictionary


def sigmoid(x: np.ndarray, k: float, x0: float) -> np.ndarray:
    return 1 / (1 + np.exp(-k * (x - x0)))


def run_ants_apply_registration(
    ref_space: str,
    input_img: str,
    output_filename: str,
    transforms: list[str],
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
        transforms (list) : List of transforms filenames. (First apply last)
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


""" Preparation step algorithms """


def send_to_RefSpace(
    input_filename: str | Path,
    output_filename: str | Path,
    json_filename: str | Path,
    name_fov: str,
    interpolator="Linear",
) -> None:
    dictionnary = load_jsonfile(json_filename)
    refspace_filename = dictionnary["RefSpace"]
    transforms_pipeline = dictionnary[name_fov]["tparams"]
    run_ants_apply_registration(
        ref_space=refspace_filename,
        input_img=input_filename,
        output_filename=output_filename,
        transforms=transforms_pipeline,
        interpolation=interpolator,
    )


def create_mask_FOV(inputfilename: str | Path, outputfilename: str | Path) -> None:
    meta_input = ni.load(str(inputfilename))
    res = np.ones(meta_input.shape)
    ni_im = ni.Nifti1Image(res, meta_input.affine, meta_input.header)
    ni.save(ni_im, outputfilename)


def create_geometric_penalty_from_fov_mapping(
    mask_fov_in_respace_filename: str | Path, geo_penalty_filename: str | Path
) -> None:
    meta = ni.load(mask_fov_in_respace_filename)
    arr = meta.get_fdata()
    newarr = nd.distance_transform_edt(arr, meta.header["pixdim"][1])
    newarr = sigmoid(newarr, 1.5, 9)
    ni.save(ni.Nifti1Image(newarr, meta.affine, meta.header), geo_penalty_filename)


def run_FOV_mapping_and_geometric_penaty_creation(
    liste_fovs: list[str],
    init_space_path: str | Path,
    ref_space_path: str | Path,
    working_path: str | Path,
) -> None:
    logger.info("========= Creating binary FOV masks =========")
    # This algortihsm is based on the t2star map
    for fov in liste_fovs:
        logger.info(f"========= Process FOV : {fov} =========")
        output_mask_nifti_filename_initspace = init_space_path / "".join(
            ["Mask_", fov, "_FOV.nii.gz"]
        )
        t2star_filename = init_space_path / "".join(["T2star_", fov, ".nii.gz"])

        output_mask_nifti_filename_refspace = ref_space_path / "".join(
            ["Mask_", fov, "_FOV.nii.gz"]
        )
        geometric_penalty_filename = ref_space_path / "".join(
            [fov, "_geometric_penalty.nii.gz"]
        )

        if not (output_mask_nifti_filename_initspace.exists()):
            create_mask_FOV(t2star_filename, output_mask_nifti_filename_initspace)

        logger.info(f"========= Send to refspace process : {fov} =========")
        if not (output_mask_nifti_filename_refspace.exists()):
            sendtorefspace_json_filename = (
                Path(working_path) / "SendToRefSpace_T2star.json"
            )
            send_to_RefSpace(
                input_filename=output_mask_nifti_filename_initspace,
                output_filename=output_mask_nifti_filename_refspace,
                json_filename=sendtorefspace_json_filename,
                name_fov=fov,
                interpolator="NearestNeighbor",
            )

        logger.info(f"========= Geometric penalty creation process : {fov} =========")
        if not (geometric_penalty_filename.exists()):
            create_geometric_penalty_from_fov_mapping(
                output_mask_nifti_filename_refspace, geometric_penalty_filename
            )
        logger.info("Done.")


def convert_gis_files_in_nifti_files(
    init_space_path: str | Path,
) -> None:
    logger.info("========= Convert GIS files into Nifti =========")

    ima_generator = init_space_path.glob("*.ima")
    for gis_filename in ima_generator:
        nifti_filename = gis_filename.with_name(gis_filename.stem + ".nii.gz")
        if not (nifti_filename.exists()):
            gkg_convert_gis_to_nifti(str(gis_filename), str(nifti_filename))


def send_nifti_files_in_init_to_refspace(
    init_space_path: str | Path,
    ref_space_path: str | Path,
    working_path: str | Path,
) -> None:
    logger.info("========= Send nifti file in initspace to refspace =========")

    nifti_generator = (
        filename
        for filename in init_space_path.glob("*.nii.gz")
        if not filename.name.startswith("Mask")
    )
    for nifti_filename_in_initspace in nifti_generator:
        modality, fov = (nifti_filename_in_initspace.name.split(".")[0].split("_"))[:2]

        if modality in DMAP_VALUES:
            modality = "DWI"

        JsonFilename = "".join(["SendToRefSpace_", modality, ".json"])

        output_refspace_filename = ref_space_path / nifti_filename_in_initspace.name
        if not (output_refspace_filename.exists()):
            logger.info(
                f"========= Send to refspace process :{nifti_filename_in_initspace.name} ========="
            )
            send_to_RefSpace(
                input_filename=nifti_filename_in_initspace,
                output_filename=output_refspace_filename,
                json_filename=Path(working_path) / JsonFilename,
                name_fov=fov,
            )


""" Fusion modality algorithms """


def create_weight(cm_filename1: str, cm_filename2: str) -> np.ndarray:
    CM_1 = ni.load(cm_filename1)
    CM_2 = ni.load(cm_filename2)

    CM_arr1 = proc.smooth_image(CM_1, proc.sigma2fwhm(0.2)).get_fdata()
    CM_arr2 = proc.smooth_image(CM_2, proc.sigma2fwhm(0.2)).get_fdata()

    return sigmoid(np.log(CM_arr1 / CM_arr2), 4, 0)


def load_geometric_penalty(filename_list: list[str], indice: int) -> np.ndarray:
    fov = Path(filename_list[indice]).name.split("_")[1]
    parent = Path(filename_list[indice]).parent
    meta = ni.load(parent / f"{fov}_geometric_penalty.nii.gz")
    return meta.get_fdata()


def load_mask_from_geometric_penalty(
    filename_list: list[str], indice: int
) -> np.ndarray:
    fov = Path(filename_list[indice]).name.split("_")[1]
    parent = Path(filename_list[indice]).parent
    meta = ni.load(parent / f"{fov}_geometric_penalty.nii.gz")
    arr = meta.get_fdata()
    arr = np.where(arr < 0.001, 0, 1)
    return arr


def fill_Qvalues_Overlapping_relaxometry(
    cm_filenames_list: list[str],
    qmap_filenames_list: list[str],
    mask_filenames_list: list[str],
    overlap_mask: np.ndarray,
    i: int,
) -> np.ndarray:
    Qvalues_Overlapping = np.zeros(ni.load(cm_filenames_list[0]).shape)
    arr_mask_weight = load_geometric_penalty(mask_filenames_list, indice=i + 1)
    arr_mask_weight_i = load_geometric_penalty(mask_filenames_list, indice=i)

    weight_2 = (
        create_weight(cm_filenames_list[i], cm_filenames_list[i + 1]) * arr_mask_weight
    )
    weight_1 = (1 - weight_2) * arr_mask_weight_i
    weight_2 = 1 - weight_1
    Qvalues_Overlapping += (
        weight_1 * ni.load(qmap_filenames_list[i]).get_fdata()
        + weight_2 * ni.load(qmap_filenames_list[i + 1]).get_fdata()
    ) * overlap_mask
    return Qvalues_Overlapping


def fill_Qvalues_Overlapping_OnlySpatialRegularization(
    qmap_filenames_list: list[str],
    mask_filenames_list: list[str],
    overlap_mask: np.ndarray,
    i: int,
) -> np.ndarray:
    Qvalues_Overlapping = np.zeros(ni.load(qmap_filenames_list[0]).shape)
    arr_mask_weight = load_geometric_penalty(mask_filenames_list, indice=i + 1)

    arr_mask_weight_i = load_geometric_penalty(mask_filenames_list, indice=i)

    weight_2 = arr_mask_weight
    weight_1 = (1 - weight_2) * arr_mask_weight_i
    weight_2 = 1 - weight_1

    arr_mask = load_mask_from_geometric_penalty(mask_filenames_list, indice=i + 1)
    arr_mask_i = load_mask_from_geometric_penalty(mask_filenames_list, indice=i)

    weight_2 = weight_2 * arr_mask
    weight_1 = weight_1 * arr_mask_i

    Qvalues_Overlapping += (
        weight_1 * ni.load(qmap_filenames_list[i]).get_fdata()
        + weight_2 * ni.load(qmap_filenames_list[i + 1]).get_fdata()
    ) * overlap_mask
    return Qvalues_Overlapping


def reconstruction(
    cm_filenames_list: list[str],
    qmap_filenames_list: list[str],
    mask_filenames_list: list[str],
    modality_in_DMAP_VALUES: bool,
) -> np.ndarray:
    res = (
        load_mask_from_geometric_penalty(mask_filenames_list, indice=0)
        * ni.load(qmap_filenames_list[0]).get_fdata()
    )
    for i in range(len(qmap_filenames_list) - 1):
        arr_mask = load_mask_from_geometric_penalty(mask_filenames_list, indice=i + 1)
        arr_mask_i = load_mask_from_geometric_penalty(mask_filenames_list, indice=i)

        if i == 0:
            overlap_mask = arr_mask * arr_mask_i
            res = (
                (arr_mask * ni.load(qmap_filenames_list[i + 1]).get_fdata()) + res
            ) * np.logical_not(overlap_mask)
        else:
            arr_mask_previous = load_mask_from_geometric_penalty(
                mask_filenames_list, indice=i - 1
            )
            overlap_mask = (
                arr_mask * arr_mask_i * np.int8(np.logical_not(arr_mask_previous))
            )
            res = (
                (arr_mask * ni.load(qmap_filenames_list[i + 1]).get_fdata())
                * np.logical_not(arr_mask * arr_mask_previous)
                + res
            ) * np.logical_not(overlap_mask)

        if modality_in_DMAP_VALUES:
            res += np.nan_to_num(
                fill_Qvalues_Overlapping_OnlySpatialRegularization(
                    qmap_filenames_list, mask_filenames_list, overlap_mask, i
                )
            )
        else:
            res += np.nan_to_num(
                fill_Qvalues_Overlapping_relaxometry(
                    cm_filenames_list,
                    qmap_filenames_list,
                    mask_filenames_list,
                    overlap_mask,
                    i,
                )
            )
    return res


def merge_QMRI_FOVs(
    json_filename: str | Path, modality_to_merge: str
) -> ni.Nifti1Image:
    fovs_filenames_to_merge_dict = load_jsonfile(json_filename)

    mask_filename_list = fovs_filenames_to_merge_dict["mask"]
    qmap_filename_list = fovs_filenames_to_merge_dict[modality_to_merge]

    is_in_DMAP_VALUES = modality_to_merge in DMAP_VALUES

    if not (is_in_DMAP_VALUES):
        confidence_map_filename_list = fovs_filenames_to_merge_dict[
            f"CM_{modality_to_merge[1:]}"
        ]
    else:
        confidence_map_filename_list = []

    qmap_block = reconstruction(
        confidence_map_filename_list,
        qmap_filename_list,
        mask_filename_list,
        is_in_DMAP_VALUES,
    )

    qmap_block = np.nan_to_num(qmap_block)
    meta = ni.load(mask_filename_list[0])
    return ni.Nifti1Image(qmap_block, meta.affine, meta.header)


def merge_QMRI_blocks(
    blocks_path: str | Path, modality: str, nbrBlocks: int
) -> ni.Nifti1Image:
    for i in range(nbrBlocks):
        block_filename = blocks_path / f"block_{str(i + 1)}_{modality}.nii.gz"
        if i < 1:
            meta_first_block = ni.load(block_filename)
            res = meta_first_block.get_fdata()
        else:
            res = np.maximum(res, ni.load(block_filename).get_fdata())
    return ni.Nifti1Image(res, meta_first_block.affine, meta_first_block.header)


def fovs_to_be_merged_are_specified(
    fov_material_filename: dict[str, list[str]], modality: str
) -> bool:
    return len(fov_material_filename.get(modality, [])) > 0


def run_Fusion_QMRI(working_directory: str | Path, nbr_of_blocks: int) -> None:
    blocks_path = Path(working_directory) / "03-Blocks"
    blocks_path.mkdir(parents=False, exist_ok=True)

    whole_hemisphere_path = Path(working_directory) / "04-Reconstruction"
    whole_hemisphere_path.mkdir(parents=False, exist_ok=True)

    for modality in QMAP_VALUES:
        for i in range(nbr_of_blocks):
            logger.info(
                f"========= Merge {modality} FOVs from block n° {i + 1}  ========="
            )

            block_json_filename = Path(working_directory) / f"block_{str(i + 1)}.json"
            fov_material_filename = load_jsonfile(block_json_filename)
            output_filename = blocks_path / f"block_{str(i + 1)}_{modality}.nii.gz"

            if fovs_to_be_merged_are_specified(fov_material_filename, modality):
                if not (output_filename.exists()):
                    block_reconstructed_nifti = merge_QMRI_FOVs(
                        block_json_filename, modality
                    )
                    ni.save(block_reconstructed_nifti, output_filename)
                else:
                    logger.info("Reconstructed block already exist.")
            else:
                logger.info(f"No filees found for {modality} reconstruction.")

        logger.info(f"Merge {modality} block {str(i + 1)}")
        reconstructed_hemisphere_filename = (
            whole_hemisphere_path / f"Reconstructed_{modality}.nii.gz"
        )
        if not (
            reconstructed_hemisphere_filename.exists()
        ) and fovs_to_be_merged_are_specified(fov_material_filename, modality):
            reconstructed_hemisphere_nifti = merge_QMRI_blocks(
                blocks_path, modality, nbr_of_blocks
            )
            ni.save(reconstructed_hemisphere_nifti, reconstructed_hemisphere_filename)
        else:
            logger.info("Whole hemisphere reconstruction already exist")


def RunPipeline(
    working_path: str | Path, nbrBlocks: int, run_merger_flag: bool
) -> None:
    if not (run_merger_flag):
        logger.info("========= Prepare the materials in the refspace =========")
        init_space_path = Path(working_path) / "01-InitSpace"
        ref_space_path = Path(working_path) / "02-RefSpace"

        fov_list = [
            file.name.split(".")[0].split("_")[-1]
            for file in init_space_path.glob("T2star_*.nii.gz")
            if not file.name.endswith("ConfidenceMap.nii.gz")
        ]
        if len(fov_list) != len(list(ref_space_path.glob("*geometric_penalty.nii.gz"))):
            run_FOV_mapping_and_geometric_penaty_creation(
                fov_list, init_space_path, ref_space_path, working_path
            )

        convert_gis_files_in_nifti_files(init_space_path)

        send_nifti_files_in_init_to_refspace(
            init_space_path, ref_space_path, working_path
        )

    else:
        logger.info("========= Merge FOVs in the refspace =========")
        run_Fusion_QMRI(working_path, nbrBlocks)


def parse_command_line(argv):
    """Parse the script's command line."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Reconstruction Script: Merge Fields of View (FOVs) into Final Space.\n\n"
        "Usage Modes:\n"
        "  - Default (no --run): prepares materials in 02-RefSpace\n"
        "  - With --run and -n: merges preprocessed blocks into final space\n\n"
        "Required Files:\n"
        "  - 01-InitSpace with modality_{FOV}.nii.gz files\n"
        "  - SendToRefSpace_{modality}.json files for each modality\n"
        "  - block_{i}.json files for each block to merge (if using --run)\n"
        "Please refer to the README.md file in the github for more details.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-p",
        "--path",
        required=True,
        help="Directory path of the working folder",
    )
    parser.add_argument(
        "-n",
        "--nbrBlocks",
        required=False,
        type=int,
        help="Number of blocks to be merged",
    )
    parser.add_argument(
        "-r",
        "--run",
        action="store_true",
        help="Run the second part of this pipeline if this flag is present",
    )

    args = parser.parse_args()
    return args


def main(argv=sys.argv):
    """The script's entry point."""
    logging.basicConfig(level=logging.INFO)
    args = parse_command_line(argv)
    return RunPipeline(args.path, args.nbrBlocks, args.run) or 0


if __name__ == "__main__":
    sys.exit(main())
