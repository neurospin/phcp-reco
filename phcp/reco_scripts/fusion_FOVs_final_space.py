"""
Created on Wed Jan 15 15:59:26 2025

@author: la272118
"""

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
    "T2w",
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
]
# QMAP_VALUES = ['QT1', 'QT2', 'QT2star', 'proton-density']
# CM_QMAP = ['CM_T1', 'CM_T2', 'CM_T2star', 'CM_T2star']
CM_QMAP = ["CM_T1", "CM_T2", "CM_T2star"]
DMAP_VALUES = [
    "T2w",
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
]  #'T2w',

""" Miscellaneous algorithms """


def print_message(message):
    print("==================================")
    print(message)
    print("==================================")


def GetMaskGkg(InputGisFilename, OutputGisFilename):
    subprocess.run(
        [
            "singularity",
            "exec",
            "--bind",
            "/neurospin:/neurospin:rw",
            "/neurospin/phcp/code/gkg/2022-12-20_gkg/2022-12-20_gkg.sif",
            "GkgExecuteCommand",
            "GetMask",
            "-i",
            InputGisFilename,
            "-o",
            OutputGisFilename,
            "-a",
            "1",
            "-format",
            "Gis",
            "-verbose",
        ]
    )


# Utilisation de cette fonction sur GIS file geomean pour avoir un masque propre
def create_ArrayFromGisFilename_16bit(GisFilename):
    DimFilename = GisFilename[:-3] + "dim"
    with open(DimFilename) as f:
        TxtInDimFile = f.read()
    TxtFirstLine = TxtInDimFile.split("\n")[0]
    shape = tuple(np.int16(TxtFirstLine.split(" ")))[:-1]
    arr = np.fromfile(GisFilename, np.float16())
    arr = np.reshape(arr, shape, order="F")
    return arr


def create_ArrayFromGisFilename_32bit(GisFilename):
    DimFilename = GisFilename[:-3] + "dim"
    with open(DimFilename) as f:
        TxtInDimFile = f.read()
    TxtFirstLine = TxtInDimFile.split("\n")[0]
    shape = tuple(np.int16(TxtFirstLine.split(" ")))[:-1]
    arr = np.fromfile(GisFilename, np.float32())
    arr = np.reshape(arr, shape, order="F")
    return arr


def arrange_ArrayToFlipHeaderFormat(arr):
    newarr = np.swapaxes(arr, 1, 2)
    newarr = np.flip(newarr, 1)
    newarr = np.flip(newarr, 0)
    return newarr


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
        nifti_filename = gis_filename.with_name(gis_filename.stem + "nii.gz")
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


def create_weight(CM_filename1, CM_filename2):
    CM_1 = ni.load(CM_filename1)
    CM_2 = ni.load(CM_filename2)

    CM_arr1 = proc.smooth_image(CM_1, proc.sigma2fwhm(0.2)).get_fdata()
    CM_arr2 = proc.smooth_image(CM_2, proc.sigma2fwhm(0.2)).get_fdata()

    return sigmoid(np.log(CM_arr1 / CM_arr2), 4, 0)


def create_mask_overlap_block(Mask_FilenameList):
    overlap_mask = np.zeros(ni.load(Mask_FilenameList[0]).shape, dtype=np.int8())
    for i in range(len(Mask_FilenameList) - 1):
        overlap_mask |= np.int8(ni.load(Mask_FilenameList[i]).get_fdata()) & np.int8(
            ni.load(Mask_FilenameList[i + 1]).get_fdata()
        )
    return overlap_mask


def fill_Qvalues_nonOverlapping(Mask_FilenameList, Qmap_FilenameList, overlap_mask):
    Qvalues_nonOverlapping = np.zeros(ni.load(Mask_FilenameList[0]).shape)
    for mask, Qmap in zip(Mask_FilenameList, Qmap_FilenameList, strict=False):
        Qvalues_nonOverlapping += (
            ni.load(mask).get_fdata() * ni.load(Qmap).get_fdata()
        ) * np.logical_not(overlap_mask)
    return Qvalues_nonOverlapping


def fill_Qvalues_Overlapping(CM_FilenameList, Qmap_FilenameList, overlap_mask):
    Qvalues_Overlapping = np.zeros(ni.load(CM_FilenameList[0]).shape)
    for i in range(len(CM_FilenameList) - 1):
        weight_2 = create_weight(CM_FilenameList[i], CM_FilenameList[i + 1])
        weight_1 = 1 - weight_2
        Qvalues_Overlapping += (
            weight_1 * ni.load(Qmap_FilenameList[i]).get_fdata()
            + weight_2 * ni.load(Qmap_FilenameList[i + 1]).get_fdata()
        ) * overlap_mask
    return Qvalues_Overlapping


def fill_Qvalues_Overlapping2(
    CM_FilenameList, Qmap_FilenameList, Mask_FilenameList, overlap_mask, i
):
    Qvalues_Overlapping = np.zeros(ni.load(CM_FilenameList[0]).shape)
    mask_weight_filename = Mask_FilenameList[i + 1].split(".")[0] + "_weight.nii.gz"
    meta_mask_weight = ni.load(mask_weight_filename)
    arr_mask_weight = meta_mask_weight.get_fdata()
    mask_weight_filename_i = Mask_FilenameList[i].split(".")[0] + "_weight.nii.gz"
    meta_mask_weight_i = ni.load(mask_weight_filename_i)
    arr_mask_weight_i = meta_mask_weight_i.get_fdata()

    weight_2 = (
        create_weight(CM_FilenameList[i], CM_FilenameList[i + 1]) * arr_mask_weight
    )
    weight_1 = (1 - weight_2) * arr_mask_weight_i
    weight_2 = 1 - weight_1
    Qvalues_Overlapping += (
        weight_1 * ni.load(Qmap_FilenameList[i]).get_fdata()
        + weight_2 * ni.load(Qmap_FilenameList[i + 1]).get_fdata()
    ) * overlap_mask
    return Qvalues_Overlapping


def fill_Qvalues_Overlapping_OnlySpatialRegularization(
    Qmap_FilenameList, Mask_FilenameList, overlap_mask, i
):
    Qvalues_Overlapping = np.zeros(ni.load(Qmap_FilenameList[0]).shape)
    mask_weight_filename = (
        Mask_FilenameList[i + 1].split(".")[0] + "_weight.nii.gz"
    )  #'_FOV_weight.nii.gz'
    meta_mask_weight = ni.load(mask_weight_filename)
    arr_mask_weight = meta_mask_weight.get_fdata()
    mask_weight_filename_i = (
        Mask_FilenameList[i].split(".")[0] + "_weight.nii.gz"
    )  #'_FOV_weight.nii.gz'
    meta_mask_weight_i = ni.load(mask_weight_filename_i)
    arr_mask_weight_i = meta_mask_weight_i.get_fdata()

    weight_2 = arr_mask_weight
    weight_1 = (1 - weight_2) * arr_mask_weight_i
    weight_2 = 1 - weight_1

    mask_weight_filename = Mask_FilenameList[i + 1]
    meta_mask_weight = ni.load(mask_weight_filename)
    arr_mask_weight = meta_mask_weight.get_fdata()
    arr_mask_weight = np.where(arr_mask_weight < 0.001, 0, 1)
    mask_weight_filename_i = Mask_FilenameList[i]
    meta_mask_weight_i = ni.load(mask_weight_filename_i)
    arr_mask_weight_i = meta_mask_weight_i.get_fdata()
    arr_mask_weight_i = np.where(arr_mask_weight_i < 0.001, 0, 1)

    weight_2 = weight_2 * arr_mask_weight
    weight_1 = weight_1 * arr_mask_weight_i

    Qvalues_Overlapping += (
        weight_1 * ni.load(Qmap_FilenameList[i]).get_fdata()
        + weight_2 * ni.load(Qmap_FilenameList[i + 1]).get_fdata()
    ) * overlap_mask
    return Qvalues_Overlapping


def reconstructionv2(CM_FilenameList, Qmap_FilenameList, Mask_FilenameList):
    res = (
        ni.load(Mask_FilenameList[0]).get_fdata()
        * ni.load(Qmap_FilenameList[0]).get_fdata()
    )
    for i in range(len(CM_FilenameList) - 1):
        meta_mask = ni.load(Mask_FilenameList[i + 1])
        arr_mask = meta_mask.get_fdata()
        meta_mask_i = ni.load(Mask_FilenameList[i])
        arr_mask_i = meta_mask_i.get_fdata()

        if i < 1:
            overlap_mask = arr_mask * arr_mask_i
            res = (
                (
                    (
                        ni.load(Mask_FilenameList[i + 1]).get_fdata()
                        * ni.load(Qmap_FilenameList[i + 1]).get_fdata()
                    )
                    + res
                )
                * np.logical_not(overlap_mask)
            )  # fill_Qvalues_nonOverlapping2(Mask_FilenameList, Qmap_FilenameList, res, overlap_mask, i)

        else:
            meta_mask_previous = ni.load(Mask_FilenameList[i - 1])
            arr_mask_previous = meta_mask_previous.get_fdata()
            overlap_mask = (
                arr_mask * arr_mask_i * np.int8(np.logical_not(arr_mask_previous))
            )  # case where mask[0] overlap with mask[1] and mask[2]
            res = (
                (
                    (
                        ni.load(Mask_FilenameList[i + 1]).get_fdata()
                        * ni.load(Qmap_FilenameList[i + 1]).get_fdata()
                    )
                    * np.logical_not(arr_mask * arr_mask_previous)
                    + res
                )
                * np.logical_not(overlap_mask)
            )  # fill_Qvalues_nonOverlapping2(Mask_FilenameList, Qmap_FilenameList, res, overlap_mask, i)

        res += np.nan_to_num(
            fill_Qvalues_Overlapping2(
                CM_FilenameList, Qmap_FilenameList, Mask_FilenameList, overlap_mask, i
            )
        )
    return res


def reconstructionv2_OnlySpatialRegularization(Qmap_FilenameList, Mask_FilenameList):
    meta_mask = ni.load(Mask_FilenameList[0].split(".")[0] + "_weight.nii.gz")
    arr_mask = meta_mask.get_fdata()
    arr_mask = np.where(arr_mask < 0.001, 0, 1)
    res = arr_mask * ni.load(Qmap_FilenameList[0]).get_fdata()
    for i in range(len(Qmap_FilenameList) - 1):
        meta_mask = ni.load(Mask_FilenameList[i + 1].split(".")[0] + "_weight.nii.gz")
        arr_mask = meta_mask.get_fdata()
        arr_mask = np.where(arr_mask < 0.001, 0, 1)
        meta_mask_i = ni.load(Mask_FilenameList[i].split(".")[0] + "_weight.nii.gz")
        arr_mask_i = meta_mask_i.get_fdata()
        arr_mask_i = np.where(arr_mask_i < 0.001, 0, 1)

        if i < 1:
            overlap_mask = arr_mask * arr_mask_i
            res = (
                ((arr_mask * ni.load(Qmap_FilenameList[i + 1]).get_fdata()) + res)
                * np.logical_not(overlap_mask)
            )  # fill_Qvalues_nonOverlapping2(Mask_FilenameList, Qmap_FilenameList, res, overlap_mask, i)
        else:
            meta_mask_previous = ni.load(
                Mask_FilenameList[i - 1].split(".")[0] + "_weight.nii.gz"
            )
            arr_mask_previous = meta_mask_previous.get_fdata()
            arr_mask_previous = np.where(arr_mask_previous < 0.001, 0, 1)
            overlap_mask = (
                arr_mask * arr_mask_i * np.int8(np.logical_not(arr_mask_previous))
            )  # case where mask[0] overlap with mask[1] and mask[2]
            res = (
                (
                    (arr_mask * ni.load(Qmap_FilenameList[i + 1]).get_fdata())
                    * np.logical_not(arr_mask * arr_mask_previous)
                    + res
                )
                * np.logical_not(overlap_mask)
            )  # fill_Qvalues_nonOverlapping2(Mask_FilenameList, Qmap_FilenameList, res, overlap_mask, i)

        res += np.nan_to_num(
            fill_Qvalues_Overlapping_OnlySpatialRegularization(
                Qmap_FilenameList, Mask_FilenameList, overlap_mask, i
            )
        )
    return res


def merge_QMRI_FOVs(
    json_filename: str | Path, modality_to_merge: str
) -> ni.Nifti1Image:
    fovs_filenames_to_merge_dict = load_jsonfile(json_filename)

    mask_filename_list = fovs_filenames_to_merge_dict["mask"]
    qmap_filename_list = fovs_filenames_to_merge_dict[modality_to_merge]

    if modality_to_merge not in DMAP_VALUES:
        confidence_map_filename_list = fovs_filenames_to_merge_dict[
            CM_QMAP[QMAP_VALUES.index(modality_to_merge)]
        ]

    meta = ni.load(mask_filename_list[0])

    if modality_to_merge not in DMAP_VALUES:
        qmap_block = reconstructionv2(
            confidence_map_filename_list, qmap_filename_list, mask_filename_list
        )
    else:
        qmap_block = reconstructionv2_OnlySpatialRegularization(
            qmap_filename_list, mask_filename_list
        )

    qmap_block = np.nan_to_num(qmap_block)
    return ni.Nifti1Image(qmap_block, meta.affine, meta.header)


def merge_QMRI_blocks(
    blocks_path: str | Path, modality: str, nbrBlocks: int
) -> ni.Nifti1Image:
    for i in range(nbrBlocks):
        block_filename = blocks_path / f"bloc_str(i + 1)_{modality}.nii.gz"
        if i < 1:
            meta_first_block = ni.load(block_filename)
            res = meta_first_block.get_fdata()
        else:
            res = np.maximum(res, ni.load(block_filename).get_fdata())
    return ni.Nifti1Image(res, meta_first_block.affine, meta_first_block.header)


def fovs_to_be_merged_are_specified(
    fov_material_filename: str | Path, modality: str
) -> bool:
    return len(fov_material_filename[modality]) > 0


def run_Fusion_QMRI(working_directory: str | Path, nbr_of_blocks: int) -> None:
    blocks_path = Path(working_directory) / "03-Blocks"
    blocks_path.mkdir(parents=False, exist_ok=True)

    whole_hemisphere_path = Path(working_directory) / "04-Reconstruction"
    whole_hemisphere_path.mkdir(parents=False, exists_ok=True)

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
        description="Merge Field of views in the final space.",
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
