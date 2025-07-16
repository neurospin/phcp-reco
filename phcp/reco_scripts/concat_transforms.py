import json
import logging
import subprocess
import sys
from pathlib import Path

import ants
import nibabel as ni
import numpy as np
import SimpleITK as sitk
from scipy.ndimage import laplace

logger = logging.getLogger(__name__)


"""Create Reference Init Space Files"""


def create_ReferenceInitSpace(filename: str, OutputDirectory: str | Path) -> None:
    meta = ni.load(filename)
    xdim, ydim, zdim = meta.shape
    newarr = np.meshgrid(
        np.arange(xdim), np.arange(ydim), np.arange(zdim), indexing="ij"
    )
    ni.save(
        ni.Nifti1Image(newarr[0], meta.affine, meta.header),
        Path(OutputDirectory) / "ReferenceInitSpace_x.nii.gz",
    )
    ni.save(
        ni.Nifti1Image(newarr[1], meta.affine, meta.header),
        Path(OutputDirectory) / "ReferenceInitSpace_y.nii.gz",
    )
    ni.save(
        ni.Nifti1Image(newarr[2], meta.affine, meta.header),
        Path(OutputDirectory) / "ReferenceInitSpace_z.nii.gz",
    )

    # Generates FOV mapping. Will be used to restrict the distorted FOV.
    ni.save(
        ni.Nifti1Image(np.ones(meta.shape), meta.affine, meta.header),
        Path(OutputDirectory) / "mask.nii.gz",
    )


"""Load & Create Transform Files"""


def Create_transformationFilesList_from_JsonFilename(JsonFilename: str) -> list[str]:
    with open(JsonFilename) as f:
        res = json.load(f)
    return res["tparams"]


def is_mat_file(path: str | Path) -> bool:
    return Path(path).suffix == ".mat"


def is_warp_file(path: str | Path) -> bool:
    return Path(path).stem.lower().endswith("warp")


def make_complete_transformationlist_with_invertflags(
    TransformationFilesList: list[str],
) -> (list[str], list[bool]):
    """
    Constructs a complete list of transformation filenames and corresponding inversion flags,
    ordered to allow composition of all transformations into the initial reference space.

    The strategy applied follows this sequence:

        Initial space
             |
        Linear_1
             |
        non_linear_1
             |
        Linear_2
             |
        non_linear_2
             |
            ...
             |
        Last_transform

    To bring everything back to the initial space, we apply the **inverse** of each linear transform in reverse order,
    skipping all the non-linear transformation.

    This results in the following application order:

        Inverse Linear_1
               |
        Inverse Linear_2
               |
            ...
               |
        Last_transform (if linear, then inverted, otherwise not added in this part)

    Returns:
        complete_transformationlist (list[str]): Ordered list of transformation file paths.
        invert_flagslist (list[bool]): Boolean flags indicating whether each corresponding transformation should be inverted.
    """
    complete_transformationlist = []
    invert_flagslist = []
    for el in TransformationFilesList:
        if ".nii" not in Path(el).suffixes:
            complete_transformationlist.append(el)
            invert_flagslist.append(True)
    TransformationFilesList.reverse()
    for el in TransformationFilesList:
        complete_transformationlist.append(el)
        invert_flagslist.append(False)
    TransformationFilesList.reverse()
    return complete_transformationlist, invert_flagslist


"""Apply Transform Files"""


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


def apply_transformationlist_To_ReferenceInitSpace_Files(
    InputFilename: str,
    OutputDirectory: str | Path,
    TransformationFilesList: list[str],
    invert_flag_list: list[bool],
) -> None:
    movingFilename = (Path(OutputDirectory) / "ReferenceInitSpace_x.nii.gz").as_posix()
    OutputRefWarpedxFilename = (
        Path(OutputDirectory) / "ReferenceInitSpaceWarped_x.nii.gz"
    ).as_posix()

    run_ants_apply_registration(
        ref_space=InputFilename,
        input_img=movingFilename,
        output_filename=OutputRefWarpedxFilename,
        transforms=TransformationFilesList,
        dimensionality=3,
        invert_flags=invert_flag_list,
    )

    movingFilename = (Path(OutputDirectory) / "ReferenceInitSpace_y.nii.gz").as_posix()
    OutputRefWarpedxFilename = (
        Path(OutputDirectory) / "ReferenceInitSpaceWarped_y.nii.gz"
    ).as_posix()

    run_ants_apply_registration(
        InputFilename,
        movingFilename,
        OutputRefWarpedxFilename,
        TransformationFilesList,
        dimensionality=3,
        invert_flags=invert_flag_list,
    )

    movingFilename = (Path(OutputDirectory) / "ReferenceInitSpace_z.nii.gz").as_posix()
    OutputRefWarpedxFilename = (
        Path(OutputDirectory) / "ReferenceInitSpaceWarped_z.nii.gz"
    ).as_posix()

    run_ants_apply_registration(
        InputFilename,
        movingFilename,
        OutputRefWarpedxFilename,
        TransformationFilesList,
        dimensionality=3,
        invert_flags=invert_flag_list,
    )

    movingFilename = (Path(OutputDirectory) / "mask.nii.gz").as_posix()
    OutputRefWarpedxFilename = (Path(OutputDirectory) / "MaskWarped.nii.gz").as_posix()
    run_ants_apply_registration(
        InputFilename,
        movingFilename,
        OutputRefWarpedxFilename,
        TransformationFilesList,
        interpolation="NearestNeighbor",
        dimensionality=3,
        invert_flags=invert_flag_list,
    )


"""Create Deformation Field"""


def change_header_to_deformation_fieldheader(ni_im: ni.Nifti1Image) -> ni.Nifti1Image:
    ni_im.header.set_intent(1007)
    ni_im.header["regular"] = b"r"
    ni_im.header.set_xyzt_units(2)
    sform = ni_im.header.get_sform()
    ni_im.header.set_qform(sform, code=1)
    ni_im.header.set_sform(sform, code=1)
    return ni_im


def create_deformation_field(Filename: str | Path, OutputDirectory: str | Path) -> None:
    meta = ni.load(Filename)
    xdim, ydim, zdim = meta.shape
    resx, resy, resz = meta.header["pixdim"][1:4]
    metax = ni.load(Path(OutputDirectory) / "ReferenceInitSpaceWarped_x.nii.gz")
    metay = ni.load(Path(OutputDirectory) / "ReferenceInitSpaceWarped_y.nii.gz")
    metaz = ni.load(Path(OutputDirectory) / "ReferenceInitSpaceWarped_z.nii.gz")
    meta_mask = ni.load(Path(OutputDirectory) / "MaskWarped.nii.gz")

    newarr = np.meshgrid(
        np.arange(xdim), np.arange(ydim), np.arange(zdim), indexing="ij"
    )
    arr_mask = meta_mask.get_fdata()

    arr_Ref = np.concatenate(
        (
            -resx * np.reshape(arr_mask * newarr[0], meta.shape + (1, 1)),
            resz * np.reshape(arr_mask * newarr[2], meta.shape + (1, 1)),
            resy * np.reshape(arr_mask * newarr[1], meta.shape + (1, 1)),
        ),
        -1,
    )
    arr_SyN = np.concatenate(
        (
            -resx * np.reshape(arr_mask * metax.get_fdata(), (metax.shape + (1, 1))),
            resz * np.reshape(arr_mask * metaz.get_fdata(), (metaz.shape + (1, 1))),
            resy * np.reshape(arr_mask * metay.get_fdata(), (metay.shape + (1, 1))),
        ),
        -1,
    )

    res = arr_SyN - arr_Ref

    ni_im = ni.Nifti1Image(res, metax.affine, metax.header)
    ni_im = change_header_to_deformation_fieldheader(ni_im)
    ni.save(ni_im, Path(OutputDirectory) / "Deformation_field_SyN.nii.gz")


""" Create Jacobian files """


def create_jacobian_files(Filename: str | Path, OutputDirectory: str | Path) -> None:
    ants_ref = ants.image_read(Filename)
    jacobian = ants.create_jacobian_determinant_image(
        ants_ref,
        Path(OutputDirectory) / "Deformation_field_SyN.nii.gz",
        do_log=False,
    )
    ants.image_write(
        jacobian,
        (Path(OutputDirectory) / "Jacobian_Deformation_field_SyN.nii.gz").as_posix(),
    )
    jacobian = ants.create_jacobian_determinant_image(
        ants_ref,
        Path(OutputDirectory) / "Deformation_field_SyN.nii.gz",
        do_log=True,
    )
    ants.image_write(
        jacobian,
        (Path(OutputDirectory) / "JacobianLog_Deformation_field_SyN.nii.gz").as_posix(),
    )


"""Create TotalAffineTransform File"""


def is_nifiti_file(path: str | Path) -> bool:
    return ".nii" in Path(path).suffixes


def ExtractMatrixTransforms(TransformationFilesList: list[str]) -> list[str]:
    MatrixTransformsFilenameList = list()
    for el in TransformationFilesList:
        if not (is_nifiti_file(el)):
            MatrixTransformsFilenameList.append(el)
    return MatrixTransformsFilenameList


def create_matrix_matFile(matrix: tuple, trans: tuple) -> np.ndarray:
    matrix, trans = (
        np.asarray(matrix).reshape((3, 3)),
        np.asarray(trans).reshape((3, 1)),
    )
    ligne = np.asarray([0, 0, 0, 1]).reshape((1, 4))
    hst = np.hstack((matrix, trans.reshape((3, 1))))
    return np.vstack((hst, ligne))


def create_matrix_txtFile(parameters: tuple) -> np.ndarray:
    matrix, trans = (
        np.asarray(parameters)[:9].reshape((3, 3)),
        np.asarray(parameters)[9:].reshape((3, 1)),
    )
    ligne = np.asarray([0, 0, 0, 1]).reshape((1, 4))
    hst = np.hstack((matrix, trans.reshape((3, 1))))
    return np.vstack((hst, ligne))


def write_total_affine_transform_txtformat(
    output_filename: str | Path, matrix_array: np.ndarray
) -> None:
    with open(output_filename, "w") as f:
        f.write("#Insight Transform File V1.0\n")
        f.write("#Transform 0\n")
        f.write("Transform: MatrixOffsetTransformBase_double_3_3\n")
        f.write("Parameters: ")
        f.write(
            " ".join(
                map(
                    str,
                    np.concatenate(
                        (
                            np.asarray(matrix_array[:3, :3]).flatten(),
                            np.asarray(matrix_array[:3, 3].reshape((3,))),
                        )
                    ),
                )
            )
        )
        f.write("\n")
        f.write("FixedParameters: 0 0 0")  # +" ".join(map(str, center)))


def compose_transformations(
    MatrixTransformsFilenameList: list[str], OutputDirectory: str | Path
) -> None:
    res = np.eye(4)
    center = np.asarray((0, 0, 0))
    for trans in MatrixTransformsFilenameList:
        information = sitk.ReadTransform(trans)

        if is_mat_file(trans):
            matrix = information.GetMatrix()
            trans = information.GetTranslation()
            if information.GetCenter() != list((0, 0, 0)):
                center = np.asarray(information.GetCenter())
            Euler = create_matrix_matFile(matrix, trans)
        else:
            parameters = information.GetParameters()
            Euler = create_matrix_txtFile(parameters)

        res = res @ Euler
    T_center, T_neg_center = (
        np.eye(4),
        np.eye(4),
    )
    T_center[:3, 3], T_neg_center[:3, 3] = center, -center
    res = (
        T_center @ res @ T_neg_center
    )  # Including the center of rotation in the affine matrix -> implies that the center of rotation is zero

    OutputFilename = Path(OutputDirectory) / "TotalAffineTransform.txt"
    write_total_affine_transform_txtformat(OutputFilename, res)


def regularize_laplacian_smoothing(
    deformationField_array: np.ndarray, iteration: int, Lambda: float, component: int
) -> (np.ndarray, np.ndarray):
    arr = deformationField_array[:, :, :, 0, component]
    for _ in range(iteration):
        arr += Lambda * laplace(arr)
    return arr, laplace(arr)


def apply_laplacian_smoothing(OutputDirectory: str | Path) -> None:
    meta_deformationField = ni.load(
        Path(OutputDirectory) / "Deformation_field_SyN.nii.gz"
    )
    arr_deformationField = meta_deformationField.get_fdata()

    arrsmoothed_deformationfield_component_x, Laplacian_deformationfield_component_x = (
        regularize_laplacian_smoothing(
            arr_deformationField, iteration=3, Lambda=0.2, component=0
        )
    )
    ni_im = ni.Nifti1Image(
        Laplacian_deformationfield_component_x, meta_deformationField.affine
    )
    ni.save(ni_im, Path(OutputDirectory) / "Laplacian_component_x.nii.gz")

    arrsmoothed_deformationfield_component_y, Laplacian_deformationfield_component_y = (
        regularize_laplacian_smoothing(
            arr_deformationField, iteration=3, Lambda=0.2, component=1
        )
    )
    ni_im = ni.Nifti1Image(
        Laplacian_deformationfield_component_y, meta_deformationField.affine
    )
    ni.save(ni_im, Path(OutputDirectory) / "Laplacian_component_y.nii.gz")

    arrsmoothed_deformationfield_component_z, Laplacian_deformationfield_component_z = (
        regularize_laplacian_smoothing(
            arr_deformationField, iteration=3, Lambda=0.2, component=2
        )
    )
    ni_im = ni.Nifti1Image(
        Laplacian_deformationfield_component_z, meta_deformationField.affine
    )
    ni.save(ni_im, Path(OutputDirectory) / "Laplacian_component_z.nii.gz")

    SmoothedDeformationField = np.stack(
        [
            arrsmoothed_deformationfield_component_x,
            arrsmoothed_deformationfield_component_y,
            arrsmoothed_deformationfield_component_z,
        ],
        axis=-1,
    )

    ni_im = ni.Nifti1Image(
        SmoothedDeformationField,
        meta_deformationField.affine,
        meta_deformationField.header,
    )
    ni.save(ni_im, Path(OutputDirectory) / "Deformation_field_SyN_Smoothed.nii.gz")


def concat_transforms(InputFilename, JsonFilename, OutputDirectory):
    # logger.info("========= Create Reference Init Space Files =========")
    # create_ReferenceInitSpace(InputFilename, OutputDirectory)

    # logger.info("========= Load & Create Transform Files =========")
    # TransformationFilesList = Create_transformationFilesList_from_JsonFilename(
    #     JsonFilename
    # )
    # MatrixTransformsFilenameList = ExtractMatrixTransforms(TransformationFilesList)

    # logger.info("========= Create TotalAffineTransform File =========")
    # compose_transformations(MatrixTransformsFilenameList, OutputDirectory)

    # NewTransformationFilesList, InvertFlagList = (
    #     make_complete_transformationlist_with_invertflags(TransformationFilesList)
    # )

    # logger.info("========= Apply Transform Files =========")
    # apply_transformationlist_To_ReferenceInitSpace_Files(
    #     InputFilename, OutputDirectory, NewTransformationFilesList, InvertFlagList
    # )

    # logger.info("========= Create Deformation Field =========")
    # create_deformation_field(InputFilename, OutputDirectory)

    logger.info("========= Create Jacobian Determinant files =========")
    create_jacobian_files(InputFilename, OutputDirectory)

    logger.info("========= Deformation Field - Laplacian Smoothing =========")
    apply_laplacian_smoothing(OutputDirectory)
    return None


def parse_command_line(argv):
    """Parse the script's command line."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Concatenates linear and non-linear transformations in ANTS format. \n"
        "This script creates the following main files in the output directory :\n"
        "    + Deformation_field_SyN.nii.gz, Concatenated raw non-linear transformations ;\n"
        "    + Deformation_field_SyN_Smoothed.nii.gz, SyN_deformation_field smoothed with Laplacian model (can be used in future scripts) ;\n"
        "    + Jacobian_Deformation_field_SyN.nii.gz, Jacobian determinant of Deformation_field_SyN.nii.gz ;\n"
        "    + JacobianLog_Deformation_field_SyN.nii.gz, log of Jacobian_Deformation_field_SyN.nii.gz;\n"
        "    + TotalAffineTransform.txt, Concatenated linear transformations.\n \n"
        "USAGE : \n    Transformation files should be applied in the following order : \n"
        "        1. Deformation_field_SyN.nii.gz or Deformation_field_SyN_Smoothed.nii.gz \n"
        "        2. TotalAffineTransform.txt.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Filename of the FOV in the initial space (ex: sub-{sub}_"
        "ses-_{ses}_T2w.nii.gz)",
    )
    parser.add_argument(
        "-j",
        "--json",
        required=True,
        help="Json filename containing path of transformation files from "
        "initial space to the final space. The json filename should followed "
        "the following architecture {  'tparams':  ["
        "transform_file_1_path, transform_file_2_path, ...,"
        " transform_file_n_path]   }",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="output directory",
    )

    args = parser.parse_args()
    return args


def main(argv=sys.argv):
    """The script's entry point."""
    logging.basicConfig(level=logging.INFO)
    args = parse_command_line(argv)
    return concat_transforms(args.input, args.json, args.output) or 0


if __name__ == "__main__":
    sys.exit(main())
