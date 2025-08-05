import logging
import subprocess

logger = logging.getLogger(__name__)


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
