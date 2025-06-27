#! /usr/bin/env python3

import json
import logging
import os
import pathlib
import re
import sys

import nibabel
import numpy
import yaml

logger = logging.getLogger(__name__)


def nibabel_orient_as_RAS(
    img: nibabel.spatialimages.SpatialImage,
) -> nibabel.spatialimages.SpatialImage:
    img_ornt = nibabel.orientations.io_orientation(img.affine)
    return img.as_reoriented(
        nibabel.orientations.ornt_transform(img_ornt, [[0, 1], [1, 1], [2, 1]])
    )


def convert_to_scaled_encoding(img, dtype, slope, inter, reset_scaling=False):
    arr = img.get_fdata(caching="unchanged")
    iinfo = numpy.iinfo(dtype)

    # The header may not be able to store full-precision values (e.g. Nifti1 is
    # float32 vs float64). This round-trip applies the rounding, so that the
    # actual stored values are used when encoding.
    tmp_header = img.header.copy()
    tmp_header.set_slope_inter(slope, inter)
    slope, inter = tmp_header.get_slope_inter()

    min_encoded = (iinfo.min * slope) + inter
    max_encoded = (iinfo.max * slope) + inter
    print(f"Data range min-max: {min_encoded} to {max_encoded}")
    print(f"Resolution: {slope}")
    new_arr = arr.copy()
    mask = numpy.isnan(arr)
    new_arr[mask] = 0
    print(
        f"{numpy.count_nonzero(mask)} NaN values encoded by zero (scaled value {inter})"
    )
    mask = arr < min_encoded
    new_arr[mask] = min_encoded
    print(f"{numpy.count_nonzero(mask)} values < {min_encoded} replaced by {iinfo.min}")
    mask = arr > max_encoded
    new_arr[mask] = max_encoded
    print(f"{numpy.count_nonzero(mask)} values > {max_encoded} replaced by {iinfo.max}")
    new_arr = numpy.rint((new_arr - inter) / slope).astype(dtype)
    new_header = img.header.copy()
    new_header.set_data_dtype(dtype)
    new_img = nibabel.Nifti1Image(new_arr, img.affine, new_header)
    # Must be set AFTER the creation of the Nifti1Image (slope/inter are RESET
    # by the constructor of Nifti1Image)
    if not reset_scaling:
        new_img.header.set_slope_inter(slope, inter)
    return new_img


def deidentify_json(json_dict: dict):
    # Remove unimportant or possibly identifying information
    json_dict.pop("DeviceSerialNumber", None)
    json_dict.pop("StationName", None)
    json_dict.pop("ScanDate", None)
    json_dict.pop("BidsGuess", None)

    # Remove the software version of the DICOM deidentification software, which
    # is not crucial but could leak timing information.
    if "DeidentificationMethod" in json_dict:
        json_dict["DeidentificationMethod"] = [
            (
                "CATI DEIDENTIFICATION - https://github.com/cati-neuroimaging/deidentification"
                if "CATI DEIDENTIFICATION" in row
                else row
            )
            for row in json_dict["DeidentificationMethod"]
        ]
    return json_dict


def encode_image(
    image_dict: dict,
    config: dict,
    input_dir: os.PathLike,
    output_dir: os.PathLike,
    *,
    dry_run=False,
):
    source_sub = config["source_sub"]
    target_sub = config["target_sub"]
    output_suffix = image_dict["suffix"]
    if not image_dict.get("source"):
        logger.error("A source is missing for %s", output_suffix)
        return
    input_subpath = pathlib.PurePosixPath(image_dict["source"])
    assert not input_subpath.is_absolute()
    output_datatype = image_dict["datatype"]
    input_full_path = input_dir / f"sub-{source_sub}" / input_subpath
    input_json_path = input_full_path.with_name(
        re.sub(r"\.nii(\.gz)?$", ".json", input_full_path.name)
    )
    output_basename = f"sub-{target_sub}_{output_suffix}"
    output_full_path = (
        output_dir
        / f"sub-{target_sub}"
        / output_datatype
        / (output_basename + ".nii.gz")
    )
    output_json_path = output_full_path.with_name(output_basename + ".json")
    if not input_full_path.is_file():
        logger.error("%s does not exist or is not a file", input_full_path)
        return

    if input_json_path.exists():
        with input_json_path.open() as f:
            json_dict = json.load(f)
        json_dict = deidentify_json(json_dict)
        if output_json_path.exists() or output_json_path.is_symlink():
            logger.error(
                "Target file already exists, not overwriting: %s", output_json_path
            )
        elif dry_run:
            logger.info("would copy %s to %s", input_json_path, output_json_path)
        else:
            logger.info("copying %s to %s", input_json_path, output_json_path)
            with output_json_path.open("w") as f:
                json.dump(json_dict, f, indent=4)
                f.write("\n")

    if output_full_path.exists() or output_full_path.is_symlink():
        logger.error(
            "Target file already exists, not overwriting: %s", output_full_path
        )
    elif dry_run:
        logger.info("would encode %s into %s", input_full_path, output_full_path)
    else:
        logger.info("encoding %s into %s", input_full_path, output_full_path)
        output_full_path.parent.mkdir(parents=True, exist_ok=True)

        encoding_dict = config["encodings"].get(image_dict.get("encoding"), {})
        geometry_dict = config["geometries"].get(
            image_dict.get("geometry", "default"), {}
        )
        if geometry_dict is None:
            geometry_dict = {}

        img = nibabel_orient_as_RAS(nibabel.load(input_full_path))

        if "expected_shape" in geometry_dict:
            expected_shape = geometry_dict["expected_shape"]
            shape = img.header.get_data_shape()
            if not numpy.array_equiv(shape, expected_shape):
                logger.error(
                    "Skipping image with unexpected shape %s (expected %s)",
                    shape,
                    expected_shape,
                )
                return

        if "affine_corr" in geometry_dict:
            affine = numpy.asarray(geometry_dict["affine_corr"]) @ img.affine
            img = nibabel.Nifti1Image(img.dataobj, affine, img.header)

        if "crop" in geometry_dict:
            crop_slicing = tuple(
                slice(min, max, None) for min, max in geometry_dict["crop"]
            )
            img = img.slicer[crop_slicing]
        else:
            img = img

        if "affine_precision_override" in geometry_dict:
            affine = numpy.asarray(geometry_dict["affine_precision_override"])
            if numpy.allclose(
                affine[:3, :3], img.affine[:3, :3], rtol=1e-4
            ) and numpy.allclose(affine[:3, 3], img.affine[:3, 3], atol=0.05):
                img = nibabel.Nifti1Image(img.dataobj, affine, img.header)
            else:
                logger.error(
                    "Refusing to override the affine matrix, "
                    "because it does not match the recorded transform:\n%s",
                    img.affine,
                )
                return
            voxel_size = numpy.linalg.norm(affine[:3, :3], axis=0)
            if not numpy.allclose(voxel_size, img.header.get_zooms(), rtol=1e-4):
                logger.error(
                    "Refusing to override the voxel size, "
                    "because it does not match the recorded voxel size:\n%s",
                    img.header.get_zooms(),
                )
                return
            img.header.set_zooms(voxel_size)
            img.header.set_qform(affine, code="aligned")
            img.header.set_sform(affine, code="aligned")

        if (
            "dtype" in encoding_dict
            or "slope" in encoding_dict
            or "inter" in encoding_dict
        ):
            img = convert_to_scaled_encoding(
                img,
                numpy.dtype(encoding_dict["dtype"]),
                float(encoding_dict["slope"]),
                float(encoding_dict["inter"]),
                reset_scaling=encoding_dict.get("reset_scaling", False),
            )
        image_dict.setdefault("intent", encoding_dict.get("intent", 0))
        image_dict.setdefault("intent_name", encoding_dict.get("intent_name", ""))
        img.header.set_intent(image_dict["intent"], name=image_dict["intent_name"])
        image_dict.setdefault("cal_min", encoding_dict.get("cal_min"))
        image_dict.setdefault("cal_max", encoding_dict.get("cal_max"))
        image_dict.setdefault("cal_max_quantile", encoding_dict.get("cal_max_quantile"))
        img.header["cal_min"] = image_dict["cal_min"]
        if image_dict["cal_max"] == "max":
            img.header["cal_max"] = img.get_fdata().max()
        elif image_dict["cal_max"] is not None:
            img.header["cal_max"] = image_dict["cal_max"]
        elif image_dict["cal_max_quantile"] is not None:
            img.header["cal_max"] = numpy.quantile(
                img.get_fdata(), image_dict["cal_max_quantile"]
            )
        full_descrip = (output_basename + " " + config.get("descrip_post", "")).encode()
        if len(full_descrip) > 79:
            descrip_crop_len = len(full_descrip) - 79 + 3
            img.header["descrip"] = (
                output_basename[:-descrip_crop_len]
                + "... "
                + config.get("descrip_post", "")
            ).encode()
        else:
            img.header["descrip"] = full_descrip
        img.header.set_xyzt_units(xyz="mm", t=None)
        img.header.set_dim_info(None, None, None)
        # img.header.set_slice_duration(None)
        # img.header.set_slice_times(None)
        img.header["toffset"] = 0.0
        img.header["aux_file"] = b""

        nibabel.save(img, output_full_path)


def prepare_published_dataset(
    input_dir: os.PathLike | str,
    output_dir: os.PathLike | str,
    config_yaml: os.PathLike | str,
    dry_run=False,
):
    with open(config_yaml) as f:
        config = yaml.safe_load(f)

    for image_dict in config["images"]:
        encode_image(image_dict, config, input_dir, output_dir, dry_run=dry_run)

    # Naming of the diffusion maps is based on this version of BIDS Enhancement Proposal (BEP16):
    # https://github.com/bids-standard/bids-bep016/blob/edab0fd74b0530ffd4d35258c2e7dc3911262fce/src/derivatives/05-diffusion-derivatives.md

    # with (output_dir / "dataset_description.json").open("w") as f:
    #     json.dump({"Name": "p-HCP", "BIDSVersion": "1.10.0"}, f, indent=4)


def parse_command_line(argv):
    """Parse the script's command line."""
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=pathlib.Path)
    parser.add_argument("output_dir", type=pathlib.Path)
    parser.add_argument("config_yaml", type=pathlib.Path)
    parser.add_argument("--dry-run", action="store_true")

    args = parser.parse_args()
    return args


def main(argv=sys.argv):
    """The script's entry point."""
    logging.basicConfig(level=logging.INFO)
    args = parse_command_line(argv)
    return (
        prepare_published_dataset(
            args.input_dir, args.output_dir, args.config_yaml, dry_run=args.dry_run
        )
        or 0
    )


if __name__ == "__main__":
    sys.exit(main())
