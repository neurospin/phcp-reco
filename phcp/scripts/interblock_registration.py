#!/usr/bin/env python3

import logging
import sys
from collections.abc import Sequence
from pathlib import Path

import nibabel as ni
import numpy as np
import scipy.ndimage as nd
from scipy.ndimage import map_coordinates

logger = logging.getLogger(__name__)


def spherical_structure(radius):
    size = 2 * radius + 1
    center = radius
    X, Y, Z = np.ogrid[:size, :size, :size]
    sphere = (X - center) ** 2 + (Y - center) ** 2 + (Z - center) ** 2
    return sphere.astype(np.uint8)


def extract_hull(
    segmentation_filename: str | Path, output_filename: str | Path
) -> None:
    meta = ni.load(segmentation_filename)
    arr = meta.get_fdata()
    structure = spherical_structure(3)
    erode = nd.binary_erosion(arr, structure=structure, iterations=2)
    dilate = nd.binary_dilation(arr, structure=structure, iterations=2)
    newarr = dilate ^ erode
    ni_im = ni.Nifti1Image(newarr, meta.affine, meta.header)
    ni.save(ni_im, output_filename)


def plane_extractor_SVD(filename: str | Path) -> None:
    meta = ni.load(filename)
    arr = meta.get_fdata()
    points = np.argwhere(arr > 0)

    logger.info("Downsampled part")
    sampling_ratio = 0.005
    indices = np.arange(len(points))
    np.random.shuffle(indices)
    sample_size = int(len(points) * sampling_ratio)
    sampled_indices = indices[:sample_size]
    downsampled_points = points[sampled_indices]

    centroid = np.mean(downsampled_points, axis=0)
    centered_points = downsampled_points - centroid

    logger.info("SVD part")
    _, _, vh = np.linalg.svd(centered_points)
    normal_vector = vh[-1, :]

    logger.info("Centroid =" + str(centroid))
    logger.info("Normal vector =" + str(normal_vector))


def get_2dplane(
    block_filename: str | Path,
    output_filename: str | Path,
    centroid: Sequence[float],
    normal_vector: Sequence[float],
    usexz: bool,
    usexy: bool,
    useyz: bool,
) -> None:
    meta = ni.load(block_filename)
    arr = meta.get_fdata()

    centroid = np.asarray(centroid, dtype=float)
    normal_vector = np.asarray(normal_vector, dtype=float)

    d = -centroid.dot(normal_vector)

    x_dim, y_dim, z_dim = arr.shape

    if usexz:
        x = np.linspace(0, x_dim - 1, x_dim)
        z = np.linspace(0, z_dim - 1, z_dim)

        xx, zz = np.meshgrid(x, z)

        yy = (-normal_vector[0] * xx - normal_vector[2] * zz - d) / normal_vector[1]

    if usexy:
        x = np.linspace(0, x_dim - 1, x_dim)
        y = np.linspace(0, y_dim - 1, y_dim)

        xx, yy = np.meshgrid(x, y)

        zz = (-normal_vector[0] * xx - normal_vector[1] * yy - d) / normal_vector[2]

    if useyz:
        z = np.linspace(0, z_dim - 1, z_dim)
        y = np.linspace(0, y_dim - 1, y_dim)

        yy, zz = np.meshgrid(y, z)

        xx = (-normal_vector[1] - normal_vector[2] * zz * yy - d) / normal_vector[0]

    # Stack the coordinates for interpolation
    coords = np.stack([xx, yy, zz], axis=-1)

    # Interpolate the values at the grid points
    plane_2d_image = map_coordinates(
        arr, [coords[..., 0], coords[..., 1], coords[..., 2]], order=3, mode="nearest"
    )

    affine = np.eye(4) * meta.header["pixdim"][1]
    affine[3][3] = 1

    ni.save(ni.Nifti1Image(plane_2d_image, affine), output_filename)


def back_proj(
    block_filename: str | Path,
    slice_to_project_filename: str | Path,
    output_filename: str | Path,
    centroid: Sequence[float],
    normal_vector: Sequence[float],
    usexz: bool,
    usexy: bool,
    useyz: bool,
) -> None:
    meta = ni.load(block_filename)
    arr = meta.get_fdata()

    meta_seg = ni.load(slice_to_project_filename)
    modified_plane_2d_image = meta_seg.get_fdata()

    centroid = np.asarray(centroid, dtype=float)
    normal_vector = np.asarray(normal_vector, dtype=float)

    d = -centroid.dot(normal_vector)

    # Define the bounding box for the grid points
    x_dim, y_dim, z_dim = arr.shape

    # Map the modified 2D plane data back to 3D space
    modified_3d_plane = np.zeros_like(arr)

    if usexz:
        x = np.linspace(0, x_dim - 1, x_dim)
        z = np.linspace(0, z_dim - 1, z_dim)

        xx, zz = np.meshgrid(x, z)

        yy = (-normal_vector[0] * xx - normal_vector[2] * zz - d) / normal_vector[1]

        # Iterate over the grid and update the 3D image
        for i in range(z_dim - 1):
            for j in range(x_dim - 1):
                x, y, z = int(xx[i, j]), int(yy[i, j]), int(zz[i, j])
                if 0 <= x < x_dim and 0 <= y < y_dim and 0 <= z < z_dim:
                    modified_3d_plane[x, y, z] = modified_plane_2d_image[i, j]
    if usexy:
        x = np.linspace(0, x_dim - 1, x_dim)
        y = np.linspace(0, y_dim - 1, y_dim)

        xx, yy = np.meshgrid(x, y)

        zz = (-normal_vector[0] * xx - normal_vector[1] * yy - d) / normal_vector[2]

        for i in range(y_dim - 1):
            for j in range(x_dim - 1):
                x, y, z = int(xx[i, j]), int(yy[i, j]), int(zz[i, j])
                if 0 <= x < x_dim and 0 <= y < y_dim and 0 <= z < z_dim:
                    modified_3d_plane[x, y, z] = modified_plane_2d_image[i, j]

    if useyz:
        z = np.linspace(0, z_dim - 1, z_dim)
        y = np.linspace(0, y_dim - 1, y_dim)

        yy, zz = np.meshgrid(y, z)

        xx = (-normal_vector[1] - normal_vector[2] * zz * yy - d) / normal_vector[0]

        for i in range(z_dim - 1):
            for j in range(y_dim - 1):
                x, y, z = int(xx[i, j]), int(yy[i, j]), int(zz[i, j])
                if 0 <= x < x_dim and 0 <= y < y_dim and 0 <= z < z_dim:
                    modified_3d_plane[x, y, z] = modified_plane_2d_image[i, j]

    ni.save(
        ni.Nifti1Image(modified_3d_plane, meta.affine, meta.header), output_filename
    )


def distance_projection(normal_vector, centroid, xa, ya, za):
    d = -centroid.dot(normal_vector)
    a, b, c = normal_vector[0], normal_vector[1], normal_vector[2]
    return np.abs(a * xa + b * ya + c * za + d) / (np.linalg.norm(normal_vector))


def lambda_proj(normal_vector, centroid, xa, ya, za):
    d = -centroid.dot(normal_vector)
    a, b, c = normal_vector[0], normal_vector[1], normal_vector[2]
    return (a * xa + b * ya + c * za + d) / (np.linalg.norm(normal_vector) ** 2)


def projection(
    cutting_surface_filename: str | Path,
    output_filename: str | Path,
    centroid: Sequence[float],
    normal_vector: Sequence[float],
) -> None:
    meta = ni.load(cutting_surface_filename)
    arr = meta.get_fdata()

    centroid = np.asarray(centroid, dtype=float)
    normal_vector = np.asarray(normal_vector, dtype=float)
    size = arr.shape
    res = np.zeros(size)
    a, b, c = normal_vector[0], normal_vector[1], normal_vector[2]
    for xa in range(size[0]):
        for ya in range(size[1]):
            for za in range(size[2]):
                if (
                    arr[xa, ya, za] == 1
                    and distance_projection(normal_vector, centroid, xa, ya, za) < 30
                ):  # np.abs(value)<5 and
                    const = lambda_proj(normal_vector, centroid, xa, ya, za)
                    xh, yh, zh = (
                        int(xa - const * a),
                        int(ya - const * b),
                        int(za - const * c),
                    )

                    if xh > 0 and yh > 0 and zh > 0:
                        res[xh][yh][zh] = 1

    ni.save(ni.Nifti1Image(res, meta.affine, meta.header), output_filename)


def parse_command_line(argv):
    """Parse the script's command line."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Inter-bloc registration toolbox",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # sub parser skull extraction
    p1 = subparsers.add_parser(
        "extract_hull", help="Extract hull of a binary segmentation"
    )
    p1.add_argument(
        "binary_segmentation_filename", help="Filename of the binary segmentation"
    )
    p1.add_argument("output_filename", help="Ouput filename of the hull")
    p1.set_defaults(
        func=lambda args: extract_hull(
            args.binary_segmentation_filename, args.output_filename
        )
    )

    # sub parser plane_extractor_SVD
    p2 = subparsers.add_parser(
        "run_svd",
        help="Extraction of the centroid and the normal vector to the plane from the SVD",
    )
    p2.add_argument("filename", help="Name of the file on which to perform the SVD")
    p2.set_defaults(func=lambda args: plane_extractor_SVD(args.filename))

    # sub parser intensity extraction
    p3 = subparsers.add_parser(
        "extract_intensity", help="Extract intensity within the equation plan"
    )
    p3.add_argument("input_filename", help="Filename of the reconstructed block")
    p3.add_argument("output_filename", help="Filename of the 2D slice extracted")
    p3.add_argument(
        "--centroid",
        nargs=3,
        type=float,
        required=True,
        help="Centroid coordinates (x y z)",
    )
    p3.add_argument(
        "--normal",
        nargs=3,
        type=float,
        required=True,
        help='Normal vector (nx ny nz). WARNING: if you have negative values, you must enclose the value in quotation marks and leave a space at the beginning (e.g. " -2.8456e-09")',
    )
    p3.add_argument(
        "--usexz", action="store_true", help="Use the algorithm to extract xz plane"
    )
    p3.add_argument(
        "--usexy", action="store_true", help="Use the algorithm to extract xy plane"
    )
    p3.add_argument(
        "--useyz", action="store_true", help="Use the algorithm to extract yz plane"
    )
    p3.set_defaults(
        func=lambda args: get_2dplane(
            args.input_filename,
            args.output_filename,
            args.centroid,
            args.normal,
            args.usexz,
            args.usexy,
            args.useyz,
        )
    )

    # sub parser back projection
    p4 = subparsers.add_parser(
        "backproj", help="Project back the intensity in the 3D space"
    )
    p4.add_argument("input_filename", help="Filename of the reconstructed block")
    p4.add_argument("seg_filename", help="Filename of the manual segmentation")
    p4.add_argument("output_filename", help="Filename of the 2D slice extracted")
    p4.add_argument(
        "--centroid",
        nargs=3,
        type=float,
        required=True,
        help="Centroid coordinates (x y z)",
    )
    p4.add_argument(
        "--normal",
        nargs=3,
        type=float,
        required=True,
        help='Normal vector (nx ny nz). WARNING: if you have negative values, you must enclose the value in quotation marks and leave a space at the beginning (e.g. " -2.8456e-09")',
    )
    p4.add_argument(
        "--usexz", action="store_true", help="Use the algorithm to extract xz plane"
    )
    p4.add_argument(
        "--usexy", action="store_true", help="Use the algorithm to extract xy plane"
    )
    p4.add_argument(
        "--useyz", action="store_true", help="Use the algorithm to extract yz plane"
    )
    p4.set_defaults(
        func=lambda args: back_proj(
            args.input_filename,
            args.seg_filename,
            args.output_filename,
            args.centroid,
            args.normal,
            args.usexz,
            args.usexy,
            args.useyz,
        )
    )

    # sub parser intensity extraction
    p5 = subparsers.add_parser(
        "project",
        help="Project in the space of filename the binary cutting surface segmentation in the estimated plane",
    )
    p5.add_argument(
        "input_filename",
        help="Filename of the manual segmentation of the cutting surface",
    )
    p5.add_argument(
        "output_filename",
        help="Filename of the cutting surface projected in the estimated plane",
    )
    p5.add_argument(
        "--centroid",
        nargs=3,
        type=float,
        required=True,
        help="Centroid coordinates (x y z)",
    )
    p5.add_argument(
        "--normal",
        nargs=3,
        type=float,
        required=True,
        help='Normal vector (nx ny nz). WARNING: if you have negative values, you must enclose the value in quotation marks and leave a space at the beginning (e.g. " -2.8456e-09")',
    )
    p4.set_defaults(
        func=lambda args: projection(
            args.input_filename, args.output_filename, args.centroid, args.normal
        )
    )

    args = parser.parse_args()
    return args


def main(argv=sys.argv):
    """The script's entry point."""
    logging.basicConfig(level=logging.INFO)
    args = parse_command_line(argv)
    return args.func(args) or 0


if __name__ == "__main__":
    sys.exit(main())
