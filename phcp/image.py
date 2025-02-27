"""Utilities for working with images"""

import nibabel


def nibabel_orient_like(
    img: nibabel.spatialimages.SpatialImage, ref: nibabel.spatialimages.SpatialImage
) -> nibabel.spatialimages.SpatialImage:
    img_ornt = nibabel.orientations.io_orientation(img.affine)
    ref_ornt = nibabel.orientations.io_orientation(ref.affine)
    return img.as_reoriented(nibabel.orientations.ornt_transform(img_ornt, ref_ornt))
