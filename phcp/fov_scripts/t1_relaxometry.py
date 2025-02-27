"""Single Compartment T1 Relaxometry Mapper"""

import glob
import json
import logging
import math
import os
import sys

import numpy as np
import nibabel as ni
import nibabel.processing as proc
import ants
from sklearn.mixture import GaussianMixture

from phcp.fsl import run_fsl_command
from phcp.gkg import (
    run_gkg_command,
    run_gkg_GetMask,
    run_gkg_SubVolume,
)
from phcp.image import nibabel_orient_like


logger = logging.getLogger(__name__)


def VFAregister(img_filename, output_filename, verbose):
    logger.info("========= Loading data =========")

    img_4D = ni.load(img_filename)
    arr_4D = img_4D.get_fdata()
    volume_number = arr_4D.shape[-1]
    ANTS_fixed_image = ants.from_numpy(arr_4D[:, :, :, 0])
    registered_volumes = []
    registered_volumes.append(ANTS_fixed_image)

    for i in range(1, volume_number):
        logger.info(
            "========= Registration volume n°%d / %d =========", i + 1, volume_number
        )
        ANTS_moving_image = ants.from_numpy(arr_4D[:, :, :, i])
        registration = ants.registration(
            ANTS_fixed_image,
            ANTS_moving_image,
            type_of_transform="TRSAA",
            aff_iterations=(200, 200, 200),
            aff_shrink_factors=(8, 4, 2),
            aff_smoothing_sigmas=(3, 2, 1),
            verbose=False,
        )
        warped_image = registration["warpedmovout"]
        registered_volumes.append(warped_image)

    logger.info("========= Concatenating and writing the registered VFA =========")
    registered_data = np.stack([img.numpy() for img in registered_volumes], axis=-1)
    ni_im = ni.Nifti1Image(registered_data, img_4D.affine)
    ni.save(ni_im, output_filename)


def createCommand_T1SingleCompartmentRelaxometryMapper(
    fileNameVFACat,
    fileNameMask,
    fileNameTRValues,
    fileNameFAValues,
    fileNameB1,
    fileNameProtonDensity,
    fileNameOutputT1,
    fileNameFittedVFA,
    verbose,
) -> list[str]:
    command = [
        "GkgExecuteCommand",
        "SingleCompartmentRelaxometryMapper",
        "-i",
        fileNameVFACat,
        "-m",
        fileNameMask,
        "-t",
        "t1-mapping-vfa-spgr",
        "-optimizerParameters",
        "50000",
        "-optimizerParameters",
        "0.001",
        "-optimizerParameters",
        "1",
        "-optimizerParameters",
        "200",
        "-optimizerParameters",
        "100",
        "-optimizerParameters",
        "10",
        "-optimizerParameters",
        "1000",
        "-scalarParameters",
        "5",
        "-scalarParameters",
        "50",
        "-scalarParameters",
        "0",
        "-scalarParameters",
        "5",
        "-scalarParameters",
        "50",
        "-scalarParameters",
        "500",
        "-scalarParameters",
        "1",
        "-scalarParameters",
        "20",
        "-scalarParameters",
        "0.002",
        "-stringParameters",
        str(fileNameTRValues),
        "-stringParameters",
        str(fileNameFAValues),
        "-stringParameters",
        str(fileNameB1),
        "-op",
        fileNameProtonDensity,
        "-ot",
        fileNameOutputT1,
        "-f",
        fileNameFittedVFA,
        "-ascii",
        "False",
        "-format",
        "nifti",
        "-verbose",
        str(verbose),
    ]
    return command


def write_TR_FA(fileNameFAValues, fileNameTRValues, FAFilenames, fileNameB1):
    with open(fileNameFAValues, "w") as file:
        with open(fileNameTRValues, "w") as file2:
            for metadata_filename in sorted(glob.iglob(FAFilenames)):
                with open(metadata_filename) as metadata_file:
                    metadata = json.load(metadata_file)
                if os.path.exists(fileNameB1):
                    file.write(str(metadata["FlipAngle"]) + "\n")
                else:
                    file.write(
                        str(math.radians(metadata["FlipAngle"])) + "\n"
                    )  # Si pas de carte  de B1 --> convertir les angles en radians
                file2.write(str(metadata["RepetitionTime"] * 1000) + "\n")

    return None


def create_ConfidenceMap(fileNameVFACat, fileNameFittedVFA, fileNameMask):
    meta_fittedvfa = ni.load(fileNameFittedVFA)
    arr_fittedvfa = meta_fittedvfa.get_fdata()

    meta_vfa = nibabel_orient_like(ni.load(fileNameVFACat), meta_fittedvfa)
    arr_vfa = meta_vfa.get_fdata()

    meta_mask = nibabel_orient_like(ni.load(fileNameMask), meta_fittedvfa)
    arr_mask = meta_mask.get_fdata()

    diff = arr_fittedvfa - arr_vfa
    std_diff = np.std(diff, axis=-1) * arr_mask

    ni_im = ni.Nifti1Image(
        std_diff,
        meta_fittedvfa.affine,
        meta_fittedvfa.header,
    )
    return ni_im


""" Correction Part """


def extract_BiasField(ants_T1img, ants_T1img_mask, verbose):
    return ants.n4_bias_field_correction(
        ants_T1img,
        ants_T1img_mask,
        spline_param=20,
        rescale_intensities=True,
        return_bias_field=True,
        verbose=verbose,
    )


def create_NewBiasField(biasField):
    flat_arr = biasField[biasField > 0].reshape(-1, 1)
    gmm = GaussianMixture(n_components=4)
    gmm.fit(flat_arr)

    means = gmm.means_.flatten().tolist()
    logger.debug("means: %s", means)
    moy_inconsistency, index_inconsistency = min(means), means.index(min(means))
    means.remove(moy_inconsistency)
    moy_normal = min(means)  # index_normal = means.index(min(means))

    stds = np.sqrt(gmm.covariances_.flatten()).tolist()
    std_inconsistency = stds[index_inconsistency]

    biasField_numpy = biasField.numpy()
    return np.where(
        biasField_numpy < (moy_inconsistency + 2 * std_inconsistency),
        moy_normal / biasField_numpy,
        1,
    )


def correct_qt1(T1NiftiFilename, verbose=True):
    logger.info("Loading data")

    T1_basename = T1NiftiFilename.split(".")[0]

    ants_T1img = ants.image_read(T1NiftiFilename)
    ants_T1img_mask = ants_T1img.get_mask()

    logger.info("Bias Field Processing")

    biasField = extract_BiasField(ants_T1img, ants_T1img_mask, verbose)

    ants.image_write(biasField, T1_basename + "BiasField1.nii.gz")

    logger.info("New Bias Field Generation")

    newbiasField = create_NewBiasField(biasField)

    meta = ni.load(T1NiftiFilename)
    fwhm = proc.sigma2fwhm(6)
    ni_im = ni.Nifti1Image(np.float32(newbiasField), meta.affine, meta.header)

    logger.info("Smoothing part")

    newbiasField_smooth = proc.smooth_image(ni_im, fwhm)
    newbiasField_smooth_Filename = T1_basename + "_BiasField.nii.gz"
    ni.save(newbiasField_smooth, newbiasField_smooth_Filename)

    newbiasField_smooth = ants.from_nibabel(newbiasField_smooth)

    unbiased_T1 = ants_T1img * newbiasField_smooth

    logger.info("Saving Part")

    T1Nifti_unbiased_Filename = T1_basename + "_rec-unbiased.nii.gz"

    ants.image_write(unbiased_T1, T1Nifti_unbiased_Filename)
    return None


def runT1RelaxometryMapper(
    subjectDirectoryGisConversion, FAFilenames, outputDirectory, verbose
):
    fileNameVFACat0 = os.path.join(subjectDirectoryGisConversion, "t1map.nii.gz")
    fileNameVFACat = os.path.join(subjectDirectoryGisConversion, "t1map-TRSAA.nii.gz")

    if not (os.path.exists(fileNameVFACat)):
        logger.info("T1 RELAXOMETRY : VFA REGISTRATION")
        VFAregister(fileNameVFACat0, fileNameVFACat, verbose)

    logger.info("SINGLE COMPARTMENT T1 RELAXOMETRY MAPPER")

    fileNameMask = os.path.join(subjectDirectoryGisConversion, "mask.nii.gz")
    t1FileName = os.path.join(subjectDirectoryGisConversion, "t1_extracted.nii.gz")
    fileNameB1 = os.path.join(subjectDirectoryGisConversion, "b1.nii.gz")
    fileNameB1inVFAspace = os.path.join(
        subjectDirectoryGisConversion, "b1_in_vfa_space.nii.gz"
    )
    fileNameTRValues = os.path.join(outputDirectory, "tr.txt")
    fileNameFAValues = os.path.join(outputDirectory, "fa.txt")

    if not (os.path.exists(t1FileName)):
        run_gkg_SubVolume(
            [
                "-i",
                fileNameVFACat,
                "-o",
                t1FileName,
                "-tIndices",
                "0",
                "-verbose",
                "true",
            ],
            output_dirs=[subjectDirectoryGisConversion],
        )

    if not (os.path.exists(fileNameMask)):
        run_gkg_GetMask(
            ["-i", t1FileName, "-o", fileNameMask, "-a", "2", "-verbose"],
            output_dirs=[subjectDirectoryGisConversion],
        )

    if not os.path.exists(fileNameB1inVFAspace):
        run_fsl_command(
            [
                "flirt",
                "-applyxfm",
                "-usesqform",
                "-interp",
                "trilinear",
                "-in",
                fileNameB1,
                "-ref",
                fileNameVFACat,
                "-out",
                fileNameB1inVFAspace,
            ]
        )

    fileNameProtonDensity = os.path.join(outputDirectory, "proton-density.nii.gz")
    fileNameFittedVFA = os.path.join(outputDirectory, "fitted-vfa.nii.gz")
    fileNameOutputT1 = os.path.join(outputDirectory, "T1.nii.gz")

    if not (os.path.exists(fileNameOutputT1)):
        write_TR_FA(
            fileNameFAValues, fileNameTRValues, FAFilenames, fileNameB1inVFAspace
        )

        command = createCommand_T1SingleCompartmentRelaxometryMapper(
            fileNameVFACat,
            fileNameMask,
            fileNameTRValues,
            fileNameFAValues,
            fileNameB1inVFAspace,
            fileNameProtonDensity,
            fileNameOutputT1,
            fileNameFittedVFA,
            verbose,
        )
        run_gkg_command(
            command,
            gkg_container_version="2022-12-20",
            input_dirs=[subjectDirectoryGisConversion],
            output_dirs=[outputDirectory],
        )
        # FIXME: outputs have incorrect qform (no translation despite seemingly
        # correct rotation) and grossly incorrect sform (no scale, no translation)

    fileNameConfidenceMap = os.path.join(outputDirectory, "T1ConfidenceMap.nii.gz")

    if not (os.path.exists(fileNameConfidenceMap)):
        ni_im = create_ConfidenceMap(fileNameVFACat, fileNameFittedVFA, fileNameMask)
        ni.save(ni_im, fileNameConfidenceMap)

    logger.info("Correction Part")

    correct_qt1(fileNameOutputT1, verbose)
    return None


def parse_command_line(argv):
    """Parse the script's command line."""
    import argparse

    parser = argparse.ArgumentParser(
        description="reconstruct T1 relaxometry based on the Variable Flip Angle method",
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Directory with the inputs t1map.nii.gz and b1.nii.gz, "
        + "typically fov/derivatives/T1mapping/sub-${sub}/ses-${ses}/01-Materials",
    )
    parser.add_argument(
        "-m",
        "--vfa",
        dest="VFAFilenames",
        help="String for glob search of VFA volume metadata, typically fov/rawdata/sub-${sub}/ses-${ses}/anat/'*_VFA.json'",
    )
    parser.add_argument(
        "-o",
        "--outputDirectory",
        dest="outputDirectory",
        help="Output directory, typically fov/derivatives/T1mapping/sub-${sub}/ses-${ses}/02-Results",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=True,
        help="print detailed information during the fit",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        dest="verbose",
        action="store_false",
        help="do not print detailed information during the fit",
    )

    args = parser.parse_args()
    return args


def main(argv=sys.argv):
    """The script's entry point."""
    logging.basicConfig(level=logging.INFO)
    args = parse_command_line(argv)
    return (
        runT1RelaxometryMapper(
            args.input, args.VFAFilenames, args.outputDirectory, args.verbose
        )
        or 0
    )


if __name__ == "__main__":
    sys.exit(main())
