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

from phcp.gkg import (
    gkg_convert_gis_to_nifti,
    run_gkg_command,
    run_gkg_GetMask,
    run_gkg_SubVolume,
)


logger = logging.getLogger(__name__)


def create_ArrayFromGisFilename(GisFilename):
    DimFilename = GisFilename[:-3] + "dim"
    TxtInDimFile = open(DimFilename, "r").read()
    TxtFirstLine = TxtInDimFile.split("\n")[0]
    shape = tuple(np.int16(TxtFirstLine.split(" ")))  # [:-1]
    arr = np.fromfile(GisFilename, np.float32())
    arr = np.reshape(arr, shape, order="F")
    return arr


def arrange_ArrayToFlipHeaderFormat(arr):
    newarr = np.swapaxes(arr, 1, 2)
    newarr = np.flip(newarr, 1)
    newarr = np.flip(newarr, 0)
    return newarr


def dearrange_ArrayToFlipHeaderFormat(arr):
    newarr = np.flip(arr, 0)
    newarr = np.flip(newarr, 1)
    newarr = np.swapaxes(newarr, 1, 2)
    return newarr


def VFAregister(img_filename, output_filename, verbose):
    logger.info("========= Loading data =========")

    img_4D = ni.load(img_filename)
    arr_4D = img_4D.get_fdata()
    volume_number = arr_4D.shape[-1]
    ANTS_fixed_image = ants.from_numpy(arr_4D[:, :, :, 0])
    registered_volumes = []
    registered_volumes.append(ANTS_fixed_image)

    for i in range(1, volume_number):
        logger.info("========= Registration volume n°%d =========", i)
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
        "gis",
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


def create_ConfidenceMap(
    fileNameVFACat, fileNameFittedVFA, fileNameMask, fileNameT1nifti
):
    arr_fittedvfa = create_ArrayFromGisFilename(fileNameFittedVFA)
    arr_fittedvfa = arrange_ArrayToFlipHeaderFormat(arr_fittedvfa)

    meta_vfa = ni.load(fileNameVFACat)
    arr_vfa = meta_vfa.get_fdata()

    meta_mask = ni.load(fileNameMask)
    arr_mask = meta_mask.get_fdata()
    arr_mask = arrange_ArrayToFlipHeaderFormat(arr_mask)

    diff = arr_fittedvfa - arr_vfa
    std_diff = np.std(diff, axis=-1) * arr_mask

    metaT1nifti = ni.load(fileNameT1nifti)

    ni_im = ni.Nifti1Image(
        dearrange_ArrayToFlipHeaderFormat(std_diff),
        metaT1nifti.affine,
        metaT1nifti.header,
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


def correct_qt1(T1GisFilename, verbose=True):
    logger.info("Loading data")

    T1NiftiFilename = T1GisFilename.split(".")[0] + ".nii.gz"

    if not (os.path.exists(T1NiftiFilename)):
        gkg_convert_gis_to_nifti(T1GisFilename, T1NiftiFilename, verbose=True)

    ants_T1img = ants.image_read(T1NiftiFilename)
    ants_T1img_mask = ants_T1img.get_mask()

    logger.info("Bias Field Processing")

    biasField = extract_BiasField(ants_T1img, ants_T1img_mask, verbose)

    ants.image_write(biasField, T1GisFilename.split(".")[0] + "BiasField1.nii.gz")

    logger.info("New Bias Field Generation")

    newbiasField = create_NewBiasField(biasField)

    meta = ni.load(T1NiftiFilename)
    fwhm = proc.sigma2fwhm(6)
    ni_im = ni.Nifti1Image(np.float32(newbiasField), meta.affine, meta.header)

    logger.info("Smoothing part")

    newbiasField_smooth = proc.smooth_image(ni_im, fwhm)
    newbiasField_smooth_Filename = T1GisFilename.split(".")[0] + "_BiasField.nii.gz"
    ni.save(newbiasField_smooth, newbiasField_smooth_Filename)

    newbiasField_smooth = ants.from_nibabel(newbiasField_smooth)

    unbiased_T1 = ants_T1img * newbiasField_smooth

    logger.info("Saving Part")

    T1Nifti_unbiased_Filename = T1GisFilename.split(".")[0] + "_rec-unbiased.nii.gz"

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
    fileNameB1 = os.path.join(subjectDirectoryGisConversion, "b1_registred.nii.gz")
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

    fileNameT1nifti = os.path.join(outputDirectory, "T1.nii.gz")
    fileNameProtonDensity = os.path.join(outputDirectory, "proton-density.ima")
    fileNameFittedVFA = os.path.join(outputDirectory, "fitted-vfa.ima")
    fileNameOutputT1 = os.path.join(outputDirectory, "T1.ima")

    if not (os.path.exists(fileNameT1nifti)):
        write_TR_FA(fileNameFAValues, fileNameTRValues, FAFilenames, fileNameB1)

        command = createCommand_T1SingleCompartmentRelaxometryMapper(
            fileNameVFACat,
            fileNameMask,
            fileNameTRValues,
            fileNameFAValues,
            fileNameB1,
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

        gkg_convert_gis_to_nifti(fileNameOutputT1, fileNameT1nifti, verbose=True)

    fileNameConfidenceMap = os.path.join(outputDirectory, "T1ConfidenceMap.nii.gz")

    if not (os.path.exists(fileNameConfidenceMap)):
        ni_im = create_ConfidenceMap(
            fileNameVFACat, fileNameFittedVFA, fileNameMask, fileNameT1nifti
        )
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
    parser.add_argument("-i", "--input", dest="vfadirectory", help="VFA Directory")
    parser.add_argument(
        "-m",
        "--vfa",
        dest="VFAFilenames",
        help="String for glob search of VFA volumes. Ex: /phcp/rawdata/sub/ses/anat/sub_ses_flip*_VFA.json",
    )
    parser.add_argument(
        "-o", "--outputDirectory", dest="outputDirectory", help="Output directory"
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
            args.vfadirectory, args.VFAFilenames, args.outputDirectory, args.verbose
        )
        or 0
    )


if __name__ == "__main__":
    sys.exit(main())
