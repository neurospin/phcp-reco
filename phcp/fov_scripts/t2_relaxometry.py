import glob
import json
import logging
import os
import subprocess
import sys

import numpy as np
import nibabel as ni


logger = logging.getLogger(__name__)


def ConversionGisToNifti(InputGisFilename, OutputNiftiFilename):
    subprocess.run(
        [
            "singularity",
            "exec",
            "--bind",
            "/neurospin:/neurospin:rw",
            "/neurospin/phcp/code/gkg/2022-12-20_gkg/2022-12-20_gkg.sif",
            "GkgExecuteCommand",
            "Gis2NiftiConverter",
            "-i",
            InputGisFilename,
            "-o",
            OutputNiftiFilename,
            "-verbose",
        ]
    )


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


""" Single Compartment T2 Relaxometry Mapper """


def createCommand_SubVolume_gkgMethod(fileNameMSME, t2FileName):
    command = (
        "GkgExecuteCommand SubVolume"
        + " -i "
        + fileNameMSME
        + " -o "
        + t2FileName
        + " -tIndices "
        + str(0)
        + " -verbose true"
    )
    return command


def createCommand_Mask_gkgMethod(fileNameMSME, fileNameMask):
    command = (
        "GkgExecuteCommand GetMask"
        + " -i "
        + fileNameMSME
        + " -o "
        + fileNameMask
        + " -a 2 "
        + " -verbose "
    )
    return command


def createCommand_SingleCompartmentRelaxometryMapper(
    fileNameMSME,
    fileNameMask,
    fileNameTEValues,
    fileNameProtonDensity,
    fileNameT2,
    fileNameFittedMSME,
    verbose,
):
    # optimizerParameters :: <P1> : NLP maximum iteration count  <P2> : NLP test size <P3> : 0->do not apply MCMC 1->apply MCMC <P4> : MCMC burnin count <P5> : MCMC sample count <P6> : MCMC interval count <P7> : MCMC maximum iteration coun
    # ' -optimizerParameters 50000 0.001 0 2000 300 50 10000' +\  ' -scalarParameters 0.75 50 0 0 8 180 1 20 5'  +\
    command = (
        "GkgExecuteCommand SingleCompartmentRelaxometryMapper"
        + " -i %s" % fileNameMSME
        + " -m %s" % fileNameMask
        + " -t t2-mapping-msme"
        + " -optimizerParameters %s " % 50000
        + " -optimizerParameters %s " % 0.001
        + " -optimizerParameters %s" % 0
        + " -optimizerParameters %s" % 2000
        + " -optimizerParameters %s" % 300
        + " -optimizerParameters %s" % 50
        + " -optimizerParameters %s" % 10000
        + " -scalarParameters %s" % 0.75
        + " -scalarParameters %s" % 50
        + " -scalarParameters %s" % 0
        + " -scalarParameters %s" % 0
        + " -scalarParameters %s" % 8
        + " -scalarParameters %s" % 180
        + " -scalarParameters %s" % 1
        + " -scalarParameters %s" % 20
        + " -scalarParameters %s" % 5
        + " -stringParameters %s" % fileNameTEValues
        + " -op %s" % fileNameProtonDensity
        + " -ot %s" % fileNameT2
        + " -f %s" % fileNameFittedMSME
        + " -ascii False"
        + " -format gis"
        + " -verbose %s" % verbose
    )
    return command


def run_command(command):
    os.system(
        "singularity exec --bind /neurospin:/neurospin:rw "
        + " /neurospin/phcp/code/gkg/2022-12-20_gkg/2022-12-20_gkg.sif "
        + command
    )
    return None


def write_EchotTime(fileNameTEValues, MSMEFilenames):
    with open(fileNameTEValues, "w") as file:
        for metadata_filename in sorted(glob.iglob(MSMEFilenames)):
            with open(metadata_filename) as metadata_file:
                metadata = json.load(metadata_file)
            file.write(str(metadata["EchoTime"] * 1000) + "\n")
    return None


def create_ConfidenceMap(
    fileNameMSME, fileNameFittedMSME, fileNameMask, fileNameT2nifti
):
    arr_fittedmsme = create_ArrayFromGisFilename(fileNameFittedMSME)
    arr_fittedmsme = arrange_ArrayToFlipHeaderFormat(arr_fittedmsme)

    meta_t2msme = ni.load(fileNameMSME)
    arr_t2msme = meta_t2msme.get_fdata()

    meta_mask_t2 = ni.load(fileNameMask)
    arr_maskt2 = meta_mask_t2.get_fdata()
    arr_maskt2 = arrange_ArrayToFlipHeaderFormat(arr_maskt2)

    diff = arr_fittedmsme - arr_t2msme
    std_diff = np.std(diff[:, :, :, :10], axis=-1) * arr_maskt2

    meta_t2nifti = ni.load(fileNameT2nifti)
    ni_im = ni.Nifti1Image(
        dearrange_ArrayToFlipHeaderFormat(std_diff),
        meta_t2nifti.affine,
        meta_t2nifti.header,
    )
    return ni_im


def runT2RelaxometryMapper(
    subjectDirectoryGisConversion, MSMEFilenames, outputDirectory, verbose
):
    logger.info("SINGLE COMPARTMENT T2 RELAXOMETRY MAPPER")

    fileNameMSME = os.path.join(subjectDirectoryGisConversion, "t2-msme.nii.gz")
    fileNameMask = os.path.join(subjectDirectoryGisConversion, "mask.nii.gz")
    fileNameTEValues = os.path.join(outputDirectory, "te.txt")
    t2FileName = os.path.join(subjectDirectoryGisConversion, "t2_extracted.nii.gz")

    if not (os.path.exists(t2FileName)):
        command = createCommand_SubVolume_gkgMethod(fileNameMSME, t2FileName)
        run_command(command)

    if not (os.path.exists(fileNameMask)):
        command = createCommand_Mask_gkgMethod(fileNameMSME, fileNameMask)
        run_command(command)

    write_EchotTime(fileNameTEValues, MSMEFilenames)

    fileNameProtonDensity = os.path.join(outputDirectory, "proton-density.ima")
    fileNameFittedMSME = os.path.join(outputDirectory, "fitted-msme.ima")
    fileNameT2 = os.path.join(outputDirectory, "T2.ima")

    command = createCommand_SingleCompartmentRelaxometryMapper(
        fileNameMSME,
        fileNameMask,
        fileNameTEValues,
        fileNameProtonDensity,
        fileNameT2,
        fileNameFittedMSME,
        verbose,
    )
    run_command(command)

    fileNameT2nifti = os.path.join(outputDirectory, "T2.nii.gz")
    ConversionGisToNifti(fileNameT2, fileNameT2nifti)

    ni_im = create_ConfidenceMap(
        fileNameMSME, fileNameFittedMSME, fileNameMask, fileNameT2nifti
    )
    ni.save(ni_im, os.path.join(outputDirectory, "T2ConfidenceMap.nii.gz"))
    return None


def parse_command_line(argv):
    """Parse the script's command line."""
    import argparse

    parser = argparse.ArgumentParser(
        description="reconstruct T2 relaxometry based on the Multi Echo Spin Echo method",
    )
    parser.add_argument(
        "-i", "--input", dest="t2msmedirectory", help="T2-MSME Directory"
    )
    parser.add_argument(
        "-m",
        "--msme",
        dest="MSMEFilenames",
        help="String for glob search of MSME volumes. Ex: /phcp/rawdata/sub/ses/anat/sub_ses_echo*_MESE.json",
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
        runT2RelaxometryMapper(
            args.t2msmedirectory,
            args.MSMEFilenames,
            args.outputDirectory,
            args.verbose,
        )
        or 0
    )


if __name__ == "__main__":
    sys.exit(main())
