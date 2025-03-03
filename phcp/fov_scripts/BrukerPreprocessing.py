"""
Script :
    PreprocessingPipeline.py
Description :
    Based on JSON files :
        --> Write in rawdata directory .bvec and .bval files for each run of DWI
        --> Merge .bvec and .bval files and write them in Derivatives/Gkg-Pipeline/sub/ses/b___/01-Materials/ folder
        --> Correct DWI data for Receiver Gain (RG).
        --> Write DWI data corrected for RG in Derivatives/Gkg-Pipeline/sub/ses/b___/01-Materials/

        --> Correct FAs for RGs
        --> Merge FAs and create VFA file
        --> Write VFA file in Derivatives/RGCorrection/sub/ses/
Needs :
    FILE ARCHITECTURE
        ->../sourcedata/
            -sub/
                -ses/
        ->../rawdata/
            #BIDS norm
            -sub/
                -ses/
                    -anat/
                    -dwi/
                    -fmap/
        ->../derivatives/
            -gkg-Pipeline/
                -sub/
                    -ses/
                        -b___/
                            -01-Materials/
            -RGCorrection/
                -sub/
                    -ses/
                        #These folders will store corrected volumes for RGs (if needed)
                        -anat/
                        -dwi/
        ->../PreprocessingJSONdirectory/
            #This file describes the modalities that will be processed.
            -RG_description.json            (ex: {"sub_-ses_":["VFA","b4500"]})
            #This file
            -sub-_ses-_description.json     (ex: {"b4500": [102, 103,..]})
Usage :
    python code/code_gitHub/gkg-preprocessing/PreprocessingPipeline.py
        -s sourcedata/
        -r headerfliprawdata/
        -d derivatives/
        -q '/../../Subject.json'
        -p /../../PreprocessingJSONdirectory/
Authors :
    Lucas Arcamone
"""

import glob
import json
import logging
import nibabel
import os
import re
import sys

import numpy

from phcp import BrukerInformationParser


logger = logging.getLogger(__name__)


def print_message(message):
    print("==================================")
    print(message)
    print("==================================")


# Creates bvec and bvac files in 'rawdata' storage space
# from files output from MRI console (DICOM)
def runBvecAndBvalWriter(
    SourceDataDirectory, RawDataDirecctory, Subject, Session, DescriptionFilename
):
    with open(DescriptionFilename, "r") as p:
        DescriptionFile = json.load(p)

    for bvalue in DescriptionFile.keys():
        for compteur in range(len(DescriptionFile[bvalue])):
            # Dicom data are stored with a BIDS-like names (sub / ses)
            DiffDirectory = os.path.join(
                SourceDataDirectory,
                Subject,
                Session,
                str(DescriptionFile[bvalue][compteur]),
            )
            DiffMethod = os.path.join(DiffDirectory, "method")
            DiffAcqp = os.path.join(DiffDirectory, "acqp")
            DiffVisu = os.path.join(DiffDirectory, "visu_pars")

            Parser = BrukerInformationParser.BrukerInformationParser(
                DiffMethod, DiffAcqp, DiffVisu, False
            )

            RawdataDirectorySaver = os.path.join(
                RawDataDirecctory, Subject, Session, "dwi"
            )
            FileName = (
                Subject
                + "_"
                + Session
                + "_acq-"
                + bvalue
                + "_run-"
                + str(compteur + 1)
                + "_dwi"
            )

            Parser.saveBvecAndBval(
                os.path.join(RawdataDirectorySaver, FileName + ".bvec"),
                os.path.join(RawdataDirectorySaver, FileName + ".bval"),
            )
    return None


def runBvalMerger(DwiDirectory, Subject, Session, bvalue, OutputDirectory):
    prefix_filename = Subject + "_" + Session + "_acq-" + bvalue
    BvalFiles = sorted(
        glob.glob(os.path.join(DwiDirectory, prefix_filename + "*.bval"))
    )
    values = ""
    for bvalfile in BvalFiles:
        with open(bvalfile, "r") as f:
            values += f.readlines()[0][:-1]

    path = os.path.join(OutputDirectory, prefix_filename + "_dwi.bval")
    with open(path, "w") as p:
        p.write(values)
    return None


def runBvecMerger(DwiDirectory, Subject, Session, bvalue, OutputDirectory):
    prefix_filename = Subject + "_" + Session + "_acq-" + bvalue

    BvecFiles = sorted(
        glob.glob(os.path.join(DwiDirectory, prefix_filename + "*.bvec"))
    )
    valuex = ""
    valuey = ""
    valuez = ""
    for bvecfile in BvecFiles:
        with open(bvecfile, "r") as f:
            value = f.readlines()
            valuex += value[0][:-1]
            valuey += value[1][:-1]
            valuez += value[2][:-1]

    path = os.path.join(OutputDirectory, prefix_filename + "_dwi.bvec")
    result = valuex + "\n" + valuey + "\n" + valuez
    with open(path, "w") as p:
        p.write(result)
    return None


def is_bvalue(keyname):
    return ("b8000" in keyname) or ("b4500" in keyname) or ("b1500" in keyname)


def read_rg(acqp_filename):
    with open(acqp_filename) as f:
        for line in f:
            match = re.match(r"##\$RG=(.*)$", line)
            if match:
                rg = float(match.group(1))
                break
    return rg


def runRGWriter(
    RawDataDirectory,
    SourceDataDirectory,
    SaveDirectory,
    Subject,
    Session,
    SubSesKeyModalityValue_json,
):
    with open(SubSesKeyModalityValue_json, "r") as p:
        SubSesKeyModalityValue_dict = json.load(p)

    DictionnaireRg = {}

    if not os.path.exists(os.path.join(SaveDirectory, Subject)):
        os.mkdir(os.path.join(SaveDirectory, Subject))

    if not os.path.exists(os.path.join(SaveDirectory, Subject, Session)):
        os.mkdir(os.path.join(SaveDirectory, Subject, Session))

    key = Subject + "_" + Session

    for modality in SubSesKeyModalityValue_dict[key]:
        if modality == "VFA":
            RawDataSubjectDirectory = os.path.join(
                RawDataDirectory, Subject, Session, "anat"
            )

            FileName = Subject + "_" + Session + "_flip-*" + "_" + modality + ".json"

        else:  # Other value available for 'modality' variable are : "b1500", "b4500", "b8000"
            RawDataSubjectDirectory = os.path.join(
                RawDataDirectory, Subject, Session, "dwi"
            )

            FileName = (
                Subject + "_" + Session + "_acq-" + modality + "_run-*" + "_dwi.json"
            )

        for File in glob.iglob(os.path.join(RawDataSubjectDirectory, FileName)):
            with open(File, "r") as r:
                JsonFile = json.load(r)

            number = str(JsonFile["SeriesNumber"])[:-4]

            DirectoryFile = os.path.join(
                SourceDataDirectory, Subject, Session, number, "acqp"
            )

            DictionnaireRg[File.split("/")[-1][:-4] + "nii.gz"] = read_rg(DirectoryFile)

        jsonObject = json.dumps(DictionnaireRg)

        SaveFile = os.path.join(SaveDirectory, Subject, Session)

        with open(SaveFile + "/RG_values.json", "w") as f:
            f.write(jsonObject)
    return None


def correct_RG(
    RGCorrectionDirectory,
    RGDictionary,
    RawdataDirectory,
    Subject,
    Session,
    FolderName,
    GlobSearch,
):
    Files_iterator = sorted(
        glob.iglob(
            os.path.join(RawdataDirectory, Subject, Session, FolderName, GlobSearch)
        )
    )

    path_str = os.path.join(RGCorrectionDirectory, Subject, Session, FolderName)

    if not os.path.exists(path_str):
        os.mkdir(path_str)

    ref_filename = Files_iterator[0]
    ref = RGDictionary[ref_filename.split("/")[-1]]
    ni_volume_ref = nibabel.load(ref_filename)
    nibabel.save(
        ni_volume_ref,
        os.path.join(
            RGCorrectionDirectory,
            Subject,
            Session,
            FolderName,
            ref_filename.split("/")[-1],
        ),
    )
    Files_iterator.remove(ref_filename)

    for filename in Files_iterator:
        RGValue = RGDictionary[filename.split("/")[-1]]
        volume = nibabel.load(
            os.path.join(
                RawdataDirectory, Subject, Session, FolderName, filename.split("/")[-1]
            )
        )
        volume_arr = volume.get_fdata()
        volume_arr_corrected = volume_arr * ref / RGValue
        volume_nii = nibabel.Nifti1Image(
            volume_arr_corrected, volume.affine, volume.header
        )

        nibabel.save(
            volume_nii,
            os.path.join(
                RGCorrectionDirectory,
                Subject,
                Session,
                FolderName,
                filename.split("/")[-1],
            ),
        )
    return None


def merge_volume(prefix, inputFilenames, outputFilename, isVFA):
    Iterator = sorted(glob.iglob(inputFilenames))

    meta_FirstVolume = nibabel.load(Iterator[0])
    arr_FirstVolume = meta_FirstVolume.get_fdata()

    if isVFA:
        Dimension = (arr_FirstVolume.shape) + (1,)
        res = numpy.reshape([arr_FirstVolume], Dimension)
    else:
        res = arr_FirstVolume

    Iterator.remove(Iterator[0])

    for volume in Iterator:
        meta = nibabel.load(volume)
        arr = meta.get_fdata()
        if isVFA:
            res = numpy.concatenate(
                (res, numpy.reshape([arr], meta.shape + (1,))), axis=-1
            )
        else:
            res = numpy.concatenate((res, arr), axis=-1)

    ni_im = nibabel.Nifti1Image(res, meta_FirstVolume.affine, meta_FirstVolume.header)
    nibabel.save(ni_im, outputFilename)
    return None


def NiiMerger(FilesDirectory, DerivativesDirectory, modality, Subject, Session):
    SessionDirectory = os.path.join(FilesDirectory, Subject, Session)
    SaveDwiDirectory = os.path.join(
        DerivativesDirectory, "gkg-Pipeline", Subject, Session
    )

    isVFA = False
    if "VFA" == modality:
        print_message("VFA Merger")
        SaveVFADirectory = os.path.join(
            DerivativesDirectory, "RGCorrection", Subject, Session
        )

        prefix = Subject + "_" + Session + "_" + "flip-*" + "_VFA.nii.gz"
        inputFilenames = os.path.join(SessionDirectory, "anat", prefix)
        outputFilename = os.path.join(
            SaveVFADirectory, Subject + "_" + Session + "_VFA.nii.gz"
        )
        isVFA = True

    elif "b8000" == modality:
        print_message("b8000 Merger")

        prefix = Subject + "_" + Session + "_" + "acq-b8000_run-*" + "_dwi.nii.gz"
        inputFilenames = os.path.join(SessionDirectory, "dwi", prefix)
        outputFilename = os.path.join(
            SaveDwiDirectory,
            "b8000",
            "01-Materials",
            Subject + "_" + Session + "_acq-b8000_dwi.nii.gz",
        )

    elif "b4500" == modality:
        print_message("b4500 Merger")

        prefix = Subject + "_" + Session + "_" + "acq-b4500_run-*" + "_dwi.nii.gz"
        inputFilenames = os.path.join(SessionDirectory, "dwi", prefix)
        outputFilename = os.path.join(
            SaveDwiDirectory,
            "b4500",
            "01-Materials",
            Subject + "_" + Session + "_acq-b4500_dwi.nii.gz",
        )
    elif "b1500" == modality:
        print_message("b1500 Merger")

        prefix = Subject + "_" + Session + "_" + "acq-b1500_run-*" + "_dwi.nii.gz"
        inputFilenames = os.path.join(SessionDirectory, "dwi", prefix)
        outputFilename = os.path.join(
            SaveDwiDirectory,
            "b1500",
            "01-Materials",
            Subject + "_" + Session + "_acq-b1500_dwi.nii.gz",
        )

    merge_volume(prefix, inputFilenames, outputFilename, isVFA)
    # FIXME: delete non-concatenated volumes, because they are never used again
    return None


def runVolumesRGCorrection(
    RGValuesJsonFilename,
    SubSesKeyModalityValue_json,
    Subject,
    Session,
    RawdataDirectory,
    DerivativesDirectory,
    RGCorrectionDirectory,
):
    with open(RGValuesJsonFilename, "r") as f:
        RGDictionary = json.load(f)

    with open(SubSesKeyModalityValue_json, "r") as p:
        SubSesKeyModalityValue_dict = json.load(p)

    Id = Subject + "_" + Session

    for modality in SubSesKeyModalityValue_dict[Id]:
        if modality[0] == "b":  # if key is a bvalue
            prefix = Subject + "_" + Session + "_acq-" + modality
            filtered_RGDictionary = {
                k: v for k, v in RGDictionary.items() if k.startswith(prefix)
            }
            if len(set(filtered_RGDictionary.values())) > 1:
                FolderName = "dwi"
                GlobSearch = (
                    Subject
                    + "_"
                    + Session
                    + "_acq-"
                    + modality
                    + "_run-*"
                    + "_dwi.nii.gz"
                )
                FilesDirectory = RGCorrectionDirectory
                correct_RG(
                    RGCorrectionDirectory,
                    RGDictionary,
                    RawdataDirectory,
                    Subject,
                    Session,
                    FolderName,
                    GlobSearch,
                )
            else:
                FilesDirectory = RawdataDirectory
                print("All RG values are the same for the " + modality + " modality")

        elif modality == "VFA":
            prefix = Subject + "_" + Session + "_flip-"
            filtered_RGDictionary = {
                k: v for k, v in RGDictionary.items() if k.startswith(prefix)
            }

            if len(set(filtered_RGDictionary.values())) > 1:
                FolderName = "anat"
                GlobSearch = (
                    Subject + "_" + Session + "_flip-*" + "_" + modality + ".nii.gz"
                )
                FilesDirectory = RGCorrectionDirectory

                correct_RG(
                    RGCorrectionDirectory,
                    RGDictionary,
                    RawdataDirectory,
                    Subject,
                    Session,
                    FolderName,
                    GlobSearch,
                )
            else:
                FilesDirectory = RawdataDirectory
                print("All RG values are the same for the " + modality + " modality")
        NiiMerger(FilesDirectory, DerivativesDirectory, modality, Subject, Session)
    return None


def bruker_preprocessing(
    SourcedataDirectory,
    RawdataDirectory,
    DerivativesDirectory,
    SubjectKeySessionValue_json,
    DescriptionsDirectory,
):
    SubSesKeyModalityValue_json = os.path.join(
        DescriptionsDirectory, "RG_description.json"
    )

    with open(SubjectKeySessionValue_json, "r") as f:
        SubjectKeySessionValue_dict = json.load(f)

    with open(SubSesKeyModalityValue_json, "r") as q:
        SubSesKeyModalityValue_dict = json.load(q)

    for subject in SubjectKeySessionValue_dict.keys():
        for session in SubjectKeySessionValue_dict[subject]:
            print_message(subject + " / " + session)

            Id = subject + "_" + session
            if is_bvalue(SubSesKeyModalityValue_dict[Id]):
                BvalueKeyNumberFileValue_json = os.path.join(
                    DescriptionsDirectory, subject + "_" + session + "_description.json"
                )

                with open(BvalueKeyNumberFileValue_json, "r") as p:
                    BvalueKeyNumberFileValue_dict = json.load(p)

                print_message("Run BvecAndBvalWriter")

                runBvecAndBvalWriter(
                    SourcedataDirectory,
                    RawdataDirectory,
                    subject,
                    session,
                    BvalueKeyNumberFileValue_json,
                )

                for bvalue in BvalueKeyNumberFileValue_dict.keys():
                    SessionDirectory = os.path.join(RawdataDirectory, subject, session)
                    DwiDirectory = os.path.join(SessionDirectory, "dwi")

                    if not os.path.exists(
                        os.path.join(
                            DerivativesDirectory,
                            "gkg-Pipeline",
                            subject,
                            session,
                            bvalue,
                            "01-Materials",
                        )
                    ):
                        os.makedirs(
                            os.path.join(
                                DerivativesDirectory,
                                "gkg-Pipeline",
                                subject,
                                session,
                                bvalue,
                                "01-Materials",
                            )
                        )

                    MaterialsDirectory = os.path.join(
                        DerivativesDirectory,
                        "gkg-Pipeline",
                        subject,
                        session,
                        bvalue,
                        "01-Materials",
                    )

                    print_message("Run BvecAndBvalMerger for " + bvalue)

                    runBvalMerger(
                        DwiDirectory, subject, session, bvalue, MaterialsDirectory
                    )

                    runBvecMerger(
                        DwiDirectory, subject, session, bvalue, MaterialsDirectory
                    )

            RGCorrectionDirectory = os.path.join(DerivativesDirectory, "RGCorrection")
            os.makedirs(RGCorrectionDirectory, exist_ok=True)

            print_message("Run RGWriter")

            runRGWriter(
                RawdataDirectory,
                SourcedataDirectory,
                RGCorrectionDirectory,
                subject,
                session,
                SubSesKeyModalityValue_json,
            )

            RGSessionDirectory = os.path.join(RGCorrectionDirectory, subject, session)
            RGValuesJsonFilename = os.path.join(RGSessionDirectory, "RG_values.json")

            print_message("Run VolumesRGCorrection")

            runVolumesRGCorrection(
                RGValuesJsonFilename,
                SubSesKeyModalityValue_json,
                subject,
                session,
                RawdataDirectory,
                DerivativesDirectory,
                RGCorrectionDirectory,
            )

    return None


def parse_command_line(argv):
    """Parse the script's command line."""
    import argparse

    parser = argparse.ArgumentParser(
        "Early pre-processings that require access to the Bruker raw data: extraction "
        "of b-values and b-vectors, and receivergain correction. Refer to the "
        "repository's README for detailed usage instructions."
    )
    parser.add_argument(
        "-s",
        "--sourcedata",
        dest="sourcedata",
        help="Sourcedata directory, e.g. fov/sourcedata",
    )
    parser.add_argument(
        "-r",
        "--rawdata",
        dest="rawdata",
        help="Rawdata directory, e.g. fov/headerfliprawdata",
    )
    parser.add_argument(
        "-d",
        "--derivatives",
        dest="derivatives",
        help="Derivatives directory, e.g. fov/derivatives",
    )
    parser.add_argument(
        "-q",
        "--subject",
        dest="subject",
        help="Subject JSON filename, e.g. derivatives/gkg-Pipeline/PreprocessingDescriptions/sub-${sub}.json",
    )
    parser.add_argument(
        "-p",
        "--descriptions",
        dest="descriptions",
        help="Descriptions directory, e.g. derivatives/gkg-Pipeline/PreprocessingDescriptions/sub-${sub}",
    )

    args = parser.parse_args()
    return args


def main(argv=sys.argv):
    """The script's entry point."""
    logging.basicConfig(level=logging.INFO)
    args = parse_command_line(argv)
    return (
        bruker_preprocessing(
            args.sourcedata,
            args.rawdata,
            args.derivatives,
            args.subject,
            args.descriptions,
        )
        or 0
    )


if __name__ == "__main__":
    sys.exit(main())
