#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 10:32:37 2023

@author: la272118
"""

import glob
import optparse
import os
import shutil
import nibabel as ni
import numpy as np


def ValueExtractor(sesDirectory, outputDirectory):
    SessionOutputDir = os.path.join(outputDirectory, sesDirectory.split("/")[-2])
    if not os.path.exists(SessionOutputDir):
        os.mkdir(SessionOutputDir)

    modalitiesDirectory = glob.iglob(os.path.join(sesDirectory, "*"))
    session = sesDirectory.split("/")[-2]
    print("================================")
    print("  Début de la session " + session)
    print("================================")
    for modality in modalitiesDirectory:
        ModOutputDir = os.path.join(outputDirectory, session, modality.split("/")[-1])
        if not os.path.exists(ModOutputDir):
            os.mkdir(ModOutputDir)

        filenamesDirectory = glob.iglob(os.path.join(modality, "*.nii.gz"))
        print("================================")
        print("  Début de la modalité " + modality.split("/")[-1])
        print("================================")

        for volume in filenamesDirectory:
            print("================================")
            print("  Début de la conversion de " + volume.split("/")[-1])
            print("================================")
            if not (
                os.path.exists(
                    os.path.join(
                        outputDirectory,
                        session,
                        modality.split("/")[-1],
                        volume.split("/")[-1],
                    )
                )
            ):
                outputFilename = os.path.join(
                    outputDirectory,
                    session,
                    modality.split("/")[-1],
                    volume.split("/")[-1],
                )
                meta = ni.load(volume)
                matRotation = np.array(
                    [[1, 0, 0, 0], [0, 0, -1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]
                )
                ni_im = ni.Nifti1Image(
                    meta.get_fdata(), np.dot(matRotation, meta.affine), meta.header
                )
                # meta.affine=np.dot(matRotation, meta.affine) #Serait plus rapide mais ne change pas les axes dans anatomist
                # ni.save(meta, outputFilename)
                ni.save(ni_im, outputFilename)

            JsonFilename = os.path.join(
                sesDirectory,
                modality.split("/")[-1],
                (volume.split("/")[-1]).split(".")[0] + ".json",
            )

            if os.path.exists(JsonFilename) and not (
                os.path.exists(
                    os.path.join(
                        outputDirectory,
                        session,
                        modality.split("/")[-1],
                        (volume.split("/")[-1]).split(".")[0] + ".json",
                    )
                )
            ):
                outputJsonFilename = os.path.join(
                    outputDirectory,
                    session,
                    modality.split("/")[-1],
                    (volume.split("/")[-1]).split(".")[0] + ".json",
                )
                shutil.copy(JsonFilename, outputJsonFilename)


################################################################################
# parser to get option(s)
################################################################################

parser = optparse.OptionParser()
parser.add_option(
    "-i", "--subDirectory", dest="subDirectory", help="Subject file to flip directory"
)
parser.add_option(
    "-o", "--outputDirectory", dest="outputDirectory", help="Output directory"
)

(options, args) = parser.parse_args()


################################################################################
# 1) Value Extractor
################################################################################

ValueExtractor(options.subDirectory, options.outputDirectory)
