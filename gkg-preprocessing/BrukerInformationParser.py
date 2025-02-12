"""
Script :
    BrukerInformationParser.py
Description :
    Extract MRI sequence information such as receiver gain, b-vectors and b-values
    from acqp, visu and method files. These files come from the source data of the MRI console.
Needs :

Usage :
    Used in PreprocessingPipeline.py
Authors :
    Yann Leprince
"""

import math


class BrukerInformationParser:
    def __init__(
        self,
        methodFileName,
        acqpFileName,
        visuParsFileName,
        applyVisuCoreTransformation,
    ):
        f = open(methodFileName, "r")
        methodLines = f.readlines()
        f.close()

        f = open(acqpFileName, "r")
        acqpLines = f.readlines()
        f.close()

        f = open(visuParsFileName, "r")
        visuParsLines = f.readlines()
        f.close()

        self._acqGradMatrix = self.getMatrix(acqpLines, "ACQ_grad_matrix")
        self._acqPatientPos = self.getString(acqpLines, "ACQ_patient_pos")

        self._sliceOrientation = self.getString(methodLines, "PVM_SPackArrSliceOrient")
        self._readOrientation = self.getString(methodLines, "PVM_SPackArrReadOrient")

        # initializing the rotation matrix to identity
        self._rotation = [
            [float(1), float(0), float(0)],
            [float(0), float(1), float(0)],
            [float(0), float(0), float(1)],
        ]
        # and customizing it according to slice and read orientations
        if (self._sliceOrientation == "sagittal") and (self._readOrientation == "H_F"):
            self._rotation = [
                [float(0), float(0), float(1)],
                [float(0), float(-1), float(0)],
                [float(1), float(0), float(0)],
            ]

        elif (self._sliceOrientation == "coronal") and (self._readOrientation == "H_F"):
            self._rotation = [
                [float(0), float(1), float(0)],
                [float(0), float(0), float(1)],
                [float(1), float(0), float(0)],
            ]

        elif (self._sliceOrientation == "axial") and (self._readOrientation == "L_R"):
            self._rotation = [
                [float(0), float(1), float(0)],
                [float(0), float(0), float(1)],
                [float(1), float(0), float(0)],
            ]

        positionMatrices = {
            "Head_Supine": [[-1.0, 0.0, 0.0], [0.0, +1.0, 0.0], [0.0, 0.0, -1.0]],
            "Head_Prone": [[+1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
            "Feet_Supine": [[+1.0, 0.0, 0.0], [0.0, +1.0, 0.0], [0.0, 0.0, +1.0]],
            "Feet_Prone": [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, +1.0]],
        }
        self._positionMatrix = positionMatrices[self._acqPatientPos]

        self._paravisionVersion = self.getStringVector(
            visuParsLines, "VisuCreatorVersion"
        )[0]
        self._isParavision6 = False
        if self._paravisionVersion[1] == "6":
            self._isParavision6 = True

        self._visuCoreCoordinateSystems = ["L->R", "D->V", "F->H"]
        if self._isParavision6:
            self._visuCoreCoordinateSystems = ["R->L", "V->D", "F->H"]

        self._visuCoreOrientationMatrix = self.getMatrix(
            visuParsLines, "VisuCoreOrientation"
        )
        self._dwEffBValues = self.getFloatVector(methodLines, "PVM_DwEffBval")
        self._dwDir = self.getMatrix(methodLines, "PVM_DwDir")

        self._rdhToLdhTransformation3d = [
            [-1.0, 0.0, 0.0],
            [0.0, +1.0, 0.0],
            [0.0, 0.0, +1.0],
        ]

        self._ldhToRvhTransformation3d = [
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, +1.0],
        ]

        self._niftiOrientationMatrix = self.computeMatrixMultiplication(
            self._ldhToRvhTransformation3d, self._visuCoreOrientationMatrix
        )
        self._niftiOrientationMatrixDeterminant = self.computeDeterminant(
            self._niftiOrientationMatrix
        )

        self._b0Count = self.getInteger(methodLines, "PVM_DwAoImages")

        self._diffusionGradientOrientations = []

        for b0 in range(self._b0Count):
            self._diffusionGradientOrientations.append(
                [0.5773502691896257, 0.5773502691896257, 0.5773502691896257]
            )

        for o in self._dwDir:
            orientation = o
            if applyVisuCoreTransformation:
                # convert from physical frame to the subject coordinate system
                orientation = self.computeMatrixDotVector(
                    self._visuCoreOrientationMatrix, orientation
                )
                # convert from RDH frame to LDH frame
                orientation = self.computeMatrixDotVector(
                    self._rdhToLdhTransformation3d, orientation
                )
                # convert from LDH frame to RVH frame
                orientation = self.computeMatrixDotVector(
                    self._ldhToRvhTransformation3d, orientation
                )

            else:
                orientation = self.computeMatrixDotVector(self._rotation, orientation)

            self._diffusionGradientOrientations.append(orientation)

        self._receiverGain = self.getFloat(acqpLines, "RG")

    def getB0Count(self):
        return self._b0Count

    def getDiffusionGradientOrientations(self):
        return self._diffusionGradientOrientations

    def getBValues(self):
        return self._dwEffBValues

    def getReceiverGain(self):
        return self._receiverGain

    def saveBvecAndBval(self, fileNameBVec, fileNameBVal):
        f = open(fileNameBVec, "w")

        for o in self._diffusionGradientOrientations:
            f.write(str(o[0]) + " ")

        f.write("\n")

        for o in self._diffusionGradientOrientations:
            f.write(str(o[1]) + " ")

        f.write("\n")

        for o in self._diffusionGradientOrientations:
            f.write(str(o[2]) + " ")

        f.write("\n")

        f.close()

        f = open(fileNameBVal, "w")

        for b in self._dwEffBValues:
            f.write(str(b) + " ")

        f.write("\n")

        f.close()

    def getLines(self, fileNameLines, fieldName):
        startIndex = 0
        index = 0
        for line in fileNameLines:
            if "##$" + fieldName + "=" in line:
                startIndex = index
                break

            index += 1

        lines = []
        lines.append(fileNameLines[startIndex][len("##$" + fieldName + "=") : -1])

        endIndex = startIndex + 1
        while (
            "##$" not in fileNameLines[endIndex] and "$$" not in fileNameLines[endIndex]
        ):
            lines.append(fileNameLines[endIndex][:-1])
            endIndex += 1

        return lines

    def getMatrix(self, fileNameLines, fieldName):
        lines = self.getLines(fileNameLines, fieldName)

        sizeItems = lines[0][1:-1].split(",")

        items = []
        for line in lines[1:]:
            items += line.split()

        values = []
        for item in items:
            values.append(float(item))

        sizes = []
        if len(sizeItems) == 3:
            sizes = [int(sizeItems[0]), int(sizeItems[1]), int(sizeItems[2])]

        elif len(sizeItems) == 2:
            if int(sizeItems[0]) == 1:
                sizes = [
                    int(sizeItems[0]),
                    int(math.sqrt(float(sizeItems[1]))),
                    int(math.sqrt(float(sizeItems[1]))),
                ]

            else:
                sizes = [1, int(sizeItems[0]), int(sizeItems[1])]

        result = []
        if sizes[0] == 1:
            index = 0
            for i in range(0, sizes[1]):
                result.append([])
                for j in range(0, sizes[2]):
                    result[i].append(values[index])
                    index += 1

        else:
            print("getMatrix(): sizes[ 0 ] > 1 not managed")

        return result

    def getFloatVector(self, fileNameLines, fieldName):
        lines = self.getLines(fileNameLines, fieldName)

        size = int(lines[0][1:-1])

        items = []
        for line in lines[1:]:
            items += line.split()

        values = []
        for item in items:
            values.append(float(item))

        result = []
        for i in range(0, size):
            result.append(values[i])

        return result

    def getIntegerVector(self, fileNameLines, fieldName):
        lines = self.getLines(fileNameLines, fieldName)

        size = int(lines[0][1:-1])

        items = []
        for line in lines[1:]:
            items += line.split()

        values = []
        for item in items:
            values.append(int(item))

        result = []
        for i in range(0, size):
            result.append(values[i])

        return result

    def getStringVector(self, fileNameLines, fieldName):
        lines = self.getLines(fileNameLines, fieldName)

        items = []
        for line in lines[1:]:
            items += line.split()

        return items

    def getFloat(self, fileNameLines, fieldName):
        lines = self.getLines(fileNameLines, fieldName)

        return float(lines[0])

    def getInteger(self, fileNameLines, fieldName):
        lines = self.getLines(fileNameLines, fieldName)

        return int(lines[0])

    def getString(self, fileNameLines, fieldName):
        lines = self.getLines(fileNameLines, fieldName)
        if lines[0][0] == "(":
            return lines[1]

        return lines[0]

    def computeMatrixMultiplication(self, A, B):
        R = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

        for i in range(0, 3):
            for j in range(0, 3):
                for k in range(0, 3):
                    R[i][j] += A[i][k] * B[k][j]

        return R

    def computeMatrixDotVector(self, A, b):
        r = [0.0, 0.0, 0.0]

        for i in range(0, 3):
            for j in range(0, 3):
                r[i] += A[i][j] * b[j]

        return r

    def computeDeterminant(self, M):
        return (
            +M[0][0] * (+M[1][1] * M[2][2] - M[1][2] * M[2][1])
            - M[1][1] * (+M[0][0] * M[2][2] - M[0][2] * M[2][0])
            + M[0][0] * (+M[0][0] * M[1][1] - M[1][0] * M[0][1])
        )
