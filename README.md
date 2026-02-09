[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/neurospin/phcp-reco/main.svg)](https://results.pre-commit.ci/latest/github/neurospin/phcp-reco/main)
[![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/neurospin/phcp-reco/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/neurospin/phcp-reco)

# phcp-reco
Reconstruction scripts for post-mortem multi-quantitative fetal brain MRI at 11.7 teslas

## The p-HCP dataset

This repository contains the scripts that were used for preparing the p-HCP dataset:

> 🗃️ Arcamone, L. et al. Multimodal imaging of human fetal brain development at the mesoscopic scale using 11.7 T ex vivo MRI (v1). EBRAINS https://doi.org/10.25493/X3K8-Y9W (2025).

These scripts implement the data processing described in the following article, please refer to it for a description of the methods implemented here. Also, please cite it if you use these scripts for an academic publication:

> 📄 Arcamone, L. et al. Multimodal imaging of human fetal brain development at the mesoscopic scale using 11.7 T ex vivo MRI. bioRxiv. Preprint at https://doi.org/10.1101/2025.09.08.669657 (2025).


## Installation

Dependencies:

- Python 3.10
- [FSL](https://fsl.fmrib.ox.ac.uk/) (tested with version 6.0.7.13)
- [ANTS](https://github.com/ANTsX/ANTs) (tested with version 2.5.4)
- [Gkg](https://framagit.org/cpoupon/gkg) (locally distributed Singularity containers)
- [dcm2bids](https://github.com/unfmontreal/Dcm2Bids) (tested with version 3.2.0)
- [dcm2niix](https://github.com/rordenlab/dcm2niix) (tested with version 1.0.20250505)

```shell
git clone https://github.com/neurospin/phcp-reco.git

# Installation in a virtual environment
cd phcp-reco
python3 -m venv venv/
. venv/bin/activate
pip install -e .
```

You may need to adapt the contents of [phcp/config.py](phcp/config.py) to suit your local set-up.


## Data layout

### Naming of input data

Each block should be assigned a short `BlockName` to identify it, based on its anatomical location, e.g. `Inf`, `Mid`, `Sup`, `InfLat`, `SupMed`, etc.

Each FOV within a block should be assigned a `FovName` that starts with `BlockName` and further specifies the position of the FOV within the block, e.g. `InfAnt`, `InfPos`, etc.

Each imaging session corresponds to a single FOV and should be named with a concatenation of the acquisition date in `YYYYMMDD` format, and `FovName`, e.g. `20250101InfPos`.


### Source data

The images of each field of view are first exported from the console in raw Bruker format and DICOM. The DICOM data is used for mostly everything, except that some acquisition parameters cannot be found in the DICOM metadata, namely the diffusion b-values and b-vectors, as well as the receiver gain. In those cases, we resort to reading them from the Bruker files `acqp`, `method` and `visu_pars`.

```
fov/sourcedata
└── sub-{subjectID}
    └── ses-{sessionID}
        └── {ScanNumber}
            ├── acqp
            ├── method
            └── pdata
                └── 1
                    └── dicom
                        └── MRIm*.dcm
```

The majority of the DICOM files are converted to BIDS format using `dcm2bids`. Configuration files are included in the `dcm2bids` directory for the 11.7 T, 7 T, and in-situ scans.

```shell
dcm2bids \
    -p ${sub} \
    -s ${ses} \
    -c code/phcp-reco/dcm2bids/dcm2bids_11.7T.json \
    -o rawdata/ \
    -d sourcedata/sub-${sub}/ses-${ses}/*/pdata/1/dicom/
```



### Raw data

The raw data is organized using the following BIDS layout, where each BIDS “session” corresponds to a given field of view:

```
fov/rawdata
└── sub-{subjectID}
    └── ses-{sessionID}
        ├── anat
        │   ├── sub-{subjectID}_ses-{sessionID}_echo-{echo}_MEGRE.json
        │   ├── sub-{subjectID}_ses-{sessionID}_echo-{echo}_MEGRE.nii.gz
        │   ├── sub-{subjectID}_ses-{sessionID}_echo-{echo}_MESE.json
        │   ├── sub-{subjectID}_ses-{sessionID}_echo-{echo}_MESE.nii.gz
        │   ├── sub-{subjectID}_ses-{sessionID}_flip-{fa}_VFA.json
        │   ├── sub-{subjectID}_ses-{sessionID}_flip-{fa}_VFA.nii.gz
        │   ├── sub-{subjectID}_ses-{sessionID}_T2w.json
        │   └── sub-{subjectID}_ses-{sessionID}_T2w.nii.gz
        ├── dwi
        │   ├── sub-{subjectID}_ses-{sessionID}_acq-{bval}_run-{run}_dwi.json
        │   ├── sub-{subjectID}_ses-{sessionID}_acq-{bval}_run-{run}_dwi.nii.gz
        └── fmap
            ├── sub-{subjectID}_ses-{sessionID}_TB1AFI.json
            └── sub-{subjectID}_ses-{sessionID}_TB1AFI.nii.gz
```

### Header-flipped raw data

The spatial transformation contained in the DICOM metadata, and thus Nifti `rawdata`, does not match the anatomical axes. The reason is that the ParaVision software does not allow us to input the real orientation of the tissue block in the scanner, so we default to choosing "head first supine" (HFS). Therefore, prior to any further processing we flip the the Nifti header transforms, so that the axes of the new "scanner-based" referential correspond to sensible anatomical orientations. This is done so that visualizing the images the usual viewers will display the brain in a sensible orientation.

The `phcp-header-flip` script will mirror the BIDS structure of the provided `rawdata` in a new `headerfliprawdata` directory.

#### Usage

```shell
phcp-header-flip \
    --sessionDirectory fov/rawdata/sub-${sub}/ses-${ses} \
    --outputDirectory fov/headerfliprawdata
```

#### Outputs


```
fov/headerfliprawdata
└── sub-{subjectID}
    └── ses-{sessionID}
        ├── anat
        │   ├── sub-{subjectID}_ses-{sessionID}_echo-{echo}_MEGRE.json
        │   ├── sub-{subjectID}_ses-{sessionID}_echo-{echo}_MEGRE.nii.gz
        │   ├── sub-{subjectID}_ses-{sessionID}_echo-{echo}_MESE.json
        │   ├── sub-{subjectID}_ses-{sessionID}_echo-{echo}_MESE.nii.gz
        │   ├── sub-{subjectID}_ses-{sessionID}_flip-{fa}_VFA.json
        │   ├── sub-{subjectID}_ses-{sessionID}_flip-{fa}_VFA.nii.gz
        │   ├── sub-{subjectID}_ses-{sessionID}_T2w.json
        │   └── sub-{subjectID}_ses-{sessionID}_T2w.nii.gz
        ├── dwi
        │   ├── sub-{subjectID}_ses-{sessionID}_acq-{bval}_run-{run}_dwi.json
        │   ├── sub-{subjectID}_ses-{sessionID}_acq-{bval}_run-{run}_dwi.nii.gz
        └── fmap
            ├── sub-{subjectID}_ses-{sessionID}_TB1AFI.json
            └── sub-{subjectID}_ses-{sessionID}_TB1AFI.nii.gz
```


## The field-of-view processing pipeline

### Early pre-processing

Early pre-processing of the raw data is performed with the `phcp-bruker-preprocessing` script, which includes:

- extraction of the b-vectors and b-values from the Bruker source data (`acqp` and `method`, `visu_pars` files) into standard BIDS-format `.bvec` and `.bval` files;
- extraction of receiver gain values from the Bruker source data (`acqp` file) for the diffusion and VFA images, and homogenization of the values to correct for potential receiver gain differences;
- concatenation of the segmented diffusion scans and of the individual VFA scans.

#### Inputs

The `phcp-bruker-preprocessing` script uses a few JSON files as configuration. You have to create a `description` directory, which can be stored under `derivatives/gkg-Pipeline/PreprocessingDescriptions`:

```
derivatives/gkg-Pipeline/PreprocessingDescriptions
├── sub-{subjectID}
│   └── RG_description.json
└── sub-{subjectID}.json
```

- `RG_description.json` contains a list of sessions with the list of modalities to be processed for each:
```json
{
  "sub-{subjectID}_ses-{sessionID}": [
    "b1500",
    "b4500",
    "b8000",
    "VFA"
  ]
}
```

Additionally, the “subject file”, whose name is free but suggested to be in the form `sub-{subjectID}.json`, contains the list of subjects/sessions to be processed in the following form:

```json
{
  "sub-{subjectID}": ["ses-{sessionID}"]
}
```

#### Usage

```shell
phcp-brucker-preprocessing \
    --sourcedata fov/sourcedata \
    --rawdata fov/headerfliprawdata \
    --derivatives fov/derivatives \
    --subject derivatives/gkg-Pipeline/PreprocessingDescriptions/sub-${sub}.json \
    --descriptions derivatives/gkg-Pipeline/PreprocessingDescriptions/sub-${sub}
```

#### Outputs

```
fov/headerfliprawdata
└── sub-{subjectID}
    └── ses-{sessionID}
        └── dwi
            ├── sub-{subjectID}_ses-{sessionID}_acq-{bval}_run-{run}_dwi.bvec
            └── sub-{subjectID}_ses-{sessionID}_acq-{bval}_run-{run}_dwi.bval
fov/derivatives/RGCorrection
└── sub-{subjectID}
    └── ses-{sessionID}
        ├── RG_values.json
        └── sub-{subjectID}_ses-{sessionID}_VFA.nii.gz
fov/derivatives/gkg-Pipeline/
└── sub-{subjectID}
    └── ses-{sessionID}
        ├── b1500
        │   └── 01-Materials
        │       ├── sub-{subjectID}_ses-{sessionID}_acq-b1500_dwi.bval
        │       ├── sub-{subjectID}_ses-{sessionID}_acq-b1500_dwi.bvec
        │       └── sub-{subjectID}_ses-{sessionID}_acq-b1500_dwi.nii.gz
        ├── b4500
        │   └── 01-Materials
        │       ├── sub-{subjectID}_ses-{sessionID}_acq-b4500_dwi.bval
        │       ├── sub-{subjectID}_ses-{sessionID}_acq-b4500_dwi.bvec
        │       └── sub-{subjectID}_ses-{sessionID}_acq-b4500_dwi.nii.gz
        └── b8000
            └── 01-Materials
                ├── sub-{subjectID}_ses-{sessionID}_acq-b8000_dwi.bval
                ├── sub-{subjectID}_ses-{sessionID}_acq-b8000_dwi.bvec
                └── sub-{subjectID}_ses-{sessionID}_acq-b8000_dwi.nii.gz
```

### B1 mapping

The transmit B1 field inhomogeneities are mapped using an Actual Flip Angle sequence. The following script reconstructs the B1 map from the AFI data. Parameters (TR1 and TR2) are not in the JSON metadata, they have to be passed manually but default values match the protocol.

#### Inputs

```
fov/headerfliprawdata
└── sub-{subjectID}
    └── ses-{sessionID}
        └── fmap
            ├── sub-{subjectID}_ses-{sessionID}_TB1AFI.json
            └── sub-{subjectID}_ses-{sessionID}_TB1AFI.nii.gz
```

#### Usage

```shell
mkdir -p fov/derivatives/fov-reconstructed/sub-${sub}/ses-${ses}/fmap
phcp-fov-afi-b1mapping \
    fov/headerfliprawdata/sub-${sub}/ses-${ses}/fmap/sub-${sub}_ses-${ses}_TB1AFI.nii.gz \
    fov/derivatives/fov-reconstructed/sub-${sub}/ses-${ses}/fmap/sub-${sub}_ses-${ses}_TB1map.nii.gz
```

#### Outputs

```
fov/derivatives/fov-reconstructed
└── sub-{subjectID}
    └── ses-{sessionID}
        └── fmap
            └── sub-{subjectID}_ses-{sessionID}_TB1map.nii.gz
```

### T2* mapping

#### Inputs

```
fov/headerfliprawdata
└── sub-{subjectID}
    └── ses-{sessionID}
        └── anat
            ├── sub-{subjectID}_ses-{sessionID}_echo-{echo}_MEGRE.json
            └── sub-{subjectID}_ses-{sessionID}_echo-{echo}_MEGRE.nii.gz
```

#### Usage

```shell
mkdir -p fov/derivatives/T2starmapping/sub-${sub}/ses-${ses}/01-Materials
fslmerge -t fov/derivatives/T2starmapping/sub-${sub}/ses-${ses}/01-Materials/t2star-mge.nii.gz fov/headerfliprawdata/sub-${sub}/ses-${ses}/anat/*_MEGRE.nii.gz
mkdir -p fov/derivatives/T2starmapping/sub-${sub}/ses-${ses}/02-Results
phcp-t2star-relaxometry \
    --input fov/derivatives/T2starmapping/sub-${sub}/ses-${ses}/01-Materials \
    --mge fov/headerfliprawdata/sub-${sub}/ses-${ses}/anat/'*_MEGRE.json' \
    --outputDirectory fov/derivatives/T2starmapping/sub-${sub}/ses-${ses}/02-Results
```

#### Outputs

```
fov/derivatives/T2starmapping
└── sub-{subjectID}
    └── ses-{sessionID}
        └── 02-Results
            ├── proton-density.nii.gz
            ├── fitted-mge.nii.gz
            ├── T2star.nii.gz
            └── T2starConfidenceMap.nii.gz
```

### T2 mapping

#### Inputs

```
fov/headerfliprawdata
└── sub-{subjectID}
    └── ses-{sessionID}
        └── anat
            ├── sub-{subjectID}_ses-{sessionID}_echo-{echo}_MESE.json
            └── sub-{subjectID}_ses-{sessionID}_echo-{echo}_MESE.nii.gz
```

#### Usage

```shell
mkdir -p fov/derivatives/T2mapping/sub-${sub}/ses-${ses}/01-Materials
fslmerge -t fov/derivatives/T2mapping/sub-${sub}/ses-${ses}/01-Materials/t2-msme.nii.gz fov/headerfliprawdata/sub-${sub}/ses-${ses}/anat/*_MESE.nii.gz
mkdir -p fov/derivatives/T2mapping/sub-${sub}/ses-${ses}/02-Results
phcp-t2-relaxometry \
    --input fov/derivatives/T2mapping/sub-${sub}/ses-${ses}/01-Materials \
    --msme fov/headerfliprawdata/sub-${sub}/ses-${ses}/anat/'*_MESE.json' \
    --outputDirectory fov/derivatives/T2mapping/sub-${sub}/ses-${ses}/02-Results
```

#### Outputs

```
fov/derivatives/T2mapping
└── sub-{subjectID}
    └── ses-{sessionID}
        └── 02-Results
            ├── proton-density.nii.gz
            ├── fitted-msme.nii.gz
            ├── T2.nii.gz
            └── T2ConfidenceMap.nii.gz
```

### T1 mapping

#### Inputs

```
fov/headerfliprawdata
└── sub-{subjectID}
    └── ses-{sessionID}
        └── anat
            └── sub-{subjectID}_ses-{sessionID}_flip-{fa}_VFA.json
fov/derivatives/RGCorrection
└── sub-{subjectID}
    └── ses-{sessionID}
        └── sub-{subjectID}_ses-{sessionID}_VFA.nii.gz
fov/derivatives/fov-reconstructed
└── sub-{subjectID}
    └── ses-{sessionID}
        └── fmap
            └── sub-{subjectID}_ses-{sessionID}_TB1map.nii.gz
```

#### Usage

```shell
mkdir -p fov/derivatives/T1mapping/sub-${sub}/ses-${ses}/01-Materials
cp fov/derivatives/RGCorrection/sub-${sub}/ses-${ses}/sub-${sub}_ses-${ses}_VFA.nii.gz \
   fov/derivatives/T1mapping/sub-${sub}/ses-${ses}/01-Materials/t1map.nii.gz
cp fov/derivatives/fov-reconstructed/sub-${sub}/ses-${ses}/fmap/sub-${sub}_ses-${ses}_TB1map.nii.gz \
    fov/derivatives/T1mapping/sub-${sub}/ses-${ses}/01-Materials/b1.nii.gz
mkdir -p fov/derivatives/T1mapping/sub-${sub}/ses-${ses}/02-Results
phcp-t1-relaxometry \
    --input fov/derivatives/T1mapping/sub-${sub}/ses-${ses}/01-Materials \
    --vfa fov/headerfliprawdata/sub-${sub}/ses-${ses}/anat/'*_VFA.json' \
    --outputDirectory fov/derivatives/T1mapping/sub-${sub}/ses-${ses}/02-Results
```

#### Outputs

```
fov/derivatives/T1mapping
└── sub-{subjectID}
    └── ses-{sessionID}
        └── 02-Results
            ├── proton-density.nii.gz
            ├── fitted-vfa.nii.gz
            ├── T1.nii.gz
            ├── T1_rec-unbiased.nii.gz
            └── T1ConfidenceMap.nii.gz
```

### Diffusion imaging pipeline

The diffusion imaging pipeline takes care of the whole chain, from the pre-processing (artefact correction) to the diffusion models (DTI, NODDI, and tractography). JSON files are used to configure which steps to run and pass parameters.

#### Inputs

The diffusion data itself, corrected for receiver gain, is read from this hierarchy (see `phcp-bruker-preprocessing` above):

```
fov/derivatives/gkg-Pipeline/
└── sub-{subjectID}
    └── ses-{sessionID}
        ├── b1500
        │   └── 01-Materials
        │       ├── sub-{subjectID}_ses-{sessionID}_acq-b1500_dwi.bval
        │       ├── sub-{subjectID}_ses-{sessionID}_acq-b1500_dwi.bvec
        │       └── sub-{subjectID}_ses-{sessionID}_acq-b1500_dwi.nii.gz
        ├── b4500
        │   └── 01-Materials
        │       ├── sub-{subjectID}_ses-{sessionID}_acq-b4500_dwi.bval
        │       ├── sub-{subjectID}_ses-{sessionID}_acq-b4500_dwi.bvec
        │       └── sub-{subjectID}_ses-{sessionID}_acq-b4500_dwi.nii.gz
        └── b8000
            └── 01-Materials
                ├── sub-{subjectID}_ses-{sessionID}_acq-b8000_dwi.bval
                ├── sub-{subjectID}_ses-{sessionID}_acq-b8000_dwi.bvec
                └── sub-{subjectID}_ses-{sessionID}_acq-b8000_dwi.nii.gz
```

The pipeline is configured using three JSON files:

```
derivatives/gkg-Pipeline/
├── PreprocessingDescriptions
│   └── sub-{subjectID}.json
└── GkgPipelineDescriptions
    ├── GkgPipeline_description_{subjectID}.json
    └── tasks.json
```

- `sub-{subjectID}.json` is the “subject file” described above in the preprocessing section.

- `GkgPipeline_description_${subjectID}.json` contains a list of sessions with the list of modalities to be processed for each:
```json
{
  "sub-{subjectID}_ses-{sessionID}": [
    "b1500",
    "b4500",
    "b8000"
  ]
}
```
- `tasks.json` contains the list of steps that should be run:

```json
{
  "OrientationAndBValueFileDecoding": 1,
  "NonLocalMeansFiltering": 1,
  "QSpaceSamplingAddition": 1,
  "Morphologist": 1,
  "OutlierCorrection": 1,
  "EddyCurrentCorrection": 1,
  "LocalModelingDTI": 1,
  "LocalModelingQBI": 0
}
```

#### Usage

```shell
phcp-dffusion-pipeline \
    --subjectJsonFileName fov/derivatives/gkg-Pipeline/PreprocessingDescriptions/sub-${sub}.json \
    --taskJsonFileName fov/derivatives/gkg-Pipeline/GkgPipelineDescriptions/tasks.json \
    --gkgpipelineJsonFilename fov/derivatives/GkgPipelineDescriptions/GkgPipeline_description_${sub}.json\
    --outputDirectory fov/derivatives/gkg-Pipeline/
```

#### Outputs

```
fov/derivatives/gkg-Pipeline/
└── sub-{subjectID}
    └── ses-{sessionID}
        └── b1500 / b4500 / b8000 / multishell
            ├── 02-OrientationAndBValueFileDecoding
            ├── 03-NonLocalMeansFiltering
            ├── 04-QSpaceSamplingAddition
            ├── 05-MorphologistPipeline
            ├── 06-OutlierCorrection
            ├── 07-EddyCurrentAndMotionCorrection
            ├── 08-LocalModeling-DTI
            ├── 09-LocalModeling-QBI
            └── ...
```

### Outputs of the field-of-view processing pipeline

```
fov/derivatives/fov-reconstructed
└── sub-{subjectID}
    └── ses-{sessionID}
        └── fmap
            └── sub-{subjectID}_ses-{sessionID}_TB1map.nii.gz
```


## The reconstruction pipeline

The reconstruction pipeline consists of:
1. Obtaining the deformation fields that bring all fields of view into spatial correspondence. The registration strategy consists of three steps described below: intra-FOV registration, inter-FOV registration, and inter-block registration.
2. Applying these deformation fields to each field of view, and merging the data into a single image.


### Intra-FOV registration

Intra-FOV registration aims to align each modality to the modality of reference (T2-weighted). The main steps are:

1. Placing the modality of reference (T2-weighted) in 'fixed' folder and all other modalities inside the 'moving' folder;
2. Running the `phcp-intrafov-registration` script, preforming format conversion when needed (GIS to NIfTI), bias correction and intra-FOV registration.

#### Inputs of intra-FOV registration

The inputs are the two main folders ('fixed' and 'moving') containing modalities to register.
```
derivatives/Registration/sub-{subjectID}/IntraFOV_Registration/{FovName}/
├── fixed
│   └── sub-{subjectID}_ses-{sessionIDFOV1}_T2w.nii.gz
└── moving
    ├── {modality-prefix}_{FovName1}.nii.gz or .ima
    ├── {modality-prefix}_{FovName2}.nii.gz or .ima
    └── ...
```

Rename your moving using the name of the modality as the first prefix (ex: `VFA_flip-90.nii.gz`, `MSME_echo-1.nii.gz`, ...).
#### Usage

```shell
phcp-intrafov-registration --verbose \
    -i derivatives/Registration/sub-${sub}/IntraFOV_Registration/${FovName}
```

#### Outputs of intra-FOV registration
Intra-FOV creates transformation files in the 'transform_files' folder, which is created automatically. Transformation files will have the following name architecture : `{modality-prefix}_To_T2w_{object}.{extension}`. `object` refers to the type of file created and can take the following values: `0GenericAffine`, `1Warp`, `1InverseWarp`, `Regis`, `INV_Regis`. And `extension` can be `mat` if it is a linear transformation file and `nii.gz` if it is an NIfTI image such as the warped image or the deformation field.

<!-- TODO add a word on how to perform quality control on the intra-FOV registration -->

### Inter-FOV registration

Inter-FOV registration aims to reconstruct a whole-block image from its constituent FOVs. The main steps are:

1. Placing the input data and creating a JSON file describing the input FOVs and the approximate translation between them;
2. Running the `phcp-interfov-registration` script, which performs bias correction, denoising, gradient non-linearity correction, the inter-FOV translation, runs the actual registration, and initializes the geometry of the whole-block space;
3. Creating a transformation file for placing the reference (most posterior) FOV into the whole-block space, and optionally fine-tuning the fusion parameters;
4. Running `phcp-interfov-fusion` to apply the registration fuse all FOVs, yielding a whole-block T2-weighted image.

#### Inputs of inter-FOV registration

The inputs are the pre-calculated deformation field for correction of gradient non-linearity (`WarpHeaderFile_FWHM1.nii.gz`), a `Description.json` file described below, and the T2-weighted images of all FOVs. They should be prepared according to this layout:

```
derivatives/Registration/sub-{subjectID}/InterFOV_Registration/{BlockName}/
├── 01-Materials
│   ├── sub-{subjectID}_ses-{sessionIDFOV1}_T2w.nii.gz
│   ├── sub-{subjectID}_ses-{sessionIDFOV2}_T2w.nii.gz
│   ├── sub-{subjectID}_ses-{sessionIDFOV3}_T2w.nii.gz
│   ├── ...
│   └── WarpHeaderFile_FWHM1.nii.gz
└── Description.json
```

A `Description.json` file should be prepared, listing the FOVs of the block and their relative translation values. In the example below, a block was acquired in 3 FOVs with 30 mm translation between each FOV. Note that the reference FOV (the most posterior FOV) needs to be identified with a **0.0** translation, and **other FOVs must have non-zero translation**.

```json
{
    "sub-{subjectID}_ses-{sessionIDFOV1}_T2w": 0.0,
    "sub-{subjectID}_ses-{sessionIDFOV2}_T2w": 30.0,
    "sub-{subjectID}_ses-{sessionIDFOV3}_T2w": 60.0,
}
```

#### Usage

```shell
phcp-interfov-registration --verbose \
    -i derivatives/Registration/sub-{subjectID}/InterFOV_Registration/${BlockName} \
    -j derivatives/Registration/sub-{subjectID}/InterFOV_Registration/${BlockName}/Description.json
```

#### Outputs of inter-FOV registration

Each FOV is denoised with `NLMF` (non-local means filtering) and debiased with `N4`, corrected for gradient non-linearities, and the translation is applied to the header, yielding a preprocessed output with `_NLMF_N4_rec-unwarp` suffix. Note that the deformation fields are estimated and saved to the `03-TransformFiles` directory, but they are not applied to the data at this stage: therefore, quality control takes place *after* the next step of inter-FOV fusion.

```
derivatives/Registration/sub-{subjectID}/InterFOV_Registration/{BlockName}/
├── 02-Preparation
│   ├── sub-{subjectID}_ses-{sessionIDFOV1}_T2w_NLMF_N4_rec-unwarp.nii.gz
│   ├── sub-{subjectID}_ses-{sessionIDFOV2}_T2w_NLMF_N4_rec-unwarp.nii.gz
│   ├── sub-{subjectID}_ses-{sessionIDFOV3}_T2w_NLMF_N4_rec-unwarp.nii.gz
│   ├── sub-{subjectID}_ses-{sessionIDFOV1}_T2w_sub-{subjectID}_ses-{sessionIDFOV1}_T2w.json
│   ├── sub-{subjectID}_ses-{sessionIDFOV2}_T2w_sub-{subjectID}_ses-{sessionIDFOV3}_T2w.json
│   └── ...
├── 03-TransformFiles
│   └── ...
└── IntermediateSpace.nii.gz
```

#### Inputs of inter-FOV fusion

The previous process creates a reference geometry for the whole-block image, `IntermediateSpace.nii.gz`, but that reference is not derived from the FOV geometry. Therefore, the image of the reference FOV (normally the most posterior FOV) must be manually registered into that intermediate space. This can be done by:
1. Loading `IntermediateSpace.nii.gz` in ITK-SNAP as a main image, using the Jet colormap so that the background of that image appears blue on the black background (beware that this uniform image causes a crash of ITK-SNAP if Contrast inspector or Auto adjust constrast are used).
2. Loading `02-Preparation/sub-{subjectID}_ses-{sessionIDRefFOV}_T2w_NLMF_N4_rec-unwarp.nii.gz` as an additional image in overlay mode.
3. Using the Registration module of ITK-SNAP in Rigid mode, to perform manual alignment of the reference FOV at the extremity of the intermediate space. Make sure to hit *Zoom to fit* so that you can visualize the full whole-block space.
4. Saving the resulting transformation matrix into a file named `refFOV_to_intermediateSpace.txt`.

Inter-FOV fusion uses the following files:

```
derivatives/Registration/sub-{subjectID}/InterFOV_Registration/{BlockName}/
├── 01-Materials
│   ├── sub-{subjectID}_ses-{sessionIDFOV1}_T2w_NLMF_N4_rec-unwarp.nii.gz
│   ├── sub-{subjectID}_ses-{sessionIDFOV2}_T2w_NLMF_N4_rec-unwarp.nii.gz
│   ├── sub-{subjectID}_ses-{sessionIDFOV3}_T2w_NLMF_N4_rec-unwarp.nii.gz
│   ├── ...
├── 03-TransformFiles
│   └── ...
├── Description.json
├── IntermediateSpace.nii.gz
└── refFOV_to_intermediateSpace.txt
```

#### Usage

```shell
phcp-interfov-fusion \
    -i derivatives/Registration/sub-{subjectID}/InterFOV_Registration/${BlockName} \
    -j derivatives/Registration/sub-{subjectID}/InterFOV_Registration/${BlockName}/Description.json
```

#### Outputs of inter-FOV fusion

```
derivatives/Registration/sub-{subjectID}/InterFOV_Registration/{BlockName}/
├── 04-FovsDistribution
│   └── ...
└── Reconstructed_block.nii.gz
```

You should check the quality of the reconstructed block at this stage, and make adjustments to the fusion parameters if necessary (currently this requires editing the scripts themselves).


### Inter-block registration

The registration strategy is described in the data paper, but the details of registration can differ for different specimens. As a result, no registration scripts are distributed at the moment.


### Concatenation of transformations (optional)

The `phcp-concat-transforms` script concatenates a series of linear and non-linear transformations applied to a single field of view (FOV), producing two key output files:
- `total_affine_transform.txt`
- `total_deformation_field.nii.gz`

#### Example
```shell
phcp-concat-transforms \
    --input sub-${sub}_ses-${ses}_T2w.nii.gz \
    --json transform_filenames_sorted.json \
    --output /output_directory_path/
```

#### Required Inputs
- `sub-{sub}_ses-{ses}_T2w.nii.gz` defines the **initial (native) space**. The output `total_deformation_field.nii.gz` space will be the same.
- `transform_filenames_sorted.json` lists the successive transformation files, sorted in the order they should be applied (from first to last). The file should follow this structure:
```json
{
  "tparams": [
    "First_transform_filename(.txt or .nii.gz)",
    "Second_transform_filename(.txt or .nii.gz)",
    "...",
    "Last_transform_filename(.txt or .nii.gz)"
  ]
}
```


#### Output Files and Usage
The script generates the following **main outputs**:
1. `total_deformation_field.nii.gz` – A non-linear deformation field in ANTs format, usable directly with `antsApplyTransforms`.
2. `total_affine_transform.txt` – A text file representing the combined affine transformation.

To apply the transformations correctly using `antsApplyTransforms`, use the following order:
1. `total_deformation_field.nii.gz`
2. `total_affine_transform.txt`

#### Additional Files
- **Secondary files** (not discussed in this repository but included for reference): `jacobian_deformation_field.nii.gz`, `jacobianLog_deformation_field.nii.gz`, `total_deformation_field_smoothed.nii.gz`
- **Tertiary files**: Intermediate files generated during processing; not required for downstream analysis.


### Transformation into the final space and fusion

The `phcp-transform-and-fuse` script merges multiple fields of view (FOVs) from the **initial space** into the **final space**. It consists of two main stages, controlled via the `--run` flag.

#### Example Workflow

```shell
# Step 1: Prepare final-space materials (default mode)
phcp-transform-and-fuse -p fov/derivatives/Registration/sub-${sub}/fusion/

# Step 2: Merge using 3 blocks
phcp-transform-and-fuse -p fov/derivatives/Registration/sub-${sub}/fusion/ --run -n 3
```

#### Stage 1: Prepare Materials (Default, without `--run`)

This stage generates all necessary materials for reconstruction in the final space from T2star modality in the initial space, including:

- Geometric penalty maps
- Each modality transformed into the final space

All output will be stored in a folder created automatically and named `02-RefSpace`.

##### Required Inputs
Considering `Fusion` folder as the working directory :
```
derivatives/Registration
└── sub-{subjectID}
    └── Fusion
        ├── 01-InitSpace/
        │   ├── {modality}_{FovName}.nii.gz
        │   └── ...
        ├── SendToRefSpace_{modality}.json
        └── SendToRefSpace_...
```

###### 1. `01-InitSpace/` folder

Must contain your FOV data in the initial space files, for example:

```
T2star_InfPos.nii.gz
T2star_InfMid.nii.gz
...
```

- File naming format: `{modality}_{FovName}.nii.gz`, where `{modality}` in [ `QT1`, `QT2`, `QT2star`, `ADC1500`, `ADC4500`, `ADC8000`, `FA1500`, `FA4500`, `FA8000`, `TransDiff1500`, `TransDiff4500`, `TransDiff8000`,  `ParaDiff1500`, `ParaDiff4500`, `ParaDiff8000`, `GFA1500`, `GFA4500`, `GFA8000`, `T2w`] and `{FovName}` is the FOV name as described above.
- Supported formats: **GIS** or **NIfTI**
- **Must contain at least T2star modality**

###### 2. `SendToRefSpace_{modality}.json` files

You must provide one JSON file per modality (`QT1`, `QT2`, `QT2star`, `DWI`-one key for all DWI modalities-, `T2W`), with the following structure:

```json
{
  "RefSpace": "refspace_filename_path",
  "InfPos": {
    "tparams": ["last_transformation_path", "...", "first_transformation_path"]
  },
  "InfMid": {
    "tparams": ["..."]
  }
}
```
- `RefSpace` points to a Nifti file that specifies the geometry of the whole-brain space that was used for the registration process

- `tparams` is a list of filenames that describe each step of the transformation chain. Each filename can point to either a linear transformation (`.txt` or `.mat` or a non-linear deformation field (`.nii.gz`). The files must be sorted from the last transformation file to the first.

##### Stage 2: Merge Blocks (with `--run` and `-n`)

**This stage must be performed after the stage 1.** Merges the transformed data located in `02-RefSpace` using:

- **Geometric penalties** for all modalities
- **Geometric penalties** & **Goodness-of-fit weighting** for relaxometric maps

Stores the reconstructed blocks in `03-Blocks` and the final reconstruction in `04-Reconstruction`.

###### Required Inputs

```
derivatives/Registration
└── sub-{subjectID}
    └── Fusion
        ├── 02-RefSpace/
        │   └── ...
        ├── block_{i}.json
        └── block...
```

###### 1. `block_{i}.json` files

For `n` blocks, create `n` files named:

```
block_1.json
block_2.json
...
block_n.json
```

Each file should follow this structure:

```json
{
  "mask": ["mask_file_1", "...", "mask_file_n"],
  "modality_1": ["modality1_file_1", "...", "modality1_file_n"],
  "modality_2": ["..."]
}
```

###### Important Notes:

- File paths must be **ordered continuously** from one anatomical end to the other (e.g., `InfPos` > `InfMid` > `InfAnt`).
- Reverse order (e.g., `InfAnt` > `InfMid` > `InfPos`) is also valid, as long as the sequence is consistent.
- Do **not skip any FOVs**.


## Data publication

Prepare the dataset for publication (reidentification, Nifti encoding, adjustment of headers):

```shell
phcp-prepare-published-dataset \
    input_directory \
    output_bids_directory \
    publication.yaml
```


## Contributing

This repository uses [pre-commit](https://pre-commit.com/) to ensure that all committed code follows basic quality standards. Please install it and configure it to run as part of ``git commit`` by running ``pre-commit install`` in your local repository:

```shell
pipx install pre-commit
pre-commit install  # set up the hook (to run pre-commit during 'git commit')
```
