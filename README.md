# phcp-reco
Reconstruction scripts for post-mortem multi-quantitative fetal brain MRI at 11.7 teslas

## Installation

Dependencies:

- Python 3.10
- [FSL](https://fsl.fmrib.ox.ac.uk/) (tested with version 6.0.7.13)
- [ANTS](https://github.com/ANTsX/ANTs) (tested with version 2.5.4)
- [Gkg](https://framagit.org/cpoupon/gkg) (locally distributed Singularity containers)

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

The DICOM files are converted to BIDS format. TODO: document the conversion process

- dcm2niix
- Parse Bruker raw data and write bvec, bval, ReceiverGain

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
│   ├── RG_description.json
│   └── sub-{subjectID}_ses-{sessionID}_description.json
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
- `sub-{subjectID}_ses-{sessionID}_description.json` contains the list of Bruker scan numbers associated with each run of the segmented diffusion, for example::

```json
{
  "b1500": [97, 98],
  "b4500": [91, 92, 93, 94, 95, 96],
  "b8000": [82, 83, 84, 85, 86, 87, 88, 89, 90]
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
    --rawdata fov/rawdata \
    --derivatices fov/derivatives \
    --subject descriptions/sub-{subjectID}.json \
    --descriptions descriptions/sub-{subjectID}
```

#### Outputs

```
fov/rawdata
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
fov/rawdata
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
    fov/rawdata/sub-${sub}/ses-${ses}/fmap/sub-${sub}_ses-${ses}_TB1AFI.nii.gz \
    fov/derivatives/fov-reconstructed/sub-${sub}/ses-${ses}/fmap/sub-${sub}_ses-${ses}_TB1map.nii.gz
```

The T1 mapping script expects to find the B1 map in a specific location, so you should copy the map over:

```shell
mkdir -p fov/derivatives/T1mapping/sub-${sub}/ses-${ses}/01-Materials
cp fov/derivatives/fov-reconstructed/sub-${sub}/ses-${ses}/fmap/sub-${sub}_ses-${ses}_TB1map.nii.gz \
    fov/derivatives/T1mapping/sub-${sub}/ses-${ses}/01-Materials/b1.nii.gz
```

#### Outputs

```
fov/derivatives/fov-reconstructed
└── sub-{subjectID}
    └── ses-{sessionID}
        └── fmap
            └── sub-{subjectID}_ses-{sessionID}_TB1map.nii.gz
fov/derivatives/T1mapping
└── sub-{subjectID}
    └── ses-{sessionID}
        └── 01-Materials
            └── b1.nii.gz
```

### T2* mapping

#### Inputs

```
fov/rawdata
└── sub-{subjectID}
    └── ses-{sessionID}
        └── anat
            ├── sub-{subjectID}_ses-{sessionID}_echo-{echo}_MEGRE.json
            └── sub-{subjectID}_ses-{sessionID}_echo-{echo}_MEGRE.nii.gz
```

#### Usage

```shell
fslmerge -t fov/derivatives/T2starmapping/sub-${sub}/ses-${ses}/01-Materials/t2star-mge.nii.gz fov/rawdata/sub-${sub}/ses-${ses}/anat/*_MEGRE.nii.gz
mkdir -p fov/derivatives/T2starmapping/sub-${sub}/ses-${ses}/02-Results
phcp-t2star-relaxometry \
    --input fov/derivatives/T2starmapping/sub-${sub}/ses-${ses}/01-Materials \
    --mge fov/rawdata/sub-${sub}/ses-${ses}/anat/'*_MEGRE.json' \
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
fov/rawdata
└── sub-{subjectID}
    └── ses-{sessionID}
        └── anat
            ├── sub-{subjectID}_ses-{sessionID}_echo-{echo}_MESE.json
            └── sub-{subjectID}_ses-{sessionID}_echo-{echo}_MESE.nii.gz
```

#### Usage

```shell
fslmerge -t fov/derivatives/T2mapping/sub-${sub}/ses-${ses}/01-Materials/t2-msme.nii.gz fov/rawdata/sub-${sub}/ses-${ses}/anat/*_MESE.nii.gz
mkdir -p fov/derivatives/T2mapping/sub-${sub}/ses-${ses}/02-Results
phcp-t2-relaxometry \
    --input fov/derivatives/T2mapping/sub-${sub}/ses-${ses}/01-Materials \
    --msme fov/rawdata/sub-${sub}/ses-${ses}/anat/'*_MESE.json' \
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

### Outputs of the field-of-view processing pipeline

```
fov/derivatives/fov-reconstructed
└── sub-{subjectID}
    └── ses-{sessionID}
        └── fmap
            └── sub-{subjectID}_ses-{sessionID}_TB1map.nii.gz
```

## Contributing

This repository uses [pre-commit](https://pre-commit.com/) to ensure that all committed code follows basic quality standards. Please install it and configure it to run as part of ``git commit`` by running ``pre-commit install`` in your local repository:

```shell
pipx install pre-commit
pre-commit install  # set up the hook (to run pre-commit during 'git commit')
```
