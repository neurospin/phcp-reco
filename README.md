# phcp-reco
Reconstruction scripts for post-mortem multi-quantitative fetal brain MRI at 11.7 teslas


## Contributing

This repository uses [pre-commit](https://pre-commit.com/) to ensure that all committed code follows basic quality standards. Please install it and configure it to run as part of ``git commit`` by running ``pre-commit install`` in your local repository:

```shell
pipx install pre-commit
pre-commit install  # set up the hook (to run pre-commit during 'git commit')
```

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

### Reconstructed FoV data

```
fov/derivatives/fov-reconstructed
└── sub-{subjectID}
    └── ses-{sessionID}
        └── fmap
            └── sub-{subjectID}_ses-{sessionID}_TB1map.nii.gz
```

Reconstruction of the FoV data from the rawdata goes though the following steps:

```shell
phcp-fov-afi-b1mapping fov/rawdata/sub-${sub}/ses-${ses}/fmap/sub-${sub}_ses-${ses}_TB1AFI.nii.gz fov/derivatives/fov-reconstructed/sub-${sub}/ses-${ses}/fmap/sub-${sub}_ses-${ses}_TB1map.nii.gz
```
