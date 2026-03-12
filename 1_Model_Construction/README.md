# Model Construction

This directory contains the model build pipeline.

## Step 1

`01_Python_Annotations/Model_Annotations.py`

- reads the base SBML and protein annotation inputs
- writes reaction-attribute tables and annotated SBML outputs under `03_MAFBA_Models_Outputs/`

## Step 2

`02_MATLAB_Build/example_build_script.m`

- loads the annotation outputs
- prepares MATLAB model files
- writes model variants and attribute spreadsheets under `03_MAFBA_Models_Outputs/`

## Key Files

- `01_Python_Annotations/README.md`
- `02_MATLAB_Build/example_build_script.m`
- `02_MATLAB_Build/loadEcModel.m`
- `02_MATLAB_Build/generateModelExcel.m`

## Model Codes

- `iML1515_MAFBA`
- `iML1515_Homogeneous`
