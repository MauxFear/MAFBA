# Basic Analysis

This directory contains the main analysis workflows for the MAFBA model.

## Contents

- **`1_Glucose_Sweep/`**: Parameterized and homogeneous glucose-sweep entry points.
- **`2_Sectors/`**: Derives the homogeneous sector bundle in `data/outputs/Sectors/` from the staged glucose-sweep outputs.
- **`3_Sensitivity/`**: Regenerates the parameterized sensitivity dataset used by Figures 4 and 5 and stages outputs in `data/outputs/iML1515_Sensitivity/`.
- **`4_FVA/`**: Regenerates the staged MAFBA FVA dataset used by Figure 6 and writes outputs in `data/outputs/FVA_Results/`.

Helper functions shared across the analysis scripts (`load_mafba_model`, `ensure_dir`, `get_default_uptake_grid`, `detect_acetate_threshold`) are located in `utils/`.

## Usage

Run the workflows in order if you are rebuilding the staged outputs:

1. `1_Glucose_Sweep/run_glucose_sweep_homogeneous.m` for Figure 3 outputs, or `1_Glucose_Sweep/run_glucose_sweep_parametrized.m` for the default parameterized sweep
2. `2_Sectors/run_sectors.m`
3. `3_Sensitivity/run_sensitivity.m`
4. `4_FVA/run_FVA.m`

All four workflows load inputs from `data/input/`, write staged datasets under `data/outputs/`, and keep timestamped raw intermediates under each analysis folder's `output/raw/`.

Ensure the COBRA Toolbox is initialized in MATLAB before running these scripts.
