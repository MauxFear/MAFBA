# Sectors

Run `run_sectors.m` after `../1_Glucose_Sweep/run_glucose_sweep_homogeneous.m` to perform the homogeneous sector analysis.

Inputs:
- `../../data/outputs/iML1515_Homogeneous/dataFlux_Homogeneous.csv`
- `../../data/outputs/iML1515_Homogeneous/sectors_Homogeneous.csv`

Staged outputs:
- `../../data/outputs/Sectors/data_modelV0_kcatscomb_Homog.csv`
- `../../data/outputs/Sectors/sectors_modelV0_kcatscomb_Homog.csv`
- `../../data/outputs/Sectors/aceTData_modelV0_kcatscomb_Homog.csv`

Figure outputs:
- `output/Homogeneous_sectors_g*.svg`
- `output/Homogeneous_sectors_g*.png`
- `output/Homogeneous_sectors_g*_hatch.png`

What it does:
- stages the homogeneous sector bundle under `data/outputs/Sectors/`
- detects the acetate threshold from the homogeneous flux data
- generates stacked sector-allocation plots using the same visual logic as Supplementary Figure S1
