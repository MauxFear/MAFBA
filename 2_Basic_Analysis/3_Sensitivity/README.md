# Sensitivity

Run `run_sensitivity.m` to regenerate the parameterized baseline dataset used by Figures 4 and 5.

Inputs:
- `../../data/input/mafba_iML1515_ecV1_g1_MAFBA.xml`

Outputs:
- Raw timestamped run: `output/raw/sensitivity/<timestamp>/`
- Staged files:
  - `../../data/outputs/iML1515_Sensitivity/dataFlux_Sensitivity.csv`
  - `../../data/outputs/iML1515_Sensitivity/sectors_Sensitivity.csv`
  - `../../data/outputs/iML1515_Sensitivity/aceT_Sensitivity.csv`

Default gamma values: `1.0`, `0.3`, `0.25`, `0.2`, `0.15`

Default change factors: `1`, `2`, `5`, `10`, `20`, `50`, `0.5`, `0.2`, `0.1`, `0.05`, `0.02`
