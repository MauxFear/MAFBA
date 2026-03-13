# FVA

Run `run_FVA.m` to regenerate the bundled MAFBA FVA dataset used by Figure 6.

Inputs:
- `../../data/input/mafba_iML1515_ecV1_g1_MAFBA.xml`
- `../../data/input/attributes_model_ecModel_iML1515_MAFBA_g_1_v3.xlsx`

Outputs:
- Raw timestamped run: `output/raw/fva/<timestamp>/`
- Staged files:
  - `../../data/outputs/FVA_Results/FVA_results_MAFBA.csv`
  - `../../data/outputs/FVA_Results/iML1515_rxns_attrs.csv`

Default settings:
- Gamma values: `1.0`, `0.3`, `0.25`, `0.2`
- Uptake bound: `1000`
- Stored `FVA_Percent`: `0` for the exact-optimum workflow used by Figure 6
