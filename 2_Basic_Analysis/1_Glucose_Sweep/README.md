# Glucose Sweep

There are two entry points:
- `run_glucose_sweep_parametrized.m`: default sweep using `sigma_mode = 'parametrized_totalweight'`
- `run_glucose_sweep_homogeneous.m`: figure-compatible homogeneous sweep using `sigma_mode = 'homogeneous'`

The shared loader builds the working ec-model from the prepared model in `../../1_Model_Construction/02_MATLAB_Build/output/`. The returned model includes `protGroup`, allocation weights, and the irreversible exchange reactions expected by `get_pFBAsolution_ec_model`.

Both scripts follow the revision glucose-sweep logic used in A03/A04:
- canonical 40-point grid: 15 coarse points below 10 plus 25 dense points above 10
- glucose clamp on `EX_glc__D_e_b` with `lower_bound = 0`, `upper_bound = G`
- growth maximization followed by `Prot_Const` minimization
- `lambda_ac` defined at the first point where `EX_ac_e_f >= 0.01`
- incremental raw CSV/log writing during the sweep

Outputs:
- Parametrized raw run: `output/raw/glucose_sweep_parametrized/<timestamp>/`
- Homogeneous raw run: `output/raw/glucose_sweep_homogeneous/<timestamp>/`
  Each raw run contains:
  - `flux_data.csv`
  - `sectors_data.csv`
  - `aceL_data.csv`
  - `kcats_data.csv`
  - `coeffs_data.csv`
  - `run_metadata.json`
  - `logs/uptake_grid.log`
- Parametrized staged files:
  - `../../data/outputs/Parametrized/data_modelV0_kcatscomb_Parametrized.csv`
  - `../../data/outputs/Parametrized/sectors_modelV0_kcatscomb_Parametrized.csv`
  - `../../data/outputs/Parametrized/aceTData_modelV0_kcatscomb_Parametrized.csv`
- Homogeneous staged files:
  - `../../data/outputs/iML1515_Homogeneous/dataFlux_Homogeneous.csv`
  - `../../data/outputs/iML1515_Homogeneous/sectors_Homogeneous.csv`
  - `../../data/outputs/iML1515_Homogeneous/aceT_Homogeneous.csv`

Default gamma values: `1.0`, `0.3`, `0.25`, `0.2`

Use `run_glucose_sweep_homogeneous.m` for Figure 3 / Supplementary homogeneous outputs. Use `run_glucose_sweep_parametrized.m` when you want the parameterized basic-analysis sweep staged under `data/outputs/Parametrized/`.
