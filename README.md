# MAFBA

Code, data, and figure scripts for Membrane-Allocation Flux Balance Analysis (MAFBA).

## Layout

- `1_Model_Construction/`: build and export model files
- `2_Basic_Analysis/`: glucose sweeps, sectors, sensitivity, and FVA
- `Figures/`: main-figure scripts
- `Supplementary_Figures/`: supplementary-figure scripts
- `data/`: inputs and staged outputs
- `utils/`: shared helpers

## Requirements

- MATLAB with the COBRA Toolbox
- Python with `cobra`, `pandas`, `numpy`, `matplotlib`, `scipy`, `openpyxl`, `pyyaml`, and `tqdm`
- Gurobi for the main optimization workflow

## Entry Points

- Rebuild model assets from `1_Model_Construction/README.md`
- Regenerate staged outputs from `2_Basic_Analysis/README.md`
- Recreate figures from the scripts under `Figures/` and `Supplementary_Figures/`

## Main Model Input

- `data/input/mafba_iML1515_ecV1_g1_MAFBA.xml`
