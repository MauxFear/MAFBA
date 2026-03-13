# 03_MAFBA_Models_Outputs

This folder stores model-generation outputs.

## Structure
- `<model_code>/models/g<gamma>/`
  - MAFBA models saved by gamma (e.g., `g1`, `g0.25`, `g0.23`, `g0.2`).
- `<model_code>/attributes_spreadsheets/`
  - Excel exports of model attributes. Some scripts write to per-script subfolders here.
- `<model_code>/models_InputData/`
  - Step 1 reaction-attribute CSVs and kcat tables.
- `<model_code>/model_OutputData/`
  - Step 1 annotated SBML outputs.

## Notes
- Existing files may be renamed with `_backup` suffixes before new outputs are written.
