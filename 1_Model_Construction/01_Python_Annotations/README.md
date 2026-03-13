# Python Annotation (Step 1)

This step generates the annotated SBML model and reaction-attribute tables used by the MATLAB build.

## Main script
- `Model_Annotations.py`

## Inputs (from `data/input/`)
- Base model (SBML): `data/input/mafba_iML1515_ecV1_g1_MAFBA.xml`
- STEPdb processed data: `data/input/STEPdb/STEPdb_iML1515_prepared_corrected.csv`
- sMOMENT mapping: `data/input/iML1515_star_reactions_kcat_mapping_combined.json`
- Combined kcat DB: `data/input/combined_brenda_sabio_rk_iML_irr.json`
- Alter PAM data: `data/input/model_construction/PAM_alter_data.csv`

## Outputs (written under Step 3)
- `1_Model_Construction/03_MAFBA_Models_Outputs/<model_code>/models_InputData/`
  - `iML1515_irreversible_kcatDataBaseV<modelV>_<dateVersion>.csv`
  - `iML1515_irreversible_reactionDataV<modelV>_<dateVersion>.csv`
- `1_Model_Construction/03_MAFBA_Models_Outputs/<model_code>/model_OutputData/`
  - `annotated_iML1515_irreversible_V<modelV>_<dateVersion>.xml`

## How to run
1. From this folder, run:
   - `python Model_Annotations.py`
2. Edit the `modelV` / `dateVersion` values inside the script if needed.
3. `update_data` is optional. Leave it unset unless you specifically want to apply a manual curation sheet.
4. Outputs will be written to the Step 3 folder above.

## Notes
- This step expects the processed STEPdb CSV only.
