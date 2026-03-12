# %%
from pathlib import Path
import importlib.util
import cobra

base_dir = Path(__file__).resolve().parent
repo_root = base_dir.parents[1]

module_path = repo_root / "1_Model_Construction" / "01_Python_Annotations" / "Model_Annotations.py"
module_spec = importlib.util.spec_from_file_location("Model_Annotations", module_path)
Model_Annotations = importlib.util.module_from_spec(module_spec)
module_spec.loader.exec_module(Model_Annotations)

model_id = "iML1515"
modelV = 1
dateVersion = "MAFBA"
model_code = f"{model_id}_{dateVersion}"

data_input = repo_root / "data" / "input"
output_root = repo_root / "1_Model_Construction" / "03_MAFBA_Models_Outputs" / model_code
(output_root / "models_InputData").mkdir(parents=True, exist_ok=True)
(output_root / "model_OutputData").mkdir(parents=True, exist_ok=True)

model_path = data_input / f"mafba_{model_id}_ecV{modelV}_g1_{dateVersion}.xml"
update_data = None  # optional; set to the update Excel path before running if needed

# Inputs and outputs
paths = {
    "model_path": str(model_path),
    "sMOMENT_reactions_kcat_map": str(data_input / "iML1515_star_reactions_kcat_mapping_combined.json"),
    "alter_data": str(data_input / "model_construction" / "PAM_alter_data.csv"),
    "STEPdb_data": str(data_input / "STEPdb" / "STEPdb_iML1515_prepared_corrected.csv"),
    "combined_kcatdb": str(data_input / "combined_brenda_sabio_rk_iML_irr.json"),
    "input_folder": str(data_input) + "/",
    "output_folder": str(output_root) + "/",
    "update_data": str(update_data) if update_data else None,
}

suffixes_PAM = ["f", "b"]
optionIrreversible = False
additionalKcats = False
optionDefaultKcat = "mean"

annotated_model = Model_Annotations.Annotate_Model(
    paths,
    suffixes_PAM,
    modelV,
    dateVersion,
    False,
    optionIrreversible,
    additionalKcats,
    optionDefaultKcat,
)

annotated_model.save_model_as_xml()
