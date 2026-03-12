# Supplementary Figure S9 – FVA, Homogeneous Approach

## Overview

Figure S9 compares flux variability of the proteome-constrained model at different
γ values using the **homogeneous** cost approach.  The structure is identical to
Figure 6 of the main document, substituting the parametrized (MAFBA) model
with the homogeneous model.

**Panels:**
- **A** – Cumulative distribution of flux variability (all reactions)
- **B** – Cumulative distribution of flux variability (membrane-associated reactions)
- **C** – Cumulative distribution of flux variability (non-membrane reactions)
- **D** – Split-violin plots comparing γ = 0.25 vs γ = 1.0 for the top 5 subsystems

## Script

`Figure_S9_FVA_cumulative.m`

**Changes vs Figure 6:**
- Model `Homog` replaces `MAFBA`
- Legend labels read "Homogeneous approach, γ = …"
- Output prefix is `iML1515_Homog_FVA`

## Data

Located in `../../data/outputs/FVA_Results/`:

| File | Description |
|------|-------------|
| `FVA_results_Homog.csv` | FVA results for the homogeneous-cost model |
| `FVA_results_irriML1515.csv` | FVA results for irreversible iML1515 (FBA only, reference) |
| `iML1515_rxns_attrs.csv` | Reaction → subsystem mapping for Panel D (3682 reactions, 100% FVA coverage) |

**Data structure:**
```
Model, FVA_Percent, Analysis_Type, Gamma_fMSA, ATPM, Uptake_Bound,
Reaction, Original_Value, Min_Value, Max_Value, FVA_Difference
```

**Filters applied:** `Uptake_Bound = 1000`, `FVA_Percent = 0`

## How to Run

```matlab
cd Supplementary_Figures/Figure_S9/
Figure_S9_FVA_cumulative
```

## Output Files

Saved to `output/`:

| File | Panel |
|------|-------|
| `iML1515_Homog_FVA_panel_A.svg/.png` | All reactions – cumulative |
| `iML1515_Homog_FVA_panel_B.svg/.png` | Membrane reactions – cumulative |
| `iML1515_Homog_FVA_panel_C.svg/.png` | Non-membrane reactions – cumulative |
| `iML1515_Homog_FVA_panel_D.svg/.png` | Subsystem split-violin (γ=0.25 vs γ=1.0) |
