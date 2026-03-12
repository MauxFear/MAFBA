# Figure 6 - Flux Variability Analysis (FVA)

## Overview

Figure 6 presents Flux Variability Analysis (FVA) results comparing flux variability between:
- **MAFBA**: Parametrized model with both enzymatic and membrane constraints at multiple γ values
- **irreversible iML1515**: Model without enzymatic constraints (FBA only)

The comparison reveals how resource allocation constraints affect metabolic flexibility.

**Panels generated:**
- **Panel A**: All reactions - cumulative distribution
- **Panel B**: Membrane reactions only - cumulative distribution
- **Panel C**: Non-membrane reactions only - cumulative distribution
- **Panel D**: Top 5 subsystems - split-violin plot (γ=0.25 vs γ=1.00)

**Panels A-C**: Cumulative distributions on logarithmic scale  
**Panel D**: Split-violin distributions by subsystem

**Key findings:**
- Quantifies metabolic flexibility across the network
- Compares constrained vs unconstrained models
- Distinguishes variability patterns between membrane and non-membrane reactions
- Shows impact of membrane constraints (γ) on flux variability

## Script

### `Figure_6_FVA_cumulative.m` - FVA Cumulative Distribution Plots

**Purpose**: Generate cumulative distribution plots comparing flux variability between constrained and unconstrained models.

**What it does:**
- Loads FVA results from `../../data/outputs/FVA_Results/`
- Compares two models:
  - MAFBA: γ = 0.2, 0.25, 0.3, 1.0
  - irriML1515: γ = 1.0 (FBA only, no enzymatic constraints)
- Filters for uptake bound = 1000, FVA_Percent = 0
- For each analysis type (all, membrane, non-membrane):
  - Calculates cumulative distribution of FVA differences for each γ
  - Plots all conditions on same axes (logarithmic x-axis)
  - Color-codes by constraint level

**Outputs** (saved to `output/`):
- `iML1515_MAFBA_FVA_panel_A.svg/.png` (All reactions - cumulative)
- `iML1515_MAFBA_FVA_panel_B.svg/.png` (Membrane reactions - cumulative)
- `iML1515_MAFBA_FVA_panel_C.svg/.png` (Non-membrane reactions - cumulative)
- `iML1515_MAFBA_FVA_panel_D.svg/.png` (Top 5 subsystems - split-violin)

**Runtime**: ~2-3 minutes (large datasets: 519K + 324K rows + subsystem mapping)

## Data Sources

### Input Data

Located in: `../../data/outputs/FVA_Results/`

| File | Description | Size |
|------|-------------|------|
| `FVA_results_MAFBA.csv` | FVA results for MAFBA (parametrized model) | ~37 MB |
| `FVA_results_irriML1515.csv` | FVA results for irreversible iML1515 | ~24 MB |
| `iML1515_rxns_attrs.csv` | Reaction attributes (for Panel D subsystem mapping) | ~499 KB |

**Data structure:**
```
Model,FVA_Percent,Analysis_Type,Gamma_fMSA,ATPM,Uptake_Bound,Reaction,Original_Value,Min_Value,Max_Value,FVA_Difference
```

**Key columns:**
- `Analysis_Type`: all_reactions, membrane, nonmembrane
- `Gamma_fMSA`: Membrane surface area constraint parameter (γ)
- `Uptake_Bound`: Glucose uptake bound (mmol gDW⁻¹ h⁻¹)
- `FVA_Difference`: Max_Value - Min_Value (flux variability)

### Data Generation

FVA data can be regenerated using:
- **Script**: `2_Basic_Analysis/4_FVA/run_FVA.m`
- **Models**:
  - MAFBA: Both enzymatic and membrane constraints
  - irriML1515: No enzymatic constraints (standard FBA)
- **Method**: Flux Variability Analysis at optimal growth
- **Analysis types**:
  - `all_reactions`: All network reactions (~3,500)
  - `membrane`: Reactions localized to membranes (~600)
  - `nonmembrane`: Cytoplasmic and periplasmic reactions (~2,900)

## How to Run

### Requirements
- MATLAB (tested with R2020b or later)
- Data files in `../../data/outputs/FVA_Results/`

### Steps

1. **Navigate to this directory in MATLAB:**
   ```matlab
   cd Figures/Figure_6/
   ```

2. **Run the script:**
   ```matlab
   Figure_6_FVA_cumulative
   ```

3. **Check outputs:**
   ```matlab
   ls output/
   ```

## Understanding the Analysis

### Flux Variability Analysis (FVA)

FVA quantifies the **range of possible flux values** for each reaction while maintaining optimal growth:
- **Min Value**: Minimum flux consistent with optimal growth
- **Max Value**: Maximum flux consistent with optimal growth  
- **FVA Difference**: Max - Min (variability)

**High variability** indicates metabolic flexibility (multiple equivalent pathways).  
**Low variability** indicates metabolic rigidity (essential pathway).

### Analysis Types

| Type | Description | Reaction Count |
|------|-------------|----------------|
| **all_reactions** | All network reactions | ~3,500 |
| **membrane** | Inner/outer membrane reactions | ~600 |
| **nonmembrane** | Cytoplasmic + periplasmic | ~2,900 |

### Parameters

**Figure 6 configuration:**

| Parameter | Value | Meaning |
|-----------|-------|---------|
| Model | `MAFBA` | Parametrized model (both enzymatic + membrane constraints) |
| γ (Gamma_fMSA) | 0.2, 0.25, 0.3, 1.0 | Range of membrane allocation constraints |
| Uptake_Bound | 1000 | High glucose availability |
| FVA_Percent | 0 | Standard FVA (exact optimum) |

## Expected Results

- **Runtime**: 2-3 minutes (large datasets: ~843K total rows + subsystem mapping)
- **Figure outputs**: 4 panels × 2 formats = 8 files
  - Panels A-C: Cumulative distributions (all, membrane, non-membrane)
  - Panel D: Split-violin plot (top 5 subsystems)
  - Both SVG and PNG for each

### Panels A-C (Cumulative Distributions)
- **X-axis**: Logarithmic scale (10⁻³ to 10³)
- **Y-axis**: Cumulative count (linear scale)
- **Plot features**:
  - 5 curves per panel (4 γ values for mpac-iML + 1 for FBA-only)
  - Color-coded by constraint level
  - Legend shows all conditions
- **Interpretation**: 
  - Steep curves = many reactions with low variability
  - Long tails = few reactions with high variability
  - Tighter γ constraints → reduced variability (less flexibility)
  - FBA-only (irriML1515) shows highest variability
  - Membrane vs non-membrane differences reveal compartment-specific constraints

### Panel D (Split-Violin by Subsystem)
- **X-axis**: Flux variability (linear scale, 0 to max)
- **Y-axis**: Top 5 subsystems (ranked by absolute difference)
- **Plot features**:
  - Upper half (orange): γ = 0.25 (both constraints)
  - Lower half (blue): γ = 1.00 (no membrane constraint)
  - Kernel density estimation + scatter points + median lines
- **Interpretation**:
  - Shows which metabolic subsystems are most affected by membrane constraints
  - Wider distributions = more variability (flexibility)
  - Differences between halves reveal constraint impact per subsystem

## Output Naming Convention

Format: `iML1515_MAFBA_FVA_panel_{A|B|C}.svg/.png`

Examples:
- `iML1515_MAFBA_FVA_panel_A.svg` → All reactions
- `iML1515_MAFBA_FVA_panel_B.svg` → Membrane reactions
- `iML1515_MAFBA_FVA_panel_C.svg` → Non-membrane reactions

## Troubleshooting

**"Data file not found"**
- Check that FVA data exists in `../../data/outputs/FVA_Results/`
- Run `ls ../../data/outputs/FVA_Results/` to verify

**"No data found for the specified filters"**
- Script filters for:
  - Models = 'MAFBA', 'irriML1515'
  - Uptake_Bound = 1000
  - FVA_Percent = 0
  - Gamma_fMSA: MAFBA (0.2, 0.25, 0.3, 1.0), irriML1515 (1.0)
- Verify these values exist in your FVA results

**Script runs slowly**
- FVA datasets are large (~61 MB total, 843K rows)
  - MAFBA: 37 MB, 519K rows
  - irriML1515: 24 MB, 324K rows
- Expected runtime: 2-3 minutes
- Use `fprintf` messages to track progress

**Figures don't match publication**
- Check that correct models ('MAFBA', 'irriML1515') are selected
- Verify uptake = 1000 and FVA_Percent = 0
- Ensure all gamma values are present in data
- Check analysis types are correctly filtered

## Comparison with Other Figures

**Figure 3 (Flux predictions)**:
- Shows *predicted flux values* at different growth rates
- Deterministic predictions from FBA/MAFBA
- Single model, multiple growth conditions

**Figure 6 (Flux variability)**:
- Shows *ranges of possible flux values* at optimal growth
- Quantifies metabolic flexibility (FVA)
- Compares constrained (mpac-iML) vs unconstrained (FBA-only) models
- Multiple constraint levels (γ values)

Key differences:
- Figure 3: "What fluxes occur at each growth rate?" (deterministic)
- Figure 6: "How flexible are fluxes at optimal growth?" (variability)
- Figure 6 also compares resource allocation constraints vs no constraints

## Directory Structure

```
Figure_6/
├── README.md                               (this file)
├── Figure_6_FVA_cumulative.m               (plotting script)
└── output/                                 (generated figures)
    ├── iML1515_MAFBA_FVA_panel_A.svg/.png  (cumulative - all)
    ├── iML1515_MAFBA_FVA_panel_B.svg/.png  (cumulative - membrane)
    ├── iML1515_MAFBA_FVA_panel_C.svg/.png  (cumulative - non-membrane)
    └── iML1515_MAFBA_FVA_panel_D.svg/.png  (split-violin - subsystems)
```

## Data Dependencies

```
../../data/outputs/FVA_Results/
├── FVA_results_MAFBA.csv       (~37 MB, 519K rows)
├── FVA_results_irriML1515.csv  (~24 MB, 324K rows)
└── iML1515_rxns_attrs.csv      (~499 KB, reaction→subsystem mapping)
```

## Related Scripts

- **FVA generation**: `2_Basic_Analysis/4_FVA/run_FVA.m`
- **Supplementary version (homogeneous model)**: `Supplementary_Figures/Figure_S9/Figure_S9_FVA_cumulative.m`

## Notes

- Compares two model types: constrained (mpac-iML) vs unconstrained (FBA-only)
- High uptake bound (1000) ensures glucose is not limiting
- Multiple γ values for mpac-iML show impact of membrane constraints
- irriML1515 (FBA-only) serves as reference for maximum metabolic flexibility
- Analysis focuses on metabolic network flexibility under different resource allocation constraints
