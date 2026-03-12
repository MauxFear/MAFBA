# Figure 5: Sensitivity Analysis - Metabolic Flux Predictions

This directory contains the script to generate Figure 5, which demonstrates how metabolic fluxes respond to the staged parameterized baseline dataset across different gamma values.

## Overview

Figure 5 shows metabolic flux predictions at baseline conditions (CF = 1.0, no parameter perturbation) for key gamma values, providing a focused comparison of metabolic behavior under different membrane constraints.

**Panels generated:**
- γ = 1.00 (no membrane constraint)
- γ = 0.30 (moderate constraint)
- γ = 0.25 (tighter constraint)
- γ = 0.23 (tight constraint near acetate-transition regime)

**All at Change Factor = 1.0 (baseline parameters)**

**Key findings:**
- Metabolic regime transitions at different membrane constraints
- Acetate threshold shifts with membrane limitations
- Comparison of flux profiles across key gamma values

## Script

### `Figure_5_fluxes_fromdata.m` - Flux Plot Generation

**Purpose**: Generate flux plots showing sensitivity to parameter changes.

**What it does:**
- Loads parameterized baseline data from `../../data/outputs/Parametrized/`
- Filters for specific gamma values: 1.0, 0.3, 0.25, 0.23
- Filters for Change Factor = 1.0 (baseline, no perturbation)
- For each gamma value:
  - Detects acetate threshold (λ_ac)
  - Plots key metabolic fluxes vs growth rate
  - Marks acetate threshold with vertical line

**Outputs** (saved to `output/`):
- `iML1515_MAFBA_flux_g100_cf100.svg/.png` (γ = 1.00)
- `iML1515_MAFBA_flux_g030_cf100.svg/.png` (γ = 0.30)
- `iML1515_MAFBA_flux_g025_cf100.svg/.png` (γ = 0.25)
- `iML1515_MAFBA_flux_g023_cf100.svg/.png` (γ = 0.23)

**Runtime**: < 1 minute

## How to Run

### Prerequisites
- MATLAB with plotting toolbox
- Parameterized baseline data (already provided in `../../data/outputs/Parametrized/`)

### Step-by-Step

1. **Navigate to this directory in MATLAB:**
   ```matlab
   cd Figures/Figure_5/
   ```

2. **Run the script:**
   ```matlab
   Figure_5_fluxes_fromdata
   ```
   Creates SVG and PNG files in `output/`

**Quick start**: Data already provided, just run the script.

## Data Files

**Location**: `../../data/outputs/Parametrized/`

| File | Description | Size |
|------|-------------|------|
| `data_modelV0_kcatscomb_Parametrized.csv` | All reaction fluxes | varies |
| `sectors_modelV0_kcatscomb_Parametrized.csv` | Proteome sector allocations | varies |
| `aceTData_modelV0_kcatscomb_Parametrized.csv` | Acetate thresholds | varies |

**Data source**: Parameterized baseline dataset staged in `data/outputs/Parametrized/`

**Data structure:**
- Columns include `Gamma`, `Uptake_Bound`, and reaction IDs
- If a `Change_Factor` column is present the script uses `CF = 1.0`; otherwise it uses the full dataset directly

## Parameters

### Gamma (γ) - Membrane Constraint
Fraction of total proteome allocated to membrane proteins.

**Figure 5 panels:**

| γ | Constraint Level | Panel |
|---|------------------|-------|
| 1.00 | No membrane constraint | Baseline |
| 0.30 | Moderate constraint | Medium |
| 0.25 | Tighter constraint | Tight |
| 0.23 | Tighter constraint | Transition |

### Change Factor (CF) - Parameter Perturbation

**Figure 5 uses CF = 1.0 only** when that column exists (baseline, no parameter perturbation)

This focuses the analysis on membrane constraint effects without confounding parameter sensitivity.

## Key Reactions Tracked

| Reaction | Description | Plot Color |
|----------|-------------|------------|
| `EX_glc__D_e_b` | Glucose uptake | Green |
| `EX_ac_e_f` | Acetate excretion | Red |
| `AKGDH` | TCA cycle (α-KG DH) | Purple |
| `ATPS4rpp_f` | ATP synthase | Blue |
| `MALS` | Glyoxylate shunt | Pink |
| `EDD` | ED pathway | Yellow |
| `EX_co2_e_f` | CO2 excretion | Gray |

## Expected Results

- **Runtime**: < 1 minute
- **Figure outputs**: 4 gamma values × 2 formats = 8 files
  - γ = 1.00, 0.30, 0.25, 0.23
  - Both SVG and PNG for each
- **Acetate threshold**: Visible as vertical dashed line (λ_ac)
- **Comparison**: Shows how membrane constraints shift metabolic regimes

## Output Naming Convention

Format: `iML1515_MAFBA_flux_g{gamma}_cf{changefactor}.svg/.png`

Examples:
- `iML1515_MAFBA_flux_g100_cf100.svg` → γ = 1.00, CF = 1.00
- `iML1515_MAFBA_flux_g030_cf100.svg` → γ = 0.30, CF = 1.00
- `iML1515_MAFBA_flux_g025_cf100.svg` → γ = 0.25, CF = 1.00
- `iML1515_MAFBA_flux_g023_cf100.svg` → γ = 0.23, CF = 1.00

## Troubleshooting

**"Data file not found"**
- Check that parameterized data exists in `../../data/outputs/Parametrized/`
- Verify `data_modelV0_kcatscomb_Parametrized.csv` is present

**"No data for this combination"**
- Script filters for specific gamma values (1.0, 0.3, 0.25, 0.23)
- If a gamma value is missing from data, it will be skipped with a warning

**Script runs slowly**
- Large dataset (~896 MB) takes time to load
- Typically completes in < 1 minute after loading
- Only generates 4 plots (one per gamma value)

**Out of memory**
- Script closes figures after saving to conserve memory
- If still problematic, reduce number of combinations by filtering data

## Directory Structure

```
Figure_5/
├── README.md                               (this file)
├── Figure_5_fluxes_fromdata.m              (plotting script)
└── output/                                 (generated figures)
    ├── iML1515_MAFBA_flux_g100_cf100.svg/.png
    ├── iML1515_MAFBA_flux_g030_cf100.svg/.png
    └── iML1515_MAFBA_flux_g025_cf100.svg/.png
```

## Comparison with Figure 3

**Figure 3 (Homogeneous)**:
- Data generated specifically for homogeneous analysis
- All 11 gamma values shown
- Cleaner dataset from targeted simulation

**Figure 5 (Parameterized baseline)**:
- Uses the staged baseline parameterized dataset
- Selected gamma values: 1.0, 0.3, 0.25, 0.23
- Same metabolic behavior as Figure 3, different data source
- Larger dataset but filtered to 4 key conditions

---

**Language**: MATLAB  
**Dependencies**: MATLAB plotting toolbox  
**Data source**: Parameterized baseline data staged in `data/outputs/Parametrized/`  
**Status**: Complete ✓
