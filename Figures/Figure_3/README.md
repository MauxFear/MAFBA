# Figure 3: Metabolic Flux Predictions vs Growth Rate

This directory contains scripts to generate Figure 3, which demonstrates metabolic regime transitions (respiratory → mixed-acid fermentation) as glucose uptake increases under different membrane proteome allocation constraints.

## Overview

Figure 3 illustrates how tighter membrane constraints (lower γ values) shift metabolic regimes and the emergence of acetate overflow metabolism above a critical growth rate threshold (λ_ac).

**Key findings:**
- Acetate overflow emerges above critical growth rate
- Membrane constraints shift metabolic regimes
- TCA cycle, glyoxylate shunt, and ED pathway respond to constraints

## Scripts

### 1. `v1_sensAll_DataCollection_Homogeneous.m` - Data Generation

**Purpose**: Glucose uptake sweep to generate flux data across gamma values.

**What it does:**
- Loads the bundled `iML1515_MAFBA` model
- For each γ ∈ {1.0, 0.5, 0.4, 0.35, 0.3, 0.25, 0.23, 0.22, 0.2, 0.15, 0.1}:
  - Sets membrane constraint (γ × total proteome)
  - Sweeps glucose uptake: 0.2–35 mmol/gDW/h (100 points)
  - Records reaction fluxes, proteome sectors, acetate threshold

**Direct script outputs** (saved under `simResults/MAFBA_homogeneous_all_reactions/<timestamp>/simData/`):
- `data_modelV0_all_reactions_MAFBA_homogeneous_final.csv` - All reaction fluxes
- `sectors_modelV0_all_reactions_MAFBA_homogeneous_final.csv` - Proteome sector allocations
- `aceL_modelV0_all_reactions_MAFBA_homogeneous_final.csv` - Acetate thresholds per γ
- `model_modelV0_all_reactions_MAFBA_homogeneous_final.csv` - Model metadata

**Expected staged filenames** (copy or rename into `../../data/outputs/iML1515_Homogeneous/`):
- `dataFlux_Homogeneous.csv`
- `sectors_Homogeneous.csv`
- `aceT_Homogeneous.csv`
- `model_Homogeneous.csv`

**Runtime**: ~5-15 minutes

### 2. `Figure_3_fluxes_fromdata.m` - Figure Generation

**Purpose**: Creates publication-quality plots from CSV data.

**What it does:**
- Loads flux data from `../../data/outputs/iML1515_Homogeneous/`
- For each γ value:
  - Detects acetate threshold using respiration slope analysis
  - Plots key metabolic fluxes vs growth rate
  - Marks λ_ac with vertical dashed line

**Outputs** (saved to `output/`):
- `Homogeneous_flux_g100.svg/.png` (γ = 1.00)
- `Homogeneous_flux_g050.svg/.png` (γ = 0.50)
- ... (one pair per γ value)

## How to Run

### Prerequisites
- MATLAB with COBRA Toolbox
- CPLEX or Gurobi solver
- Model files in `../../1_Model_Construction/03_MAFBA_Models_Outputs/iML1515_Homogeneous/`
- Helper functions from `../../1_Model_Construction/02_MATLAB_Build/`

### Step-by-Step

1. **Navigate to this directory:**
   ```matlab
   cd Figures/Figure_3/
   ```

2. **Generate raw homogeneous sweep outputs**:
   ```matlab
   v1_sensAll_DataCollection_Homogeneous
   ```
   This creates timestamped CSV files under:
   ```text
   simResults/MAFBA_homogeneous_all_reactions/<timestamp>/simData/
   ```

3. **Stage the final files for plotting** by copying or renaming the `*_final.csv` outputs into:
   ```text
   ../../data/outputs/iML1515_Homogeneous/
   ```
   using this mapping:
   ```text
   data_modelV0_all_reactions_MAFBA_homogeneous_final.csv    -> dataFlux_Homogeneous.csv
   sectors_modelV0_all_reactions_MAFBA_homogeneous_final.csv -> sectors_Homogeneous.csv
   aceL_modelV0_all_reactions_MAFBA_homogeneous_final.csv    -> aceT_Homogeneous.csv
   model_modelV0_all_reactions_MAFBA_homogeneous_final.csv   -> model_Homogeneous.csv
   ```

4. **Generate figures:**
   ```matlab
   Figure_3_fluxes_fromdata
   ```
   Creates SVG/PNG files in `output/`

### Quick Start
If the staged files already exist in `../../data/outputs/iML1515_Homogeneous/`, skip steps 2-3:
```matlab
Figure_3_fluxes_fromdata
```

## Data Files

Located in `../../data/outputs/iML1515_Homogeneous/`:

| File | Description | Size |
|------|-------------|------|
| `dataFlux_Homogeneous.csv` | All reaction fluxes | 1100 × ~2700 |
| `sectors_Homogeneous.csv` | Proteome sector allocations | 1100 × 17 |
| `aceT_Homogeneous.csv` | Acetate thresholds | 11 × 2 |
| `model_Homogeneous.csv` | Model metadata | Variable |
| `avg_*.csv` | Mean values (optional) | Variable |
| `stdv_*.csv` | Standard deviations (optional) | Variable |

**Data structure:**
- Columns: `Gamma`, `Uptake_Bound`, `Growth`, reaction IDs, sector allocations
- Rows: 11 γ values × 100 glucose uptake points = 1100 total

## Model

**iML1515_MAFBA**
- Base: iML1515 (E. coli genome-scale model)
- ~2,700 reactions (irreversible: `_f` forward, `_b` backward)
- ~1,800 metabolites
- Protein allocation constraints (total + membrane)
- Kinetic parameters from Alter Database v2

## Gamma (γ) Parameter

Fraction of total proteome allocated to membrane proteins:

| γ | Constraint | Meaning |
|---|------------|---------|
| 1.0 | No membrane limit | Baseline |
| 0.3 | 30% to membrane | Moderate limitation |
| 0.25 | 25% to membrane | Tighter constraint |
| 0.2 | 20% to membrane | Strong limitation |

## Key Reactions

| Reaction | Description | Color |
|----------|-------------|-------|
| `EX_glc__D_e_b` | Glucose uptake | Green |
| `EX_ac_e_f` | Acetate excretion | Red |
| `AKGDH` | TCA cycle (α-KG DH) | Purple |
| `ATPS4rpp_f` | ATP synthase | Blue |
| `MALS` | Glyoxylate shunt | Pink |
| `EDD` | ED pathway | Yellow |
| `EX_co2_e_f` | CO2 excretion | Gray |

## Expected Results

- **Data generation**: 5-15 minutes
- **Success rate**: >95% optimal solutions
- **Figure outputs**: 11 plots (22 files: SVG + PNG)
- **Acetate threshold**: Decreases with lower γ

## Troubleshooting

**"Model file not found"**
- Verify model exists in `../../1_Model_Construction/03_MAFBA_Models_Outputs/iML1515_Homogeneous/`

**"Solver not available"**
- Install CPLEX or Gurobi
- Change solver in script (line 13): `changeCobraSolver('cplex')`

**"Data file not found" (plotting)**
- Run `v1_sensAll_DataCollection_Homogeneous` first
- Copy the final `simResults/.../simData/*_final.csv` files into `../../data/outputs/iML1515_Homogeneous/` using the mapping above
- Check the staged files exist in `../../data/outputs/iML1515_Homogeneous/`

## Directory Structure

```
Figure_3/
├── README.md                                    (this file)
├── v1_sensAll_DataCollection_Homogeneous.m      (data generation)
├── Figure_3_fluxes_fromdata.m                   (plotting)
└── output/                                      (generated figures)
    ├── Homogeneous_flux_g100.svg/.png
    ├── Homogeneous_flux_g050.svg/.png
    └── ...
```

---

**Language**: MATLAB  
**Dependencies**: COBRA Toolbox, CPLEX/Gurobi, iML1515_MAFBA model  
**Status**: Complete ✓
