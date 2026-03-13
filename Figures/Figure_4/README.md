# Figure 4: Parameterized Proteome Sector Allocation vs Growth Rate

This directory contains scripts to generate Figure 4 using the parameterized baseline dataset staged in `data/outputs/Parametrized/`.

## Overview

Figure 4 shows stacked area plots of proteome sector allocations for the parameterized approach, illustrating how cells redistribute protein resources between different functional categories as growth rate increases under selected membrane constraints.

**Key findings:**
- Uses the same baseline sensitivity slice as Figure 5 (`Change_Factor = 1.0`)
- Compares proteome sectors for γ = 1.00, 0.30, and 0.25
- Preserves the stacked-area sector view while switching to the parameterized dataset
- Acetate threshold (λ_ac) marks the transition point in each panel

## Scripts

### 1. `Figure_4_sectors_fromdata.m` - MATLAB Implementation

**Purpose**: Generate stacked area plots of proteome sectors from the baseline sensitivity data.

**What it does:**
- Loads sectors and flux data from `../../data/outputs/Parametrized/`
- Uses `Change_Factor = 1.0` when that column exists; otherwise it uses the full dataset in the file
- For γ = 1.00, 0.30, and 0.25:
  - Detects acetate threshold (λ_ac)
  - Creates stacked area plots of proteome sectors vs growth rate
  - Adds hatching patterns for distinction
  - Marks protein limits (φP_max, φM_max)

**Outputs** (saved to `output/`):
- `iML1515_MAFBA_parametrized_sectors_g100_cf100.svg/.png`
- `iML1515_MAFBA_parametrized_sectors_g100_cf100_hatch.png`
- ... (three files per γ value)

**Runtime**: < 1 minute

### 2. `Figure_4_sectors_python.py` - Python Implementation (Alternative)

**Purpose**: Alternative Python implementation for the same sector visualization.

The Python version follows the same visual logic as the MATLAB version and uses tolerance-based filtering for floating-point `Gamma` values.

## How to Run

### Prerequisites

- **MATLAB version**
  - MATLAB with plotting toolbox
  - `hatchfill2.m` utility in `../../utils/`
- **Python version**
  - Python 3.7+
  - pandas, numpy, matplotlib

### Step-by-Step

**Option A: MATLAB (recommended)**

1. Navigate to this directory in MATLAB:
   ```matlab
   cd Figures/Figure_4/
   ```

2. Run the script:
   ```matlab
   Figure_4_sectors_fromdata
   ```

**Option B: Python**

1. Navigate to this directory:
   ```bash
   cd Figures/Figure_4/
   ```

2. Run the script:
   ```bash
   python Figure_4_sectors_python.py
   ```

## Data Files

**Location**: `../../data/outputs/Parametrized/`

| File | Description | Size |
|------|-------------|------|
| `sectors_modelV0_kcatscomb_Parametrized.csv` | Proteome sector allocations | varies |
| `data_modelV0_kcatscomb_Parametrized.csv` | Flux data used for acetate-threshold detection | varies |
| `aceTData_modelV0_kcatscomb_Parametrized.csv` | Precomputed acetate thresholds (reference) | varies |

**Data source**: Parameterized baseline dataset (same staged dataset used by Figure 5)

**Data structure:**
- Columns include `Gamma`, `Growth`, sector fractions (`phiCm`, `phiCc`, `phiR`, `phiEc`, `phiEr`, `phiEm`), and limits (`phiPmax`, `phiMmax`)
- If a `Change_Factor` column is present it is filtered to `1.0`; otherwise the full CSV is used as-is

## Figure 4 Parameters

| γ | Constraint | Included |
|---|------------|----------|
| 1.00 | No membrane constraint | Yes |
| 0.30 | Moderate membrane limitation | Yes |
| 0.25 | Tighter membrane limitation | Yes |

**Change factor:** `1.0` only (baseline parameterization)

## Expected Results

- **Runtime**: < 1 minute
- **Figure outputs**: 3 γ values × 3 formats = 9 files (MATLAB)
- **Figure outputs**: 3 γ values × 2 formats = 6 files (Python)
- **Acetate threshold**: visible as vertical dashed line (λ_ac)
- **Protein limits**: horizontal dashed lines show φP_max and φM_max

## Troubleshooting

**"Data file not found"**
- Check that `../../data/outputs/Parametrized/` exists
- Verify both `sectors_modelV0_kcatscomb_Parametrized.csv` and `data_modelV0_kcatscomb_Parametrized.csv` are present

**"hatchfill2 not found" (MATLAB)**
- Check that `hatchfill2.m` exists in `../../utils/`

**Python module not found**
```bash
pip install pandas numpy matplotlib
```

## Directory Structure

```text
Figure_4/
├── README.md
├── Figure_4_sectors_fromdata.m
├── Figure_4_sectors_python.py
└── output/
    ├── iML1515_MAFBA_parametrized_sectors_g100_cf100.svg/.png
    ├── iML1515_MAFBA_parametrized_sectors_g100_cf100_hatch.png
    └── ...
```

---

**Language**: MATLAB (primary), Python (alternative)  
**Dependencies**: MATLAB plotting + `hatchfill2.m`; pandas, numpy, matplotlib  
**Data source**: Baseline parameterized sensitivity analysis  
**Status**: Updated manuscript mapping
