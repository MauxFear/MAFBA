# Supplementary Figure S1: Homogeneous Proteome Sector Allocation vs Growth Rate

This directory contains the homogeneous sector-allocation workflow for Supplementary Figure S1.

## Overview

Supplementary Figure S1 shows stacked area plots of proteome sector allocations for the homogeneous workflow, illustrating how cells redistribute protein resources between different functional categories as growth rate increases under different membrane constraints.

**Key findings:**
- Metabolic sectors (Cm, Em, Ec) dominate at high growth rates
- Ribosomal sector (R) increases with growth rate
- Membrane constraints alter sector allocations
- Acetate threshold (λ_ac) marks metabolic regime shift

## Scripts

### 1. `Figure_S1_sectors_fromdata.m` - MATLAB Implementation

**Purpose**: Generate stacked area plots of proteome sectors from CSV data.

**What it does:**
- Loads sectors and flux data from `../../data/outputs/Sectors/`
- For each γ value:
  - Detects acetate threshold (λ_ac)
  - Creates stacked area plot of proteome sectors vs growth rate
  - Adds hatching patterns for distinction
  - Marks protein limits (φP_max, φM_max)

**Outputs** (saved to `output/`):
- `Homogeneous_sectors_g100.svg/.png` (γ = 1.00)
- `Homogeneous_sectors_g100_hatch.png` (with hatching)
- ... (three files per γ value)

**Runtime**: < 1 minute

### 2. `Figure_S1_sectors_python.py` - Python Implementation (Alternative)

**Purpose**: Alternative Python implementation for sector visualization.

**What it does:**
- Same functionality as MATLAB version
- Uses matplotlib for plotting
- Outputs PNG and SVG formats

**Outputs** (saved to `output/`):
- `Homogeneous_sectors_g*.png/svg`

**Runtime**: < 1 minute

## How to Run

### Prerequisites
- **MATLAB version**:
  - MATLAB with plotting toolbox
  - `hatchfill2.m` utility (located in `../../utils/`)
  
- **Python version**:
  - Python 3.7+
  - pandas, numpy, matplotlib

### Step-by-Step

**Option A: MATLAB (recommended)**
1. Navigate to this directory in MATLAB:
   ```matlab
   cd Supplementary_Figures/Figure_S1/
   ```

2. Run the script:
   ```matlab
   Figure_S1_sectors_fromdata
   ```

**Option B: Python**
1. Navigate to this directory:
   ```bash
   cd Supplementary_Figures/Figure_S1/
   ```

2. Run the script:
   ```bash
   python Figure_S1_sectors_python.py
   ```

## Data Files

**Location**: `../../data/outputs/Sectors/`

| File | Description | Size |
|------|-------------|------|
| `sectors_modelV0_kcatscomb_Homog.csv` | Proteome sector allocations | ~133 KB |
| `data_modelV0_kcatscomb_Homog.csv` | Flux data (for acetate detection) | ~7.8 MB |

**Data source**: `data/outputs/Sectors/`

**Data structure:**
- Columns: `Gamma`, `Growth`, sector fractions (φCm, φCc, φR, φEc, φEr, φEm, φP, φM), limits (φPmax, φMmax)
- Rows: currently staged for γ = 1.0, 0.3, 0.25, and 0.2 (100 glucose uptake points each)

## Proteome Sectors

| Sector | Symbol | Description | Color |
|--------|--------|-------------|-------|
| Metabolic (cytoplasmic) | φCm | Cytoplasmic metabolic enzymes | Green |
| Metabolic (membrane) | φEm | Membrane metabolic enzymes | Yellow |
| Respiration | φEr | Respiratory chain enzymes | Blue |
| Ribosomes | φR | Translation machinery | Purple |
| Catabolism | φEc | Catabolic enzymes | Red |
| Core | φCc | Core cellular processes | Orange |

**Protein limits:**
- φP_max: Maximum proteome fraction (total)
- φM_max: Maximum membrane proteome (γ × φP_max)

## Gamma (γ) Parameter

Fraction of total proteome allocated to membrane proteins:

| γ | Constraint | Availability |
|---|------------|--------------|
| 1.0 | No membrane limit | Present |
| 0.3 | 30% to membrane | Present |
| 0.25 | 25% to membrane | Present |
| 0.2 | 20% to membrane | Present |

## Expected Results

- **Runtime**: < 1 minute
- **Figure outputs**: one stacked-area plot per staged γ value
- **Acetate threshold**: Visible as vertical dashed line (λ_ac)
- **Protein limits**: Horizontal dashed lines show φP_max and φM_max

## Troubleshooting

**"Data file not found"**
- Ensure the staged homogeneous sector bundle exists
- Check files exist in `../../data/outputs/Sectors/`
- If needed, regenerate them with `../../2_Basic_Analysis/2_Sectors/run_sectors.m`

**"hatchfill2 not found" (MATLAB)**
- Check that `hatchfill2.m` exists in `../../utils/`
- Script automatically adds utils to path

**MATLAB hatching not working**
- Uncomment line 20 to disable gamma filtering: `% gamma_values = [1];`
- Check that graphics renderer supports hatching

**Python module not found**
- Install required packages:
  ```bash
  pip install pandas numpy matplotlib
  ```

## Directory Structure

```
Figure_S1/
├── README.md                            (this file)
├── Figure_S1_sectors_fromdata.m         (MATLAB plotting)
├── Figure_S1_sectors_python.py          (Python plotting - alternative)
└── output/                              (generated figures)
    ├── Homogeneous_sectors_g100.svg/.png
    ├── Homogeneous_sectors_g100_hatch.png
    └── ...
```

---

**Language**: MATLAB (primary), Python (alternative)  
**Dependencies**: MATLAB plotting, hatchfill2.m (MATLAB); pandas, matplotlib (Python)  
**Data source**: Staged homogeneous outputs  
**Status**: Supplementary workflow preserved
