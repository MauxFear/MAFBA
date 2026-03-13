# Figure 8 - Normalized Protein Expression Analysis

## Overview

Figure 8 analyzes the effect of heterologous protein expression on metabolic regime transitions in *E. coli* by examining how cytosolic (φ_Y) and membrane (φ_Z) protein allocations affect the acetate overflow threshold across different membrane surface area constraints (γ values).

Key insight: By normalizing protein allocations relative to maximum capacity (φ_max) or membrane-constrained capacity (φ_Mmax), the analysis reveals fundamental trade-offs between protein allocation and metabolic efficiency.

## Figure Panels

### Panel B: Membrane Protein Expression (φ_Z / φ_max^P)
Shows how allocating protein budget to membrane proteins affects the acetate threshold. Membrane proteins increase respiratory capacity but compete for total protein budget.

### Panel C: Cytosolic Protein Expression (φ_Y / φ_max^P)
Shows how allocating protein budget to cytosolic proteins affects the acetate threshold. Cytosolic proteins increase metabolic capacity but compete for total protein budget.

### Panel D: Membrane-Normalized Cytosolic Expression (φ_Y / φ_max^M)
Shows cytosolic protein allocation normalized by membrane-constrained maximum (φ_max^M = γ × φ_max). This reveals how membrane constraints modulate the impact of cytosolic protein expression.

### Combined Plot (reference)
Overlays both Y and Z protein effects for direct comparison across all gamma values. This comprehensive view is provided for analysis but does not appear in the main figure.

## Running the Analysis

### Prerequisites
- MATLAB R2020b or later
- Filtered dataset: `../../data/outputs/Protein_Expression/data_protein_expression_figure8.csv` (~605 MB)

### Execution
```matlab
cd Figures/Figure_8/
Figure_8_protein_normalized
```

Runtime: ~30 seconds

### Output
Four figures generated in `output/` directory (SVG and PNG formats):
- `iML1515_MAFBA_protNorm_panel_B.svg/.png` - Membrane protein (Panel B)
- `iML1515_MAFBA_protNorm_panel_C.svg/.png` - Cytosolic protein (Panel C)
- `iML1515_MAFBA_protNorm_panel_D.svg/.png` - Membrane-normalized cytosolic (Panel D)
- `iML1515_MAFBA_protNorm_combined.svg/.png` - Combined reference plot

## Data Source

**Input**: Filtered protein expression sensitivity analysis data containing systematic variations of:
- **Membrane constraints**: γ = 1.0, 0.35, 0.3, 0.25, 0.22, 0.2
- **Protein allocations**: All Y_Level (cytosolic) and Z_Level (membrane) combinations
- **Reaction fluxes**: Complete flux profiles for acetate threshold detection

**Filtering**: Only gamma values are filtered; all protein allocation levels and reaction columns are retained to ensure compatibility with the acetate threshold detection algorithm.

**Size**: 41,400 rows × 3,697 columns (~605 MB)

### Regenerating Filtered Data
If the filtered dataset needs to be regenerated, filter the full protein expression sweep output
to include only the required gamma values (`1.0, 0.35, 0.3, 0.25, 0.22, 0.2`) and save
to `data/outputs/Protein_Expression/data_protein_expression_figure8.csv`.

## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| φ_max (φ_max^P) | 0.484 | Maximum total protein allocation (48.4% of cell mass) |
| φ_max^M | γ × φ_max | Membrane-constrained maximum protein allocation |
| λ_ac | Variable | Acetate overflow threshold (growth rate where overflow begins) |

## Methodology

### Normalization Approaches

1. **φ_max^P normalization** (Panels B, C): Protein allocation divided by theoretical maximum (0.484)
   - Shows fraction of total cellular protein capacity allocated
   - Independent of membrane constraints

2. **φ_max^M normalization** (Panel D): Protein allocation divided by membrane-constrained maximum (γ × 0.484)
   - Shows fraction of membrane-available protein capacity
   - Reveals membrane-limited regime behavior

### Acetate Threshold Detection

The acetate overflow threshold (λ_ac) is detected by identifying when membrane flux (ATP synthase) saturates, triggering the transition from pure respiration to mixed respiration + fermentation metabolism. The detection algorithm:

1. Calculates slope of membrane flux vs. growth rate
2. Identifies saturation point (slope drops below 1% of maximum)
3. Finds first point where acetate flux exceeds threshold (0.01 mmol/gDW/h)
4. Returns growth rate at this transition point

This method requires complete flux profiles for accurate detection.

## Scientific Context

**Biological question**: How does expressing heterologous proteins (Y for cytosolic, Z for membrane) affect cellular metabolic regime transitions?

**Key findings**:
- Cytosolic protein expression generally decreases acetate threshold (earlier overflow)
- Membrane protein expression has complex effects depending on membrane constraints
- Membrane-normalized view (Panel D) reveals that γ values collapse curves differently than absolute normalization
- Trade-offs between protein allocation and metabolic efficiency are γ-dependent

## File Structure

```
Figure_8/
├── Figure_8_protein_normalized.m    # Main analysis script
├── README.md                         # This file
└── output/                           # Generated figures (8 files)
```

## Related Scripts

- **Protein expression sensitivity (1D sweep)**: `Figures/Figure_7/Figure_7_protein_expression.m`
- **Contour view**: `Supplementary_Figures/Figure_S10/Figure_S10_contour.m`

## Notes

- All columns from source data are retained to ensure algorithm compatibility
- Panel A is not used; panels are numbered B, C, D per figure specification
- Combined plot serves as comprehensive reference but is not part of main figure panels

## Citation

If using this analysis, cite:
- MAFBA framework paper (protein allocation modeling)
- φ_max value: Schmidt et al. (2016) - proteomic measurements
- iML1515 metabolic model: Monk et al. (2017)
