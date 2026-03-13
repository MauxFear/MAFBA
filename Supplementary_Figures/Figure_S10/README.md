# Supplementary Figure S10

**Correlation between acetate threshold and normalized heterologous protein level.**

## Description

Filled contour plots showing how the acetate threshold (λ_ac) depends on the
level of a heterologous protein and the membrane constraint parameter γ.
White space indicates an infeasible region (no growth solution).

| Panel | Protein type | Normalisation | Key observation |
|-------|-------------|---------------|-----------------|
| A | Cytosolic (φ_Z) | φ_P_max | For γ > 0.25, λ_ac depends only on protein level; below γ = 0.25, it correlates with both γ and protein level |
| B | Membrane (φ_Y) | φ_P_max | λ_ac is always correlated with both γ and protein level |
| C | Membrane (φ_Y) | φ_M_max (= γ · φ_P_max) | Nonlinear effects; normalising by the membrane capacity reveals a γ-independent region |
| D | Cytosolic (φ_Z) | φ_P_max | Annotated version of Panel A showing the two linear-correlation regions |

## Files

| File | Description |
|------|-------------|
| `Figure_S10_contour.m` | MATLAB script that generates all four panels |
| `Figure_S10_contour_reference.m` | Reference contour script using MAT-file inputs |
| `output/` | Generated SVG figures |

## Data

| File | Location | Description |
|------|----------|-------------|
| `aceL_protein_expression_combined.csv` | `data/outputs/Protein_Expression/` | Pre-computed acetate thresholds for heterologous protein sweeps (81 KB) |

### Data columns

| Column | Description |
|--------|-------------|
| `Gamma` | Membrane constraint parameter γ |
| `Z_Level` | Cytosolic heterologous protein allocation (φ_Z) |
| `Y_Level` | Membrane heterologous protein allocation (φ_Y) |
| `Acetate_threshold` | Acetate threshold λ_ac (h⁻¹) |
| `Change_Factor` | kcat scaling factor (1 = baseline) |

## How to Run

Open MATLAB, navigate to this folder, and run:

```matlab
Figure_S10_contour
```

Output SVGs are written to `output/`:
- `panel_A_cytosolic_phiPmax.svg`
- `panel_B_membrane_phiPmax.svg`
- `panel_C_membrane_phiMmax.svg`
- `panel_D_annotated.svg`

## Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `phi_max` | 0.484 | Total proteome capacity (φ_P_max, g/g_DW) |
| `Change_Factor` filter | 1 | Baseline kcat values only |
