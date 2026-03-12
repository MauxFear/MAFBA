# Supplementary Figure S2

**Flux simulations using the heterogeneous approach.**

## Description

Average metabolic fluxes vs growth rate for γ = 0.25, derived from 1000 simulations
in which the protein cost parameter *w* was drawn independently at random for each
reaction (heterogeneous approach). Key metabolic fluxes are shown together with the
acetate threshold (λ_ac).

## Files

| File | Description |
|------|-------------|
| `Figure_S2_heterogeneous.m` | MATLAB plotting script — loads pre-computed data for γ = 0.25 and generates Figure S2 |
| `Heterogeneous_plot.m` | MATLAB script that generates plots for all available γ values |
| `output/` | Generated figure files |

## Data

Pre-computed averaged flux tables (MATLAB `.mat` format, ~12 KB each):

| File | Location | γ |
|------|----------|---|
| `wc_iML_rand_table_g_025.mat` | `data/outputs/Heterogeneous/` | 0.25 (primary — Figure S2) |
| `wc_iML_rand_table_g_1.mat` | `data/outputs/Heterogeneous/` | 1.00 |
| `wc_iML_rand_table_g_023.mat` | `data/outputs/Heterogeneous/` | 0.23 |
| `wc_iML_rand_table_g_02.mat` | `data/outputs/Heterogeneous/` | 0.20 |
| `wc_iML_rand_table_g_015.mat` | `data/outputs/Heterogeneous/` | 0.15 |

Each `.mat` file contains a MATLAB `table` variable named `T` with columns:
`Var1` (wCc sweep values), `gr`, `ac`, `glc`, `memb`, `akgdh`, `mals`, `edd`, `co2`,
`phiCm`, `phiCc`, `phiR`, `phiEc`, `phiEr`, `phiEm`, `phiM`, `phi_max`.

> These files contain averages over 1000 LP solutions per growth-rate point,
> generated using random protein cost assignments to represent cellular heterogeneity.

## How to Run

Open MATLAB, navigate to this folder, and run:

```matlab
Figure_S2_heterogeneous
```

Output: `output/Figure_S2_heterogeneous_g025.png`

To generate plots for all gamma values:

```matlab
Heterogeneous_plot
```

Output: `output/heterogeneous_flux_g*.png` (one file per γ)

## Simulation Details

| Parameter | Value |
|-----------|-------|
| γ (membrane ratio) | 0.25 |
| Simulations per point | 1000 |
| Growth-rate points | 100 |
| Weight distribution | Exponential box (`expbox`), mean *ŵ* = 0.000822 |
| Variable parameter | *w_Cc* (chaperone group weight, swept 0 → 1) |
| Carbon source | Glucose (iML1515 model) |
