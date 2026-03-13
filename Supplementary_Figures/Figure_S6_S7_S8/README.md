# Supplementary Figures S6, S7, and S8

## Figure S6 – Acetate secretion profiles under translocation machinery costs

Acetate secretion profiles across γ values when the proteome cost of the
translocation machinery ranges from 5% to 30% of the ribosomal proteome
fraction.  The curves illustrate how increasing the translocation cost changes
the acetate secretion onset and profile as γ varies.  Vertical dashed lines
indicate the acetate threshold (λac).

**Output:** `output/Figure_S6.{png,pdf}`

## Figure S7 – Proteome sectors with added translocation machinery costs, γ = 0.30

Proteome sector allocation at γ = 0.30 for each export fraction.  Stacked area
plots show cumulative proteome sector fractions versus growth rate.  Reference
lines indicate φMmax (red dashed) and φPmax (black dashed); the vertical dotted
line marks λac.

**Output:** `output/Figure_S7.{png,pdf}`

## Figure S8 – Proteome sectors with added translocation machinery costs, γ = 0.25

Same layout as Figure S7 at γ = 0.25.  Incorporating the translocation
machinery costs tightens the membrane constraint and reduces the acetate
threshold at this γ value.

**Output:** `output/Figure_S8.{png,pdf}`

---

## Scripts

| File | Description |
|------|-------------|
| `run_figure_s6_s7_s8.py` | Main entry point – plot from existing data or re-run analysis |
| `figure_s6_s7_s8_plot.py` | Plotting functions for S6, S7, and S8 |
| `protein_export_sweep.py` | Protein export gamma sweep simulation class |

## Data

Pre-computed simulation results are in `data/outputs/Figure_S6_S7_S8/`:

| File | Contents |
|------|----------|
| `flux_data.csv` | Per-point metabolic flux solutions |
| `sectors_data_phiT.csv` | Per-point proteome sector fractions (with translocation sector) |
| `aceL_data.csv` | Detected acetate threshold (λac) for each (f, γ) combination |
| `analysis_metadata.json` | Parameters (export fractions, γ values, φmax) |

## Usage

**Plot from existing data (default):**

```bash
python run_figure_s6_s7_s8.py
```

**Re-run the simulation and plot:**

```bash
python run_figure_s6_s7_s8.py --run-analysis
```

The simulation requires `data/input/mafba_iML1515_ecV1_g1_MAFBA.xml` and
`data/input/attributes_model_ecModel_iML1515_MAFBA_g_1_v3.xlsx`.

## Parameters

| Parameter | Value |
|-----------|-------|
| Export fractions (f) | 0.05, 0.10, 0.20, 0.30 |
| γ values | 0.15, 0.20, 0.25, 0.30, 0.35, 1.00 |
| γ for Figure S7 | 0.30 |
| γ for Figure S8 | 0.25 |
| Acetate detection threshold | 0.01 mmol/gDW/h |
| Solver | Gurobi (falls back to GLPK) |
