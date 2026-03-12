# Supplementary Figure S5 (panels A and B)

**Manuscript files:** `Revised_Submission_BJ_MAFBA/Supplementary_Information/figures/Figure_S5a.pdf`, `Figure_S5b.pdf`

## Caption (from manuscript)

**Robustness of the acetate threshold to kcat uncertainty.** Acetate secretion threshold λac as a function of the membrane proteome ratio γ = φM,max/φP,max for the baseline parametrization (original kcat values) and for Monte Carlo perturbations of kcat with different coefficients of variation (CV). **(A)** Uncertainty applied to **all reactions**. **(B)** Uncertainty applied to **membrane-related reactions only**. Symbols show the mean acetate threshold over Monte Carlo samples and error bars denote one standard deviation.

## Plot type

- **S5a:** λac vs γ — baseline + CV 10%, 20%, 50%, 100% (all reactions).
- **S5b:** Same layout — membrane reactions only.

## Scripts

| Script | Purpose |
|--------|--------|
| `figure_s5_plot.py` | Reads kcat sampling run outputs (all-reactions and membrane-only), draws S5a and S5b (acetate threshold vs γ, mean ± SD), saves to `output/Figure_S5a.pdf` and `output/Figure_S5b.pdf`. |
| `run_figure_s5.py` | Runs the plot with default paths. Use `--run-analysis` to run kcat sampling for CV 10%, 20%, 50%, 100% first, then plot. |

## Data

All inputs and outputs use paths under the project `data/` folder only.

- **S5a (all reactions):** Run directories containing `analysis_metadata.json` (with `cv`) and `aceT_data.csv`. Default: `data/outputs/Figure_S5/all_reactions/` (timestamped subdirs).
- **S5b (membrane only):** Same structure. Default: `data/outputs/Figure_S5/membrane_only/`.
- **Inputs for sampling:** `data/input/mafba_iML1515_ecV1_g1_MAFBA.xml`, `data/input/attributes_model_ecModel_iML1515_MAFBA_g_1_v3.xlsx`.

## How to run

**Plot only (existing data):**

```bash
cd Supplementary_Figures/Figure_S5
python run_figure_s5.py
```

With custom paths (must be under the project):

```bash
python figure_s5_plot.py --all-reactions-dir ../../data/outputs/Figure_S5/all_reactions \
  --membrane-only-dir ../../data/outputs/Figure_S5/membrane_only \
  --output-dir output
```

**Run sampling for CV 10%, 20%, 50%, 100%, then plot:**

```bash
python run_figure_s5.py --run-analysis
```

Requires COBRApy, a solver (e.g. Gurobi), and the model at `data/input/mafba_iML1515_ecV1_g1_MAFBA.xml` and attributes at `data/input/attributes_model_ecModel_iML1515_MAFBA_g_1_v3.xlsx`.

## Outputs

- `output/Figure_S5a.pdf`, `output/Figure_S5a.png` — panel A (all reactions).
- `output/Figure_S5b.pdf`, `output/Figure_S5b.png` — panel B (membrane only).
