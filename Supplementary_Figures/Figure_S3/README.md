# Supplementary Figure S3

**Manuscript file:** `Revised_Submission_BJ_MAFBA/Supplementary_Information/figures/Figure_S3.pdf`

## Caption (from manuscript)

**Sensitivity analysis of the acetate secretion threshold to variations in proteome cost coefficients.** A sensitivity analysis was conducted to evaluate the effect of varying the proteome cost parameter (w), scaled by a change factor, and to estimate the impact on the acetate threshold (λac) at different γ values. Changes in the proteome cost coefficients significantly affect the impact on the acetate threshold when γ is reduced.

## Plot type

Line or contour: **λac (acetate threshold)** vs **γ** (or ratio φM,max/φP,max) for different **change factors** — proteome cost w scaled (all reactions).

## Scripts that can help recreate

- **Data:** `data/outputs/Acetate_Threshold_Sensitivity/acetate_threshold_sensitivity_all_reactions.csv`
- **Visualization:** `Figure_S3_acetate_threshold_sensitivity.m`

## Data

- `data/outputs/Acetate_Threshold_Sensitivity/acetate_threshold_sensitivity_all_reactions.csv`

This table contains the full `Gamma x Change_Factor` sweep and the
`Consensus_Threshold` column required for the line and contour plots.

## Status

Run
`Figure_S3_acetate_threshold_sensitivity.m` from this folder; outputs are
saved to `output/`.
