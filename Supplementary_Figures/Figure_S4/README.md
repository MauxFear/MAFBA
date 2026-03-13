# Supplementary Figure S4

**Manuscript file:** `Revised_Submission_BJ_MAFBA/Supplementary_Information/figures/Figure_S4.pdf`

## Caption (from manuscript)

**Sensitivity of the acetate secretion threshold to the protein cost of transporter reactions.** A sensitivity analysis was conducted to evaluate the effect of varying the proteome cost parameter (w), scaled by a change factor, and to estimate the impact on the acetate threshold (λac) at different γ values. A reduction in the proteome cost coefficients increases the acetate threshold up to ≈ 0.78. On the other hand, an increase in the proteome cost coefficients significantly reduces the acetate threshold.

## Plot type

Line or contour: **λac vs γ** for different change factors — **transporter reactions only** (not all reactions).

## Scripts that can help recreate

- **Data:** `data/outputs/Acetate_Threshold_Transport_Sensitivity/acetate_threshold_transport_sensitivity.csv`
- **Visualization:** `Figure_S4_transport_acetate_threshold_sensitivity.m`

## Data

- `data/outputs/Acetate_Threshold_Transport_Sensitivity/acetate_threshold_transport_sensitivity.csv`

This table contains the full `Gamma x Change_Factor` sweep and the threshold
data needed for the main and contour plots.

## Status

Run
`Figure_S4_transport_acetate_threshold_sensitivity.m` from this folder; outputs
are saved to `output/`.

The script generates:
- `acetate_threshold_main_vs_param.*`
- `acetate_threshold_contour_plot.*`
