#!/usr/bin/env python3
"""
Supplementary Figure S5: Robustness of acetate threshold to kcat uncertainty.

Generates two panels:
- S5a: λac vs γ for baseline + Monte Carlo kcat (multiple CVs), all reactions.
- S5b: Same, membrane-related reactions only.

Input: run directories containing analysis_metadata.json (with 'cv') and aceT_data.csv.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

try:
    import cobra
except ImportError:
    cobra = None

# Style
sns.set_style("whitegrid")
plt.rcParams["font.size"] = 11
plt.rcParams["font.family"] = "Arial"
plt.rcParams["figure.dpi"] = 300
plt.rcParams["axes.linewidth"] = 1.2

# Error bars use 95% CI (1.96 * SD / sqrt(n))
USE_STD = False  # False = 95% CI; True = ±1 SD


def project_root() -> Path:
    """Project root (three levels up from this script)."""
    return Path(__file__).resolve().parent.parent.parent


def _ace_file(run_path: Path) -> Optional[Path]:
    """Return the acetate-threshold CSV in a run dir (aceT_data.csv or aceL_data.csv)."""
    for name in ("aceT_data.csv", "aceL_data.csv"):
        p = run_path / name
        if p.exists():
            return p
    return None


def discover_runs(base_dir: Path) -> Dict[float, Path]:
    """Discover timestamped run dirs and map cv (from metadata) -> run Path."""
    run_dirs = {}
    if not base_dir.exists():
        return run_dirs
    for item in base_dir.iterdir():
        if not item.is_dir():
            continue
        if "comparison" in item.name or "TEST" in item.name:
            continue
        meta_path = item / "analysis_metadata.json"
        if not meta_path.exists():
            continue
        try:
            with open(meta_path) as f:
                meta = json.load(f)
            cv = meta.get("cv")
            if cv is not None and _ace_file(item) is not None:
                run_dirs[cv] = item
        except Exception:
            continue
    return run_dirs


def load_and_combine_aceT(run_dirs: Dict[float, Path]) -> Optional[pd.DataFrame]:
    """Load acetate-threshold CSV from each run, add 'cv' column, concatenate."""
    if not run_dirs:
        return None
    dfs = []
    for cv, run_path in sorted(run_dirs.items()):
        ace_path = _ace_file(run_path)
        df = pd.read_csv(ace_path)
        df["cv"] = cv
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


def compute_baseline(
    gamma_values: List[float],
    model_path: Path,
    solver: str = "gurobi",
) -> Optional[pd.DataFrame]:
    """Compute baseline acetate threshold (original kcat) for given gammas."""
    if cobra is None:
        return None
    if not model_path.exists():
        return None
    try:
        model = cobra.io.read_sbml_model(str(model_path))
        try:
            model.solver = solver
        except Exception:
            pass
    except Exception:
        return None

    membrane_constraint_id = "Memb_Const"
    if membrane_constraint_id not in model.reactions:
        return None
    membrane_rxn = model.reactions.get_by_id(membrane_constraint_id)
    baseline_membrane_bound = membrane_rxn.upper_bound

    points_below_10 = np.linspace(0.2, 9.8, 15)
    points_above_10 = 10 + 25 * (np.linspace(0, 1, 25) ** 2)
    glucose_uptake_values = np.concatenate([points_below_10, points_above_10])

    glc_id = "EX_glc__D_e_b"
    ac_id = "EX_ac_e_f"
    if glc_id not in model.reactions:
        return None

    baseline_data = []
    for gamma in gamma_values:
        membrane_rxn.upper_bound = gamma * baseline_membrane_bound
        growth_rates = []
        acetate_fluxes = []
        for uptake in glucose_uptake_values:
            with model:
                model.reactions.get_by_id(glc_id).upper_bound = uptake
                try:
                    sol = model.optimize()
                    if sol.status == "optimal":
                        growth = sol.objective_value
                        acetate_flux = sol.fluxes.get(ac_id, 0.0)
                    else:
                        growth = 0.0
                        acetate_flux = 0.0
                except Exception:
                    growth = 0.0
                    acetate_flux = 0.0
                growth_rates.append(max(growth, 0))
                acetate_fluxes.append(max(acetate_flux, 0))
        acetate_threshold = 0.0
        for idx, ac in enumerate(acetate_fluxes):
            if ac > 1e-6:
                acetate_threshold = growth_rates[idx]
                break
        baseline_data.append({"gamma": gamma, "acetate_threshold": acetate_threshold})
    return pd.DataFrame(baseline_data)


def plot_acetate_vs_gamma_panel(
    combined_aceT: pd.DataFrame,
    baseline_df: Optional[pd.DataFrame],
    title: str,
    output_dir: Path,
    filename_stem: str,
) -> None:
    """
    Single-panel: acetate threshold vs gamma, baseline + CV series (mean ± SD).
    """
    # Exclude CV=5% if present
    cv_values = sorted(combined_aceT["cv"].unique())
    cv_values = [c for c in cv_values if c != 0.05]
    if not cv_values:
        cv_values = sorted(combined_aceT["cv"].unique())

    base_pastel = ["#6BA3D6", "#7FCD91", "#E57373", "#DDA0DD", "#F4C2C2"]
    colors = base_pastel[: len(cv_values)]

    fig, ax = plt.subplots(figsize=(10, 6))

    for idx, cv in enumerate(cv_values):
        cv_data = combined_aceT[combined_aceT["cv"] == cv]
        summary = (
            cv_data.groupby("gamma")["acetate_threshold"]
            .agg(["mean", "std", "count"])
            .reset_index()
        )
        if USE_STD:
            yerr = summary["std"].values
        else:
            yerr = 1.96 * summary["std"].values / np.sqrt(summary["count"].values)
        ax.errorbar(
            summary["gamma"],
            summary["mean"],
            yerr=yerr,
            marker="o",
            capsize=5,
            capthick=2,
            linewidth=2.5,
            markersize=7,
            color=colors[idx],
            label=f"CV={cv:.1%}",
            alpha=0.9,
        )

    if baseline_df is not None:
        ax.plot(
            baseline_df["gamma"],
            baseline_df["acetate_threshold"],
            marker="s",
            linestyle="--",
            linewidth=2.5,
            markersize=8,
            color="#FF9F5A",
            label="Baseline (Original kcat values)",
            zorder=10,
        )

    ax.set_xlabel(r"$\gamma$ ratio ($\phi M_{max}/\phi P_{max}$)", fontsize=13)
    ax.set_ylabel(r"Acetate threshold ($\lambda_{ac}$, 1/hr)", fontsize=13)
    ax.set_title(title, fontsize=13, pad=15)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=11, framealpha=0.95)
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        out_path = output_dir / f"{filename_stem}.{ext}"
        plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"  Saved {filename_stem}.pdf/.png")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate Supplementary Figure S5 (S5a: all reactions, S5b: membrane only)."
    )
    root = project_root()
    default_all = root / "data" / "outputs" / "Figure_S5" / "all_reactions"
    default_membrane = root / "data" / "outputs" / "Figure_S5" / "membrane_only"
    default_out = Path(__file__).resolve().parent / "output"

    parser.add_argument(
        "--all-reactions-dir",
        type=Path,
        default=default_all,
        dest="all_reactions_dir",
        help="Directory with run subdirs for all-reactions kcat sampling (each: analysis_metadata.json, aceT_data.csv)",
    )
    parser.add_argument(
        "--membrane-only-dir",
        type=Path,
        default=default_membrane,
        dest="membrane_only_dir",
        help="Directory with run subdirs for membrane-only kcat sampling",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=default_out,
        help="Output directory for Figure_S5a and Figure_S5b",
    )
    parser.add_argument(
        "--model-path",
        type=Path,
        default=None,
        help="Path to SBML model for baseline (default: data/input/mafba_iML1515_ecV1_g1_MAFBA.xml)",
    )
    parser.add_argument(
        "--no-baseline",
        action="store_true",
        help="Skip baseline computation and plot only Monte Carlo series",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    model_path = args.model_path or (root / "data" / "input" / "mafba_iML1515_ecV1_g1_MAFBA.xml")

    print("Figure S5: kcat uncertainty")
    print("=" * 60)

    # Panel A: all reactions
    all_runs = discover_runs(args.all_reactions_dir)
    if not all_runs:
        print(f"No all-reactions data found under {args.all_reactions_dir}")
        print("  Provide run subdirs with analysis_metadata.json and aceT_data.csv (multiple CVs, e.g. 0.1, 0.2, 0.5, 1.0).")
    else:
        print(f"All reactions: {len(all_runs)} runs (CVs: {sorted(all_runs.keys())})")
        combined_all = load_and_combine_aceT(all_runs)
        gamma_values = sorted(combined_all["gamma"].unique())
        baseline_df = None
        if not args.no_baseline:
            baseline_df = compute_baseline(gamma_values, model_path)
            if baseline_df is not None:
                print("  Baseline (original kcat) computed.")
        plot_acetate_vs_gamma_panel(
            combined_all,
            baseline_df,
            title="Sensitivity analysis of kcat uncertainty by Monte Carlo sampling",
            output_dir=args.output_dir,
            filename_stem="Figure_S5a",
        )

    # Panel B: membrane only
    membrane_runs = discover_runs(args.membrane_only_dir)
    if not membrane_runs:
        print(f"No membrane-only data found under {args.membrane_only_dir}")
        print("  Provide run subdirs with analysis_metadata.json and aceT_data.csv.")
    else:
        print(f"Membrane only: {len(membrane_runs)} runs (CVs: {sorted(membrane_runs.keys())})")
        combined_membrane = load_and_combine_aceT(membrane_runs)
        gamma_values = sorted(combined_membrane["gamma"].unique())
        baseline_df = None
        if not args.no_baseline:
            baseline_df = compute_baseline(gamma_values, model_path)
        plot_acetate_vs_gamma_panel(
            combined_membrane,
            baseline_df,
            title="Sensitivity analysis of membrane kcat uncertainty by Monte Carlo sampling",
            output_dir=args.output_dir,
            filename_stem="Figure_S5b",
        )

    print("=" * 60)
    print(f"Outputs written to: {args.output_dir}")


if __name__ == "__main__":
    main()
