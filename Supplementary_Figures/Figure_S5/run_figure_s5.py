#!/usr/bin/env python3
"""
Run Supplementary Figure S5: plot S5a and S5b from kcat sampling data.

Usage:
  python run_figure_s5.py                    # Plot only (uses existing data)
  python run_figure_s5.py --run-analysis     # Run kcat sampling for CV 10%, 20%, 50%, 100%, then plot

Data: timestamped run dirs with analysis_metadata.json and aceT_data.csv.
"""

import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

# CV set for figure (10%, 20%, 50%, 100%)
S5_CV_VALUES = [0.10, 0.20, 0.50, 1.0]
S5_GAMMA_VALUES = [1.0, 0.35, 0.30, 0.25, 0.20, 0.15]


def run_kcat_sampling():
    """Run kcat Monte Carlo sampling (all reactions and membrane-only) for each CV; writes to data/outputs/Figure_S5/."""
    out_all = PROJECT_ROOT / "data" / "outputs" / "Figure_S5" / "all_reactions"
    out_membrane = PROJECT_ROOT / "data" / "outputs" / "Figure_S5" / "membrane_only"
    model_path = PROJECT_ROOT / "data" / "input" / "mafba_iML1515_ecV1_g1_MAFBA.xml"
    kcat_path = PROJECT_ROOT / "data" / "input" / "attributes_model_ecModel_iML1515_MAFBA_g_1_v3.xlsx"
    if not model_path.exists() or not kcat_path.exists():
        print("Missing model or kcat file under data/input/.")
        return False
    out_all.mkdir(parents=True, exist_ok=True)
    out_membrane.mkdir(parents=True, exist_ok=True)
    solver = "gurobi"
    n_iter = 500
    kcat_threshold = 1e7

    try:
        from kcat_sampling_all import KcatMonteCarloAnalysis
    except ImportError as e:
        print(f"Cannot import kcat_sampling_all: {e}")
        return False
    print("Running kcat sampling (all reactions) for CV =", S5_CV_VALUES)
    for cv in S5_CV_VALUES:
        a = KcatMonteCarloAnalysis(
            model_path=model_path,
            kcat_data_path=kcat_path,
            output_dir=out_all,
            solver=solver,
        )
        a.run_full_analysis(
            n_iterations=n_iter,
            kcat_threshold=kcat_threshold,
            cv=cv,
            gamma_values=S5_GAMMA_VALUES,
        )
        print(f"  All reactions CV={cv:.0%} -> {a.output_dir}")

    try:
        from kcat_sampling_membrane import KcatMonteCarloMembraneAnalysis
    except ImportError as e:
        print(f"Cannot import kcat_sampling_membrane: {e}")
        return False
    print("Running kcat sampling (membrane only) for CV =", S5_CV_VALUES)
    for cv in S5_CV_VALUES:
        a = KcatMonteCarloMembraneAnalysis(
            model_path=model_path,
            kcat_data_path=kcat_path,
            output_dir=out_membrane,
            solver=solver,
        )
        a.run_full_analysis(
            n_iterations=n_iter,
            kcat_threshold=kcat_threshold,
            cv=cv,
            gamma_values=S5_GAMMA_VALUES,
        )
        print(f"  Membrane only CV={cv:.0%} -> {a.output_dir}")
    return True


def main():
    parser = argparse.ArgumentParser(description="Generate Figure S5 (S5a + S5b)")
    parser.add_argument(
        "--run-analysis",
        action="store_true",
        help="Run kcat sampling for CV 10%%, 20%%, 50%%, 100%% then plot",
    )
    args = parser.parse_args()

    if args.run_analysis:
        if not run_kcat_sampling():
            sys.exit(1)
    # Call the plot with an empty argument list so figure_s5_plot's argparse
    # uses its defaults (data/outputs/Figure_S5/... paths).
    old_argv = sys.argv
    sys.argv = [old_argv[0]]
    try:
        from figure_s5_plot import main as plot_main
        plot_main()
    finally:
        sys.argv = old_argv


if __name__ == "__main__":
    main()
