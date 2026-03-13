#!/usr/bin/env python3
"""
Runner for Supplementary Figures S6, S7, and S8.

Default behaviour (no flags): plot from existing data in data/outputs/Figure_S6_S7_S8/.
Pass --run-analysis to re-run the protein export gamma sweep first.

Usage:
    python run_figure_s6_s7_s8.py                  # plot only
    python run_figure_s6_s7_s8.py --run-analysis   # run simulation then plot
"""

import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))


def project_root() -> Path:
    for p in [SCRIPT_DIR, *SCRIPT_DIR.parents]:
        if (p / 'data').is_dir() and (p / 'Supplementary_Figures').is_dir():
            return p
    return SCRIPT_DIR.parent.parent.parent


PROJECT_ROOT = project_root()

# Input data
MODEL_PATH = PROJECT_ROOT / 'data' / 'input' / 'mafba_iML1515_ecV1_g1_MAFBA.xml'
ATTRIBUTES_PATH = PROJECT_ROOT / 'data' / 'input' / 'attributes_model_ecModel_iML1515_MAFBA_g_1_v3.xlsx'

# Output data dir (flat: one analysis result per run)
DATA_OUTPUT_DIR = PROJECT_ROOT / 'data' / 'outputs' / 'Figure_S6_S7_S8'

# Analysis parameters
EXPORT_FRACTIONS = [0.05, 0.10, 0.20, 0.30]
GAMMA_VALUES = [1.0, 0.35, 0.30, 0.25, 0.20, 0.15]
ACETATE_THRESHOLD = 0.01
SOLVER = 'gurobi'


def run_analysis() -> bool:
    """Run the protein export gamma sweep and write results to DATA_OUTPUT_DIR."""
    from protein_export_sweep import ProteinExportGammaSweep

    if not MODEL_PATH.exists():
        print(f"Model not found: {MODEL_PATH}")
        return False
    if not ATTRIBUTES_PATH.exists():
        print(f"Attributes file not found: {ATTRIBUTES_PATH}")
        return False

    DATA_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Running protein export gamma sweep...")
    analysis = ProteinExportGammaSweep(
        model_path=MODEL_PATH,
        attributes_path=ATTRIBUTES_PATH,
        output_dir=DATA_OUTPUT_DIR,
        solver=SOLVER,
    )
    analysis.load_model()
    analysis.run_analysis(
        export_fractions=EXPORT_FRACTIONS,
        gamma_values=GAMMA_VALUES,
        acetate_threshold=ACETATE_THRESHOLD,
    )
    print(f"  Analysis results written to: {analysis.output_dir}")
    return True


def run_plot(data_dir: Path) -> None:
    """Generate S6 and S7 from data in data_dir."""
    from figure_s6_s7_s8_plot import main as plot_main

    old_argv = sys.argv
    sys.argv = [old_argv[0],
                '--data-dir', str(data_dir),
                '--output-dir', str(SCRIPT_DIR / 'output')]
    try:
        plot_main()
    finally:
        sys.argv = old_argv


def find_latest_run(base_dir: Path) -> Path:
    """Return the most recently created timestamped subdirectory, or base_dir itself."""
    if (base_dir / 'analysis_metadata.json').exists():
        return base_dir
    subdirs = sorted([d for d in base_dir.iterdir() if d.is_dir()], reverse=True)
    for d in subdirs:
        if (d / 'analysis_metadata.json').exists():
            return d
    return base_dir


def main() -> None:
    parser = argparse.ArgumentParser(description='Runner for Supplementary Figures S6, S7, and S8.')
    parser.add_argument('--run-analysis', action='store_true',
                        help='Re-run the protein export gamma sweep before plotting.')
    args = parser.parse_args()

    print("Supplementary Figures S6 / S7 / S8")
    print("=" * 60)

    if args.run_analysis:
        ok = run_analysis()
        if not ok:
            print("Analysis failed; aborting.")
            sys.exit(1)

    data_dir = find_latest_run(DATA_OUTPUT_DIR)
    if not (data_dir / 'analysis_metadata.json').exists():
        print(f"No analysis data found in {DATA_OUTPUT_DIR}.")
        print("Run with --run-analysis to generate data first.")
        sys.exit(1)

    print(f"Plotting from: {data_dir}")
    run_plot(data_dir)
    print("=" * 60)
    print(f"Outputs written to: {SCRIPT_DIR / 'output'}")


if __name__ == '__main__':
    main()
