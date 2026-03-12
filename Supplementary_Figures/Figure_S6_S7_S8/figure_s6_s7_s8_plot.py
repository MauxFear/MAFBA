#!/usr/bin/env python3
"""
Generate Supplementary Figures S6, S7, and S8.

Reads output data from a protein export gamma sweep run and produces:
  Figure S6 – Acetate secretion profiles across γ values (one panel per γ,
               curves for each export fraction f, with vertical λac markers).
  Figure S7 – Proteome sector allocation at γ = 0.30 (one panel per f value,
               stacked area showing cumulative sector fractions vs growth rate).
  Figure S8 – Proteome sector allocation at γ = 0.25 (same layout as S7).

Expected input directory layout:
    <data_dir>/
        flux_data.csv
        sectors_data.csv   (or sectors_data_phiT.csv)
        aceL_data.csv
        analysis_metadata.json

Outputs are written to <output_dir>/Figure_S6.{png,pdf}, Figure_S7.{png,pdf},
and Figure_S8.{png,pdf}.
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_style("whitegrid")
plt.rcParams['font.size'] = 11
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['figure.dpi'] = 300
plt.rcParams['axes.linewidth'] = 1.2

# Sector colors (matching manuscript Figure 3b palette)
SECTOR_COLORS = {
    'phiCm': '#65AF65',
    'phiEm': '#D2C725',
    'phiEr': '#266DD7',
    'phiT':  '#E78AC3',
    'phiR':  '#9065D1',
    'phiEc': '#C74848',
    'phiCc': '#EC8327',
}

SECTOR_ORDER = ['phiCm', 'phiEm', 'phiEr', 'phiT', 'phiR', 'phiEc', 'phiCc']

# γ values for sector-allocation figures
S7_GAMMA = 0.30
S8_GAMMA = 0.25


def load_data(data_dir: Path):
    flux_data = pd.read_csv(data_dir / 'flux_data.csv')
    sectors_file = data_dir / 'sectors_data_phiT.csv'
    if not sectors_file.exists():
        sectors_file = data_dir / 'sectors_data.csv'
    sectors_data = pd.read_csv(sectors_file)
    aceL_data = pd.read_csv(data_dir / 'aceL_data.csv')
    with open(data_dir / 'analysis_metadata.json') as f:
        metadata = json.load(f)
    return flux_data, sectors_data, aceL_data, metadata


def plot_figure_s6(flux_data: pd.DataFrame, aceL_data: pd.DataFrame,
                   metadata: dict, output_dir: Path) -> None:
    """
    Figure S6: Acetate flux vs growth rate, one subplot per γ.
    Each subplot shows curves for each export fraction f, with a vertical
    dashed line at the corresponding acetate threshold λac.
    """
    export_fractions = metadata['parameters']['export_fractions']
    gamma_values = metadata['parameters']['gamma_values']

    n_gamma = len(gamma_values)
    n_cols = 3
    n_rows = int(np.ceil(n_gamma / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 5 * n_rows))
    axes = np.array(axes).flatten()

    colors = sns.color_palette("husl", len(export_fractions))

    for gamma_idx, gamma in enumerate(gamma_values):
        ax = axes[gamma_idx]

        for f_idx, export_frac in enumerate(export_fractions):
            subset = flux_data[
                (flux_data['gamma'] == gamma) &
                (flux_data['export_fraction'] == export_frac)
            ].copy().sort_values('mu')

            aceL_sub = aceL_data[
                (aceL_data['gamma'] == gamma) &
                (aceL_data['export_fraction'] == export_frac)
            ]

            if 'EX_ac_e_f' in subset.columns and len(subset) > 0:
                ax.plot(subset['mu'], subset['EX_ac_e_f'],
                        color=colors[f_idx], linewidth=2.5,
                        label=f'f = {export_frac:.2f}',
                        marker='o', markersize=4, alpha=0.8)

                if not aceL_sub.empty:
                    lam = aceL_sub['lambda_ac_mu'].values[0]
                    if not np.isnan(lam):
                        ax.axvline(lam, color=colors[f_idx],
                                   linestyle='--', linewidth=2, alpha=0.7, zorder=1)

        ax.set_xlabel('Growth rate \u03bc (1/h)', fontsize=11, fontweight='bold')
        ax.set_ylabel('Acetate flux (mmol/gDW/h)', fontsize=11, fontweight='bold')
        ax.set_title(f'\u03b3 = {gamma:.2f}', fontsize=12, fontweight='bold')
        ax.legend(loc='upper left', framealpha=0.9, fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)

    for idx in range(n_gamma, len(axes)):
        axes[idx].axis('off')

    plt.suptitle('Acetate Production vs Growth Rate (by \u03b3 value)',
                 fontsize=15, fontweight='bold', y=0.995)
    plt.tight_layout()

    for ext in ('png', 'pdf'):
        fig.savefig(output_dir / f'Figure_S6.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved Figure_S6.pdf/.png")


def plot_figure_s7(sectors_data: pd.DataFrame, aceL_data: pd.DataFrame,
                   metadata: dict, output_dir: Path) -> None:
    """Figure S7: proteome sector allocation at γ = 0.30, one panel per f."""
    _plot_sectors_at_gamma(sectors_data, aceL_data, metadata, S7_GAMMA, output_dir, 'Figure_S7')


def plot_figure_s8(sectors_data: pd.DataFrame, aceL_data: pd.DataFrame,
                   metadata: dict, output_dir: Path) -> None:
    """Figure S8: proteome sector allocation at γ = 0.25, one panel per f."""
    _plot_sectors_at_gamma(sectors_data, aceL_data, metadata, S8_GAMMA, output_dir, 'Figure_S8')


def _plot_sectors_at_gamma(sectors_data: pd.DataFrame, aceL_data: pd.DataFrame,
                            metadata: dict, gamma: float,
                            output_dir: Path, filename_stem: str) -> None:
    """
    Shared helper: proteome sector allocation at a fixed γ, one panel per f.
    Produces <filename_stem>.{png,pdf} in output_dir.
    """
    export_fractions = metadata['parameters']['export_fractions']
    phi_max = metadata['model_info']['phi_max']
    phiMmax = gamma * phi_max

    for col in ('phiT', 'phiCc'):
        if col not in sectors_data.columns:
            sectors_data = sectors_data.copy()
            sectors_data[col] = 0.0

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    legend_items = {}

    for idx, export_frac in enumerate(export_fractions):
        ax = axes[idx]

        subset = sectors_data[
            (sectors_data['export_fraction'] == export_frac) &
            (sectors_data['gamma'] == gamma)
        ].copy().sort_values('mu')

        aceL_sub = aceL_data[
            (aceL_data['export_fraction'] == export_frac) &
            (aceL_data['gamma'] == gamma)
        ]

        growth_rates = subset['mu'].values
        if len(growth_rates) == 0:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'f = {export_frac:.2f}', fontsize=12, fontweight='bold')
            continue

        slack_base = ['phiEc', 'phiEm', 'phiEr', 'phiT', 'phiCm', 'phiR']
        stacked = np.zeros(len(growth_rates))
        for s in slack_base:
            if s in subset.columns:
                stacked += subset[s].values
        subset = subset.copy()
        subset['phiCc'] = np.maximum(phi_max - stacked, 0.0)

        cumulative = np.zeros(len(growth_rates))
        for sector in SECTOR_ORDER:
            vals = subset[sector].values if sector in subset.columns else np.zeros(len(growth_rates))
            label = f'\u03c6{sector[3:]}' if sector != 'phiT' else '\u03c6T (translocation)'
            poly = ax.fill_between(
                growth_rates, cumulative, cumulative + vals,
                color=SECTOR_COLORS[sector], label=label,
                alpha=0.9, edgecolor='black', linewidth=0.5
            )
            if label not in legend_items:
                legend_items[label] = poly
            cumulative += vals

        phi_line = ax.axhline(phiMmax, color='red', linestyle='--', linewidth=2, alpha=0.8, zorder=10)
        if '\u03c6$M_{max}$' not in legend_items:
            phi_line.set_label('\u03c6$M_{max}$')
            legend_items['\u03c6$M_{max}$'] = phi_line

        phiP_line = ax.axhline(phi_max, color='black', linestyle='--', linewidth=1.5, alpha=0.9, zorder=9)
        if '\u03c6$P_{max}$' not in legend_items:
            phiP_line.set_label('\u03c6$P_{max}$')
            legend_items['\u03c6$P_{max}$'] = phiP_line

        if not aceL_sub.empty:
            lam = aceL_sub['lambda_ac_mu'].values[0]
            if not np.isnan(lam):
                lam_line = ax.axvline(lam, color='blue', linestyle=':', linewidth=2, alpha=0.7, zorder=10)
                if '$\u03bb_{ac}$' not in legend_items:
                    lam_line.set_label('$\u03bb_{ac}$')
                    legend_items['$\u03bb_{ac}$'] = lam_line

        ax.set_xlabel('Growth rate (1/h)', fontsize=11, fontweight='bold')
        ax.set_ylabel('Cumulative proteome fraction', fontsize=11, fontweight='bold')
        ax.set_title(f'f = {export_frac:.2f}', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(left=0, right=1.0)
        ax.set_ylim(bottom=0, top=phi_max * 1.1)

    plt.suptitle(f'Proteome Sector Allocation \u2013 \u03b3 = {gamma:.2f}',
                 fontsize=15, fontweight='bold', y=0.995)

    handles = list(legend_items.values())
    labels = list(legend_items.keys())
    fig.legend(handles, labels, loc='center left', bbox_to_anchor=(0.78, 0.5),
               framealpha=0.95, title='Proteome sectors', fontsize=10)

    plt.tight_layout(rect=[0, 0, 0.78, 0.96])
    for ext in ('png', 'pdf'):
        fig.savefig(output_dir / f'{filename_stem}.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {filename_stem}.pdf/.png")


def project_root() -> Path:
    here = Path(__file__).resolve().parent
    for p in [here, *here.parents]:
        if (p / 'data').is_dir() and (p / 'Supplementary_Figures').is_dir():
            return p
    return here.parent.parent.parent


def main() -> None:
    root = project_root()
    default_data = root / 'data' / 'outputs' / 'Figure_S6_S7_S8'
    default_out = Path(__file__).resolve().parent / 'output'

    parser = argparse.ArgumentParser(
        description='Generate Supplementary Figures S6, S7, and S8.'
    )
    parser.add_argument('--data-dir', type=Path, default=default_data,
                        dest='data_dir',
                        help='Directory containing analysis outputs (CSVs + metadata.json)')
    parser.add_argument('--output-dir', type=Path, default=default_out,
                        dest='output_dir',
                        help='Directory for generated figures')
    args = parser.parse_args()

    if not args.data_dir.exists():
        print(f"Data directory not found: {args.data_dir}")
        print("Run with --run-analysis flag via run_figure_s6_s7.py, or provide a valid --data-dir.")
        return

    required = ['flux_data.csv', 'aceL_data.csv', 'analysis_metadata.json']
    missing = [f for f in required if not (args.data_dir / f).exists()]
    if missing:
        print(f"Missing required files in {args.data_dir}: {missing}")
        return

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("Figure S6 / S7: loading data...")
    flux_data, sectors_data, aceL_data, metadata = load_data(args.data_dir)
    print(f"  flux_data: {len(flux_data)} rows")
    print(f"  sectors_data: {len(sectors_data)} rows")
    print(f"  aceL_data: {len(aceL_data)} rows")
    print(f"  γ values: {metadata['parameters']['gamma_values']}")
    print(f"  export fractions: {metadata['parameters']['export_fractions']}")

    print("Generating figures...")
    plot_figure_s6(flux_data, aceL_data, metadata, args.output_dir)
    plot_figure_s7(sectors_data, aceL_data, metadata, args.output_dir)
    plot_figure_s8(sectors_data, aceL_data, metadata, args.output_dir)

    print(f"Done. Outputs in: {args.output_dir}")


if __name__ == '__main__':
    main()
