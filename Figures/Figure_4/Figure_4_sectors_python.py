"""
Figure 4 - Parameterized Proteome Sector Allocation (Python version)
===================================================================

Alternative Python implementation for generating Figure 4 sector plots
from the same baseline sensitivity dataset used by Figure 5.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


def select_float_rows(table_data, column, target, atol=1e-9):
    return table_data[np.isclose(table_data[column], target, atol=atol)].copy()


def resolve_column(table_data, *candidates):
    for candidate in candidates:
        if candidate in table_data.columns:
            return candidate
    raise KeyError(f"None of the candidate columns were found: {candidates}")


def detect_acetate_threshold(table_data):
    if table_data.empty:
        return 0.0

    growth_col = resolve_column(table_data, "rx_BIOMASS_Ec_iML1515_WT_75p37M", "BIOMASS_Ec_iML1515_WT_75p37M")
    memb_col = resolve_column(table_data, "rx_ATPS4rpp_f", "ATPS4rpp_f")
    acetate_col = resolve_column(table_data, "rx_EX_ac_e_f", "EX_ac_e_f")

    table_data = table_data.sort_values(by=growth_col).reset_index(drop=True)
    gr = table_data[growth_col].to_numpy()
    memb = table_data[memb_col].to_numpy()
    ac = table_data[acetate_col].to_numpy()

    if len(gr) < 2:
        return float(gr[0]) if len(gr) == 1 else 0.0

    with np.errstate(divide="ignore", invalid="ignore"):
        slope = np.diff(memb) / np.diff(gr)

    valid = slope[np.isfinite(slope)]
    if valid.size == 0:
        ac_idx = np.where(ac >= 0.01)[0]
        return float(gr[ac_idx[0]]) if ac_idx.size else float(gr[0])

    max_slope = valid.max()
    for threshold, acetate_cutoff in [(0.01, 0.01), (0.05, 0.0), (0.2, 0.0)]:
        memb_sat = np.where(slope <= max_slope * threshold)[0]
        if memb_sat.size == 0:
            continue
        acetate_idx = np.where(ac[memb_sat] >= acetate_cutoff)[0]
        if acetate_idx.size:
            return float(gr[memb_sat[acetate_idx[0]]])

    return float(gr[0])


def plot_sectors(sectors_filename, data_filename, figure_prefix, output_folder):
    try:
        table_sectors = pd.read_csv(sectors_filename)
        table_data = pd.read_csv(data_filename)
    except FileNotFoundError as exc:
        print(f"Error: {exc}")
        return

    gamma_values = [1.0, 0.3, 0.25]
    change_factor = 1.0
    print(f"Plotting for gamma values: {gamma_values}")
    if "Change_Factor" in table_data.columns:
        print(f"Using Change_Factor = {change_factor}")
    else:
        print("No Change_Factor column found. Using the full dataset.")

    colors = {
        "Cm": "#65AF65",
        "Em": "#D2C725",
        "Er": "#266DD7",
        "R": "#9065D1",
        "Ec": "#C74848",
        "Cc": "#EC8327",
        "ac_line": "lightseagreen",
    }

    os.makedirs(output_folder, exist_ok=True)

    for gamma_val in gamma_values:
        flux_slice = select_float_rows(table_data, "Gamma", gamma_val)
        if "Change_Factor" in flux_slice.columns:
            flux_slice = select_float_rows(flux_slice, "Change_Factor", change_factor)
        if flux_slice.empty:
            print(f"Skipping gamma {gamma_val}: no baseline flux data")
            continue

        sectors_slice = select_float_rows(table_sectors, "Gamma", gamma_val)
        if "Change_Factor" in sectors_slice.columns:
            sectors_slice = select_float_rows(sectors_slice, "Change_Factor", change_factor)
        if sectors_slice.empty:
            print(f"Skipping gamma {gamma_val}: no baseline sector data")
            continue

        sectors_slice = sectors_slice.sort_values(by="Growth").reset_index(drop=True)
        gr_ac = detect_acetate_threshold(flux_slice)

        phi_p_max = (
            sectors_slice["phiR"]
            + sectors_slice["phiEc"]
            + sectors_slice["phiEr"]
            + sectors_slice["phiCm"]
            + sectors_slice["phiEm"]
        ).max()
        phi_m_max = gamma_val * phi_p_max
        gr = sectors_slice["Growth"]
        phi_cc = phi_p_max - (
            sectors_slice["phiR"]
            + sectors_slice["phiEc"]
            + sectors_slice["phiEr"]
            + sectors_slice["phiCm"]
            + sectors_slice["phiEm"]
        )
        plot_data = {
            r"$\phi$ C$_{m}$": sectors_slice["phiCm"],
            r"$\phi$ E$_{m}$": sectors_slice["phiEm"],
            r"$\phi$ E$_{r}$": sectors_slice["phiEr"],
            r"$\phi$ R": sectors_slice["phiR"],
            r"$\phi$ E$_{c}$": sectors_slice["phiEc"],
            r"$\phi$ C$_{c}$": phi_cc,
        }

        labels = list(plot_data.keys())
        data_values = list(plot_data.values())
        plot_colors = [colors["Cm"], colors["Em"], colors["Er"], colors["R"], colors["Ec"], colors["Cc"]]

        fig, ax = plt.subplots(figsize=(9, 8))
        plt.rcParams["hatch.linewidth"] = 2.0
        ax.stackplot(gr, data_values, labels=labels, colors=plot_colors, edgecolor="black", linewidth=1)

        hatch_map = {
            r"$\phi$ C$_{m}$": "||",
            r"$\phi$ E$_{m}$": "\\\\",
            r"$\phi$ E$_{r}$": "//",
            r"$\phi$ R": "",
            r"$\phi$ E$_{c}$": "",
            r"$\phi$ C$_{c}$": "",
        }

        handles = []
        y_bottom = np.zeros(len(gr))
        for i, label in enumerate(labels):
            y_top = y_bottom + data_values[i]
            hatch = hatch_map.get(label)
            if hatch:
                ax.fill_between(gr, y_bottom, y_top, hatch=hatch, facecolor="none", edgecolor="white", linewidth=0)
                handles.append(Patch(facecolor=plot_colors[i], hatch=hatch, label=label, edgecolor="white"))
            else:
                handles.append(Patch(facecolor=plot_colors[i], label=label, edgecolor="white"))
            y_bottom = y_top

        cumulative_data = np.cumsum(data_values, axis=0)
        for data in cumulative_data:
            ax.plot(gr, data, color="black", linewidth=1)

        ac_line = Line2D([0], [0], color=colors["ac_line"], linestyle="--", linewidth=3, label=r"$\lambda_{ac}$")
        handles.append(ac_line)

        ax.set_title(f"$\\gamma = {gamma_val}$", fontsize=22)
        ax.set_xlabel("Growth rate (h$^{-1}$)", fontsize=20)
        ax.set_ylabel("Cumulative proteome fraction", fontsize=20)

        ax.axhline(y=phi_p_max, color="black", linestyle="--", linewidth=1)
        ax.text(gr.max() * 1.01, phi_p_max, f"$\\phi^P_{{max}} = {phi_p_max:.2f}$", va="center", ha="left", fontsize=20)

        ax.axhline(y=phi_m_max, color="black", linestyle="--", linewidth=2)
        ax.text(gr.max() * 1.01, phi_m_max, f"$\\phi^M_{{max}} = {phi_m_max:.2f}$", va="center", ha="left", fontsize=20)

        ax.axvline(x=gr_ac - 0.0005, color="black", linestyle="-", linewidth=3.5)
        ax.axvline(x=gr_ac, color=colors["ac_line"], linestyle="--", linewidth=3)

        ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.45), ncol=4, frameon=False, fontsize=20)
        ax.set_xlim(0, float(gr.max()))
        ax.set_ylim(0, phi_p_max)
        ax.tick_params(axis="both", which="major", labelsize=18)
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.35)

        base_name = f"{figure_prefix}_sectors_g{int(gamma_val * 100):03d}_cf100"
        filename_png = os.path.join(output_folder, f"{base_name}.png")
        plt.savefig(filename_png, dpi=300, bbox_inches="tight")
        filename_svg = os.path.join(output_folder, f"{base_name}.svg")
        plt.savefig(filename_svg, dpi=300, bbox_inches="tight")
        print(f"  Saved: {filename_png}")
        print(f"  Saved: {filename_svg}")
        plt.show()


if __name__ == "__main__":
    print("=" * 80)
    print("FIGURE 4 - PARAMETERIZED PROTEOME SECTOR ALLOCATION (Python)")
    print("=" * 80)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(os.path.dirname(script_dir))
    data_folder = os.path.join(repo_root, "data", "outputs", "Parametrized")

    print(f"\nData folder: {data_folder}")

    sectors_file = os.path.join(data_folder, "sectors_modelV0_kcatscomb_Parametrized.csv")
    data_file = os.path.join(data_folder, "data_modelV0_kcatscomb_Parametrized.csv")

    if not os.path.exists(sectors_file):
        print(f"ERROR: Sectors file not found: {sectors_file}")
        raise SystemExit(1)
    if not os.path.exists(data_file):
        print(f"ERROR: Data file not found: {data_file}")
        raise SystemExit(1)

    print(f"Sectors file: {os.path.basename(sectors_file)}")
    print(f"Data file: {os.path.basename(data_file)}")
    print("\nGenerating figures...")

    output_folder = os.path.join(script_dir, "output")
    print(f"Output folder: {output_folder}")

    plot_sectors(sectors_file, data_file, "iML1515_MAFBA_parametrized", output_folder)

    print("\n✓ Figure 4 generation complete!")
