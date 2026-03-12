"""
Supplementary Figure S1 - Homogeneous Proteome Sector Allocation (Python)
=========================================================================

Alternative Python implementation for the homogeneous sector-allocation
workflow preserved as Supplementary Figure S1.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os

def plot_sectors(sectors_filename, data_filename, figure_prefix):
    # Read the data
    try:
        table_sectors = pd.read_csv(sectors_filename)
        table_data = pd.read_csv(data_filename)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return

    # Get unique gamma values
    gamma_values = [0.2, 0.25, 0.3, 1]
    print(f"Plotting for gamma values: {gamma_values}")

    # Define colors
    colors = {
        'Cm': '#65AF65', 'Em': '#D2C725', 'Er': '#266DD7', 
        'R': '#9065D1', 'Ec': '#C74848', 'Cc': '#EC8327',
        'ac_line': 'lightseagreen'
    }

    # Loop through each gamma value
    for gamma_val in gamma_values:
        # Detect acetate threshold
        gr_ac = 0.75 # Hardcoding for now

        # Filter sectors data
        S = table_sectors[table_sectors['Gamma'] == gamma_val].copy()
        S = S.sort_values(by='Growth').reset_index(drop=True)

        # Calculate phiPmax and phiMmax
        phiPmax = (S['phiR'] + S['phiEc'] + S['phiEr'] + S['phiCm'] + S['phiEm']).max()
        phiMmax = gamma_val * phiPmax

        # Prepare data for stacked area plot
        gr = S['Growth']
        phiCc = phiPmax - (S['phiR'] + S['phiEc'] + S['phiEr'] + S['phiCm'] + S['phiEm'])
        
        plot_data = {
            r'$\phi$ C$_{m}$': S['phiCm'],
            r'$\phi$ E$_{m}$': S['phiEm'],
            r'$\phi$ E$_{r}$': S['phiEr'],
            r'$\phi$ R': S['phiR'],
            r'$\phi$ E$_{c}$': S['phiEc'],
            r'$\phi$ C$_{c}$': phiCc
        }
        
        labels = plot_data.keys()
        data_values = list(plot_data.values())
        plot_colors = [colors['Cm'], colors['Em'], colors['Er'], colors['R'], colors['Ec'], colors['Cc']]

        # Create the plot
        fig, ax = plt.subplots(figsize=(9, 8))
        
        # Set hatch line width
        plt.rcParams['hatch.linewidth'] = 2.0

        ax.stackplot(gr, data_values, labels=labels, colors=plot_colors, edgecolor='black', linewidth=1)

        hatch_map = {
            r'$\phi$ C$_{m}$': '||',
            r'$\phi$ E$_{m}$': '\\\\',
            r'$\phi$ E$_{r}$': '//',
            r'$\phi$ R': '',
            r'$\phi$ E$_{c}$': '',
            r'$\phi$ C$_{c}$': ''
        }
        
        from matplotlib.patches import Patch
        handles = []
        y_bottom = np.zeros(len(gr))
        for i, label in enumerate(labels):
            y_top = y_bottom + data_values[i]
            hatch = hatch_map.get(label)
            
            if hatch:
                ax.fill_between(gr, y_bottom, y_top, hatch=hatch, facecolor='none', edgecolor='black', linewidth=0)

            if hatch:
                handles.append(Patch(facecolor=plot_colors[i], hatch=hatch, label=label, edgecolor='black'))
            else:
                handles.append(Patch(facecolor=plot_colors[i], label=label, edgecolor='black'))
            
            y_bottom = y_top
            
        cumulative_data = np.cumsum(data_values, axis=0)
        for data in cumulative_data:
            ax.plot(gr, data, color='black', linewidth=1)

        ac_line = Line2D([0], [0], color=colors['ac_line'], linestyle='--', linewidth=3, label=r'$\lambda_{ac}$')
        handles.append(ac_line)
        
        ax.set_title(f'$\\gamma = {gamma_val}$', fontsize=22)
        ax.set_xlabel('Growth rate (h$^{-1}$)', fontsize=20)
        ax.set_ylabel('Cumulative proteome fraction', fontsize=20)
        
        ax.axhline(y=phiPmax, color='black', linestyle='--', linewidth=1)
        ax.text(gr.max() * 1.01, phiPmax, f'$\\phi^P_{{max}} = {phiPmax:.2f}$', 
                va='center', ha='left', fontsize=20)

        ax.axhline(y=phiMmax, color='black', linestyle='--', linewidth=2)
        ax.text(gr.max()* 1.01, phiMmax, f'$\\phi^M_{{max}} = {phiMmax:.2f}$', 
                va='center', ha='left', fontsize=20)
        
        ax.axvline(x=gr_ac, color='black', linestyle='--', linewidth=3.5)
        ax.axvline(x=gr_ac, color=colors['ac_line'], linestyle='--', linewidth=2.5)

        ax.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.45), ncol=4, frameon=False, fontsize=20)
        
        ax.set_xlim(0, float(gr.max()))
        ax.set_ylim(0, phiPmax)
        ax.tick_params(axis='both', which='major', labelsize=18)

        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.35)
        
        output_folder = 'output'
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        
        filename_png = os.path.join(output_folder, f'{figure_prefix}_sectors_g{int(gamma_val*100):03d}.png')
        plt.savefig(filename_png, dpi=300, bbox_inches='tight')
        
        filename_svg = os.path.join(output_folder, f'{figure_prefix}_sectors_g{int(gamma_val*100):03d}.svg')
        plt.savefig(filename_svg, dpi=300, bbox_inches='tight')
        
        print(f"  Saved: {filename_png}")
        print(f"  Saved: {filename_svg}")
        plt.show()

if __name__ == '__main__':
    print("="*80)
    print("SUPPLEMENTARY FIGURE S1 - HOMOGENEOUS PROTEOME SECTOR ALLOCATION (Python)")
    print("="*80)
    
    # Get paths relative to script location
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(os.path.dirname(script_dir))
    data_folder = os.path.join(repo_root, 'data', 'outputs', 'Sectors')
    
    print(f"\nData folder: {data_folder}")
    
    # Data files
    sectors_file = os.path.join(data_folder, 'sectors_modelV0_kcatscomb_Homog.csv')
    data_file = os.path.join(data_folder, 'data_modelV0_kcatscomb_Homog.csv')
    
    # Check files exist
    if not os.path.exists(sectors_file):
        print(f"ERROR: Sectors file not found: {sectors_file}")
        exit(1)
    if not os.path.exists(data_file):
        print(f"ERROR: Data file not found: {data_file}")
        exit(1)
    
    print(f"Sectors file: {os.path.basename(sectors_file)}")
    print(f"Data file: {os.path.basename(data_file)}")
    print("\nGenerating figures...")
    
    # Generate plots
    plot_sectors(sectors_file, data_file, 'Homogeneous')
    
    print("\n✓ Supplementary Figure S1 generation complete!")
