# cspell: disable
"""
Figure 1 - Protein Concentration Analysis
==========================================

This script generates Figure 1B and Figure 1C from the manuscript:
- Figure 1B: Total protein concentration (g/gDCW) vs growth rate
- Figure 1C: Inner membrane protein concentration (g/gDCW) vs growth rate

Data sources:
- Schmidt et al., 2016 (Proteomics data - Tables S11, S13, S23)
  
Output files:
- Figure_1B_schmidt2016_ProtConc.png/svg
- Figure_1C_schmidt2016_InnerMembProt_Loc_Conc.png/svg
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import os

# ==========================================
# CONFIGURATION & CONSTANTS
# ==========================================

# Output directory
OUTPUT_DIR = 'output'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Input data paths
SCHMIDT_DATA_PATH = '../../data/input/Proteomics/Schmidt2016_DataToPlot.xlsx'

# Physical constants
CELL_DRY_WEIGHT = 268.36  # g/L (molar cell weight)

# Plot configuration
PLOT_DPI_PNG = 150
PLOT_DPI_SVG = 300
FONT_SIZE = 24
FONT_FAMILY = 'serif'
FIGSIZE = (10, 8)

# Color palette
COLORS = {
    'Navy Blue': '#216DB7',
    'Steel Blue': '#4682B4',
    'Dark Red': '#8B0000',
    'Crimson': '#DC143C',
    'Forest Green': '#228B22',
    'Olive': '#6B8E23'
}

# ==========================================
# HELPER FUNCTIONS
# ==========================================

def load_excel_all_sheets(file_path):
    """Load all sheets from an Excel file."""
    try:
        sheets = pd.read_excel(file_path, sheet_name=None)
        print(f"Successfully loaded: {file_path}")
        return sheets
    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return {}
    except Exception as e:
        print(f"Error loading file: {e}")
        return {}


def prepare_schmidt_table(df, header_row_idx, start_row_idx):
    """
    Prepare a Schmidt data table by setting proper column names and extracting data rows.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Raw dataframe from Excel sheet
    header_row_idx : int
        Row index containing column headers
    start_row_idx : int
        Row index where data starts
        
    Returns:
    --------
    pd.DataFrame
        Cleaned dataframe with proper headers
    """
    df.columns = df.iloc[header_row_idx, :]
    return df.iloc[start_row_idx:, :].copy(deep=True)


def create_condition_code(df, condition_col, strain_col):
    """Create a unique condition code by combining condition and strain."""
    return df[condition_col] + '_' + df[strain_col]


def plot_concentration_vs_growth_rate(df, x_col, y_col, y_err_col, 
                                     ylabel, filename, color, ylim=None):
    """
    Create a publication-quality concentration vs growth rate plot.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Data to plot
    x_col : str
        Column name for x-axis (growth rate)
    y_col : str
        Column name for y-axis (concentration)
    y_err_col : str
        Column name for y error bars
    ylabel : str
        Y-axis label
    filename : str
        Output filename (without extension)
    color : str
        Color for data points
    ylim : tuple, optional
        Y-axis limits (min, max)
    """
    fig, ax = plt.subplots(figsize=FIGSIZE)
    plt.rcParams.update({'font.size': 18})
    
    # Extract data
    x_data = df[x_col]
    y_data = df[y_col]
    y_err = df[y_err_col]
    
    # Calculate statistics
    mean_val = np.mean(y_data)
    std_err = stats.sem(y_data)
    conf_interval = std_err * stats.t.ppf((1 + 0.95) / 2., len(y_data) - 1)
    
    # Plot data with error bars
    ax.errorbar(x_data, y_data, yerr=y_err, fmt='o', 
                color=color, markersize=12, capsize=5, alpha=0.8,
                label='Protein Concentration')
    
    # Add mean line
    ax.axhline(y=mean_val, color=COLORS['Crimson'], 
               linestyle='--', linewidth=2)
    
    # Add confidence interval shading
    x_min, x_max = ax.get_xlim()
    fill_x = np.insert(x_data.astype(float).values, 0, x_min)
    fill_x = np.append(fill_x, x_max)
    ax.fill_between(fill_x, 
                     [mean_val - conf_interval] * len(fill_x),
                     [mean_val + conf_interval] * len(fill_x),
                     color=color, alpha=0.4, label='95% Confidence Interval')
    
    # Add mean annotation
    ax.text(0.95, mean_val - 0.003 * (ax.get_ylim()[1] - ax.get_ylim()[0]), 
            f'Mean: {mean_val:.3f}', 
            ha='right', va='center', fontsize=FONT_SIZE, 
            fontfamily=FONT_FAMILY, transform=ax.get_yaxis_transform())
    
    # Customize plot
    ax.set_xlabel(r'Growth rate ($h^{-1}$)', fontsize=FONT_SIZE, fontfamily=FONT_FAMILY)
    ax.set_ylabel(ylabel, fontsize=FONT_SIZE, fontfamily=FONT_FAMILY)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=FONT_SIZE - 2)
    ax.grid(False)
    ax.set_facecolor('white')
    ax.set_xlim(x_min, x_max)
    
    if ylim:
        ax.set_ylim(ylim)
    
    # Save plot
    output_path = os.path.join(OUTPUT_DIR, filename)
    fig.savefig(f'{output_path}.png', format='png', bbox_inches='tight', 
                transparent=False, dpi=PLOT_DPI_PNG)
    fig.savefig(f'{output_path}.svg', format='svg', bbox_inches='tight', 
                transparent=False, dpi=PLOT_DPI_SVG)
    
    plt.show()
    print(f"Saved: {filename}")


# ==========================================
# DATA PROCESSING FUNCTIONS
# ==========================================

def load_and_prepare_metadata(all_sheets):
    """
    Load and prepare sample metadata and condition information.
    
    Returns:
    --------
    tuple : (samples_df, conditions_df)
        samples_df: DataFrame with sample metadata indexed by File Name
        conditions_df: DataFrame with condition metadata indexed by ConditionCode
    """
    # Process SamplesIds sheet
    samples = all_sheets['SamplesIds']
    samples.columns = samples.iloc[1, :]
    samples = samples.iloc[2:, :].copy(deep=True)
    samples.set_index('File Name', inplace=True)
    samples['ConditionCode'] = create_condition_code(
        samples, 
        'Growth Condition (in biological triplicates)', 
        'Strain'
    )
    
    # Process Table S23 (conditions)
    conditions = all_sheets['Table S23']
    conditions.columns = conditions.iloc[1, :]
    conditions = conditions.iloc[2:, :].copy(deep=True)
    conditions['ConditionCode'] = create_condition_code(
        conditions,
        'Growth condition',
        'Strain'
    )
    conditions.set_index('ConditionCode', inplace=True)
    
    return samples, conditions


def process_table_s11_total_protein(all_sheets, samples, conditions):
    """
    Process Table S11 to calculate total protein concentration.
    
    Returns:
    --------
    pd.DataFrame
        Summary table with growth rate and total protein concentration
    """
    table_s11 = all_sheets['Table S11']
    table_s11 = prepare_schmidt_table(table_s11, header_row_idx=2, start_row_idx=3)
    
    # Sum protein masses for each sample - USE ONLY FIRST 22 ROWS (protein groups)
    subtable_s11 = table_s11.iloc[:22].copy(deep=True)
    
    # Convert all columns to numeric and sum
    subtable_s11_numeric = subtable_s11.apply(pd.to_numeric, errors='coerce')
    column_sums = subtable_s11_numeric.sum(axis=0, numeric_only=True)
    
    # Create DataFrame with sums
    sum_table = column_sums.reset_index()
    sum_table.columns = ['Column', 'Sum']
    
    # Add experimental growth rate from row 23 (index 20 after header removal)
    exp_growth_rate_row = table_s11.iloc[20, :]
    sum_table['Experimental_Growth_rate'] = exp_growth_rate_row.values
    
    # Drop first row (typically a metadata column like 'COG Group')
    sum_table.drop(0, inplace=True)
    sum_table.set_index('Column', inplace=True)
    
    # Merge with sample metadata to get ConditionCode
    sum_table = pd.merge(sum_table, samples, left_index=True, right_index=True, how='left')
    
    # Add volume data from Table S23
    volume_dict = conditions['Single cell volume [fl]1'].to_dict()
    sum_table['Volume'] = [volume_dict.get(key, np.nan) for key in sum_table['ConditionCode']]
    sum_table['Concentration'] = sum_table['Sum'] / (sum_table['Volume'] * CELL_DRY_WEIGHT)
    
    # Group by condition
    grouped = sum_table.groupby('ConditionCode').agg({
        'Sum': ['mean', 'std'],
        'Concentration': ['mean', 'std'],
        'Experimental_Growth_rate': ['mean']
    })
    grouped.columns = ['Sum', 'Sum_std', 'Concentration', 'Conc_std', 'growthrate']
    
    # Merge with condition metadata to get volume and protein count
    grouped = pd.merge(grouped, conditions[['Single cell volume [fl]1', 'Number of Proteins Identified (FDR 1%)2']], 
                      left_index=True, right_index=True, how='left')
    
    # Create summary table
    summary = pd.DataFrame()
    summary['Total Protein Mass (fg/cell)'] = grouped['Sum']
    summary['Total Protein Mass Stdev'] = grouped['Sum_std']
    summary['Growth rate (h-1)'] = grouped['growthrate']
    summary['Single cell volume (fL)'] = grouped['Single cell volume [fl]1']
    summary['Number of Proteins Identified'] = grouped['Number of Proteins Identified (FDR 1%)2']
    summary['Total Protein Concentration (gProtein / L)'] = (
        summary['Total Protein Mass (fg/cell)'] / summary['Single cell volume (fL)']
    )
    summary['Total Protein Concentration (gProtein / gDCW)'] = (
        summary['Total Protein Concentration (gProtein / L)'] / CELL_DRY_WEIGHT
    )
    summary['Total Protein Mass per DCW Stdev'] = (
        summary['Total Protein Mass Stdev'] / summary['Single cell volume (fL)']
    ) / CELL_DRY_WEIGHT
    
    return summary.sort_values(by='Growth rate (h-1)')


def process_table_s13_membrane_protein(all_sheets, samples, conditions):
    """
    Process Table S13 to calculate inner membrane protein concentration.
    
    Returns:
    --------
    pd.DataFrame
        Summary table with growth rate and membrane protein concentration by location
    """
    table_s13 = all_sheets['Table S13']
    table_s13 = prepare_schmidt_table(table_s13, header_row_idx=2, start_row_idx=3)
    table_s13.set_index('Uniprot Accession', inplace=True)
    
    # Standardize location names
    location_map = {
        'Cell Inner Membrane': 'Cell Inner Membrane',
        'Cell inner membrane': 'Cell Inner Membrane',
        'Cell membrane': 'Cell Inner Membrane',
        'Cell outer membrane': 'Cell Outer Membrane',
        'Cytoplasm': 'Cytoplasm',
        'Membrane': 'Cell Inner Membrane',
        'Outer Membrane': 'Cell Outer Membrane',
        'Periplasm': 'Periplasm',
        'Secreted': 'Secreted'
    }
    
    table_s13['Cellular Protein Location'] = [
        location_map.get(x, x) 
        for x in table_s13['Cellular protein location (according to www.uniprot.org)']
    ]
    
    # Get sample columns (excluding metadata columns)
    sample_columns = [col for col in table_s13.columns[3:-1] if col in samples.index]
    
    # Convert sample columns to numeric
    for col in sample_columns:
        table_s13[col] = pd.to_numeric(table_s13[col], errors='coerce')
    
    # Group by location and sum protein masses for each sample
    grouped_by_location = table_s13.groupby('Cellular Protein Location')[sample_columns].agg('sum')
    
    # Transpose to get samples as rows, locations as columns
    transposed = grouped_by_location.T
    transposed.index.name = 'Sample Id'
    
    # Merge with sample metadata
    location_data = pd.merge(transposed, samples, left_index=True, right_index=True, how='left')
    
    # Get numeric columns (location columns)
    location_columns = transposed.columns.tolist()
    
    # Ensure numeric types
    for col in location_columns:
        location_data[col] = pd.to_numeric(location_data[col], errors='coerce')
    
    # Group by condition and calculate mean/std for each location
    agg_dict = {col: ['mean', 'std'] for col in location_columns}
    location_summary = location_data.dropna(subset=['ConditionCode']).groupby('ConditionCode').agg(agg_dict)
    location_summary.columns = ['_'.join(col).strip() for col in location_summary.columns]
    
    # Store location column names before merge
    mean_columns = [col for col in location_summary.columns if '_mean' in str(col)]
    std_columns = [col for col in location_summary.columns if '_std' in str(col)]
    
    # Merge with condition metadata
    location_summary = pd.merge(location_summary, conditions, left_index=True, right_index=True, how='left')
    
    # Calculate protein concentrations (divide by volume * cell weight)
    
    for col in mean_columns + std_columns:
        location_summary[col] = location_summary[col].div(
            location_summary['Single cell volume [fl]1'] * CELL_DRY_WEIGHT, 
            axis=0
        )
    
    # Calculate total protein concentration
    location_summary['Total Protein Concentration'] = location_summary[mean_columns].sum(axis=1)
    
    return location_summary


# ==========================================
# MAIN ANALYSIS
# ==========================================

def main():
    """Execute the complete Figure 1 analysis pipeline."""
    
    print("="*80)
    print("FIGURE 1 GENERATION - PROTEIN CONCENTRATION ANALYSIS")
    print("="*80)
    
    # -------------------------------------------------------------------------
    # 1. Load Data
    # -------------------------------------------------------------------------
    print("\n1. Loading Schmidt et al. 2016 data...")
    all_sheets = load_excel_all_sheets(SCHMIDT_DATA_PATH)
    
    if not all_sheets:
        print("Error: Failed to load data. Exiting.")
        return
    
    # -------------------------------------------------------------------------
    # 2. Prepare Metadata
    # -------------------------------------------------------------------------
    print("\n2. Processing metadata (sample IDs and conditions)...")
    samples, conditions = load_and_prepare_metadata(all_sheets)
    print(f"   - {len(samples)} samples loaded")
    print(f"   - {len(conditions)} conditions loaded")
    
    # -------------------------------------------------------------------------
    # 3. Figure 1B: Total Protein Concentration
    # -------------------------------------------------------------------------
    print("\n3. Processing Figure 1B: Total Protein Concentration...")
    summary_total = process_table_s11_total_protein(all_sheets, samples, conditions)
    
    print("\n   Creating Figure 1B plot...")
    plot_concentration_vs_growth_rate(
        df=summary_total,
        x_col='Growth rate (h-1)',
        y_col='Total Protein Concentration (gProtein / gDCW)',
        y_err_col='Total Protein Mass per DCW Stdev',
        ylabel=r'Protein Concentration ($g/g_{DCW}$)',
        filename='Figure_1B_schmidt2016_ProtConc',
        color=COLORS['Forest Green'],
        ylim=(0.312, 0.326)
    )
    
    # -------------------------------------------------------------------------
    # 4. Figure 1C: Inner Membrane Protein Concentration
    # -------------------------------------------------------------------------
    print("\n4. Processing Figure 1C: Inner Membrane Protein Concentration...")
    location_summary = process_table_s13_membrane_protein(all_sheets, samples, conditions)
    
    print("\n   Creating Figure 1C plot...")
    
    # Prepare data for plotting
    plot_data = location_summary[['Growth rate (h-1)', 
                                   'Cell Inner Membrane_mean', 
                                   'Cell Inner Membrane_std']].copy()
    
    # Create figure with dual y-axis
    fig, ax = plt.subplots(figsize=FIGSIZE)
    plt.rcParams.update({'font.size': 18})
    
    x_data = plot_data['Growth rate (h-1)']
    y_data = plot_data['Cell Inner Membrane_mean']
    y_err = plot_data['Cell Inner Membrane_std']
    
    # Calculate statistics
    mean_val = np.mean(y_data)
    std_err = stats.sem(y_data)
    conf_interval = std_err * stats.t.ppf((1 + 0.90) / 2., len(y_data) - 1)
    
    # Plot data
    ax.errorbar(x_data, y_data, yerr=y_err, fmt='o',
                color=COLORS['Navy Blue'], markersize=10, capsize=5, alpha=0.7,
                label='Inner Membrane Protein Concentration')
    
    # Add mean and confidence interval
    ax.axhline(y=mean_val, color=COLORS['Navy Blue'], linestyle='--', linewidth=2)
    x_min, x_max = ax.get_xlim()
    fill_x = np.insert(x_data.astype(float).values, 0, x_min)
    fill_x = np.append(fill_x, x_max)
    ax.fill_between(fill_x,
                     [mean_val - conf_interval] * len(fill_x),
                     [mean_val + conf_interval] * len(fill_x),
                     color=COLORS['Navy Blue'], alpha=0.3, label='_nolegend_')
    
    # Customize plot
    ax.set_xlabel(r'Growth rate ($h^{-1}$)', fontsize=FONT_SIZE, fontfamily=FONT_FAMILY)
    ax.set_ylabel(r'Protein Concentration ($g/g_{DCW}$)', fontsize=FONT_SIZE, fontfamily=FONT_FAMILY)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=FONT_SIZE - 2)
    ax.grid(False)
    ax.set_facecolor('white')
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(0, y_data.max() * 1.5)
    
    # Add secondary y-axis showing percentage of total protein
    ax2 = ax.twinx()
    mean_total_concentration = location_summary['Total Protein Concentration'].mean()
    y1_lim = ax.get_ylim()
    ax2.set_ylim(y1_lim[0] / mean_total_concentration * 100, 
                 y1_lim[1] / mean_total_concentration * 100)
    ax2.set_ylabel(r'Protein Fraction of Total (%)', fontsize=FONT_SIZE, fontfamily=FONT_FAMILY)
    
    # Save plot
    output_path = os.path.join(OUTPUT_DIR, 'Figure_1C_schmidt2016_InnerMembProt_Loc_Conc')
    fig.savefig(f'{output_path}.png', format='png', bbox_inches='tight', 
                transparent=False, dpi=PLOT_DPI_PNG)
    fig.savefig(f'{output_path}.svg', format='svg', bbox_inches='tight', 
                transparent=False, dpi=PLOT_DPI_SVG)
    
    plt.show()
    print(f"Saved: Figure_1C_schmidt2016_InnerMembProt_Loc_Conc")
    

if __name__ == "__main__":
    main()
