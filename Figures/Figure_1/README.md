# Figure 1: Proteome Constraints

This directory contains the scripts and data processing logic to reproduce the experimental data plots shown in Figure 1 (Panels B and C) of the manuscript.

## Contents

-   **`Figure1_Protein_Concentrations.py`**: A Python script that:
    1.  Loads proteomics data from Schmidt et al. (2016).
    2.  Calculates the total protein concentration ($g/g_{DCW}$) at different growth rates (Panel B).
    3.  Calculates the inner membrane protein concentration ($g/g_{DCW}$) at different growth rates (Panel C).
    4.  Generates the plots and saves them to the `output/` directory.

-   **`output/`**: (Generated) Contains the resulting plots:
    -   `Figure_1B_Total_Protein_Concentration.png/svg`
    -   `Figure_1C_InnerMembrane_Concentration.png/svg`

## Data Sources

The script relies on the following input data, located in `../../data/input/`:

1.  **Proteomics Data**: `data/input/Proteomics/Schmidt2016_DataToPlot.xlsx`
    -   Source: Schmidt et al., *Nature Biotechnology*, 2016.
    -   Contains absolute protein quantification (Table S6, S11, S13).

## How to Run

```bash
# From the root of the repository
cd Figures/Figure_1
python Figure1_Protein_Concentrations.py
```
