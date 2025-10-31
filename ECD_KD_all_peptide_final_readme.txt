HDX-MS/MS (ECD) Data Processing and Kd determination in Signal Amino Acid level

by De Lin, October 2025

Overview

This Python script processes **ECD (Electron Capture Dissociation)** data to determine the apparent dissociation constants (Kd) for individual amino acid positions within peptides or proteins.  
It automates data extraction, differential uptake analysis, Michaelis–Menten fitting, visualization, and statistical evaluation for both `c` and `z` fragment ions.

The results are exported as structured Excel files and high-quality fit plots.

Features

Automatic peptide identification and data filtering
Calculation of deuterium protection relative to the parent peptide
Concentration determination based on experimental dilution
Weighted Michaelis–Menten curve fitting with propagated standard deviation
Separate analysis of parent, c, and z fragment ions
Automatic generation of fitted curves (.jpg) and parameter tables (.xlsx)
Combined PDF summary of all fitted figures for reporting

Input Requirements

1. Source Data: A CSV file (e.g., 250414_GCN5_4250_ECD_Kd_stateData_202507.csv -- state data exported from NewMarket analysis) containing columns such as:

Column Name	Description
Sequence	Peptide sequence
Start / End	Residue numbering for each peptide
Fragment	Fragment ion label (e.g., c3, z8)
Uptake	        Deuterium uptake value
Uptake SD	Standard deviation of uptake
State	        Experimental condition index (e.g., 0–8)


2. User Inputs

When the script starts, it will prompt for:

please input source file: 250414_GCN5_4250_ECD_Kd_stateData_202507
please input compound name: DDD02444250_202507
please input top concentration in uM: 128


source file: CSV filename (without extension)
compound name: Used to name output folders and result files
top concentration: The highest ligand concentration (μM)

Output Files

All results are saved in a new folder named after your compound (e.g., DDD02444250_202507/).

Output File	Description
*_data_filtered.xlsx	Cleaned peptide or fragment data
*_Parent_data_filtered.xlsx	Processed parent peptide data
*_parent_Kd_fit.jpg	Michaelis–Menten fit plot for parent peptide
*_parent_Kd_fit_params.xlsx	Fitted parameters for parent peptide
_c_output.xlsx / _z_output.xlsx	Intermediate c or z fragment calculations
*_c_Kd_fit_params.xlsx / *_z_Kd_fit_params.xlsx	Fitted Kd parameters for each fragment ion
*_Kd_fit_params.xlsx	Combined c/z fragment Kd results
Kd_figures.pdf	All fit plots combined into a single PDF

Dependencies

Install required packages before running the script: numpy pandas matplotlib scipy openpyxl xlsxwriter

Output Interpretation

Parameter	Description
Dmax	        Maximum deuterium protection (%)
Dmax_std	Standard deviation of Dmax
Kd	        Apparent dissociation constant (µM)
Kd_std	        Standard deviation of Kd
R²	        Goodness of fit
p_vs_Frag_Parent	Significance vs. parent peptide Kd (t-test)

The combined output (*_Kd_fit_params.xlsx) summarizes residue-level Kd trends derived from both c and z ions.

Notes
The reference (State 0) is treated as the unbound or control condition.
Only fragments with sufficient uptake difference (>0.1) are included in Kd fitting.
If a fitting fails for a fragment, it is skipped and reported in the terminal.
The script automatically combines all .jpg plots into a single PDF for review.


