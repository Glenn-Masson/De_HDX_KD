HDX-MS Data Processing and Kd determination in peptide level

by De Lin, October 2025

Overview

This Python script automates data preprocessing, analysis, and Michaelis–Menten fitting for HDX-MS (Hydrogen–Deuterium Exchange Mass Spectrometry) experiments.
It takes a source CSV file (all result table) exported from HDEaminer analsysis containing peptide-level data, cleans and organizes it, calculates deuterium protection percentages, and performs a weighted Michaelis–Menten fit to estimate Dmax and Kd parameters for each peptide.

The results are exported as structured Excel files and high-quality fit plots.

Features

1. Automated filtering, renaming, and data cleaning
2. Peptide-wise data separation and summary
3. Calculation of deuterium uptake and protection percentage
4. Automatic concentration calculation based on experimental setup
5. Weighted Michaelis–Menten curve fitting using standard deviation
6. Generation of fitted plots (.jpg) for each peptide
7. Export of all intermediate and final results to Excel files

Input Requirements

1. Source Data: A CSV file (e.g., GCN5_4250_300s_202510.csv) containing columns such as:
Column Name	Description
Sequence	Peptide sequence
Charge	        Charge state
Deut Time	Exposure time (e.g., 0s, 300.00s)
Protein State	Experimental state number
Exp Cent	Experimental centroid (m/z)

2. User Inputs

When the script starts, it will prompt for:

please input source file: GCN5_4250_300s_202510
please input compound name: GCN5_4250_300s
please input top concentration in uM: 128


source file: CSV filename (without extension)
compound name: Used to name output folders and result files
top concentration: The highest ligand concentration (μM)

Output Files

All results are saved in a new folder named after your compound (e.g., MyCompound/).

Output File	Description
*_peptide_select.xlsx	Each peptide saved in a separate sheet
*_peptide_D_protection.xlsx	Deuterium protection calculation per peptide
*_fitting_results.xlsx	Michaelis–Menten fit results and plots
fit_*_*.jpg	Curve-fitting plots for each peptide
*_merged_output.xlsx	Combined summary of all peptides
*_merged_output_final.xlsx	Final filtered summary (Dmax > 10, R² > 0.8)

Dependencies

Install required packages before running the script: numpy pandas matplotlib scipy openpyxl xlsxwriter

Output Interpretation

Dmax: Estimated maximum deuterium protection
Dmax_std: Standard deviation of Dmax
Kd: Dissociation constant (μM)
Kd_std: Standard deviation of Kd
R²: Goodness of fit
Final filtered results (*_merged_output_final.xlsx) contain only peptides with:
Dmax > 10
R² > 0.8

Notes
The script assumes exposure times “0s” correspond to undeuterated deuterated conditions.
The concentration calculation is based on dilution steps and experimental constants.


