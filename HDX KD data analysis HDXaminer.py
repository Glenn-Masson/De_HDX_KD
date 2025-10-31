# by de Lin May 2025

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize as opt

source_file = input("please input source file: ")
Comp_name = input("please input compound name: ")
Top_con = int(input("please input top concentration in uM: "))
output_folder = f"{Comp_name}"
os.makedirs(output_folder, exist_ok=True)


df = pd.read_csv(f'{source_file}.csv')
#df = df[~df["Confidence"].str.contains("Low", case=False, na=False)]
df = df.rename(columns={"Charge": "z", "Deut Time": "Exposure", "Protein State": "State", "Exp Cent": "Center"})
# df["Exposure"] = df["Exposure"].astype(str).replace({"0s": "0", "300.00s": "300"}).astype(int)
df["Exposure"] = (
    df["Exposure"]
    .astype(str)
    .str.replace(r"s$", "", regex=True)        # remove s
    .str.replace(r"\.00$", "", regex=True)     # remove .00
    .astype(float)
    .astype(int)
)

df = df[df["Sequence"].str.len() <= 29]

print(df)
df['Sequence_z'] = df['Sequence'].astype(str) + '_' + df['z'].astype(str)
# select column
column_name = "Sequence_z"

# get unique value
unique_values = df[column_name].unique()
output1_path = os.path.join(output_folder, f"{Comp_name}_peptide_select.xlsx")
with pd.ExcelWriter(output1_path, engine="xlsxwriter") as writer:
    # Search the unique values and filter out the corresponding rows
    for value in unique_values:
        filtered_df = df[df[column_name] == value]  # Filter all rows for a specific value
        filtered_df.to_excel(writer, sheet_name=str(value), index=False)  # save to different sheet

# Read all sheets of an Excel file
input_file = output1_path
sheets = pd.read_excel(input_file, sheet_name=None)
output_file = os.path.join(output_folder, f"{Comp_name}_peptide_D_protection.xlsx") # result saved to Excel file

# Read all sheets of an Excel file
xls = pd.ExcelFile(input_file)

# Create an ExcelWriter object
with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
    for sheet_name in xls.sheet_names:
        # Read current sheet
        df = pd.read_excel(xls, sheet_name=sheet_name)

        # Setting parameters
        State_number = 0
        Exposure_number = 0

        # Filter data
        filtered_df_1 = df[~((df['State'] != State_number) & (df['Exposure'] == Exposure_number))]
        filtered_df_2 = filtered_df_1[['Sequence', 'State', 'Exposure', 'z', 'Center']].copy()

        # calculate mz
        filtered_df_2.loc[:, 'mz'] = filtered_df_2['z'] * filtered_df_2['Center']

        # Generate State_expo column
        filtered_df_2.loc[:, 'State_expo'] = filtered_df_2['State'].astype(str) + '_' + filtered_df_2[
            'Exposure'].astype(str)
        result = []

        for state in np.unique(filtered_df_2.State_expo):  # List all the state expo
            result.append(
                [state] + filtered_df_2.loc[
                    filtered_df_2.State_expo == state].mz.tolist())  # find out mz with a specified state, make it as a list

        tb2 = pd.DataFrame(result, columns=["State_expo", "mz1", "mz2", "mz3"])
        tb2.loc[:, 'Sequence'] = sheet_name

        tb2['mz_average'] = tb2[['mz1', 'mz2', 'mz3']].mean(axis=1)
        print(tb2)
        average_value_0_0 = tb2.iloc[0, 1:4].mean()  # Non deteurium
        average_value_0_5 = tb2.iloc[1, 1:4].mean()  # max deteurium
        D_Prot_max = average_value_0_5 - average_value_0_0  # max D uptake
        tb2["D_mz1"] = tb2["mz1"] - average_value_0_0  # D uptake
        tb2["D_mz2"] = tb2["mz2"] - average_value_0_0
        tb2["D_mz3"] = tb2["mz3"] - average_value_0_0
        tb2["D_mz1_prot"] = ((D_Prot_max - tb2["D_mz1"]) / D_Prot_max) * 100  # D uptake protection percentage effect
        tb2["D_mz2_prot"] = ((D_Prot_max - tb2["D_mz2"]) / D_Prot_max) * 100
        tb2["D_mz3_prot"] = ((D_Prot_max - tb2["D_mz3"]) / D_Prot_max) * 100
        tb2['D_prot_average'] = tb2[['D_mz1_prot', 'D_mz2_prot', 'D_mz3_prot']].mean(axis=1)  # average D uptake protection percentage effect
        tb2["D_Prot_std"] = tb2[['D_mz1_prot', 'D_mz2_prot', 'D_mz3_prot']].std(axis=1)
        tb3 = tb2[~tb2["State_expo"].isin(["0_0", "0_5"])]

        tb3 = tb3.copy()  # Make it a true copy if it was derived from slicing
        tb3["stage"] = tb3["State_expo"].str.split("_").str[0]
        tb3["stage"] = pd.to_numeric(tb3["stage"], errors='coerce')
        tb3.loc[:, f"{Comp_name}_conc_uM"] = 2 ** tb3["stage"] * (Top_con * 57 / 60) / 2 ** 8
        tb3 = tb3[["Sequence", f"{Comp_name}_conc_uM", "D_prot_average", "D_Prot_std"]]  # final tabl
        # Save to a new Excel file
        tb3.to_excel(writer, sheet_name=sheet_name, index=False)
print(f"The processed data has been saved to {output_file}")

def michaelis_menten_weighted_fit(df, x_col, y_col, std_col):
    """Fitting a DataFrame to a Michaelis-Menten"""

    X = df[x_col].values
    Y = df[y_col].values
    Y_err = df[std_col].values  # Standard Deviation

    # Defining the Michaelis-Menten formula
    def model(X, Dmax, Kd):
        return (Dmax * X) / (Kd + X)
    try:
        # Initial parameter guess
        initial_guess = [max(Y), np.median(X)]

        # Weighted fit using standard deviation
        params, cov = opt.curve_fit(model, X, Y, p0=initial_guess, sigma=Y_err, absolute_sigma=True)
        Dmax_opt, Kd_opt = params

        # Calculate standard deviations (sqrt of diagonal of covariance matrix)
        Dmax_std, Kd_std = np.sqrt(np.diag(cov))

        # Calculating the R² value
        residuals = Y - model(X, Dmax_opt, Kd_opt)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((Y - np.mean(Y)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)

        return round(Dmax_opt, 0), round(Dmax_std, 0), round(Kd_opt, 1), round(Kd_std, 1), round(r_squared, 3)
    except (RuntimeError, ValueError) as e:
        print(f"Fit failed: {e}")
        return None, None, None, None, None
def process_excel_sheets(input_file1, output_file1, x_col, y_col, std_col):
    """Perform Michaelis-Menten fitting on multiple sheets of an Excel file and save the results"""

    # read all sheet
    sheet_data = pd.read_excel(input_file1, sheet_name=None)

    results = {}

    with pd.ExcelWriter(output_file1, engine="openpyxl") as writer:
        for sheet_name, df in sheet_data.items():
            try:
                # Perform the fitting
                Dmax_opt, Dmax_std, Kd_opt, Kd_std, r_squared = michaelis_menten_weighted_fit(df, x_col, y_col, std_col)

                # Recording Results
                results[sheet_name] = [Dmax_opt, Dmax_std, Kd_opt, Kd_std, r_squared]

                # Generate a fitted curve
                X_fit = np.linspace(df[x_col].min(), df[x_col].max(), 100)
                Y_fit = (Dmax_opt * X_fit) / (Kd_opt + X_fit)

                # Draw a scatter plot + fit a curve
                plt.figure(figsize=(8, 6))
                plt.errorbar(df[x_col], df[y_col], yerr=df[std_col], fmt='o', color='blue', alpha=0.6, label="Data")
                plt.plot(X_fit, Y_fit, color='red', label=f"Fit: Dmax={Dmax_opt:.2f}, Kd={Kd_opt:.2f}")
                plt.xlabel(x_col)
                plt.ylabel(y_col)
                plt.title(f"Michaelis-Menten Fit ({sheet_name})")
                plt.legend()
                plt.grid()
                plt.savefig(os.path.join(output_folder, f"fit_{Comp_name}_{sheet_name}.jpg"))
                plt.close()

                # Save the resulting DataFrame to Excel
                df_results = pd.DataFrame({
                    "Dmax":[Dmax_opt],
                    "Dmax_std": [Dmax_std],
                    "Kd":[Kd_opt],
                    "Kd_std": [Kd_std],
                    "R²":[r_squared]
                })
                df_results.to_excel(writer, sheet_name=f"{sheet_name}_fit", index=False)

            except Exception as e:
                print(f"⚠️ Error processing {sheet_name}: {e}")

    print("✅ Processing completed, fitting results saved!")

# 📌 Calling a function：
input_file1 = output_file  # input file
output_file1 = os.path.join(output_folder, f"{Comp_name}_fitting_results.xlsx") # Excel file with saved results
x_column = f"{Comp_name}_conc_uM"  # X-axis (ligand concentration)
y_column = "D_prot_average"  # Y-axis (binding affinity)
std_column = "D_Prot_std"  # Y-axis standard deviation (error)

process_excel_sheets(input_file1, output_file1, x_column, y_column, std_column)
output_file_1 = os.path.join(output_folder, f"{Comp_name}_merged_output.xlsx") # Excel file with saved results

# read all sheet
xls = pd.ExcelFile(output_file1)
all_sheets = []

# Traverse each sheet and add a sheet name column
for sheet_name in xls.sheet_names:
    df = pd.read_excel(xls, sheet_name=sheet_name)

    # Add the Sheet name in the first column
    df.insert(0, "Sequence", sheet_name)

    # Add DataFrame to List
    all_sheets.append(df)

# Merge all data
merged_df = pd.concat(all_sheets, ignore_index=True)
merged_df_sorted = merged_df.sort_values(by="R²", ascending=False)
# Save to a new Excel file
merged_df_sorted.to_excel(output_file_1, index=False)

df_final = merged_df_sorted[(merged_df_sorted["Dmax"] > 10) & (merged_df_sorted["R²"] > 0.8)]
output_file_2 = os.path.join(output_folder, f"{Comp_name}_merged_output_final.xlsx")
df_final.to_excel(output_file_2, index=False)
print(f"✅ Data merge completed and saved to {output_file_2}")




