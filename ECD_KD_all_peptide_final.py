# For signal amino acid KD determination from ECD data
# By De Lin 15th May 2025

import os
from PIL import Image
import numpy as np
import pandas as pd
import re
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind_from_stats

source_file = input("please input source_stateData file: ")
Comp_name = input("please input compound name: ")
text = "GHKEVQLKDQILGVLDYLEKQQSAWPFLKPVSLSEAPDYYDIIKEPTDILTMRRKARHGDYKTKEDFGIELKRMFDNCRLYNAPTTIYFKYANELQTLIWPKYEAI"
Top_con = int(input("please input top concentration in uM: "))
output_folder = f"{Comp_name}"
os.makedirs(output_folder, exist_ok=True)
#
df = pd.read_csv(f'{source_file}.csv')
#peptides
subset_df = df[df['Fragment'].isna()][['Sequence', 'Start', 'End']]
peptide_df = subset_df.drop_duplicates()

def michaelis_menten(S, Dmax, Kd):
    return (Dmax * S) / (Kd + S)

for idx, pep_row in peptide_df.iterrows():
    params_df_z = pd.DataFrame()
    params_df_c = pd.DataFrame()
    seq = pep_row['Sequence']
    start = pep_row['Start']
    end = pep_row['End']
    # df['Fragment'] = df['Fragment'].str.split('/').str[1]
    df_filtered_parent = df[(df['Uptake'] != 0) & (df['Fragment'].isna()) & (df['Start'] >= start) & (
                df['Start'] <= end)]
    df['Fragment'] = df['Fragment'].apply(lambda x: x.split('/')[1] if isinstance(x, str) and '/' in x else x)
    df_filtered_parent.loc[:, 'Sequence_Frag'] = df_filtered_parent['Sequence'].astype(str) + '_parent'
    df['Sequence_Frag'] = df['Sequence'].astype(str) + '_' + df['Fragment'].astype(str)
    df_filtered = df[(df['Uptake'] != 0) & (df['Fragment'] != 'null') & (df['Start'] >= start) & (
                df['Start'] <= end)]
    # select column
    column_name = "Sequence_Frag"
    df_filtered = df_filtered.rename(columns={
        'Uptake SD': 'UptakeSD',
    })
    df_filtered_parent = df_filtered_parent.rename(columns={
        'Uptake SD': 'UptakeSD',
    })
    df_filtered_parent['conc_uM'] = 2 ** df_filtered_parent['State'] * (Top_con * 57 / 60) / 2 ** 8
    ref_row = df_filtered_parent[df_filtered_parent['State'] == 0].iloc[0]
    uptake_0 = ref_row['Uptake']
    uptake_0_SD = ref_row['UptakeSD']

    df_filtered_parent['Frag_Parent'] = (uptake_0 - df_filtered_parent['Uptake']) / uptake_0 * 100
    df_filtered_parent['Frag_Parent_SD'] = 100 * np.sqrt(
                    (df_filtered_parent['UptakeSD'] / uptake_0 ** 2 * uptake_0_SD) ** 2 +
                    (df_filtered_parent['UptakeSD'] / uptake_0) ** 2
                )
    df_filtered_parent = df_filtered_parent[df_filtered_parent['State'] != 0]
    df_filtered = df_filtered.dropna(subset=['Fragment'])
    print(df_filtered)


    output_filename = os.path.join(output_folder, f"{Comp_name}_{seq}_data_filtered.xlsx")
    df_filtered.to_excel(output_filename, index=False)
    output_filename_parent = os.path.join(output_folder, f"{Comp_name}_{seq}_Parent_data_filtered.xlsx")
    df_filtered_parent.to_excel(output_filename_parent, index=False)


    def round_sig(x, sig=3):
        if pd.isnull(x):
            return x
        if x == 0:
            return 0
        return round(x, sig - int(np.floor(np.log10(abs(x)))) - 1)


    # Customize the number of significant digits retained in each column
    sig_figs = {
        'Dmax': 2,
        'Dmax_std': 1,
        'Kd': 2,
        'Kd_std': 1,
        'R2': 3
    }

    xdata = df_filtered_parent['conc_uM']
    ydata = df_filtered_parent['Frag_Parent']
    yerr = df_filtered_parent['Frag_Parent_SD']

    # Fitting
    # params, _ = curve_fit(michaelis_menten, xdata, ydata, sigma=yerr, absolute_sigma=True,
    #                       p0=[1, 0.1])
    popt, pcov = curve_fit(michaelis_menten, xdata, ydata, sigma=yerr, absolute_sigma=True, bounds=(0, np.inf))
    Dmax, Kd = popt

    # Calculate the standard deviation (by taking the square root of the diagonal elements of the covariance matrix).
    perr = np.sqrt(np.diag(pcov))
    Dmax_std, Kd_std = perr

    # Calculate R² (coefficient of determination)
    y_pred = michaelis_menten(xdata, Dmax, Kd)
    ss_res = np.sum((ydata - y_pred) ** 2)
    ss_tot = np.sum((ydata - np.mean(ydata)) ** 2)
    r2 = 1 - (ss_res / ss_tot)

    params_dict = {}

    params_dict['Frag_Parent'] = {'Dmax': Dmax, 'Dmax_std': Dmax_std, 'Kd': Kd, 'Kd_std': Kd_std, 'R2': r2}

    # Curve fitting
    S_fit = np.linspace(0, xdata.max(), 100)
    v_fit = michaelis_menten(S_fit, Dmax, Kd)

    # Plotting
    plt.errorbar(
        xdata,
        ydata,
        yerr=yerr,
        fmt='o',
        label=f'Frag_Parent_data'
    )

    plt.plot(
        S_fit,
        v_fit,
        label=f'Frag_Parent_fit (Dmax={Dmax:.2f}, Kd={Kd:.2f}, R2={r2:.2f})'
    )

    # 5. Improve the appearance of the plot
    plt.xlabel('Substrate Concentration [S] (µM)')
    plt.ylabel('D_Protection (%)')
    plt.title(f"{Comp_name}_{seq}_parent_Fit")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=6)
    plt.tight_layout()
    plot_path1 = os.path.join(output_folder, f"{Comp_name}_{seq}_parent_Kd_fit.jpg")
    plt.savefig(plot_path1, dpi=300, bbox_inches='tight')
    # plt.show()
    plt.close()

    # 6. (Optional) Save parameter results as a DataFrame.
    params_df = pd.DataFrame.from_dict(params_dict, orient='index').reset_index()
    params_df.rename(columns={'index': 'Fragment'}, inplace=True)

    for col, sig in sig_figs.items():
        params_df.loc[:, col] = params_df[col].apply(lambda x: round_sig(x, sig))

    plot_path_parent = os.path.join(output_folder, f"{Comp_name}_{seq}_parent_Kd_fit_params.xlsx")
    params_df.to_excel(plot_path_parent, index=False)



    # 1. Pivot table operation: set 'Sequence_Frag' as rows, 'State' as columns, and 'Uptake' and 'Uptake SD' as values
    pivot_df_filtered = df_filtered.pivot(index='Sequence_Frag', columns='State', values=['Uptake', 'UptakeSD'])
    pivot_df_filtered.columns = [
        f"{v}_{s}" if v != '' else s
        for v, s in pivot_df_filtered.columns
    ]
    # 2. Reset the index to turn 'Sequence_Frag' back into a regular column
    pivot_df_filtered = pivot_df_filtered.reset_index()

    # 3. Create a copy df1, and extract the 'Fragment' part from the 'Sequence_Frag' column.
    df1 = pivot_df_filtered.copy()
    df1['Fragment'] = df1['Sequence_Frag'].str.split('_').str[1]

    # # 4. To reorder the columns so that 'Fragment' appears right after 'Sequence_Frag'
    cols = list(df1.columns)
    cols.insert(cols.index('Sequence_Frag') + 1, cols.pop(cols.index('Fragment')))
    df1 = df1[cols]

    target_col = 'Fragment'

    df1['Position'] = df1[target_col].astype(str).apply(
        lambda x: ' '.join(re.findall(r'\d+', x))
    )
    df1['Position'] = pd.to_numeric(df1['Position'], errors='coerce')
    df1.columns = df1.columns.astype(str)
    chars = ['c', 'z']
    try:
        output_filename1 = os.path.join(output_folder, f"{Comp_name}_{seq}_pivot_df_filtered.xlsx")
        df1.to_excel(output_filename1, index=False)
    except Exception as e:
        print(f"Error writing Excel file: {e}")


    for char in chars:
        count = df1[df1['Fragment'].astype(str).str.contains(char)].shape[0]
        if count > 1:
            df2 = df1[df1[target_col].str.contains(char, case=False, na=False)]
            df2_sorted = df2.sort_values(by='Position', ascending=True).reset_index(drop=True)
            df2_sorted['Frag'] = 'Frag' + '_' + df2_sorted['Fragment'] + '_' + df2_sorted['Fragment'].shift(1)

            # remove leading underscores from the first row
            df2_sorted['Frag'] = df2_sorted['Frag'].str.strip('_')

            columns_to_calculate = ['Uptake_0', 'Uptake_1', 'Uptake_2', 'Uptake_3', 'Uptake_4', 'Uptake_5', 'Uptake_6',
                                    'Uptake_7', 'Uptake_8']
            columns_to_calculateSD = ['UptakeSD_0', 'UptakeSD_1', 'UptakeSD_2', 'UptakeSD_3', 'UptakeSD_4',
                                      'UptakeSD_5',
                                      'UptakeSD_6',
                                      'UptakeSD_7', 'UptakeSD_8']

            # To compute the difference between the current row and the previous row starting from the third row
            for col in columns_to_calculate:
                # To name the new column as diff_ plus the original column name
                diff_col = f'diff_{col}'

                # To calculate the row-wise difference starting from the third row (index 2)
                df2_sorted[diff_col] = df2_sorted[col] - df2_sorted[col].shift(1)

                # First two rows as NaN
                df2_sorted.loc[0, diff_col] = None

            for col in columns_to_calculateSD:
                # To name the new column as diff_ plus the original column name
                diffSD_col = f'diff_{col}'

                #  calculate the row-wise difference starting from the third row (index 2)
                df2_sorted[diffSD_col] = np.sqrt(df2_sorted[col] ** 2 + df2_sorted[col].shift(1) ** 2)

                # First two rows as NaN
                df2_sorted.loc[0, diffSD_col] = None

            output_filename2 = os.path.join(output_folder, f"{Comp_name}_{seq}_rows_with_{char}_Frag.xlsx")
            df2_sorted.to_excel(output_filename2, index=False)


            # def michaelis_menten(S, Dmax, Kd):
            #     return (Dmax * S) / (Kd + S)
            print(f"\n🔍 Processing: {Comp_name}_{seq}_{char}")
            df3 = df2_sorted.copy()
            # difference calculation
            columns_to_subtract = ['diff_Uptake_1', 'diff_Uptake_2', 'diff_Uptake_3', 'diff_Uptake_4', 'diff_Uptake_5',
                                   'diff_Uptake_6', 'diff_Uptake_7', 'diff_Uptake_8']
            columns_to_subtractSD = ['diff_UptakeSD_1', 'diff_UptakeSD_2', 'diff_UptakeSD_3', 'diff_UptakeSD_4',
                                     'diff_UptakeSD_5',
                                     'diff_UptakeSD_6', 'diff_UptakeSD_7', 'diff_UptakeSD_8']
            base_col = 'diff_Uptake_0'
            base_colSD = 'diff_UptakeSD_0'

            df3['diff_Uptake_8_0'] = df3[base_col] - df3['diff_Uptake_8']

            for col in columns_to_subtract:
                new_col = f'Delta_{col}'
                df3[new_col] = (df3[base_col] - df3[col]) / df3[base_col] * 100

            for col in columns_to_subtractSD:
                new_col = f'Delta_{col}'
                col1 = col.replace('SD', '')
                df3[new_col] = 100 * np.sqrt(
                    (df3[col1] / df3[base_col] ** 2 * df3[base_colSD]) ** 2 +
                    (df3[col] / df3[base_col]) ** 2
                )

            # save difference result
            output_filename3 = os.path.join(output_folder, f"{Comp_name}_{seq}_{char}_output1.xlsx")
            df3.to_excel(output_filename3, index=False)

            # Filter + Transpose
            filtered_df = df3[df3['diff_Uptake_8_0'] > 0.1].copy()
            # filtered_df = df3.copy()
            columns_to_transpose = [f"Delta_diff_Uptake_{i}" for i in range(1, 9)] + [f"Delta_diff_UptakeSD_{i}" for i
                                                                                      in range(1, 9)]
            transposed = filtered_df[columns_to_transpose].transpose()
            transposed.columns = filtered_df['Frag'].values
            transposed = transposed.reset_index().rename(columns={'index': 'Sample'})

            transposed['State'] = transposed['Sample'].apply(
                lambda x: pd.to_numeric(re.findall(r'\d+', x)[0] if re.findall(r'\d+', x) else None, errors='coerce'))
            transposed['conc_uM'] = 2 ** transposed['State'] * (Top_con * 57 / 60) / 2 ** 8

            # Save Transposed Table
            output_filename4 = os.path.join(output_folder, f"{Comp_name}_{seq}_{char}_output2.xlsx")
            transposed.to_excel(output_filename4, index=False)

            # Separate the data of Delta_diff_Uptake and Delta_diff_UptakeSD
            df_uptake = transposed.iloc[:8, :].copy()  # The first 8 rows are Delta_diff_Uptake
            df_sd = transposed.iloc[8:, :].copy()  # The last 8 rows are Delta_diff_Uptake
            df_sd = df_sd.rename(columns={col: f"{col}_SD" for col in df_sd.columns})
            df_sd = df_sd.rename(columns={'State_SD': 'State', 'conc_uM_SD': 'conc_uM'})

            df_merged = pd.merge(df_uptake, df_sd, on=['State', 'conc_uM'], how='inner')
            df_merged = df_merged.drop(columns=['Sample_SD'])
            # Specify columns to put in front
            priority_cols = ['Sample', 'State', 'conc_uM']

            # the remaining columns sorted alphabetically
            other_cols = sorted([col for col in df_merged.columns if col not in priority_cols])

            # Reorder the columns
            df_merged = df_merged[priority_cols + other_cols]

            df_merged_parent = df_merged.merge(df_filtered_parent[['State', 'Frag_Parent', 'Frag_Parent_SD']], on='State', how='left')

            df_merged = df_merged_parent.copy()

            output_filename_combined = os.path.join(output_folder, f"{Comp_name}_{seq}_{char}_output3.xlsx")
            df_merged.to_excel(output_filename_combined, index=False)

            required_cols = {'Sample', 'State', 'conc_uM'}
            other_cols = set(df_merged.columns) - required_cols

            if other_cols:
                # 1. Select data column and standard deviation column
                value_cols = [col for col in df_merged.columns if col.startswith('Frag') and not col.endswith('_SD')]
                params_dict = {}

                for col in value_cols:
                    sd_col = f"{col}_SD"

                    # Check whether the corresponding SD column exists
                    if sd_col not in df_merged.columns:
                        print(f"⚠️ Missing standard deviation column: {sd_col}")
                        continue

                    try:
                        xdata = df_merged['conc_uM']
                        ydata = df_merged[col]
                        yerr = df_merged[sd_col]

                        # Fitting
                        # params, _ = curve_fit(michaelis_menten, xdata, ydata, sigma=yerr, absolute_sigma=True,
                        #                       p0=[1, 0.1])
                        popt, pcov = curve_fit(michaelis_menten, xdata, ydata, sigma=yerr, absolute_sigma=True, bounds=(0, np.inf))
                        Dmax, Kd = popt

                        # Calculate the standard deviation (by taking the square root of the diagonal elements of the covariance matrix).
                        perr = np.sqrt(np.diag(pcov))
                        Dmax_std, Kd_std = perr

                        # Calculate R² (coefficient of determination)
                        y_pred = michaelis_menten(xdata, Dmax, Kd)
                        ss_res = np.sum((ydata - y_pred) ** 2)
                        ss_tot = np.sum((ydata - np.mean(ydata)) ** 2)
                        r2 = 1 - (ss_res / ss_tot)

                        params_dict[col] = {'Dmax': Dmax, 'Dmax_std': Dmax_std,'Kd': Kd, 'Kd_std': Kd_std, 'R2': r2}

                        # Curve fitting
                        S_fit = np.linspace(0, xdata.max(), 100)
                        v_fit = michaelis_menten(S_fit, Dmax, Kd)

                        # Plotting
                        plt.errorbar(
                            xdata,
                            ydata,
                            yerr=yerr,
                            fmt='o',
                            label=f'{col} data'
                        )

                        plt.plot(
                            S_fit,
                            v_fit,
                            label=f'{col} fit (Dmax={Dmax:.2f}, Kd={Kd:.2f}, R2={r2:.2f})'
                        )

                    except Exception as e:
                        print(f"❌ Fitting failed for column {col}: {e}")

                # 5. Improve the appearance of the plot
                plt.xlabel('Substrate Concentration [S] (µM)')
                plt.ylabel('D_Protection (%)')
                plt.title(f"{Comp_name}_{seq}_{char}_Frag_Fit")
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=6)
                plt.tight_layout()
                plot_path1 = os.path.join(output_folder, f"{Comp_name}_{seq}_{char}_Kd_fit.jpg")
                plt.savefig(plot_path1, dpi=300, bbox_inches='tight')
                # plt.show()
                plt.close()

                # 6. (Optional) Save parameter results as a DataFrame.
                params_df = pd.DataFrame.from_dict(params_dict, orient='index').reset_index()
                params_df.rename(columns={'index': 'Fragment'}, inplace=True)

                # 1. Get reference row（Fragment == 'Frag_Parent'）
                ref_row = params_df[params_df['Fragment'] == 'Frag_Parent'].iloc[0]
                ref_mean = ref_row['Kd']
                ref_std = ref_row['Kd_std']
                ref_n = ref_row.get('n', 3)  # 假设重复数为 3

                # 2. t-test on all remaining rows
                p_values = []

                for _, row in params_df.iterrows():
                    if row['Fragment'] == 'Frag_Parent':
                        p_values.append(np.nan)  # no statistical test on self vs. self
                        continue

                    mean = row['Kd']
                    std = row['Kd_std']
                    n = row.get('n', 3)

                    _, p = ttest_ind_from_stats(
                        mean1=ref_mean, std1=ref_std, nobs1=ref_n,
                        mean2=mean, std2=std, nobs2=n,
                        equal_var=False
                    )
                    p_values.append(p)

                params_df['p_vs_Frag_Parent'] = p_values
                params_df['Dmax_std'] = params_df['Dmax_std'].round(1)
                params_df['R2'] = params_df['R2'].round(2)

                pos = text.index(seq)
                length = len(seq)

                if params_df['Fragment'].astype(str).str.startswith('Frag_c').any():
                    dfc = params_df.copy()
                    dfc = dfc[dfc['Fragment'] != 'Frag_Parent']
                    # Extract the number after the second single quote
                    dfc['End'] = dfc['Fragment'].str.extract(r"^[^']*'[^']*'(\d+)")

                    # Extract the number after the fourth single quote
                    dfc['Start'] = dfc['Fragment'].str.extract(r"^[^']*'[^']*'[^']*'[^']*'(\d+)")

                    # Convert to integer type
                    dfc['Start'] = dfc['Start'].astype(int)
                    dfc['End'] = dfc['End'].astype(int)

                    # Add two columns n (start AA) and m (number of AA)
                    dfc['n'] = dfc['Start'] + pos + 2
                    dfc['m'] = dfc['End'] - dfc['Start']


                    dfc['extracted'] = dfc.apply(lambda row: text[row['n'] - 1: row['n'] - 1 + row['m']], axis=1)
                    dfc['AA'] = dfc['extracted'].str[0] + dfc['n'].astype(str) + dfc['extracted'].str[1:]
                    cols = [col for col in dfc.columns if col != 'Kd'] + ['Kd']
                    dfc = dfc[cols]

                    params_df = params_df[params_df['Fragment'] == 'Frag_Parent']
                    params_df = params_df[params_df.columns]  # Rearrange the columns in the order of df1

                    # Append to the end of df1 (do not change df2)
                    df_combined = pd.concat([dfc, params_df], ignore_index=True)
                    df_combined['Dmax'] = df_combined['Dmax'].astype(int)
                    df_combined['Kd'] = df_combined['Kd'].round(1)
                    df_combined['Kd_std'] = df_combined['Kd_std'].round(1)
                    df_combined['p_vs_Frag_Parent'] = df_combined['p_vs_Frag_Parent'].round(3)
                    desired_order = ['Kd', 'Kd_std', 'p_vs_Frag_Parent']
                    other_cols = [col for col in dfc.columns if col not in desired_order]
                    df_combined = df_combined[other_cols + desired_order]
                    params_df_c = df_combined.sort_values(by='n', ascending=True)
                    plot_path2 = os.path.join(output_folder, f"{Comp_name}_{seq}_{char}_Kd_fit_params.xlsx")
                    params_df_c.to_excel(plot_path2, index=False)

                if params_df['Fragment'].astype(str).str.startswith('Frag_z').any():
                    dfz = params_df.copy()
                    dfz = dfz[dfz['Fragment'] != 'Frag_Parent']
                    # Extract the number after the first single quote
                    dfz['End'] = dfz['Fragment'].str.extract(r"^[^']*'(\d+)")

                    # Extract the number after the second single quote
                    dfz['Start'] = dfz['Fragment'].str.extract(r"^[^']*'[^']*'(\d+)")

                    # Convert to integer type
                    dfz['Start'] = dfz['Start'].astype(int)
                    dfz['End'] = dfz['End'].astype(int)

                    # Add two columns n (start AA) and m (number of AA)
                    dfz['n'] = pos + length + 1 - dfz['Start']
                    dfz['m'] = dfz['End'] - dfz['Start']

                    dfz['extracted'] = dfz.apply(lambda row: text[row['n'] - 1: row['n'] - 1 + row['m']], axis=1)
                    dfz['AA'] = dfz['extracted'].str[0] + dfz['n'].astype(str) + dfz['extracted'].str[1:]
                    cols = [col for col in dfz.columns if col != 'Kd'] + ['Kd']
                    dfz = dfz[cols]
                    dfz['Dmax'] = dfz['Dmax'].astype(int)
                    dfz['Kd'] = dfz['Kd'].round(1)

                    params_df = params_df[params_df['Fragment'] == 'Frag_Parent']
                    params_df = params_df[params_df.columns]

                    # Append to the end of df1 (do not change df2)
                    df_combined = pd.concat([dfz, params_df], ignore_index=True)

                    df_combined['Dmax'] = df_combined['Dmax'].astype(int)
                    df_combined['Kd'] = df_combined['Kd'].round(1)
                    df_combined['Kd_std'] = df_combined['Kd_std'].round(1)
                    df_combined['p_vs_Frag_Parent'] = df_combined['p_vs_Frag_Parent'].round(3)
                    desired_order = ['Kd', 'Kd_std', 'p_vs_Frag_Parent']
                    other_cols = [col for col in dfz.columns if col not in desired_order]
                    df_combined = df_combined[other_cols + desired_order]
                    params_df_z = df_combined.sort_values(by='n', ascending=True)
                    plot_path2 = os.path.join(output_folder, f"{Comp_name}_{seq}_{char}_Kd_fit_params.xlsx")
                    params_df_z.to_excel(plot_path2, index=False)

                print(f"✅ Finish processing {Comp_name}_{seq}_{char}")
            else:
                print(f"✅ Finish processing {Comp_name}_{seq}_{char}")
    if isinstance(params_df_c, pd.DataFrame) and not params_df_c.empty and isinstance(params_df_z, pd.DataFrame) and not params_df_z.empty:
        combined_df_cz = pd.concat([params_df_c, params_df_z], ignore_index=True)
        plot_path2 = os.path.join(output_folder, f"{Comp_name}_{seq}_Kd_fit_params.xlsx")
        combined_df_cz.to_excel(plot_path2, index=False)
jpg_files = [f for f in os.listdir(output_folder) if f.lower().endswith('.jpg')]
jpg_files.sort()  # Optional: sort by filename.

# Open all images and convert them to RGB mode
images = [Image.open(os.path.join(output_folder, f)).convert('RGB') for f in jpg_files]

# Save as a PDF file
if images:
    output_path = os.path.join(output_folder, 'Kd_figures.pdf')
    images[0].save(output_path, save_all=True, append_images=images[1:])

print("PDF finish combination：", output_path)
