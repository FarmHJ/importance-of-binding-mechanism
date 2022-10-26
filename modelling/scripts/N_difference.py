# To make sure that the hERG current and action potential are not affected by
# the change in Hill's coefficient

import numpy as np
import os
import pandas as pd

saved_data_dir = '../../simulation_data/sensitivity_analysis/N/'

# Read data for drugs
saved_data_dir = '../../simulation_data/'
filename = 'SA_alldrugs.csv'
drug_df = pd.read_csv(saved_data_dir + filename,
                      header=[0, 1], index_col=[0],
                      skipinitialspace=True)

# drug_list = drug_df[('drug', 'drug')].values
drug_list = ['bepridil', 'cisapride', 'dofetilide', 'quinidine',
             'sotalol', 'terfenadine', 'verapamil']

saved_data_dir = '../../simulation_data/sensitivity_analysis/N/'
percentage_diff_filename = 'N_percentage_diff.csv'
for drug in drug_list:
    # Read data
    filepath = saved_data_dir + 'SA_' + drug + '_N.csv'
    df = pd.read_csv(filepath,
                     header=[0, 1], index_col=[0],
                     skipinitialspace=True)

    RMSD_arr = df['RMSE']['RMSE'].values
    drug_RMSD = drug_df.loc[drug_df[('drug', 'drug')] == drug][('RMSE', 'RMSE')]
    delta_RMSD = (RMSD_arr - drug_RMSD) / drug_RMSD
    df[("RMSD", "deltaRMSD_ratio")] = delta_RMSD
    df.to_csv(filepath)

    delta_RMSD_stats = [min(delta_RMSD), max(delta_RMSD), np.mean(delta_RMSD)]

    delta_RMSD_df = pd.DataFrame(
        [drug] + delta_RMSD_stats,
        index=['drug', 'deltaRMSD_min', 'deltaRMSD_max', 'deltaRMSD_mean'])

    if os.path.exists(saved_data_dir + percentage_diff_filename):
        combined_df = pd.read_csv(saved_data_dir + percentage_diff_filename,
                                  header=[0, 1], index_col=[0],
                                  skipinitialspace=True)
        combined_df = pd.concat([combined_df, delta_RMSD_df.T])
    else:
        combined_df = delta_RMSD_df.T

    combined_df.to_csv(saved_data_dir + percentage_diff_filename)
