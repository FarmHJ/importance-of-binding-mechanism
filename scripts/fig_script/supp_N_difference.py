# Plot the range of RMSD for each synthetic drug when varying the Hill
# coefficient

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

# Read APD90 differences for all synthetic drug
data_dir = '../../simulation_data/'
filename = 'SA_alldrugs_opt.csv'
drug_df = pd.read_csv(data_dir + filename,
                      header=[0, 1], index_col=[0],
                      skipinitialspace=True)
drug_list = drug_df[('drug', 'drug')].values
drug_list = drug_list[:-1]

# Define directories and variables
data_dir = '../../simulation_data/supp_mat/APD90diff_N/'
percentage_diff_filename = 'N_percentage_diff.csv'
first_iter = True
RMSD_boxplot = []
delta_RMSD_boxplot = []

for drug in drug_list:
    # Read RMSD of each synthetic drug when the Hill coefficient varies
    filepath = data_dir + 'SA_' + drug + '_N.csv'
    df = pd.read_csv(filepath,
                     header=[0, 1], index_col=[0],
                     skipinitialspace=True)

    N_range = np.array(df['param_values']['N'].values)
    RMSD_arr = np.array(df['RMSE']['RMSE'].values)

    RMSD_arr_boxplot = RMSD_arr[~np.isnan(RMSD_arr)]
    RMSD_boxplot.append(RMSD_arr_boxplot)

    # Calculate the ratio of the difference in RMSD
    drug_RMSD = drug_df.loc[drug_df[('drug', 'drug')] == drug][
        ('RMSE', 'RMSE')].values
    delta_RMSD = np.abs(RMSD_arr - drug_RMSD)
    delta_RMSD_ratio = delta_RMSD / drug_RMSD
    df[("RMSD", "deltaRMSD_ratio")] = delta_RMSD_ratio

    # Calculate the summary statistics of the changes in the RMSD
    delta_RMSD = delta_RMSD[~np.isnan(delta_RMSD)]
    delta_RMSD_boxplot.append(delta_RMSD)
    delta_RMSD_stats = [min(delta_RMSD), max(delta_RMSD), np.mean(delta_RMSD)]

    delta_RMSD_ratio_stats = [i / drug_RMSD for i in delta_RMSD_stats]

    # Create dataframe to save results
    delta_RMSD_df = pd.DataFrame(
        [drug] + delta_RMSD_stats + delta_RMSD_ratio_stats,
        index=['drug', 'deltaRMSD_min', 'deltaRMSD_max', 'deltaRMSD_mean',
               'deltaRMSD_ratio_min', 'deltaRMSD_ratio_max',
               'deltaRMSD_ratio_mean'])

    if first_iter:
        combined_df = delta_RMSD_df.T
        first_iter = False
    else:
        combined_df = pd.concat([combined_df, delta_RMSD_df.T])

combined_df.to_csv(data_dir + percentage_diff_filename)

RMSD_mean = np.mean(combined_df['deltaRMSD_mean'].values)

# Plot the distribution of the RMSD for each synthetic drug
fig = plt.figure(figsize=(10, 3))
plt.rcParams.update({'font.size': 9})
ax = fig.add_subplot(1, 2, 1)
ax.boxplot(RMSD_boxplot)
ax.set_xticks(np.arange(12) + 1, labels=drug_list)
ax.set_ylabel('RMSD (ms)')

# Plot a histogram of the changes in RMSD for all synthetic drugs
ax2 = fig.add_subplot(1, 2, 2)
plt.setp(ax.get_xticklabels(), rotation=45, ha='right',
         rotation_mode='anchor')

for i in range(12):
    if i == 0:
        arr = delta_RMSD_boxplot[i]
    else:
        arr = np.concatenate((arr, delta_RMSD_boxplot[i]))
ax2.hist(arr, bins=25)
ax2.set_xlabel('RMSD difference (ms)')
ax2.set_ylabel('Number of virtual drugs')

# Add panel label
fig.text(0.075, 0.9, '(A)', fontsize=11)
fig.text(0.5, 0.9, '(B)', fontsize=11)

# Save figure
fig_dir = '../../figures/supp_mat/APD90diff_N/'
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)
plt.savefig(fig_dir + 'RMSD_N.pdf',
            bbox_inches='tight')

# Show mean and standard deviation of the histogram
print(np.mean(arr))
print(np.std(arr))
