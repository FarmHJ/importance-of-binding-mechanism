# To make sure that the hERG current and action potential are not affected by
# the change in Hill's coefficient

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

saved_data_dir = '../../simulation_data/sensitivity_analysis/N/'

# Read data for drugs
saved_data_dir = '../../simulation_data/'
filename = 'SA_alldrugs.csv'
drug_df = pd.read_csv(saved_data_dir + filename,
                      header=[0, 1], index_col=[0],
                      skipinitialspace=True)

drug_list = drug_df[('drug', 'drug')].values
# drug_list = ['cisapride', 'dofetilide', 'quinidine',
#              'sotalol', 'terfenadine', 'verapamil']

saved_data_dir = '../../simulation_data/sensitivity_analysis/N/'
percentage_diff_filename = 'N_percentage_diff.csv'
first_iter = True
RMSD_boxplot = []
delta_RMSD_boxplot = []

# plt.figure()
for drug in drug_list:
    # print(drug)
    # Read data
    filepath = saved_data_dir + 'SA_' + drug + '_N_2.csv'
    df = pd.read_csv(filepath,
                     header=[0, 1], index_col=[0],
                     skipinitialspace=True)

    N_range = np.array(df['param_values']['N'].values)
    RMSD_arr = np.array(df['RMSE']['RMSE'].values)

    RMSD_arr_boxplot = RMSD_arr[~np.isnan(RMSD_arr)]
    RMSD_boxplot.append(RMSD_arr_boxplot)
#     sort_ind = [i[0] for i in sorted(enumerate(N_range), key=lambda x:x[1])]
#     N_range = sorted(N_range)
#     RMSD_arr = [RMSD_arr[i] for i in sort_ind]
#     plt.plot(N_range, RMSD_arr)

    drug_RMSD = drug_df.loc[drug_df[('drug', 'drug')] == drug][
        ('RMSE', 'RMSE')].values
    delta_RMSD = np.abs(RMSD_arr - drug_RMSD)
    delta_RMSD_ratio = delta_RMSD / drug_RMSD
    df[("RMSD", "deltaRMSD_ratio")] = delta_RMSD_ratio
    # df.to_csv(filepath)

    delta_RMSD = delta_RMSD[~np.isnan(delta_RMSD)]
    delta_RMSD_boxplot.append(delta_RMSD)
    delta_RMSD_stats = [min(delta_RMSD), max(delta_RMSD), np.mean(delta_RMSD)]

    delta_RMSD_ratio_stats = [i / drug_RMSD for i in delta_RMSD_stats]

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

combined_df.to_csv(saved_data_dir + 'SA_N_stats.csv')
combined_df.to_csv(saved_data_dir + percentage_diff_filename)

RMSD_mean = np.mean(combined_df['deltaRMSD_mean'].values)

fig = plt.figure(figsize=(10, 3))
plt.rcParams.update({'font.size': 9})
ax = fig.add_subplot(1, 2, 1)
ax.boxplot(RMSD_boxplot)
ax.set_xticks(np.arange(12) + 1, labels=drug_list)
ax.set_ylabel('RMSD (ms)')

ax2 = fig.add_subplot(1, 2, 2)
# ax2.boxplot(delta_RMSD_boxplot)
# ax2.set_xticks(np.arange(12) + 1, labels=drug_list)
# ax2.set_ylabel('RMSD difference (ms)')

plt.setp(ax.get_xticklabels(), rotation=45, ha='right',
         rotation_mode='anchor')
# plt.setp(ax2.get_xticklabels(), rotation=45, ha='right',
#          rotation_mode='anchor')
# plt.savefig('../../figures/testing/APD90_diff.pdf',
#             bbox_inches='tight')

for i in range(12):
    if i == 0:
        arr = delta_RMSD_boxplot[i]
    else:
        arr = np.concatenate((arr, delta_RMSD_boxplot[i]))
ax2.hist(arr, bins=25)
ax2.set_xlabel('RMSD difference (ms)')
ax2.set_ylabel('Number of virtual drugs')

fig.text(0.075, 0.9, '(A)', fontsize=11)
fig.text(0.5, 0.9, '(B)', fontsize=11)

plt.savefig('../../figures/sensitivity_analysis/N/RMSD_N.pdf',
            bbox_inches='tight')

print(np.mean(arr))
print(np.std(arr))
