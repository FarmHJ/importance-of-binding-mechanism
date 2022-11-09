import matplotlib
import numpy as np
import os
import pandas as pd

import modelling

SA_model = modelling.SensitivityAnalysis()

# Read data for space
saved_data_dir = '../../simulation_results/'
file_prefix = 'SA_APD'
result_files = [saved_data_dir + f for f in os.listdir(saved_data_dir)
                if f.startswith(file_prefix)]

first_iter = True
for file in result_files:
    df = pd.read_csv(file,
                     header=[0, 1], index_col=[0],
                     skipinitialspace=True)
    if first_iter:
        combined_df = df
        first_iter = False
    else:
        combined_df = pd.concat([combined_df, df])

combined_df = combined_df.sort_values(
    by=[('param_values', 'Ku'), ('param_values', 'Kmax'),
        ('param_values', 'Vhalf')], ascending=[False, True, True])

Vhalf_range = combined_df['param_values']['Vhalf'].values
Kmax_range = combined_df['param_values']['Kmax'].values
Ku_range = combined_df['param_values']['Ku'].values

RMSError = combined_df['RMSE']['RMSE'].values
MError = combined_df['ME']['ME'].values

nan_ind = [i for i in range(len(RMSError)) if np.isnan(RMSError[i]) or
           np.isnan(MError[i])]
Error_space = RMSError * MError / np.abs(MError)

fig = modelling.figures.FigureStructure(figsize=(10, 3),
                                        gridspec=(1, 3),
                                        wspace=0.3,)
plot = modelling.figures.FigurePlot()
cmap_SD = matplotlib.cm.get_cmap('Oranges')
cmap_CS = matplotlib.cm.get_cmap('Blues')

previous_i = 0

# Choose a random Ku and Kmax
changing_Kmax_ids = [0]
for i in range(1, int(len(Kmax_range))):
    if Kmax_range[i] != Kmax_range[i - 1]:
        changing_Kmax_ids.append(i)
print(changing_Kmax_ids[:7])
chosen_Kmax_id = []
# for choice in [10, 200, 400]:
for choice in [0, 2, 4]:
    chosen_Kmax_id.append((changing_Kmax_ids[choice - 1], changing_Kmax_ids[choice]))

for previous_i, i in chosen_Kmax_id:
    print(previous_i, i)

#         norm = matplotlib.colors.Normalize(0, i - previous_i)

#         for r in range(previous_i, i):
#             APD_trapping = combined_df.iloc[[r]]['APD_trapping'].values[0]
#             APD_conductance = combined_df.iloc[[r]]['APD_conductance'].values[0]

#             fig.axs[0][0].plot(np.arange(len(APD_trapping)),
#                                APD_trapping, 'o', color=cmap_SD(norm(r)))
#             fig.axs[0][0].plot(np.arange(len(APD_conductance)),
#                                APD_conductance, 'o', color=cmap_CS(norm(r)))

#         previous_i = i
#         break

# saved_fig_dir = '../../figures/testing/'
# fig.savefig(saved_fig_dir + 'test_SI.pdf')
