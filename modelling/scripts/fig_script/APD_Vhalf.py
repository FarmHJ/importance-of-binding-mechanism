import matplotlib
import numpy as np
import os
import pandas as pd

import modelling

SA_model = modelling.SensitivityAnalysis()

# Read data for space
saved_data_dir = '../../simulation_results/SA_space/'
file_prefix = 'SA_allparam'
result_files = [saved_data_dir + f for f in os.listdir(saved_data_dir) if f.startswith(file_prefix)]
# fields = ['param_values', 'RMSE', 'drug_conc_AP', 'APD_trapping', 'APD_conductance']

first_iter = True
for file in result_files:
    df = pd.read_csv(file,
                     header=[0, 1], index_col=[0],
                     skipinitialspace=True)  # , usecols=fields)
    if first_iter:
        combined_df = df
        first_iter = False
    else:
        combined_df = pd.concat([combined_df, df])
        
combined_df = combined_df.sort_values(by=[('param_values', 'Ku'), ('param_values', 'Kmax'), ('param_values', 'Vhalf')],
                                      ascending=[False, True, True])

Vhalf_range = combined_df['param_values']['Vhalf'].values
Kmax_range = combined_df['param_values']['Kmax'].values
Ku_range = combined_df['param_values']['Ku'].values


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
chosen_Kmax_id = []
for choice in [10, 200, 400]:
    chosen_Kmax_id.append((changing_Kmax_ids[choice - 1],
                           changing_Kmax_ids[choice]))

for e, ids in enumerate(chosen_Kmax_id):
    previous_i, i = ids

    norm = matplotlib.colors.Normalize(0, i - previous_i)

    for r in range(previous_i, i):
        drug_conc = combined_df.iloc[[r]]['drug_conc_AP'].values[0]
        APD_trapping = combined_df.iloc[[r]]['APD_trapping'].values[0]
        APD_conductance = combined_df.iloc[[r]]['APD_conductance'].values[0]

        fig.axs[0][e].plot(drug_conc, APD_trapping, 'o',
                           color=cmap_SD(norm(r - previous_i)))
        fig.axs[0][e].plot(drug_conc, APD_conductance, 'o',
                           color=cmap_CS(norm(r - previous_i)))

    fig.axs[0][e].set_title('Kmax = ' + str(Kmax_range[previous_i]))
    fig.axs[0][e].set_xscale("log", nonpositive='clip')
    fig.axs[0][e].set_xlabel('Normalised drug concentration (nM)')
    fig.axs[0][e].set_ylabel(r'APD$_{90}$ (ms)')

saved_fig_dir = '../../figures/testing/'
fig.savefig(saved_fig_dir + 'test_SI.pdf')
