import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

import modelling

SA_model = modelling.SensitivityAnalysis()

# # Read data for space
# saved_data_dir = '../../simulation_results/SA_space/'
# file_prefix = 'SA_allparam'
# result_files = [saved_data_dir + f for f in os.listdir(saved_data_dir)
#                 if f.startswith(file_prefix)]
# # fields = ['param_values', 'RMSE', 'drug_conc_AP',
# #           'APD_trapping', 'APD_conductance']

# first_iter = True
# for file in result_files:
#     df = pd.read_csv(file,
#                      header=[0, 1], index_col=[0],
#                      skipinitialspace=True)  # , usecols=fields)
#     if first_iter:
#         combined_df = df
#         first_iter = False
#     else:
#         combined_df = pd.concat([combined_df, df])

# combined_df = combined_df.sort_values(by=[('param_values', 'Ku'),
#                                           ('param_values', 'Kmax'),
#                                           ('param_values', 'Vhalf')],
#                                       ascending=[False, True, True])

# Vhalf_range = combined_df['param_values']['Vhalf'].values
# Kmax_range = combined_df['param_values']['Kmax'].values
# Ku_range = combined_df['param_values']['Ku'].values
# RMSE_range = combined_df['RMSE']['RMSE'].values

# fig = modelling.figures.FigureStructure(figsize=(10, 3),
#                                         gridspec=(1, 3),
#                                         wspace=0.3,)
# fig2 = modelling.figures.FigureStructure(figsize=(10, 2.5),
#                                          gridspec=(1, 3),
#                                          wspace=0.3,)
# plot = modelling.figures.FigurePlot()
# cmap_SD = matplotlib.cm.get_cmap('Oranges')
# cmap_CS = matplotlib.cm.get_cmap('Blues')

# previous_i = 0

# # Choose a random Ku and Kmax
# changing_Kmax_ids = [0]
# for i in range(1, int(len(Kmax_range))):
#     if Kmax_range[i] != Kmax_range[i - 1]:
#         changing_Kmax_ids.append(i)
# chosen_Kmax_id = []
# for choice in [10, 200, 400]:
#     chosen_Kmax_id.append((changing_Kmax_ids[choice - 1],
#                            changing_Kmax_ids[choice]))

# for e, ids in enumerate(chosen_Kmax_id):
#     previous_i, i = ids

#     norm = matplotlib.colors.Normalize(0, i - previous_i)

#     Vhalf_chosen_id = [previous_i, (i - previous_i) / 2, i - 1]

#     for r in range(previous_i, i):
#         drug_conc = combined_df.iloc[[r]]['drug_conc_AP'].values[0]
#         APD_trapping = combined_df.iloc[[r]]['APD_trapping'].values[0]
#         APD_conductance = combined_df.iloc[[r]]['APD_conductance'].values[0]

#         fig.axs[0][e].plot(drug_conc, APD_trapping, 'o',
#                            color=cmap_SD(norm(r - previous_i)))
#         fig.axs[0][e].plot(drug_conc, APD_conductance, 'o',
#                            color=cmap_CS(norm(r - previous_i)))

#     fig.axs[0][e].set_title(r'$(K_u, K_{max}) = $' + "(%.2e, %.2e)" %
#                             (Ku_range[previous_i], Kmax_range[previous_i]))
#     fig.axs[0][e].set_xscale("log", nonpositive='clip')
#     fig.axs[0][e].set_xlabel('Normalised drug concentration')
#     fig.axs[0][e].set_ylabel(r'APD$_{90}$ (ms)')

# saved_fig_dir = '../../figures/testing/'
# fig.savefig(saved_fig_dir + 'test_SI.pdf')

# previous_i, i = chosen_Kmax_id[1]
# norm = matplotlib.colors.Normalize(0, i - previous_i)
# Vhalf_chosen_id = [previous_i, int((i - previous_i) / 2 + previous_i), i - 1]

# for i in range(len(Vhalf_chosen_id)):
#     for num, r in enumerate(Vhalf_chosen_id):
#         drug_conc = combined_df.iloc[[r]]['drug_conc_AP'].values[0]
#         APD_trapping = combined_df.iloc[[r]]['APD_trapping'].values[0]
#         APD_conductance = combined_df.iloc[[r]]['APD_conductance'].values[0]
#         EAD_marker = [1050 if (i >= 1000 or j >= 1000) else None for (i, j)
#                       in zip(APD_trapping, APD_conductance)]

#         if num != i:
#             fig2.axs[0][i].plot(drug_conc, APD_trapping, 'o',
#                                 color='grey', alpha=0.5, zorder=-10)
#             fig2.axs[0][i].plot(drug_conc, APD_conductance, '^',
#                                 color='grey', alpha=0.5, zorder=-10)
#         else:
#             fig2.axs[0][i].plot(drug_conc, APD_trapping, 'o',
#                                 color='orange', label='SD model')
#             fig2.axs[0][i].plot(drug_conc, APD_conductance, '^',
#                                 color='blue', label='CS model')
#             fig2.axs[0][i].scatter(drug_conc, EAD_marker, color='k',
#                                    marker=(5, 2), label='EAD-like AP')

#         print('RMSE = ', RMSE_range[r])
#         fig2.axs[0][i].set_title(r'$V_{half-trap} = $' + "%.2e" %
#                                  (Vhalf_range[Vhalf_chosen_id[i]]))
#         fig2.axs[0][i].set_xscale("log", nonpositive='clip')
#         fig2.axs[0][i].set_xlabel('Normalised drug concentration (nM)')
#         fig2.axs[0][i].set_ylabel(r'APD$_{90}$ (ms)')

# fig2.axs[0][0].legend(handlelength=1)
# fig2.fig.text(0.095, 0.9, '(A)', fontsize=11)
# fig2.fig.text(0.38, 0.9, '(B)', fontsize=11)
# fig2.fig.text(0.66, 0.9, '(C)', fontsize=11)

# saved_fig_dir = '../../figures/sensitivity_analysis/'
# fig2.savefig(saved_fig_dir + 'APD_Vhalf.pdf')

#################################

# Read data
saved_data_dir = '../../simulation_results/'
file_prefix = 'SA_APD'
result_files = [saved_data_dir + f for f in os.listdir(saved_data_dir)
                if f.startswith(file_prefix)]
# fields = [('param_values', 'Vhalf'), ('param_values', 'Kmax'),
#           ('param_values', 'Ku'), ('RMSE', 'RMSE'), ('ME', 'ME')]
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

# combined_df = combined_df[fields]
# print(combined_df.keys())
# combined_df = combined_df.sort_values(by=[('param_values', 'Vhalf')],
#                                       ascending=[True])

nan_df = pd.read_csv(saved_data_dir + 'SA_space/filling_nan.csv',
                     header=[0, 1], index_col=[0],
                     skipinitialspace=True)
combined_df = pd.concat([combined_df, nan_df])
# print(combined_df.keys())

RMSError = combined_df['RMSE']['RMSE'].values
MError = combined_df['ME']['ME'].values

nan_ind = [i for i in range(len(RMSError)) if np.isnan(RMSError[i])
           or np.isnan(MError[i])]
Error_fullspace = RMSError * MError / np.abs(MError)
cmin = min(Error_fullspace)
cmax = max(Error_fullspace)

fig = modelling.figures.FigureStructure(figsize=(10, 3),
                                        gridspec=(1, 3),
                                        wspace=0.1,)
plot = modelling.figures.FigurePlot()

cmap = plt.get_cmap('RdYlBu_r')

chosen_Vhalf_value = [-219.45, -100, -1.147]
for i in range(3):
    chosen_Vhalf_df = combined_df[np.abs(combined_df[('param_values', 'Vhalf')]
                                  - chosen_Vhalf_value[i]) < 0.01]

    # Vhalf_range = chosen_Vhalf_df['param_values']['Vhalf'].values
    Kmax_range = chosen_Vhalf_df['param_values']['Kmax'].values
    Ku_range = chosen_Vhalf_df['param_values']['Ku'].values

    RMSError = chosen_Vhalf_df['RMSE']['RMSE'].values
    MError = chosen_Vhalf_df['ME']['ME'].values

    nan_ind = [i for i in range(len(RMSError)) if np.isnan(RMSError[i])
               or np.isnan(MError[i])]
    Error_space = RMSError * MError / np.abs(MError)

    # should be the same normalisation as all data, have to change
    cmap_norm = matplotlib.colors.Normalize(cmin, cmax)
    scale_map = matplotlib.cm.ScalarMappable(norm=cmap_norm, cmap=cmap)

    fig.axs[0][i].scatter(Kmax_range, Ku_range,
                          c=scale_map.to_rgba(Error_space),
                          s=5, marker='o', zorder=-10)
    fig.axs[0][i].set_xscale('log')
    fig.axs[0][i].set_yscale('log')
    fig.axs[0][i].set_title(r'$V_{half-trap} = $' + "%.2e" %
                            (chosen_Vhalf_value[i]))

fig.sharex([r"$K_\mathrm{max}$"] * 3)
fig.sharey([r"$K_u$"])

cax = fig.axs[0][2].inset_axes([1.03, 0, 0.03, 1])
scale_map.set_array(Error_fullspace)
fig.fig.colorbar(scale_map, orientation='vertical', ax=fig.axs, cax=cax)

fig.fig.text(0.105, 0.9, '(A)', fontsize=11)
fig.fig.text(0.37, 0.9, '(B)', fontsize=11)
fig.fig.text(0.64, 0.9, '(C)', fontsize=11)

saved_fig_dir = '../../figures/sensitivity_analysis/'
fig.savefig(saved_fig_dir + 'SA_Vhalf_2D.pdf')
