import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import os
import pandas as pd

import modelling

saved_fig_dir = '../../figures/sensitivity_analysis/'

# 3D plot in parameter space
# Plot for known drugs
param_lib = modelling.BindingParameters()
drug_list = param_lib.drug_compounds

SA_model = modelling.SensitivityAnalysis()
param_names = SA_model.param_names

starting_param_df = pd.DataFrame([1] * 5, index=param_names).T
ComparisonController = modelling.ModelComparison(starting_param_df)

# Read data for drugs
saved_data_dir = '../../simulation_data/'
filename = 'SA_alldrugs.csv'
df = pd.read_csv(saved_data_dir + filename,
                 header=[0, 1], index_col=[0],
                 skipinitialspace=True)

Vhalf_list = df['param_values']['Vhalf'].values
Kmax_list = df['param_values']['Kmax'].values
Ku_list = df['param_values']['Ku'].values
drug_list = df['drug']['drug'].values

RMSError_drug = df['RMSE']['RMSE'].values
MError_drug = df['ME']['ME'].values

Error_drug = np.array(RMSError_drug) * np.array(MError_drug) / \
    np.abs(np.array(MError_drug))

# Read data for space
saved_data_dir = '../../simulation_results/'
file_prefix = 'SA_APD'
result_files = [saved_data_dir + f for f in os.listdir(saved_data_dir)
                if f.startswith(file_prefix)]

error_range = 30

first_iter = True
for file in result_files:
    df = pd.read_csv(file,
                     header=[0, 1], index_col=[0],
                     skipinitialspace=True)
    chosen_df = df.loc[df['RMSE']['RMSE'] < error_range]

    if first_iter:
        combined_df = df
        combined_chosen_df = chosen_df
        first_iter = False
    else:
        combined_chosen_df = pd.concat([combined_chosen_df, chosen_df])
        combined_df = pd.concat([combined_df, df])

Vhalf_range = combined_df['param_values']['Vhalf'].values
Kmax_range = combined_df['param_values']['Kmax'].values
Ku_range = combined_df['param_values']['Ku'].values

RMSError = combined_df['RMSE']['RMSE'].values
MError = combined_df['ME']['ME'].values

nan_ind = [i for i in range(len(RMSError)) if np.isnan(RMSError[i]) or
           np.isnan(MError[i])]
Error_space = RMSError * MError / np.abs(MError)

cmin = min(min(Error_drug), min(Error_space))
cmax = max(max(Error_drug), max(Error_space))

Vhalf_range = [Vhalf_range[i] for i in range(len(Vhalf_range))
               if i not in nan_ind]
Kmax_range = [Kmax_range[i] for i in range(len(Kmax_range))
              if i not in nan_ind]
Ku_range = [Ku_range[i] for i in range(len(Ku_range))
            if i not in nan_ind]
Error_space = [Error_space[i] for i in range(len(Error_space))
               if i not in nan_ind]

Vhalf_chosen = combined_chosen_df['param_values']['Vhalf'].values
Kmax_chosen = combined_chosen_df['param_values']['Kmax'].values
Ku_chosen = combined_chosen_df['param_values']['Ku'].values

# Read data for curve
saved_data_dir = '../../simulation_results/SA_curve/'
file_prefix = 'SA_curve'
result_files2 = [saved_data_dir + f for f in os.listdir(saved_data_dir)
                 if f.startswith(file_prefix)]

first_iter = True
for file in result_files2:
    df = pd.read_csv(file,
                     header=[0, 1], index_col=[0],
                     skipinitialspace=True)
    chosen_df = df.loc[df['RMSE']['RMSE'] < error_range]

    if first_iter:
        curve_chosen_df = chosen_df
        first_iter = False
    else:
        curve_chosen_df = pd.concat([curve_chosen_df, chosen_df])

Vhalf_curve = curve_chosen_df['param_values']['Vhalf'].values
Kmax_curve = curve_chosen_df['param_values']['Kmax'].values
Ku_curve = curve_chosen_df['param_values']['Ku'].values


def log_tick_formatter(val, pos=None):
    return f"$10^{{{int(val)}}}$"


fig = plt.figure(figsize=(10, 5))

gs = fig.add_gridspec(1, 2, wspace=0.1)
axs = [fig.add_subplot(gs[0, j], projection='3d') for j in range(2)]

cmap = plt.get_cmap('RdYlBu_r')
cmap_norm = matplotlib.colors.Normalize(cmin, cmax)
scale_map = matplotlib.cm.ScalarMappable(norm=cmap_norm, cmap=cmap)

axs[0].scatter(Vhalf_range, np.log10(Kmax_range), np.log10(Ku_range),
               c=scale_map.to_rgba(Error_space),
               s=5, marker='o', zorder=-10, alpha=0.5)
axs[0].view_init(20, 40)

axs[1].scatter(Vhalf_chosen, np.log10(Kmax_chosen), np.log10(Ku_chosen),
               c='dimgrey', s=10, marker='o', zorder=-10, alpha=0.5)
axs[1].scatter(Vhalf_curve, np.log10(Kmax_curve), np.log10(Ku_curve),
               c='dimgrey', s=10, marker='o', zorder=-10, alpha=0.5)
axs[1].scatter(Vhalf_list, np.log10(Kmax_list), np.log10(Ku_list),
               c=scale_map.to_rgba(Error_drug),
               s=100, marker='^', zorder=-1)
axs[1].view_init(10, 40)

for i in range(2):
    axs[i].set_xlabel(r"$V_\mathrm{half-trap}$")
    axs[i].set_ylabel(r"$K_\mathrm{max}$")
    axs[i].set_zlabel(r"$K_u$")

    axs[i].set_xlim(min(Vhalf_range), max(Vhalf_range))
    axs[i].set_ylim(min(np.log10(Kmax_range)), max(np.log10(Kmax_range)))
    axs[i].set_zlim(min(np.log10(Ku_range)), max(np.log10(Ku_range)))

    axs[i].yaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
    axs[i].yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    axs[i].zaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
    axs[i].zaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    axs[i].set_rasterization_zorder(0)

cax = axs[0].inset_axes([0.5, -0.08, 1, 0.03])
scale_map.set_array(Error_space)
fig.colorbar(scale_map, orientation='horizontal', ax=axs, cax=cax)

fig.text(0.075, 0.75, '(A)', fontsize=11)
fig.text(0.5, 0.75, '(B)', fontsize=11)

plt.subplots_adjust(hspace=0)

plt.savefig(saved_fig_dir + 'Fig_SA_3D_test.png', bbox_inches='tight')

############################################

fig = plt.figure(figsize=(10, 5))

gs = fig.add_gridspec(1, 2, wspace=0.1)
axs = [fig.add_subplot(gs[0, j], projection='3d') for j in range(2)]

cmap = plt.get_cmap('RdYlBu_r')
cmap_norm = matplotlib.colors.Normalize(cmin, cmax)
scale_map = matplotlib.cm.ScalarMappable(norm=cmap_norm, cmap=cmap)

for i in range(2):
    axs[i].scatter(Vhalf_chosen, np.log10(Kmax_chosen), np.log10(Ku_chosen),
                   c='dimgrey', s=10, marker='o', zorder=-10, alpha=0.5)
    axs[i].scatter(Vhalf_curve, np.log10(Kmax_curve), np.log10(Ku_curve),
                   c='dimgrey', s=10, marker='o', zorder=-10, alpha=0.5)
    axs[i].scatter(Vhalf_list, np.log10(Kmax_list), np.log10(Ku_list),
                   c=scale_map.to_rgba(Error_drug), s=100, marker='^')

axs[0].view_init(35, 60)
axs[1].view_init(10, 10)


for i in range(2):
    axs[i].set_xlabel(r"$V_\mathrm{half-trap}$")
    axs[i].set_ylabel(r"$K_\mathrm{max}$")
    axs[i].set_zlabel(r"$K_u$")

    axs[i].set_xlim(min(Vhalf_range), max(Vhalf_range))
    axs[i].set_ylim(min(np.log10(Kmax_range)), max(np.log10(Kmax_range)))
    axs[i].set_zlim(min(np.log10(Ku_range)), max(np.log10(Ku_range)))

    axs[i].yaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
    axs[i].yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    axs[i].zaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
    axs[i].zaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    axs[i].set_rasterization_zorder(0)

cax = axs[0].inset_axes([0.5, -0.08, 1, 0.03])
scale_map.set_array(Error_space)
fig.colorbar(scale_map, orientation='horizontal', ax=axs, cax=cax)

fig.text(0.075, 0.75, '(A)', fontsize=11)
fig.text(0.5, 0.75, '(B)', fontsize=11)

plt.subplots_adjust(hspace=0)

plt.savefig(saved_fig_dir + 'FigS_SA_3D.png', bbox_inches='tight')
