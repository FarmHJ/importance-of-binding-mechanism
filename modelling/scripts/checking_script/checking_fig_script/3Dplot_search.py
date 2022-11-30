import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import os
import pandas as pd

import modelling

# 3D plot in parameter space
# Plot for known drugs
param_lib = modelling.BindingParameters()
drug_list = param_lib.drug_compounds

SA_model = modelling.SensitivityAnalysis()
param_names = SA_model.param_names

starting_param_df = pd.DataFrame([1] * 5, index=param_names).T
ComparisonController = modelling.ModelComparison(starting_param_df)

# Read data for drugs
saved_data_dir = '../../../simulation_data/'
filename = 'SA_alldrugs.csv'
df = pd.read_csv(saved_data_dir + filename,
                 header=[0, 1], index_col=[0],
                 skipinitialspace=True)

Vhalf_list = df['param_values']['Vhalf'].values
Kmax_list = df['param_values']['Kmax'].values
Ku_list = df['param_values']['Ku'].values
drug_list = df['drug']['drug'].values

RMSError_drug = df['RMSE']['RMSE'].values
MAError_drug = df['ME']['ME'].values

Error_drug = np.array(RMSError_drug) * np.array(MAError_drug) / np.abs(np.array(MAError_drug))

# Read data for space
saved_data_dir = '../../../simulation_results/'
file_prefix = 'SA_APD'
result_files = [saved_data_dir + f for f in os.listdir(saved_data_dir) if f.startswith(file_prefix)]

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
        

# combined_df = combined_df.sort_values(by=[('param_values', 'Ku'), ('param_values', 'Kmax'), ('param_values', 'Vhalf')],
#                                       ascending=[False, True, True])

Vhalf_range = combined_df['param_values']['Vhalf'].values
Kmax_range = combined_df['param_values']['Kmax'].values
Ku_range = combined_df['param_values']['Ku'].values

RMSError = combined_df['RMSE']['RMSE'].values
MError = combined_df['ME']['ME'].values

nan_ind = [i for i in range(len(RMSError)) if np.isnan(RMSError[i]) or np.isnan(MError[i])]
Error_space = RMSError * MError / np.abs(MError)

cmin = min(min(Error_drug), min(Error_space))
cmax = max(max(Error_drug), max(Error_space))

# Read data for space
saved_data_dir = '../../../simulation_results/SA_space/'
file_prefix = 'filling_nan'
result_files = [saved_data_dir + f for f in os.listdir(saved_data_dir) if f.startswith(file_prefix)]

first_iter = True
for file in result_files:
    nan_df = pd.read_csv(file,
                     header=[0, 1], index_col=[0],
                     skipinitialspace=True)
    
    if first_iter:
        combined_nan_df = nan_df
        first_iter = False
    else:
        combined_nan_df = pd.concat([combined_nan_df, nan_df])

Vhalf_nan = combined_nan_df['param_values']['Vhalf'].values
Kmax_nan = combined_nan_df['param_values']['Kmax'].values
Ku_nan = combined_nan_df['param_values']['Ku'].values

RMSError_nan = combined_nan_df['RMSE']['RMSE'].values
MError_nan = combined_nan_df['ME']['ME'].values

# nan_ind = [i for i in range(len(RMSError)) if np.isnan(RMSError[i]) or np.isnan(MError[i])]
Error_nan = RMSError_nan * MError_nan / np.abs(MError_nan)

Vhalf_range = combined_df['param_values']['Vhalf'].values
Kmax_range = combined_df['param_values']['Kmax'].values
Ku_range = combined_df['param_values']['Ku'].values

Error_space = RMSError * MError / np.abs(MError)

Vhalf_range = [Vhalf_range[i] for i in range(len(Vhalf_range)) if i not in nan_ind]
Kmax_range = [Kmax_range[i] for i in range(len(Kmax_range)) if i not in nan_ind]
Ku_range = [Ku_range[i] for i in range(len(Ku_range)) if i not in nan_ind]
Error_space = [Error_space[i] for i in range(len(Error_space)) if i not in nan_ind]

Vhalf_chosen = combined_chosen_df['param_values']['Vhalf'].values
Kmax_chosen = combined_chosen_df['param_values']['Kmax'].values
Ku_chosen = combined_chosen_df['param_values']['Ku'].values

# Read data for curve
saved_data_dir = '../../../simulation_results/SA_curve/'
file_prefix = 'SA_curve'
result_files2 = [saved_data_dir + f for f in os.listdir(saved_data_dir) if f.startswith(file_prefix)]

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
           s=10, marker='o', zorder=-10, alpha=0.5)
axs[0].scatter(Vhalf_nan, np.log10(Kmax_nan), np.log10(Ku_nan),
           c=scale_map.to_rgba(Error_nan),
           s=10, marker='o', zorder=-10, alpha=0.5)
axs[0].view_init(20, 40)


axs[1].scatter(Vhalf_chosen, np.log10(Kmax_chosen), np.log10(Ku_chosen),
           c='dimgrey',
           s=10, marker='o', zorder=-10, alpha=0.5)
axs[1].scatter(Vhalf_curve, np.log10(Kmax_curve), np.log10(Ku_curve),
           c='dimgrey',
           s=10, marker='o', zorder=-10, alpha=0.5)
axs[1].scatter(Vhalf_list, np.log10(Kmax_list), np.log10(Ku_list),
           c=scale_map.to_rgba(Error_drug),
           s=100, marker='^', zorder=-1)
# axs[1].plot(X, np.log10(Y), np.log10(Z))
# axs[1].contour(Vhalf_chosen, np.log(Kmax_chosen), np.log(Ku_chosen),
#                zdir='x', offset=min(Vhalf_range), cmap=matplotlib.cm.coolwarm)
axs[1].view_init(30, 30)
# handles, labels = ax.get_legend_handles_labels()
# unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if
#           l not in labels[:i]]
# ax.legend(*zip(*unique), loc='upper left', bbox_to_anchor=(1.0, 1.0))
# ax.set_facecolor('silver')

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
fig.colorbar(scale_map, orientation='horizontal', ax=axs, cax=cax, ) # , shrink=0.5, )

fig.text(0.075, 0.75, '(A)', fontsize=11)
fig.text(0.5, 0.75, '(B)', fontsize=11)

# plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.subplots_adjust(hspace=0)

saved_fig_dir = '../../../figures/testing/'
plt.savefig(saved_fig_dir + 'test.png', bbox_inches='tight')

# Read data for space
saved_data_dir = '../../../simulation_results/parameter_space/'
file_prefix = 'parameter_space_curve'
result_files = [saved_data_dir + f for f in os.listdir(saved_data_dir) if f.startswith(file_prefix)]

# error_range = 20

first_iter = True
for file in result_files:
    df = pd.read_csv(file,
                     header=[0], index_col=[0],
                     skipinitialspace=True)
#     chosen_df = df.loc[df['RMSE']['RMSE'] < error_range]
    
    if first_iter:
        combined_df = df
#         combined_chosen_df = chosen_df
        first_iter = False
    else:
#         combined_chosen_df = pd.concat([combined_chosen_df, chosen_df])
        combined_df = pd.concat([combined_df, df])

Vhalf_curve = combined_df['Vhalf'].values
Kmax_curve = combined_df['Kmax'].values
Ku_curve = combined_df['Ku'].values

unique_pair = [i for i in range(len(Ku_curve)) if (Kmax_curve[i], Ku_curve[i])
               not in zip(Kmax_curve[:i], Ku_curve[:i])]
Y_curve = [Kmax_curve[i] for i in unique_pair]
Z_curve = [Ku_curve[i] for i in unique_pair]
X_curve = [min(Vhalf_curve)] * len(unique_pair)

unique_pair_end = unique_pair[1:] + [len(Ku_curve)]
color = np.array(unique_pair_end) - np.array(unique_pair)
# color[-1] = 25

sort_ind = [i[0] for i in sorted(enumerate(Y_curve), key=lambda x:x[1])]
Y_curve = sorted(Y_curve)
Z_curve = [Z_curve[i] for i in sort_ind]
# color = [color[i] for i in sort_ind]

cmap = plt.get_cmap('RdYlBu_r')
cmap_norm = matplotlib.colors.Normalize(1, 25)

plt.figure()
# for i in range(len(Y)):
plt.scatter(Y_curve, Z_curve)
plt.scatter(Y, Z)  # , color=cmap(cmap_norm(color[i])))
plt.xscale('log')
plt.yscale('log')
plt.show()
