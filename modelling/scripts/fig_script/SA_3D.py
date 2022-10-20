import matplotlib
import matplotlib.pyplot as plt
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
MAError_drug = df['MAE']['MAE'].values

# Read data for space
saved_data_dir = '../../simulation_data/sensitivity_analysis/'
file_prefix = 'SA_allparam_'
result_files = [saved_data_dir + f for f in os.listdir(saved_data_dir) if f.startswith(file_prefix)]

saved_data_dir = '../../simulation_results/'
file_prefix = 'SA_allparam_gaps_'
result_files2 = [saved_data_dir + f for f in os.listdir(saved_data_dir) if f.startswith(file_prefix)]

result_files.extend(result_files2)

Vhalf_range = np.array([])
Kmax_range = np.array([])
Ku_range = np.array([])

RMSError = np.array([])
MAError = np.array([])

param_id = np.array([])
RMSE_range = 10

for file in result_files:
    df = pd.read_csv(file,
                     header=[0, 1], index_col=[0],
                     skipinitialspace=True)

    df = df.loc[df['RMSE']['RMSE'] < RMSE_range]
    print(df)
    # Vhalf_range = np.concatenate((Vhalf_range, df['param_values']['Vhalf'].values))
    # Kmax_range = np.concatenate((Kmax_range, df['param_values']['Kmax'].values))
    # Ku_range = np.concatenate((Ku_range, df['param_values']['Ku'].values))

    # RMSError = np.concatenate((RMSError, df['RMSE']['RMSE'].values))
    # MAError = np.concatenate((MAError, df['MAE']['MAE'].values))

    # param_id = np.concatenate((param_id, df['param_id']['param_id'].values))

RMSError_drug = np.array(RMSError_drug) * np.array(MAError_drug) / np.abs(np.array(MAError_drug))
RMSError_space = RMSError * MAError / np.abs(MAError)

cmin = min(min(RMSError_drug), min(RMSError_space))
cmax = max(max(RMSError_drug), max(RMSError_space))

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

cmap = plt.get_cmap('jet')
cmap_norm = matplotlib.colors.Normalize(cmin, cmax)
scale_map = matplotlib.cm.ScalarMappable(norm=cmap_norm, cmap=cmap)

# ax.scatter(Vhalf_list, np.log(Kmax_list), np.log(Ku_list),
#            c=scale_map.to_rgba(RMSError_drug),
#            s=100, marker='^', zorder=-10)
# ax.scatter(Vhalf_range, np.log(Kmax_range), np.log(Ku_range),
#            c=scale_map.to_rgba(RMSError_space),
#            s=10, marker='o', zorder=-10, alpha=0.5)

chosen_ind = [i for i, e in enumerate(RMSError_space) if e < 20 and e > -20]
Vhalf_chosen = np.array([Vhalf_range[i] for i in chosen_ind])
Kmax_chosen = np.array([Kmax_range[i] for i in chosen_ind])
Ku_chosen = np.array([Ku_range[i] for i in chosen_ind])
RMSE_chosen = np.array([RMSError_space[i] for i in chosen_ind])

# # Filling in gaps
# res = 5
# Vhalf_range = SA_model.param_explore('Vhalf', res)
# Kmax_range = SA_model.param_explore('Kmax', res)
# Ku_range = SA_model.param_explore('Ku', res)

# Vhalf_fullrange = SA_model.param_explore_gaps(Vhalf_range, 3, 'Vhalf')
# Kmax_fullrange = SA_model.param_explore_gaps(Kmax_range, 3, 'Kmax')
# Ku_fullrange = SA_model.param_explore_gaps(Ku_range, 3, 'Ku')

# # Do I have to make the value absolute
# Vhalf_min_diff = min(np.array(sorted(Vhalf_fullrange)[1:]) -
#                      np.array(sorted(Vhalf_fullrange)[:-1]))
# Kmax_min_diff = min(np.array(sorted(Kmax_fullrange)[1:]) -
#                     np.array(sorted(Kmax_fullrange)[:-1]))
# Ku_min_diff = min(np.array(sorted(Ku_fullrange)[1:]) -
#                   np.array(sorted(Ku_fullrange)[:-1]))

# X, Y, Z = np.meshgrid(Vhalf_fullrange, np.log(Kmax_fullrange),
#                       np.log(Ku_fullrange))
Kmax_chosen = np.log(Kmax_chosen)
Ku_chosen = np.log(Ku_chosen)

print(sorted(Ku_chosen))

half_len = int(np.ceil(len(Vhalf_chosen) / 2))
Vhalf_2D = np.array([Vhalf_chosen[:half_len], Vhalf_chosen[half_len-1:]])
Kmax_2D = np.array([Kmax_chosen[:half_len], Kmax_chosen[half_len-1:]])
Ku_2D = np.array([Ku_chosen[:half_len], Ku_chosen[half_len-1:]])
print(len(Vhalf_chosen[:half_len]))
print(len(Vhalf_chosen[half_len:]))
# Plot contour surfaces
_ = ax.scatter(
    Vhalf_2D, Ku_2D, Kmax_2D)
# _ = ax.contourf(
#     X[0, :, :], data[0, :, :], Z[0, :, :],
#     zdir='y', offset=0, **kw
# )
# C = ax.contourf(
#     data[:, -1, :], Y[:, -1, :], Z[:, -1, :],
#     zdir='x', offset=X.max(), **kw
# )
# scale_map.set_array(RMSError_drug)
# fig.colorbar(scale_map)

# handles, labels = ax.get_legend_handles_labels()
# unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if
#           l not in labels[:i]]
# ax.legend(*zip(*unique), loc='upper left', bbox_to_anchor=(1.0, 1.0))
# ax.set_facecolor('silver')
ax.set_xlabel('Vhalf')
ax.set_ylabel('Kmax')
ax.set_zlabel('Ku')
ax.set_rasterization_zorder(0)

saved_fig_dir = '../../figures/testing/'
plt.savefig(saved_fig_dir + 'test.pdf')
