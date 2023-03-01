"""
To plot the APD90 difference in the domain of the 3 interested
parameters.
"""

import numpy as np
import os
import pandas as pd
import plotly.graph_objects as go

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
drug_df = pd.read_csv(saved_data_dir + filename,
                 header=[0, 1], index_col=[0],
                 skipinitialspace=True)

Vhalf_list = drug_df['param_values']['Vhalf'].values
Kmax_list = drug_df['param_values']['Kmax'].values
Ku_list = drug_df['param_values']['Ku'].values
drug_list = drug_df['drug']['drug'].values

# Read data for space
saved_data_dir = '../../../simulation_results/'
file_prefix = 'SA_APD'
result_files = [saved_data_dir + f for f in os.listdir(saved_data_dir) if f.startswith(file_prefix)]

error_choice = 'RMSE'
error_choice2 = 'ME'

Error_drug = drug_df[error_choice][error_choice].values
Error_drug2 = drug_df[error_choice2][error_choice2].values

Vhalf_range = np.array([])
Kmax_range = np.array([])
Ku_range = np.array([])

RMSError = np.array([])
Error_space = np.array([])
Error_space2 = np.array([])

param_id = np.array([])

for file in result_files:
    df = pd.read_csv(file,
                     header=[0, 1], index_col=[0],
                     skipinitialspace=True)

    Vhalf_range = np.concatenate((Vhalf_range, df['param_values']['Vhalf'].values))
    Kmax_range = np.concatenate((Kmax_range, df['param_values']['Kmax'].values))
    Ku_range = np.concatenate((Ku_range, df['param_values']['Ku'].values))

    Error_space = np.concatenate((Error_space, df[error_choice][error_choice].values))
    Error_space2 = np.concatenate((Error_space2, df[error_choice2][error_choice2].values))

    param_id = np.concatenate((param_id, df['param_id']['param_id'].values))

Error_drug = np.array(Error_drug) * np.array(Error_drug2) / np.abs(np.array(Error_drug2))
Error_space = Error_space * Error_space2 / np.abs(Error_space2)

cmin = min(min(Error_drug), min(Error_space))
cmax = max(max(Error_drug), max(Error_space))

fig = go.Figure()

hovertext = np.empty(shape=(12,2,1), dtype='object')
hovertext[:,0] = np.array(drug_list).reshape(-1,1)
hovertext[:,1] = np.array(Error_drug).reshape(-1,1)
# hovertext[:,2] = np.array(Error_drug).reshape(-1,1)
fig.add_trace(
    go.Scatter3d(
        x=Vhalf_list,
        y=Kmax_list,
        z=Ku_list,
        mode='markers',
        marker_symbol='diamond',
        name='',
        customdata=hovertext,
        hovertemplate='<b>%{customdata[0]}</b> <br>Error = %{customdata[1]:.2e}',
        marker=dict(
            color=Error_drug,
            colorscale='Portland',
            colorbar=dict(thickness=20),
            cmin=cmin,
            cmax=cmax
        )
    )
)

bin_arr = np.arange(np.floor(min(Error_space)), max(Error_space) + 4, 5)
bins = [(bin_arr[i], bin_arr[i + 1]) for i in range(len(bin_arr) - 1)]

for lb, ub in bins:
    chosen_ind = [i for i, e in enumerate(Error_space) if e < ub and e > lb]
    Vhalf_chosen = np.array([Vhalf_range[i] for i in chosen_ind])
    Kmax_chosen = np.array([Kmax_range[i] for i in chosen_ind])
    Ku_chosen = np.array([Ku_range[i] for i in chosen_ind])
    paramid_chosen = np.array([param_id[i] for i in chosen_ind])
    Error_chosen = np.array([Error_space[i] for i in chosen_ind])

    Vhalf_bg = np.array([Vhalf_range[i] for i in range(len(Error_space)) if i not in chosen_ind])
    Kmax_bg = np.array([Kmax_range[i] for i in range(len(Error_space)) if i not in chosen_ind])
    Ku_bg = np.array([Ku_range[i] for i in range(len(Error_space)) if i not in chosen_ind])
    paramid_bg = np.array([param_id[i] for i in range(len(Error_space)) if i not in chosen_ind])

    hovertext = np.empty(shape=(len(paramid_chosen),2,1), dtype='object')
    hovertext[:,0] = np.array(paramid_chosen).reshape(-1,1)
    hovertext[:,1] = np.array(Error_chosen).reshape(-1,1)

    fig.add_trace(
        go.Scatter3d(
            visible=True,
            x=Vhalf_chosen,
            y=Kmax_chosen,
            z=Ku_chosen,
            mode='markers',
            name='',
            customdata=hovertext,
            hovertemplate='<b>id: %{customdata[0]}</b> <br>RMSD = %{customdata[1]}',
            marker=dict(
                color=Error_chosen,
                colorscale='Portland',
                opacity=0.7,
                size=3,
                colorbar=dict(thickness=20),
                cmin=cmin,
                cmax=cmax
            )
        )
    )
#     fig.add_trace(
#         go.Scatter3d(
#             visible=False,
#             x=Vhalf_bg,
#             y=Kmax_bg,
#             z=Ku_bg,
#             mode='markers',
#             name='',
# #             customdata=hovertext,
# #             hovertemplate='<b>id: %{customdata[0]}</b> <br>RMSD = %{customdata[1]} <br>MD = %{customdata[2]}',
#             marker=dict(
#                 color='lightgray',
#                 opacity=0.3,
#                 size=3,
#             )
#         )
#     )

sets = []
for i in range(int((len(fig.data) - 1))):
    param_set = dict(
        method="update",
        args=[{"visible": [False] * len(fig.data)},
              {"title": error_choice + " at range " + "%d" % bins[i][0] + " to " +
               "%d" % bins[i][1]}],
        label= "(" + "%d" % bins[i][0] + ", " + "%d" % bins[i][1] + ")"
    )
    param_set["args"][0]["visible"][0] = True
    param_set["args"][0]["visible"][:i + 1] = [True] * i
#     param_set["args"][0]["visible"][2 * i + 1] = True
#     param_set["args"][0]["visible"][2 * i + 2] = True
    sets.append(param_set)

sliders = [dict(
    active=5,
    currentvalue={"prefix": "Bins at "},
    pad={"t": 5},
    steps=sets
)]
        
fig.update_layout(sliders=sliders,
                  scene = dict(
                    xaxis_title='Vhalf',
                    yaxis_title='Kmax',
                    zaxis_title='Ku',
                    xaxis = dict(range=[min(Vhalf_range), max(Vhalf_range)]),
                    yaxis = dict(dtick=1,
                                 type='log',
                                 range=[np.log10(min(Kmax_range)), np.log10(max(Kmax_range))]),
                    zaxis = dict(dtick=1,
                                 type='log',
                                 range=[np.log10(min(Ku_range)), np.log10(max(Ku_range))])),
                  scene_aspectmode='manual',
                  scene_aspectratio=dict(x=1, y=1.2, z=1))
#                   margin=dict(r=20, l=10, b=10, t=10))

fig.show()

res = 5
Kmax_temp = SA_model.param_explore('Vhalf', res)

# Filling in gaps
Kmax_fullrange = SA_model.param_explore_gaps(Kmax_temp, 3, 'Vhalf')
print(len(Kmax_fullrange))
# print(np.log10(sorted(Kmax_fullrange)))