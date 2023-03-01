import numpy as np
import os
import pandas as pd
import plotly.graph_objects as go

import modelling

category = False

# 3D plot in parameter space
# Plot for known drugs
param_lib = modelling.BindingParameters()
drug_list = param_lib.drug_compounds

SA_model = modelling.SensitivityAnalysis()
param_names = SA_model.param_names

starting_param_df = pd.DataFrame([1] * 5, index=param_names).T
ComparisonController = modelling.ModelComparison(starting_param_df)

discrete_colors = ['red', 'blue', 'black']
APD_diff_label = ['similar', 'SD higher', 'CS higher']

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

# Read data for space
saved_data_dir = '../../../simulation_data/parameter_space_exploration/SA_space/'
file_prefix = 'SA_allparam_uniform_opt_'
result_files = [saved_data_dir + f for f in os.listdir(saved_data_dir) if f.startswith(file_prefix)]

# saved_data_dir = '../../../simulation_results/'
# file_prefix = 'SA_allparam_gaps_'
# result_files2 = [saved_data_dir + f for f in os.listdir(saved_data_dir) if f.startswith(file_prefix)]

# result_files.extend(result_files2)

# # Read data for space
# saved_data_dir = '../../../simulation_results/'
# file_prefix = 'SA_APD'
# result_files = [saved_data_dir + f for f in os.listdir(saved_data_dir) if f.startswith(file_prefix)]

# Read data for curve
saved_data_dir = '../../../simulation_data/parameter_space_exploration/SA_curve/'
file_prefix = 'SA_curve_uniform_opt'
result_files.extend([saved_data_dir + f for f in os.listdir(saved_data_dir) if f.startswith(file_prefix)])

fig = go.Figure()

Vhalf_range = np.array([])
Kmax_range = np.array([])
Ku_range = np.array([])

RMSError = np.array([])
MAError = np.array([])

param_id = np.array([])

for file in result_files:
    df = pd.read_csv(file,
                     header=[0, 1], index_col=[0],
                     skipinitialspace=True)

    Vhalf_range = np.concatenate((Vhalf_range, df['param_values']['Vhalf'].values))
    Kmax_range = np.concatenate((Kmax_range, df['param_values']['Kmax'].values))
    Ku_range = np.concatenate((Ku_range, df['param_values']['Ku'].values))

    RMSError = np.concatenate((RMSError, df['RMSE']['RMSE'].values))
    MAError = np.concatenate((MAError, df['ME']['ME'].values))

    param_id = np.concatenate((param_id, df['param_id']['param_id'].values))

if category:
    color_code_list = []
    for i in range(len(drug_list)):
        if RMSError_drug[i] < 100:
            color_code = 0
        elif MAError_drug[i] > 0:
            color_code = 1
        else:
            color_code = 2

        color_code_list.append(color_code)

    fig.add_trace(
        go.Scatter3d(
            x=Vhalf_list,
            y=Kmax_list,
            z=Ku_list,
            mode='markers',
            marker_symbol='diamond',
            name='drugs',
            marker=dict(
                color=color_code_list,
                colorscale=discrete_colors
            )
        )
    )

    color_code_list = []
    for i in range(len(param_id)):
        if RMSError[i] < 100:
            color_code = 0
        elif MAError[i] > 0:
            color_code = 1  # SD higher
        else:
            color_code = 2  # CS higher
        color_code_list.append(color_code)

    fig.add_trace(
        go.Scatter3d(
            x=Vhalf_range,
            y=Kmax_range,
            z=Ku_range,
            mode='markers',
            name='space',
            marker=dict(
                color=color_code_list,
                colorscale=discrete_colors,
                opacity=0.5,
                size=5
            )
        )
    )
else:
    RMSError_drug = np.array(RMSError_drug)  #  * np.array(MAError_drug) / np.abs(np.array(MAError_drug))
    RMSError_space = RMSError  #  * MAError / np.abs(MAError)
    
    cmin = min(min(RMSError_drug), min(RMSError_space))
    cmax = max(max(RMSError_drug), max(RMSError_space))
    
    hovertext = np.empty(shape=(12,3,1), dtype='object')
    hovertext[:,0] = np.array(drug_list).reshape(-1,1)
    hovertext[:,1] = np.array(RMSError_drug).reshape(-1,1)
    hovertext[:,2] = np.array(MAError_drug).reshape(-1,1)
    fig.add_trace(
        go.Scatter3d(
            x=Vhalf_list,
            y=Kmax_list,
            z=Ku_list,
            mode='markers',
            marker_symbol='diamond',
            name='',
            customdata=hovertext,
            hovertemplate='<b>%{customdata[0]}</b> <br>RMSD = %{customdata[1]:.2e}',  # <br>MD = %{customdata[2]:.2e}',
            marker=dict(
                color=RMSError_drug,
                colorscale='Portland',
                colorbar=dict(thickness=20),
                cmin=cmin,
                cmax=cmax
            )
        )
    )
    
#     bin_arr = np.arange(-160, 130, 20)
    bin_arr = np.arange(0, 160, 30)
    bins = [(bin_arr[i], bin_arr[i + 1]) for i in range(len(bin_arr) - 1)]
    
    for lb, ub in bins:
        chosen_ind = [i for i, e in enumerate(RMSError_space) if e < ub and e > lb]
        Vhalf_chosen = np.array([Vhalf_range[i] for i in chosen_ind])
        Kmax_chosen = np.array([Kmax_range[i] for i in chosen_ind])
        Ku_chosen = np.array([Ku_range[i] for i in chosen_ind])
        RMSE_chosen = np.array([RMSError_space[i] for i in chosen_ind])
        paramid_chosen = np.array([param_id[i] for i in chosen_ind])
        MAE_chosen = np.array([MAError[i] for i in chosen_ind])
        
        Vhalf_bg = np.array([Vhalf_range[i] for i in range(len(RMSError_space)) if i not in chosen_ind])
        Kmax_bg = np.array([Kmax_range[i] for i in range(len(RMSError_space)) if i not in chosen_ind])
        Ku_bg = np.array([Ku_range[i] for i in range(len(RMSError_space)) if i not in chosen_ind])
        paramid_bg = np.array([param_id[i] for i in range(len(RMSError_space)) if i not in chosen_ind])

        hovertext = np.empty(shape=(len(paramid_chosen),3,1), dtype='object')
        hovertext[:,0] = np.array(paramid_chosen).reshape(-1,1)
        hovertext[:,1] = np.array(RMSE_chosen).reshape(-1,1)
        hovertext[:,2] = np.array(MAE_chosen).reshape(-1,1)

        fig.add_trace(
            go.Scatter3d(
                visible=True,
                x=Vhalf_chosen,
                y=Kmax_chosen,
                z=Ku_chosen,
                mode='markers',
                name='',
                customdata=hovertext,
                hovertemplate='<b>id: %{customdata[0]}</b> <br>RMSD = %{customdata[1]} <br>MD = %{customdata[2]}',
                marker=dict(
                    color=RMSE_chosen,
                    colorscale='Portland',
                    opacity=0.7,
                    size=3,
                    colorbar=dict(thickness=20),
                    cmin=cmin,
                    cmax=cmax
                )
            )
        )
#         fig.add_trace(
#             go.Scatter3d(
#                 visible=False,
#                 x=Vhalf_bg,
#                 y=Kmax_bg,
#                 z=Ku_bg,
#                 mode='markers',
#                 name='',
#                 customdata=hovertext,
#                 hovertemplate='<b>id: %{customdata[0]}</b> <br>RMSD = %{customdata[1]} <br>MD = %{customdata[2]}',
#                 marker=dict(
#                     color='#cccccc',
#                     opacity=0.4,
#                     size=3,
#                 )
#             )
#         )

sets = []
for i in range(int((len(fig.data) - 1))):
    param_set = dict(
        method="update",
        args=[{"visible": [False] * len(fig.data)},
              {"title": "APD90 difference at range " + "%d" % bins[i][0] + " to " +
               "%d" % bins[i][1]}],
        label= "(" + "%d" % bins[i][0] + ", " + "%d" % bins[i][1] + ")"
    )
    param_set["args"][0]["visible"][0] = True
    param_set["args"][0]["visible"][i + 1] = True
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

# Read data for space
saved_data_dir = '../../../simulation_results/'
file_prefix = 'SA_APD'
result_files = [saved_data_dir + f for f in os.listdir(saved_data_dir) if f.startswith(file_prefix)]

# Read data for curve
saved_data_dir = '../../../simulation_results/SA_curve/'
file_prefix = 'SA_curve'
result_files2 = [saved_data_dir + f for f in os.listdir(saved_data_dir) if f.startswith(file_prefix)]

# # Read data for curve
# saved_data_dir = '../../../simulation_data/sensitivity_analysis/'
# file_prefix = 'parameter_space_curve3'
# result_files2 = [saved_data_dir + f for f in os.listdir(saved_data_dir) if f.startswith(file_prefix)]

fig = go.Figure()

Vhalf_range = np.array([])
Kmax_range = np.array([])
Ku_range = np.array([])

RMSError_space = np.array([])

param_id = np.array([])

for file in result_files:
    df = pd.read_csv(file,
                     header=[0, 1], index_col=[0],
                     skipinitialspace=True)

    Vhalf_range = np.concatenate((Vhalf_range, df['param_values']['Vhalf'].values))
    Kmax_range = np.concatenate((Kmax_range, df['param_values']['Kmax'].values))
    Ku_range = np.concatenate((Ku_range, df['param_values']['Ku'].values))

    RMSError_space = np.concatenate((RMSError_space, df['RMSE']['RMSE'].values))

    param_id = np.concatenate((param_id, df['param_id']['param_id'].values))

cmin = min(min(RMSError_drug), min(RMSError_space))
cmax = max(max(RMSError_drug), max(RMSError_space))

lb, ub = 0, 30
chosen_ind = [i for i, e in enumerate(RMSError_space) if e < ub and e > lb]
Vhalf_chosen = np.array([Vhalf_range[i] for i in chosen_ind])
Kmax_chosen = np.array([Kmax_range[i] for i in chosen_ind])
Ku_chosen = np.array([Ku_range[i] for i in chosen_ind])
RMSE_chosen = np.array([RMSError_space[i] for i in chosen_ind])
paramid_chosen = np.array([param_id[i] for i in chosen_ind])

hovertext = np.empty(shape=(len(paramid_chosen),3,1), dtype='object')
hovertext[:,0] = np.array(paramid_chosen).reshape(-1,1)
hovertext[:,1] = np.array(RMSE_chosen).reshape(-1,1)

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
#             color='blue',
            color=RMSE_chosen,
            colorscale='Portland',
            opacity=0.7,
            size=3,
            colorbar=dict(thickness=20),
            cmin=cmin,
            cmax=cmax
        )
    )
)

Vhalf_curve = np.array([])
Kmax_curve = np.array([])
Ku_curve = np.array([])

RMSError_curve = np.array([])

paramid_curve = np.array([])

for file in result_files2:
    df = pd.read_csv(file,
                     header=[0, 1], index_col=[0],
                     skipinitialspace=True)

    Vhalf_curve = np.concatenate((Vhalf_curve, df['param_values']['Vhalf'].values))
    Kmax_curve = np.concatenate((Kmax_curve, df['param_values']['Kmax'].values))
    Ku_curve = np.concatenate((Ku_curve, df['param_values']['Ku'].values))

    RMSError_curve = np.concatenate((RMSError_curve, df['RMSE']['RMSE'].values))

    paramid_curve = np.concatenate((paramid_curve, df['param_id']['param_id'].values))

# for file in result_files2:
#     df = pd.read_csv(file,
#                      header=[0], index_col=[0],
#                      skipinitialspace=True)

#     Vhalf_curve = np.concatenate((Vhalf_curve, df['Vhalf'].values))
#     Kmax_curve = np.concatenate((Kmax_curve, df['Kmax'].values))
#     Ku_curve = np.concatenate((Ku_curve, df['Ku'].values))

#     paramid_curve = np.concatenate((paramid_curve, df['param_id'].values))
    
lb, ub = 0, 30
curve_chosen_ind = [i for i, e in enumerate(RMSError_curve) if e < ub and e > lb]
# curve_chosen_ind = np.arange(len(paramid_curve))
Vhalf_curve_chosen = np.array([Vhalf_curve[i] for i in curve_chosen_ind])
Kmax_curve_chosen = np.array([Kmax_curve[i] for i in curve_chosen_ind])
Ku_curve_chosen = np.array([Ku_curve[i] for i in curve_chosen_ind])
RMSE_curve_chosen = np.array([RMSError_curve[i] for i in curve_chosen_ind])
paramid_curve_chosen = np.array([paramid_curve[i] for i in curve_chosen_ind])

hovertext = np.empty(shape=(len(paramid_curve_chosen),3,1), dtype='object')
hovertext[:,0] = np.array(paramid_curve_chosen).reshape(-1,1)
# hovertext[:,1] = np.array(RMSError_curve).reshape(-1,1)

fig.add_trace(
    go.Scatter3d(
        visible=True,
        x=Vhalf_curve_chosen,
        y=Kmax_curve_chosen,
        z=Ku_curve_chosen,
        mode='markers',
        name='',
        customdata=hovertext,
        hovertemplate='<b>id: %{customdata[0]}</b>',  # <br>RMSD = %{customdata[1]}',
        marker=dict(
#             color='red',
            color=RMSE_curve_chosen,
            colorscale='Portland',
            opacity=0.7,
            size=3,
            colorbar=dict(thickness=20),
            cmin=cmin,
            cmax=cmax
        )
    )
)

fig.update_layout(
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
                                 range=[np.log10(min(Ku_curve)), np.log10(max(Ku_range))])),
                  scene_aspectmode='manual',
                  scene_aspectratio=dict(x=1, y=1.2, z=1))
#                   margin=dict(r=20, l=10, b=10, t=10))

fig.show()