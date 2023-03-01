import numpy as np
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import modelling

SA_model = modelling.SensitivityAnalysis()

colors = px.colors.sequential.Plasma
print(colors)

# # Read data for space
# saved_data_dir = '../../../simulation_results/SA_space/'
# file_prefix = 'SA_allparam'
# result_files = [saved_data_dir + f for f in os.listdir(saved_data_dir) if f.startswith(file_prefix)]
# # fields = ['param_values', 'RMSE', 'drug_conc_AP', 'APD_trapping', 'APD_conductance']

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
        
# combined_df = combined_df.sort_values(by=[('param_values', 'Ku'), ('param_values', 'Kmax'), ('param_values', 'Vhalf')],
#                                       ascending=[False, True, True])

# Vhalf_range = combined_df['param_values']['Vhalf'].values
# Kmax_range = combined_df['param_values']['Kmax'].values
# Ku_range = combined_df['param_values']['Ku'].values

# RMSError = combined_df['RMSE']['RMSE'].values
# MError = combined_df['ME']['ME'].values

# nan_ind = [i for i in range(len(RMSError)) if np.isnan(RMSError[i]) or np.isnan(MError[i])]
# Error_space = RMSError * MError / np.abs(MError)

# # cmin = min(min(Error_drug), min(Error_space))
# # cmax = max(max(Error_drug), max(Error_space))

# Read data for space
saved_data_dir = '../../../simulation_results/'
file_prefix = 'SA_APD'
result_files = [saved_data_dir + f for f in os.listdir(saved_data_dir) if f.startswith(file_prefix)]

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

RMSError = combined_df['RMSE']['RMSE'].values
MError = combined_df['ME']['ME'].values

nan_ind = [i for i in range(len(RMSError)) if np.isnan(RMSError[i]) or np.isnan(MError[i])]
Error_space = RMSError * MError / np.abs(MError)

# Kmax_fullrange = SA_model.param_explore('Kmax', 5)
# Kmax_fullrange = SA_model.param_explore_gaps(Kmax_fullrange, 3, 'Kmax')
# Kmax_length = len(Kmax_fullrange)

# Ku_fullrange = SA_model.param_explore('Ku', 5)
# Ku_fullrange = SA_model.param_explore_gaps(Ku_fullrange, 3, 'Ku')
# Ku_length = len(Ku_fullrange)

fig = make_subplots(rows=2, cols=1, row_heights=[0.7, 0.3])

previous_i = 0
count = 0
plot_ids = []
Kmax_values = [Kmax_range[0]]
Ku_values = [Ku_range[0]]
for i in range(1, int(len(Kmax_range) / 10)):
    if Kmax_range[i] != Kmax_range[i - 1]:
        
        param_range = Vhalf_range[previous_i:i]
        RMSError_plot = RMSError[previous_i:i]
        fig.add_trace(
            go.Scatter(
                visible=True,
                x=param_range,
                y=RMSError_plot,
                mode='lines+markers',
                name='root mean square difference'
            ), row=2, col=1
        )
        plot_id1 = count * 1
        count += 1

        for r in range(previous_i, i):
            APD_trapping = combined_df.iloc[[r]]['APD_trapping'].values[0]
            APD_conductance = combined_df.iloc[[r]]['APD_conductance'].values[0]

            fig.add_trace(
                go.Scatter(
                    visible=True,
                    x=np.arange(len(APD_trapping)),
                    y=APD_trapping,
                    mode='lines+markers',
                    name='APD_trapping',
                    line=dict(  # color=colors[min(r, 9)],)
                        color="#ffa500",)
#                              colorscale='OrRd')
                ), row=1, col=1,
            )
            fig.add_trace(
                go.Scatter(
                    visible=True,
                    x=np.arange(len(APD_conductance)),
                    y=APD_conductance,
                    mode='lines+markers',
                    name='APD_conductance',
                    marker=dict(  # color=r,
                        color="#0000ff",)
#                              colorscale='GnBu')
                ), row=1, col=1,
            )
            count += 2
        plot_id2 = count * 1
        
        plot_ids.append((plot_id1, plot_id2))
        Kmax_values.append(Kmax_range[i])
        Ku_values.append(Ku_range[i])
        previous_i = i
        
sets = []
for i in range(len(plot_ids)):
    param_set = dict(
        method="update",
        args=[{"visible": [False] * len(fig.data)},],
#               {"title": "Drug " + drug + " at parameter " + param_interest + " = " +
#                "%.3f" % param_range[i]}],
        label="(%.2e, %.2e)" % (Ku_values[i], Kmax_values[i])
    )

    param_set["args"][0]["visible"][plot_ids[i][0]] = True
    length = plot_ids[i][1] - plot_ids[i][0]
    param_set["args"][0]["visible"][plot_ids[i][0]:plot_ids[i][1]] = [True] * length
    sets.append(param_set)

sliders = [dict(
    active=5,
    currentvalue={"prefix": "(Ku, Kmax) = "},
    pad={"t": 10},
    steps=sets
)]

fig.update_layout(sliders=sliders)  # , yaxis1=dict(range=[0, 1050]),
#                   xaxis1=dict(range=[np.log10(min_drug_conc), np.log10(max_drug_conc)]))

fig.update_xaxes(title_text="Normalised drug concentration", row=1, col=1)
fig.update_xaxes(title_text="Parameter value", row=2, col=1)

fig.update_yaxes(title_text="APD90", row=1, col=1)
fig.update_yaxes(title_text="APD difference", row=2, col=1)
fig.show()