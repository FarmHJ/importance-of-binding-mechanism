"""
The APD90s are saved but not the action potentials themselves.
So we have to choose the parameter sets where the RMS difference and
mean difference look weird, then simulate the action potentials again.
"""

import myokit
import numpy as np
import os
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import modelling

## Set up models and classes
saved_data_dir = '../../../simulation_data/sensitivity_analysis/'

model_dir = '../../../model/'
AP_model_filename = 'ohara-cipa-v1-2017.mmt'

# Load AP model
APmodel, _, x = myokit.load(model_dir + AP_model_filename)
AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')
pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)
base_conductance = APmodel.get('ikr.gKr').value()

steady_state_pulse = 1000
save_signal = 2
norm_constant = 1

Hill_model = modelling.HillsModel()
param_lib = modelling.BindingParameters()

## Set up the choices of parameter sets and their Hill's curve coefficients
choosing_method = 'param_id' # 'RMSError' or 'param_id'
if choosing_method == 'RMSError':
    choice_label = ['previous', 'next']
    param_interest = 'N'
    
    drug_list = ['dofetilide', 'verapamil', 'terfenadine',
             'cisapride', 'quinidine', 'sotalol']
    drug = drug_list[1]
    
    saved_data_dir = saved_data_dir + param_interest + '/'
    filename = 'SA_' + drug + '_' + param_interest + '.csv'
    
    slider_prefix = param_interest + " = "

elif choosing_method == 'param_id':
    
    id_interest = [2768, 2769]
    choice_label = ['%d' % i for i in id_interest]
    id_file = 2
    
    filename = 'SA_allparam_' + str(id_file) + '.csv'
    
    slider_prefix = 'Parameter id: '
df = pd.read_csv(saved_data_dir + filename,
                    header=[0, 1], index_col=[0],
                    skipinitialspace=True)
# data included: drug_conc_Hill, peak_current, Hill_curve, param_values,
# drug_conc_AP, APD_trapping, APD_conductance and MSE

if choosing_method == 'RMSError':
    
    df = df.sort_values(by=[('param_values', param_interest)])
    
    # Choosing interesting or weird output of RMS difference and mean difference
    RMSError = df['RMSE']['RMSE'].values
    RMSError_change = [RMSError[i + 1] - RMSError[i] for i in range(len(RMSError) - 1)]
    RMS_bigchange_ind = [(i, i + 1) for i in range(len(RMSError_change)) if RMSError_change[i] < 0]

    paramid_choice = RMS_bigchange_ind[0]
    
elif choosing_method == 'param_id':
    df = df.reset_index(drop=True)
    paramid_choice = []
    # Choosing parameter id
    for chosen_id in id_interest:
        paramid_choice.append(df[df[('param_id', 'param_id')] == chosen_id].index[0])

param_set_list = []
Hills_coef_list = []

drugconc_list = []
APDtrap_list = []
APDconduct_list = []

paramrange_list = []
RMSError_list = []
MError_list = []
for id_i in paramid_choice:
    param_set = df.iloc[[id_i]]['param_values']
    Hills_coef = df.iloc[[id_i]]['Hill_curve'].values[0]
    
    param_set_list.append(param_set)
    Hills_coef_list.append(Hills_coef)
    
    drug_conc_AP = df.iloc[[id_i]]['drug_conc_AP'].values[0]
    APD_trapping = df.iloc[[id_i]]['APD_trapping'].values[0]
    APD_conductance = df.iloc[[id_i]]['APD_conductance'].values[0]
    
    drugconc_list.append(drug_conc_AP)
    APDtrap_list.append(APD_trapping)
    APDconduct_list.append(APD_conductance)
    
    if choosing_method == 'RMSError':
        paramrange_list.append(df['param_values'][param_interest].values[id_i])
    RMSError_list.append(df['RMSE']['RMSE'].values[id_i])
    MError_list.append(df['MAE']['MAE'].values[id_i])
    
if choosing_method == 'param_id':
    paramrange_list = id_interest

## Run simulations for chosen parameter sets
logtrap_short = []
logconduct_short = []
drugconc_short = []
for point_n in range(len(param_set_list)):
    param_values = param_set_list[point_n]
    param_values['EC50'] = 1

    Hills_coef = Hills_coef_list[point_n]
    drug_conc_AP  = drugconc_list[point_n]
    APD_trapping = APDtrap_list[point_n]
    APD_conductance = APDconduct_list[point_n]
    
    EAD_trapping = [ind for ind, i in enumerate(APD_trapping) if i >= 1000]
    EAD_conductance = [ind for ind, i in enumerate(APD_conductance) if i >= 1000]
    EAD_conc_max = int(max(EAD_trapping[0], EAD_conductance[0]) + 1)
    EAD_conc_min = int(min(EAD_trapping[0], EAD_conductance[0]) - 1) 

    log_conc_trap = []
    log_conc_conduct = []
    for i in range(len(drug_conc_AP[EAD_conc_min:EAD_conc_max])):
        # Run simulation for trapping model
        log = AP_model.custom_simulation(
            param_values, drug_conc_AP[i + EAD_conc_min], steady_state_pulse,
            timestep=0.1, save_signal=save_signal,
#             abs_tol=1e-7, rel_tol=1e-8,
            log_var=['engine.time', 'membrane.V'])
        log_conc_trap.append(log)

        # Run simulation for conductance model
        reduction_scale = Hill_model.simulate(
            Hills_coef, drug_conc_AP[i + EAD_conc_min] * norm_constant)
        d2 = AP_model.conductance_simulation(
            base_conductance * reduction_scale, steady_state_pulse,
            timestep=0.1, save_signal=save_signal,
#             abs_tol=1e-7, rel_tol=1e-8,
            log_var=['engine.time', 'membrane.V'])
        log_conc_conduct.append(d2)

    drugconc_short.append(drug_conc_AP[EAD_conc_min:EAD_conc_max])
    logtrap_short.append(log_conc_trap)
    logconduct_short.append(log_conc_conduct)

## Plot figure
fig = make_subplots(rows=2, cols=2, row_heights=[0.6, 0.4],
                    subplot_titles=("SD model", "CS model", "APD difference", "APD"))

for i in range(len(paramid_choice)):
    drug_conc = drugconc_short[i]
    colorscale = ['hsl(0,100%,' + str(int(c)) + '%)' for c in np.linspace(80, 20, len(drug_conc))]
    for j in range(len(logtrap_short[i])):
        log = logtrap_short[i][j]
        log_time = np.concatenate((log.time(), log.time() + max(log.time())))
        log_voltage = np.concatenate((log['membrane.V', 0], log['membrane.V', 1]))

        fig.add_trace(
            go.Scatter(
                x=log_time,
                y=log_voltage,
                mode='lines',
                name=choice_label[i] + ' at ' + '{:.2e}'.format(drug_conc[j]),
                line=dict(color=colorscale[j])
            ), row=1, col=1)

        d2 = logconduct_short[i][j]
        log_time = np.concatenate((d2.time(), d2.time() + max(d2.time())))
        log_voltage = np.concatenate((d2['membrane.V', 0], d2['membrane.V', 1]))
        
        fig.add_trace(
            go.Scatter(
                x=log_time,
                y=log_voltage,
                mode='lines',
                name=choice_label[i] + ' at ' + '{:.2e}'.format(drug_conc[j]),
                line=dict(color=colorscale[j]),
                showlegend=False,
            ), row=1, col=2)

    fig.add_trace(
        go.Scatter(
            x=drugconc_list[i],
            y=APDtrap_list[i],
            mode='lines+markers',
            name='APD_trapping',
            line=dict(color='orange'),
        ), row=2, col=2,
    )
    fig.add_trace(
        go.Scatter(
            x=drugconc_list[i],
            y=APDconduct_list[i],
            mode='lines+markers',
            name='APD_conductance',
            line=dict(color='blue'),
        ), row=2, col=2,
        )

fig.add_trace(
    go.Scatter(
        visible=True,
        x=paramrange_list,
        y=RMSError_list,
        mode='lines+markers',
        name='root mean square difference',
        line=dict(color="#ff0000"),
    ), row=2, col=1
)
fig.add_trace(
    go.Scatter(
        visible=True,
        x=paramrange_list,
        y=MError_list,
        mode='lines+markers',
        name='mean difference',
        line=dict(color="brown"),
    ), row=2, col=1
)

sets = []
for i in range(len(paramid_choice)):
    if choosing_method == 'RMSError':
        fig_title = "Drug " + drug + " at parameter " + param_interest + " = " + \
                    "%.3f" % param_range[i]
    elif choosing_method == 'param_id':
        fig_title = "Parameter values (Vhalf, Kmax, Ku) = (" + "%.2e" % param_set_list[i]['Vhalf'] + \
                    ', ' + "%.2e" % param_set_list[i]['Kmax'] + ', ' + "%.2e" % param_set_list[i]['Ku'] + ')'
    param_set = dict(
        method="update",
        args=[{"visible": [False] * len(fig.data)},
              {"title": fig_title}
             ],
        label=choice_label[i]
    )
    param_set["args"][0]["visible"][-1] = True
    param_set["args"][0]["visible"][-2] = True
    param_set["args"][0]["visible"][
        len(logtrap_short[i]) * 2 * i:len(logtrap_short[i]) * 2 * (i + 1) + 2] = \
            [True] * (len(logtrap_short[i]) * 2 + 2)
    sets.append(param_set)



sliders = [dict(
    active=0,
    currentvalue={"prefix": slider_prefix},
    pad={"t": 50},
    steps=sets
)]

fig.update_layout(sliders=sliders)

fig.update_xaxes(title_text="Normalised drug concentration", type="log", row=2, col=2)
fig.update_xaxes(title_text="Parameter value", row=2, col=1)

fig.update_yaxes(title_text="APD90", row=2, col=2)
fig.update_yaxes(title_text="APD difference", row=2, col=1)
fig.show()