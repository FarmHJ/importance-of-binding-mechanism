# To plot the range of parameter values
# Drug binding-related parameters

import matplotlib.pyplot as plt
import pandas as pd

import modelling


testing_fig_dir = '../../figures/testing/'
final_fig_dir = '../../figures/binding_kinetics_comparison/' + \
    'OHaraCiPA_model/sensitivity_analysis/'

check_plot = False
saved_fig_dir = final_fig_dir

saved_data_dir = '../../simulation_data/sensitivity_analysis/'

df = pd.read_csv(saved_data_dir + 'SA_cisapride_Vhalf.csv',
                 header=[0, 1], index_col=[0],
                 skipinitialspace=True)
# data included: drug_conc_Hill, peak_current, Hill_curve, param_values,
# drug_conc_AP, APD_trapping, APD_conductance and MSE

param_lib = modelling.BindingParameters()
drug = 'cisapride'
param_true = param_lib.binding_parameters[drug]['Vhalf']

# Plot Hill curve
if check_plot:
    row_ind = 0
    drug_conc_Hill = df.iloc[[row_ind]]['drug_conc_Hill'].values[0]
    peak_current = df.iloc[[row_ind]]['peak_current'].values[0]
    Hill_curve = df.iloc[[row_ind]]['Hill_curve'].values[0]

runs = len(df.index)
interest_param_values = df['param_values']['Vhalf'].values
MSEs = df['MSE']['MSE'].values

plt.figure()
plt.plot(interest_param_values, MSEs, 'o-')
plt.vlines(param_true, min(MSEs), max(MSEs), colors='red')
# plt.xscale('log')
plt.savefig(saved_fig_dir + 'cisapride_Vhalf_MSE.pdf')
