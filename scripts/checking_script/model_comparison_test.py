import matplotlib.pyplot as plt
import myokit
import numpy as np
import os
import pandas as pd
import pints

import modelling

drug = 'droperidol'
protocol_name = 'Milnes'

# Define directories to save simulated data
data_dir = '../../simulation_data/model_comparison/' + \
    drug + '/' + protocol_name + '/'
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)
result_filename = 'Hill_curve.txt'

# Model directory
current_model_filepath = '../../math_model/ohara-cipa-v1-2017-IKr-opt.mmt'
AP_model_filepath = '../../math_model/ohara-cipa-v1-2017-opt.mmt'

# Load current model and set Milnes' protocol
model, _, x = myokit.load(current_model_filepath)
current_model = modelling.BindingKinetics(model)

protocol_params = modelling.ProtocolParameters()
protocol = protocol_params.protocol_parameters[protocol_name]['function']
current_model.protocol = protocol

# Load AP model and set current protocol
APmodel, _, x = myokit.load(AP_model_filepath)
AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')
pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)

# Parameters used in simulations
offset = 50
save_signal = 2
repeats = 1000
APD_points = 20

# Get name of parameters
SA_model = modelling.SensitivityAnalysis()
param_names = SA_model.param_names

# Define parameter values of virtual drug
param_lib = modelling.BindingParameters()
param_values = param_lib.binding_parameters[drug]
param_values = pd.DataFrame(param_values, index=[0])
ComparisonController = modelling.ModelComparison(param_values)

# Compute Hill curve of the virtual drug with the SD model
Hill_curve_coefs, drug_conc_Hill, peaks_norm = \
    ComparisonController.compute_Hill(current_model)
# parameters of Hill curve are based on normalised drug concentration

with open(data_dir + result_filename, 'w') as f:
    for x in Hill_curve_coefs:
        f.write(pints.strfloat(x) + '\n')

plt.figure()
plt.plot(drug_conc_Hill, peaks_norm)
plt.xscale('log')
fig_dir = '../../figures/model_comparison/' + drug + '/' + \
    protocol_name + '/'
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)
plt.savefig(fig_dir + 'Hill_curve.pdf')

print(drug_conc_Hill)

# Define drug concentration range similar to the drug concentration used
# to infer Hill curve
drug_conc_AP = 10**np.linspace(np.log10(drug_conc_Hill[1]),
                               np.log10(max(drug_conc_Hill)),
                               APD_points)

# Simulate APs and APD90s of the ORd-SD model and the ORd-CS model
APD_trapping, APD_conductance, drug_conc_AP = \
    ComparisonController.APD_sim(
        AP_model, Hill_curve_coefs, drug_conc=drug_conc_AP,
        EAD=True)

# Calculate RMSD and MD of simulated APD90 of the two models
RMSError = ComparisonController.compute_RMSE(APD_trapping,
                                             APD_conductance)
MAError = ComparisonController.compute_ME(APD_trapping,
                                          APD_conductance)

plt.figure()
plt.plot(drug_conc_AP, APD_trapping)
plt.plot(drug_conc_AP, APD_conductance)
plt.xscale('log')
plt.savefig(fig_dir + 'APDs.pdf')

print(drug_conc_AP)
# # Create dataframe to save results
# conc_Hill_ind = ['conc_' + str(i) for i, _ in
#                     enumerate(drug_conc_Hill)]
# conc_AP_ind = ['conc_' + str(i) for i, _ in enumerate(drug_conc_AP)]
# index_dict = {'param_id': ['param_id'],
#                 'drug_conc_Hill': conc_Hill_ind,
#                 'peak_current': conc_Hill_ind,
#                 'Hill_curve': ['Hill_coef', 'IC50'],
#                 'param_values': param_names, 'drug_conc_AP': conc_AP_ind,
#                 'APD_trapping': conc_AP_ind,
#                 'APD_conductance': conc_AP_ind, 'RMSE': ['RMSE'],
#                 'ME': ['ME']}
# all_index = [(i, j) for i in index_dict.keys() for j in index_dict[i]]
# index = pd.MultiIndex.from_tuples(all_index)

# big_df = pd.DataFrame(
#     [param_id] + drug_conc_Hill + list(peaks_norm) +
#     list(Hill_curve_coefs) + list(param_values.values[0]) +
#     list(drug_conc_AP) + APD_trapping + APD_conductance +
#     [RMSError] + [MAError], index=index)
