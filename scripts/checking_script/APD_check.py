# Plot the signed RMSD of the explored parameter space and the region where
# the RMSD between the APD90 of the ORd-SD model and the ORd-CS model is small
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import myokit
import numpy as np
import os
import pandas as pd

import modelling

# Define directory to save figure
fig_dir = '../../figures/parameter_exploration/'
fig_dir = '../../testing_figures/APD_check/'
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

# Get list of synthetic drugs' names
param_lib = modelling.BindingParameters()
drug_list = param_lib.drug_compounds

# Get the name of parameters
SA_model = modelling.SensitivityAnalysis()
param_names = SA_model.param_names

root_dir = '../../simulation_data/'

# Read simulated data of virtual drugs in the parameter space
# data_dir = root_dir + 'parameter_space_exploration/'
# file_prefix = 'SA_APD'
data_dir = root_dir + 'parameter_space_exploration/SA_space/'
file_prefix = 'SA_allparam_uniform_opt'
# data_dir = root_dir + 'parameter_space_exploration/SA_curve/'
# file_prefix = 'SA_curve_uniform_opt'
# file_prefix = 'SA_allparam'
result_files = [data_dir + f for f in os.listdir(data_dir)
                if f.startswith(file_prefix)]

# Define the range where the RMSD between the APD90s of the two models are
# small
error_range = 30

# Load results and extract points where the RMSD value is within the defined
# range
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

conc_AP_arr = combined_df['drug_conc_AP']
APD_SD_arr = combined_df['APD_trapping']
APD_CS_arr = combined_df['APD_conductance']

data_points, conc_points = np.shape(APD_SD_arr)

plt.figure()
count = 0
check_id = []
for i in range(data_points):
    conc = conc_AP_arr.iloc[i]
    APD = APD_SD_arr.iloc[i]
    conc = conc[~np.isnan(APD)]
    APD = APD[~np.isnan(APD)]
    res = all(i <= j for i, j in zip(APD, APD[1:]))
    diff = APD.values[1:] - APD.values[:-1]
    if not res:
        plt.plot(conc.values[:-1], diff, label=str(count))
        print(APD.values)
        check_id.append(i)
        count += 1
plt.legend()
plt.xscale('log')
plt.savefig(fig_dir + 'SD_space_APD.pdf')
print(count)

## Simulate AP
# Model directory
AP_model_filepath = '../../math_model/ohara-cipa-v1-2017-opt.mmt'

# Load AP model and set current protocol
APmodel, _, x = myokit.load(AP_model_filepath)
AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')
pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)

APD_points = 20
abs_tol = 1e-7
rel_tol = 1e-8
save_signal = 2

check_id = [check_id[i] for i in [3, 4]]
plt.figure()
for id in check_id:
    param_values = combined_df.iloc[[id]]['param_values']
    # Run simulation for trapping model
    conc = conc_AP_arr.iloc[id]
    drug_conc = [conc[i] for i in [-7, -5]]
    for conc_i in drug_conc:
        log = AP_model.custom_simulation(
            param_values, conc_i, 1000,
            timestep=0.1, save_signal=save_signal, abs_tol=abs_tol,
            rel_tol=rel_tol, log_var=['engine.time', 'membrane.V'])

        plt.plot(log.time(), log['membrane.V', 0],
                 label=str(id) + ' ' + str(conc_i))
        # Compute APD90
        APD_trapping_pulse = []
        for pulse in range(save_signal):
            apd90, index = AP_model.APD90(log['membrane.V', pulse], 50, 0.1)
            APD_trapping_pulse.append([apd90, index])

        print(APD_trapping_pulse)

plt.legend()
plt.savefig(fig_dir + 'SD_space_AP.pdf')

plt.figure()
count = 0
check_id = []
for i in range(data_points):
    conc = conc_AP_arr.iloc[i]
    APD = APD_CS_arr.iloc[i]
    conc = conc[~np.isnan(APD)]
    APD = APD[~np.isnan(APD)]
    res = all(i <= j for i, j in zip(APD, APD[1:]))
    diff = APD.values[1:] - APD.values[:-1]
    if not res:
        plt.plot(conc.values[:-1], diff, label=str(count))
        print(APD.values)
        check_id.append(i)
        count += 1
plt.legend()
plt.xscale('log')
plt.savefig(fig_dir + 'CS_space_APD.pdf')
print(count)

# check_id = [check_id[i] for i in [3, 4]]
Hill_model = modelling.HillsModel()
base_conductance = AP_model.original_constants['gKr']

# plt.figure()
# for id in check_id:
#     Hill_curve_coefs = combined_df.iloc[[id]]['Hill_curve']
#     # Run simulation for trapping model
#     conc = conc_AP_arr.iloc[id]
#     # drug_conc = [conc[i] for i in [-7, -5]]
#     for conc_i in drug_conc:
#         reduction_scale = Hill_model.simulate(
#             Hill_curve_coefs, drug_conc[i])
#         d2 = AP_model.conductance_simulation(
#             base_conductance * reduction_scale, 1000,
#             timestep=0.1, save_signal=save_signal, abs_tol=abs_tol,
#             rel_tol=rel_tol, log_var=['engine.time', 'membrane.V'])

#         plt.plot(log.time(), log['membrane.V', 0],
#                  label=str(id) + ' ' + str(conc_i))
#         # Compute APD90
#         APD_trapping_pulse = []
#         for pulse in range(save_signal):
#             apd90, index = AP_model.APD90(log['membrane.V', pulse], 50, 0.1)
#             APD_trapping_pulse.append([apd90, index])

#         print(APD_trapping_pulse)

# plt.legend()
# plt.savefig(fig_dir + 'SD_space_AP.pdf')