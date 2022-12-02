# Simulate action potentials at transient phase
import myokit
import numpy as np
import os
import pandas as pd
import pints
import sys

import modelling

# Define drug and protocol
drug = sys.argv[1]
protocol_name = 'Milnes'
protocol_params = modelling.ProtocolParameters()
pulse_time = protocol_params.protocol_parameters[protocol_name]['pulse_time']
protocol = protocol_params.protocol_parameters[protocol_name]['function']

# Define drug concentration range for each drug of interest
if drug == 'dofetilide':
    drug_conc = [0, 100, 200]  # nM
elif drug == 'verapamil':
    drug_conc = [0, 1000, 10000]  # nM

# Define directories to save simulated data
root_dir = '../../simulation_data/model_comparison/'
data_dir = root_dir + drug + '/' + protocol_name + '/'
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)
result_filename = 'Hill_curve.txt'

# Load IKr model
model = '../../math_model/ohara-cipa-v1-2017-IKr.mmt'
model, _, x = myokit.load(model)

current_model = modelling.BindingKinetics(model)
current_model.protocol = protocol

# Load AP model
APmodel = '../../math_model/ohara-cipa-v1-2017.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')

# Set up current protocol
pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)
base_conductance = APmodel.get('ikr.gKr').value()

# Load Hill curve parameters or compute them if not previously done
if not os.path.isfile(data_dir + result_filename):
    param_lib = modelling.BindingParameters()
    param_values = param_lib.binding_parameters[drug]
    param_values = pd.DataFrame(data=param_values, index=[0])
    ComparisonController = modelling.ModelComparison(param_values)
    estimates, _, _ = ComparisonController.compute_Hill(current_model)
    with open(data_dir + result_filename, 'w') as f:
        for x in estimates:
            f.write(pints.strfloat(x) + '\n')
else:
    estimates = np.loadtxt(data_dir + result_filename, unpack=True)
    estimates = np.array(estimates)
Hill_model = modelling.HillsModel()

offset = 50
repeats = 7
save_signal = repeats

# Get steady state AP for control condition
if not os.path.isfile(root_dir + 'steady_state_control.csv'):
    control_log = AP_model.drug_simulation(drug, 0, 1000)
    control_log.save_csv(root_dir + 'steady_state_control.csv')
else:
    control_log = myokit.DataLog.load_csv(root_dir +
                                          'steady_state_control.csv')

# Simulate AP after adding drugs for both the ORd-SD model and the ORd-CS model
# Repeated for 7 pulses
AP_conductance = []
AP_trapping = []

for i in range(len(drug_conc)):
    print('simulating for drug concentration: ' + str(drug_conc[i]))
    log = AP_model.drug_simulation(
        drug, drug_conc[i], repeats, save_signal=save_signal,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'],
        set_state=control_log)
    log.save_csv(data_dir + 'SD_AP_transient_pulses' + str(repeats) +
                 '_conc' + str(drug_conc[i]) + '_paced.csv')

    reduction_scale = Hill_model.simulate(estimates[:2], drug_conc[i])
    d2 = AP_model.conductance_simulation(
        base_conductance * reduction_scale, repeats,
        save_signal=save_signal,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'],
        set_state=control_log)
    d2.save_csv(data_dir + 'CS_AP_transient_pulses' + str(repeats) +
                '_conc' + str(drug_conc[i]) + '_paced.csv')

# Remove drug free condition
drug_conc = drug_conc[1:]

APD_conductance = []
APD_trapping = []

# Simulate the AP models for 300 pulses to show transition of APD to steady state
repeats = 300
save_signal = repeats

for i in range(len(drug_conc)):
    print('simulating concentration: ' + str(drug_conc[i]))
    log = AP_model.drug_simulation(
        drug, drug_conc[i], repeats, save_signal=save_signal,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'],
        set_state=control_log)

    APD_trapping_pulse = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(log['membrane.V', pulse], offset, 0.1)
        APD_trapping_pulse.append(apd90)

    APD_trapping.append(APD_trapping_pulse)

    reduction_scale = Hill_model.simulate(estimates[:2], drug_conc[i])
    d2 = AP_model.conductance_simulation(
        base_conductance * reduction_scale, repeats,
        save_signal=save_signal,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'],
        set_state=control_log)

    APD_conductance_pulse = []

    for pulse in range(save_signal):
        apd90 = AP_model.APD90(d2['membrane.V', pulse], offset, 0.1)
        APD_conductance_pulse.append(apd90)

    APD_conductance.append(APD_conductance_pulse)

# Save simulated APD
APD_trapping_df = pd.DataFrame(np.array(APD_trapping))
APD_trapping_df['drug concentration'] = drug_conc
APD_trapping_df.to_csv(data_dir + 'SD_APD_transient_paced.csv')
APD_conductance_df = pd.DataFrame(np.array(APD_conductance))
APD_conductance_df['drug concentration'] = drug_conc
APD_conductance_df.to_csv(data_dir + 'CS_APD_transient_paced.csv')
