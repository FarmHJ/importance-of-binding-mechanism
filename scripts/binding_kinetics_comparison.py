#
# Calibrate the ionic conductance of the CS model from the SD model and
# compare the APD90 at steady state.
# Output:
# 1. IKr at various drug concentration simulated from the SD model.
# 2. Fitted Hill curve over peak IKr from the SD model.
# 3. IKr at various drug concentration simulated from the corresponding
#    CS model.
# 4. 2 pulses of action potential and their APD90s simulated from both models
#    (steady state).
# 5. APD90 simulated at various drug concentrations for both models (steady
#    state).
#

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

# Define the range of drug concentration for a given drug
drug_conc_lib = modelling.DrugConcentrations()
drug_conc = drug_conc_lib.drug_concentrations[drug]['coarse']
repeats = 1000

# Define directories to save simulated data
data_dir = '../simulation_data/model_comparison/' + \
    drug + '/' + protocol_name + '/'
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)
result_filename = 'Hill_curve.txt'

# Load IKr model
model = '../math_model/ohara-cipa-v1-2017-IKr-opt.mmt'
model, _, x = myokit.load(model)

current_model = modelling.BindingKinetics(model)
current_model.protocol = protocol

# Define tolerance value
abs_tol = 1e-7
rel_tol = 1e-8

# Simulate IKr of the SD model for a range of drug concentrations
# Extract the peak of IKr
peaks = []
for i in range(len(drug_conc)):
    log = current_model.drug_simulation(
        drug, drug_conc[i], repeats,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'], abs_tol=abs_tol,
        rel_tol=rel_tol)
    peak, _ = current_model.extract_peak(log, 'ikr.IKr')
    peaks.append(peak[-1])

    log.save_csv(data_dir + 'SD_current_' + str(drug_conc[i]) + '.csv')

# Normalise drug response (peak current)
peaks = (peaks - min(peaks)) / (max(peaks) - min(peaks))

# Fit drug response to Hill curve
Hill_model = modelling.HillsModel()
optimiser = modelling.HillsModelOpt(Hill_model)
if not os.path.isfile(data_dir + result_filename):
    estimates, _ = optimiser.optimise(drug_conc, peaks)
    with open(data_dir + result_filename, 'w') as f:
        for x in estimates:
            f.write(pints.strfloat(x) + '\n')
else:
    estimates = np.loadtxt(data_dir + result_filename, unpack=True)
    estimates = np.array(estimates)

# Compare peak current
base_conductance = model.get('ikr.gKr').value()
current_model.current_head = current_model.model.get('ikr')
for i in range(len(drug_conc)):

    reduction_scale = Hill_model.simulate(estimates[:2], drug_conc[i])
    d2 = current_model.conductance_simulation(
        base_conductance * reduction_scale, repeats,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'],
        abs_tol=abs_tol, rel_tol=rel_tol)

    d2.save_csv(data_dir + 'CS_current_' + str(drug_conc[i]) + '.csv')

#
# Propagate to action potential
#

# Set AP model
APmodel = '../math_model/ohara-cipa-v1-2017-opt.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')

# Define current protocol
pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)
base_conductance = APmodel.get('ikr.gKr').value()
protocol_duration = AP_model.protocol.characteristic_time()

offset = 50
save_signal = 2
# Use different repeats for plotting purpose - so that EAD-like behaviour
# happens on the same pulse
if drug == 'dofetilide':
    repeats_SD = 1000 + 1
    repeats_CS = 1000
else:
    repeats_SD = 1000
    repeats_CS = 1000

AP_conductance = []
AP_trapping = []
APD_conductance = []
APD_trapping = []

# Simulate AP of the AP-SD model and the AP-CS model
# Compute APD90
for i in range(len(drug_conc)):
    print('simulating concentration: ' + str(drug_conc[i]))
    log = AP_model.drug_simulation(
        drug, drug_conc[i], repeats_SD, save_signal=save_signal,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'], abs_tol=abs_tol,
        rel_tol=rel_tol)
    log.save_csv(data_dir + 'SD_AP_' + str(drug_conc[i]) + '.csv')

    # APD_trapping_pulse = []
    # for pulse in range(save_signal):
    #     apd90 = AP_model.APD90(log['membrane.V', pulse], offset, 0.1)
    #     APD_trapping_pulse.append(apd90)

    Vm_signal = list(log['membrane.V', 0])
    for pulse in range(1, save_signal):
        Vm_signal += list(log['membrane.V', pulse])
    APD_trapping_pulse = AP_model.APD90_update(
        log.time(), Vm_signal, offset, protocol_duration)

    AP_trapping.append(log)
    APD_trapping.append(APD_trapping_pulse)

    reduction_scale = Hill_model.simulate(estimates[:2], drug_conc[i])
    d2 = AP_model.conductance_simulation(
        base_conductance * reduction_scale, repeats_CS,
        save_signal=save_signal,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'], abs_tol=abs_tol,
        rel_tol=rel_tol)
    d2.save_csv(data_dir + 'CS_AP_' + str(drug_conc[i]) + '.csv')

    # APD_conductance_pulse = []
    # for pulse in range(save_signal):
    #     apd90 = AP_model.APD90(d2['membrane.V', pulse], offset, 0.1)
    #     APD_conductance_pulse.append(apd90)

    Vm_signal = list(d2['membrane.V', 0])
    for pulse in range(1, save_signal):
        Vm_signal += list(d2['membrane.V', pulse])
    APD_conductance_pulse = AP_model.APD90_update(
        d2.time(), Vm_signal, offset, protocol_duration)

    AP_conductance.append(d2)
    APD_conductance.append(APD_conductance_pulse)

    print('done concentration: ' + str(drug_conc[i]))

# Save simulated APD90 of both the AP-SD model and the AP-CS model
column_name = ['pulse ' + str(i) for i in range(save_signal)]
APD_trapping_df = pd.DataFrame(APD_trapping, columns=column_name)
APD_trapping_df['drug concentration'] = drug_conc
APD_trapping_df.to_csv(data_dir + 'SD_APD_pulses' +
                       str(int(save_signal)) + '.csv')
APD_conductance_df = pd.DataFrame(APD_conductance, columns=column_name)
APD_conductance_df['drug concentration'] = drug_conc
APD_conductance_df.to_csv(data_dir + 'CS_APD_pulses' +
                          str(int(save_signal)) + '.csv')

# Define drug concentration range for steady state APD90 comparison between
# models
drug_conc = drug_conc_lib.drug_concentrations[drug]['fine']
repeats = 1000
save_signal = 2

APD_conductance = []
APD_trapping = []

for i in range(len(drug_conc)):
    print('simulating concentration: ' + str(drug_conc[i]))

    # Run simulation for the AP-SD model till steady state
    log = AP_model.drug_simulation(
        drug, drug_conc[i], repeats, save_signal=save_signal,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'], abs_tol=abs_tol,
        rel_tol=rel_tol)

    # Compute APD90 of simulated AP
    # APD_trapping_pulse = []
    # for pulse in range(save_signal):
    #     apd90 = AP_model.APD90(log['membrane.V', pulse], offset, 0.1)
    #     APD_trapping_pulse.append(apd90)

    Vm_signal = list(log['membrane.V', 0])
    for pulse in range(1, save_signal):
        Vm_signal += list(log['membrane.V', pulse])
    APD_trapping_pulse = AP_model.APD90_update(
        log.time(), Vm_signal, offset, protocol_duration)

    APD_trapping.append(APD_trapping_pulse)

    # Run simulation for the AP-CS model till steady state
    reduction_scale = Hill_model.simulate(estimates[:2], drug_conc[i])
    d2 = AP_model.conductance_simulation(
        base_conductance * reduction_scale, repeats,
        save_signal=save_signal,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'], abs_tol=abs_tol,
        rel_tol=rel_tol)

    # Compute APD90 of simulated AP
    # APD_conductance_pulse = []
    # for pulse in range(save_signal):
    #     apd90 = AP_model.APD90(d2['membrane.V', pulse], offset, 0.1)
    #     APD_conductance_pulse.append(apd90)

    Vm_signal = list(d2['membrane.V', 0])
    for pulse in range(1, save_signal):
        Vm_signal += list(d2['membrane.V', pulse])
    APD_conductance_pulse = AP_model.APD90_update(
        d2.time(), Vm_signal, offset, protocol_duration)

    APD_conductance.append(APD_conductance_pulse)

    print('done concentration: ' + str(drug_conc[i]))

# Compute APD90 with AP behaviour in alternating cycles
# APD_trapping = [max(i) for i in APD_trapping]
# APD_conductance = [max(i) for i in APD_conductance]

APD_trapping = [float('nan') if np.isnan(i).any() else max(i)
                for i in APD_trapping]
APD_conductance = [float('nan') if np.isnan(i).any() else max(i)
                   for i in APD_conductance]

# Save APD90 data
APD_trapping_df = pd.DataFrame(np.array(APD_trapping), columns=['APD'])
APD_trapping_df['drug concentration'] = drug_conc
APD_trapping_df.to_csv(data_dir + 'SD_APD_fine.csv')
APD_conductance_df = pd.DataFrame(np.array(APD_conductance), columns=['APD'])
APD_conductance_df['drug concentration'] = drug_conc
APD_conductance_df.to_csv(data_dir + 'CS_APD_fine.csv')
