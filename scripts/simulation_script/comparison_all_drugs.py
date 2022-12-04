# To compare the Hill curve of drugs for different protocols.
# Generate Hill curve of drug simulations under different protocols.
# Compare the APD90 of the ORd-SD model and the ORd-CS model for each
# synthetic drug
# Output:
# 1. Hill curve of drugs from the SD model under Milnes' protocol stimulation
# 2. The APD90 of the ORd-SD model and the ORd-CS model for a given drug

import myokit
import numpy as np
import os
import pandas as pd
import sys

import modelling

# Define drug
drug = sys.argv[1]

# Define directory to save data
data_dir = '../../simulation_data/model_comparison/' + drug + '/'
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)

# Load IKr model
model = '../../math_model/ohara-cipa-v1-2017-IKr.mmt'
model, _, x = myokit.load(model)
current_model = modelling.BindingKinetics(model)
base_conductance = current_model.model.get('ikr.gKr').value()
current_model.current_head = current_model.model.get('ikr')

# Set up the range of drug concentration
drug_conc_lib = modelling.DrugConcentrations()
drug_conc = drug_conc_lib.drug_concentrations[drug]['coarse']

# Set up protocols
protocol_params = modelling.ProtocolParameters()
protocol_list = protocol_params.protocols
protocols = []
for prot in protocol_list:
    protocols.append(protocol_params.protocol_parameters[prot]['function'])

# Define Hill equation
Hill_model = modelling.HillsModel()
optimiser = modelling.HillsModelOpt(Hill_model)

Hill_coef_df = pd.DataFrame(columns=['Hill coefficient', 'IC50', 'protocol'])
repeats = 1000

for p in range(len(protocols)):
    current_model.protocol = protocols[p]

    # Simulate IKr and compute the peak current
    peaks = []
    for i in range(len(drug_conc)):
        log = current_model.drug_simulation(drug, drug_conc[i], repeats)
        peak, _ = current_model.extract_peak(log, 'ikr.IKr')
        peaks.append(peak[-1])

    # Normalise the peak currents and fit the Hill curve
    peaks = (peaks - min(peaks)) / (max(peaks) - min(peaks))
    estimates, _ = optimiser.optimise(drug_conc, peaks)

    Hill_df = pd.DataFrame({'Hill coefficient': [estimates[0]],
                            'IC50': [estimates[1]],
                            'protocol': [protocol_list[p]]})
    Hill_coef_df = pd.concat([Hill_coef_df, Hill_df])

# Save parameters of Hill curves 
Hill_coef_df.to_csv(data_dir + 'Hill_curves.csv')

# Get Hill curve of the drug stimulated with Milnes' protocol
Hill_eq = Hill_coef_df.loc[Hill_coef_df['protocol'] == protocol_list[0]]
Hill_eq = Hill_eq.values.tolist()[0][:-1]

# Load AP model and set current protocol
APmodel = '../../math_model/ohara-cipa-v1-2017.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')
pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)

# Define constants
base_conductance = APmodel.get('ikr.gKr').value()
offset = 50
save_signal = 2
drug_conc = drug_conc_lib.drug_concentrations[drug]['fine']

APD_conductance = []
APD_trapping = []

abs_tol = 1e-7
rel_tol = 1e-10

for i, conc in enumerate(drug_conc):
    # Simulate AP of the ORd-SD model
    log = AP_model.drug_simulation(
        drug, conc, repeats, save_signal=save_signal,
        abs_tol=abs_tol, rel_tol=rel_tol,
        log_var=['engine.time', 'membrane.V'])

    # Calculate the APD90
    APD_trapping_pulse = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(log['membrane.V', pulse], offset, 0.1)
        APD_trapping_pulse.append(apd90)

    APD_trapping.append(APD_trapping_pulse)

    # Simulate AP of the ORd-CS model
    reduction_scale = Hill_model.simulate(Hill_eq, conc)
    d2 = AP_model.conductance_simulation(
        base_conductance * reduction_scale, repeats,
        save_signal=save_signal, abs_tol=abs_tol, rel_tol=rel_tol,
        log_var=['engine.time', 'membrane.V'])

    # Calculate the APD90
    APD_conductance_pulse = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(d2['membrane.V', pulse], offset, 0.1)
        APD_conductance_pulse.append(apd90)

    APD_conductance.append(APD_conductance_pulse)

# Save APD of the two models
data_dir = data_dir + protocol_list[0] + '/'
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)

APD_trapping_df = pd.DataFrame(np.array(APD_trapping))
APD_trapping_df['drug concentration'] = drug_conc
APD_trapping_df.to_csv(data_dir + 'SD_APD_pulses' + str(save_signal) + '.csv')
APD_conductance_df = pd.DataFrame(np.array(APD_conductance))
APD_conductance_df['drug concentration'] = drug_conc
APD_conductance_df.to_csv(data_dir + 'CS_APD_pulses' + str(save_signal) +
                          '.csv')
