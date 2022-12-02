# To compare the Hill curve of the SD model and the CS model for different
# protocols.
# Generate Hill curve of drug simulations under different protocols.
# Comparison of the SD model and the CS model.
# Output:
# 1. IKr under effect of drug, simulated with Milnes' protocol, Pneg80, P0 and
# P40 protocol.
# 2. Hill curve of drugs from the SD model and the CS model, calibrated with
# Hill curve from Milnes' protocol.

import myokit
import os
import pandas as pd

import modelling

# Load IKr model
model = '../../math_model/ohara-cipa-v1-2017-IKr.mmt'
model, _, x = myokit.load(model)
current_model = modelling.BindingKinetics(model)
current_model.current_head = current_model.model.get('ikr')

# Set up protocols
protocol_params = modelling.ProtocolParameters()
protocol_name = ['Milnes', 'Pneg80', 'P0', 'P40']
protocols = []
for prot in protocol_name:
    protocols.append(protocol_params.protocol_parameters[prot]['function'])
color = ['orange', 'blue', 'red', 'green']

# Set up library of parameters
param_lib = modelling.BindingParameters()
Hill_model = modelling.HillsModel()
drug_conc_lib = modelling.DrugConcentrations()

# Define constants
abs_tol = 1e-7
rel_tol = 1e-8
repeats = 1000

drugs = ['dofetilide', 'verapamil']

for drug in drugs:
    data_dir = '../../simulation_data/model_comparison/' + drug + '/'

    # Set up drug concentration
    drug_conc = drug_conc_lib.drug_concentrations[drug]['coarse']

    Hill_coef_df = pd.DataFrame(columns=['Hill coefficient', 'IC50',
                                         'protocol'])

    # Compute the Hill curves for all protocols given
    for p in range(len(protocols)):
        current_model.protocol = protocols[p]

        # Define the parameter values of the drug
        param_values = param_lib.binding_parameters[drug]
        param_values = pd.DataFrame(data=param_values, index=[0])
        ComparisonController = modelling.ModelComparison(param_values)

        # Generate the Hill curve of the SD model with different protocols
        Hill_coef, _, _ = ComparisonController.compute_Hill(current_model)

        Hill_df = pd.DataFrame({'Hill coefficient': [Hill_coef[0]],
                                'IC50': [Hill_coef[1]],
                                'protocol': [protocol_name[p]]})
        Hill_coef_df = pd.concat([Hill_coef_df, Hill_df])

    # Save parameters of Hill curve -- Hill coefficient and IC50
    Hill_coef_df.to_csv(data_dir + 'Hill_curves.csv')

# Load AP model
APmodel = '../../math_model/ohara-cipa-v1-2017.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')

# Set current protocol
pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)
base_conductance = APmodel.get('ikr.gKr').value()

# Define constants and variables for APD90 simulation of ORd-CS model
drug = 'verapamil'
repeats = 1000
save_signal = 2
offset = 50
drug_conc = drug_conc_lib.drug_concentrations[drug]['fine']

# Define directories to save data
prot_dir = '../../simulation_data/model_comparison/' + drug + '/protocols/'
if not os.path.isdir(prot_dir):
    os.makedirs(prot_dir)

# Simulate APD90 of ORd-CS model whose ionic conductance is calibrated with
# different protocols
for p in protocol_name:

    APD_conductance = []
    # Read Hill curve of given protocol
    Hill_eq = Hill_coef_df.loc[Hill_coef_df['protocol'] == p]
    Hill_eq = Hill_eq.values.tolist()[0][:-1]

    # Simulate AP and calculate APD90
    for i in range(len(drug_conc)):
        print('simulating for drug concentration: ' + str(drug_conc[i]))

        reduction_scale = Hill_model.simulate(Hill_eq, drug_conc[i])
        d2 = AP_model.conductance_simulation(
            base_conductance * reduction_scale, repeats,
            save_signal=save_signal, abs_tol=abs_tol, rel_tol=rel_tol,
            log_var=['engine.time', 'membrane.V'])

        APD_conductance_pulse = []
        for pulse in range(save_signal):
            apd90 = AP_model.APD90(d2['membrane.V', pulse], offset, 0.1)
            APD_conductance_pulse.append(apd90)

        APD_conductance.append(APD_conductance_pulse)

    # Save APD90s
    column_name = ['pulse ' + str(i) for i in range(save_signal)]
    APD_conductance_df = pd.DataFrame(APD_conductance, columns=column_name)
    APD_conductance_df['drug concentration'] = drug_conc
    APD_conductance_df.to_csv(prot_dir + 'CS_APD_' + p + '.csv')
