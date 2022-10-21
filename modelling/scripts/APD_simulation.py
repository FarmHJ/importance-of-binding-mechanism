# Simulate action potentials and calculate their APD90
import myokit
import numpy as np
import pandas as pd
import sys

import modelling


drug = sys.argv[1]
protocol_name = sys.argv[2]
protocol_params = modelling.ProtocolParameters()
pulse_time = protocol_params.protocol_parameters[protocol_name]['pulse_time']
protocol = protocol_params.protocol_parameters[protocol_name]['function']
##########################
# Drug concentration can be set referencing to hERG_peak figure
if drug == 'dofetilide':
    drug_conc = 10.0**np.linspace(-1, 2.5, 20)
    # drug_conc_fine = 10**np.linspace(2.132, 2.316, 5)
    # drug_conc = np.append(drug_conc, drug_conc_fine)
elif drug == 'verapamil':
    drug_conc = 10.0**np.linspace(-1, 5, 20)
repeats = 1000
drug_labels = [str(i) + ' nM' for i in drug_conc]

saved_data_dir = '../../simulation_data/binding_kinetics_comparison/' + \
    drug + '/' + protocol_name + '/'
result_filename = 'OHaraCiPA-conductance-fit.txt'

# Set AP model
APmodel = '../../model/ohara-cipa-v1-2017.mmt'

# Load Hill curve parameters
estimates = np.loadtxt(saved_data_dir + result_filename, unpack=True)
estimates = np.array(estimates)
Hill_model = modelling.HillsModel()

# Check action potential
# Load model
APmodel, _, x = myokit.load(APmodel)

AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')
pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)
base_conductance = APmodel.get('ikr.gKr').value()

offset = 50
if drug == 'verapamil' and protocol_name == 'P0':
    repeats = 83
else:
    repeats = 1000
save_signal = 2

APD_conductance = []
APD_trapping = []

for i in range(len(drug_conc)):
    print('simulating concentration: ' + str(drug_conc[i]))
    log = AP_model.drug_simulation(
        drug, drug_conc[i], repeats, timestep=0.1, save_signal=save_signal,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'])

    APD_trapping_pulse = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(log['membrane.V', pulse], offset, 0.1)
        APD_trapping_pulse.append(apd90)

    APD_trapping.append(APD_trapping_pulse)

    reduction_scale = Hill_model.simulate(estimates[:2], drug_conc[i])
    d2 = AP_model.conductance_simulation(
        base_conductance * reduction_scale, repeats, timestep=0.1,
        save_signal=save_signal,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'])

    APD_conductance_pulse = []

    for pulse in range(save_signal):
        apd90 = AP_model.APD90(d2['membrane.V', pulse], offset, 0.1)
        APD_conductance_pulse.append(apd90)

    APD_conductance.append(APD_conductance_pulse)

    print('done concentration: ' + str(drug_conc[i]))

print(APD_trapping)
print(APD_conductance)
APD_trapping = [max(i) for i in APD_trapping]
APD_conductance = [max(i) for i in APD_conductance]

APD_trapping_df = pd.DataFrame(np.array(APD_trapping), columns=['APD'])
APD_trapping_df['drug concentration'] = drug_conc
APD_trapping_df.to_csv(saved_data_dir + 'CiPA_APD_fine.csv')
APD_conductance_df = pd.DataFrame(np.array(APD_conductance), columns=['APD'])
APD_conductance_df['drug concentration'] = drug_conc
APD_conductance_df.to_csv(saved_data_dir + 'conductance_APD_fine.csv')
