# To compare the Hill equation of models with and without trapping kinetics
# for different protocols.
# Generate Hill equation of drug simulations under different protocols.
# Comparison of CiPA hERG model and hERG model without trapping kinetics.
# Output:
# 1. hERG current under effect of drug, simulated with the Milnes protocol,
# P-80, P0 and P40 protocol.
# 2. Hill curve of drugs from CiPA hERG model and model without trapping
# kinetics, calibrated with Hill curve from Milnes protocol.

import matplotlib
import matplotlib.pyplot as plt
import myokit
import numpy as np
import pandas as pd
import sys

import modelling

drug = sys.argv[1]

data_dir = '../../testing_data/model_comparison/' + drug + '/'

fig_dir = '../../testing_figures/model_comparison/' + drug + '/'

# Load IKr model
model = '../../math_model/ohara-cipa-v1-2017-IKr.mmt'
model, _, x = myokit.load(model)
current_model = modelling.BindingKinetics(model)
base_conductance = current_model.model.get('ikr.gKr').value()
current_model.current_head = current_model.model.get('ikr')

# Set up variables
drug_conc_lib = modelling.DrugConcentrations()
drug_conc = drug_conc_lib.drug_concentrations[drug]['coarse']

protocol_params = modelling.ProtocolParameters()
protocol_list = protocol_params.protocols
protocols = []
for prot in protocol_list:
    protocols.append(protocol_params.protocol_parameters[prot]['function'])
color = ['orange', 'blue', 'red', 'green']

Hill_model = modelling.HillsModel()
optimiser = modelling.HillsModelOpt(Hill_model)

Hill_coef_df = pd.DataFrame(columns=['Hill coefficient', 'IC50', 'protocol'])
repeats = 1000

cmap = matplotlib.cm.get_cmap('viridis')
norm = matplotlib.colors.Normalize(0, len(drug_conc))

for p in range(len(protocols)):
    current_model.protocol = protocols[p]

    # Simulate hERG current
    peaks = []

    for i in range(len(drug_conc)):
        log = current_model.drug_simulation(drug, drug_conc[i], repeats)
        peak, _ = current_model.extract_peak(log, 'ikr.IKr')
        peaks.append(peak[-1])

    plt.figure()
    plt.plot(drug_conc, peaks, 'o')
    plt.xscale("log")
    plt.title(protocol_list[p])
    plt.savefig(fig_dir + protocol_list[p] + ".pdf")
    plt.close()

    peaks = (peaks - min(peaks)) / (max(peaks) - min(peaks))
    estimates, _ = optimiser.optimise(drug_conc, peaks)

    Hill_df = pd.DataFrame({'Hill coefficient': [estimates[0]],
                            'IC50': [estimates[1]],
                            'protocol': [protocol_list[i]]})
    Hill_coef_df = pd.concat([Hill_coef_df, Hill_df])

Hill_coef_df.to_csv(data_dir + 'Hill_curves.csv')

# Simulate and save hERG current under stimulant of Milnes protocol
Hill_eq = Hill_coef_df.loc[Hill_coef_df['protocol'] == protocol_list[0]]
Hill_eq = Hill_eq.values.tolist()[0][:-1]

save_signal = 2
current_model.protocol = protocols[0]
for i, conc in enumerate(drug_conc):
    log = current_model.drug_simulation(
        drug, conc, repeats,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'])
    log.save_csv(data_dir + protocol_list[0] + '/SD_hERG_' + str(conc) +
                 '.csv')

    reduction_scale = Hill_model.simulate(Hill_eq, conc)

    log = current_model.conductance_simulation(
        base_conductance * reduction_scale, repeats, timestep=0.1,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'])
    log.save_csv(data_dir + protocol_list[0] + '/CS_hERG_' + str(conc) +
                 '.csv')

# Set AP model
# Action potential
APmodel = '../../math_model/ohara-cipa-v1-2017.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')
pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)
base_conductance = APmodel.get('ikr.gKr').value()
offset = 50

# Transient phase
repeats = 10
save_signal = repeats

AP_conductance = []
AP_trapping = []

for i, conc in enumerate(drug_conc):
    log = AP_model.drug_simulation(
        drug, conc, repeats, save_signal=save_signal,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'])
    log.save_csv(data_dir + protocol_list[0] + '/SD_AP_transient_pulses' +
                 str(repeats) + '_' + str(conc) + '.csv')

    reduction_scale = Hill_model.simulate(Hill_eq, conc)
    d2 = AP_model.conductance_simulation(
        base_conductance * reduction_scale, repeats,
        save_signal=save_signal,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'])
    d2.save_csv(data_dir + protocol_list[0] + '/CS_AP_transient_pulses' +
                str(repeats) + '_conc' + str(conc) + '.csv')

repeats = 1000
save_signal = 2

for i, conc in enumerate(drug_conc):
    log = AP_model.drug_simulation(
        drug, conc, repeats, save_signal=save_signal,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'])

    log.save_csv(data_dir + protocol_list[0] + '/SD_AP_conc' + str(conc) +
                 '.csv')

    reduction_scale = Hill_model.simulate(Hill_eq, conc)

    d2 = AP_model.conductance_simulation(
        base_conductance * reduction_scale, repeats,
        save_signal=save_signal,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'])

    d2.save_csv(data_dir + protocol_list[0] + '/CS_AP_conc' +
                str(drug_conc[i]) + '.csv')

save_signal = repeats
APD_conductance = []
APD_trapping = []
# No need to run simulation when drug concentration is zero
drug_conc = drug_conc_lib.drug_concentrations[drug]['fine']

for i, conc in enumerate(drug_conc):
    log = AP_model.drug_simulation(
        drug, conc, repeats, save_signal=save_signal,
        log_var=['engine.time', 'membrane.V'])

    APD_trapping_pulse = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(log['membrane.V', pulse], offset, 0.1)
        APD_trapping_pulse.append(apd90)

    APD_trapping.append(APD_trapping_pulse)

    reduction_scale = Hill_model.simulate(Hill_eq, conc)

    d2 = AP_model.conductance_simulation(
        base_conductance * reduction_scale, repeats,
        save_signal=save_signal,
        log_var=['engine.time', 'membrane.V'])

    APD_conductance_pulse = []

    for pulse in range(save_signal):
        apd90 = AP_model.APD90(d2['membrane.V', pulse], offset, 0.1)
        APD_conductance_pulse.append(apd90)

    APD_conductance.append(APD_conductance_pulse)

APD_trapping_df = pd.DataFrame(np.array(APD_trapping))
APD_trapping_df['drug concentration'] = drug_conc
APD_trapping_df.to_csv(data_dir + protocol_list[0] + '/SD_APD_pulses' +
                       str(save_signal) + '.csv')
APD_conductance_df = pd.DataFrame(np.array(APD_conductance))
APD_conductance_df['drug concentration'] = drug_conc
APD_conductance_df.to_csv(data_dir + protocol_list[0] + '/CS_APD_pulses' +
                          str(save_signal) + '.csv')
