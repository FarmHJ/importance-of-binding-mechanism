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

saved_data_dir = '../../simulation_data/binding_kinetics_comparison/' + \
    drug + '/'

testing_fig_dir = '../../figures/testing/'
final_fig_dir = '../../figures/binding_kinetics_comparison/' + drug + '/'

saved_fig_dir = final_fig_dir

# Load IKr model
model = '../../model/ohara-cipa-v1-2017-IKr.mmt'
model, _, x = myokit.load(model)

# Set up variables
drug_conc_lib = modelling.DrugConcentrations()
drug_conc = drug_conc_lib.drug_concentrations[drug]['coarse']

repeats = 1000

protocol_params = modelling.ProtocolParameters()
protocol_list = protocol_params.protocols
protocols = []
pulse_times = []
for prot in protocol_list:
    pulse_times.append(protocol_params.protocol_parameters[prot]['pulse_time'])
    protocols.append(protocol_params.protocol_parameters[prot]['function'])
color = ['orange', 'blue', 'red', 'green']

drug_model = modelling.BindingKinetics(model)

Hill_model = modelling.HillsModel()
optimiser = modelling.HillsModelOpt(Hill_model)

base_conductance = drug_model.model.get('ikr.gKr').value()
drug_model.current_head = drug_model.model.get('ikr')

cmap = matplotlib.cm.get_cmap('viridis')
norm = matplotlib.colors.Normalize(0, len(drug_conc))

peak_combine = []
model_name = []
for p in range(len(protocols)):
    drug_model.protocol = protocols[p]

    # Simulate hERG current
    peaks = []

    for i in range(len(drug_conc)):
        log = drug_model.drug_simulation(drug, drug_conc[i], repeats)
        peak, _ = drug_model.extract_peak(log, 'ikr.IKr')
        peaks.append(peak[-1])

    plt.figure()
    plt.plot(drug_conc, peaks, 'o')
    plt.xscale("log")
    plt.title(protocol_list[p])
    plt.savefig(saved_fig_dir + protocol_list[p] + ".pdf")
    plt.close()
    peak_combine.append(peaks)

Hill_coef_df = pd.DataFrame(columns=['Hill coefficient', 'IC50', 'protocol'])

for i in range(len(peak_combine)):
    peaks = peak_combine[i]
    peaks = (peaks - min(peaks)) / (max(peaks) - min(peaks))
    estimates, _ = optimiser.optimise(drug_conc, peaks)

    Hill_df = pd.DataFrame({'Hill coefficient': [estimates[0]],
                            'IC50': [estimates[1]],
                            'protocol': [protocol_list[i]]})
    Hill_coef_df = pd.concat([Hill_coef_df, Hill_df])

Hill_coef_df.to_csv(saved_data_dir + 'Hill_curves.csv')

# If Hill curves are computed previously
# Hill_coef_df = pd.read_csv(saved_data_dir + 'Hill_curves.csv',
#                            usecols=['Hill coefficient', 'IC50', 'protocol'])

# Simulate and save hERG current under stimulant of Milnes protocol
Hill_eq = Hill_coef_df.loc[Hill_coef_df['protocol'] == protocol_list[0]]
Hill_eq = Hill_eq.values.tolist()[0][:-1]

save_signal = 2
drug_model.protocol = protocols[0]
for i, conc in enumerate(drug_conc):
    log = drug_model.drug_simulation(
        drug, conc, repeats,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'])
    log.save_csv(saved_data_dir + protocol_list[0] + '/CiPA_hERG_' +
                 str(conc) + '.csv')

    reduction_scale = Hill_model.simulate(Hill_eq, conc)

    log = drug_model.conductance_simulation(
        base_conductance * reduction_scale, repeats, timestep=0.1,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'])
    log.save_csv(saved_data_dir + protocol_list[0] + '/conductance_hERG_' +
                 str(conc) + '.csv')

# Set AP model
APmodel = '../../model/ohara-cipa-v1-2017.mmt'

# Action potential
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
        drug, conc, repeats, timestep=0.1, save_signal=save_signal,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'])
    log.save_csv(saved_data_dir + protocol_list[0] +
                 '/CiPA_AP_transient_pulses' + str(repeats) +
                 '_' + str(conc) + '.csv')

    reduction_scale = Hill_model.simulate(Hill_eq, conc)
    d2 = AP_model.conductance_simulation(
        base_conductance * reduction_scale, repeats, timestep=0.1,
        save_signal=save_signal,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'])
    d2.save_csv(saved_data_dir + protocol_list[0] +
                '/conductance_AP_transient_pulses' +
                str(repeats) + '_' + str(conc) + '.csv')

repeats = 1000
save_signal = 2

if drug == 'dofetilide' or 'ranolazine':
    drug_conc = drug_conc[:-1]

for i, conc in enumerate(drug_conc):
    log = AP_model.drug_simulation(
        drug, conc, repeats, timestep=0.1, save_signal=save_signal,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'])

    log.save_csv(saved_data_dir + protocol_list[0] + '/CiPA_AP_' + str(conc) +
                 '.csv')

    reduction_scale = Hill_model.simulate(Hill_eq, conc)

    if drug == 'ranolazine' and conc >= 1e6:
        d2 = AP_model.conductance_simulation(
            base_conductance * reduction_scale, repeats, timestep=0.1,
            save_signal=save_signal, abs_tol=1e-6, rel_tol=1e-5,
            log_var=['engine.time', 'membrane.V', 'ikr.IKr'])
    else:
        d2 = AP_model.conductance_simulation(
            base_conductance * reduction_scale, repeats, timestep=0.1,
            save_signal=save_signal,
            log_var=['engine.time', 'membrane.V', 'ikr.IKr'])

    d2.save_csv(saved_data_dir + protocol_list[0] + '/conductance_AP_' +
                str(drug_conc[i]) + '.csv')

save_signal = repeats
APD_conductance = []
APD_trapping = []
# No need to run simulation when drug concentration is zero
drug_conc = drug_conc_lib.drug_concentrations[drug]['fine']

for i, conc in enumerate(drug_conc):
    log = AP_model.drug_simulation(
        drug, conc, repeats, timestep=0.1, save_signal=save_signal,
        log_var=['engine.time', 'membrane.V'])

    APD_trapping_pulse = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(log['membrane.V', pulse], offset, 0.1)
        APD_trapping_pulse.append(apd90)

    APD_trapping.append(APD_trapping_pulse)

    reduction_scale = Hill_model.simulate(Hill_eq, conc)

    d2 = AP_model.conductance_simulation(
        base_conductance * reduction_scale, repeats, timestep=0.1,
        save_signal=save_signal,
        log_var=['engine.time', 'membrane.V'])

    APD_conductance_pulse = []

    for pulse in range(save_signal):
        apd90 = AP_model.APD90(d2['membrane.V', pulse], offset, 0.1)
        APD_conductance_pulse.append(apd90)

    APD_conductance.append(APD_conductance_pulse)

APD_trapping_df = pd.DataFrame(np.array(APD_trapping))
APD_trapping_df['drug concentration'] = drug_conc
APD_trapping_df.to_csv(saved_data_dir + protocol_list[0] +
                       '/CiPA_APD_pulses' + str(repeats) + '.csv')
APD_conductance_df = pd.DataFrame(np.array(APD_conductance))
APD_conductance_df['drug concentration'] = drug_conc
APD_conductance_df.to_csv(saved_data_dir + protocol_list[0] +
                          '/conductance_APD_pulses' + str(repeats) + '.csv')
