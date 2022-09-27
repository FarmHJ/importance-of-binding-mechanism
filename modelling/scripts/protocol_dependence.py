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
import myokit
import numpy as np
import pandas as pd

import modelling

drug = 'verapamil'

saved_data_dir = '../../simulation_data/binding_kinetics_comparison/' + \
    drug + '/'

testing_fig_dir = '../../figures/testing/'
final_fig_dir = \
    '../../figures/binding_kinetics_comparison/OHaraCiPA_model/protocol/' + \
    drug + '/'

saved_fig_dir = testing_fig_dir

# Load IKr model
model = '../../model/ohara-cipa-v1-2017-IKr.mmt'
model, _, x = myokit.load(model)

# Set AP model
APmodel = '../../model/ohara-cipa-v1-2017.mmt'

# Set up variables
# pulse_time = 25e3
pulse_times = [25e3, 5400, 5400, 5400]
if drug == 'dofetilide':
    drug_conc = [0, 0.1, 1, 30, 100, 300, 500, 1000]  # nM
elif drug == 'verapamil':
    drug_conc = [0, 0.1, 1, 30, 300, 500, 1000, 10000, 1e5]  # nM

repeats = 1000

protocol_params = modelling.ProtocolParameters()
PMilnes = modelling.ProtocolLibrary().Milnes(pulse_times[0])
Pneg80 = modelling.ProtocolLibrary().Pneg80(pulse_times[1])
P0 = modelling.ProtocolLibrary().P0(pulse_times[1])
P40 = modelling.ProtocolLibrary().P40(pulse_times[1])
protocols = [PMilnes, Pneg80, P0, P40]
protocol_name = ['Milnes', 'Pneg80', 'P0', 'P40']
color = ['orange', 'blue', 'red', 'green']

drug_model = modelling.BindingKinetics(model)

Hill_model = modelling.HillsModel()
optimiser = modelling.HillsModelOpt(Hill_model)
max_grid = np.ceil(np.log(drug_conc[-1])) + 1
conc_grid = np.arange(-3, max_grid, 1)  # for plotting

conductance = drug_model.model.get('ikr.gKr').value()
drug_model.current_head = drug_model.model.get('ikr')

cmap = matplotlib.cm.get_cmap('viridis')
norm = matplotlib.colors.Normalize(0, len(drug_conc))

peak_combine = []
model_name = []
for p in range(len(protocols)):
    drug_model.protocol = protocols[p]

    fig_protocol = modelling.figures.FigureStructure(figsize=(2, 0.7))
    plot_protocol = modelling.figures.FigurePlot()
    # Simulate hERG current
    peaks = []

    for i in range(len(drug_conc)):
        log = drug_model.drug_simulation(drug, drug_conc[i], repeats)
        peak, _ = drug_model.extract_peak(log, 'ikr.IKr')
        peaks.append(peak[-1])

    plot_protocol.add_single(fig_protocol.axs[0][0], log, 'membrane.V',
                             color='k', label=protocol_name[p])
    fig_protocol.sharex([' '], [(0, pulse_times[p])])
    fig_protocol.axs[0][0].set_yticks(protocol_params.protocol_parameters[
        protocol_name[p]]['voltage_points'])
    fig_protocol.axs[0][0].set_title(protocol_name[p] + " protocol")
    fig_protocol.axs[0][0].spines['top'].set_visible(False)
    fig_protocol.axs[0][0].spines['right'].set_visible(False)
    fig_protocol.savefig(saved_fig_dir + "../" + protocol_name[p] + ".pdf")

    peak_combine.append(peaks)

fig_Hill = modelling.figures.FigureStructure(figsize=(4, 3))
plot_Hill = modelling.figures.FigurePlot()

Hill_coef_df = pd.DataFrame(columns=['Hill coefficient', 'IC50', 'protocol'])

for i in range(len(peak_combine)):
    peaks = peak_combine[i]
    peaks = (peaks - min(peaks)) / (max(peaks) - min(peaks))
    estimates, _ = optimiser.optimise(drug_conc, peaks)

    Hill_df = pd.DataFrame({'Hill coefficient': [estimates[0]],
                            'IC50': [estimates[1]],
                            'protocol': [protocol_name[i]]})
    Hill_coef_df = pd.concat([Hill_coef_df, Hill_df])

    fig_Hill.axs[0][0].plot(np.log(drug_conc[1:]), peaks[1:], 'o',
                            color=color[i])
    fig_Hill.axs[0][0].plot(conc_grid, Hill_model.simulate(
        estimates[:2], np.exp(conc_grid)), '-', color=color[i],
        label=protocol_name[i])
    # fig_Hill.axs[0][0].text(
    #     -3, 0.5 - i * 0.13,
    #     model_name[i] + r': $n=$%.2f' % (estimates[0]) + "\n" +
    #     r'IC$_{50}=$%.2f' % (estimates[1]), fontsize=7)  # , wrap=True)

fig_Hill.axs[0][0].set_xlabel('Drug concentration (log)')
fig_Hill.axs[0][0].set_ylabel('Normalised peak current')
fig_Hill.axs[0][0].legend()
fig_Hill.savefig(saved_fig_dir + "peak_hERG_Hill_protocols.pdf")

Hill_coef_df.to_csv(saved_data_dir + 'Hill_curves.csv')

# Action potential
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')
pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)
base_conductance = APmodel.get('ikr.gKr').value()

repeats = 1000
save_signal = 2
offset = 50
if drug == 'dofetilide':
    drug_conc = 10.0**np.linspace(-1, 2.5, 20)
elif drug == 'verapamil':
    drug_conc = 10.0**np.linspace(-1, 5, 20)
    drug_conc = drug_conc[-2:]
    protocol_name = ['P0', 'P40']

for p in protocol_name:
    saved_data_dir = '../../simulation_data/binding_kinetics_comparison/' + \
        drug + '/' + p + '/'
    APD_conductance = []
    Hill_eq = Hill_coef_df.loc[Hill_coef_df['protocol'] == p]
    Hill_eq = Hill_eq.values.tolist()[0][:-1]
    print(p)
    if p == 'P0':
        repeats = 150
    elif p == 'P40':
        repeats = 175
    else:
        repeats = 5

    for i in range(len(drug_conc)):
        print(drug_conc[i])

        reduction_scale = Hill_model.simulate(Hill_eq, drug_conc[i])
        d2 = AP_model.conductance_simulation(
            base_conductance * reduction_scale, repeats, timestep=0.1,
            save_signal=save_signal,
            log_var=['engine.time', 'membrane.V'])

        APD_conductance_pulse = []
        for pulse in range(save_signal):
            apd90 = AP_model.APD90(d2['membrane.V', pulse], offset, 0.1)
            APD_conductance_pulse.append(apd90)

        APD_conductance.append(APD_conductance_pulse)

    column_name = ['pulse ' + str(i) for i in range(save_signal)]
    APD_conductance_df = pd.DataFrame(APD_conductance, columns=column_name)
    APD_conductance_df['drug concentration'] = drug_conc
    APD_conductance_df.to_csv(saved_data_dir + 'conductance_APD_pulses' +
                              str(int(save_signal)) + '.csv')
