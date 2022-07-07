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

run_sim = True

drug = 'dofetilide'
protocol_name = 'Milnes'
# pulse_time = 25e3
# protocol = modelling.ProtocolLibrary().Milnes(pulse_time)
drug_conc = [0, 0.1, 1, 30, 100, 300, 500, 1000]

saved_data_dir = '../../simulation_data/binding_kinetics_comparison/' + \
    drug + '/' + protocol_name + '/'
result_filename = 'OHaraCiPA-conductance-fit.txt'

testing_fig_dir = '../../figures/testing/'
final_fig_dir = '../../figures/binding_kinetics_comparison/OHaraCiPA_model/\
    protocol/' + drug + '/'

saved_fig_dir = testing_fig_dir

# Load IKr model
model = '../../model/ohara-cipa-v1-2017-IKr.mmt'
model, _, x = myokit.load(model)

# Set AP model
APmodel = '../../model/ohara-cipa-v1-2017.mmt'

# Set up variables
# pulse_time = 25e3
pulse_times = [25e3, 5400, 5400, 5400]
# drug = 'dofetilide'
# drug_conc = [0, 0.1, 1, 30, 300, 500, 1000] #nM #dofetilide
# drug = 'verapamil'
# drug_conc = [0, 0.1, 1, 30, 300, 500, 1000, 1e4, 1e5] #nM #verapamil
# drug_conc = [30]
repeats = 1000

PMilnes = modelling.ProtocolLibrary().Milnes(pulse_times[0])
Pneg80 = modelling.ProtocolLibrary().Pneg80(pulse_times[1])
P0 = modelling.ProtocolLibrary().P0(pulse_times[1])
P40 = modelling.ProtocolLibrary().P40(pulse_times[1])
protocols = [PMilnes, Pneg80, P0, P40]
protocol_name = ['Milnes', 'Pneg80', 'P0', 'P40']
color = ['orange', 'blue', 'red', 'green', 'black']

drug_model = modelling.BindingKinetics(model)

Hill_model = modelling.HillsModel()
optimiser = modelling.HillsModelOpt(Hill_model)
max_grid = np.ceil(np.log(drug_conc[-1])) + 1
conc_grid = np.arange(-3, max_grid, 1)  # for plotting

params = np.loadtxt(saved_data_dir + result_filename, unpack=True)
params = np.array(params)
Hillcoef = params[0]
IC50 = params[1]

conductance_scale_est = []
for i in range(len(drug_conc)):
    scale = Hill_model.simulate([Hillcoef, IC50], drug_conc[i])
    conductance_scale_est.append(scale)

conductance_scale_df = pd.DataFrame(np.transpose(conductance_scale_est),
                                    columns=[drug])
conductance_scale_df['drug_concentration'] = drug_conc
conductance_scale_df = conductance_scale_df.set_index('drug_concentration')

conductance = drug_model.model.get('ikr.gKr').value()
drug_model.current_head = drug_model.model.get('ikr')

cmap = matplotlib.cm.get_cmap('viridis')
norm = matplotlib.colors.Normalize(0, len(drug_conc))

for p in range(len(protocols)):
    drug_model.protocol = protocols[p]

    fig_log = modelling.figures.FigureStructure(figsize=(5, 4),
                                                gridspec=(3, 1))
    plot_log = modelling.figures.FigurePlot()
    # Simulate hERG current
#     total_log = []
    peaks = []

#     total_log_conduct = []
    peaks_conduct = []
    for i in range(len(drug_conc)):
        log = drug_model.drug_simulation(drug, drug_conc[i], repeats)
        peak, _ = drug_model.extract_peak(log, 'ikr.IKr')
        peaks.append(peak[-1])
        # total_log.append(log)
        plot_log.add_single(fig_log.axs[1][0], log, 'ikr.IKr',
                            color=cmap(norm(i)),
                            label=str(drug_conc[i]) + ' nM')

        scale = conductance_scale_df.iloc[i][drug]
        log_conduct = drug_model.conductance_simulation(conductance * scale,
                                                        repeats)
        peak, _ = drug_model.extract_peak(log_conduct, 'ikr.IKr')
        peaks_conduct.append(peak[-1])
        # total_log_conduct.append(log_conduct)
        plot_log.add_single(fig_log.axs[2][0], log_conduct, 'ikr.IKr',
                            color=cmap(norm(i)),
                            label=str(drug_conc[i]) + ' nM')

    plot_log.add_single(fig_log.axs[0][0], log, 'membrane.V', color='k')
    fig_log.axs[0][0].set_title(protocol_name[p] + ' protocol')
    fig_log.sharex(['Time (s)'], [(0, pulse_times[p])])
    fig_log.sharey(['Voltage', 'Current\n(trapping)',
                    'Current\n(conductance)'])
    fig_log.adjust_ticks(fig_log.axs[2][0], pulse_times[p])
    fig_log.axs[2][0].legend(loc='lower right', bbox_to_anchor=(1.4, 0))
    fig_log.savefig(saved_fig_dir + "hERG_" + protocol_name[p] + "_concs.svg",
                    format='svg')

    if p == 0:
        peak_combine = [peaks_conduct, peaks]
        model_name = ['conductance', 'OHaraCiPA - ' + protocol_name[p]]
    else:
        peak_combine += [peaks]
        model_name += ['OHaraCiPA - ' + protocol_name[p]]

fig_Hill = modelling.figures.FigureStructure(figsize=(4, 3))
plot_Hill = modelling.figures.FigurePlot()

for i in range(len(peak_combine)):
    peaks = peak_combine[i]
    peaks = (peaks - min(peaks)) / (max(peaks) - min(peaks))
    estimates, _ = optimiser.optimise(drug_conc, peaks)

    fig_Hill.axs[0][0].plot(np.log(drug_conc[1:]), peaks[1:], 'o',
                            color=color[i])
    fig_Hill.axs[0][0].plot(conc_grid, Hill_model.simulate(
        estimates[:2], np.exp(conc_grid)), '-', color=color[i],
        label=model_name[i])
    fig_Hill.axs[0][0].text(
        -3, 0.35 - i * 0.1,
        model_name[i] + r': $n=$%.2f' % (estimates[0]) + "\n" +
        r'IC$_{50}=$%.2f' % (estimates[1]), fontsize=7)  # , wrap=True)

fig_Hill.axs[0][0].set_xlabel('Drug concentration (log)')
fig_Hill.axs[0][0].set_ylabel('Normalised peak current')
fig_Hill.fig.legend()
# fig_Hill.fig.tight_layout(pad=0.4)
fig_Hill.savefig(saved_fig_dir + "peak_hERG_Hill_protocols.pdf")

# Action potential
# APmodel, _, x = myokit.load(APmodel)
# AP_model = modelling.BindingKinetics(APmodel)
# pulse_time = 1000
# AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)
# conductance = APmodel.get('ikr.gKr').value()

# repeats = 1000
# plotting_pulse_time = 750

# # Simulate AP
# log2 = AP_model.drug_simulation(drug, drug_conc, repeats,
#         current_head='ikr', timestep=0.1)

# # Plot AP
# fig = modelling.figures.CurrentPlot(AP_model)
# fig.add_plot_AP(log2, plotting_pulse_time)
# plt.savefig(saved_fig_dir + "AP_OHaraCiPA_Milnes.pdf")
