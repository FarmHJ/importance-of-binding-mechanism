# Calibrate the conductance reduction model from synthetic experimental data
# (CiPA hERG model) and compare the APD90 at steady state and transient phase.
# Output:
# 1. hERG current at various drug concentration simulated from CiPA hERG model.
# 2. Normalised peak hERG current vs drug concentration in log scale.
# 3. Fitted Hill curve over peak hERG current from CiPA hERG model.
# 4. hERG current at various drug concentration simulated from conductance
#    calibrated hERG model.
# 5. Comparison of peak hERG current between both models, i.e. CiPA hERG model
#    and conductance hERG model.
# 6. Comparison of hERG current between both models.
# 7. 2 pulses of action potential simulated from both models (steady state).
# 8. APD90 of both pulses for both models (steady state).
# 9. 10 pulses of action potentials simulated from both models
#    (transient phase).
# 10. APD90 of all pulses for both models (transient phase).

import matplotlib
import matplotlib.pyplot as plt
import myokit
import numpy as np
import pandas as pd
import pints

import modelling

run_sim = True
steady_state = False
plot_fig = True

drug = 'dofetilide'
protocol_name = 'P0'
protocol_params = modelling.ProtocolParameters()
pulse_time = protocol_params.protocol_parameters[protocol_name]['pulse_time']
protocol = protocol_params.protocol_parameters[protocol_name]['function']
if drug == 'dofetilide':
    drug_conc = [0, 0.1, 1, 30, 100, 300, 500, 1000]  # nM
elif drug == 'verapamil':
    # drug_conc = [0, 0.1, 1, 30, 300, 500, 1000, 10000, 1e5]  # nM
    drug_conc = [0, 0.1, 1, 30, 300, 1000, 10000, 1e5]
repeats = 1000
drug_labels = [str(i) + ' nM' for i in drug_conc]

saved_data_dir = '../../simulation_data/binding_kinetics_comparison/' + \
    drug + '/' + protocol_name + '/'
result_filename = 'OHaraCiPA-conductance-fit.txt'

testing_fig_dir = '../../figures/testing/'
final_fig_dir = '../../figures/binding_kinetics_comparison/' + drug + '/' + \
    protocol_name + '/'

saved_fig_dir = final_fig_dir

# Load IKr model
model = '../../model/ohara-cipa-v1-2017-IKr.mmt'
model, _, x = myokit.load(model)

# Set AP model
APmodel = '../../model/ohara-cipa-v1-2017.mmt'

drug_model = modelling.BindingKinetics(model)
drug_model.protocol = protocol

# Run simulation
total_log = []
peaks = []
for i in range(len(drug_conc)):
    log = drug_model.drug_simulation(drug, drug_conc[i], repeats)
    peak, _ = drug_model.extract_peak(log, 'ikr.IKr')
    peaks.append(peak[-1])
    total_log.append(log)

# Plot hERG currents for different drug concentrations
if plot_fig:
    fig_plot = modelling.figures.ReferenceStructure()
    fig = fig_plot.current_concs(total_log, pulse_time, drug_conc)
    fig.savefig(saved_fig_dir + "hERG_trapping_" + drug + "_concs.pdf")

# Plot drug response against drug concentration curve
peaks = (peaks - min(peaks)) / (max(peaks) - min(peaks))

plt.rcParams.update({'font.size': 9})

if plot_fig:
    plt.figure(figsize=(4, 3))
    plt.plot(np.log(drug_conc), peaks, 'o')
    plt.xlabel('Drug concentration (log)')
    plt.ylabel('Normalised peak current')
    plt.tight_layout()
    plt.savefig(saved_fig_dir + "peak_hERG_trapping_" + drug + "_concs.pdf")

# Fit Hill curve
Hill_model = modelling.HillsModel()
optimiser = modelling.HillsModelOpt(Hill_model)
if run_sim:
    estimates, _ = optimiser.optimise(drug_conc, peaks)
    with open(saved_data_dir + result_filename, 'w') as f:
        for x in estimates:
            f.write(pints.strfloat(x) + '\n')
else:
    estimates = np.loadtxt(saved_data_dir + result_filename, unpack=True)
    estimates = np.array(estimates)

# Plot fitting result of Hill curve
max_grid = np.ceil(np.log(drug_conc[-1]))
conc_grid = np.arange(-3, max_grid, 1)

if plot_fig:
    plt.figure(figsize=(4, 3))
    plt.plot(np.log(drug_conc[1:]), peaks[1:], 'o', label='peak current')
    plt.plot(conc_grid, Hill_model.simulate(estimates[:2], np.exp(conc_grid)),
             'k', label='fitted Hill eq')
    plt.xlabel('Drug concentration (log)')
    plt.ylabel('Normalised peak current')
    plt.tight_layout()
    plt.legend()
    plt.savefig(saved_fig_dir + "fitted_peak_hERG_trapping_" + drug
                + "_concs.pdf")

# Scale peak current for simple conductance model (without trapping)
conductance_scale_est = []
for i in range(len(drug_conc)):
    scale = Hill_model.simulate(estimates[:2], drug_conc[i])
    conductance_scale_est.append(scale)

# Reformat to dataframe
conductance_scale_df = pd.DataFrame(np.transpose(conductance_scale_est),
                                    columns=[drug])
conductance_scale_df['drug_concentration'] = drug_conc
conductance_scale_df = conductance_scale_df.set_index('drug_concentration')

# Compare peak current
conductance = model.get('ikr.gKr').value()

peaks_conductance = []
log_conductance = []

peaks_trapping = []
log_trapping = []

drug_model.current_head = drug_model.model.get('ikr')
for i in range(len(drug_conc)):
    log = drug_model.drug_simulation(drug, drug_conc[i], repeats)
    log_trapping.append(log)

    peaks, _ = drug_model.extract_peak(log, 'ikr.IKr')
    peaks_trapping.append(peaks[-1])

    scale = conductance_scale_df.iloc[i][drug]
    d2 = drug_model.conductance_simulation(conductance * scale, repeats)
    # , current_head='ikr')
    log_conductance.append(d2)

    peaks, _ = drug_model.extract_peak(d2, 'ikr.IKr')
    peaks_conductance.append(peaks[-1])

if not plot_fig:
    # fig = modelling.figures.CurrentPlot(drug_model)
    # fig.add_plot_current_various(log_conductance, drug_conc, pulse_time)
    # fig.adjust_ticks(fig.axs[1], pulse_time)
    # plt.savefig(saved_fig_dir + "hERG_conductance_" + drug + "_concs.pdf")

    fig = modelling.figures.FigureStructure(figsize=(5, 4),
                                            gridspec=(3, 1))
    plot = modelling.figures.FigurePlot()
    cmap = matplotlib.cm.get_cmap('viridis')

    labels = [str(i) + ' nM' for i in drug_conc]
    plot.add_single(fig.axs[0][0], log_trapping[0], 'membrane.V', color='k')
    plot.add_multiple(fig.axs[1][0], log_trapping, 'ikr.IKr', labels=labels,
                      color=cmap)
    plot.add_multiple(fig.axs[2][0], log_conductance, 'ikr.IKr', labels=labels,
                      color=cmap)

    fig.axs[2][0].legend()
    fig.sharex(['Time (s)'], [(0, pulse_time)])
    fig.sharey(['Voltage (mV)', 'Current (A/F)\ntrapped model',
               'Current (A/F)\nnontrapped model'])
    fig.adjust_ticks(fig.axs[2][0], pulse_time)
    fig.savefig(saved_fig_dir + 'combine_comparison_viridis.pdf')

# Check hERG current
if plot_fig:
    row, col = 3, 3
    fig = modelling.figures.ReferenceStructure()
    fig.hERG_compare(log_trapping, log_conductance, drug_conc, pulse_time,
                     grid=(row, col))
    plt.savefig(saved_fig_dir + "hERG_compare_" + drug + "_concs.pdf")

# Check hERG peak current
peaks_trapping_norm = (peaks_trapping - min(peaks_trapping)) / (
    max(peaks_trapping) - min(peaks_trapping))
peaks_conductance_norm = (peaks_conductance - min(peaks_conductance)) / (
    max(peaks_conductance) - min(peaks_conductance))

if plot_fig:
    plt.figure(figsize=(4, 3))
    plt.plot(np.log(drug_conc[1:]), peaks_trapping_norm[1:],
             'o', label='trapping')
    plt.plot(np.log(drug_conc[1:]), peaks_conductance_norm[1:],
             'o', label='w/o trapping')
    plt.xlabel('Drug concentration (log)')
    plt.ylabel('Normalised peak current')
    plt.legend()
    plt.tight_layout()
    plt.savefig(saved_fig_dir + "norm_peak_compare_hERG_" + drug
                + "_concs.pdf")

# Check action potential
# Load model
APmodel, _, x = myokit.load(APmodel)

AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')
pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)
conductance = APmodel.get('ikr.gKr').value()

offset = 50
if steady_state:
    if drug == 'verapamil' and protocol_name == 'P0':
        repeats = 83
    else:
        repeats = 150
    save_signal = 2
else:
    repeats = 10
    save_signal = repeats

AP_conductance = []
APD_conductance = []
hERG_peak_conductance = []

AP_trapping = []
APD_trapping = []
hERG_peak_trapping = []

for i in range(len(drug_conc)):
    log = AP_model.drug_simulation(drug, drug_conc[i], repeats,  # + 1
                                   timestep=0.1, save_signal=save_signal)
    APD_trapping_pulse = []
    hERG_trapping_pulse = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(log['membrane.V', pulse], offset, 0.1)
        APD_trapping_pulse.append(apd90)

        hERG_peak = np.max(log['ikr.IKr', pulse])
        hERG_trapping_pulse.append(hERG_peak)

    AP_trapping.append(log)
    APD_trapping.append(APD_trapping_pulse)
    hERG_peak_trapping.append(hERG_trapping_pulse)

    scale = conductance_scale_df.iloc[i][drug]
    d2 = AP_model.conductance_simulation(conductance * scale, repeats,
                                         timestep=0.1, save_signal=save_signal)
    APD_conductance_pulse = []
    hERG_conductance_pulse = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(d2['membrane.V', pulse], offset, 0.1)
        APD_conductance_pulse.append(apd90)

        hERG_peak = np.max(d2['ikr.IKr', pulse])
        hERG_conductance_pulse.append(hERG_peak)

    AP_conductance.append(d2)
    APD_conductance.append(APD_conductance_pulse)
    hERG_peak_conductance.append(hERG_conductance_pulse)

if plot_fig:
    if steady_state:

        # Plot action potentials - steady state
        plotting_pulse_time = pulse_time * save_signal

        fig = modelling.figures.FigureStructure(figsize=(8, 4),
                                                gridspec=(3, 2),
                                                wspace=0.08)
        plot = modelling.figures.FigurePlot()
        cmap = matplotlib.cm.get_cmap('viridis')

        plot.add_continuous(fig.axs[0][0], AP_trapping[0], 'stimulus.i_stim')
        plot.add_multiple_continuous(fig.axs[1][0], AP_trapping, 'membrane.V',
                                     cmap=cmap, labels=drug_labels)
        plot.add_multiple_continuous(fig.axs[2][0], AP_trapping, 'ikr.IKr',
                                     cmap=cmap, labels=drug_labels)
        plot.add_continuous(fig.axs[0][1], AP_conductance[0],
                            'stimulus.i_stim')
        plot.add_multiple_continuous(fig.axs[1][1], AP_conductance,
                                     'membrane.V', cmap=cmap,
                                     labels=drug_labels)
        plot.add_multiple_continuous(fig.axs[2][1], AP_conductance, 'ikr.IKr',
                                     cmap=cmap, labels=drug_labels)
        fig.axs[0][0].set_title('with trapping', fontsize=8)
        fig.axs[0][1].set_title('without trapping', fontsize=8)

        unique = fig.legend_without_duplicate_labels(fig.axs[2][1])
        fig.axs[2][1].legend(*zip(*unique), loc='right',
                             bbox_to_anchor=(1.4, 0.6))
        fig.sharex(['Time (ms)'] * 2, [(0, plotting_pulse_time)] * 2)
        fig.sharey(['Current\nstimulus', 'Voltage (mV)', 'Current (A/F)'])
        fig.savefig(saved_fig_dir + "2AP_steady_trapping_nontrapping.pdf")

        # APD plot at steady state
        APD_trap_plot = []
        APD_conduct_plot = []
        APD_trap_plot_previous = []
        APD_conduct_plot_previous = []
        for i in range(len(APD_trapping)):
            APD_trap_plot.append(APD_trapping[i][-1])
            APD_conduct_plot.append(APD_conductance[i][-1])
            APD_trap_plot_previous.append(APD_trapping[i][-2])
            APD_conduct_plot_previous.append(APD_conductance[i][-2])

        # Plot APD90 for 2 pulses
        plt.figure(figsize=(4, 3))
        plt.plot(np.log(drug_conc[1:]), APD_trap_plot[1:],
                 'o-', color='orange', label='trapping - last pulse')
        plt.plot(np.log(drug_conc[1:]), APD_trap_plot_previous[1:],
                 '^-', color='orange', label='trapping - 2nd last pulse')
        plt.plot(np.log(drug_conc[1:]), APD_conduct_plot[1:],
                 'o--', color='blue', label='w/o trapping - last pulse')
        plt.plot(np.log(drug_conc[1:]), APD_conduct_plot_previous[1:],
                 '^--', color='blue', label='w/o trapping - 2nd last pulse')
        plt.xlabel('Drug concentration (log)')
        plt.ylabel(r'APD$_{90}$')
        plt.legend()
        plt.tight_layout()
        plt.savefig(saved_fig_dir + "APD90_compare_2pulses_hERG_" + drug
                    + "_concs.pdf")
        plt.close()

    else:
        # Plot hERG peak - transient phase
        cmap = matplotlib.cm.get_cmap('tab10')
        norm = matplotlib.colors.Normalize(0, len(drug_conc))

    #     hERG_reduc_trap = []
    #     hERG_reduc_conduct = []

        fig, ax = plt.subplots(1, 1, figsize=(4, 3))
        for i in range(len(drug_conc)):
            ax.plot(np.arange(save_signal), np.array(hERG_peak_trapping[i]),
                    'o-', label=str(drug_conc[i]) + 'nM - trapping',
                    color=cmap(norm(i)))
            ax.plot(np.arange(save_signal), np.array(hERG_peak_conductance[i]),
                    '^--', label=str(drug_conc[i]) + 'nM - w/o trapping',
                    color=cmap(norm(i)))

            # # peak reduction
            # hERG_peak_conc_trap = hERG_peak_trapping[i]
            # hERG_reduc_trap.append((hERG_peak_conc_trap[0]
            #                         - hERG_peak_conc_trap[-1])
            #                         / hERG_peak_conc_trap[0])

            # hERG_peak_conc_conduct = hERG_peak_conductance[i]
            # hERG_reduc_conduct.append((hERG_peak_conc_conduct[0]
            #                            - hERG_peak_conc_conduct[-1])
            #                            / hERG_peak_conc_conduct[0])

        ax.legend(loc='upper right', bbox_to_anchor=(1.75, 1.5))
        plt.xlabel('Sweeps')
        plt.ylabel(r'APD$_{90}$')
        plt.savefig(saved_fig_dir + "peak_compare_hERG_" + drug + "_concs.pdf",
                    bbox_inches='tight')
        plt.close()

    #     plt.figure()
    #     plt.plot(hERG_reduc_trap)
    #     plt.plot(hERG_reduc_conduct)
    #     plt.savefig(saved_fig_dir + 'peak_reduction.pdf')

        # Plot action potentials - transient phase
        labels = [str(i) + ' nM' for i in drug_conc]
        cmap = matplotlib.cm.get_cmap('viridis')
        plotting_pulse_time = pulse_time * save_signal

        fig = modelling.figures.FigureStructure(figsize=(10, 4),
                                                gridspec=(2, 1),
                                                height_ratios=[1, 1])
        plot = modelling.figures.FigurePlot()

        plot.add_multiple_continuous(fig.axs[0][0], AP_trapping, 'membrane.V',
                                     cmap=cmap)
        plot.add_multiple_continuous(fig.axs[1][0], AP_trapping, 'ikr.IKr',
                                     labels=labels, cmap=cmap)

        fig.legend_without_duplicate_labels(fig.axs[1][0])
        fig.sharex(['Time (s)'], [(0, plotting_pulse_time)])
        fig.sharey(['Voltage (mV)', 'hERG current'])
        fig.adjust_ticks(fig.axs[1][0], plotting_pulse_time)

        fig.savefig(saved_fig_dir + "10AP_transient_hERG_trapping_" + drug
                    + "_concs.pdf")
        plt.close()

        fig = modelling.figures.FigureStructure(figsize=(10, 4),
                                                gridspec=(2, 1),
                                                height_ratios=[1, 1])
        plot = modelling.figures.FigurePlot()

        plot.add_multiple_continuous(fig.axs[0][0], AP_conductance,
                                     'membrane.V', cmap=cmap)
        plot.add_multiple_continuous(fig.axs[1][0], AP_conductance, 'ikr.IKr',
                                     labels=labels, cmap=cmap)

        fig.legend_without_duplicate_labels(fig.axs[1][0])
        fig.sharex(['Time (s)'], [(0, plotting_pulse_time)])
        fig.sharey(['Voltage (mV)', 'hERG current'])
        fig.adjust_ticks(fig.axs[1][0], plotting_pulse_time)

        fig.savefig(saved_fig_dir + "10AP_transient_hERG_conductance_" + drug
                    + "_concs.pdf")
        plt.close('all')

        # Plot phase plane
#         norm = matplotlib.colors.Normalize(0, repeats)
#         for i in range(len(AP_conductance)):
#             fig, ax = plt.subplots(1, 1, figsize=(4, 3))
#             for pulse in range(repeats):
#                 ax.scatter(AP_conductance[i]['membrane.V', pulse],
#                            AP_conductance[i]['ikr.IKr', pulse],
#                            color=cmap(norm(pulse)), zorder=-10)
#             ax.set_rasterization_zorder(0)
#             plt.savefig(saved_fig_dir + "phase_plane_" + str(drug_conc[i])
#                         + ".pdf")

        # Plot APD90 of 10 pulses
        cmap_trap = matplotlib.cm.get_cmap('Blues')
        cmap_conduct = matplotlib.cm.get_cmap('YlOrBr')
        cmap = matplotlib.cm.get_cmap('tab10')
        norm = matplotlib.colors.Normalize(0, len(drug_conc))

        fig, ax = plt.subplots(1, 1, figsize=(4, 3))
        for i in range(len(drug_conc)):
            ax.plot(np.arange(save_signal), np.array(APD_trapping[i]), 'o-',
                    label=str(drug_conc[i]) + 'nM - trapping',
                    color=cmap(norm(i)))
            ax.plot(np.arange(save_signal), np.array(APD_conductance[i]),
                    '^--', label=str(drug_conc[i]) + 'nM - w/o trapping',
                    color=cmap(norm(i)))

        ax.legend(loc='upper right', bbox_to_anchor=(1.75, 1.5))
        plt.xlabel('Sweeps')
        plt.ylabel(r'APD$_{90}$')
        plt.savefig(saved_fig_dir + "10APD90_compare_hERG_" + drug
                    + "_concs.pdf", bbox_inches='tight')
        plt.close()

        # if drug == 'verapamil':
        #     fig = modelling.figures.FigureStructure(figsize=(8, 6),
        #                                             gridspec=(3, 3),
        #                                             height_ratios=[1, 1, 1],
        #                                             hspace=0.2, wspace=0.08)
        #     for i in range(len(drug_conc)):
        #         for j in range(len(drug_conc)):
        #             if i == j:
        #                 fig.axs[i // 3][i % 3].plot(
        #                     np.arange(save_signal), np.array(APD_trapping[j]),
        #                     'o-', label='with trapping',
        #                     color='orange', zorder=5)
        #                 fig.axs[i // 3][i % 3].plot(
        #                     np.arange(save_signal),
        #                     np.array(APD_conductance[j]), '^--',
        #                     label='w/o trapping', color='blue', zorder=5)
        #             else:
        #                 fig.axs[i // 3][i % 3].plot(
        #                     np.arange(save_signal), np.array(APD_trapping[j]),
        #                     'o-', label=str(drug_conc[i]) + 'nM - trapping',
        #                     color='grey', alpha=0.2, zorder=1)
        #                 fig.axs[i // 3][i % 3].plot(
        #                     np.arange(save_signal),
        #                     np.array(APD_conductance[j]), '^--',
        #                     label=str(drug_conc[i]) + 'nM - w/o trapping',
        #                     color='grey', alpha=0.2, zorder=1)
        #         fig.axs[i // 3][i % 3].set_title(str(drug_conc[i]) + 'nM')

        #     handles, labels = fig.axs[0][0].get_legend_handles_labels()
        #     lgds = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if
        #             i <= 1]
        #     fig.axs[0][2].legend(*zip(*lgds), loc='right',
        #                          bbox_to_anchor=(1, 1.4))
        #     fig.sharex(['Sweeps'] * 3)
        #     fig.sharey([r"APD$_{90}$"] * 3)

        # if drug == 'dofetilide':
        fig = modelling.figures.FigureStructure(figsize=(10, 4),
                                                gridspec=(2, 4),
                                                height_ratios=[1, 1],
                                                hspace=0.2, wspace=0.08)
        for i in range(len(drug_conc)):
            for j in range(len(drug_conc)):
                if i == j:
                    fig.axs[int(i / 4)][i % 4].plot(
                        np.arange(save_signal), np.array(APD_trapping[j]),
                        'o-', label='with trapping',
                        color='orange', zorder=5)
                    fig.axs[int(i / 4)][i % 4].plot(
                        np.arange(save_signal),
                        np.array(APD_conductance[j]), '^--',
                        label='w/o trapping', color='blue', zorder=5)
                else:
                    fig.axs[int(i / 4)][i % 4].plot(
                        np.arange(save_signal), np.array(APD_trapping[j]),
                        'o-', label=str(drug_conc[i]) + 'nM - trapping',
                        color='grey', alpha=0.2, zorder=1)
                    fig.axs[int(i / 4)][i % 4].plot(
                        np.arange(save_signal),
                        np.array(APD_conductance[j]), '^--',
                        label=str(drug_conc[i]) + 'nM - w/o trapping',
                        color='grey', alpha=0.2, zorder=1)
            fig.axs[int(i / 4)][i % 4].set_title(str(drug_conc[i]) + 'nM')

        handles, labels = fig.axs[0][0].get_legend_handles_labels()
        lgds = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if
                i <= 1]
        fig.axs[0][3].legend(*zip(*lgds), loc='right',
                             bbox_to_anchor=(1, 1.4))
        fig.sharex(['Sweeps'] * 4)
        fig.sharey([r"APD$_{90}$"] * 2)

        fig.savefig(saved_fig_dir + "APD90_compare_transient_multiples.pdf")

        # Plot APD90 of 10 pulses - adjustments for specific case (verapamil,
        # Milnes)
        # fig, (ax, ax2) = plt.subplots(2, 1, figsize=(4, 3), sharex=True,
        #         gridspec_kw={'height_ratios': [1, 3]})
        # for i in range(len(drug_conc)):
        #     # mask = ~np.isnan(APD_trapping[i])
        #     ax2.plot(np.arange(save_signal), np.array(APD_trapping[i]), 'o-',
        #             label=str(int(drug_conc[i])) + 'nM - trapping',
        #             # color=cmap_trap(norm(i)))
        #             color=cmap(norm(i)))
        #     ax.plot(np.arange(save_signal), np.array(APD_trapping[i]), 'o-',
        #             color=cmap(norm(i)))
        #     ax2.plot(np.arange(save_signal), np.array(APD_conductance[i]),
        #              '^--',
        #              label=str(int(drug_conc[i])) + 'nM - w/o trapping',
        #             # color=cmap_conduct(norm(i)))
        #             color=cmap(norm(i)))
        #     ax.plot(np.arange(save_signal), np.array(APD_conductance[i]),
        #             '^--', color=cmap(norm(i)))
        #
        # ax.set_ylim(950, 1000)
        # ax2.set_ylim(200, 650)
        #
        # ax.spines['bottom'].set_visible(False)
        # ax2.spines['top'].set_visible(False)
        # ax.xaxis.tick_top()
        # ax.tick_params(labeltop=False)
        # ax2.xaxis.tick_bottom()
        #
        # ax2.legend(loc='upper right', bbox_to_anchor=(1.75, 1.5))
        # plt.xlabel('Sweeps')
        # plt.ylabel(r'APD$_{90}$')
        # plt.savefig(saved_fig_dir + "10APD90_compare_hERG_" + drug
        #             + "_concs.pdf", bbox_inches='tight')
        # plt.close()
