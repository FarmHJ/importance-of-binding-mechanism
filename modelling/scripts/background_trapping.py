# Introduces the idea of trapping and justifies the use of the Milnes protocol
import myokit
import numpy as np

import modelling

steady_state = False
hERG_model = False

drugs = ['dofetilide', 'verapamil']
drug_concs = [30, 1000]  # nM
# drug_label = ['drug free'] + [str(drug_concs[i]) + ' nM ' + drugs[i]
#                               for i in range(len(drugs))]
drug_label = ['drug free', 'trapped drug', 'nontrapped drug']

testing_fig_dir = '../../figures/testing/'
final_fig_dir = '../../figures/background/trapping/'
saved_fig_dir = testing_fig_dir

if hERG_model:
    model = '../../model/ohara-cipa-v1-2017-IKr.mmt'
    model, _, x = myokit.load(model)

    pulse_time = 25e3
    protocol = modelling.ProtocolLibrary().Milnes(pulse_time)

    drug_model = modelling.BindingKinetics(model)
    drug_model.protocol = protocol

if not hERG_model:
    APmodel = '../../model/ohara-cipa-v1-2017.mmt'
    APmodel, _, x = myokit.load(APmodel)

    pulse_time = 1000
    protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)

    AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')
    AP_model.protocol = protocol

if steady_state and hERG_model:
    repeats = 1000

    fig = modelling.figures.FigureStructure(figsize=(8, 4), gridspec=(3, 3),
                                            wspace=0.05)
    plot = modelling.figures.FigurePlot()

    log_all = []
    log_control = drug_model.drug_simulation(drugs[0], 0, repeats)
    log_all.append(log_control)

    plot.add_single(fig.axs[0][0], log_control, 'membrane.V')
    plot.state_occupancy_plot(fig.axs[2][0], log_control, model)

    for d in range(len(drugs)):
        log = drug_model.drug_simulation(drugs[d], drug_concs[d], repeats)
        log_all.append(log)

        plot.add_single(fig.axs[0][d + 1], log, 'membrane.V')
        plot.state_occupancy_plot(fig.axs[2][d + 1], log, model, legend=False)

    for i in range(len(log_all)):
        for j in range(len(log_all)):
            if i == j:
                plot.add_single(fig.axs[1][i], log_all[j], 'ikr.IKr')
            else:
                plot.add_single(fig.axs[1][i], log_all[j], 'ikr.IKr',
                                color='grey', alpha=0.5)
        fig.axs[1][i].text(24500, 0.65, drug_label[i], fontsize=8, ha='right')

    fig.sharex(['Time (s)'] * (len(drugs) + 1),
               [(0, pulse_time)] * (len(drugs) + 1))
    fig.sharey(['Voltage (mV)', 'Current (A/F)', 'State occupancy'])
    for i in range(len(log_all)):
        fig.adjust_ticks(fig.axs[2][i], pulse_time)
    fig.savefig(saved_fig_dir + "hERG_drug_state_occupancy.pdf")

elif steady_state and not hERG_model:
    repeats = 1000

    fig = modelling.figures.FigureStructure(figsize=(8, 6), gridspec=(4, 3),
                                            wspace=0.08)
    plot = modelling.figures.FigurePlot()

    log_all = []
    log_control = AP_model.drug_simulation(drugs[0], 0, repeats, timestep=0.1)
    log_all.append(log_control)

    plot.add_single(fig.axs[0][0], log_control, 'stimulus.i_stim')
    plot.state_occupancy_plot(fig.axs[3][0], log_control, APmodel)

    for d in range(len(drugs)):
        log = AP_model.drug_simulation(drugs[d], drug_concs[d], repeats)
        log_all.append(log)

        plot.add_single(fig.axs[0][d + 1], log, 'stimulus.i_stim')
        plot.state_occupancy_plot(fig.axs[3][d + 1], log,
                                  APmodel, legend=False)

    for i in range(len(log_all)):
        for j in range(len(log_all)):
            if i == j:
                plot.add_single(fig.axs[1][i], log_all[j], 'membrane.V')
                plot.add_single(fig.axs[2][i], log_all[j], 'ikr.IKr')
            else:
                plot.add_single(fig.axs[1][i], log_all[j], 'membrane.V',
                                color='grey', alpha=0.5)
                plot.add_single(fig.axs[2][i], log_all[j], 'ikr.IKr',
                                color='grey', alpha=0.5)
        fig.axs[1][i].text(900, 40, drug_label[i], fontsize=8, ha='right')

    fig.sharex(['Time (ms)'] * (len(drugs) + 1),
               [(0, pulse_time)] * (len(drugs) + 1))
    fig.sharey(['Current\nstimulus', 'Voltage (mV)',
                'Current (A/F)', 'State occupancy'])
    fig.savefig(saved_fig_dir + "AP_drug_state_occupancy.pdf")

else:
    repeats = 10

    fig = modelling.figures.FigureStructure(figsize=(8, 4), gridspec=(3, 1))
    plot = modelling.figures.FigurePlot()

    for d in range(len(drugs)):
        log = AP_model.drug_simulation(drugs[d], drug_concs[d], repeats,
                                       save_signal=repeats)
        if d == 0:
            max_hERG = np.max(log['ikr.IKr', 0])
        plot.add_continuous(fig.axs[d + 1][0], log, 'ikr.IKr')

    plot.add_continuous(fig.axs[0][0], log, 'membrane.V')
    fig.sharex(['Time (s)'], [(0, pulse_time * repeats)])
    fig.sharey(['Voltage (mV)', 'Current (A/F)\ntrapped drug',
                'Current (A/F)\nnontrapped drug'],
               [None, (0, max_hERG + 0.05), (0, max_hERG + 0.05)])
    fig.adjust_ticks(fig.axs[2][0], pulse_time * repeats)
    fig.savefig(saved_fig_dir + "hERG_trapped_nontrapped.pdf")
