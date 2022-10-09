# Compares the hERG current and action potential of state dependen drug block
# and conductance scaling
import matplotlib
import myokit
import numpy as np
import os
import pandas as pd

import modelling

# Set up directory for figures
testing_fig_dir = '../../figures/testing/'
final_fig_dir = '../../figures/binding_kinetics_comparison/'
# final_fig_dir = '../../figures/conferences/'

saved_fig_dir = final_fig_dir

# Set up figure's main grid
fig = modelling.figures.FigureStructure(
    # specs for conference poster
    # figsize=(10, 9),
    # gridspec=(3, 2), hspace=0.6,
    # wspace=0.55,
    figsize=(10, 8),
    gridspec=(3, 2), hspace=0.5,
    wspace=0.45,
    height_ratios=[1, 1, 1],
    plot_in_subgrid=True)
plot = modelling.figures.FigurePlot()

subgridspecs = [(2, 1)] * 4 + [(1, 2)] * 2
subgs = []
for i in range(3 * 2):
    subgs.append(fig.gs[i].subgridspec(*subgridspecs[i], wspace=0.1,
                                       hspace=0.05))
axs = [[[fig.fig.add_subplot(subgs[k][i, j]) for j in range(
    subgridspecs[k][1])] for i in range(subgridspecs[k][0])] for
    k in range(len(subgs))]

# Load data
drug = 'dofetilide'
protocol_name = 'Milnes'
saved_data_dir = '../../simulation_data/binding_kinetics_comparison/' + \
    drug + '/' + protocol_name + '/'

trapping_data_files = [f for f in os.listdir(saved_data_dir) if
                       f.startswith('CiPA_AP_transient_')]
conductance_data_files = [f for f in os.listdir(saved_data_dir) if
                          f.startswith('conductance_AP_transient_')]
conc_label = [fname[26:-4] for fname in trapping_data_files]
drug_conc = [float(fname[26:-4]) for fname in trapping_data_files]

# Sort in increasing order of drug concentration
sort_ind = [i[0] for i in sorted(enumerate(drug_conc), key=lambda x:x[1])]
drug_conc = sorted(drug_conc)
trapping_data_files = [trapping_data_files[i] for i in sort_ind]
conductance_data_files = [conductance_data_files[i] for i in sort_ind]
conc_label = [conc_label[i] for i in sort_ind]

CiPA_AP_log = []
conductance_AP_log = []
for i in range(len(trapping_data_files)):
    CiPA_AP_log.append(myokit.DataLog.load_csv(
        saved_data_dir + trapping_data_files[i]))
    conductance_AP_log.append(myokit.DataLog.load_csv(
        saved_data_dir + conductance_data_files[i]))

# Initiate constants and variables
labels = [i + ' nM' for i in conc_label]
plotting_pulse_time = 1000 * 7
cmap = matplotlib.cm.get_cmap('viridis')

# Top left panel
panel1 = axs[0]
plot.add_multiple_continuous(panel1[0][0], CiPA_AP_log, 'membrane.V',
                             cmap=cmap)
plot.add_multiple_continuous(panel1[1][0], CiPA_AP_log, 'ikr.IKr',
                             labels=labels, cmap=cmap)

unique = fig.legend_without_duplicate_labels(panel1[1][0])
panel1[1][0].legend(*zip(*unique), loc='lower left',
                    bbox_to_anchor=(1.0, 0), handlelength=1)
fig.sharex(['Time (s)'], [(0, plotting_pulse_time)],
           axs=panel1, subgridspec=subgridspecs[0])
fig.sharey(['Voltage\n(mV)', 'hERG current'],
           axs=panel1, subgridspec=subgridspecs[0])
fig.adjust_ticks(panel1[1][0], plotting_pulse_time)
panel1[0][0].set_title(drug + "-like drug, state dependent drug block")

# Second row, left panel
panel3 = axs[2]
plot.add_multiple_continuous(panel3[0][0], conductance_AP_log, 'membrane.V',
                             cmap=cmap)
plot.add_multiple_continuous(panel3[1][0], conductance_AP_log, 'ikr.IKr',
                             labels=labels, cmap=cmap)

unique = fig.legend_without_duplicate_labels(panel3[1][0])
panel3[1][0].legend(*zip(*unique), loc='lower left',
                    bbox_to_anchor=(1.0, 0), handlelength=1)
fig.sharex(['Time (s)'], [(0, plotting_pulse_time)],
           axs=panel3, subgridspec=subgridspecs[2])
fig.sharey(['Voltage\n(mV)', 'hERG current'],
           axs=panel3, subgridspec=subgridspecs[2])
fig.adjust_ticks(panel3[1][0], plotting_pulse_time)
panel3[0][0].set_title(drug + "-like drug, conductance scaling")

# Load APD values
APD_trapping = pd.read_csv(saved_data_dir + 'CiPA_APD_transient.csv')
APD_conductance = pd.read_csv(saved_data_dir + 'conductance_APD_transient.csv')

APD_trapping_concs = []
APD_conductance_concs = []
drug_conc = APD_trapping['drug concentration'].values.tolist()
for i in range(len(drug_conc)):
    APD_conc = APD_trapping.loc[
        APD_trapping['drug concentration'] == drug_conc[i]]
    APD_conc = APD_conc.drop(
        columns='drug concentration').values.tolist()[0][1:]
    saved_signal = len(APD_conc)
    APD_trapping_concs.append(APD_conc)

    APD_conc = APD_conductance.loc[
        APD_conductance['drug concentration'] == drug_conc[i]]
    APD_conc = APD_conc.drop(
        columns='drug concentration').values.tolist()[0][1:]
    APD_conductance_concs.append(APD_conc)

# Third row, left panel
panel5 = axs[4]
for i in range(len(drug_conc)):
    APD_plot = [APD_conductance_concs[i][ind] for ind in
                range(len(APD_conductance_concs[i])) if ind % 2 == 0]
    panel5[0][i].plot(np.arange(saved_signal / 2) * 2, APD_plot, 'o', ms=0.9,
                      label='conductance\nscaling', color='blue')

    APD_plot = [APD_trapping_concs[i][ind] for ind in
                range(len(APD_trapping_concs[i])) if ind % 2 == 0]
    panel5[0][i].plot(np.arange(saved_signal / 2) * 2, APD_plot, 'o', ms=0.9,
                      label='state dependent\ndrug block', color='orange')
    panel5[0][i].set_title(str(drug_conc[i]) + ' nM')

panel5[0][0].plot(0, 1050, 'o', color='k', marker=(5, 2))
panel5[0][1].plot(150, 1050, 'o', color='k', marker=(5, 2),
                  label='EAD-like\nbehaviour')
# panel5[0][1].legend(loc='lower left', bbox_to_anchor=(1.0, 0),
#                     handlelength=1)
fig.sharex(['Sweeps'] * 2,
           axs=panel5, subgridspec=subgridspecs[4])
fig.sharey([r"APD$_{90}$ (ms)"],
           axs=panel5, subgridspec=subgridspecs[4])

# Load data for verapamil
drug = 'verapamil'
protocol_name = 'Milnes'
saved_data_dir = '../../simulation_data/binding_kinetics_comparison/' + \
    drug + '/' + protocol_name + '/'

trapping_data_files = [f for f in os.listdir(saved_data_dir) if
                       f.startswith('CiPA_AP_transient_')]
conductance_data_files = [f for f in os.listdir(saved_data_dir) if
                          f.startswith('conductance_AP_transient_')]
conc_label = [fname[26:-4] for fname in trapping_data_files]
drug_conc = [float(fname[26:-4]) for fname in trapping_data_files]

# Sort in increasing order of drug concentration
sort_ind = [i[0] for i in sorted(enumerate(drug_conc), key=lambda x:x[1])]
drug_conc = sorted(drug_conc)
trapping_data_files = [trapping_data_files[i] for i in sort_ind]
conductance_data_files = [conductance_data_files[i] for i in sort_ind]
conc_label = [conc_label[i] for i in sort_ind]

CiPA_AP_log = []
conductance_AP_log = []
for i in range(len(trapping_data_files)):
    CiPA_AP_log.append(myokit.DataLog.load_csv(
        saved_data_dir + trapping_data_files[i]))
    conductance_AP_log.append(myokit.DataLog.load_csv(
        saved_data_dir + conductance_data_files[i]))

# Initiate constants and variables
labels = [i + ' nM' for i in conc_label]
plotting_pulse_time = 1000 * 7
cmap = matplotlib.cm.get_cmap('viridis')

# Top right panel
panel2 = axs[1]
plot.add_multiple_continuous(panel2[0][0], CiPA_AP_log, 'membrane.V',
                             cmap=cmap)
plot.add_multiple_continuous(panel2[1][0], CiPA_AP_log, 'ikr.IKr',
                             labels=labels, cmap=cmap)

unique = fig.legend_without_duplicate_labels(panel2[1][0])
panel2[1][0].legend(*zip(*unique), loc='lower left',
                    bbox_to_anchor=(1.0, 0), handlelength=1)
fig.sharex(['Time (s)'], [(0, plotting_pulse_time)],
           axs=panel2, subgridspec=subgridspecs[1])
fig.sharey(['Voltage\n(mV)', 'hERG current'],
           axs=panel2, subgridspec=subgridspecs[1])
fig.adjust_ticks(panel2[1][0], plotting_pulse_time)
panel2[0][0].set_title(drug + "-like drug, state dependent drug block")

# Second row, right panel
panel4 = axs[3]
plot.add_multiple_continuous(panel4[0][0], conductance_AP_log, 'membrane.V',
                             cmap=cmap)
plot.add_multiple_continuous(panel4[1][0], conductance_AP_log, 'ikr.IKr',
                             labels=labels, cmap=cmap)

unique = fig.legend_without_duplicate_labels(panel4[1][0])
panel4[1][0].legend(*zip(*unique), loc='lower left',
                    bbox_to_anchor=(1.0, 0), handlelength=1)
fig.sharex(['Time (s)'], [(0, plotting_pulse_time)],
           axs=panel4, subgridspec=subgridspecs[3])
fig.sharey(['Voltage\n(mV)', 'hERG current'],
           axs=panel4, subgridspec=subgridspecs[3])
fig.adjust_ticks(panel4[1][0], plotting_pulse_time)
panel4[0][0].set_title(drug + "-like drug, conductance scaling")

# Load APD values
APD_trapping = pd.read_csv(saved_data_dir + 'CiPA_APD_transient.csv')
APD_conductance = pd.read_csv(saved_data_dir + 'conductance_APD_transient.csv')

APD_trapping_concs = []
APD_conductance_concs = []
drug_conc = APD_trapping['drug concentration'].values.tolist()
for i in range(len(drug_conc)):
    APD_conc = APD_trapping.loc[
        APD_trapping['drug concentration'] == drug_conc[i]]
    APD_conc = APD_conc.drop(
        columns='drug concentration').values.tolist()[0][1:]
    saved_signal = len(APD_conc)
    APD_trapping_concs.append(APD_conc)

    APD_conc = APD_conductance.loc[
        APD_conductance['drug concentration'] == drug_conc[i]]
    APD_conc = APD_conc.drop(
        columns='drug concentration').values.tolist()[0][1:]
    APD_conductance_concs.append(APD_conc)

# Third row, right panel
panel6 = axs[5]
for i in range(len(drug_conc)):
    APD_plot = [APD_trapping_concs[i][ind] for ind in
                range(len(APD_trapping_concs[i])) if ind % 2 == 0]
    panel6[0][i].plot(np.arange(saved_signal / 2) * 2, APD_plot, 'o', ms=0.9,
                      label='state dependent\ndrug block', color='orange')
    panel6[0][i].plot(np.arange(saved_signal), APD_conductance_concs[i], 'o',
                      ms=0.9, label='conductance\nscaling', color='blue')
    panel6[0][i].set_title(str(drug_conc[i]) + ' nM')

panel6[0][1].plot(150, 1050, 'o', color='k', marker=(5, 2),
                  label='EAD-like\nbehaviour')
lgnd = panel6[0][1].legend(loc='lower left', bbox_to_anchor=(1.0, 0),
                           handlelength=1)
for handle in lgnd.legendHandles:
    handle.set_markersize(6)
fig.sharex(['Sweeps'] * 2,
           axs=panel6, subgridspec=subgridspecs[5])
fig.sharey([r"APD$_{90}$ (ms)"],
           axs=panel6, subgridspec=subgridspecs[5])

fig.fig.set_size_inches(10, 8)
fig.fig.text(0.075, 0.905, '(a)', fontsize=11)
fig.fig.text(0.535, 0.905, '(b)', fontsize=11)
fig.fig.text(0.075, 0.625, '(c)', fontsize=11)
fig.fig.text(0.535, 0.625, '(d)', fontsize=11)
fig.fig.text(0.075, 0.325, '(e)', fontsize=11)
fig.fig.text(0.535, 0.325, '(f)', fontsize=11)

fig.savefig(saved_fig_dir + "AP_transient.svg", format='svg')
