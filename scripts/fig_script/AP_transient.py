# Compares the IKr and AP of the ORd-SD model and the ORd-CS model at first
# few pulses after addition of drugs
import matplotlib
import myokit
import numpy as np
import os
import pandas as pd

import modelling

# Set up directory for figures
fig_dir = '../../figures/model_comparison/'
root_dir = '../../simulation_data/model_comparison/'

# Set up figure's main grid
fig = modelling.figures.FigureStructure(
    figsize=(10, 8), gridspec=(3, 2), hspace=0.55,
    wspace=0.25, height_ratios=[1, 1, 1],
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

control_log = myokit.DataLog.load_csv(root_dir + 'steady_state_control.csv')

# Read file names of action potential data for both AP models with
# dofetilide-like drug
drug = 'dofetilide'
protocol_name = 'Milnes'
data_dir = root_dir + drug + '/' + protocol_name + '/'

SD_fileprefix = 'SD_AP_transient_pulses7_conc'
CS_fileprefix = 'CS_AP_transient_pulses7_conc'
trapping_data_files = [f for f in os.listdir(data_dir) if
                       f.startswith(SD_fileprefix)]
conductance_data_files = [f for f in os.listdir(data_dir) if
                          f.startswith(CS_fileprefix)]
conc_label = [fname[len(SD_fileprefix):-10] for fname in trapping_data_files]
drug_conc = [float(fname[len(SD_fileprefix):-10])
             for fname in trapping_data_files]

# Sort drug concentrations in increasing order
sort_ind = [i[0] for i in sorted(enumerate(drug_conc), key=lambda x:x[1])]
drug_conc = sorted(drug_conc)
trapping_data_files = [trapping_data_files[i] for i in sort_ind]
conductance_data_files = [conductance_data_files[i] for i in sort_ind]
conc_label = [conc_label[i] for i in sort_ind]

# Load action potential data
trapping_AP_log = []
conductance_AP_log = []
for i in range(len(trapping_data_files)):
    trapping_AP_log.append(myokit.DataLog.load_csv(
        data_dir + trapping_data_files[i]))
    conductance_AP_log.append(myokit.DataLog.load_csv(
        data_dir + conductance_data_files[i]))

# Initiate constants and variables
labels = [i + ' nM' for i in conc_label]
plotting_pulse_time = 1000 * 7
cmap = matplotlib.colors.ListedColormap(
    matplotlib.cm.get_cmap('tab10')(range(3)))

# Top left panel
# Plot the action potential and IKr of the ORd-SD model with various
# concentrations of a dofetilide-like drug
panel1 = axs[0]
# Plot AP and IKr at steady state of control condition
plot.add_single(panel1[0][0], control_log, 'membrane.V', color=cmap(0))
plot.add_single(panel1[1][0], control_log, 'ikr.IKr', color=cmap(0))
# Indicate time point of drug addition
panel1[0][0].axvline(999, color='k', linestyle='--')
panel1[1][0].axvline(999, color='k', linestyle='--')

# Plot AP and IKr after addition of drugs
plot.add_multiple_continuous(panel1[0][0], trapping_AP_log, 'membrane.V',
                             cmap=cmap, starting_pos=1)
plot.add_multiple_continuous(panel1[1][0], trapping_AP_log, 'ikr.IKr',
                             labels=labels, cmap=cmap, starting_pos=1)

# Adjust figure details
unique = fig.legend_without_duplicate_labels(panel1[1][0])
panel1[1][0].legend(*zip(*unique), loc='lower left',
                    bbox_to_anchor=(0, 2.0), handlelength=1, ncol=3,
                    columnspacing=1)
fig.sharex(['Time (s)'], [(0, plotting_pulse_time)],
           axs=panel1, subgridspec=subgridspecs[0])
fig.sharey(['Voltage\n(mV)', 'Current (A/F)'],
           axs=panel1, subgridspec=subgridspecs[0])
fig.adjust_ticks(panel1[1][0], plotting_pulse_time)
panel1[0][0].set_title(drug + "-like drug, ORd-SD model", y=1.25)

# Second row, left panel
# Plot the action potential and IKr of the ORd-CS model with various
# concentrations of a dofetilide-like drug
panel3 = axs[2]
# Plot AP and IKr at steady state of control condition
plot.add_single(panel3[0][0], control_log, 'membrane.V', color=cmap(0))
plot.add_single(panel3[1][0], control_log, 'ikr.IKr', color=cmap(0))
# Indicate location of drug addition
panel3[0][0].axvline(999, color='k', linestyle='--')
panel3[1][0].axvline(999, color='k', linestyle='--')

# Plot AP and IKr after addition of drugs
plot.add_multiple_continuous(panel3[0][0], conductance_AP_log, 'membrane.V',
                             cmap=cmap, starting_pos=1)
plot.add_multiple_continuous(panel3[1][0], conductance_AP_log, 'ikr.IKr',
                             labels=labels, cmap=cmap, starting_pos=1)

# Adjust figure details
unique = fig.legend_without_duplicate_labels(panel3[1][0])
panel3[1][0].legend(*zip(*unique), loc='lower left',
                    bbox_to_anchor=(0, 2.0), handlelength=1, ncol=3,
                    columnspacing=1)
fig.sharex(['Time (s)'], [(0, plotting_pulse_time)],
           axs=panel3, subgridspec=subgridspecs[2])
fig.sharey(['Voltage\n(mV)', 'Current (A/F)'],
           axs=panel3, subgridspec=subgridspecs[2])
fig.adjust_ticks(panel3[1][0], plotting_pulse_time)
panel3[0][0].set_title(drug + "-like drug, ORd-CS model", y=1.25)

# Third row, left panel
# Plot the APD90 of the AP models with dofetilide-like drug for 300 pulses
panel5 = axs[4]

# Load APD values of the first 300 pulses for dofetilide-like drug
APD_trapping = pd.read_csv(data_dir + 'SD_APD_transient_paced.csv')
APD_conductance = pd.read_csv(data_dir + 'CS_APD_transient_paced.csv')

# Organise data for plotting
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

# Plot APD of both models for the first 300 pulses (dofetilide-like drug)
for i in range(len(drug_conc)):
    APD_plot = [APD_conductance_concs[i][ind] for ind in
                range(len(APD_conductance_concs[i])) if ind % 2 == 0]
    panel5[0][i].plot(np.arange(saved_signal / 2) * 2, APD_plot, 'o', ms=0.9,
                      label='ORd-CS model', color='blue')

    APD_plot = [APD_trapping_concs[i][ind] for ind in
                range(len(APD_trapping_concs[i])) if ind % 2 == 1]
    panel5[0][i].plot(np.arange(saved_signal / 2) * 2, APD_plot, 'o', ms=0.9,
                      label='ORd-SD model', color='orange')
    panel5[0][i].set_title(str(drug_conc[i]) + ' nM')

# Adjust figure details
panel5[0][0].plot(0, 1050, 'o', marker=(5, 2), color='k')
panel5[0][1].plot(150, 1050, 'o', marker=(5, 2), color='k',
                  label='EAD-like AP')
handles, labels = panel5[0][1].get_legend_handles_labels()
lgd_order = [1, 0, 2]
lgnd = panel5[0][1].legend([handles[idx] for idx in lgd_order],
                           [labels[idx] for idx in lgd_order],
                           loc='lower right', bbox_to_anchor=(1.0, 0),
                           handlelength=1)
for handle in lgnd.legendHandles:
    handle.set_markersize(6)
fig.sharex(['Pulses'] * 2,
           axs=panel5, subgridspec=subgridspecs[4])
fig.sharey([r"APD$_{90}$ (ms)"],
           axs=panel5, subgridspec=subgridspecs[4])

# Read file names of AP data for both AP models with verapamil-like drug
drug = 'verapamil'
data_dir = root_dir + drug + '/' + protocol_name + '/'

trapping_data_files = [f for f in os.listdir(data_dir) if
                       f.startswith(SD_fileprefix)]
conductance_data_files = [f for f in os.listdir(data_dir) if
                          f.startswith(CS_fileprefix)]
conc_label = [fname[len(SD_fileprefix):-10] for fname in trapping_data_files]
drug_conc = [float(fname[len(SD_fileprefix):-10])
             for fname in trapping_data_files]

# Sort drug concentration in increasing order
sort_ind = [i[0] for i in sorted(enumerate(drug_conc), key=lambda x:x[1])]
drug_conc = sorted(drug_conc)
trapping_data_files = [trapping_data_files[i] for i in sort_ind]
conductance_data_files = [conductance_data_files[i] for i in sort_ind]
conc_label = [conc_label[i] for i in sort_ind]

# Load AP data
trapping_AP_log = []
conductance_AP_log = []
for i in range(len(trapping_data_files)):
    trapping_AP_log.append(myokit.DataLog.load_csv(
        data_dir + trapping_data_files[i]))
    conductance_AP_log.append(myokit.DataLog.load_csv(
        data_dir + conductance_data_files[i]))

# Initiate constants and variables
labels = [i + ' nM' for i in conc_label]
labels[-1] = r"$10^4$ nM"
plotting_pulse_time = 1000 * 7

# Top right panel
# Plot the action potential and IKr of the ORd-SD model with various
# concentrations of a verapamil-like drug
panel2 = axs[1]
# Plot AP and IKr at steady state of control condition
plot.add_single(panel2[0][0], control_log, 'membrane.V', color=cmap(0))
plot.add_single(panel2[1][0], control_log, 'ikr.IKr', color=cmap(0))
# Indicate time point of drug addition
panel2[0][0].axvline(999, color='k', linestyle='--')
panel2[1][0].axvline(999, color='k', linestyle='--')

# Plot AP and IKr after addition of drugs
plot.add_multiple_continuous(panel2[0][0], trapping_AP_log, 'membrane.V',
                             cmap=cmap, starting_pos=1)
plot.add_multiple_continuous(panel2[1][0], trapping_AP_log, 'ikr.IKr',
                             labels=labels, cmap=cmap, starting_pos=1)

# Adjust figure details
unique = fig.legend_without_duplicate_labels(panel2[1][0])
panel2[1][0].legend(*zip(*unique), loc='lower left',
                    bbox_to_anchor=(0, 2.0), handlelength=1, ncol=3,
                    columnspacing=1)
fig.sharex(['Time (s)'], [(0, plotting_pulse_time)],
           axs=panel2, subgridspec=subgridspecs[1])
fig.sharey(['Voltage\n(mV)', 'Current (A/F)'],
           axs=panel2, subgridspec=subgridspecs[1])
fig.adjust_ticks(panel2[1][0], plotting_pulse_time)
panel2[0][0].set_title(drug + "-like drug, ORd-SD model", y=1.25)

# Second row, right panel
# Plot the action potential and IKr of the ORd-CS model with various
# concentrations of a verapamikl-like drug
panel4 = axs[3]
# Plot AP and IKr at steady state of control condition
plot.add_single(panel4[0][0], control_log, 'membrane.V', color=cmap(0))
plot.add_single(panel4[1][0], control_log, 'ikr.IKr', color=cmap(0))
# Indicate location of drug addition
panel4[0][0].axvline(999, color='k', linestyle='--')
panel4[1][0].axvline(999, color='k', linestyle='--')

# Plot AP and IKr after addition of drugs
plot.add_multiple_continuous(panel4[0][0], conductance_AP_log, 'membrane.V',
                             cmap=cmap, starting_pos=1)
plot.add_multiple_continuous(panel4[1][0], conductance_AP_log, 'ikr.IKr',
                             labels=labels, cmap=cmap, starting_pos=1)

# Adjust figure details
unique = fig.legend_without_duplicate_labels(panel4[1][0])
panel4[1][0].legend(*zip(*unique), loc='lower left',
                    bbox_to_anchor=(0, 2.0), handlelength=1, ncol=3,
                    columnspacing=1)
fig.sharex(['Time (s)'], [(0, plotting_pulse_time)],
           axs=panel4, subgridspec=subgridspecs[3])
fig.sharey(['Voltage\n(mV)', 'Current (A/F)'],
           axs=panel4, subgridspec=subgridspecs[3])
fig.adjust_ticks(panel4[1][0], plotting_pulse_time)
panel4[0][0].set_title(drug + "-like drug, ORd-CS model", y=1.25)

# Third row, right panel
# Plot the APD90 of the AP models with verapamil-like drug for 300 pulses
panel6 = axs[5]

# Load APD values of the first 300 pulses for verapamil-like drug
APD_trapping = pd.read_csv(data_dir + 'SD_APD_transient_paced.csv')
APD_conductance = pd.read_csv(data_dir + 'CS_APD_transient_paced.csv')

# Organise data for plotting
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

# Plot APD of both models for the first 300 pulses (verapamil-like drug)
for i in range(len(drug_conc)):
    APD_plot = [APD_trapping_concs[i][ind] for ind in
                range(len(APD_trapping_concs[i])) if ind % 2 == 0]
    panel6[0][i].plot(np.arange(saved_signal / 2) * 2, APD_plot, 'o', ms=0.9,
                      label='ORd-SD model', color='orange')
    panel6[0][i].plot(np.arange(saved_signal), APD_conductance_concs[i], 'o',
                      ms=0.9, label='ORd-CS model', color='blue')
    panel6[0][i].set_title(str(drug_conc[i]) + ' nM')

# Adjust figure details
panel6[0][1].plot(150, 1050, 'o', color='k', marker=(5, 2),
                  label='EAD-like AP')
lgnd = panel6[0][1].legend(loc='upper left', bbox_to_anchor=(-1.1, 1.0),
                           handlelength=1)
for handle in lgnd.legendHandles:
    handle.set_markersize(6)
fig.sharex(['Pulses'] * 2,
           axs=panel6, subgridspec=subgridspecs[5])
fig.sharey([r"APD$_{90}$ (ms)"],
           axs=panel6, subgridspec=subgridspecs[5])

# Add panel labels
fig.fig.set_size_inches(10, 8)
fig.fig.text(0.1, 0.905, '(A)', fontsize=11)
fig.fig.text(0.535, 0.905, '(B)', fontsize=11)
fig.fig.text(0.1, 0.625, '(C)', fontsize=11)
fig.fig.text(0.535, 0.625, '(D)', fontsize=11)
fig.fig.text(0.1, 0.325, '(E)', fontsize=11)
fig.fig.text(0.535, 0.325, '(F)', fontsize=11)

# Save figure
fig.savefig(fig_dir + "AP_transient_paced.pdf")
