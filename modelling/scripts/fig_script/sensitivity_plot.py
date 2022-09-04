# Plots figures to see how APD90 changes with different drugs
import matplotlib
import myokit
import numpy as np
import os
import pandas as pd
import sys

import modelling

drug = sys.argv[1]
protocol_name = 'Milnes'

# Set up directory for figures
testing_fig_dir = '../../figures/testing/'
final_fig_dir = '../../figures/binding_kinetics_comparison/' + drug + '/' + \
    protocol_name + '/'

saved_fig_dir = final_fig_dir

fig = modelling.figures.FigureStructure(
    figsize=(10, 8),
    gridspec=(3, 2), hspace=0.5,
    wspace=0.3,
    height_ratios=[1] * 3,
    plot_in_subgrid=True)
plot = modelling.figures.FigurePlot()

subgridspecs = [(2, 1)] * 4 + [(1, 1)] * 2
height_ratios = [[1, 2]] * 2 + [None] * 4
subgs = []
for i in range(3 * 2):
    subgs.append(fig.gs[i].subgridspec(*subgridspecs[i], wspace=0.1,
                                       hspace=0.05,
                                       height_ratios=height_ratios[i]))
axs = [[[fig.fig.add_subplot(subgs[k][i, j]) for j in range(
    subgridspecs[k][1])] for i in range(subgridspecs[k][0])] for
    k in range(len(subgs))]


saved_data_dir = '../../simulation_data/binding_kinetics_comparison/' + \
    drug + '/' + protocol_name + '/'

# Top panel
panel1 = axs[0]
panel2 = axs[1]
# Load hERG current data
trapping_data_files = [f for f in os.listdir(saved_data_dir) if
                       f.startswith('CiPA_hERG_')]
conductance_data_files = [f for f in os.listdir(saved_data_dir) if
                          f.startswith('conductance_hERG_')]
conc_label = [fname[10:-4] for fname in trapping_data_files]
drug_conc = [float(fname[10:-4]) for fname in trapping_data_files]

# Sort in increasing order of drug concentration
sort_ind = [i[0] for i in sorted(enumerate(drug_conc), key=lambda x:x[1])]
drug_conc = sorted(drug_conc)
trapping_data_files = [trapping_data_files[i] for i in sort_ind]
conductance_data_files = [conductance_data_files[i] for i in sort_ind]
conc_label = [conc_label[i] for i in sort_ind]

CiPA_hERG_log = []
conductance_hERG_log = []
for i in range(len(trapping_data_files)):
    CiPA_hERG_log.append(myokit.DataLog.load_csv(
        saved_data_dir + trapping_data_files[i]))
    conductance_hERG_log.append(myokit.DataLog.load_csv(
        saved_data_dir + conductance_data_files[i]))

# Initiate constants and variables
labels = [i + ' nM' for i in conc_label]
pulse_time = 25e3
cmap = matplotlib.cm.get_cmap('viridis')

# Plot figure
plot.add_single(panel1[0][0], CiPA_hERG_log[0], 'membrane.V', color='k')
plot.add_multiple(panel1[1][0], CiPA_hERG_log, 'ikr.IKr', labels=labels,
                  color=cmap)
plot.add_single(panel2[0][0], conductance_hERG_log[0], 'membrane.V', color='k')
plot.add_multiple(panel2[1][0], conductance_hERG_log, 'ikr.IKr', labels=labels,
                  color=cmap)

# panel1[1][0].legend()
panel1[0][0].set_title('State dependent drug block')
panel2[0][0].set_title('Drug block by conductance scaling')
for i in range(2):
    fig.sharex(['Time (s)'], [(0, pulse_time)],
               axs=axs[i], subgridspec=subgridspecs[i])
    fig.sharey(['Voltage\n(mV)', 'Current (A/F)'],
               axs=axs[i], subgridspec=subgridspecs[i])
fig.adjust_ticks(panel1[1][0], pulse_time)
fig.adjust_ticks(panel2[1][0], pulse_time)

# Second row panels
panel3 = axs[2]
panel4 = axs[3]

trapping_data_files = [f for f in os.listdir(saved_data_dir) if
                       f.startswith('CiPA_AP_') and not
                       f.startswith('CiPA_AP_transient')]
conductance_data_files = [f for f in os.listdir(saved_data_dir) if
                          f.startswith('conductance_AP_') and not
                          f.startswith('conductance_AP_transient')]
conc_label = [fname[8:-4] for fname in trapping_data_files]
drug_conc = [float(fname[8:-4]) for fname in trapping_data_files]

# Sort in increasing order of drug concentration
sort_ind = [i[0] for i in sorted(enumerate(drug_conc), key=lambda x:x[1])]
drug_conc = sorted(drug_conc)
trapping_data_files = [trapping_data_files[i] for i in sort_ind]
conductance_data_files = [conductance_data_files[i] for i in sort_ind]
conc_label = [conc_label[i] for i in sort_ind]
labels = [i + ' nM' for i in conc_label]

CiPA_AP_log = []
conductance_AP_log = []
for i in range(len(trapping_data_files)):
    CiPA_AP_log.append(myokit.DataLog.load_csv(
        saved_data_dir + trapping_data_files[i]))
    conductance_AP_log.append(myokit.DataLog.load_csv(
        saved_data_dir + conductance_data_files[i]))

APD_trapping = pd.read_csv(saved_data_dir + 'CiPA_APD_pulses1000.csv')
APD_conductance = pd.read_csv(saved_data_dir + 'conductance_APD_pulses1000.csv')

APD_trapping = [max(APD_trapping.loc[i].values.tolist()[(1 + 998):-1]) for i in
                range(APD_trapping.shape[0])]
APD_conductance = [max(APD_conductance.loc[i].values.tolist()[(1 + 998):-1]) for i in
                   range(APD_conductance.shape[0])]

# Initiate constants and variables
plotting_pulse_time = 1000 * 2

# Remove repeated signals at high concentrations
second_EAD_trap = [i for i, e in enumerate(APD_trapping)
                   if e == 1000][1:]
second_EAD_conduct = [i for i, e in enumerate(APD_conductance)
                      if e == 1000][1:]
if len(second_EAD_trap) <= len(second_EAD_conduct):
    AP_trapping_plot = [e for i, e in enumerate(CiPA_AP_log)
                        if i not in second_EAD_trap]
    AP_conductance_plot = [e for i, e in enumerate(conductance_AP_log)
                           if i not in second_EAD_trap]
else:
    AP_trapping_plot = [e for i, e in enumerate(CiPA_AP_log)
                        if i not in second_EAD_conduct]
    AP_conductance_plot = [e for i, e in enumerate(conductance_AP_log)
                           if i not in second_EAD_conduct]

# Plot figure
plot.add_multiple_continuous(panel3[0][0], AP_trapping_plot,
                             'membrane.V', cmap=cmap,
                             labels=labels)
plot.add_multiple_continuous(panel3[1][0], AP_trapping_plot,
                             'ikr.IKr', cmap=cmap, labels=labels)
plot.add_multiple_continuous(panel4[0][0], AP_conductance_plot,
                             'membrane.V', cmap=cmap,
                             labels=labels)
plot.add_multiple_continuous(panel4[1][0], AP_conductance_plot,
                             'ikr.IKr', cmap=cmap, labels=labels)
panel3[0][0].set_title('State dependent drug block')
panel4[0][0].set_title('Drug block by conductance scaling')

unique = fig.legend_without_duplicate_labels(panel4[1][0])
panel4[1][0].legend(*zip(*unique), loc='lower left',
                    bbox_to_anchor=(1.0, 0))
for i in range(2):
    fig.sharex(['Time (s)'], [(0, plotting_pulse_time)],
               axs=axs[i + 2], subgridspec=subgridspecs[i + 2])
    fig.sharey(['Voltage\n(mV)', 'Current (A/F)'],
               axs=axs[i + 2], subgridspec=subgridspecs[i + 2])

# Bottom left panel
panel5 = axs[4]

APD_trapping = pd.read_csv(saved_data_dir + 'CiPA_APD_pulses1000.csv')
APD_conductance = pd.read_csv(saved_data_dir + 'conductance_APD_pulses1000.csv')

drug_conc = APD_trapping['drug concentration'].values.tolist()
APD_trapping = [max(APD_trapping.loc[i].values.tolist()[1 + 998:-1]) for i in
                range(APD_trapping.shape[0])]
APD_conductance = [max(APD_conductance.loc[i].values.tolist()[1 + 998:-1])
                   for i in range(APD_conductance.shape[0])]
EAD_marker = [1050 if (i >= 1000 or j >= 1000) else None for (i, j)
              in zip(APD_trapping[1:], APD_conductance[1:])]

####################
# Only for dofetilide
# drug_conc = drug_conc[:-5]
# APD_trapping = APD_trapping[:-5]
# APD_conductance = APD_conductance[:-5]
# EAD_marker = EAD_marker[:-5]

# Left panel
panel5[0][0].plot(drug_conc[1:], APD_trapping[1:],
                  'o', color='orange', label='state dependent drug block')
panel5[0][0].plot(drug_conc[1:], APD_conductance[1:],
                  '^', color='blue', label='conductance scaling', alpha=0.8)
panel5[0][0].plot(drug_conc[1:], EAD_marker, marker=(5, 2), color='k',
                  label='EAD-like AP')
panel5[0][0].set_xscale("log", nonpositive='clip')
panel5[0][0].set_xlabel('Drug concentration (nM)')
panel5[0][0].set_ylabel(r'APD$_{90}$ (ms)')
panel5[0][0].legend(handlelength=1)

# Bottom right panel
panel6 = axs[5]

Hill_model = modelling.HillsModel()
color = ['orange', 'blue', 'red', 'green']

Hill_coef_df = pd.read_csv(saved_data_dir + '../Hill_curves.csv')

drug_conc_lib = modelling.DrugConcentrations()
conc_grid = drug_conc_lib.drug_concentrations[drug]['fine']

protocol_list = modelling.ProtocolParameters().protocols

for p, prot in enumerate(protocol_list):
    Hill_eq = Hill_coef_df.loc[
        Hill_coef_df['protocol'] == prot]
    Hill_eq = Hill_eq.values.tolist()[0][1:-1]

    panel6[0][0].plot(conc_grid, Hill_model.simulate(
        Hill_eq, conc_grid), '-', color=color[p],
        label=prot)

panel6[0][0].set_xscale("log", nonpositive='clip')
panel6[0][0].set_xlabel('Drug concentration (nM)')
panel6[0][0].set_ylabel('Normalised\npeak current')
panel6[0][0].legend()

fig.savefig(saved_fig_dir + 'sensitivity_summary.pdf')

# Figure 2
APD_trapping = pd.read_csv(saved_data_dir + 'CiPA_APD_pulses1000.csv')
APD_conductance = pd.read_csv(
    saved_data_dir + 'conductance_APD_pulses1000.csv')

APD_trapping_concs = []
APD_conductance_concs = []
drug_conc = APD_trapping['drug concentration'].values.tolist()
drug_conc = [i for i in drug_conc if i != 0]

for _, conc in enumerate(drug_conc):
    APD_conc = APD_trapping.loc[
        APD_trapping['drug concentration'] == conc]
    APD_conc = APD_conc.drop(
        columns='drug concentration').values.tolist()[0][1:]
    saved_signal = len(APD_conc)
    APD_trapping_concs.append(APD_conc)

    APD_conc = APD_conductance.loc[
        APD_conductance['drug concentration'] == conc]
    APD_conc = APD_conc.drop(
        columns='drug concentration').values.tolist()[0][1:]
    APD_conductance_concs.append(APD_conc)

nrow = int(np.ceil(len(drug_conc) / 4))
fig = modelling.figures.FigureStructure(
    figsize=(10, 2.5 * nrow),
    gridspec=(nrow, 4), hspace=0.3,
    wspace=0.1,
    height_ratios=[1] * nrow)

for i in range(len(drug_conc)):
    # APD_plot = [APD_trapping_concs[i][ind] for ind in
    #             range(len(APD_trapping_concs[i])) if ind % 2 == 0]
    APD_plot = APD_trapping_concs[i]
    pulse_num = len(APD_plot)
    fig.axs[int(i / 4)][i % 4].plot(np.arange(pulse_num), APD_plot, 'o',
                                    ms=0.9,
                                    label='state dependent\ndrug block',
                                    color='orange', zorder=-10)
    # APD_plot = [APD_conductance_concs[i][ind] for ind in
    #             range(len(APD_conductance_concs[i])) if ind % 2 == 0]
    APD_plot = APD_conductance_concs[i]
    fig.axs[int(i / 4)][i % 4].plot(np.arange(pulse_num), APD_plot, 'o',
                                    ms=0.9, label='conductance\nscaling',
                                    color='blue', zorder=-10)
    fig.axs[int(i / 4)][i % 4].set_title(str(drug_conc[i]) + 'nM')
    fig.axs[int(i / 4)][i % 4].set_rasterization_zorder(0)

lim_APD = []
for r in range(nrow):
    combine_APD = []
    for i in range(4):
        combine_APD += APD_trapping_concs[i + 4 * r]
        combine_APD += APD_conductance_concs[i + 4 * r]
    lim_APD.append((min(combine_APD) - 30, max(combine_APD) + 30))
lgnd = fig.axs[0][3].legend(loc='lower left', bbox_to_anchor=(1.0, 0),
                            handlelength=1)
for handle in lgnd.legendHandles:
    handle.set_markersize(6)
fig.sharex(['Sweeps'] * 4)
fig.sharey([r"APD$_{90}$ (ms)"] * nrow, lim_APD)

fig.savefig(saved_fig_dir + 'sensitivity_transientAPD.pdf')
