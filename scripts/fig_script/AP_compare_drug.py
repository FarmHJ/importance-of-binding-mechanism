# Compares the hERG current and action potential of state dependent drug block
# and conductance scaling
import pandas as pd
import matplotlib
import myokit
import numpy as np
import os

import modelling

drug = 'dofetilide'
protocol_name = 'Milnes'

saved_data_dir = '../../simulation_data/binding_kinetics_comparison/' + \
    drug + '/' + protocol_name + '/'

testing_fig_dir = '../../figures/testing/'
final_fig_dir = '../../figures/binding_kinetics_comparison/' + drug + '/' + \
    protocol_name + '/'

saved_fig_dir = final_fig_dir

fig = modelling.figures.FigureStructure(figsize=(10, 5), gridspec=(2, 2),
                                        height_ratios=[3, 4], hspace=0.4,
                                        width_ratios=[2.5, 1.2],
                                        plot_in_subgrid=True)
plot = modelling.figures.FigurePlot()

subgridspecs = [(2, 2), (1, 1)] * 2
subgs = []
# subgs.append(fig.gs[0].subgridspec(*subgridspecs[0], wspace=0.05,
#                                    hspace=0.4))
# subgs.append(fig.gs[1].subgridspec(*subgridspecs[1], wspace=0.05,
#                                    hspace=0.05))
# axs = [[[fig.fig.add_subplot(subgs[0][i, 0])] for
#         i in range(subgridspecs[0][0])]]
# axs.append([[fig.fig.add_subplot(subgs[0][:, 1])]])
# axs.append([[fig.fig.add_subplot(subgs[1][i, j]) for j in range(
#     subgridspecs[1][1])] for i in range(subgridspecs[1][0])])

subgs.append(fig.gs[0].subgridspec(*subgridspecs[0], wspace=0.08,
                                   hspace=0.08, height_ratios=[1, 2]))
for i in range(1, 2 * 2):
    subgs.append(fig.gs[i].subgridspec(*subgridspecs[i], wspace=0.08,
                                       hspace=0.08))
# axs = [[[fig.fig.add_subplot(subgs[k][i, j]) for j in range(
#          subgridspecs[k][1])] for i in range(subgridspecs[k][0])]
#        for k in [0, 2]]
# axs.append([[fig.fig.add_subplot(subgs[3][0, 0])]])
axs = [[[fig.fig.add_subplot(subgs[k][i, j]) for j in range(
         subgridspecs[k][1])] for i in range(subgridspecs[k][0])]
       for k in range(4)]

# Bottom panel
panel2 = axs[2]
# Load action potential data
trapping_data_files = [f for f in os.listdir(saved_data_dir) if
                       f.startswith('CiPA_AP_') and not
                       f.startswith('CiPA_AP_tran')]
conductance_data_files = [f for f in os.listdir(saved_data_dir) if
                          f.startswith('conductance_AP_') and not
                          f.startswith('conductance_AP_tran')]
conc_label = [fname[8:-4] for fname in trapping_data_files]
drug_conc = [float(fname[8:-4]) for fname in trapping_data_files]

# for verapamil
# removing_ind = drug_conc.index(500.0)
# drug_conc.pop(removing_ind)
# trapping_data_files.pop(removing_ind)
# conductance_data_files.pop(removing_ind)
# conc_label.pop(removing_ind)
# conc_label[-3] = r"$10^3$"
# conc_label[-2] = r"$10^4$"
# conc_label[-1] = r"$10^5$"

# Sort in increasing order of drug concentration
sort_ind = [i[0] for i in sorted(enumerate(drug_conc), key=lambda x:x[1])]
drug_conc = sorted(drug_conc)
trapping_data_files = [trapping_data_files[i] for i in sort_ind]
conductance_data_files = [conductance_data_files[i] for i in sort_ind]
conc_label = [conc_label[i] for i in sort_ind]

# Initiate constants and variables
labels = [i + ' nM' for i in conc_label]
# labels = labels[:-2]
cmap = matplotlib.cm.get_cmap('viridis')

CiPA_AP_log = []
conductance_AP_log = []
for i in range(len(trapping_data_files)):
    CiPA_AP_log.append(myokit.DataLog.load_csv(
        saved_data_dir + trapping_data_files[i]))
    conductance_AP_log.append(myokit.DataLog.load_csv(
        saved_data_dir + conductance_data_files[i]))

# APD_trapping = pd.read_csv(saved_data_dir + 'CiPA_APD_pulses2.csv')
# APD_conductance = pd.read_csv(saved_data_dir + 'conductance_APD_pulses2.csv')
APD_trapping = pd.read_csv(saved_data_dir + 'CiPA_APD_pulses1000.csv')
APD_conductance = pd.read_csv(saved_data_dir +
                              'conductance_APD_pulses1000.csv')

APD_trapping = [max(APD_trapping.loc[i].values.tolist()[1:-1]) for i in
                range(APD_trapping.shape[0])]
APD_conductance = [max(APD_conductance.loc[i].values.tolist()[1:-1]) for i in
                   range(APD_conductance.shape[0])]

# for dofetilide - because APD90 of 1000nM not run
APD_trapping.append(1000)
APD_conductance.append(1000)

# Initiate constants and variables
plotting_pulse_time = 1000 * 2

# Remove repeated signals at high concentrations
second_EAD_trap = [i for i, e in enumerate(APD_trapping)
                   if e == 1000][1:]
second_EAD_conduct = [i for i, e in enumerate(APD_conductance)
                      if e == 1000][1:]

if len(second_EAD_trap) <= len(second_EAD_conduct):
    chosen_conc_ind = second_EAD_trap
else:
    chosen_conc_ind = second_EAD_conduct

AP_trapping_plot = [e for i, e in enumerate(CiPA_AP_log)
                    if i not in chosen_conc_ind]
AP_conductance_plot = [e for i, e in enumerate(conductance_AP_log)
                       if i not in chosen_conc_ind]

# Plot figure
plot.add_multiple_continuous(panel2[0][0], AP_trapping_plot,
                             'membrane.V', cmap=cmap,
                             labels=labels)
plot.add_multiple_continuous(panel2[1][0], AP_trapping_plot,
                             'ikr.IKr', cmap=cmap, labels=labels)
plot.add_multiple_continuous(panel2[0][1], AP_conductance_plot,
                             'membrane.V', cmap=cmap,
                             labels=labels)
plot.add_multiple_continuous(panel2[1][1], AP_conductance_plot,
                             'ikr.IKr', cmap=cmap, labels=labels)
panel2[0][0].set_title('ORd-state dependent model')
panel2[0][1].set_title('ORd-conductance scaling model')

unique = fig.legend_without_duplicate_labels(panel2[1][1])
# panel2[1][1].legend(*zip(*unique), loc='lower left',
#                     bbox_to_anchor=(1.0, 0), handlelength=1)
fig.sharex(['Time (ms)'] * 2, [(0, plotting_pulse_time)] * 2,
           axs=panel2, subgridspec=subgridspecs[2])
fig.sharey(['Voltage (mV)', 'Current (A/F)'],
           axs=panel2, subgridspec=subgridspecs[2])

# Top left panel
panel1 = axs[0]
# Load hERG current data
trapping_data_files = [f for f in os.listdir(saved_data_dir) if
                       f.startswith('CiPA_hERG_')]
conductance_data_files = [f for f in os.listdir(saved_data_dir) if
                          f.startswith('conductance_hERG_')]
drug_conc = [float(fname[10:-4]) for fname in trapping_data_files]

# for verapamil
removing_ind = drug_conc.index(500.0)
trapping_data_files.pop(removing_ind)
conductance_data_files.pop(removing_ind)

# Sort in increasing order of drug concentration
trapping_data_files = [trapping_data_files[i] for i in sort_ind]
conductance_data_files = [conductance_data_files[i] for i in sort_ind]

CiPA_hERG_log = []
conductance_hERG_log = []
for i in range(len(trapping_data_files)):
    CiPA_hERG_log.append(myokit.DataLog.load_csv(
        saved_data_dir + trapping_data_files[i]))
    conductance_hERG_log.append(myokit.DataLog.load_csv(
        saved_data_dir + conductance_data_files[i]))

hERG_trapping_plot = [e for i, e in enumerate(CiPA_hERG_log)
                      if i not in chosen_conc_ind]
hERG_conductance_plot = [e for i, e in enumerate(conductance_hERG_log)
                         if i not in chosen_conc_ind]

# Initiate constants and variables
# labels = labels[:-2]
pulse_time = 25e3

# Plot figure
plot.add_single(panel1[0][0], hERG_trapping_plot[0], 'membrane.V', color='k')
plot.add_multiple(panel1[1][0], hERG_trapping_plot, 'ikr.IKr', labels=labels,
                  color=cmap)
plot.add_single(panel1[0][1], hERG_conductance_plot[0], 'membrane.V',
                color='k')
plot.add_multiple(panel1[1][1], hERG_conductance_plot, 'ikr.IKr',
                  labels=labels, color=cmap)

panel1[1][1].legend(handlelength=0.9, ncol=2, columnspacing=0.9, frameon=False)
panel1[0][0].set_title('SD model')
panel1[0][1].set_title('CS model')
fig.sharex(['Time (s)'] * 2, [(0, pulse_time)] * 2,
           axs=panel1, subgridspec=subgridspecs[0])
fig.sharey(['Voltage (mV)', 'Current (A/F)'],
           axs=panel1, subgridspec=subgridspecs[0])
for i in range(2):
    panel1[i][0].spines['top'].set_visible(False)
    panel1[i][0].spines['right'].set_visible(False)
    panel1[i][1].spines['top'].set_visible(False)
    panel1[i][1].spines['right'].set_visible(False)
    fig.adjust_ticks(panel1[i][0], pulse_time)
    fig.adjust_ticks(panel1[i][1], pulse_time)

# Bottom right panel
panel3 = axs[3]

APD_trapping = pd.read_csv(saved_data_dir + 'CiPA_APD_fine.csv')
APD_conductance = pd.read_csv(saved_data_dir + 'conductance_APD_fine.csv')
# APD_trapping = pd.read_csv(saved_data_dir + 'CiPA_APD_pulses1000.csv')
# APD_conductance = pd.read_csv(saved_data_dir +
#                               'conductance_APD_pulses1000.csv')

drug_conc = APD_trapping['drug concentration'].values.tolist()
APD_trapping = [max(APD_trapping.loc[i].values.tolist()[1:-1]) for i in
                range(APD_trapping.shape[0])]
APD_conductance = [max(APD_conductance.loc[i].values.tolist()[1:-1])
                   for i in range(APD_conductance.shape[0])]
EAD_marker = [1050 if (i >= 1000 or j >= 1000) else None for (i, j)
              in zip(APD_trapping[1:], APD_conductance[1:])]

panel3[0][0].plot(drug_conc[1:], APD_trapping[1:], 'o', color='orange',
                  label='ORd-SD model')
panel3[0][0].plot(drug_conc[1:], APD_conductance[1:], '^', color='blue',
                  label='ORd-CS model', alpha=0.8)
panel3[0][0].plot(drug_conc[1:], EAD_marker, 'o', color='k',
                  marker=(5, 2), label='EAD-like AP')
panel3[0][0].set_xscale("log", nonpositive='clip')
panel3[0][0].set_xlabel('Drug concentration (nM)')
panel3[0][0].set_ylabel(r'APD$_{90}$ (ms)')
panel3[0][0].legend(handlelength=1)

# Top right panel
panel4 = axs[1]

data_dir = '../../data/'
exp_data = pd.read_csv(data_dir + drug + '.csv',
                       header=[0], index_col=[0],
                       skipinitialspace=True)
control_log = myokit.DataLog.load_csv(data_dir + 'control_pulses10.csv')

cmap = matplotlib.colors.ListedColormap(
    matplotlib.cm.get_cmap('tab10')(range(4)))

# Load hERG current data
sim_data_files = [f for f in os.listdir(data_dir) if
                  f.startswith(drug + '_')]
drug_conc = [float(fname[len(drug) + 5:-13]) for fname in sim_data_files]
sort_ind = [i[0] for i in sorted(enumerate(drug_conc), key=lambda x:x[1])]

# Sort in increasing order of drug concentration
sim_data_files = [sim_data_files[i] for i in sort_ind]

sim_log = []
for i in range(len(sim_data_files)):
    sim_log.append(myokit.DataLog.load_csv(
        data_dir + sim_data_files[i]))

for i, conc_i in enumerate(drug_conc):

    current = exp_data.loc[exp_data['conc'] == conc_i]
    log = sim_log[i]

    max_exp = max(current['exp'].values)
    exp_repeats = []
    for exp in range(1, max_exp):
        current_exp = current.loc[current['exp'] == exp]

        max_sweep = max(current_exp['sweep'].values)
        current_temp = current_exp.loc[current_exp['sweep'] == max_sweep]
        exp_repeats.append(current_temp['frac'].values)

    min_time = min(current_temp.index)
    max_time = max(current_temp.index)
    log_range_min = np.argmin(np.abs(np.array(log.time()) - min_time))
    log_range_max = np.argmin(np.abs(np.array(log.time()) - max_time))

    log_plot = log['ikr.IKr'][log_range_min:log_range_max + 1]
    control_log_plot = control_log['ikr.IKr'][log_range_min:log_range_max + 1]

    panel4[0][0].plot(current_temp.index - current_temp.index[0],
             np.mean(exp_repeats, axis=0), 'o', ms=1.2,
             zorder=-10, color=cmap(i), label=str(int(conc_i)) + ' nM')
    panel4[0][0].fill_between(
        current_temp.index - current_temp.index[0], np.mean(exp_repeats, axis=0) -
        np.std(exp_repeats, axis=0), np.mean(exp_repeats, axis=0) +
        np.std(exp_repeats, axis=0), color=cmap(i), alpha=.3)

    panel4[0][0].plot(np.array(log.time()[log_range_min:log_range_max + 1]) - log.time()[log_range_min],
             np.array(log_plot) / np.array(control_log_plot), zorder=1, color='k')

panel4[0][0].set_xlim(min_time  - current_temp.index[0], max_time - current_temp.index[0])
panel4[0][0].set_ylim(top=1.3)
panel4[0][0].set_xlabel('Time (s)')
panel4[0][0].set_ylabel('Fractional block')
fig.adjust_ticks(panel4[0][0], 1e4)
lgnd = panel4[0][0].legend(ncol=2, columnspacing=1, handlelength=1)
for handle in lgnd.legendHandles:
    handle.set_markersize(4)

# Add panel letter
fig.fig.set_size_inches(10, 5.5)
fig.fig.text(0.1, 0.905, '(A)', fontsize=11)
fig.fig.text(0.1, 0.495, '(B)', fontsize=11)
fig.fig.text(0.64, 0.905, '(C)', fontsize=11)
fig.fig.text(0.64, 0.495, '(D)', fontsize=11)

fig.savefig(saved_fig_dir + "model_compare_test.svg", format='svg')
# fig.savefig(saved_fig_dir + "grid_test.pdf")
