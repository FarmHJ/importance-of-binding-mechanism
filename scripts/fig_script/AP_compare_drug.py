# Compares the IKr, AP and APD90 of the SD model and the CS model
import matplotlib
import myokit
import os
import pandas as pd
import sys

import modelling

# Define drug and protocol
drug = sys.argv[1]
protocol_name = 'Milnes'

# Define directories to read data and save plotted figures
data_dir = '../../simulation_data/model_comparison/' + \
    drug + '/' + protocol_name + '/'
fig_dir = '../../figures/model_comparison/' + drug + '/' + \
    protocol_name + '/'
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

# Set up structure of the figure
fig = modelling.figures.FigureStructure(figsize=(10, 5), gridspec=(2, 2),
                                        height_ratios=[3, 4], hspace=0.4,
                                        width_ratios=[2.5, 1.2],
                                        plot_in_subgrid=True)
plot = modelling.figures.FigurePlot()

subgridspecs = [(2, 2), (1, 1)] * 2
subgs = []
subgs.append(fig.gs[0].subgridspec(*subgridspecs[0], wspace=0.08,
                                   hspace=0.08, height_ratios=[1, 2]))
for i in range(1, 2 * 2):
    subgs.append(fig.gs[i].subgridspec(*subgridspecs[i], wspace=0.08,
                                       hspace=0.08))
axs = [[[fig.fig.add_subplot(subgs[k][i, j]) for j in range(
         subgridspecs[k][1])] for i in range(subgridspecs[k][0])]
       for k in [0, 2]]
axs.append([[fig.fig.add_subplot(subgs[3][0, 0])]])

# Bottom left panel
# Plot action potentials and the corresponding IKr of the ORd-SD model and the
# ORd-CS model stimulated to steady state
panel2 = axs[1]

# Read files name of action potential data
trapping_data_files = [f for f in os.listdir(data_dir) if
                       f.startswith('SD_AP_') and not
                       f.startswith('SD_AP_tran')]
conductance_data_files = [f for f in os.listdir(data_dir) if
                          f.startswith('CS_AP_') and not
                          f.startswith('CS_AP_tran')]
conc_label = [fname[6:-4] for fname in trapping_data_files]
drug_conc = [float(fname[6:-4]) for fname in trapping_data_files]

# Sort in increasing order of drug concentration
sort_ind = [i[0] for i in sorted(enumerate(drug_conc), key=lambda x:x[1])]
drug_conc = sorted(drug_conc)
trapping_data_files = [trapping_data_files[i] for i in sort_ind]
conductance_data_files = [conductance_data_files[i] for i in sort_ind]
conc_label = [conc_label[i] for i in sort_ind]

if drug == 'verapamil':
    conc_label[-2] = r"$10^4$"
    conc_label[-1] = r"$10^5$"

# Initiate constants and variables
labels = [i + ' nM' for i in conc_label]
cmap = matplotlib.cm.get_cmap('viridis')

# Load action potential and APD data
trapping_AP_log = []
conductance_AP_log = []
for i in range(len(trapping_data_files)):
    trapping_AP_log.append(myokit.DataLog.load_csv(
        data_dir + trapping_data_files[i]))
    conductance_AP_log.append(myokit.DataLog.load_csv(
        data_dir + conductance_data_files[i]))

APD_trapping = pd.read_csv(data_dir + 'SD_APD_pulses2.csv')
APD_conductance = pd.read_csv(data_dir + 'CS_APD_pulses2.csv')

APD_trapping = [max(APD_trapping.loc[i].values.tolist()[1:-1]) for i in
                range(APD_trapping.shape[0])]
APD_conductance = [max(APD_conductance.loc[i].values.tolist()[1:-1]) for i in
                   range(APD_conductance.shape[0])]

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

AP_trapping_plot = [e for i, e in enumerate(trapping_AP_log)
                    if i not in chosen_conc_ind]
AP_conductance_plot = [e for i, e in enumerate(conductance_AP_log)
                       if i not in chosen_conc_ind]

# Plot AP and IKr at various drug concentrations
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

# Adjust axes
fig.sharex(['Time (ms)'] * 2, [(0, plotting_pulse_time)] * 2,
           axs=panel2, subgridspec=subgridspecs[2])
fig.sharey(['Voltage (mV)', 'Current (A/F)'],
           axs=panel2, subgridspec=subgridspecs[2])

# Top left panel
# Shows Milnes' protocol and the IKr stimulated by it at various drug
# concentration for both the SD model and the CS model
panel1 = axs[0]
# Read files name of IKr data
trapping_data_files = [f for f in os.listdir(data_dir) if
                       f.startswith('SD_current_')]
conductance_data_files = [f for f in os.listdir(data_dir) if
                          f.startswith('CS_current_')]
drug_conc = [float(fname[11:-4]) for fname in trapping_data_files]

# Sort in increasing order of drug concentration
trapping_data_files = [trapping_data_files[i] for i in sort_ind]
conductance_data_files = [conductance_data_files[i] for i in sort_ind]

# Load IKr data
trapping_hERG_log = []
conductance_hERG_log = []
for i in range(len(trapping_data_files)):
    trapping_hERG_log.append(myokit.DataLog.load_csv(
        data_dir + trapping_data_files[i]))
    conductance_hERG_log.append(myokit.DataLog.load_csv(
        data_dir + conductance_data_files[i]))

hERG_trapping_plot = [e for i, e in enumerate(trapping_hERG_log)
                      if i not in chosen_conc_ind]
hERG_conductance_plot = [e for i, e in enumerate(conductance_hERG_log)
                         if i not in chosen_conc_ind]

# Initiate constants and variables
pulse_time = 25e3

# Plot Milnes' protocol and Ikr
plot.add_single(panel1[0][0], hERG_trapping_plot[0], 'membrane.V', color='k')
plot.add_multiple(panel1[1][0], hERG_trapping_plot, 'ikr.IKr', labels=labels,
                  color=cmap)
plot.add_single(panel1[0][1], hERG_conductance_plot[0], 'membrane.V',
                color='k')
plot.add_multiple(panel1[1][1], hERG_conductance_plot, 'ikr.IKr',
                  labels=labels, color=cmap)

# Adjust figure details
panel1[1][1].legend(handlelength=0.9, ncol=2, columnspacing=0.9)
panel1[0][0].set_title('SD model')
panel1[0][1].set_title('CS model')
fig.sharex(['Time (s)'] * 2, [(0, pulse_time)] * 2,
           axs=panel1, subgridspec=subgridspecs[0])
fig.sharey(['Voltage (mV)', 'Current (A/F)'],
           axs=panel1, subgridspec=subgridspecs[0])
for i in range(2):
    fig.adjust_ticks(panel1[i][0], pulse_time)
    fig.adjust_ticks(panel1[i][1], pulse_time)

# Bottom right panel
# Plots the APD90 calculated for both models with drugs at various drug
# concentrations
panel3 = axs[2]

# Load APD data
APD_trapping = pd.read_csv(data_dir + 'SD_APD_fine.csv')
APD_conductance = pd.read_csv(data_dir + 'CS_APD_fine.csv')

# Identify EAD-like behaviour
drug_conc = APD_trapping['drug concentration'].values.tolist()
APD_trapping = APD_trapping['APD'].values.tolist()
APD_conductance = APD_conductance['APD'].values.tolist()
EAD_marker = [1050 if (i >= 1000 or j >= 1000) else None for (i, j)
              in zip(APD_trapping, APD_conductance)]

# Plot APD90 of both models
panel3[0][0].plot(drug_conc, APD_trapping, 'o', color='orange',
                  label='ORd-SD model')
panel3[0][0].plot(drug_conc, APD_conductance, '^', color='blue',
                  label='ORd-CS model', alpha=0.8)
panel3[0][0].scatter(drug_conc, EAD_marker, marker=(5, 2),
                     color='k', label='EAD-like AP')
panel3[0][0].set_xscale("log", nonpositive='clip')
panel3[0][0].set_xlabel('Drug concentration (nM)')
panel3[0][0].set_ylabel(r'APD$_{90}$ (ms)')
panel3[0][0].legend(handlelength=1)

# Add panel letter
fig.fig.set_size_inches(10, 5.5)
fig.fig.text(0.1, 0.905, '(A)', fontsize=11)
fig.fig.text(0.1, 0.495, '(B)', fontsize=11)
fig.fig.text(0.64, 0.495, '(C)', fontsize=11)

fig.savefig(fig_dir + "model_compare.svg", format='svg')
