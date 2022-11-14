# Plots figures to see how APD90 changes with different drugs
import matplotlib
import myokit
import numpy as np
import os
import pandas as pd
import sys

import modelling


# Set up directory for figures
testing_fig_dir = '../../figures/testing/'
final_fig_dir = '../../figures/binding_kinetics_comparison/'

saved_fig_dir = final_fig_dir

fig = modelling.figures.FigureStructure(
    figsize=(10, 5),
    gridspec=(2, 3), hspace=1,
    wspace=0.9,
    height_ratios=[1] * 2,
    plot_in_subgrid=True)
plot = modelling.figures.FigurePlot()

subgridspecs = [(1, 1)] * 5
subgs = []
for i in range(5):
    subgs.append(fig.gs[i + 1].subgridspec(*subgridspecs[i], wspace=0.05,
                                           hspace=0.05))
axs = [[[fig.fig.add_subplot(subgs[k][i, j]) for j in range(
    subgridspecs[k][1])] for i in range(subgridspecs[k][0])] for
    k in range(len(subgs))]

drug = 'verapamil'
protocol_name = 'Milnes'
saved_data_dir = '../../simulation_data/binding_kinetics_comparison/' + \
    drug + '/' + protocol_name + '/'

# Showing hERG currents
panel1 = axs[0]
panel2 = axs[3]
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

peaks = []
peaks_pos = []
peaks_firstpulse = []
peaks_pos_firstpulse = []
for log in CiPA_hERG_log:
    peaks_firstpulse.append(np.max(log['ikr.IKr'][:250000]))
    peaks_pos_firstpulse.append(log.time()[np.argmax(log['ikr.IKr'][:250000])])
    peaks.append(np.max(log['ikr.IKr']))
    peaks_pos.append(log.time()[np.argmax(log['ikr.IKr'])])

# Initiate constants and variables
labels = [i + ' nM' for i in conc_label]
pulse_time = 25e3
cmap = matplotlib.cm.get_cmap('viridis')

# Plot figure
plot.add_multiple(panel1[0][0], CiPA_hERG_log, 'ikr.IKr', labels=labels,
                  color=cmap)
panel1[0][0].plot(peaks_pos_firstpulse[1:], peaks_firstpulse[1:], 'kx')
plot.add_multiple(panel2[0][0], conductance_hERG_log, 'ikr.IKr',
                  labels=labels, color=cmap)

# panel1[1][0].legend()
panel1[0][0].set_title('State-dependent drug block')
panel2[0][0].set_title('Conductance scaling drug block')
fig.sharex(['Time (s)'], [(0, pulse_time)],
           axs=panel1, subgridspec=subgridspecs[0])
fig.sharex(['Time (s)'], [(0, pulse_time)],
           axs=panel2, subgridspec=subgridspecs[3])
fig.sharey(['Current (A/F)'],
           axs=panel1, subgridspec=subgridspecs[0])
fig.sharey(['Current (A/F)'],
           axs=panel2, subgridspec=subgridspecs[3])
fig.adjust_ticks(panel1[0][0], pulse_time)
fig.adjust_ticks(panel2[0][0], pulse_time)

# Showing Hill curve
panel3 = axs[1]
panel4 = axs[4]

peaks = (peaks - min(peaks)) / (max(peaks) - min(peaks))
panel3[0][0].plot(drug_conc[1:], peaks[1:], 'x', color='k')

panel3[0][0].set_xscale("log", nonpositive='clip')
panel3[0][0].set_xlabel('Drug concentration (nM)')
panel3[0][0].set_ylabel('Normalised peak current')

Hill_model = modelling.HillsModel()
Hill_coef_df = pd.read_csv(saved_data_dir + '../Hill_curves.csv')

drug_conc_lib = modelling.DrugConcentrations()
conc_grid = drug_conc_lib.drug_concentrations[drug]['fine']

Hill_eq = Hill_coef_df.loc[
    Hill_coef_df['protocol'] == protocol_name]
Hill_eq = Hill_eq.values.tolist()[0][1:-1]

panel4[0][0].plot(conc_grid, Hill_model.simulate(
    Hill_eq, conc_grid), linestyle='-', color='k')

panel4[0][0].set_xscale("log", nonpositive='clip')
panel4[0][0].set_xlabel('Drug concentration (nM)')
panel4[0][0].set_ylabel('Ionic conductance scale')
panel4[0][0].set_title('Fitted Hill curve')

# Show example AP
panel5 = axs[2]
conductance_data_files = [f for f in os.listdir(saved_data_dir) if
                          f.startswith('conductance_AP_') and not
                          f.startswith('conductance_AP_transient')]
conc_label = [fname[15:-4] for fname in conductance_data_files]
drug_conc = [float(fname[15:-4]) for fname in conductance_data_files]

# Sort in increasing order of drug concentration
sort_ind = [i[0] for i in sorted(enumerate(drug_conc), key=lambda x:x[1])]
drug_conc = sorted(drug_conc)
conductance_data_files = [conductance_data_files[i] for i in sort_ind]
conc_label = [conc_label[i] for i in sort_ind]
labels = [i + ' nM' for i in conc_label]

conductance_AP_log = []
for i in range(len(conductance_data_files)):
    conductance_AP_log.append(myokit.DataLog.load_csv(
        saved_data_dir + conductance_data_files[i]))

# Initiate constants and variables
plotting_pulse_time = 1000 * 2

# Plot figure
plot.add_multiple_continuous(panel5[0][0], conductance_AP_log,
                             'membrane.V', cmap=cmap,
                             labels=labels)
panel5[0][0].set_title('ORd-conductance scaling')

fig.sharex(['Time (ms)'], [(0, plotting_pulse_time)],
           axs=panel5, subgridspec=subgridspecs[2])
fig.sharey(['Voltage\n(mV)'],
           axs=panel5, subgridspec=subgridspecs[2])

for i in range(5):
    axs[i][0][0].spines['top'].set_visible(False)
    axs[i][0][0].spines['right'].set_visible(False)
fig.savefig(saved_fig_dir + 'method_diagram.svg', format='svg')
