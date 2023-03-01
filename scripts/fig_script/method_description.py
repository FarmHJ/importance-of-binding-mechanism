#
# Figure 2
# Describe the fitting procedure of the CS model to the SD model.
#

import matplotlib
import myokit
import numpy as np
import os
import pandas as pd

import modelling


# Set up directory for figures
fig_dir = '../../figures/background/'
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

# Set up structure of the figure
fig = modelling.figures.FigureStructure(figsize=(10, 5),
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

# Choose an example
drug = 'verapamil'
protocol_name = 'Milnes'
data_dir = '../../simulation_data/model_comparison/' + drug + '/' + \
    protocol_name + '/'

# Panels that show hERG currents
panel1 = axs[0]
panel2 = axs[3]

# Load hERG current data
SD_fileprefix = 'SD_current_'
CS_fileprefix = 'CS_current_'
SD_data_files = [f for f in os.listdir(data_dir) if
                 f.startswith(SD_fileprefix)]
CS_data_files = [f for f in os.listdir(data_dir) if
                 f.startswith(CS_fileprefix)]
drug_conc = [float(fname[len(SD_fileprefix):-4]) for fname in
             SD_data_files]
drug_conc_CS = [float(fname[len(CS_fileprefix):-4]) for fname in
                CS_data_files]

# Sort in increasing order of drug concentration
sort_ind = [i[0] for i in sorted(enumerate(drug_conc), key=lambda x:x[1])]
sort_ind_CS = [i[0] for i in sorted(enumerate(drug_conc_CS),
               key=lambda x:x[1])]
drug_conc = sorted(drug_conc)
trapping_data_files = [SD_data_files[i] for i in sort_ind]
conductance_data_files = [CS_data_files[i] for i in sort_ind_CS]

trapping_hERG_log = []
conductance_hERG_log = []
for i in range(len(trapping_data_files)):
    trapping_hERG_log.append(myokit.DataLog.load_csv(
        data_dir + trapping_data_files[i]))
    conductance_hERG_log.append(myokit.DataLog.load_csv(
        data_dir + conductance_data_files[i]))

# Initiate constants and variables
pulse_time = 25e3
cmap = matplotlib.cm.get_cmap('viridis')

peaks = []
peaks_pos = []
peaks_firstpulse = []
peaks_pos_firstpulse = []
for log in trapping_hERG_log:
    peaks_firstpulse.append(np.max(log['ikr.IKr'][:pulse_time * 10]))
    peaks_pos_firstpulse.append(log.time()[np.argmax(log['ikr.IKr']
                                                     [:pulse_time * 10])])
    peaks.append(np.max(log['ikr.IKr']))
    peaks_pos.append(log.time()[np.argmax(log['ikr.IKr'])])

# Plot figure
plot.add_multiple(panel1[0][0], trapping_hERG_log, 'ikr.IKr', color=cmap)
panel1[0][0].plot(peaks_pos_firstpulse[1:], peaks_firstpulse[1:], 'kx')
plot.add_multiple(panel2[0][0], conductance_hERG_log, 'ikr.IKr', color=cmap)

# Adjust panels
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

# Panels that show the Hill curve
panel3 = axs[1]
panel4 = axs[4]

# Plot peaks of hERG current
peaks = (peaks - min(peaks)) / (max(peaks) - min(peaks))
panel3[0][0].plot(drug_conc[1:], peaks[1:], 'x', color='k')

# Adjust axes
panel3[0][0].set_xscale("log", nonpositive='clip')
panel3[0][0].set_xlabel('Drug concentration (nM)')
panel3[0][0].set_ylabel('Normalised peak current')

# Load coefficients for Hill curve
Hill_model = modelling.HillsModel()
Hill_coef_df = pd.read_csv(data_dir + '../Hill_curves.csv')

Hill_eq = Hill_coef_df.loc[
    Hill_coef_df['protocol'] == protocol_name]
Hill_eq = Hill_eq.values.tolist()[0][1:-1]

# Define drug concentration
drug_conc_lib = modelling.DrugConcentrations()
conc_grid = drug_conc_lib.drug_concentrations[drug]['fine']

# Plot peaks of hERG current and the Hill curve
panel4[0][0].plot(drug_conc[1:], peaks[1:], 'x', color='k')
panel4[0][0].plot(conc_grid, Hill_model.simulate(
    Hill_eq, conc_grid), linestyle='-', color='k')

# Adjust axes
panel4[0][0].set_xscale("log", nonpositive='clip')
panel4[0][0].set_xlabel('Drug concentration (nM)')
panel4[0][0].set_ylabel('Ionic conductance scale')
panel4[0][0].set_title('Fitted Hill curve')

# Panel to show an example AP
panel5 = axs[2]

# Load AP data
CS_AP_prefix = 'CS_AP_'
conductance_data_files = [f for f in os.listdir(data_dir) if
                          f.startswith(CS_AP_prefix) and not
                          f.startswith('CS_AP_transient')]
drug_conc = [float(fname[len(CS_AP_prefix):-4]) for fname in
             conductance_data_files]

# Sort in increasing order of drug concentration
sort_ind = [i[0] for i in sorted(enumerate(drug_conc), key=lambda x:x[1])]
drug_conc = sorted(drug_conc)
conductance_data_files = [conductance_data_files[i] for i in sort_ind]

conductance_AP_log = []
for i in range(len(conductance_data_files)):
    conductance_AP_log.append(myokit.DataLog.load_csv(
        data_dir + conductance_data_files[i]))

# Initiate constants and variables
plotting_pulse_time = 1000 * 2

# Plot figure
plot.add_multiple_continuous(panel5[0][0], conductance_AP_log,
                             'membrane.V', cmap=cmap)
panel5[0][0].set_title('AP-conductance scaling')

# Adjust axes
fig.sharex(['Time (ms)'], [(0, plotting_pulse_time)],
           axs=panel5, subgridspec=subgridspecs[2])
fig.sharey(['Voltage\n(mV)'],
           axs=panel5, subgridspec=subgridspecs[2])

for i in range(5):
    axs[i][0][0].spines['top'].set_visible(False)
    axs[i][0][0].spines['right'].set_visible(False)

# Save figure
fig.savefig(fig_dir + 'method_diagram.svg', format='svg')
