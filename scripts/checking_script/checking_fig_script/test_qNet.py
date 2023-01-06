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
    drug + '/' + protocol_name + '/qNet/'
fig_dir = '../../figures/model_comparison/' + drug + '/' + \
    protocol_name + '/qNet/'
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

fig = modelling.figures.FigureStructure(figsize=(10, 5), gridspec=(1, 2))
plot = modelling.figures.FigurePlot()

# Read files name of inet data
trapping_file_prefix = 'SD_inet_'
trapping_data_files = [f for f in os.listdir(data_dir) if
                       f.startswith(trapping_file_prefix)]
conductance_data_files = [f for f in os.listdir(data_dir) if
                          f.startswith('CS_inet_')]
# conc_label = [fname[len(trapping_file_prefix):-4] for fname in
#               trapping_data_files]
drug_conc = [float(fname[len(trapping_file_prefix):-4]) for fname in
             trapping_data_files]

sort_ind = [i[0] for i in sorted(enumerate(drug_conc), key=lambda x:x[1])]
drug_conc = sorted(drug_conc)
trapping_data_files = [trapping_data_files[i] for i in sort_ind]
conductance_data_files = [conductance_data_files[i] for i in sort_ind]

cmap = matplotlib.cm.get_cmap('viridis')

# Load action potential and APD data
SD_inet_log = []
CS_inet_log = []
for i in range(len(trapping_data_files)):
    SD_inet_log.append(pd.read_csv(data_dir + trapping_data_files[i]))
    # CS_inet_log.append(myokit.DataLog.load_csv(
    #     data_dir + conductance_data_files[i]))
    CS_inet_log.append(pd.read_csv(data_dir + conductance_data_files[i]))

# Initiate constants and variables
plotting_pulse_time = 2000

# Plot AP and IKr at various drug concentrations
for i in range(len(trapping_data_files)):
    fig.axs[0][0].plot(SD_inet_log[i]['time'], SD_inet_log[i]['inet'])
    fig.axs[0][1].plot(CS_inet_log[i]['time'], CS_inet_log[i]['inet'])
# plot.add_multiple_continuous(panel2[1][1], AP_conductance_plot,
#                              'ikr.IKr', cmap=cmap, labels=labels)
fig.axs[0][0].set_title('ORd-state dependent model')
fig.axs[0][1].set_title('ORd-conductance scaling model')

# # Adjust axes
fig.sharex(['Time (ms)'] * 2, [(0, plotting_pulse_time)] * 2)
fig.sharey(['Current (A/F)'])

fig.savefig(fig_dir + 'inet.pdf')


fig = modelling.figures.FigureStructure(figsize=(5, 5), gridspec=(1, 1))
plot = modelling.figures.FigurePlot()

# Load qNets data
qNets = pd.read_csv(data_dir + 'qNets.csv')

drug_conc = qNets['drug_conc']
SD_qNet = qNets['SD']
CS_qNet = qNets['CS']
# Plot APD90 of both models
fig.axs[0][0].plot(drug_conc, SD_qNet, 'o', color='orange',
                   label='ORd-SD model')
fig.axs[0][0].plot(drug_conc, CS_qNet, '^', color='blue',
                   label='ORd-CS model', alpha=0.8)
fig.axs[0][0].set_xscale("log", nonpositive='clip')
fig.axs[0][0].set_xlabel('Drug concentration (nM)')
# panel3[0][0].set_ylabel(r'APD$_{90}$ (ms)')
# panel3[0][0].legend(handlelength=1)

fig.savefig(fig_dir + 'qnet.pdf')