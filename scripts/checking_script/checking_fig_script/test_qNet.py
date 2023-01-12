import matplotlib
import myokit
import os
import pandas as pd
import sys

import modelling

# Define drug and protocol
# drug = sys.argv[1]
param_lib = modelling.BindingParameters()
# Define drug
drug_list = param_lib.drug_compounds[:-1]
protocol_name = 'Milnes'

for drug in drug_list:
    # Define directories to read data and save plotted figures
    data_dir = '../../simulation_data/model_comparison/' + \
        drug + '/' + protocol_name + '/qNet/'
    fig_dir = '../../figures/model_comparison/' + drug + '/' + \
        protocol_name + '/qNet/'
    if not os.path.isdir(fig_dir):
        os.makedirs(fig_dir)

    # fig = modelling.figures.FigureStructure(figsize=(10, 5), gridspec=(1, 2))
    fig = modelling.figures.FigureStructure(figsize=(4, 5), gridspec=(1, 1))
    plot = modelling.figures.FigurePlot()

    # Read files name of inet data
    # trapping_file_prefix = 'SD_inet_'
    trapping_file_prefix = 'SD_AP_inetcurrents_multiion_'
    trapping_data_files = [f for f in os.listdir(data_dir) if
                        f.startswith(trapping_file_prefix)]
    # conductance_data_files = [f for f in os.listdir(data_dir) if
    #                         f.startswith('CS_AP_inetcurrents_')]
    # conc_label = [fname[len(trapping_file_prefix):-4] for fname in
    #               trapping_data_files]
    drug_conc = [float(fname[len(trapping_file_prefix):-4]) for fname in
                trapping_data_files]

    sort_ind = [i[0] for i in sorted(enumerate(drug_conc), key=lambda x:x[1])]
    drug_conc = sorted(drug_conc)
    trapping_data_files = [trapping_data_files[i] for i in sort_ind]
    # conductance_data_files = [conductance_data_files[i] for i in sort_ind]

    cmap = matplotlib.cm.get_cmap('viridis')

    # Load action potential and APD data
    SD_inet_log = []
    # CS_inet_log = []
    for i in range(len(trapping_data_files)):
        # SD_inet_log.append(pd.read_csv(data_dir + trapping_data_files[i]))
        SD_inet_log.append(myokit.DataLog.load_csv(
            data_dir + trapping_data_files[i]))
        # CS_inet_log.append(myokit.DataLog.load_csv(
        #     data_dir + conductance_data_files[i]))
        # CS_inet_log.append(pd.read_csv(data_dir + conductance_data_files[i]))

    # Initiate constants and variables
    plotting_pulse_time = 2000

    # Plot AP and IKr at various drug concentrations
    for i in range(len(trapping_data_files)):
        # fig.axs[0][0].plot(SD_inet_log[i]['time'], SD_inet_log[i]['inet'])
        # fig.axs[0][1].plot(CS_inet_log[i]['time'], CS_inet_log[i]['inet'])
        fig.axs[0][0].plot(SD_inet_log[i].time(), SD_inet_log[i]['membrane.V'])
        # fig.axs[0][1].plot(CS_inet_log[i].time(), CS_inet_log[i]['membrane.V'])
    # plot.add_multiple_continuous(panel2[1][1], AP_conductance_plot,
    #                              'ikr.IKr', cmap=cmap, labels=labels)
    fig.axs[0][0].set_title('ORd-state dependent model')
    # fig.axs[0][1].set_title('ORd-conductance scaling model')

    # # Adjust axes
    # fig.sharex(['Time (ms)'] * 2, [(0, plotting_pulse_time)] * 2)
    # fig.sharey(['Current (A/F)'])

    # fig.savefig(fig_dir + 'inet.pdf')
    fig.savefig(fig_dir + 'AP_qnet.pdf')


# fig = modelling.figures.FigureStructure(figsize=(5, 5), gridspec=(1, 1))
# plot = modelling.figures.FigurePlot()

# # Load qNets data
# qNets = pd.read_csv(data_dir + 'qNets.csv')

# param_lib = modelling.BindingParameters()
# Cmax = param_lib.binding_parameters[drug]['Cmax']

# drug_conc = qNets['drug_conc']
# drug_conc_multiple = drug_conc / Cmax
# SD_qNet = qNets['SD']
# CS_qNet = qNets['CS']
# # Plot APD90 of both models
# fig.axs[0][0].plot(drug_conc_multiple, SD_qNet, 'o', color='orange',
#                    label='ORd-SD model')
# fig.axs[0][0].plot(drug_conc_multiple, CS_qNet, '^', color='blue',
#                    label='ORd-CS model', alpha=0.8)
# # fig.axs[0][0].set_xscale("log", nonpositive='clip')
# fig.axs[0][0].set_xlabel(r"Drug concentration ($\times \mathrm{C}_\mathrm{max}$)")
# fig.axs[0][0].set_ylabel('qNet (As/F)')
# fig.axs[0][0].legend(handlelength=1)

# fig.savefig(fig_dir + 'qnet.pdf')

# Plot qNet for AP-SD model for all drugs
fig = modelling.figures.FigureStructure(figsize=(5, 5), gridspec=(1, 1))
plot = modelling.figures.FigurePlot()

fig_dir = '../../figures/model_comparison/'

marker_color_dict = {
    'quinidine': {'m': 'o', 'c': 'red'},
    'dofetilide': {'m': 'd', 'c': 'red'},
    'bepridil': {'m': 's', 'c': 'red'},
    'sotalol': {'m': '^', 'c': 'red'},
    'chlorpromazine': {'m': 'o', 'c': 'blue'},
    'cisapride': {'m': 'd', 'c': 'blue'},
    'terfenadine': {'m': 's', 'c': 'blue'},
    'ondansetron': {'m': '^', 'c': 'blue'},
    'diltiazem': {'m': 'o', 'c': 'green'},
    'mexiletine': {'m': 'd', 'c': 'green'},
    'ranolazine': {'m': 's', 'c': 'green'},
    'verapamil': {'m': '^', 'c': 'green'}}

drug_list.remove('quinidine')
# drug_list = ['dofetilide', 'bepridil', 'terfenadine', 'verapamil',
#              'ranolazine', 'mexiletine']
for drug in drug_list:
    # Load qNets data
    data_dir = '../../simulation_data/model_comparison/' + \
        drug + '/' + protocol_name + '/qNet/'
    qNets = pd.read_csv(data_dir + 'qNets_multiion.csv')

    Cmax = param_lib.binding_parameters[drug]['Cmax']

    drug_conc = qNets['drug_conc']
    drug_conc_multiple = drug_conc / Cmax
    SD_qNet = qNets['SD']
    # CS_qNet = qNets['CS']
    # Plot APD90 of both models
    marker = marker_color_dict[drug]['m']
    color = marker_color_dict[drug]['c']
    fig.axs[0][0].plot(drug_conc_multiple, SD_qNet, marker + '-', label=drug,
                       color=color)
                    # , color='orange',
                    # label='ORd-SD model')
    # fig.axs[0][0].plot(drug_conc_multiple, CS_qNet, '^', color='blue',
    #                 label='ORd-CS model', alpha=0.8)
    # fig.axs[0][0].set_xscale("log", nonpositive='clip')
fig.axs[0][0].set_xlabel(r"Drug concentration ($\times \mathrm{C}_\mathrm{max}$)")
fig.axs[0][0].set_ylabel('qNet (As/F)')
fig.axs[0][0].legend(handlelength=1)

fig.savefig(fig_dir + 'qnet_multiion.pdf')
