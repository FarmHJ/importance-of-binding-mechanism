import numpy as np
import pandas as pd

import modelling

testing_fig_dir = '../../figures/testing/'
# final_fig_dir = \
#     '../../figures/binding_kinetics_comparison/OHaraCiPA_model/protocol/'
final_fig_dir = '../../figures/conferences/'

saved_fig_dir = final_fig_dir

# Set up figure's main grid
# fig = modelling.figures.FigureStructure(figsize=(10, 6),
#                                         gridspec=(3, 1), hspace=0.6,
#                                         height_ratios=[1, 2, 2],
#                                         plot_in_subgrid=True)
fig = modelling.figures.FigureStructure(figsize=(10, 7.3),
                                        gridspec=(3, 1), hspace=0.7,
                                        height_ratios=[1, 2, 2],
                                        plot_in_subgrid=True)
plot = modelling.figures.FigurePlot()

subgridspecs = [(1, 4), (1, 2), (1, 4)]
wspaces = [0.55, 0.3, 0.1]
subgs = []
for i in range(3):
    subgs.append(fig.gs[i].subgridspec(*subgridspecs[i], wspace=wspaces[i]))
axs = [[[fig.fig.add_subplot(subgs[k][i, j]) for j in range(
    subgridspecs[k][1])] for i in range(subgridspecs[k][0])] for
    k in range(len(subgs))]

# Set up variables
protocol_params = modelling.ProtocolParameters()
protocol_name = ['Milnes', 'Pneg80', 'P0', 'P40']

# Top panel
panel1 = axs[0]

for i in range(len(protocol_name)):
    protocol = protocol_params.protocol_parameters[protocol_name[i]][
        'function']
    pulse_time = protocol_params.protocol_parameters[protocol_name[i]][
        'pulse_time']
    protocol_log = protocol.log_for_interval(0, pulse_time, for_drawing=True)

    for j, _ in enumerate(protocol_name):
        if i == j:
            panel1[0][i].plot(protocol_log['time'], protocol_log['pace'], 'k',
                              zorder=5)
        else:
            panel1[0][j].plot(protocol_log['time'], protocol_log['pace'],
                              '#cccccc', alpha=0.9)
    panel1[0][i].set_yticks(protocol_params.protocol_parameters[
        protocol_name[i]]['voltage_points'])
    panel1[0][i].set_title(protocol_name[i])
    panel1[0][i].set_ylabel('Voltage (mV)')
    panel1[0][i].spines['top'].set_visible(False)
    panel1[0][i].spines['right'].set_visible(False)
    fig.adjust_ticks(panel1[0][i], pulse_time)

fig.sharex(['Time (s)'] * 4, [(0, 25e3)] + [(0, 5400)] * 3,
           axs=panel1, subgridspec=subgridspecs[0])

# Middle panel
panel2 = axs[1]

# Set up variables
drugs = ['dofetilide', 'verapamil']
max_drug_conc = [3, 5]
Hill_model = modelling.HillsModel()
color = ['orange', 'blue', 'red', 'green']

for i in range(len(drugs)):
    saved_data_dir = '../../simulation_data/binding_kinetics_comparison/' + \
        drugs[i] + '/'
    Hill_coef_df = pd.read_csv(saved_data_dir + 'Hill_curves.csv')

    conc_grid = 10.0**np.linspace(-1, max_drug_conc[i], 20)

    for p in range(len(protocol_name)):
        Hill_eq = Hill_coef_df.loc[
            Hill_coef_df['protocol'] == protocol_name[p]]
        Hill_eq = Hill_eq.values.tolist()[0][1:-1]

        panel2[0][i].plot(conc_grid, Hill_model.simulate(
            Hill_eq, conc_grid), '-', color=color[p],
            label=protocol_name[p])

    panel2[0][i].set_xscale("log", nonpositive='clip')
    panel2[0][i].set_title(drugs[i] + '-like drug')
    panel2[0][i].set_xlabel('Drug concentration (nM)')
    panel2[0][i].set_ylabel('Normalised\npeak current')
    panel2[0][i].legend()

# Bottom panel
panel3 = axs[2]

# Set up variables
drug = 'verapamil'

for p, prot in enumerate(protocol_name):

    saved_data_dir = '../../simulation_data/binding_kinetics_comparison/' + \
        drug + '/' + prot + '/'

    APD = pd.read_csv(saved_data_dir + 'conductance_APD_pulses2.csv')

    drug_conc = APD['drug concentration'].values.tolist()
    APD = [max(APD.loc[i].values.tolist()[1:-1]) for i in range(APD.shape[0])]
    EAD_marker = [1070 if (i >= 1000) else None for i in APD[1:]]

    for q, _ in enumerate(protocol_name):
        if q == p:
            panel3[0][p].plot(drug_conc[1:], APD[1:], '^', color='blue',
                              zorder=5)
        else:
            panel3[0][q].plot(drug_conc[1:], APD[1:], '^', color='#cccccc',
                              alpha=0.9)
    panel3[0][p].plot(drug_conc[1:], EAD_marker, 'o', color='k',
                      marker=(5, 2), label='EAD-like AP')
    panel3[0][p].set_xscale("log", nonpositive='clip')
    panel3[0][p].set_xlabel('Drug concentration (nM)')
    panel3[0][p].set_title(prot)

panel3[0][0].legend(handlelength=1)
fig.sharey([r"APD$_{90}$ (ms)"],
           axs=panel3, subgridspec=subgridspecs[2])

fig.fig.text(0.075, 0.905, '(a)', fontsize=11)
fig.fig.text(0.075, 0.675, '(b)', fontsize=11)
fig.fig.text(0.535, 0.675, '(c)', fontsize=11)
fig.fig.text(0.075, 0.355, '(d)', fontsize=11)

fig.savefig(saved_fig_dir + "protocol_dependence.svg", format='svg')
