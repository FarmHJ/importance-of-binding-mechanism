# Compares the hERG current and action potential of state dependen drug block
# and conductance scaling
import pandas as pd
import matplotlib
import myokit
import os

import modelling

drug = 'verapamil'
protocol_name = 'Milnes'

saved_data_dir = '../../simulation_data/binding_kinetics_comparison/' + \
    drug + '/' + protocol_name + '/'

testing_fig_dir = '../../figures/testing/'
final_fig_dir = '../../figures/binding_kinetics_comparison/' + drug + '/' + \
    protocol_name + '/'

saved_fig_dir = final_fig_dir

fig = modelling.figures.FigureStructure(figsize=(10, 5), gridspec=(2, 1),
                                        height_ratios=[2, 3], hspace=0.5,
                                        plot_in_subgrid=True)
plot = modelling.figures.FigurePlot()

subgridspecs = [(2, 2), (2, 2)]
subgs = []
subgs.append(fig.gs[0].subgridspec(*subgridspecs[0], wspace=0.05,
                                   hspace=0.05, height_ratios=[1, 2]))
subgs.append(fig.gs[1].subgridspec(*subgridspecs[1], wspace=0.05,
                                   hspace=0.05))
axs = [[[fig.fig.add_subplot(subgs[k][i, j]) for j in range(
    subgridspecs[k][1])] for i in range(subgridspecs[k][0])] for
    k in range(len(subgs))]

# Top panel
panel1 = axs[0]
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
# labels = labels[:-2]
pulse_time = 25e3
cmap = matplotlib.cm.get_cmap('viridis')

# Plot figure
plot.add_single(panel1[0][0], CiPA_hERG_log[0], 'membrane.V', color='k')
plot.add_multiple(panel1[1][0], CiPA_hERG_log, 'ikr.IKr', labels=labels,
                  color=cmap)
plot.add_single(panel1[0][1], conductance_hERG_log[0], 'membrane.V', color='k')
plot.add_multiple(panel1[1][1], conductance_hERG_log, 'ikr.IKr', labels=labels,
                  color=cmap)

# panel1[1][0].legend()
panel1[0][0].set_title('State dependent drug block')
panel1[0][1].set_title('Drug block by conductance scaling')
fig.sharex(['Time (s)'] * 2, [(0, pulse_time)] * 2,
           axs=panel1, subgridspec=subgridspecs[0])
fig.sharey(['Voltage\n(mV)', 'Current (A/F)'],
           axs=panel1, subgridspec=subgridspecs[0])
fig.adjust_ticks(panel1[1][0], pulse_time)
fig.adjust_ticks(panel1[1][1], pulse_time)

# Bottom panel
panel2 = axs[1]
# Load action potential data
trapping_data_files = [f for f in os.listdir(saved_data_dir) if
                       f.startswith('CiPA_AP_')]
conductance_data_files = [f for f in os.listdir(saved_data_dir) if
                          f.startswith('conductance_AP_')]
trapping_data_files = [trapping_data_files[i] for i in sort_ind]
conductance_data_files = [conductance_data_files[i] for i in sort_ind]

CiPA_AP_log = []
conductance_AP_log = []
for i in range(len(trapping_data_files)):
    CiPA_AP_log.append(myokit.DataLog.load_csv(
        saved_data_dir + trapping_data_files[i]))
    conductance_AP_log.append(myokit.DataLog.load_csv(
        saved_data_dir + conductance_data_files[i]))

APD_trapping = pd.read_csv(saved_data_dir + 'CiPA_APD_pulses2.csv')
APD_conductance = pd.read_csv(saved_data_dir + 'conductance_APD_pulses2.csv')

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
panel2[0][0].set_title('State dependent drug block')
panel2[0][1].set_title('Drug block by conductance scaling')

unique = fig.legend_without_duplicate_labels(panel2[1][1])
panel2[1][1].legend(*zip(*unique), loc='lower left',
                    bbox_to_anchor=(1.0, 0))
fig.sharex(['Time (ms)'] * 2, [(0, plotting_pulse_time)] * 2,
           axs=panel2, subgridspec=subgridspecs[1])
fig.sharey(['Voltage (mV)', 'Current (A/F)'],
           axs=panel2, subgridspec=subgridspecs[1])

# Add panel letter
fig.fig.set_size_inches(10, 5)
fig.fig.text(0.075, 0.925, '(a)', fontsize=11)
fig.fig.text(0.075, 0.525, '(b)', fontsize=11)

fig.savefig(saved_fig_dir + "hERG_AP_compare.svg", format='svg')
