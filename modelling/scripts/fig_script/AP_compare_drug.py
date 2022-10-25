# Compares the hERG current and action potential of state dependent drug block
# and conductance scaling
import pandas as pd
import matplotlib
import myokit
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

# panelC_title = 'Verapamil-like drug'

fig = modelling.figures.FigureStructure(figsize=(10, 5.5), gridspec=(2, 1),
                                        height_ratios=[1, 1], hspace=0.45,
                                        plot_in_subgrid=True)
plot = modelling.figures.FigurePlot()

subgridspecs = [(2, 2), (2, 2)]
subgs = []
subgs.append(fig.gs[0].subgridspec(*subgridspecs[0], wspace=0.3,
                                   hspace=0.4))
subgs.append(fig.gs[1].subgridspec(*subgridspecs[1], wspace=0.05,
                                   hspace=0.05))
axs = [[[fig.fig.add_subplot(subgs[0][i, 0])] for
        i in range(subgridspecs[0][0])]]
axs.append([[fig.fig.add_subplot(subgs[0][:, 1])]])
axs.append([[fig.fig.add_subplot(subgs[1][i, j]) for j in range(
    subgridspecs[1][1])] for i in range(subgridspecs[1][0])])

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
print(trapping_data_files)
# for verapamil
# removing_ind = drug_conc.index(500.0)
# drug_conc.pop(removing_ind)
# trapping_data_files.pop(removing_ind)
# conductance_data_files.pop(removing_ind)
# conc_label.pop(removing_ind)
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

print(APD_trapping)
print(APD_conductance)
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
panel2[0][0].set_title('State dependent drug block')
panel2[0][1].set_title('Conductance scaling drug block')

unique = fig.legend_without_duplicate_labels(panel2[1][1])
panel2[1][1].legend(*zip(*unique), loc='lower left',
                    bbox_to_anchor=(1.0, 0), handlelength=1)
fig.sharex(['Time (ms)'] * 2, [(0, plotting_pulse_time)] * 2,
           axs=panel2, subgridspec=subgridspecs[1])
fig.sharey(['Voltage (mV)', 'Current (A/F)'],
           axs=panel2, subgridspec=subgridspecs[1])

# Top left panel
panel1 = axs[0]
# Load hERG current data
trapping_data_files = [f for f in os.listdir(saved_data_dir) if
                       f.startswith('CiPA_hERG_')]
conductance_data_files = [f for f in os.listdir(saved_data_dir) if
                          f.startswith('conductance_hERG_')]
drug_conc = [float(fname[10:-4]) for fname in trapping_data_files]
print(trapping_data_files)
# for verapamil
# removing_ind = drug_conc.index(500.0)
# trapping_data_files.pop(removing_ind)
# conductance_data_files.pop(removing_ind)

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
plot.add_multiple(panel1[0][0], hERG_trapping_plot, 'ikr.IKr', labels=labels,
                  color=cmap)
plot.add_multiple(panel1[1][0], hERG_conductance_plot, 'ikr.IKr',
                  labels=labels, color=cmap)

panel1[1][0].legend()
panel1[0][0].set_title('State dependent drug block')
panel1[1][0].set_title('Conductance scaling drug block')
fig.sharex(['Time (s)'], [(0, pulse_time)],
           axs=panel1, subgridspec=(2, 1))
fig.sharey(['Current (A/F)'] * 2,
           axs=panel1, subgridspec=(2, 1))
for i in range(2):
    panel1[i][0].spines['top'].set_visible(False)
    panel1[i][0].spines['right'].set_visible(False)
    fig.adjust_ticks(panel1[i][0], pulse_time)

# Top right panel
panel3 = axs[1]

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
print(drug_conc)
print(APD_trapping)
print(APD_conductance)

panel3[0][0].plot(drug_conc[1:], APD_trapping[1:], 'o', color='orange',
                  label='SD model')
panel3[0][0].plot(drug_conc[1:], APD_conductance[1:], '^', color='blue',
                  label='CS model', alpha=0.8)
panel3[0][0].plot(drug_conc[1:], EAD_marker, 'o', color='k',
                  marker=(5, 2), label='EAD-like AP')
panel3[0][0].set_xscale("log", nonpositive='clip')
panel3[0][0].set_xlabel('Drug concentration (nM)')
panel3[0][0].set_ylabel(r'APD$_{90}$ (ms)')
# panel3[0][0].set_title(panelC_title)
panel3[0][0].legend(handlelength=1)

# Add panel letter
fig.fig.set_size_inches(10, 5.5)
fig.fig.text(0.075, 0.925, '(A)', fontsize=11)
fig.fig.text(0.075, 0.455, '(B)', fontsize=11)
fig.fig.text(0.525, 0.925, '(C)', fontsize=11)

# fig.savefig(saved_fig_dir + "model_compare.svg", format='svg')
