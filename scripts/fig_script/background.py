#
# Figure 1
# Introduces the idea of trapping and justifies the use of the Milnes protocol.
# Plots figure that introduces the SD model and the trapping mechanism.
#

import myokit
import numpy as np
import os
import pandas as pd

import modelling

drugs = ['dofetilide', 'verapamil']
drug_conc = [30, 1000]
drug_label = ['drug free', 'dofetilide', 'verapamil']

fig_dir = '../../figures/background/'
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

expdata_dir = '../../exp_data/'
data_dir = '../../simulation_data/background/'

# Set up structure of the figure
fig = modelling.figures.FigureStructure(figsize=(12, 8), gridspec=(2, 2),
                                        height_ratios=[2, 3], hspace=0.3,
                                        wspace=0.3, plot_in_subgrid=True)
plot = modelling.figures.FigurePlot()

subgridspecs = [(1, 1), (3, 1), (3, 3), (3, 3)]
subgs = []
for i in range(4):
    height_ratio = [1] + [2] * (subgridspecs[i][0] - 1)
    subgs.append(fig.gs[i].subgridspec(*subgridspecs[i], wspace=0.1,
                                       hspace=0.05,
                                       height_ratios=height_ratio))
axs = [[[fig.fig.add_subplot(subgs[k + 1][i, j]) for j in range(
    subgridspecs[k + 1][1])] for i in range(subgridspecs[k + 1][0])] for
    k in range(len(subgs) - 1)]

# Top right panel
# Shows the fractional current of the SD model for first 10 pulses after the
# addition of dofetilide and verapamil
panel2 = axs[0]

pulse_time = 25e3
repeats = 10

control_log = myokit.DataLog.load_csv(data_dir +
                                      'control_Milnes_current_pulses10.csv')

for i, drug in enumerate(drugs):

    # Load experimental data
    df = pd.read_csv(expdata_dir + drug + '.csv',
                     header=[0], index_col=[0],
                     skipinitialspace=True)
    current = df.loc[df['conc'] == drug_conc[i]]

    # Load simulated data
    log = myokit.DataLog.load_csv(
        data_dir + drug + '_Milnes_current_pulses10.csv')

    max_sweeps = max(current['sweep'].values)
    for sweep in range(1, max_sweeps + 1):

        current_sweep = current.loc[current['sweep'] == sweep]
        max_exp = max(current['exp'].values)
        exp_repeats = []
        for exp in range(1, max_exp):
            current_exp = current_sweep.loc[current_sweep['exp'] == exp]
            exp_repeats.append(current_exp['frac'].values[11:])

        min_time = min(current_exp.index[11:])
        max_time = max(current_exp.index[11:])
        log_range_min = np.argmin(np.abs(np.array(log.time()) - min_time))
        log_range_max = np.argmin(np.abs(np.array(log.time()) - max_time))

        # Organise simulated data
        log_plot = log['ikr.IKr', sweep - 1][
            log_range_min + 1:log_range_max + 1]
        control_log_plot = control_log['ikr.IKr', sweep - 1][
            log_range_min + 1:log_range_max + 1]

        # Plot experimental data
        panel2[i + 1][0].plot(
            current_exp.index[11:] + (sweep - 1) * pulse_time,
            np.mean(exp_repeats, axis=0), 'o', ms=1.2,
            zorder=-10, color='orange', label=str(drug_conc[i]) + ' nM')
        panel2[i + 1][0].fill_between(
            current_exp.index[11:] + (sweep - 1) * pulse_time,
            np.mean(exp_repeats, axis=0) - np.std(exp_repeats, axis=0),
            np.mean(exp_repeats, axis=0) + np.std(exp_repeats, axis=0),
            color='orange', alpha=.3, zorder=-10)

        # Plot the voltage-clamp protocol
        if i == 0:
            panel2[i][0].plot(np.array(log.time()) + (sweep - 1) * pulse_time,
                              log['membrane.V', sweep - 1], zorder=1,
                              color='k')

        # Plot simulated data
        panel2[i + 1][0].plot(
            np.array(log.time()[log_range_min + 1:log_range_max + 1]) +
            (sweep - 1) * pulse_time,
            np.array(log_plot) / np.array(control_log_plot),
            zorder=1, color='b')

        panel2[i + 1][0].set_ylim(0, 1.1)
        panel2[i][0].set_rasterization_zorder(0)

# Add description text
panel2[1][0].text(248000, 0.9, 'dofetilide', fontsize=8, ha='right')
panel2[2][0].text(248000, 0.9, 'verapamil', fontsize=8, ha='right')

# Adjust axes
fig.sharex(['Time (s)'], [(0, pulse_time * repeats)],
           axs=panel2, subgridspec=(3, 1))
fig.sharey(['Voltage\n(mV)', 'Fractional\ncurrent', 'Fractional\ncurrent'],
           axs=panel2, subgridspec=(3, 1))
fig.adjust_ticks(panel2[2][0], pulse_time * repeats)

# Bottom left panel
panel3 = axs[1]

# Load steady state IKr data for drug free, addition of dofetilide
# and verapamil conditions (Milnes' protocol)
drug_free_log = myokit.DataLog.load_csv(
    data_dir + 'drug_free_Milnes_current.csv')
trapped_log = myokit.DataLog.load_csv(
    data_dir + 'dofetilide_Milnes_current.csv')
nontrapped_log = myokit.DataLog.load_csv(
    data_dir + 'verapamil_Milnes_current.csv')
log_all = [drug_free_log, trapped_log, nontrapped_log]

# Load hERG channel model
model = '../../math_model/ohara-cipa-v1-2017-IKr-opt.mmt'
model, _, x = myokit.load(model)

# Plot state occupancies of the hERG channel
color_seq = ['#7e7e7e', '#986464', '#989864', '#986496', '#988364',
             '#64986a', '#74a9cf', '#045a8d', '#2b8cbe']

for d in range(len(drugs) + 1):
    plot.add_single(panel3[0][d], log_all[d], 'membrane.V')
    if d == 2:
        plot.state_occupancy_plot(panel3[2][d], log_all[d], model,
                                  color_seq=color_seq)
    else:
        plot.state_occupancy_plot(panel3[2][d], log_all[d],
                                  model, legend=False, color_seq=color_seq)

# Plot IKr
for i in range(len(log_all)):
    for j in range(len(log_all)):
        if i == j:
            plot.add_single(panel3[1][i], log_all[j], 'ikr.IKr')
        else:
            plot.add_single(panel3[1][i], log_all[j], 'ikr.IKr',
                            color='grey', alpha=0.5)
    panel3[1][i].text(24500, 0.8, drug_label[i], fontsize=8,
                      ha='right', va='top')

# Adjust axes
for col in range(3):
    for tick in panel3[2][col].get_xticklabels():
        tick.set_ha('right')
fig.sharex(['Time (s)'] * (len(drugs) + 1),
           [(0, pulse_time)] * (len(drugs) + 1),
           axs=panel3, subgridspec=subgridspecs[3])
fig.sharey(['Voltage\n(mV)', 'Current (A/F)', 'State\noccupancy'],
           axs=panel3, subgridspec=subgridspecs[3])
for i in range(len(log_all)):
    fig.adjust_ticks(panel3[2][i], pulse_time)

# Bottom right panel
panel4 = axs[2]

# Load steady state IKr data for drug free, addition of a dofetilide-like drug
# and a verapamil-like drug conditions (AP clamp protocol)
drug_free_log = myokit.DataLog.load_csv(
    data_dir + 'drug_free_APclamp_current.csv')
trapped_log = myokit.DataLog.load_csv(
    data_dir + 'dofetilide_APclamp_current.csv')
nontrapped_log = myokit.DataLog.load_csv(
    data_dir + 'verapamil_APclamp_current.csv')
log_all = [drug_free_log, trapped_log, nontrapped_log]

# Load AP model
APmodel = '../../math_model/ohara-cipa-v1-2017-opt.mmt'
APmodel, _, x = myokit.load(APmodel)
pulse_time = 1000

# Plot state occupancies of the hERG channel
for d in range(len(drugs) + 1):
    plot.state_occupancy_plot(panel4[2][d], log_all[d],
                              APmodel, legend=False, color_seq=color_seq)

for i in range(len(log_all)):
    for j in range(len(log_all)):
        if i == j:
            plot.add_single(panel4[0][i], log_all[j], 'membrane.V')
            plot.add_single(panel4[1][i], log_all[j], 'ikr.IKr')
        else:
            plot.add_single(panel4[1][i], log_all[j], 'ikr.IKr',
                            color='grey', alpha=0.5)
    panel4[1][i].text(980, 0.9, drug_label[i], fontsize=8,
                      ha='right', va='top')

# Adjust axes
for col in range(3):
    for tick in panel4[2][col].get_xticklabels():
        tick.set_ha('right')
fig.sharex(['Time (ms)'] * (len(drugs) + 1),
           [(0, pulse_time)] * (len(drugs) + 1),
           axs=panel4, subgridspec=subgridspecs[2])
fig.sharey(['Voltage\n(mV)', 'Current\n(A/F)', 'State\noccupancy'],
           axs=panel4, subgridspec=subgridspecs[2])

# Add panel letter
fig.fig.set_size_inches(10, 6.5)
fig.fig.text(0.075, 0.925, '(A)', fontsize=11)
fig.fig.text(0.5, 0.925, '(B)', fontsize=11)
fig.fig.text(0.075, 0.525, '(C)', fontsize=11)
fig.fig.text(0.5, 0.525, '(D)', fontsize=11)

fig.savefig(fig_dir + "background.svg", format='svg')
