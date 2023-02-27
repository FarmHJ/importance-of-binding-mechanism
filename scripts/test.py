import myokit
import os
import numpy as np
import pandas as pd

import modelling

drugs = ['dofetilide', 'verapamil']
drug_concs = [30, 1000]  # nM
short_label = ['drug_free', 'dofetilide', 'verapamil']

data_dir = '../simulation_data/background/'
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)

# Simulating the SD model with Milnes' protocol for 10 pulses after addition
# of a dofetilide-like drug and verapamil-like drug

# Load hERG model
model = '../math_model/ohara-cipa-v1-2017-IKr-opt.mmt'
model, _, x = myokit.load(model)
current_model = modelling.BindingKinetics(model)

# Define Milnes' protocol
pulse_time = 25e3
protocol = myokit.Protocol()
protocol.schedule(-80, 0, 800, period=pulse_time)
protocol.schedule(-90, 800, 100, period=pulse_time)
protocol.schedule(-80, 900, 100, period=pulse_time)
protocol.schedule(-80, 11000, 14000, period=pulse_time)
current_model.protocol = protocol
repeats = 10

print(protocol.characteristic_time())

abs_tol = 1e-7
rel_tol = 1e-8

# Simulate control condition


# Simulate 10 pulses after drug addition from steady state of control condition
control_log_single = current_model.drug_simulation(drugs[0], 0, 1000,
                                                   abs_tol=abs_tol,
                                                   rel_tol=rel_tol,
                                                   protocol_period=pulse_time)

drug_label = ['drug free', 'dofetilide', 'verapamil']

fig_dir = '../testing_figures/'
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

expdata_dir = '../exp_data/'
data_dir = '../simulation_data/background/'

# Set up structure of the figure
fig = modelling.figures.FigureStructure(figsize=(9, 6), gridspec=(3, 1))
plot = modelling.figures.FigurePlot()

# Top right panel
# Shows the fractional current of the SD model for first 10 pulses after the
# addition of dofetilide and verapamil
pulse_time = 25e3
repeats = 10

control_log = current_model.drug_simulation(drugs[0], 0, 1000, save_signal=10,
                                            abs_tol=abs_tol, rel_tol=rel_tol,
                                            protocol_period=pulse_time)

for i, drug in enumerate(drugs):

    # Load experimental data
    df = pd.read_csv(expdata_dir + drug + '.csv',
                     header=[0], index_col=[0],
                     skipinitialspace=True)
    current = df.loc[df['conc'] == drug_concs[i]]

    # Load simulated data
    log = current_model.drug_simulation(
        drug, drug_concs[i], repeats, save_signal=repeats,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'],
        set_state=control_log_single, abs_tol=abs_tol, rel_tol=rel_tol,
        protocol_period=pulse_time)

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
        fig.axs[i + 1][0].plot(
            current_exp.index[11:] + (sweep - 1) * pulse_time,
            np.mean(exp_repeats, axis=0), 'o', ms=1.2,
            zorder=-10, color='orange', label=str(drug_concs[i]) + ' nM')
        fig.axs[i + 1][0].fill_between(
            current_exp.index[11:] + (sweep - 1) * pulse_time,
            np.mean(exp_repeats, axis=0) - np.std(exp_repeats, axis=0),
            np.mean(exp_repeats, axis=0) + np.std(exp_repeats, axis=0),
            color='orange', alpha=.3, zorder=-10)

        # Plot the voltage-clamp protocol
        if i == 0:
            fig.axs[i][0].plot(np.array(log.time()) + (sweep - 1) * pulse_time,
                               log['membrane.V', sweep - 1], zorder=1,
                               color='k')

        # Plot simulated data
        fig.axs[i + 1][0].plot(
            np.array(log.time()[log_range_min + 1:log_range_max + 1]) +
            (sweep - 1) * pulse_time,
            np.array(log_plot) / np.array(control_log_plot),
            zorder=1, color='b')

        fig.axs[i + 1][0].set_ylim(0, 1.1)
        fig.axs[i][0].set_rasterization_zorder(0)

# Add description text
fig.axs[1][0].text(248000, 0.9, 'dofetilide', fontsize=8, ha='right')
fig.axs[2][0].text(248000, 0.9, 'verapamil', fontsize=8, ha='right')

# Adjust axes
fig.sharex(['Time (s)'], [(0, pulse_time * repeats)])
fig.sharey(['Voltage\n(mV)', 'Fractional\ncurrent', 'Fractional\ncurrent'])
fig.adjust_ticks(fig.axs[2][0], pulse_time * repeats)

fig.savefig(fig_dir + 'test_Milnes.pdf')
