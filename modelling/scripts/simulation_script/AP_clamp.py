import matplotlib.pyplot as plt
import myokit
import numpy as np

import modelling

model = '../../model/ohara-cipa-v1-2017-IKr.mmt'
model, _, x = myokit.load(model)

# Load action potential
saved_data_dir = '../../simulation_data/background/'
APlog = myokit.DataLog.load_csv(saved_data_dir + 'nontrapped_AP.csv')
APlog2 = myokit.DataLog.load_csv(saved_data_dir + 'nontrapped_AP.csv')

APlog = APlog.npview()
APlog2 = APlog2.npview()

times = APlog['engine.time']
voltages = APlog['membrane.V']

times2 = APlog2['engine.time']
voltages2 = APlog2['membrane.V']

# plt.figure(figsize=(10, 5))
# plt.xlabel('Time (ms)')
# plt.ylabel('Voltage (mV)')
# plt.plot(times, voltages, label='trapped')
# plt.plot(times2, voltages2, label='nontrapped')
# plt.legend()
# plt.show()

sim = myokit.Simulation(model)
sim.set_fixed_form_protocol(times, voltages)

current_model = modelling.BindingKinetics(model)

tmax = APlog['engine.time'][-1] + 1
sim.pre(tmax * 1000)
log_temp = current_model.sim.run(tmax, log_times=times)
# log_temp = current_model.drug_APclamp('verapamil', 0, times, voltages, tmax, 1000)

plt.figure(figsize=(10, 6))

plt.subplot(2, 1, 1)
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')
plt.plot(times, voltages, label='CSV data')
plt.plot(log_temp['engine.time'], log_temp['membrane.V'], label='Simulation')
plt.legend()

plt.subplot(2, 1, 2)
plt.xlabel('Time (ms)')
plt.ylabel('Current (nA)')
plt.plot(log_temp['engine.time'], log_temp['ikr.IKr'], label='Simulation')
plt.legend()

plt.show()

fig = modelling.figures.FigureStructure(figsize=(8, 4), gridspec=(3, 3),
                                        wspace=0.05)
plot = modelling.figures.FigurePlot()

drugs = ['dofetilide', 'verapamil']
# drugs = ['verapamil', 'dofetilide']
short_label = ['drug_free', 'dofetilide', 'verapamil']
drug_label = ['drug free', 'dofetilide-like drug', 'verapamil-like drug']
drug_concs = [30, 1000]
repeats = 1

log_all = []
log_control = current_model.drug_APclamp(drugs[0], 0, times, voltages,
                                         tmax, repeats)
log_control.save_csv(
    saved_data_dir + short_label[0] + '_APclamp_current.csv')
log_all.append(log_control)

plot.add_single(fig.axs[0][0], log_control, 'membrane.V')
plot.state_occupancy_plot(fig.axs[2][0], log_control, model)

for d in range(len(drugs)):
    log = current_model.drug_APclamp(drugs[d], drug_concs[d], times, voltages, tmax, repeats)
    log_all.append(log)

    log.save_csv(
        saved_data_dir + short_label[d + 1] + '_APclamp_current.csv')

    plot.add_single(fig.axs[0][d + 1], log, 'membrane.V')
    plot.state_occupancy_plot(fig.axs[2][d + 1], log, model, legend=False)

for i in range(len(log_all)):
    for j in range(len(log_all)):
        if i == j:
            plot.add_single(fig.axs[1][i], log_all[j], 'ikr.IKr')
        else:
            plot.add_single(fig.axs[1][i], log_all[j], 'ikr.IKr',
                            color='grey', alpha=0.5)
    fig.axs[1][i].text(970, 0.8, drug_label[i], fontsize=8, ha='right')

fig.sharex(['Time (ms)'] * (len(drugs) + 1),
           [(0, tmax)] * (len(drugs) + 1))
fig.sharey(['Voltage (mV)', 'Current (A/F)', 'State occupancy'])
# for i in range(len(log_all)):
#     fig.adjust_ticks(fig.axs[2][i], tmax)
# fig.fig.show()
saved_fig_dir = '../../figures/testing/'
fig.savefig(saved_fig_dir + "hERG_drug_state_occupancy.pdf")
