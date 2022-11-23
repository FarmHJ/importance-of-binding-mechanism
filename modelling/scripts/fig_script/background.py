# Introduces the idea of trapping and justifies the use of the Milnes protocol
import myokit
import numpy as np

import modelling

steady_state = False
hERG_model = True

drugs = ['dofetilide', 'verapamil']
drug_label = ['drug free', 'dofetilide-like\ndrug', 'verapamil-like\ndrug']

testing_fig_dir = '../../figures/testing/'
final_fig_dir = '../../figures/background/trapping/'
# final_fig_dir = '../../figures/conferences/'
saved_fig_dir = final_fig_dir

saved_data_dir = '../../simulation_data/background/'

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
# subgs.append(fig.gs[-1].subgridspec(*subgridspecs[-1], wspace=0.1,
#                                     hspace=0.05,
#                                     height_ratios=[1] * subgridspecs[-1][0]))
axs = [[[fig.fig.add_subplot(subgs[k + 1][i, j]) for j in range(
    subgridspecs[k + 1][1])] for i in range(subgridspecs[k + 1][0])] for
    k in range(len(subgs) - 1)]

# Top right panel
trapped_log = myokit.DataLog.load_csv(
    saved_data_dir + 'trapped_current_transient.csv',
    precision=myokit.DOUBLE_PRECISION)
nontrapped_log = myokit.DataLog.load_csv(
    saved_data_dir + 'nontrapped_current_transient.csv')

pulse_time = 25e3
repeats = 10
max_hERG = np.max(trapped_log['ikr.IKr', 0])

length = int(len(trapped_log.time()) / 2)
panel2 = axs[0]
plot.add_continuous(panel2[0][0], trapped_log, 'membrane.V')
plot.add_continuous(panel2[1][0], trapped_log, 'ikr.IKr')
plot.add_continuous(panel2[2][0], nontrapped_log, 'ikr.IKr')

# Add on to cut off section
panel2[0][0].plot(
    np.array(trapped_log.time()) + max(trapped_log.time()) / 2,
    trapped_log['membrane.V', 0][length:] + trapped_log[
        'membrane.V', 1][:length], 'k')
panel2[1][0].plot(
    np.array(trapped_log.time()) + max(trapped_log.time()) / 2,
    trapped_log['ikr.IKr', 0][length:] + trapped_log['ikr.IKr', 1][:length],
    'k')
panel2[2][0].plot(
    np.array(trapped_log.time()) + max(trapped_log.time()) / 2,
    nontrapped_log['ikr.IKr', 0][length:] + nontrapped_log[
        'ikr.IKr', 1][:length], 'k')

panel2[1][0].text(248000, 0.6, 'dofetilide-like drug', fontsize=8, ha='right')
panel2[2][0].text(248000, 0.6, 'verapamil-like drug', fontsize=8, ha='right')

fig.sharex(['Time (s)'], [(0, pulse_time * repeats)],
           axs=panel2, subgridspec=(3, 1))
fig.sharey(['Voltage\n(mV)', 'Current\n(A/F)', 'Current\n(A/F)'],
           [None, (0, max_hERG + 0.05), (0, max_hERG + 0.05)],
           axs=panel2, subgridspec=(3, 1))
fig.adjust_ticks(panel2[2][0], pulse_time * repeats)

# Bottom left panel
drug_free_log = myokit.DataLog.load_csv(
    saved_data_dir + 'drug_free_current.csv')
trapped_log = myokit.DataLog.load_csv(
    saved_data_dir + 'trapped_current.csv')
nontrapped_log = myokit.DataLog.load_csv(
    saved_data_dir + 'nontrapped_current.csv')
log_all = [drug_free_log, trapped_log, nontrapped_log]

model = '../../model/ohara-cipa-v1-2017-IKr.mmt'
model, _, x = myokit.load(model)
pulse_time = 25e3

panel3 = axs[1]
for d in range(len(drugs) + 1):
    plot.add_single(panel3[0][d], log_all[d], 'membrane.V')
    if d == 2:
        plot.state_occupancy_plot(panel3[2][d], log_all[d], model)
    else:
        plot.state_occupancy_plot(panel3[2][d], log_all[d],
                                  model, legend=False)

for i in range(len(log_all)):
    for j in range(len(log_all)):
        if i == j:
            plot.add_single(panel3[1][i], log_all[j], 'ikr.IKr')
        else:
            plot.add_single(panel3[1][i], log_all[j], 'ikr.IKr',
                            color='grey', alpha=0.5)
    panel3[1][i].text(24500, 0.75, drug_label[i], fontsize=8,
                      ha='right', va='top')

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
drug_free_log = myokit.DataLog.load_csv(
    saved_data_dir + 'drug_free_APclamp_current.csv')
trapped_log = myokit.DataLog.load_csv(
    saved_data_dir + 'dofetilide_APclamp_current.csv')
nontrapped_log = myokit.DataLog.load_csv(
    saved_data_dir + 'verapamil_APclamp_current.csv')
log_all = [drug_free_log, trapped_log, nontrapped_log]

APmodel = '../../model/ohara-cipa-v1-2017.mmt'
APmodel, _, x = myokit.load(APmodel)
pulse_time = 1000

panel4 = axs[2]
for d in range(len(drugs) + 1):
    if d == 2:
        plot.state_occupancy_plot(panel4[2][d], log_all[d], APmodel)
    else:
        plot.state_occupancy_plot(panel4[2][d], log_all[d],
                                  APmodel, legend=False)

for i in range(len(log_all)):
    for j in range(len(log_all)):
        if i == j:
            plot.add_single(panel4[0][i], log_all[j], 'membrane.V')
            plot.add_single(panel4[1][i], log_all[j], 'ikr.IKr')
        else:
            plot.add_single(panel4[0][i], log_all[j], 'membrane.V',
                            color='grey', alpha=0.5)
            plot.add_single(panel4[1][i], log_all[j], 'ikr.IKr',
                            color='grey', alpha=0.5)
    panel4[1][i].text(980, 0.8, drug_label[i], fontsize=8,
                      ha='right', va='top')

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
# fig.fig.set_size_inches(20, 8)
fig.fig.text(0.075, 0.925, '(A)', fontsize=11)
fig.fig.text(0.5, 0.925, '(B)', fontsize=11)
fig.fig.text(0.075, 0.525, '(C)', fontsize=11)
fig.fig.text(0.5, 0.525, '(D)', fontsize=11)

fig.savefig(saved_fig_dir + "background_compile.svg", format='svg')
