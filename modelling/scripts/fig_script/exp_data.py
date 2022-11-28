import matplotlib
import matplotlib.pyplot as plt
import myokit
import numpy as np
import pandas as pd

import modelling

data_dir = '../../data/'
drug = 'verapamil'
repeats = 10

df = pd.read_csv(data_dir + drug + '.csv',
                 header=[0], index_col=[0],
                 skipinitialspace=True)

if drug == 'dofetilide':
    drug_conc = [1, 3, 10, 30]
elif drug == 'verapamil':
    drug_conc = [30, 100, 300, 1000]

# Load IKr model
model = '../../model/ohara-cipa-v1-2017-IKr.mmt'
model, _, x = myokit.load(model)

protocol_name = 'Milnes'
protocol_params = modelling.ProtocolParameters()
protocol = protocol_params.protocol_parameters[protocol_name]['function']

drug_model = modelling.BindingKinetics(model)
drug_model.protocol = protocol

t_max = 25e3
protocol3 = myokit.Protocol()
protocol3.schedule(-80, 0, 800, period=t_max)
protocol3.schedule(-90, 800, 100, period=t_max)
protocol3.schedule(-80, 900, 100, period=t_max)
protocol3.schedule(-80, 11000, 14000 - 1, period=t_max)
drug_model.protocol = protocol3

control_log = drug_model.drug_simulation(drug, 0, repeats)
control_log.save_csv(data_dir + 'control_pulses10.csv')
# fig = modelling.figures.FigureStructure(figsize=(10, 6), gridspec=(2, 2),
#                                         height_ratios=[1, 1], hspace=0.3)
# plot = modelling.figures.FigurePlot()
fig = plt.figure(figsize=(5, 3))
cmap = matplotlib.colors.ListedColormap(
    matplotlib.cm.get_cmap('tab10')(range(4)))

for i, conc_i in enumerate(drug_conc):
    print(conc_i)
    current = df.loc[df['conc'] == conc_i]

    log = drug_model.drug_simulation(
        drug, conc_i, repeats,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'],
        set_state=control_log)
    log.save_csv(data_dir + drug + '_conc' + str(conc_i) + '_pulses10.csv')

    max_exp = max(current['exp'].values)
    exp_repeats = []
    for exp in range(1, max_exp):
        current_exp = current.loc[current['exp'] == exp]

        max_sweep = max(current_exp['sweep'].values)
        current_temp = current_exp.loc[current_exp['sweep'] == max_sweep]
        exp_repeats.append(current_temp['frac'].values)

    min_time = min(current_temp.index)
    max_time = max(current_temp.index)
    log_range_min = np.where(log.time() == min_time)[0][0]
    log_range_max = np.where(log.time() == max_time)[0][0]

    log_plot = log['ikr.IKr'][log_range_min:log_range_max + 1]
    control_log_plot = control_log['ikr.IKr'][log_range_min:log_range_max + 1]

    plt.plot(current_temp.index - current_temp.index[0],
             np.mean(exp_repeats, axis=0), 'o', ms=1.2,
             zorder=-10, color=cmap(i), label=str(conc_i) + ' nM')
    plt.fill_between(
        current_temp.index - current_temp.index[0], np.mean(exp_repeats, axis=0) -
        np.std(exp_repeats, axis=0), np.mean(exp_repeats, axis=0) +
        np.std(exp_repeats, axis=0), color=cmap(i), alpha=.3)

    plt.plot(log.time()[log_range_min:log_range_max + 1] - log.time()[log_range_min],
             log_plot / control_log_plot, zorder=1, color='k')
    # fig.axs[int(i / 2)][i % 2].plot(log.time()[log_range_min:log_range_max + 1],
    #                                 log_plot, zorder=1)
    # fig.axs[int(i / 2)][i % 2].plot(log.time()[log_range_min:log_range_max + 1],
    #                                 control_log_plot, zorder=1)
    # print('simulation max: ')
    # print(max(log_plot / control_log_plot))
plt.xlim(min_time  - current_temp.index[0], max_time - current_temp.index[0])
# fig.axs[int(i / 2)][i % 2].set_ylim(0.9 * min(current_temp.index),
#                                     max(current_temp.index))
# plt.title(str(conc_i) + ' nM ' + drug)
plt.xlabel('Time (ms)')
plt.ylabel('Fractional block')
plt.legend()
# fig.savefig('../../data/' + drug + '_expdata_repeats' + str(repeats) + '.pdf')
plt.savefig('../../data/testing.pdf')
plt.close()
