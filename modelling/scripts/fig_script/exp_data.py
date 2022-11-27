import matplotlib.pyplot as plt
import myokit
import numpy as np
import pandas as pd

import modelling

data_dir = '../../data/'
drug = 'dofetilide'
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

fig = modelling.figures.FigureStructure(figsize=(10, 6), gridspec=(2, 2),
                                        height_ratios=[1, 1], hspace=0.3)
plot = modelling.figures.FigurePlot()

for i, conc_i in enumerate(drug_conc):
    print(conc_i)
    current = df.loc[df['conc'] == conc_i]

    log = drug_model.drug_simulation(
        drug, conc_i, repeats,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'],
        set_state=control_log)

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

    fig.axs[int(i / 2)][i % 2].plot(current_temp.index,
                                    np.mean(exp_repeats, axis=0), '-',
                                    zorder=-10)
    fig.axs[int(i / 2)][i % 2].fill_between(
        current_temp.index, np.mean(exp_repeats, axis=0) -
        np.std(exp_repeats, axis=0), np.mean(exp_repeats, axis=0) +
        np.std(exp_repeats, axis=0), color='b', alpha=.3)

    fig.axs[int(i / 2)][i % 2].plot(log.time()[log_range_min:log_range_max + 1],
                                    log_plot / control_log_plot,
                                    zorder=1)
    # fig.axs[int(i / 2)][i % 2].plot(log.time()[log_range_min:log_range_max + 1],
    #                                 log_plot, zorder=1)
    # fig.axs[int(i / 2)][i % 2].plot(log.time()[log_range_min:log_range_max + 1],
    #                                 control_log_plot, zorder=1)
    # print('simulation max: ')
    # print(max(log_plot / control_log_plot))
    fig.axs[int(i / 2)][i % 2].set_xlim(min_time, max_time)
    # fig.axs[int(i / 2)][i % 2].set_ylim(0.9 * min(current_temp.index),
    #                                     max(current_temp.index))
    fig.axs[int(i / 2)][i % 2].set_title(str(conc_i) + ' nM ' + drug)
    fig.axs[int(i / 2)][i % 2].set_xlabel('Time (ms)')
    fig.axs[int(i / 2)][i % 2].set_ylabel('Fractional block')
# fig.savefig('../../data/' + drug + '_expdata_repeats' + str(repeats) + '.pdf')
fig.savefig('../../data/testing.pdf')
plt.close()
