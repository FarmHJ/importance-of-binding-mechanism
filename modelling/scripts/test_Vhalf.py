# Introduces the idea of trapping and justifies the use of the Milnes protocol
import myokit
import numpy as np

import modelling


drugs = ['dofetilide', 'verapamil']
drug_concs = [30, 1000]  # nM

testing_fig_dir = '../../figures/testing/'
saved_fig_dir = testing_fig_dir

saved_data_dir = '../../simulation_data/background/'

model = '../../model/ohara-cipa-v1-2017-IKr.mmt'
model, _, x = myokit.load(model)

drug_model = modelling.BindingKinetics(model)
param_lib = modelling.BindingParameters()
pulse_time = 5400
repeats = 1000


def testing_protocol(Vhalf, pulse_time):
    protocol = myokit.Protocol()
    if Vhalf > -50:
        protocol.schedule(-80, 0, 100, period=pulse_time)
        protocol.schedule(Vhalf, 100, 5000, period=pulse_time)
        protocol.schedule(-60, 5100, 200, period=pulse_time)
        protocol.schedule(-80, 5300, 100, period=pulse_time)
    else:
        protocol.schedule(-80, 0, 100, period=pulse_time)
        protocol.schedule(20, 100, 500, period=pulse_time)
        protocol.schedule(-50, 600, 200, period=pulse_time)
        protocol.schedule(Vhalf, 800, 4500, period=pulse_time)
        protocol.schedule(-80, 5300, 100, period=pulse_time)

    return protocol


fig = modelling.figures.FigureStructure(figsize=(8, 4), gridspec=(3, 2),
                                        wspace=0.05)
plot = modelling.figures.FigurePlot()

for d in range(len(drugs)):

    Vhalf = param_lib.binding_parameters[drugs[d]]['Vhalf']

    print(Vhalf)
    protocol = testing_protocol(Vhalf, pulse_time)
    drug_model.protocol = protocol

    log = drug_model.drug_simulation(drugs[d], drug_concs[d], repeats)

    plot.add_single(fig.axs[0][d], log, 'membrane.V')
    plot.add_single(fig.axs[1][d], log, 'ikr.IKr')
    plot.state_occupancy_plot(fig.axs[2][d], log, model, legend=False)

fig.sharex(['Time (s)'] * (len(drugs)),
            [(0, pulse_time)] * (len(drugs)))
fig.sharey(['Voltage (mV)', 'Current (A/F)', 'State occupancy'])
for i in range(len(drugs)):
    fig.axs[0][i].set_title(drugs[i])
    fig.adjust_ticks(fig.axs[2][i], pulse_time)
fig.savefig(saved_fig_dir + "hERG_half_occupancy_Vhalf.pdf")
