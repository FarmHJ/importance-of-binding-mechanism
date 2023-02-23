# Introduces the idea of trapping and justifies the use of the Milnes protocol
import myokit
import os

import modelling

drugs = ['cisapride']
drug_concs = [1000]  # nM
short_label = ['cisapride']

data_dir = '../../simulation_data/background/'
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)

# Simulating the SD model with Milnes' protocol for 10 pulses after addition
# of a dofetilide-like drug and verapamil-like drug

# Load hERG model
model = '../../math_model/ohara-cipa-v1-2017-IKr-opt.mmt'
model, _, x = myokit.load(model)
current_model = modelling.BindingKinetics(model)

# Define Milnes' protocol
pulse_time = 25e3
protocol = myokit.Protocol()
protocol.schedule(-80, 0, 800, period=pulse_time)
protocol.schedule(-90, 800, 100, period=pulse_time)
protocol.schedule(-80, 900, 100, period=pulse_time)
protocol.schedule(-80, 11000, 14000 - 1, period=pulse_time)
current_model.protocol = protocol
repeats = 10

abs_tol = 1e-7
rel_tol = 1e-8

# # Simulate control condition
# control_log = current_model.drug_simulation(drugs[0], 0, 1000, save_signal=10,
#                                             abs_tol=abs_tol, rel_tol=rel_tol)
# control_log.save_csv(data_dir + 'control_Milnes_current_pulses10.csv')

# # Simulate 10 pulses after drug addition from steady state of control condition
# control_log_single = current_model.drug_simulation(drugs[0], 0, 1000)
# for i, conc_i in enumerate(drug_concs):

#     log = current_model.drug_simulation(
#         drugs[i], conc_i, repeats, save_signal=repeats,
#         log_var=['engine.time', 'membrane.V', 'ikr.IKr'],
#         set_state=control_log_single, abs_tol=abs_tol, rel_tol=rel_tol)
#     log.save_csv(data_dir + drugs[i] + '_Milnes_current_pulses10.csv')

# Simulating the SD model with Milnes' protocol for 1000 pulses till steady
# state under drug free, addition of dofetilide-like drug and verapamil-like
# drug conditions

# Define Milnes' protocol
Milnes_protocol = modelling.ProtocolLibrary().Milnes(pulse_time)
current_model.protocol = Milnes_protocol

# Simulate SD model under drug-free condition
repeats = 1000
# log_control = current_model.drug_simulation(drugs[0], 0, repeats,
#                                             abs_tol=abs_tol, rel_tol=rel_tol)
# log_control.save_csv(data_dir + short_label[0] + '_Milnes_current.csv')

# Simulate the SD model under addition of drugs
# for d in range(len(drugs)):
d = 0
log = current_model.drug_simulation(drugs[d], drug_concs[d], repeats,
                                    abs_tol=abs_tol, rel_tol=rel_tol)
log.save_csv(data_dir + short_label[d] + '_Milnes_current.csv')

fig = modelling.figures.FigureStructure(figsize=(5, 5), gridspec=(2, 1),
                                        height_ratios=[1, 1], hspace=0.1)
plot = modelling.figures.FigurePlot()
color_seq = ['#7e7e7e', '#986464', '#989864', '#986496', '#988364',
             '#64986a', '#74a9cf', '#045a8d', '#2b8cbe']

plot.state_occupancy_plot(fig.axs[1][0], log, model, color_seq=color_seq)
plot.add_single(fig.axs[0][0], log, 'ikr.IKr')

# Adjust axes
fig.sharex(['Time (s)'], [(0, pulse_time)])
fig.sharey(['Current (A/F)', 'State\noccupancy'])
fig.adjust_ticks(fig.axs[1][0], pulse_time)

fig.savefig('../../testing_figures/cisapride_current.pdf')

# Simulating the SD model with AP clamp protocol for 1000 pulses till steady
# state under drug free, addition of dofetilide-like drug and verapamil-like
# drug conditions

# Load AP model
APmodel = '../../math_model/ohara-cipa-v1-2017-opt.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')

# Define current pulse
pulse_time = 1000
protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)
AP_model.protocol = protocol

# Simulate AP for AP clamp protocol
# APclamp = AP_model.drug_simulation(drugs[0], drug_concs[0], repeats,
#                                    abs_tol=abs_tol, rel_tol=rel_tol)
# APclamp.save_csv(data_dir + 'APclamp.csv')
APclamp = myokit.DataLog.load_csv(data_dir + 'APclamp.csv')

# Set up AP clamp protocol
times = APclamp['engine.time']
voltages = APclamp['membrane.V']
tmax = times[-1] + 1

# Simulate SD model with AP clamp protocol under drug free condition
# log_control = current_model.drug_APclamp(drugs[0], 0, times, voltages,
#                                          tmax, repeats, abs_tol=abs_tol,
#                                          rel_tol=rel_tol)
# log_control.save_csv(data_dir + short_label[0] + '_APclamp_current.csv')

# Simulate the SD model under addition of drugs
# for d in range(len(drugs)):
log = current_model.drug_APclamp(drugs[d], drug_concs[d], times,
                                 voltages, tmax, repeats, abs_tol=abs_tol,
                                 rel_tol=rel_tol)
log.save_csv(data_dir + short_label[d] + '_APclamp_current.csv')

fig = modelling.figures.FigureStructure(figsize=(5, 6), gridspec=(3, 1),
                                        height_ratios=[1, 1, 1], hspace=0.1)

plot.state_occupancy_plot(fig.axs[2][0], log, model, color_seq=color_seq)
plot.add_single(fig.axs[1][0], log, 'ikr.IKr')
plot.add_single(fig.axs[0][0], log, 'membrane.V')

# Adjust axes
fig.sharex(['Time (s)'], [(0, pulse_time)])
fig.sharey(['Voltage\n(mV)', 'Current (A/F)', 'State\noccupancy'])
fig.adjust_ticks(fig.axs[2][0], pulse_time)

fig.savefig('../../testing_figures/cisapride_AP.pdf')
