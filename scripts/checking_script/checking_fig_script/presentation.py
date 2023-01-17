# Compares the hERG current and action potential of state dependen drug block
# and conductance scaling
import pandas as pd
import matplotlib
import myokit
import numpy as np
import os

import modelling

data_dir = '../../simulation_data/background/'
fig_dir = '../../testing_figures/presentation/'

# Load hERG channel model
model = '../../math_model/ohara-cipa-v1-2017-IKr-opt.mmt'
model, _, x = myokit.load(model)

# Plot state occupancies of the hERG channel
color_seq = ['#7e7e7e', '#986464', '#989864', '#986496', '#988364',
             '#64986a', '#74a9cf', '#045a8d', '#2b8cbe']

#########################
#
# Simulate and plot AP of CS model with ion channel currents
#
#########################
# Load AP model
APmodel = '../../math_model/ohara-cipa-v1-2017-opt.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')

# Define current pulse
pulse_time = 1000
protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)
AP_model.protocol = protocol

repeats = 1000
APlog = AP_model.drug_simulation('dofetilide', 0, repeats)

signal = APlog['membrane.V']
APA = max(signal) - min(signal)
APD90_v = min(signal) + 0.1 * APA
index = np.abs(np.array(signal) - APD90_v).argmin()

# Plot AP of AP-CS model
fig = modelling.figures.FigureStructure(figsize=(6, 7), gridspec=(2, 1),
                                        height_ratios=[2, 4], hspace=0.3,
                                        plot_in_subgrid=True)
plot = modelling.figures.FigurePlot()

subgridspecs = [(1, 1), (3, 2)]
subgs = []
for i in range(2):
    subgs.append(fig.gs[i].subgridspec(*subgridspecs[i], wspace=0.25,
                                       hspace=0.3))
axs = [[[fig.fig.add_subplot(subgs[k][i, j]) for j in range(
    subgridspecs[k][1])] for i in range(subgridspecs[k][0])] for
    k in range(len(subgs))]

top_panel = axs[0]
plot.add_single(top_panel[0][0], APlog, 'membrane.V')
top_panel[0][0].axhline(y=APD90_v, xmin=50 / 500, xmax=index / 5000,
                        color='red')
top_panel[0][0].text(150, APD90_v + 1, r'APD$_{90}$', fontsize=10,
                     ha='left', va='bottom', color='red')

top_panel[0][0].set_xlim(0, 500)
top_panel[0][0].set_xlabel('Time (ms)')
top_panel[0][0].set_ylabel(r'$V_m$ (mV)')
top_panel[0][0].set_title('Action potential')

top_panel[0][0].spines['right'].set_visible(False)
top_panel[0][0].spines['top'].set_visible(False)

current_list = ['inal.INaL', 'ical.ICaL', 'ikr.IKr', 'iks.IKs', 'ik1.IK1',
                'ito.Ito']
current_title = [r'$I_\mathrm{NaL}$', r'$I_\mathrm{CaL}$', r'$I_\mathrm{Kr}$',
                 r'$I_\mathrm{Ks}$', r'$I_\mathrm{K1}$', r'$I_\mathrm{to}$']

bottom_panel = axs[1]
count = 0
for i in range(len(current_list)):
    plot.add_single(bottom_panel[int(i / 2)][i % 2], APlog, current_list[i])
    bottom_panel[int(i / 2)][i % 2].set_title(current_title[i])
    bottom_panel[int(i / 2)][i % 2].set_xlim(0, 500)

    bottom_panel[int(i / 2)][i % 2].spines['right'].set_visible(False)
    bottom_panel[int(i / 2)][i % 2].spines['top'].set_visible(False)

bottom_panel[0][0].set_xticks([])
bottom_panel[0][1].set_xticks([])
bottom_panel[1][0].set_xticks([])
bottom_panel[1][1].set_xticks([])
bottom_panel[2][0].set_xlabel('Time (ms)')
bottom_panel[2][1].set_xlabel('Time (ms)')

bottom_panel[0][0].set_ylabel('Current (A/F)')
bottom_panel[1][0].set_ylabel('Current (A/F)')
bottom_panel[2][0].set_ylabel('Current (A/F)')

fig.savefig(fig_dir + "AP-CS_currents.svg", format='svg')

# #########################
# #
# # Plot hERG current of CS model under Milnes protocol
# #
# #########################
# # Load steady state IKr data for drug free, addition of dofetilide
# # and verapamil conditions (Milnes' protocol)
# drug_free_log = myokit.DataLog.load_csv(
#     data_dir + 'drug_free_Milnes_current.csv')

# fig = modelling.figures.FigureStructure(figsize=(5, 4), gridspec=(3, 1))

# plot.add_single(fig.axs[0][0], drug_free_log, 'membrane.V')
# plot.state_occupancy_plot(fig.axs[2][0], drug_free_log, model,
#                           color_seq=color_seq)
# plot.add_single(fig.axs[1][0], drug_free_log, 'ikr.IKr')

# l_lim, u_lim = fig.axs[1][0].get_ylim()

# pulse_time = 25e3
# fig.sharex(['Time (s)'], [(0, pulse_time)])
# fig.sharey(['Voltage\n(mV)', 'Current (A/F)', 'State\noccupancy'])
# fig.adjust_ticks(fig.axs[2][0], pulse_time)
# fig.savefig(fig_dir + "CS_drugfree.svg", format='svg')

# # #########################
# # #
# # # Plot hERG current of SD model under Milnes protocol
# # #
# # #########################
# trapped_log = myokit.DataLog.load_csv(
#     data_dir + 'dofetilide_Milnes_current.csv')

# fig = modelling.figures.FigureStructure(figsize=(5, 4), gridspec=(3, 1))

# plot.add_single(fig.axs[0][0], trapped_log, 'membrane.V')
# plot.state_occupancy_plot(fig.axs[2][0], trapped_log, model,
#                           color_seq=color_seq, legend=False)
# plot.add_single(fig.axs[1][0], trapped_log, 'ikr.IKr',
#                 label='30 nM dofetilide')
# fig.axs[1][0].set_ylim(bottom=l_lim, top=u_lim)

# fig.sharex(['Time (s)'], [(0, pulse_time)])
# fig.sharey(['Voltage\n(mV)', 'Current (A/F)', 'State\noccupancy'])
# fig.adjust_ticks(fig.axs[2][0], pulse_time)
# fig.savefig(fig_dir + "SD_dofetilide.svg", format='svg')

# #########################
# #
# # Simulate and plot AP of AP-CS model
# #
# #########################
# # Plot AP of AP-CS model
# fig = modelling.figures.FigureStructure(figsize=(5, 5), gridspec=(4, 1))

# plot.add_single(fig.axs[0][0], APlog, 'stimulus.i_stim')
# plot.add_single(fig.axs[1][0], APlog, 'membrane.V')
# plot.state_occupancy_plot(fig.axs[3][0], APlog, APmodel, color_seq=color_seq)
# plot.add_single(fig.axs[2][0], APlog, 'ikr.IKr')

# l_lim, u_lim = fig.axs[2][0].get_ylim()

# fig.sharex(['Time (s)'], [(0, pulse_time)])
# fig.sharey(['Current\nstimulus', 'Voltage\n(mV)',
#             r"$I_\mathrm{Kr}$ (A/F)", 'State\noccupancy'])
# fig.savefig(fig_dir + "AP-CS_drugfree.svg", format='svg')

# # #########################
# # #
# # # Plot AP of AP-SD model
# # #
# # #########################
# trapped_APlog = myokit.DataLog.load_csv(
#     data_dir + 'APclamp.csv')

# # Plot AP of AP-CS model
# fig = modelling.figures.FigureStructure(figsize=(5, 5), gridspec=(4, 1))

# plot.add_single(fig.axs[0][0], trapped_APlog, 'stimulus.i_stim')
# plot.add_single(fig.axs[1][0], trapped_APlog, 'membrane.V')
# plot.state_occupancy_plot(fig.axs[3][0], trapped_APlog, APmodel,
#                           color_seq=color_seq, legend=False)
# plot.add_single(fig.axs[2][0], trapped_APlog, 'ikr.IKr',
#                 label='30 nM dofetilide')
# fig.axs[2][0].set_ylim(bottom=l_lim, top=u_lim)

# fig.sharex(['Time (s)'], [(0, pulse_time)])
# fig.sharey(['Current\nstimulus', 'Voltage\n(mV)',
#             r"$I_\mathrm{Kr}$ (A/F)", 'State\noccupancy'])
# fig.savefig(fig_dir + "AP-SD_dofetilide.svg", format='svg')

# #########################
# #
# # Plot IKr and AP of AP-CS model
# #
# #########################
# base_conductance = APmodel.get('ikr.gKr').value()
# APlog_inhibit = AP_model.conductance_simulation(base_conductance * 0.5,
#                                                 repeats)

# # Plot AP of AP-CS model
# fig = modelling.figures.FigureStructure(figsize=(5, 4), gridspec=(2, 1),
#                                         height_ratios=[1, 1])

# plot.add_multiple(fig.axs[0][0], [APlog, APlog_inhibit], 'membrane.V')
# plot.add_multiple(fig.axs[1][0], [APlog, APlog_inhibit], 'ikr.IKr')

# fig.sharex(['Time (s)'], [(0, pulse_time)])
# fig.sharey(['Voltage\n(mV)', r"$I_\mathrm{Kr}$ (A/F)"])
# fig.savefig(fig_dir + "AP_prolong.svg", format='svg')
