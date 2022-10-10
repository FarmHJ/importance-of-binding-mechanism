# Probably repeating sensitivity_analysis.py
# To plot the range of parameter values
# Drug binding-related parameters

import csv
import itertools
import matplotlib
import matplotlib.pyplot as plt
import myokit
import numpy as np
import pandas as pd

import modelling

saved_data_dir = '../../simulation_data/sensitivity_analysis/'

testing_fig_dir = '../../figures/testing/'
final_fig_dir = '../../figures/binding_kinetics_comparison/' + \
    'OHaraCiPA_model/sensitivity_analysis/'

check_plot = True
saved_fig_dir = final_fig_dir

param_lib = modelling.BindingParameters()
drug_list = param_lib.drug_compounds

Vhalf = []
Kmax = []
Ku = []
N = []
EC50 = []

for drug in drug_list:
    Vhalf.append(param_lib.binding_parameters[drug]['Vhalf'])
    Kmax.append(param_lib.binding_parameters[drug]['Kmax'])
    Ku.append(param_lib.binding_parameters[drug]['Ku'])
    N.append(param_lib.binding_parameters[drug]['N'])
    EC50.append(param_lib.binding_parameters[drug]['EC50'])

all_params = [Vhalf, Kmax, Ku, N, EC50]
labels = ['Vhalf', 'Kmax', 'Ku', 'N', 'EC50']
param_ranges = [[-150, -50], [30, 1e6], [5e-3, 7e-2], [0.8, 1.2], [1e3, 1e8]]

d = {}
for i, param_name in enumerate(labels):
    d[param_name] = param_ranges[i]

param_cat_df = pd.DataFrame(data=d, index=['low_border', 'high_border'])

# Plot of categorisation of parameters
# Set up colormap
cmap = matplotlib.cm.get_cmap('Set1')
# color_list = [cmap(i) for i in range(9)] + [cmap(i) for i in range(3)]
color_list = [cmap(i) for i in range(9)] + \
    ['darkturquoise', 'darkgreen', 'magenta']
repeating_cmap = matplotlib.colors.ListedColormap(color_list)

nrow = 2
# fig = modelling.figures.FigureStructure(
#     figsize=(10, 2.7 * nrow),
#     gridspec=(nrow, 3), hspace=0.7,
#     wspace=0.2,
#     height_ratios=[1] * nrow)

APD_diff = [0, 1, 2, 3, 1, 4, 2, 2, 2, 2, 2, 3]
markers = ['*', 'o', 'X', '2', '^']
APD_rough_diff = [0, 1, 1, 2, 1, 0, 1, 1, 1, 1, 1, 2]
discrete_colors = ['red', 'blue', 'black', 'purple', 'olive']
APD_rough_diff_label = ['similar', 'SD higher', 'CS higher']
# for i, label in enumerate(labels):
#     fig.axs[int(i / 3)][i % 3].axhline(
#         param_ranges[i][0], xmin=0, xmax=len(all_params[i]),
#         c='red')
#     fig.axs[int(i / 3)][i % 3].axhline(
#         param_ranges[i][1], xmin=0, xmax=len(all_params[i]),
#         c='red')
#     for j in range(len(drug_list)):
#         fig.axs[int(i / 3)][i % 3].scatter(
#             np.arange(len(all_params[i]))[j], all_params[i][j],
#             # c=np.arange(len(drug_list))[j], cmap=repeating_cmap[j],
#             marker=markers[APD_diff[j]])
#     fig.axs[int(i / 3)][i % 3].set_title(label)
#     fig.axs[int(i / 3)][i % 3].set_xticks(np.arange(12), labels=drug_list)
#     plt.setp(fig.axs[int(i / 3)][i % 3].get_xticklabels(), rotation=45,
#              ha='right', rotation_mode='anchor')

# fig.axs[0][2].axvline(2.5, ymin=1.5e-5, ymax=np.exp(0.3), ls='--')
# fig.axs[0][2].axvline(5.5, ymin=1.5e-5, ymax=np.exp(0.3), ls='--')
# fig.axs[0][1].set_yscale("log", nonpositive='clip')
# fig.axs[0][2].set_yscale("log", nonpositive='clip')
# fig.axs[1][1].set_yscale("log", nonpositive='clip')

# fig.savefig(saved_fig_dir + 'parameter_categorisation.svg', format='svg')

# Plot the color-coded APD90 changes in a 3-dimensional plot of the 3
# interested parameters, i.e. Vhalf, Kmax and Ku
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
for j in range(len(drug_list)):
    ax.scatter(Vhalf[j], np.log(Kmax[j]), np.log(Ku[j]),
            #    c=color_list[j], label=str(drug_list[j]),
               c=discrete_colors[APD_rough_diff[j]], s=100,
               label=APD_rough_diff_label[APD_rough_diff[j]])
            #    marker=markers[APD_diff[j]], s=100)

handles, labels = ax.get_legend_handles_labels()
unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if
            l not in labels[:i]]
ax.legend(*zip(*unique), loc='upper left', bbox_to_anchor=(1.0, 1.0))
# ax.set_facecolor('silver')
ax.set_xlabel('Vhalf')
ax.set_ylabel('Kmax')
ax.set_zlabel('Ku')
plt.show()

# # Plot transition rates against voltage and drug concentration
# # During Milnes protocol and action potential


# def Kb_fn(drug_conc_norm, Kmax, n, Ku):
#     Kb = np.divide(Ku * Kmax * drug_conc_norm**n, (1 + drug_conc_norm**n))

#     return Kb


# def Kn_fn(V, Vhalf):
#     Kt = 3.5e-5
#     Kn = Kt / (1 + np.exp((Vhalf - V) / 6.789))

#     return Kn


# drug_conc_norm = 10**np.linspace(-11, 2, 28)
# voltage_arr = np.arange(-150, 50, 10)
# Kb_drugs = []
# Kn_drugs = []
# for i, drug in enumerate(drug_list):
#     Kb = Kb_fn(drug_conc_norm, Kmax[i], N[i], Ku[i])
#     Kb_drugs.append(Kb)

#     Kn = Kn_fn(voltage_arr, Vhalf[i])
#     Kn_drugs.append(Kn)

# fig = modelling.figures.FigureStructure(
#     figsize=(10, 3),
#     gridspec=(1, 3), wspace=0.2)

# for i in range(len(drug_list)):
#     fig.axs[0][0].plot(drug_conc_norm, Kb_drugs[i], c=repeating_cmap(i))
#     fig.axs[0][1].hlines(Ku[i], 0, 100, colors=repeating_cmap(i))
#     fig.axs[0][2].plot(voltage_arr, Kn_drugs[i], c=repeating_cmap(i),
#                        label=drug_list[i])

# fig.axs[0][0].set_title('Kb')
# fig.axs[0][0].set_yscale('log')
# fig.axs[0][0].set_xscale('log')
# fig.axs[0][0].set_xlabel('Normalised drug concentration')

# fig.axs[0][1].set_title('Ku')
# fig.axs[0][1].set_xlabel('Normalised drug concentration')

# fig.axs[0][2].set_title('Kn')
# fig.axs[0][2].set_xlabel('Voltage')
# fig.axs[0][2].legend(loc='lower left', bbox_to_anchor=(1.0, 0))
# fig.savefig(saved_fig_dir + 'transition_rates.pdf')

# # Simple SA - taking mean value of each category
# # Split parameter value ranges to three categories: low, medium and high
# category_means = []
# for i, param in enumerate(all_params):
#     low_cat = [k for k in param if k < param_ranges[i][0]]
#     med_cat = [k for k in param if k >= param_ranges[i][0]
#                and k <= param_ranges[i][1]]
#     high_cat = [k for k in param if k > param_ranges[i][1]]
#     category_means.append([np.mean(low_cat), np.mean(med_cat),
#                            np.mean(high_cat)])

# d = {}
# for i, param_name in enumerate(labels):
#     d[param_name] = category_means[i]
# param_cat_mean_df = pd.DataFrame(
#     data=d, index=['low_gp_mean', 'mid_gp_mean', 'high_gp_mean'])

# param_cat_df = pd.concat([param_cat_df, param_cat_mean_df])

# param_cat_df.to_csv(saved_data_dir + 'parameter_categorisation.csv')

# # Load IKr model
# model = '../../model/ohara-cipa-v1-2017-IKr.mmt'
# model, _, x = myokit.load(model)

# protocol_params = modelling.ProtocolParameters()
# pulse_time = protocol_params.protocol_parameters['Milnes']['pulse_time']
# protocol = protocol_params.protocol_parameters['Milnes']['function']

# drug_model = modelling.BindingKinetics(model)
# drug_model.protocol = protocol

# drug_conc_Hill = list(np.append(0, 10**np.linspace(-1, 10, 12)))

# # Set AP model
# APmodel = '../../model/ohara-cipa-v1-2017.mmt'

# # Load model
# APmodel, _, x = myokit.load(APmodel)

# AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')
# pulse_time = 1000
# AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)
# base_conductance = APmodel.get('ikr.gKr').value()

# offset = 50
# repeats_AP = 800
# save_signal = 2

# drug_conc_AP = 10**np.linspace(-1, 1, 3)
# # drug_conc = 10**np.linspace(-1, 7, 20)

# # Run simulation
# repeats = 1000
# APD_dict = {}

# for num, param_comb in enumerate(itertools.product(*category_means)):
#     param_values = pd.DataFrame(param_comb, index=labels)
#     param_values = param_values.T

#     peaks = []
#     for i in range(len(drug_conc_Hill)):
#         log = drug_model.custom_simulation(
#             param_values, drug_conc_Hill[i], repeats,
#             log_var=['engine.time', 'ikr.IKr'])
#         peak, _ = drug_model.extract_peak(log, 'ikr.IKr')
#         peaks.append(peak[-1])

#     peaks = (peaks - min(peaks)) / (max(peaks) - min(peaks))

#     # Fit Hill curve
#     Hill_model = modelling.HillsModel()
#     optimiser = modelling.HillsModelOpt(Hill_model)
#     Hill_curve, _ = optimiser.optimise(drug_conc_Hill, peaks)

#     if check_plot:
#         max_grid = np.ceil(np.log(drug_conc_Hill[-1]))
#         conc_grid = np.arange(-2, max_grid, 1)

#         plt.figure(figsize=(4, 3))
#         plt.plot(np.log(drug_conc_Hill[1:]), peaks[1:], 'o',
#                  label='peak current')
#         plt.plot(conc_grid,
#                  Hill_model.simulate(Hill_curve[:2], np.exp(conc_grid)),
#                  'k', label='fitted Hill eq')
#         plt.xlabel('Drug concentration (log)')
#         plt.ylabel('Normalised peak current')
#         plt.tight_layout()
#         plt.savefig(saved_fig_dir + "Hill_curve_" + str(num) + ".pdf")

#     data_dict = {'drug_conc_Hill': drug_conc_Hill, 'peak_current': peaks,
#                  'Hill_curve': Hill_curve[:2]}

#     APD_trapping = []
#     APD_conductance = []
#     for i in range(len(drug_conc_AP)):
#         # Run simulation for trapping model
#         log = AP_model.custom_simulation(
#             param_values, drug_conc_AP[i], repeats_AP, timestep=0.1,
#             save_signal=save_signal,
#             log_var=['engine.time', 'membrane.V'])

#         # Compute APD90
#         APD_trapping_pulse = []
#         for pulse in range(save_signal):
#             apd90 = AP_model.APD90(log['membrane.V', pulse], offset, 0.1)
#             APD_trapping_pulse.append(apd90)
#         APD_trapping.append(APD_trapping_pulse)

#         # Run simulation for conductance model
#         reduction_scale = Hill_model.simulate(Hill_curve[:2], drug_conc_AP[i])
#         d2 = AP_model.conductance_simulation(
#             base_conductance * reduction_scale, repeats_AP, timestep=0.1,
#             save_signal=save_signal, abs_tol=1e-6, rel_tol=1e-5,
#             log_var=['engine.time', 'membrane.V'])

#         # Compute APD90
#         APD_conductance_pulse = []
#         for pulse in range(save_signal):
#             apd90 = AP_model.APD90(d2['membrane.V', pulse], offset, 0.1)
#             APD_conductance_pulse.append(apd90)
#         APD_conductance.append(APD_conductance_pulse)

#     APD_trapping = [max(i) for i in APD_trapping]
#     APD_conductance = [max(i) for i in APD_conductance]

#     MSError = np.sum((np.array(APD_trapping) - np.array(APD_conductance))**2) \
#         / len(APD_trapping)

#     APD_result = {**data_dict, 'param_seq': labels,
#                   'param_values': param_values.values[0],
#                   'drug_conc_AP': drug_conc_AP, 'APD_trapping': APD_trapping,
#                   'APD_conductance': APD_conductance, 'MSE': MSError}

#     APD_dict[num] = APD_result

#     print(APD_dict)

# w = csv.writer(open(saved_data_dir + 'SA_param_categories.csv', 'w'))
# for key, val in APD_dict.items():
#     w.writerow([key, val])
