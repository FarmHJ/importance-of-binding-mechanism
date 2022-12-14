# Calibrate the ionic conductance of the CS model from the SD model and compare
# the APD90 at steady state and transient phase.
# Output:
# 1. hERG current at various drug concentration simulated from CiPA hERG model.
# 2. Normalised peak hERG current vs drug concentration in log scale.
# 3. Fitted Hill curve over peak hERG current from CiPA hERG model.
# 4. hERG current at various drug concentration simulated from conductance
#    calibrated hERG model.
# 5. Comparison of peak hERG current between both models, i.e. CiPA hERG model
#    and conductance hERG model.
# 6. Comparison of hERG current between both models.
# 7. 2 pulses of action potential simulated from both models (steady state).
# 8. APD90 of two pulses for both models (steady state).
# 9. 10 pulses of action potentials simulated from both models
#    (transient phase).
# 10. APD90 of 10 pulses for both models (transient phase).

import matplotlib
import matplotlib.pyplot as plt
import myokit
import numpy as np
import os
import pandas as pd
import pints
import sys

import modelling

# Define drug and protocol
drug = sys.argv[1]
protocol_name = sys.argv[2]
protocol_params = modelling.ProtocolParameters()
pulse_time = protocol_params.protocol_parameters[protocol_name]['pulse_time']
protocol = protocol_params.protocol_parameters[protocol_name]['function']

# Define drug concentration range
drug_conc_lib = modelling.DrugConcentrations()
drug_conc = drug_conc_lib.drug_concentrations[drug]['coarse']
repeats = 1000
drug_labels = [str(i) + ' nM' for i in drug_conc]

data_dir = '../../testing_data/model_comparison/' + drug + '/' + \
    protocol_name + '/'
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)
result_filename = 'Hill_curve.txt'

fig_dir = '../../testing_figures/model_comparison/' + drug + '/' + \
    protocol_name + '/'
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

# Load current model
model = '../../math_model/ohara-cipa-v1-2017-IKr.mmt'
model, _, x = myokit.load(model)
current_model = modelling.BindingKinetics(model)
current_model.protocol = protocol

abs_tol = 1e-7
rel_tol = 1e-10

# Simulate IKr of the SD model for a range of drug concentrations
# Extract the peak of IKr
total_log = []
peaks = []
for i in range(len(drug_conc)):
    log = current_model.drug_simulation(
        drug, drug_conc[i], repeats,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'],
        abs_tol=abs_tol, rel_tol=rel_tol)
    peak, _ = current_model.extract_peak(log, 'ikr.IKr')
    peaks.append(peak[-1])
    total_log.append(log)

    log.save_csv(data_dir + 'SD_hERG_' + str(drug_conc[i]) + '.csv')

# Plot IKr for different drug concentrations
fig_plot = modelling.figures.ReferenceStructure()
fig = fig_plot.current_concs(total_log, pulse_time, drug_conc)
fig.savefig(fig_dir + "hERG_trapping_" + drug + "_concs.pdf")

# Normalise drug response (peak current)
peaks = (peaks - min(peaks)) / (max(peaks) - min(peaks))

# Plot drug response (peak current) against drug concentration curve
plt.rcParams.update({'font.size': 9})

plt.figure(figsize=(4, 3))
plt.plot(np.log(drug_conc[1:]), peaks[1:], 'o')
plt.xlabel('Drug concentration (log)')
plt.ylabel('Normalised peak current')
plt.tight_layout()
plt.savefig(fig_dir + "peak_hERG_trapping_" + drug + "_concs.pdf")
plt.close()

# Fit drug response to Hill curve
Hill_model = modelling.HillsModel()
optimiser = modelling.HillsModelOpt(Hill_model)
if not os.path.isfile(data_dir + result_filename):
    estimates, _ = optimiser.optimise(drug_conc, peaks)
    with open(data_dir + result_filename, 'w') as f:
        for x in estimates:
            f.write(pints.strfloat(x) + '\n')
else:
    estimates = np.loadtxt(data_dir + result_filename, unpack=True)
    estimates = np.array(estimates)

# Plot fitting result of Hill curve
max_grid = np.ceil(np.log(drug_conc[-1]))
conc_grid = np.arange(-3, max_grid, 1)

plt.figure(figsize=(4, 3))
plt.plot(np.log(drug_conc[1:]), peaks[1:], 'o', label='peak current')
plt.plot(conc_grid, Hill_model.simulate(estimates[:2], np.exp(conc_grid)),
         'k', label='fitted Hill eq')
plt.xlabel('Drug concentration (log)')
plt.ylabel('Normalised peak current')
plt.tight_layout()
plt.legend()
plt.savefig(fig_dir + "fitted_peak_hERG_trapping_" + drug + "_concs.pdf")
plt.close()

# Compare peak current with the CS model which has its ionic conductance
# calibrated
conductance = model.get('ikr.gKr').value()

peaks_conductance = []
log_conductance = []

current_model.current_head = current_model.model.get('ikr')
for i in range(len(drug_conc)):

    reduction_scale = Hill_model.simulate(estimates[:2], drug_conc[i])
    d2 = current_model.conductance_simulation(
        conductance * reduction_scale, repeats,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'],
        abs_tol=abs_tol, rel_tol=rel_tol)
    log_conductance.append(d2)

    d2.save_csv(data_dir + 'CS_hERG_' + str(drug_conc[i]) + '.csv')

    peak, _ = current_model.extract_peak(d2, 'ikr.IKr')
    peaks_conductance.append(peak[-1])

# Plot IKr of both the SD model and the CS model
fig = modelling.figures.FigureStructure(figsize=(5, 4),
                                        gridspec=(3, 1))
plot = modelling.figures.FigurePlot()
cmap = matplotlib.cm.get_cmap('viridis')

labels = [str(i) + ' nM' for i in drug_conc]
plot.add_single(fig.axs[0][0], total_log[0], 'membrane.V', color='k')
plot.add_multiple(fig.axs[1][0], total_log, 'ikr.IKr', labels=labels,
                  color=cmap)
plot.add_multiple(fig.axs[2][0], log_conductance, 'ikr.IKr', labels=labels,
                  color=cmap)

fig.axs[2][0].legend()
fig.sharex(['Time (s)'], [(0, pulse_time)])
fig.sharey(['Voltage (mV)', 'Current (A/F)\ntrapped model',
            'Current (A/F)\nnontrapped model'])
fig.adjust_ticks(fig.axs[2][0], pulse_time)
fig.savefig(fig_dir + 'combine_comparison_viridis.pdf')

# Check if the peak IKr of both models match
row, col = 3, 3
fig = modelling.figures.ReferenceStructure()
fig.hERG_compare(total_log, log_conductance, drug_conc, pulse_time,
                 grid=(row, col))
plt.savefig(fig_dir + "hERG_compare_" + drug + "_concs.pdf")

# Check Hill curve generated by the CS model
# Make sure they are the same as the Hill curve generated by the SD model
peaks_conductance_norm = (peaks_conductance - min(peaks_conductance)) / (
    max(peaks_conductance) - min(peaks_conductance))

plt.figure(figsize=(4, 3))
plt.plot(np.log(drug_conc[1:]), peaks[1:],
         'o', label='SD model')
plt.plot(np.log(drug_conc[1:]), peaks_conductance_norm[1:],
         'o', label='CS model')
plt.xlabel('Drug concentration (log)')
plt.ylabel('Normalised peak current')
plt.legend()
plt.tight_layout()
plt.savefig(fig_dir + "norm_peak_compare_hERG_" + drug + "_concs.pdf")

#
# Propagate the effect to action potential level
#
# Load AP model and the current impulse
APmodel = '../../math_model/ohara-cipa-v1-2017.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')
pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)
conductance = APmodel.get('ikr.gKr').value()

# Define the number of repeating pulses
offset = 50
repeats = 1000
save_signal = 2

# Simulate AP of the ORd-SD model and the ORd-CS model till steady state
# Then compute their APD90s and the peak current
AP_conductance = []
APD_conductance = []
# hERG_peak_conductance = []

AP_trapping = []
APD_trapping = []
# hERG_peak_trapping = []

for i in range(len(drug_conc)):
    print('simulating concentration: ' + str(drug_conc[i]))
    log = AP_model.drug_simulation(
        drug, drug_conc[i], repeats, save_signal=save_signal,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'],
        abs_tol=abs_tol, rel_tol=rel_tol)
    log.save_csv(data_dir + 'SD_AP_conc' + str(drug_conc[i]) + '_steady.csv')
    APD_trapping_pulse = []
    hERG_trapping_pulse = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(log['membrane.V', pulse], offset, 0.1)
        APD_trapping_pulse.append(apd90)

        # hERG_peak = np.max(log['ikr.IKr', pulse])
        # hERG_trapping_pulse.append(hERG_peak)

    AP_trapping.append(log)
    APD_trapping.append(APD_trapping_pulse)
    # hERG_peak_trapping.append(hERG_trapping_pulse)

    reduction_scale = Hill_model.simulate(estimates[:2], drug_conc[i])
    d2 = AP_model.conductance_simulation(
        conductance * reduction_scale, repeats, save_signal=save_signal,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'],
        abs_tol=abs_tol, rel_tol=rel_tol)
    d2.save_csv(data_dir + 'CS_AP_conc' + str(drug_conc[i]) + '_steady.csv')
    APD_conductance_pulse = []
    hERG_conductance_pulse = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(d2['membrane.V', pulse], offset, 0.1)
        APD_conductance_pulse.append(apd90)

        # hERG_peak = np.max(d2['ikr.IKr', pulse])
        # hERG_conductance_pulse.append(hERG_peak)

    AP_conductance.append(d2)
    APD_conductance.append(APD_conductance_pulse)
    # hERG_peak_conductance.append(hERG_conductance_pulse)

    print('done concentration: ' + str(drug_conc[i]))

# Save APD90 and peak current of both the ORd-SD model and the ORd-CS model
column_name = ['pulse ' + str(i) for i in range(save_signal)]
APD_trapping_df = pd.DataFrame(APD_trapping, columns=column_name)
APD_trapping_df['drug concentration'] = drug_conc
APD_trapping_df.to_csv(data_dir + 'SD_APD_pulses' +
                       str(int(save_signal)) + '.csv')
# hERG_peak_trapping_df = pd.DataFrame(hERG_peak_trapping, columns=column_name)
# hERG_peak_trapping_df['drug concentration'] = drug_conc
# hERG_peak_trapping_df.to_csv(data_dir + 'SD_hERGpeak_pulses' +
#                              str(int(save_signal)) + '.csv')
APD_conductance_df = pd.DataFrame(APD_conductance, columns=column_name)
APD_conductance_df['drug concentration'] = drug_conc
APD_conductance_df.to_csv(data_dir + 'CS_APD_pulses' +
                          str(int(save_signal)) + '.csv')
# hERG_peak_conductance_df = pd.DataFrame(hERG_peak_conductance,
#                                         columns=column_name)
# hERG_peak_conductance_df['drug concentration'] = drug_conc
# hERG_peak_conductance_df.to_csv(
#     data_dir + 'CS_hERGpeak_pulses' + str(int(save_signal)) + '.csv')

# Take the max APD90
APD_trapping = [max(i) for i in APD_trapping]
APD_conductance = [max(i) for i in APD_conductance]

# Remove repeated signals at high concentrations
second_EAD_trap = [i for i, e in enumerate(APD_trapping)
                   if e == 1000][1:]
second_EAD_conduct = [i for i, e in enumerate(APD_conductance)
                      if e == 1000][1:]
AP_trapping_plot = [e for i, e in enumerate(AP_trapping)
                    if i not in second_EAD_trap]
AP_conductance_plot = [e for i, e in enumerate(AP_conductance)
                       if i not in second_EAD_conduct]

# Plot AP and IKr at various drug concentrations - steady state
plotting_pulse_time = pulse_time * save_signal

fig = modelling.figures.FigureStructure(figsize=(8, 3),
                                        gridspec=(2, 2),
                                        height_ratios=[1, 1],
                                        wspace=0.08)
plot = modelling.figures.FigurePlot()
cmap = matplotlib.cm.get_cmap('viridis')

plot.add_multiple_continuous(fig.axs[0][0], AP_trapping_plot,
                             'membrane.V', cmap=cmap,
                             labels=drug_labels)
plot.add_multiple_continuous(fig.axs[1][0], AP_trapping_plot,
                             'ikr.IKr', cmap=cmap, labels=drug_labels)
plot.add_multiple_continuous(fig.axs[0][1], AP_conductance_plot,
                             'membrane.V', cmap=cmap,
                             labels=drug_labels)
plot.add_multiple_continuous(fig.axs[1][1], AP_conductance_plot,
                             'ikr.IKr', cmap=cmap, labels=drug_labels)
fig.axs[0][0].set_title('ORd-state dependent model', fontsize=8)
fig.axs[0][1].set_title('ORd-conductance scaling model', fontsize=8)

unique = fig.legend_without_duplicate_labels(fig.axs[1][1])
fig.axs[1][1].legend(*zip(*unique), loc='right',
                     bbox_to_anchor=(1.45, 0.6))
fig.sharex(['Time (ms)'] * 2, [(0, plotting_pulse_time)] * 2)
fig.sharey(['Voltage (mV)', 'Current (A/F)'])
fig.savefig(fig_dir + "2AP_steady_trapping_nontrapping.pdf")

# Identify EAD-like behaviour
EAD_marker = [1050 if (i >= 1000 or j >= 1000) else None for (i, j)
              in zip(APD_trapping[1:], APD_conductance[1:])]
filename = "APD90_compare_2pulses_hERG_" + drug + "_concs.pdf"

# Plot APD90 of both models
plt.figure(figsize=(4, 3))
plt.plot(np.log(drug_conc[1:]), APD_trapping[1:],
         'o', color='orange', label='ORd-SD model')
plt.plot(np.log(drug_conc[1:]), APD_conductance[1:],
         '^', color='blue', label='ORd-CS model', alpha=0.8)
plt.plot(np.log(drug_conc[1:]), EAD_marker, 'o',
         color='k', marker=(5, 2))
plt.xlabel('Drug concentration (log)')
plt.ylabel(r'APD$_{90}$')
plt.legend()
plt.tight_layout()
plt.savefig(fig_dir + filename)
plt.close()

# Simulate AP of the ORd-SD model and the ORd-CS model for 300 pulses
# Then compute their APD90s and the peak current
repeats = 300
save_signal = repeats

AP_conductance = []
APD_conductance = []
hERG_peak_conductance = []

AP_trapping = []
APD_trapping = []
hERG_peak_trapping = []

for i in range(len(drug_conc)):
    print('simulating concentration: ' + str(drug_conc[i]))
    log = AP_model.drug_simulation(
        drug, drug_conc[i], repeats, save_signal=save_signal,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'],
        abs_tol=abs_tol, rel_tol=rel_tol)
    log.save_csv(data_dir + 'SD_AP_conc' + str(drug_conc[i]) +
                 '_transient.csv')
    APD_trapping_pulse = []
    hERG_trapping_pulse = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(log['membrane.V', pulse], offset, 0.1)
        APD_trapping_pulse.append(apd90)

        hERG_peak = np.max(log['ikr.IKr', pulse])
        hERG_trapping_pulse.append(hERG_peak)

    AP_trapping.append(log)
    APD_trapping.append(APD_trapping_pulse)
    hERG_peak_trapping.append(hERG_trapping_pulse)

    reduction_scale = Hill_model.simulate(estimates[:2], drug_conc[i])
    d2 = AP_model.conductance_simulation(
        conductance * reduction_scale, repeats, save_signal=save_signal,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'],
        abs_tol=abs_tol, rel_tol=rel_tol)
    d2.save_csv(data_dir + 'CS_AP_' + str(drug_conc[i]) + '.csv')
    APD_conductance_pulse = []
    hERG_conductance_pulse = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(d2['membrane.V', pulse], offset, 0.1)
        APD_conductance_pulse.append(apd90)

        hERG_peak = np.max(d2['ikr.IKr', pulse])
        hERG_conductance_pulse.append(hERG_peak)

    AP_conductance.append(d2)
    APD_conductance.append(APD_conductance_pulse)
    hERG_peak_conductance.append(hERG_conductance_pulse)

    print('done concentration: ' + str(drug_conc[i]))

# Save APD90 and peak current of both the ORd-SD model and the ORd-CS model
column_name = ['pulse ' + str(i) for i in range(save_signal)]
APD_trapping_df = pd.DataFrame(APD_trapping, columns=column_name)
APD_trapping_df['drug concentration'] = drug_conc
APD_trapping_df.to_csv(data_dir + 'SD_APD_pulses' +
                       str(int(save_signal)) + '.csv')
hERG_peak_trapping_df = pd.DataFrame(hERG_peak_trapping, columns=column_name)
hERG_peak_trapping_df['drug concentration'] = drug_conc
hERG_peak_trapping_df.to_csv(data_dir + 'SD_hERGpeak_pulses' +
                             str(int(save_signal)) + '.csv')
APD_conductance_df = pd.DataFrame(APD_conductance, columns=column_name)
APD_conductance_df['drug concentration'] = drug_conc
APD_conductance_df.to_csv(data_dir + 'CS_APD_pulses' +
                          str(int(save_signal)) + '.csv')
hERG_peak_conductance_df = pd.DataFrame(hERG_peak_conductance,
                                        columns=column_name)
hERG_peak_conductance_df['drug concentration'] = drug_conc
hERG_peak_conductance_df.to_csv(
    data_dir + 'CS_hERGpeak_pulses' + str(int(save_signal)) + '.csv')

# Plot peak of IKr - transient phase
cmap = matplotlib.cm.get_cmap('tab10')
norm = matplotlib.colors.Normalize(0, len(drug_conc))

hERG_reduc_trap = []
hERG_reduc_conduct = []

fig, ax = plt.subplots(1, 1, figsize=(4, 3))
for i in range(len(drug_conc)):
    ax.plot(np.arange(save_signal), np.array(hERG_peak_trapping[i]),
            'o-', label=str(drug_conc[i]) + 'nM: ORd-SD model',
            color=cmap(norm(i)))
    ax.plot(np.arange(save_signal), np.array(hERG_peak_conductance[i]),
            '^--', label=str(drug_conc[i]) + 'nM: ORd-CS model',
            color=cmap(norm(i)))

    # # peak reduction
    # hERG_peak_conc_trap = hERG_peak_trapping[i]
    # hERG_reduc_trap.append((hERG_peak_conc_trap[0]
    #                         - hERG_peak_conc_trap[-1])
    #                        / hERG_peak_conc_trap[0])

    # hERG_peak_conc_conduct = hERG_peak_conductance[i]
    # hERG_reduc_conduct.append((hERG_peak_conc_conduct[0]
    #                            - hERG_peak_conc_conduct[-1])
    #                           / hERG_peak_conc_conduct[0])

ax.legend(loc='upper right', bbox_to_anchor=(1.75, 1.5))
plt.xlabel('Sweeps')
plt.ylabel(r'APD$_{90}$')
plt.savefig(fig_dir + "peak_compare_hERG_" + drug + "_concs.pdf",
            bbox_inches='tight')
plt.close()

# plt.figure()
# plt.plot(hERG_reduc_trap)
# plt.plot(hERG_reduc_conduct)
# plt.savefig(fig_dir + 'peak_reduction_transient.pdf')
# plt.close()

# Plot APs of both models for first 10 pulses
labels = [str(i) + ' nM' for i in drug_conc]
cmap = matplotlib.cm.get_cmap('viridis')
plotting_pulse_time = pulse_time * save_signal

fig = modelling.figures.FigureStructure(figsize=(10, 4),
                                        gridspec=(2, 1),
                                        height_ratios=[1, 1])
plot = modelling.figures.FigurePlot()

plot.add_multiple_continuous(fig.axs[0][0], AP_trapping[:10], 'membrane.V',
                             cmap=cmap)
plot.add_multiple_continuous(fig.axs[1][0], AP_trapping[:10], 'ikr.IKr',
                             labels=labels, cmap=cmap)

fig.legend_without_duplicate_labels(fig.axs[1][0])
fig.sharex(['Time (s)'], [(0, plotting_pulse_time)])
fig.sharey(['Voltage (mV)', 'hERG current'])
fig.adjust_ticks(fig.axs[1][0], plotting_pulse_time)
fig.axs[0][0].set_title(drug + ", ORd-SD model")

fig.savefig(fig_dir + "10AP_transient_hERG_SD_" + drug + "_concs.pdf")
plt.close()

fig = modelling.figures.FigureStructure(figsize=(10, 4),
                                        gridspec=(2, 1),
                                        height_ratios=[1, 1])
plot = modelling.figures.FigurePlot()

plot.add_multiple_continuous(fig.axs[0][0], AP_conductance[:10],
                             'membrane.V', cmap=cmap)
plot.add_multiple_continuous(fig.axs[1][0], AP_conductance[:10], 'ikr.IKr',
                             labels=labels, cmap=cmap)

fig.legend_without_duplicate_labels(fig.axs[1][0])
fig.sharex(['Time (s)'], [(0, plotting_pulse_time)])
fig.sharey(['Voltage (mV)', 'hERG current'])
fig.adjust_ticks(fig.axs[1][0], plotting_pulse_time)
fig.axs[0][0].set_title(drug + ", ORD-CS model")

fig.savefig(fig_dir + "10AP_transient_hERG_CS_" + drug + "_concs.pdf")
plt.close('all')

# Plot APD90 of 10 pulses
cmap_trap = matplotlib.cm.get_cmap('Blues')
cmap_conduct = matplotlib.cm.get_cmap('YlOrBr')
cmap = matplotlib.cm.get_cmap('tab10')
norm = matplotlib.colors.Normalize(0, len(drug_conc))

fig, ax = plt.subplots(1, 1, figsize=(4, 3))
for i in range(len(drug_conc)):
    ax.plot(np.arange(10), np.array(APD_trapping[i][:10]), 'o-',
            label=str(drug_conc[i]) + 'nM: ORd-SD model',
            color=cmap(norm(i)))
    ax.plot(np.arange(10), np.array(APD_conductance[i][:10]),
            '^--', label=str(drug_conc[i]) + 'nM: ORd-CS model',
            color=cmap(norm(i)))

ax.legend(loc='upper right', bbox_to_anchor=(1.75, 1.5))
plt.xlabel('Sweeps')
plt.ylabel(r'APD$_{90}$')
plt.savefig(fig_dir + "10APD90_compare_hERG_" + drug + "_concs.pdf",
            bbox_inches='tight')
plt.close()

fig = modelling.figures.FigureStructure(figsize=(5, 2),
                                        gridspec=(1, 2),
                                        wspace=0.08)
for i in range(len(drug_conc)):
    for j in range(len(drug_conc)):
        if i == j:
            APD_plot = [APD_trapping[j][ind] for ind in
                        range(len(APD_trapping[j])) if ind % 2 == 0]
            fig.axs[int(i / 4)][i % 4].plot(
                np.arange(save_signal / 2) * 2, np.array(APD_plot),
                'o-', label='ORd-SD model',
                color='orange', zorder=5)
            APD_plot = [APD_conductance[j][ind] for ind in
                        range(len(APD_conductance[j])) if ind % 2 == 0]
            fig.axs[int(i / 4)][i % 4].plot(
                np.arange(save_signal / 2) * 2,
                np.array(APD_plot), '^--',
                label='ORd-CS model', color='blue', zorder=5)
    fig.axs[int(i / 4)][i % 4].set_title(str(drug_conc[i]) + 'nM')

fig.axs[0][1].plot(150, 1050, 'o', color='k', marker=(5, 2))
handles, labels = fig.axs[0][0].get_legend_handles_labels()
lgds = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if
        i <= 1]
fig.axs[0][0].legend(*zip(*lgds), loc='upper right')
fig.sharex(['Sweeps'] * 2)
fig.sharey([r"APD$_{90}$"])

fig.savefig(fig_dir + "APD90_compare_transient_multiples.pdf")
