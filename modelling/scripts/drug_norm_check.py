# To make sure that the hERG current and action potential are not affected by
# the normalisation

import matplotlib.pyplot as plt
import myokit
import numpy as np
import os
import pandas as pd

import modelling

saved_data_dir = '../../simulation_data/sensitivity_analysis/EC50/'
saved_fig_dir = '../../figures/sensitivity_analysis/EC50/'
check_plot = True

# Load IKr model
model = '../../model/ohara-cipa-v1-2017-IKr.mmt'
model, _, x = myokit.load(model)

protocol_name = 'Milnes'
protocol_params = modelling.ProtocolParameters()
pulse_time = protocol_params.protocol_parameters[protocol_name]['pulse_time']
protocol = protocol_params.protocol_parameters[protocol_name]['function']

drug_model = modelling.BindingKinetics(model)
drug_model.protocol = protocol
hERG_base_conductance = drug_model.original_constants['gKr']

repeats = 1000

# Set AP model
APmodel = '../../model/ohara-cipa-v1-2017.mmt'

# Load model
APmodel, _, x = myokit.load(APmodel)

AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')
pulse_time_AP = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time_AP)
base_conductance = APmodel.get('ikr.gKr').value()

offset = 50
repeats_AP = 800
save_signal = 2
APD_points = 20
plotting_AP_pulse_time = pulse_time_AP * save_signal

Hill_model = modelling.HillsModel()

param_lib = modelling.BindingParameters()
drug_list = param_lib.drug_compounds

SA_model = modelling.ParameterCategory()
param_names = SA_model.param_names
parameter_interest = 'EC50'

filename = 'error.csv'
if os.path.exists(saved_data_dir + filename):
    previous_df = pd.read_csv(saved_data_dir + filename,
                              header=[0], index_col=[0],
                              skipinitialspace=True)
    ran_drugs = previous_df['drug'].values
else:
    ran_drugs = []

drug_list = [i for i in drug_list if i not in ran_drugs]
print(drug_list)
for drug in drug_list:
    # Set up variables
    Vhalf = param_lib.binding_parameters[drug]['Vhalf']
    Kmax = param_lib.binding_parameters[drug]['Kmax']
    Ku = param_lib.binding_parameters[drug]['Ku']
    Hill_n = param_lib.binding_parameters[drug]['N']
    half_effect_conc = param_lib.binding_parameters[drug]['EC50']

    all_params = [Vhalf, Kmax, Ku, Hill_n, half_effect_conc]

    orig_param_values = pd.DataFrame(all_params, index=param_names)
    orig_param_values = orig_param_values.T
    ComparisonController = modelling.ModelComparison(orig_param_values)

    ComparisonController.drug_param_values = orig_param_values

    # Simulate with actual drug concentration
    Hill_curve_coefs, drug_conc_Hill, peaks_norm = \
        ComparisonController.compute_Hill(drug_model)

    total_log = []
    for i in range(len(drug_conc_Hill)):
        log = drug_model.custom_simulation(
            orig_param_values, drug_conc_Hill[i], 1000,
            log_var=['engine.time', 'ikr.IKr'])
        total_log.append(log)

    log_conductance = []
    for i in range(len(drug_conc_Hill)):
        reduction_scale = Hill_model.simulate(
            Hill_curve_coefs, drug_conc_Hill[i])
        log = drug_model.conductance_simulation(
            hERG_base_conductance * reduction_scale, 1000,
            log_var=['engine.time', 'ikr.IKr'])
        log_conductance.append(log)

    # Normalise drug concentration
    drug_conc_Hill_norm = [i / (np.power(half_effect_conc, 1 / Hill_n))
                           for i in drug_conc_Hill]

    # Change EC50 value to 1
    param_values = orig_param_values
    param_values[parameter_interest][0] = 1
    ComparisonController.drug_param_values = param_values

    # Simulate with normalised drug concentration
    Hill_curve_coefs_norm, drug_conc_Hill_norm, peaks_norm2 = \
        ComparisonController.compute_Hill(drug_model,
                                          drug_conc=drug_conc_Hill_norm)

    log_norm = []
    for i in range(len(drug_conc_Hill_norm)):
        log = drug_model.custom_simulation(
            param_values, drug_conc_Hill_norm[i], 1000,
            log_var=['engine.time', 'ikr.IKr'])
        log_norm.append(log)

    log_conductance_norm = []
    for i in range(len(drug_conc_Hill_norm)):
        reduction_scale = Hill_model.simulate(
            Hill_curve_coefs_norm, drug_conc_Hill_norm[i])
        log = drug_model.conductance_simulation(
            hERG_base_conductance * reduction_scale,
            1000, log_var=['engine.time', 'ikr.IKr'])
        log_conductance_norm.append(log)

    # Rescale normalised drug concentration for plotting purposes
    drug_conc_Hill_rescale = [i * (np.power(half_effect_conc, 1 / Hill_n))
                              for i in drug_conc_Hill_norm]

    # Plot peak current simulated with actual drug concentration and
    # normalised drug concentration
    plt.figure()
    plt.plot(drug_conc_Hill, peaks_norm, 'b',
             label='Actual drug concentration')
    plt.plot(drug_conc_Hill_rescale, peaks_norm2, 'r',
             label='Normalised drug concentration (rescaled)')
    plt.xscale('log')
    plt.xlabel('Drug concentration')
    plt.ylabel('Normalised peak current')
    plt.legend()
    plt.savefig(saved_fig_dir + 'peak_compare_' + drug + '.pdf')
    plt.close()

    # Plot hERG current simulated with actual drug concentration and
    # normalised drug concentration (same color scale)
    fig = modelling.figures.FigureStructure(figsize=(10, 3),
                                            gridspec=(2, 2),
                                            hspace=0.3,
                                            height_ratios=[1, 1])
    plot = modelling.figures.FigurePlot()

    labels = [str(i) + ' nM' for i in drug_conc_Hill]
    labels_norm = [str(i) for i in drug_conc_Hill_norm]
    plot.add_multiple(fig.axs[0][0], total_log, 'ikr.IKr', labels=labels)
    plot.add_multiple(fig.axs[1][0], log_norm, 'ikr.IKr',
                      labels=labels_norm)
    plot.add_multiple(fig.axs[0][1], log_conductance, 'ikr.IKr', labels=labels)
    plot.add_multiple(fig.axs[1][1], log_conductance_norm, 'ikr.IKr',
                      labels=labels_norm)

    fig.axs[0][0].set_title('Drug concentration - SD model')
    fig.axs[1][0].set_title('Normalised drug concentration - SD model')
    fig.axs[0][1].set_title('Drug concentration - CS model')
    fig.axs[1][1].set_title('Normalised drug concentration - CS model')
    fig.sharex(['Time (s)'] * 2, [(0, pulse_time)] * 2)
    fig.sharey(['Current (A/F)'] * 2)
    fig.adjust_ticks(fig.axs[1][0], pulse_time)
    fig.savefig(saved_fig_dir + 'hERG_comparison_' + drug + '.pdf')

    # Compare action potentials
    # Simulate AP with actual drug concentration
    drug_conc_AP = 10**np.linspace(np.log10(drug_conc_Hill[1]),
                                   np.log10(max(drug_conc_Hill)),
                                   APD_points)

    AP_log = []
    AP_conductance = []
    APD_actual = []
    APD_conductance = []
    orig_param_values[parameter_interest][0] = half_effect_conc
    for i in range(len(drug_conc_AP)):
        log = AP_model.custom_simulation(orig_param_values, drug_conc_AP[i],
                                         1000, save_signal=save_signal,
                                         abs_tol=1e-7, rel_tol=1e-8,
                                         log_var=['engine.time', 'membrane.V'])
        AP_log.append(log)

        # Compute APD90
        APD_pulse = []
        for pulse in range(save_signal):
            apd90 = AP_model.APD90(log['membrane.V', pulse], offset, 0.1)
            APD_pulse.append(apd90)
        APD_actual.append(APD_pulse)

        # Run simulation for conductance model
        reduction_scale = Hill_model.simulate(
            Hill_curve_coefs, drug_conc_AP[i])
        d2 = AP_model.conductance_simulation(
            base_conductance * reduction_scale, 1000,
            timestep=0.1, save_signal=save_signal,
            abs_tol=1e-7, rel_tol=1e-8,
            log_var=['engine.time', 'membrane.V'])
        AP_conductance.append(d2)

        # Compute APD90
        APD_conductance_pulse = []
        for pulse in range(save_signal):
            apd90 = AP_model.APD90(d2['membrane.V', pulse], offset, 0.1)
            APD_conductance_pulse.append(apd90)
        APD_conductance.append(APD_conductance_pulse)

    APD_actual = [max(i) for i in APD_actual]
    APD_conductance = [max(i) for i in APD_conductance]

    # Normalise drug concentration
    drug_conc_AP_norm = [i / (np.power(half_effect_conc, 1 / Hill_n))
                         for i in drug_conc_AP]

    # Simulate AP with normalised drug concentration
    AP_log_norm = []
    AP_conductance_norm = []
    APD_normalised = []
    APD_conductance_norm = []
    param_values[parameter_interest][0] = 1
    for i in range(len(drug_conc_AP_norm)):
        log = AP_model.custom_simulation(param_values,
                                         drug_conc_AP_norm[i],
                                         1000, save_signal=save_signal,
                                         log_var=['engine.time', 'membrane.V'])
        AP_log_norm.append(log)

        # Compute APD90
        APD_pulse = []
        for pulse in range(save_signal):
            apd90 = AP_model.APD90(log['membrane.V', pulse], offset, 0.1)
            APD_pulse.append(apd90)
        APD_normalised.append(APD_pulse)

        # Run simulation for conductance model
        reduction_scale = Hill_model.simulate(
            Hill_curve_coefs_norm, drug_conc_AP_norm[i])
        d2 = AP_model.conductance_simulation(
            base_conductance * reduction_scale, 1000,
            timestep=0.1, save_signal=save_signal, abs_tol=1e-7,
            rel_tol=1e-8, log_var=['engine.time', 'membrane.V'])
        AP_conductance_norm.append(d2)

        # Compute APD90
        APD_conductance_pulse = []
        for pulse in range(save_signal):
            apd90 = AP_model.APD90(d2['membrane.V', pulse], offset, 0.1)
            APD_conductance_pulse.append(apd90)
        APD_conductance_norm.append(APD_conductance_pulse)

    APD_normalised = [max(i) for i in APD_normalised]
    APD_conductance_norm = [max(i) for i in APD_conductance_norm]

    RMSError_trap = ComparisonController.compute_RMSE(
        APD_actual, APD_normalised)
    MAError_trap = ComparisonController.compute_MAE(
        APD_actual, APD_normalised)
    RMSError_conductance = ComparisonController.compute_RMSE(
        APD_conductance, APD_conductance_norm)
    MAError_conductance = ComparisonController.compute_MAE(
        APD_conductance, APD_conductance_norm)

    fig = modelling.figures.FigureStructure(figsize=(10, 3),
                                            gridspec=(2, 2),
                                            hspace=0.3,
                                            height_ratios=[1, 1])
    plot = modelling.figures.FigurePlot()

    labels = [str(i) + ' nM' for i in drug_conc_AP]
    labels_norm = [str(i) for i in drug_conc_AP_norm]
    plot.add_multiple_continuous(fig.axs[0][0], AP_log,
                                 'membrane.V', labels=labels)
    plot.add_multiple_continuous(fig.axs[1][0], AP_log_norm,
                                 'membrane.V', labels=labels_norm)
    plot.add_multiple_continuous(fig.axs[0][1], AP_conductance,
                                 'membrane.V', labels=labels)
    plot.add_multiple_continuous(fig.axs[1][1], AP_conductance_norm,
                                 'membrane.V', labels=labels_norm)

    fig.axs[0][0].set_title('Drug concentration - SD model', fontsize=8)
    fig.axs[1][0].set_title('Normalised drug concentration - SD model',
                            fontsize=8)
    fig.axs[0][1].set_title('Drug concentration - CS model', fontsize=8)
    fig.axs[1][1].set_title('Normalised drug concentration - CS model',
                            fontsize=8)

    # unique = fig.legend_without_duplicate_labels(fig.axs[0][0])
    # fig.axs[0][0].legend(*zip(*unique), loc='right',
    #                      bbox_to_anchor=(1.45, 0.6))
    # unique = fig.legend_without_duplicate_labels(fig.axs[1][0])
    # fig.axs[1][0].legend(*zip(*unique), loc='right',
    #                      bbox_to_anchor=(1.45, 0.6))
    fig.sharex(['Time (ms)'] * 2, [(0, plotting_AP_pulse_time)] * 2)
    fig.sharey(['Voltage (mV)', 'Voltage (mV)'])
    fig.savefig(saved_fig_dir + 'AP_comparison_' + drug + '.pdf')

    # Rescale normalised drug concentration for plotting purposes
    drug_conc_AP_rescale = [i * (np.power(half_effect_conc, 1 / Hill_n))
                            for i in drug_conc_AP_norm]

    # Plot APD90 simulated with actual drug concentration and
    # normalised drug concentration
    fig = modelling.figures.FigureStructure(figsize=(10, 2),
                                            gridspec=(1, 2),)
    plot = modelling.figures.FigurePlot()

    fig.axs[0][0].plot(drug_conc_AP, APD_actual, 'o', color='b',
                       label='Actual drug conc')
    fig.axs[0][0].plot(drug_conc_AP_rescale, APD_normalised, '^', color='r',
                       label='Norm drug conc (rescaled)')
    fig.axs[0][1].plot(drug_conc_AP, APD_conductance, 'o', color='b',
                       label='Actual drug conc')
    fig.axs[0][1].plot(drug_conc_AP_rescale, APD_conductance_norm, '^',
                       color='r', label='Norm drug conc (rescaled)')
    for i in range(2):
        fig.axs[0][i].set_xscale('log')
        fig.axs[0][i].set_xlabel('Drug concentration')
        fig.axs[0][i].set_ylabel(r'APD$_{90}$ (ms)')
        fig.axs[0][i].legend()
    fig.savefig(saved_fig_dir + 'APD_compare_' + drug + '.pdf')
    plt.close()

    filename = 'error.csv'
    error_df = pd.DataFrame([drug, RMSError_trap, MAError_trap,
                             RMSError_conductance, MAError_conductance],
                            index=['drug', 'RMSError_trap', 'MAError_trap',
                                   'RMSError_conductance',
                                   'MAError_conductance'])
    error_df = error_df.T
    if os.path.exists(saved_data_dir + filename):
        previous_df = pd.read_csv(saved_data_dir + filename,
                                  header=[0], index_col=[0],
                                  skipinitialspace=True)
        comb_df = pd.concat([previous_df, error_df])
    else:
        comb_df = error_df
    comb_df.to_csv(saved_data_dir + filename)
