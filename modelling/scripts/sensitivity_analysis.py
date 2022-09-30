# To plot the range of parameter values
# Drug binding-related parameters

import matplotlib.pyplot as plt
import myokit
import numpy as np
import os
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

drug = 'cisapride'
Vhalf = param_lib.binding_parameters[drug]['Vhalf']
Kmax = param_lib.binding_parameters[drug]['Kmax']
Ku = param_lib.binding_parameters[drug]['Ku']
Hill_n = param_lib.binding_parameters[drug]['N']
half_effect_conc = param_lib.binding_parameters[drug]['EC50']

all_params = [Vhalf, Kmax, Ku, Hill_n, half_effect_conc]
# print(all_params)
param_names = ['Vhalf', 'Kmax', 'Ku', 'N', 'EC50']

# Perturbing the parameters
# perturbing Vhalf
# parameter_interest = 'Vhalf'
# small_range = np.linspace(-220, -180, 2)
# cat_range = np.linspace(-100, -50, 10)
# cat_range2 = np.linspace(-15, -1, 10)
# param_range = np.concatenate((small_range, cat_range, cat_range2))
# param_range = small_range
# params_ranges = [param_range]

# Kmax
parameter_interest = 'Kmax'
interest_params = ['Kmax']
small_range = 10**np.linspace(0, np.log10(30), 10)
# [5.584e+01, 2.060e+05]
cat_range = 10**np.linspace(np.log10(50), np.log10(3e5), 10)
# [3.735e+07, 1.000e+08]
cat_range2 = 10**np.linspace(6, 9, 10)
param_range = np.concatenate((small_range, cat_range, cat_range2))
# params_ranges.append(param_range)
params_ranges = [cat_range]

# Load IKr model
model = '../../model/ohara-cipa-v1-2017-IKr.mmt'
model, _, x = myokit.load(model)

protocol_params = modelling.ProtocolParameters()
pulse_time = protocol_params.protocol_parameters['Milnes']['pulse_time']
protocol = protocol_params.protocol_parameters['Milnes']['function']

drug_model = modelling.BindingKinetics(model)
drug_model.protocol = protocol

# Set AP model
APmodel = '../../model/ohara-cipa-v1-2017.mmt'

# Load model
APmodel, _, x = myokit.load(APmodel)

AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')
pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)
base_conductance = APmodel.get('ikr.gKr').value()

offset = 50
repeats_AP = 800
save_signal = 2
repeats = 1000

orig_param_values = pd.DataFrame(all_params, index=param_names)
orig_param_values = orig_param_values.T

for k, parameter_interest in enumerate(interest_params):
    filename = 'SA_' + drug + '_' + parameter_interest + '.csv'
    if os.path.exists(saved_data_dir + filename):
        os.remove(saved_data_dir + filename)

    param_range = params_ranges[k]
    param_values = orig_param_values

    for num, param in enumerate(param_range):

        param_values[parameter_interest][0] = param

        drug_conc_Hill = list(np.append(0, 10**np.linspace(-1, 3, 7)))

        peaks = []
        for i in range(len(drug_conc_Hill)):
            log = drug_model.custom_simulation(
                param_values, drug_conc_Hill[i], repeats,
                log_var=['engine.time', 'ikr.IKr'])
            peak, _ = drug_model.extract_peak(log, 'ikr.IKr')
            peaks.append(peak[-1])

        peaks_norm = (peaks - min(peaks)) / (max(peaks) - min(peaks))

        # Make sure there are enough data points for the head of Hill curve
        checker_threshold = 0.95
        data_pt_checker = [True]
        while sum(data_pt_checker) < 3:
            drug_conc_Hill.insert(1, drug_conc_Hill[1] / np.sqrt(10))
            log = drug_model.custom_simulation(
                param_values, drug_conc_Hill[1], repeats,
                log_var=['engine.time', 'ikr.IKr'])
            peak, _ = drug_model.extract_peak(log, 'ikr.IKr')
            peaks.insert(1, peak[-1])
            peaks_norm = (peaks - min(peaks)) / (max(peaks) - min(peaks))
            data_pt_checker = [True if i > checker_threshold else False
                               for i in peaks_norm]

        # Make sure there are enough data points for the tail of Hill curve
        checker_threshold = 0.05
        data_pt_checker = [True]
        while sum(data_pt_checker) < 3:
            drug_conc_Hill = drug_conc_Hill + [max(drug_conc_Hill) *
                                               np.sqrt(10)]
            log = drug_model.custom_simulation(
                param_values, drug_conc_Hill[-1], repeats,
                log_var=['engine.time', 'ikr.IKr'])
            peak, _ = drug_model.extract_peak(log, 'ikr.IKr')
            peaks.append(peak[-1])
            peaks_norm = (peaks - min(peaks)) / (max(peaks) - min(peaks))
            data_pt_checker = [True if i < checker_threshold else False
                               for i in peaks_norm]

        # # Fit Hill curve
        Hill_model = modelling.HillsModel()
        optimiser = modelling.HillsModelOpt(Hill_model)
        Hill_curve, _ = optimiser.optimise(drug_conc_Hill, peaks_norm)

        if check_plot:
            max_grid = np.ceil(np.log(drug_conc_Hill[-1]))
            min_grid = np.floor(np.log(drug_conc_Hill[1]))
            conc_grid = np.linspace(min_grid, max_grid, 20)

            plt.figure(figsize=(4, 3))
            plt.plot(np.log(drug_conc_Hill[1:]), peaks_norm[1:], 'o',
                     label='peak current')
            plt.plot(conc_grid, Hill_model.simulate(Hill_curve[:2],
                     np.exp(conc_grid)), 'k', label='fitted Hill eq')
            plt.xlabel('Drug concentration (log)')
            plt.ylabel('Normalised peak current')
            plt.tight_layout()
            plt.savefig(saved_fig_dir + "Hill_" + drug + '_' +
                        parameter_interest + '_' + str(num) + ".pdf")

        APD_trapping = []
        APD_conductance = []
        drug_conc_AP = 10**np.linspace(np.log10(drug_conc_Hill[1]),
                                       np.log10(max(drug_conc_Hill)), 3)
        for i in range(len(drug_conc_AP)):
            # Run simulation for trapping model
            log = AP_model.custom_simulation(
                param_values, drug_conc_AP[i], repeats_AP, timestep=0.1,
                save_signal=save_signal,
                log_var=['engine.time', 'membrane.V'])

            # Compute APD90
            APD_trapping_pulse = []
            for pulse in range(save_signal):
                apd90 = AP_model.APD90(log['membrane.V', pulse], offset, 0.1)
                APD_trapping_pulse.append(apd90)
            APD_trapping.append(APD_trapping_pulse)

            # Run simulation for conductance model
            reduction_scale = Hill_model.simulate(Hill_curve[:2],
                                                  drug_conc_AP[i])
            d2 = AP_model.conductance_simulation(
                base_conductance * reduction_scale, repeats_AP, timestep=0.1,
                save_signal=save_signal, abs_tol=1e-6, rel_tol=1e-5,
                log_var=['engine.time', 'membrane.V'])

            # Compute APD90
            APD_conductance_pulse = []
            for pulse in range(save_signal):
                apd90 = AP_model.APD90(d2['membrane.V', pulse], offset, 0.1)
                APD_conductance_pulse.append(apd90)
            APD_conductance.append(APD_conductance_pulse)

        APD_trapping = [max(i) for i in APD_trapping]
        APD_conductance = [max(i) for i in APD_conductance]

        MSError = np.sum((np.array(APD_trapping) -
                          np.array(APD_conductance))**2) / len(APD_trapping)
        MSError = np.sqrt(MSError)

        conc_Hill_ind = ['conc_' + str(i) for i, _ in
                         enumerate(drug_conc_Hill)]
        conc_AP_ind = ['conc_' + str(i) for i, _ in enumerate(drug_conc_AP)]
        index_dict = {'drug_conc_Hill': conc_Hill_ind,
                      'peak_current': conc_Hill_ind,
                      'Hill_curve': ['Hill_coef', 'IC50'],
                      'param_values': param_names, 'drug_conc_AP': conc_AP_ind,
                      'APD_trapping': conc_AP_ind,
                      'APD_conductance': conc_AP_ind, 'MSE': ['MSE']}
        all_index = [(i, j) for i in index_dict.keys() for j in index_dict[i]]
        index = pd.MultiIndex.from_tuples(all_index)

        big_df = pd.DataFrame(
            drug_conc_Hill + list(peaks) + list(Hill_curve[:2]) +
            list(param_values.values[0]) + list(drug_conc_AP) + APD_trapping +
            APD_conductance + [MSError], index=index)

        if num != 0:
            previous_df = pd.read_csv(saved_data_dir + filename,
                                      header=[0, 1], index_col=[0],
                                      skipinitialspace=True)

            comb_df = pd.concat([previous_df, big_df.T])
        else:
            comb_df = big_df.T

        comb_df.to_csv(saved_data_dir + filename)

        os.system('cp ' + saved_data_dir + filename + ' ' +
                  saved_data_dir + filename[:-4] + '_copy.csv')
