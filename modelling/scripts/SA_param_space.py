# To plot the range of parameter values
# Drug binding-related parameters

import itertools
import myokit
import numpy as np
import os
import pandas as pd
import pints
import time

import modelling

param_space_dir = '../../simulation_results/parameter_space/'
saved_data_dir = '../../simulation_results/SA_curve/'

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

offset = 50
repeats_AP = 800
save_signal = 2
repeats = 1000
APD_points = 20

param_lib = modelling.BindingParameters()
drug_list = param_lib.drug_compounds

SA_model = modelling.SensitivityAnalysis()
param_names = SA_model.param_names
# res = 5
# Vhalf_range = SA_model.param_explore('Vhalf', res)
# Ku_range = SA_model.param_explore('Ku', res)
# Kmax_range = SA_model.param_explore('Kmax', res)
res = 5
Vhalf_fullrange = SA_model.param_explore('Vhalf', res)
Vhalf_fullrange = SA_model.param_explore_gaps(Vhalf_fullrange, 3, 'Vhalf')
# Vhalf_done = sorted(Vhalf_fullrange)[::5]
# Vhalf_done = sorted(Vhalf_fullrange)[15::2]
# Vhalf_range = sorted([Vhalf_fullrange[i] for i in range(len(Vhalf_fullrange))
#                       if Vhalf_fullrange[i] not in Vhalf_done])
Vhalf_range = sorted(Vhalf_fullrange)

# print(Vhalf_range)
# # remove -100.6775 as it is close to -100
Vhalf_range.pop(11)
# Vhalf_range = np.array(Vhalf_range)

# Kmax_range = np.concatenate((np.linspace(1.11078, 1.4771, 3 + 2)[1:-1],
#                              np.linspace(1.4771, 2.6387, 7 + 2)[1:-1]))
# Kmax_range = 10**Kmax_range

# Ku_range = np.concatenate((np.linspace(-3.924, -3.1009, res + 2)[1:-1],
#                            np.linspace(-4.7471, -3.924, 7 + 2)[1:-1]))
# Ku_range = 10**Ku_range

Kmax_fullrange = SA_model.param_explore('Kmax', res)
Kmax_fullrange = SA_model.param_explore_gaps(Kmax_fullrange, 3, 'Kmax')
Kmax_fullrange = np.log10(sorted(Kmax_fullrange))

Ku_fullrange = SA_model.param_explore('Ku', res)
Ku_fullrange = SA_model.param_explore_gaps(Ku_fullrange, 3, 'Ku')
Ku_fullrange = np.log10(sorted(Ku_fullrange))

Kmax_range1 = np.linspace(Kmax_fullrange[3], Kmax_fullrange[4], 3 + 2)[1:-1]
Kmax_range1 = 10**Kmax_range1

Ku_range1 = np.concatenate((np.linspace(Ku_fullrange[4], Ku_fullrange[7], 9),
                            Ku_fullrange[7:]))
Ku_range1 = 10**Ku_range1

Kmax_range2 = Kmax_fullrange[6:]
Kmax_range2 = 10**Kmax_range2

Ku_range2 = np.linspace(Ku_fullrange[0], Ku_fullrange[1], 3 + 2)[1:-1]
Ku_range2 = 10**Ku_range2

Kmax_range3 = Kmax_fullrange[6:]
Kmax_range3 = 10**Kmax_range3
Ku_range3 = np.linspace(-5, Ku_fullrange[0], 3 + 2)[:-1]
Ku_range3 = 10**Ku_range3

starting_param_df = pd.DataFrame([1] * 5, index=param_names).T
ComparisonController = modelling.ModelComparison(starting_param_df)


def param_evaluation(param_values):

    param_id = param_values['param_id'][0]
    param_values = param_values.drop(columns=['param_id'])

    ComparisonController.drug_param_values = param_values

    Hill_curve_coefs, drug_conc_Hill, peaks_norm = \
        ComparisonController.compute_Hill(drug_model,
                                          parallel=False)
    # parameters of Hill's curve are based on normalised drug concentration

    drug_conc_AP = 10**np.linspace(np.log10(drug_conc_Hill[1]),
                                   np.log10(max(drug_conc_Hill)),
                                   APD_points)

    if isinstance(Hill_curve_coefs, str):
        Hill_curve_coefs = [float("nan")] * 2
        APD_trapping = [float("Nan")] * APD_points
        APD_conductance = [float("Nan")] * APD_points
        RMSError = float("Nan")
        MAError = float("Nan")
    else:
        # Simulate action potentials
        try:
            APD_trapping, APD_conductance, drug_conc_AP = \
                ComparisonController.APD_sim(
                    AP_model, Hill_curve_coefs, drug_conc=drug_conc_AP,
                    EAD=True)

            RMSError = ComparisonController.compute_RMSE(APD_trapping,
                                                         APD_conductance)
            MAError = ComparisonController.compute_MAE(APD_trapping,
                                                       APD_conductance)
        except myokit.SimulationError:
            APD_trapping = [float("Nan")] * APD_points
            APD_conductance = [float("Nan")] * APD_points
            RMSError = float("Nan")
            MAError = float("Nan")

    # Create dataframe to save results
    conc_Hill_ind = ['conc_' + str(i) for i, _ in
                     enumerate(drug_conc_Hill)]
    conc_AP_ind = ['conc_' + str(i) for i, _ in enumerate(drug_conc_AP)]
    index_dict = {'param_id': ['param_id'],
                  'drug_conc_Hill': conc_Hill_ind,
                  'peak_current': conc_Hill_ind,
                  'Hill_curve': ['Hill_coef', 'IC50'],
                  'param_values': param_names, 'drug_conc_AP': conc_AP_ind,
                  'APD_trapping': conc_AP_ind,
                  'APD_conductance': conc_AP_ind, 'RMSE': ['RMSE'],
                  'MAE': ['MAE']}
    all_index = [(i, j) for i in index_dict.keys() for j in index_dict[i]]
    index = pd.MultiIndex.from_tuples(all_index)

    big_df = pd.DataFrame(
        [param_id] + drug_conc_Hill + list(peaks_norm) +
        list(Hill_curve_coefs) + list(param_values.values[0]) +
        list(drug_conc_AP) + APD_trapping + APD_conductance +
        [RMSError] + [MAError], index=index)

    return big_df


# Assuming drug concentration are all normalised, the EC50 value in the model
# becomes 1.
# Since Hill's coefficient, N, does not affect APD difference behaviour, it
# can be fixed at any value.
# For simplicity, let N = 1.

# sample_filepath = saved_data_dir + 'parameter_space_res5.csv'
# for curve
# simulations for parameters in parameter_space_curve2.csv not completed
for curve_num in ['curve2', 'curve3']:
    sample_filepath = param_space_dir + 'parameter_space_' + curve_num + '.csv'
    print('Running for parameter space file number : ', curve_num)
    param_space = []
    # if os.path.exists(sample_filepath):
    param_values_df = pd.read_csv(sample_filepath,
                                  header=[0], index_col=[0],
                                  skipinitialspace=True)
    for i in range(len(param_values_df.index)):
        param_space.append(param_values_df.iloc[[i]])
    # else:
    #     # counter = 0
    #     # for curve
    #     counter = 20000
    #     param_values_df = pd.DataFrame(columns=param_names)
    #     for Kmax_range_i, Ku_range_i in zip([Kmax_range1, Kmax_range3], [Ku_range1, Ku_range3]):
    #         for Vhalf, Kmax, Ku in itertools.product(
    #                 Vhalf_range, Kmax_range_i, Ku_range_i):

    #             param_values = pd.DataFrame([counter, Vhalf, Kmax, Ku, 1, 1],
    #                                         index=['param_id'] + param_names)
    #             param_values = param_values.T
    #             param_space.append(param_values)
    #             param_values_df = pd.concat([param_values_df, param_values])
    #             counter += 1
    #     # for curve
    #     param_values_df.to_csv(param_space_dir + 'parameter_space_curve3.csv')

    total_samples = len(param_space)
    samples_per_save = 1000
    samples_split_n = int(np.ceil(total_samples / samples_per_save))
    total_saving_file_num = np.arange(samples_split_n)

    file_prefix = 'SA_' + curve_num + '_'
    evaluation_result_files = [f for f in os.listdir(saved_data_dir) if
                               f.startswith(file_prefix)]
    if len(evaluation_result_files) == 0:
        file_id_dict = {}
        for i in range(samples_split_n):
            file_id_dict[i] = param_values_df['param_id'].values[
                i * samples_per_save: (i + 1) * samples_per_save]
        saving_file_dict = {'file_num': total_saving_file_num,
                            'sample_id_each_file': file_id_dict}
    else:
        # to include missing file
        result_files_num = [int(fname[len(file_prefix):-4]) for fname in
                            evaluation_result_files]
        file_num_to_run = []
        file_id_dict = {}
        missing_file = [i for i in total_saving_file_num
                        if i not in result_files_num]
        for i in missing_file:
            file_id_dict[i] = param_values_df['param_id'].values[
                i * samples_per_save: (i + 1) * samples_per_save]
            file_num_to_run.append(i)
        for file in evaluation_result_files:
            file_num = int(file[len(file_prefix):-4])
            saved_results_df = pd.read_csv(saved_data_dir + file,
                                           header=[0, 1], index_col=[0],
                                           skipinitialspace=True)
            ran_values = saved_results_df['param_id']['param_id'].values
            expected_ids = param_values_df['param_id'].values[
                file_num * samples_per_save: (file_num + 1) * samples_per_save]
            param_space_id = [i for i in expected_ids if i not in ran_values]
            if len(param_space) != 0:
                file_num_to_run.append(file_num)
                file_id_dict[file_num] = param_space_id
        saving_file_dict = {'file_num': sorted(file_num_to_run),
                            'sample_id_each_file': file_id_dict}

    n_workers = 8
    evaluator = pints.ParallelEvaluator(param_evaluation,
                                        n_workers=n_workers)
    for file_num in saving_file_dict['file_num']:
        print('Starting function evaluation for file number: ', file_num)
        current_time = time.strftime("%H:%M:%S", time.localtime())
        print('Starting time: ', current_time)
        samples_to_run = saving_file_dict['sample_id_each_file'][file_num]
        samples_num = len(samples_to_run)
        filename = file_prefix + str(file_num) + '.csv'

        for i in range(int(np.ceil(samples_num / n_workers))):
            subset_samples_to_run = samples_to_run[
                n_workers * i:n_workers * (i + 1)]
            print('Running samples ', int(subset_samples_to_run[0]), ' to ',
              int(subset_samples_to_run[-1]))
            subset_param_space = param_values_df.loc[
                param_values_df['param_id'].isin(subset_samples_to_run)]
            param_space = []
            for i in range(len(subset_param_space.index)):
                param_space.append(subset_param_space.iloc[[i]])

            big_df = evaluator.evaluate(param_space)

            if os.path.exists(saved_data_dir + filename):
                combined_df = pd.read_csv(saved_data_dir + filename,
                                          header=[0, 1], index_col=[0],
                                          skipinitialspace=True)
                for i in range(len(big_df)):
                    combined_df = pd.concat([combined_df, big_df[i].T])
            else:
                combined_df = big_df[0].T
                for i in range(1, len(big_df)):
                    combined_df = pd.concat([combined_df, big_df[i].T])

            combined_df.to_csv(saved_data_dir + filename)
