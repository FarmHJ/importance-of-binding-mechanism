#
# Explore the parameter space of drug-related parameters (Vhalf, Kmax and Ku).
# Compute the APD90 differences between the AP-SD model and the AP-CS model
# for a given virtual drug.
#

import itertools
import myokit
import numpy as np
import os
import pandas as pd
import pints
import time

import modelling

# Define directory to save simulation data
data_filepath = '../simulation_data/parameter_space_exploration/'
if not os.path.exists(data_filepath):
    os.makedirs(data_filepath)

# Model directory
current_model_filepath = '../math_model/ohara-cipa-v1-2017-IKr-opt.mmt'
AP_model_filepath = '../math_model/ohara-cipa-v1-2017-opt.mmt'

# Load current model and set Milnes' protocol
model, _, x = myokit.load(current_model_filepath)
current_model = modelling.BindingKinetics(model)

protocol_params = modelling.ProtocolParameters()
protocol = protocol_params.protocol_parameters['Milnes']['function']
current_model.protocol = protocol

# Load AP model and set current protocol
APmodel, _, x = myokit.load(AP_model_filepath)
AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')
pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)

# Define constants for simulations
offset = 50
save_signal = 2
repeats = 1000
APD_points = 20

# Get name of parameters
SA_model = modelling.SensitivityAnalysis()
param_names = SA_model.param_names

# Define parameter space
Vhalf_fullrange = SA_model.param_explore_uniform('Vhalf')
Kmax_fullrange = SA_model.param_explore_uniform('Kmax')
Ku_fullrange = SA_model.param_explore_uniform('Ku')

starting_param_df = pd.DataFrame([1] * 5, index=param_names).T
ComparisonController = modelling.ModelComparison(starting_param_df)


def param_evaluation(param_values):

    # Define parameter values of virtual drug
    param_id = param_values['param_id'][0]
    param_values = param_values.drop(columns=['param_id'])

    ComparisonController.drug_param_values = param_values

    # Compute Hill curve of the virtual drug with the SD model
    Hill_curve_coefs, drug_conc_Hill, peaks_norm = \
        ComparisonController.compute_Hill(current_model,
                                          parallel=False)
    # parameters of Hill curve are based on normalised drug concentration

    # Define drug concentration range similar to the drug concentration used
    # to infer Hill curve
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
        try:
            # Simulate APs and APD90s of the AP-SD model and the AP-CS model
            APD_trapping, APD_conductance, drug_conc_AP = \
                ComparisonController.APD_sim(
                    AP_model, Hill_curve_coefs, drug_conc=drug_conc_AP,
                    EAD=True)

            # Calculate RMSD and MD of simulated APD90 of the two models
            RMSError = ComparisonController.compute_RMSE(APD_trapping,
                                                         APD_conductance)
            MAError = ComparisonController.compute_ME(APD_trapping,
                                                      APD_conductance)
        except myokit.SimulationError:
            APD_trapping = [float("Nan")] * APD_points
            APD_conductance = [float("Nan")] * APD_points
            RMSError = float("Nan")
            MAError = float("Nan")
            print('simulation error')

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
                  'ME': ['ME']}
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
# Since Hill coefficient, N, does not affect APD difference behaviour, it
# can be fixed at any value.
# For simplicity, let N = 1.

# Save defined parameter space or load previously saved parameter space
if not os.path.isdir(data_filepath + 'parameter_space/'):
    os.makedirs(data_filepath + 'parameter_space/')
sample_filepath = data_filepath + 'parameter_space/parameter_space.csv'

param_space = []
if os.path.exists(sample_filepath):
    param_values_df = pd.read_csv(sample_filepath,
                                  header=[0], index_col=[0],
                                  skipinitialspace=True)
    for i in range(len(param_values_df.index)):
        param_space.append(param_values_df.iloc[[i]])
else:
    counter = 0
    param_values_df = pd.DataFrame(columns=param_names)
    for Vhalf, Kmax, Ku in itertools.product(
            Vhalf_fullrange, Kmax_fullrange, Ku_fullrange):

        param_values = pd.DataFrame([counter, Vhalf, Kmax, Ku, 1, 1],
                                    index=['param_id'] + param_names)
        param_values = param_values.T
        param_space.append(param_values)
        param_values_df = pd.concat([param_values_df, param_values])
        counter += 1
    param_values_df.to_csv(sample_filepath)

# Set up variables for data saving
total_samples = len(param_space)
samples_per_save = 1000
samples_split_n = int(np.ceil(total_samples / samples_per_save))
total_saving_file_num = np.arange(samples_split_n)

# Determine completed simulations so that it is not repeated
file_prefix = 'SA_allparam_'
data_dir = data_filepath + 'SA_space/'
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)
evaluation_result_files = [f for f in os.listdir(data_dir) if
                           f.startswith(file_prefix)]
if len(evaluation_result_files) == 0:
    file_id_dict = {}
    for i in range(samples_split_n):
        file_id_dict[i] = param_values_df['param_id'].values[
            i * samples_per_save: (i + 1) * samples_per_save]
    saving_file_dict = {'file_num': total_saving_file_num,
                        'sample_id_each_file': file_id_dict}
else:
    # Check for missing results file
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
        saved_results_df = pd.read_csv(data_dir + file,
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

# Use PINTS' parallel evaluator to evaluate the APD90 difference for each
# virtual drugs in the parameter space
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

        if os.path.exists(data_dir + filename):
            combined_df = pd.read_csv(data_dir + filename,
                                      header=[0, 1], index_col=[0],
                                      skipinitialspace=True)
            for i in range(len(big_df)):
                combined_df = pd.concat([combined_df, big_df[i].T])
        else:
            combined_df = big_df[0].T
            for i in range(1, len(big_df)):
                combined_df = pd.concat([combined_df, big_df[i].T])

        combined_df.to_csv(data_dir + filename)
