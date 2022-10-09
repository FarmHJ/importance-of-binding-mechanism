import pandas as pd
from SALib.sample import saltelli
from SALib.analyze import sobol
import myokit
import numpy as np
import os
import pints
import time

import modelling


param_names = ['Vhalf', 'Kmax', 'Ku', 'N', 'EC50']

saved_data_filepath = os.path.dirname(__file__)

# Model directory
model_dir = os.path.join(saved_data_filepath, 'model')
current_model_filename = 'ohara-cipa-v1-2017-IKr.mmt'
current_model_filepath = os.path.join(model_dir, current_model_filename)
AP_model_filename = 'ohara-cipa-v1-2017.mmt'
AP_model_filepath = os.path.join(model_dir, AP_model_filename)

# Load models
model, _, x = myokit.load(current_model_filepath)

protocol_params = modelling.ProtocolParameters()
protocol = protocol_params.protocol_parameters['Milnes']['function']

BKmodel = modelling.BindingKinetics(model)
BKmodel.protocol = protocol

APmodel, _, x = myokit.load(AP_model_filepath)
AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')
pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)

Hill_model = modelling.HillsModel()

init_param_values = pd.DataFrame([0] * 5, index=param_names)
init_param_values = init_param_values.T
ComparisonController = modelling.ModelComparison(init_param_values)


def model_comparison(param_values):

    orig_param_values = pd.DataFrame(param_values, index=param_names)
    orig_param_values = orig_param_values.T

    # Log transformation
    orig_param_values['Kmax'] = 10**orig_param_values['Kmax']
    orig_param_values['Ku'] = 10**orig_param_values['Ku']
    orig_param_values['EC50'] = 10**orig_param_values['EC50']

    start_time = time.time()
    ComparisonController.drug_param_values = orig_param_values
    Hill_curve_coefs, drug_conc_Hill, _ = \
        ComparisonController.compute_Hill(BKmodel)
    evaluation_Hill_time = time.time() - start_time

    if isinstance(Hill_curve_coefs, str):
        return float("nan"), float("nan"), np.inf, 0

    drug_conc_AP = 10**np.linspace(np.log10(drug_conc_Hill[1]),
                                   np.log10(max(drug_conc_Hill)),
                                   20)
    start_time = time.time()
    try:
        APD_trapping, APD_conductance = ComparisonController.APD_sim(
            AP_model, Hill_curve_coefs, drug_conc=drug_conc_AP)
        RMSError = ComparisonController.compute_RMSE(APD_trapping,
                                                     APD_conductance)
        MAError = ComparisonController.compute_MAE(APD_trapping,
                                                   APD_conductance)

    except myokit.SimulationError:
        RMSError = float("Nan")
        MAError = float("nan")
    evaluation_RMSE_time = time.time() - start_time

    return RMSError, MAError, evaluation_Hill_time, evaluation_RMSE_time


# Define the model inputs
problem = {
    'num_vars': 5,
    'names': param_names,
    'bounds': [[-210, -0.1],
               [0, 9],
               [-5, 0],
               [0.5, 1.5],
               [0, 10]],
}

# Generate samples
sample_filepath = os.path.join(saved_data_filepath, "param_value_samples.txt")
if os.path.exists(sample_filepath):
    param_values = np.loadtxt(sample_filepath)
else:
    samples_n = 1024
    param_values = saltelli.sample(problem, samples_n)
    np.savetxt(os.path.join(saved_data_filepath, "param_value_samples.txt"),
               param_values)

# Set up frequency to save files
total_samples = param_values.shape[0]
samples_per_save = 1000
samples_split_n = int(np.ceil(total_samples / samples_per_save))

# Find out completed and saved function evaluations
data_dir = os.getcwd()
evaluation_result_files = [f for f in os.listdir(data_dir) if
                           f.startswith('MSError_evaluations_')]
if len(evaluation_result_files) == 0:
    saving_file_num = np.arange(samples_split_n)
else:
    result_files_num = [int(fname[20:-4]) for fname in
                        evaluation_result_files]
    if len(result_files_num) == max(result_files_num) + 1:
        saving_file_num = np.arange(max(result_files_num) + 1, samples_split_n)
    else:
        missing_file_num = [x for x in range(result_files_num[-1] + 1)
                            if x not in result_files_num]
        saving_file_num = np.arange(max(result_files_num) + 1, samples_split_n)
        saving_file_num = np.concatenate((np.array(missing_file_num),
                                          saving_file_num))

# Set up parallel evaluator
n_workers = 40
evaluator = pints.ParallelEvaluator(model_comparison, n_workers=n_workers)

# Evaluate function
for i in saving_file_num:
    print('Starting function evaluation for file number: ', i)
    subset_param_values = param_values[
        samples_per_save * i:samples_per_save * (i + 1)]
    print('Samples ', samples_per_save * i, 'to', samples_per_save * (i + 1) - 1)
    Y = []
    for log_num in range(int(samples_per_save / (n_workers * 5))):
        ind_multiplier = n_workers * 5

        print('Running evaluation for sample number ',
              ind_multiplier * log_num + samples_per_save * i, ' to ',
              ind_multiplier * (log_num + 1) + samples_per_save * i - 1)
        current_time = time.strftime("%H:%M:%S", time.localtime())
        print('Starting time: ', current_time)

        subsubset_param_values = subset_param_values[
            ind_multiplier * log_num:ind_multiplier * (log_num + 1)]
        output = evaluator.evaluate(subsubset_param_values)
        Y += output
    np.savetxt(os.path.join(
        saved_data_filepath, "MSError_evaluations_" + str(i) + ".txt"),
        np.array(Y))

print('All function evaluations completed.')
print('Loading all results of function evaluations.')
# Load all function evaluations
evaluation_result_files = [f for f in os.listdir(data_dir) if
                           f.startswith('MSError_evaluations_')]
result_files_num = [int(fname[20:-4]) for fname in evaluation_result_files]
sort_ind = [i[0] for i in sorted(enumerate(result_files_num), key=lambda x:x[1])]
sorted_result_files = [evaluation_result_files[i] for i in sort_ind]

first_iter = True
for result_file in sorted_result_files:
    if first_iter:
        Y = np.loadtxt(result_file)
        first_iter = False
    else:
        loaded_result = np.loadtxt(result_file)
        Y = np.concatenate((Y, loaded_result))

print('Performing analysis.')
# Perform analysis
Si = sobol.analyze(problem, Y[:, 0], parallel=True, n_processors=n_workers)

print('Saving data.')
# Save data
total_Si, first_Si, second_Si = Si.to_df()

total_Si.to_csv(os.path.join(saved_data_filepath, 'total_Si.csv'))
first_Si.to_csv(os.path.join(saved_data_filepath, 'first_Si.csv'))
second_Si.to_csv(os.path.join(saved_data_filepath, 'second_Si.csv'))
