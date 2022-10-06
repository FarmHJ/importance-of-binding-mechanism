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
model_dir = os.path.join(os.path.dirname(__file__), 'model')
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
        return float("nan"), np.inf, 0

    drug_conc_range = (np.log10(drug_conc_Hill[1]),
                       np.log10(drug_conc_Hill[-1]))
    start_time = time.time()
    try:
        RMSError = ComparisonController.compute_RMSE(
            AP_model, Hill_curve_coefs, drug_conc=drug_conc_range)
    except myokit.SimulationError:
        RMSError = (float("Nan"), 0, 0)
    evaluation_RMSE_time = time.time() - start_time

    if np.isnan(RMSError[0]):
        return RMSError[0], evaluation_Hill_time, evaluation_RMSE_time
    else:
        return RMSError[0], evaluation_Hill_time, evaluation_RMSE_time


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
samples_n = 1024
param_values = saltelli.sample(problem, samples_n)
np.savetxt(os.path.join(saved_data_filepath, "param_value_samples.txt"), param_values)

# Evaluate function
evaluator = pints.ParallelEvaluator(model_comparison, n_workers=35)
output = evaluator.evaluate(param_values)
Y = np.array(output)
np.savetxt(os.path.join(saved_data_filepath, "MSError_evaluations.txt"), Y)

# Perform analysis
Si = sobol.analyze(problem, Y[:, 0], parallel=True, n_processors=35)

# Save data
total_Si, first_Si, second_Si = Si.to_df()

total_Si.to_csv(os.path.join(saved_data_filepath, 'total_Si.csv'))
first_Si.to_csv(os.path.join(saved_data_filepath, 'first_Si.csv'))
second_Si.to_csv(os.path.join(saved_data_filepath, 'second_Si.csv'))

