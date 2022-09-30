import pandas as pd
from SALib.sample import saltelli
from SALib.analyze import sobol
import matplotlib.pyplot as plt
import myokit
import numpy as np
import os
import pints

import modelling


param_names = ['Vhalf', 'Kmax', 'Ku', 'N', 'EC50']

data_dir = os.path.join(os.path.dirname(__file__), 'simulation_data')
saved_data_filepath = os.path.join(data_dir, 'sensitivity_analysis')

if not os.path.exists(saved_data_filepath):
    os.makedirs(saved_data_filepath)

fig_dir = os.path.join(os.path.dirname(__file__), 'figures')
saved_fig_filepath = os.path.join(fig_dir, 'sensitivity_analysis')

if not os.path.exists(saved_fig_filepath):
    os.makedirs(saved_fig_filepath)

# Model directory
model_dir = os.path.join(os.path.dirname(__file__), 'model')
current_model_filename = 'ohara-cipa-v1-2017-IKr.mmt'
current_model_filepath = os.path.join(model_dir, current_model_filename)
AP_model_filename = 'ohara-cipa-v1-2017.mmt'
AP_model_filepath = os.path.join(model_dir, AP_model_filename)

# Load models
model, _, x = myokit.load(current_model_filepath)

protocol_params = modelling.ProtocolParameters()
pulse_time = protocol_params.protocol_parameters['Milnes']['pulse_time']
protocol = protocol_params.protocol_parameters['Milnes']['function']

BKmodel = modelling.BindingKinetics(model)
BKmodel.protocol = protocol

APmodel, _, x = myokit.load(AP_model_filepath)
AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')
pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)

Hill_model = modelling.HillsModel()

# TODO: Turn off parallelisation of PINTS default
def model_comparison(param_values):
    orig_param_values = pd.DataFrame(param_values, index=param_names)
    orig_param_values = orig_param_values.T

    # Log transformation
    orig_param_values['Kmax'] = 10**orig_param_values['Kmax']
    orig_param_values['Ku'] = 10**orig_param_values['Ku']
    orig_param_values['EC50'] = 10**orig_param_values['EC50']

    ComparisonController = modelling.ModelComparison(orig_param_values)
    Hill_curve_coefs, drug_conc_Hill, _ = \
        ComparisonController.compute_Hill(BKmodel)

    drug_conc_range = (np.log10(drug_conc_Hill[1]),
                       np.log10(drug_conc_Hill[-1]))
    MSError, _, _ = ComparisonController.compute_MSE(
        AP_model, Hill_curve_coefs, drug_conc=drug_conc_range)

    return MSError


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

# TODO: Estimate time usage
# Generate samples
param_values = saltelli.sample(problem, 4)

# Run model (example)
Y = np.zeros([param_values.shape[0]])

# TODO: parallelise
for i, X in enumerate(param_values):
    print('evaluating...')
    Y[i] = model_comparison(X)
    print(Y)

# evaluator = pints.ParallelEvaluator(model_comparison)
# Y = evaluator.evaluate(param_values)

# Perform analysis
Si = sobol.analyze(problem, Y, print_to_console=True, parallel=True,)
#                    n_processors=)

# Save data
total_Si, first_Si, second_Si = Si.to_df()
total_Si.to_csv(os.path.join(saved_data_filepath, 'total_Si.csv'))
first_Si.to_csv(os.path.join(saved_data_filepath, 'first_Si.csv'))
second_Si.to_csv(os.path.join(saved_data_filepath, 'second_Si.csv'))
