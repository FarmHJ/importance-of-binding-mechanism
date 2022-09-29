import pandas as pd
from SALib.sample import saltelli
from SALib.analyze import sobol
import matplotlib.pyplot as plt
import myokit
import numpy as np
import os

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


def model_comparison(param_values):
    orig_param_values = pd.DataFrame(param_values, index=param_names)
    orig_param_values = orig_param_values.T

    ComparisonController = modelling.ModelComparison(orig_param_values)
    Hill_curve_coefs, drug_conc_Hill, peaks_norm = \
        ComparisonController.compute_Hill(BKmodel)

    max_grid = np.ceil(np.log(drug_conc_Hill[-1]))
    min_grid = np.floor(np.log(drug_conc_Hill[1]))
    conc_grid = np.linspace(min_grid, max_grid, 20)

    plt.figure(figsize=(4, 3))
    plt.plot(np.log(drug_conc_Hill[1:]), peaks_norm[1:], 'o',
             label='peak current')
    plt.plot(conc_grid, Hill_model.simulate(Hill_curve_coefs,
             np.exp(conc_grid)), 'k', label='fitted Hill eq')
    plt.xlabel('Drug concentration (log)')
    plt.ylabel('Normalised peak current')
    plt.tight_layout()
    plt.savefig(os.path.join(saved_fig_filepath, "Hill_curve_check.pdf"))

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
               [1, 1e9],
               [1e-5, 1],
               [0.5, 1.5],
               [1, 1e10]],
    # 'bounds': [[-210, -0.1],
    #            [0, 9],
    #            [-5, 0],
    #            [0.5, 1.5],
    #            [0, 10]],
    # 'dists': ['unif', 'lognorm', 'lognorm', 'unif', 'lognorm']
}

# Generate samples
param_values = saltelli.sample(problem, 8)

# Run model (example)
Y = np.zeros([param_values.shape[0]])

for i, X in enumerate(param_values):
    print('evaluating...')
    Y[i] = model_comparison(X)

# Perform analysis
Si = sobol.analyze(problem, Y, print_to_console=True, parallel=True,)
#                    n_processors=)

# Save data
total_Si, first_Si, second_Si = Si.to_df()
total_Si.to_csv(os.path.join(saved_data_filepath, 'total_Si.csv'))
first_Si.to_csv(os.path.join(saved_data_filepath, 'first_Si.csv'))
second_Si.to_csv(os.path.join(saved_data_filepath, 'second_Si.csv'))
