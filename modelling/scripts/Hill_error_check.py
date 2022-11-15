# Check why peak of hERG current does not follow sigmoid curve
# Ans: the drug concentration range I was looking at are either too high or
# too low that I was focusing only on the part close to one or zero, thus the
# fluctuate is large.

import matplotlib.pyplot as plt
import myokit
import numpy as np
import pandas as pd

import modelling

saved_data_dir = '../../simulation_data/sensitivity_analysis/N/'

# Load IKr model
model = '../../model/ohara-cipa-v1-2017-IKr.mmt'
model, _, x = myokit.load(model)

protocol_params = modelling.ProtocolParameters()
pulse_time = protocol_params.protocol_parameters['Milnes']['pulse_time']
protocol = protocol_params.protocol_parameters['Milnes']['function']

drug_model = modelling.BindingKinetics(model)
drug_model.protocol = protocol

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
# repeats = 1000
# APD_points = 20

param_lib = modelling.BindingParameters()
drug_list = ['mexiletine']
SA_model = modelling.SensitivityAnalysis()
param_names = SA_model.param_names
parameter_interest = 'N'

starting_param_df = pd.DataFrame([1] * 5, index=param_names).T
ComparisonController = modelling.ModelComparison(starting_param_df)


def param_evaluation(param_values):

    orig_half_effect_conc = param_values['EC50'][0]
    # param_values[parameter_interest][0] = param
    param_values['EC50'][0] = 1
    ComparisonController.drug_param_values = param_values

    Hill_n = param_values['N'][0]
    norm_constant = np.power(orig_half_effect_conc, 1 / Hill_n)

    Hill_curve_coefs, drug_conc_Hill, peaks_norm = \
        ComparisonController.compute_Hill(drug_model,
                                          norm_constant=norm_constant)

    return Hill_curve_coefs, drug_conc_Hill, peaks_norm


filename = 'SA_mexiletine_N_test.csv'
result_df = pd.read_csv(saved_data_dir + filename,
                        header=[0, 1], index_col=[0],
                        skipinitialspace=True)

Hill_coefs = result_df['Hill_curve']['Hill_coef'].values
nan_ind = [i for i in range(len(Hill_coefs)) if np.isnan(Hill_coefs[i])]

Vhalf = result_df['param_values']['Vhalf'].values[nan_ind[0]]
Kmax = result_df['param_values']['Kmax'].values[nan_ind[0]]
Ku = result_df['param_values']['Ku'].values[nan_ind[0]]
N = result_df['param_values']['N'].values[nan_ind[0]]
EC = result_df['param_values']['EC50'].values[nan_ind[0]]
paramid = 0
param_values = pd.DataFrame([paramid, Vhalf, Kmax, Ku, N, EC],
                            index=['param_id'] + param_names)
param_values = param_values.T
Hill_curve_coefs, drug_conc_Hill, peaks_norm = param_evaluation(param_values)

print(Hill_curve_coefs)

plt.figure()
plt.plot(drug_conc_Hill, peaks_norm, 'o')
plt.xscale('log')
plt.savefig('../../figures/testing/Hill_error.pdf')

# combined_result_df = pd.concat([combined_result_df, result_df.T])

# combined_result_df.to_csv(saved_data_dir + filename)
