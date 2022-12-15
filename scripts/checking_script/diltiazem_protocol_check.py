import matplotlib.pyplot as plt
import myokit
import numpy as np
import pandas as pd

import modelling

# saved_data_dir = '../../simulation_data/binding_kinetics_comparison/diltiazem/'

# Load IKr model
model = '../../math_model/ohara-cipa-v1-2017-IKr.mmt'
model, _, x = myokit.load(model)

protocol_params = modelling.ProtocolParameters()
# protocol_name = 'Milnes'
protocol_list = protocol_params.protocols
protocol_name = 'P40'
pulse_time = protocol_params.protocol_parameters[protocol_name]['pulse_time']
protocol = protocol_params.protocol_parameters[protocol_name]['function']

drug_model = modelling.BindingKinetics(model)
drug_model.protocol = protocol

param_lib = modelling.BindingParameters()
drug_list = ['diltiazem']

SA_model = modelling.SensitivityAnalysis()
param_names = SA_model.param_names

starting_param_df = pd.DataFrame([1] * 5, index=param_names).T
ComparisonController = modelling.ModelComparison(starting_param_df)


def param_evaluation(param_values):

    # orig_half_effect_conc = param_values['EC50'][0]
    # param_values[parameter_interest][0] = param
    # param_values['EC50'][0] = 1
    ComparisonController.drug_param_values = param_values

    # Hill_n = param_values['N'][0]
    # norm_constant = np.power(orig_half_effect_conc, 1 / Hill_n)

    Hill_curve_coefs, drug_conc_Hill, peaks_norm = \
        ComparisonController.compute_Hill(drug_model,) 
                                        #   norm_constant=norm_constant)

    return Hill_curve_coefs, drug_conc_Hill, peaks_norm

drug = drug_list[0]
Vhalf = param_lib.binding_parameters[drug]['Vhalf']
Kmax = param_lib.binding_parameters[drug]['Kmax']
Ku = param_lib.binding_parameters[drug]['Ku']
N = param_lib.binding_parameters[drug]['N']
EC50 = param_lib.binding_parameters[drug]['EC50']

paramid = 0
param_values = pd.DataFrame([paramid, Vhalf, Kmax, Ku, N, EC50],
                            index=['param_id'] + param_names)
param_values = param_values.T

plt.figure()
Hill_coef_df = pd.DataFrame(columns=['Hill coefficient', 'IC50', 'protocol'])
for prot in protocol_list:
    # pulse_time = protocol_params.protocol_parameters[prot]['pulse_time']
    protocol = protocol_params.protocol_parameters[prot]['function']
    drug_model.protocol = protocol

    Hill_curve_coefs, drug_conc_Hill, peaks_norm = \
        param_evaluation(param_values)
    print(Hill_curve_coefs)
    Hill_df = pd.DataFrame({'Hill coefficient': [Hill_curve_coefs[0]],
                            'IC50': [Hill_curve_coefs[1]],
                            'protocol': [prot]})
    Hill_coef_df = pd.concat([Hill_coef_df, Hill_df])

# Hill_coef_df.to_csv(saved_data_dir + 'Hill_curves.csv')

    # Hill_model = modelling.HillsModel()
    # reduction_scale = Hill_model.simulate(Hill_curve_coefs, drug_conc_Hill)

# Hill_coef_df = pd.read_csv(saved_data_dir + 'Hill_curves.csv',
#                            header=[0], index_col=[0])
# Hill_coef_df.iloc[[3]] = Hill_curve_coefs + ['P40']

# print(Hill_coef_df)

    plt.plot(drug_conc_Hill, peaks_norm, 'o')
# plt.plot(drug_conc_Hill, reduction_scale, '-')
plt.xscale('log')
plt.savefig('../../testing_figures/testing/Hill_error_diltiazem.pdf')

# combined_result_df = pd.concat([combined_result_df, result_df.T])

# combined_result_df.to_csv(saved_data_dir + filename)
