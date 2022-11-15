# To plot the range of parameter values
# Drug binding-related parameters

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
    # The parameters of Hill's curve are based on the normalised
    # drug concentration
    # Hill's coefficient remains the same but IC50 -> IC50/EC50

    # drug_conc_AP = 10**np.linspace(np.log10(drug_conc_Hill[1]),
    #                                np.log10(max(drug_conc_Hill)),
    #                                APD_points)

    # if isinstance(Hill_curve_coefs, str):
    #     Hill_curve_coefs = [float("nan")] * 2
    #     APD_trapping = [float("Nan")] * APD_points
    #     APD_conductance = [float("Nan")] * APD_points
    #     RMSError = float("Nan")
    #     MAError = float("Nan")
    # else:
    #     # Simulate action potentials
    #     try:
    #         APD_trapping, APD_conductance, drug_conc_AP = \
    #             ComparisonController.APD_sim(
    #                 AP_model, Hill_curve_coefs, drug_conc=drug_conc_AP,
    #                 EAD=True)

    #         RMSError = ComparisonController.compute_RMSE(APD_trapping,
    #                                                      APD_conductance)
    #         MAError = ComparisonController.compute_ME(APD_trapping,
    #                                                   APD_conductance)
    #     except myokit.SimulationError:
    #         APD_trapping = [float("Nan")] * APD_points
    #         APD_conductance = [float("Nan")] * APD_points
    #         RMSError = float("Nan")
    #         MAError = float("Nan")

    # # Create dataframe to save results
    # conc_Hill_ind = ['conc_' + str(i) for i, _ in
    #                  enumerate(drug_conc_Hill)]
    # conc_AP_ind = ['conc_' + str(i) for i, _ in enumerate(drug_conc_AP)]
    # index_dict = {'drug_conc_Hill': conc_Hill_ind,
    #               'peak_current': conc_Hill_ind,
    #               'Hill_curve': ['Hill_coef', 'IC50'],
    #               'param_values': param_names, 'drug_conc_AP': conc_AP_ind,
    #               'APD_trapping': conc_AP_ind,
    #               'APD_conductance': conc_AP_ind, 'RMSE': ['RMSE'],
    #               'ME': ['ME']}
    # all_index = [(i, j) for i in index_dict.keys() for j in index_dict[i]]
    # index = pd.MultiIndex.from_tuples(all_index)

    # param_values['EC50'][0] = orig_half_effect_conc
    # big_df = pd.DataFrame(
    #     drug_conc_Hill + list(peaks_norm) + list(Hill_curve_coefs) +
    #     list(param_values.values[0]) + list(drug_conc_AP) + APD_trapping +
    #     APD_conductance + [RMSError] + [MAError], index=index)

    # return big_df
    return Hill_curve_coefs, drug_conc_Hill, peaks_norm


filename = 'SA_mexiletine_N_test.csv'
result_df = pd.read_csv(saved_data_dir + filename,
                        header=[0, 1], index_col=[0],
                        skipinitialspace=True)

Vhalf = result_df['param_values']['Vhalf'].values[0]
Kmax = result_df['param_values']['Kmax'].values[0]
Ku = result_df['param_values']['Ku'].values[0]
N = result_df['param_values']['N'].values[0]
EC = result_df['param_values']['EC50'].values[0]
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
