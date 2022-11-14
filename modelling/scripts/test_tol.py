# To plot the range of parameter values
# Drug binding-related parameters

import matplotlib.pyplot as plt
import myokit
import numpy as np
import os
import pandas as pd

import modelling

param_space_dir = '../../simulation_results/parameter_space/'
# saved_data_dir = '../../simulation_results/SA_curve/'
saved_data_dir = '../../simulation_results/'
file_prefix = 'SA_APD'
result_files = [saved_data_dir + f for f in os.listdir(saved_data_dir)
                if f.startswith(file_prefix)]

first_iter = True
for file in result_files:
    df = pd.read_csv(file,
                     header=[0, 1], index_col=[0],
                     skipinitialspace=True)

    chosen_df = df.loc[np.isnan(df['RMSE']['RMSE'])]

    if first_iter:
        combined_chosen_df = chosen_df
        points = len(df.index)
        first_iter = False
    else:
        combined_chosen_df = pd.concat([combined_chosen_df, chosen_df])
        points += len(df.index)

param_id = combined_chosen_df['param_id']['param_id'].values
Vhalf_range = combined_chosen_df['param_values']['Vhalf'].values
Kmax_range = combined_chosen_df['param_values']['Kmax'].values
Ku_range = combined_chosen_df['param_values']['Ku'].values

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

starting_param_df = pd.DataFrame([1] * 5, index=param_names).T
ComparisonController = modelling.ModelComparison(starting_param_df)


def param_evaluation(param_values):

    param_values = param_values.drop(columns=['param_id'])

    ComparisonController.drug_param_values = param_values

    Hill_curve_coefs, drug_conc_Hill, peaks_norm = \
        ComparisonController.compute_Hill(drug_model)
    # parameters of Hill's curve are based on normalised drug concentration

    drug_conc_AP = 10**np.linspace(np.log10(drug_conc_Hill[1]),
                                   np.log10(max(drug_conc_Hill)),
                                   APD_points)

    # APD_trapping, APD_conductance, drug_conc_AP = \
    #     ComparisonController.APD_sim(
    #         AP_model, Hill_curve_coefs, drug_conc=drug_conc_AP,
    #         EAD=True, abs_tol=abs_tol, rel_tol=rel_tol)

    for i in range(17, len(drug_conc_AP)):
        conc_ind = 17
#         log = AP_model.custom_simulation(param_values, drug_conc_AP[i], 1000,
#                                          timestep=0.1, save_signal=2,
#                                          abs_tol=1e-8, rel_tol=1e-4)
#         print('large tolerance value done')

        log = AP_model.custom_simulation(param_values, drug_conc_AP[conc_ind], 5.5,
                                         timestep=0.1, save_signal=2,
                                         abs_tol=1e-9, rel_tol=1e-10)
        print('small tolerance value done')
        print(log.time()[-10:])
        plt.figure()
        plt.plot(log.time(), log['0.ikr.C1'], 'o')
        plt.plot(log.time(), log['1.ikr.C1'], 'o')
        plt.savefig('../../figures/testing/tol_check.pdf')

    return 0


missing_id = [41]
# first_iter = True
# for i in range(len(Vhalf_range)):
for i in missing_id:
    paramid = param_id[i]
    # print(paramid)
    Vhalf = Vhalf_range[i]
    Kmax = Kmax_range[i]
    Ku = Ku_range[i]
    param_values = pd.DataFrame([paramid, Vhalf, Kmax, Ku, 1, 1],
                                index=['param_id'] + param_names)
    param_values = param_values.T
    result_df = param_evaluation(param_values)
    # , abs_tol=1e-8, rel_tol=1e-4)
    # print('high tolerance value done')
    # result_df = param_evaluation(param_values, abs_tol=1e-9, rel_tol=1e-10)
    # print('low tolerance value done')
    # result_df = param_evaluation(param_values, abs_tol=1e-7, rel_tol=1e-8)
    # print('mid tolerance value done')
    # if first_iter:
    #     combined_result_df = result_df.T
    #     first_iter = False
    # else:
    # combined_result_df = pd.concat([combined_result_df, result_df.T])

    # combined_result_df.to_csv(saved_data_dir + filename)
