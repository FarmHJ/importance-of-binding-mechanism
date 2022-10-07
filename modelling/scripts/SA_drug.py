# To plot the range of parameter values
# Drug binding-related parameters

import matplotlib.pyplot as plt
import myokit
import numpy as np
import os
import pandas as pd

import modelling

saved_data_dir = '../../simulation_data/sensitivity_analysis/'

testing_fig_dir = '../../figures/testing/'
final_fig_dir = '../../figures/binding_kinetics_comparison/' + \
    'OHaraCiPA_model/sensitivity_analysis/'

check_plot = True
saved_fig_dir = final_fig_dir

param_lib = modelling.BindingParameters()
drug_list = param_lib.drug_compounds

drug = 'dofetilide'
Vhalf = param_lib.binding_parameters[drug]['Vhalf']
Kmax = param_lib.binding_parameters[drug]['Kmax']
Ku = param_lib.binding_parameters[drug]['Ku']
Hill_n = param_lib.binding_parameters[drug]['N']
half_effect_conc = param_lib.binding_parameters[drug]['EC50']

all_params = [Vhalf, Kmax, Ku, Hill_n, half_effect_conc]

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
base_conductance = APmodel.get('ikr.gKr').value()

offset = 50
repeats_AP = 800
save_signal = 2
repeats = 1000
APD_points = 20

# Perturbing the parameters
SA_model = modelling.SensitivityAnalysis()
param_ranges = []
param_names = SA_model.param_names

orig_param_values = pd.DataFrame(all_params, index=param_names)
orig_param_values = orig_param_values.T
ComparisonController = modelling.ModelComparison(orig_param_values)

for k, parameter_interest in enumerate(param_names):
    filename = 'SA_' + drug + '_' + parameter_interest + '.csv'
    if os.path.exists(saved_data_dir + filename):
        os.remove(saved_data_dir + filename)

    param_range = SA_model.param_explore_drug(drug, parameter_interest)
    param_values = orig_param_values

    for num, param in enumerate(param_range):

        param_values[parameter_interest][0] = param
        ComparisonController.drug_param_values = param_values

        Hill_curve_coefs, drug_conc_Hill, peaks_norm = \
            ComparisonController.compute_Hill(drug_model)

        drug_conc_range = (np.log10(drug_conc_Hill[1]),
                           np.log10(max(drug_conc_Hill)))
        try:
            APD_trapping, APD_conductance = ComparisonController.APD_sim(
                AP_model, Hill_curve_coefs, drug_conc=drug_conc_range)

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
        drug_conc_AP = 10**np.linspace(*drug_conc_range, APD_points)
        conc_AP_ind = ['conc_' + str(i) for i, _ in enumerate(drug_conc_AP)]
        index_dict = {'drug_conc_Hill': conc_Hill_ind,
                      'peak_current': conc_Hill_ind,
                      'Hill_curve': ['Hill_coef', 'IC50'],
                      'param_values': param_names, 'drug_conc_AP': conc_AP_ind,
                      'APD_trapping': conc_AP_ind,
                      'APD_conductance': conc_AP_ind, 'RMSE': ['RMSE'],
                      'MAE': ['MAE']}
        all_index = [(i, j) for i in index_dict.keys() for j in index_dict[i]]
        index = pd.MultiIndex.from_tuples(all_index)

        big_df = pd.DataFrame(
            drug_conc_Hill + list(peaks_norm) + list(Hill_curve_coefs) +
            list(param_values.values[0]) + list(drug_conc_AP) + APD_trapping +
            APD_conductance + [RMSError] + [MAError], index=index)

        if num != 0:
            previous_df = pd.read_csv(saved_data_dir + filename,
                                      header=[0, 1], index_col=[0],
                                      skipinitialspace=True)
            comb_df = pd.concat([previous_df, big_df.T])
        else:
            comb_df = big_df.T

        comb_df.to_csv(saved_data_dir + filename)

        os.system('cp ' + saved_data_dir + filename + ' ' +
                  saved_data_dir + filename[:-4] + '_copy.csv')
