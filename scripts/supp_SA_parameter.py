# Check that when changing Hill coefficient of each synthetic drug, the range
# of the RMSD between APD90s of the ORd-SD model and the ORd-CS model is small

import myokit
import numpy as np
import os
import pandas as pd
import pints

import modelling

# Define directory to save simulation data
data_dir = '../../simulation_data/supp_mat/APD90diff_N/'
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)

# Load IKr model and set up protocol
model = '../../math_model/ohara-cipa-v1-2017-IKr.mmt'
model, _, x = myokit.load(model)
drug_model = modelling.BindingKinetics(model)

protocol_params = modelling.ProtocolParameters()
protocol = protocol_params.protocol_parameters['Milnes']['function']
drug_model.protocol = protocol

# Load AP model
APmodel = '../../model/ohara-cipa-v1-2017.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')
pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)
base_conductance = APmodel.get('ikr.gKr').value()

# Define parameters used in simulations
offset = 50
save_signal = 2
repeats = 1000
APD_points = 20

# Getn list of synthetic drugs and the name of parameters
param_lib = modelling.BindingParameters()
drug_list = param_lib.drug_compounds
SA_model = modelling.SensitivityAnalysis()
param_names = SA_model.param_names
parameter_interest = 'N'


def param_evaluation(param, param_values):

    # Define parameter values of virtual drug
    orig_half_effect_conc = param_values['EC50'][0]
    param_values[parameter_interest][0] = param
    param_values['EC50'][0] = 1
    ComparisonController.drug_param_values = param_values

    # Calculate the normalising constant
    Hill_n = param_values['N'][0]
    norm_constant = np.power(orig_half_effect_conc, 1 / Hill_n)

    # Compute Hill curve of the synthetic drug with the SD model
    Hill_curve_coefs, drug_conc_Hill, peaks_norm = \
        ComparisonController.compute_Hill(drug_model,
                                          norm_constant=norm_constant,
                                          parallel=False)
    # The parameters of Hill's curve are based on the normalised
    # drug concentration
    # Hill's coefficient remains the same but IC50 -> IC50/EC50

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
            # Simulate APs and APD90s of the ORd-SD model and the ORd-CS model
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

    # Create dataframe to save results
    conc_Hill_ind = ['conc_' + str(i) for i, _ in
                     enumerate(drug_conc_Hill)]
    conc_AP_ind = ['conc_' + str(i) for i, _ in enumerate(drug_conc_AP)]
    index_dict = {'drug_conc_Hill': conc_Hill_ind,
                  'peak_current': conc_Hill_ind,
                  'Hill_curve': ['Hill_coef', 'IC50'],
                  'param_values': param_names, 'drug_conc_AP': conc_AP_ind,
                  'APD_trapping': conc_AP_ind,
                  'APD_conductance': conc_AP_ind, 'RMSE': ['RMSE'],
                  'ME': ['ME']}
    all_index = [(i, j) for i in index_dict.keys() for j in index_dict[i]]
    index = pd.MultiIndex.from_tuples(all_index)

    param_values['EC50'][0] = orig_half_effect_conc
    big_df = pd.DataFrame(
        drug_conc_Hill + list(peaks_norm) + list(Hill_curve_coefs) +
        list(param_values.values[0]) + list(drug_conc_AP) + APD_trapping +
        APD_conductance + [RMSError] + [MAError], index=index)

    return big_df


for drug in drug_list:
    # Get parameter values of each synthetic drug
    Vhalf = param_lib.binding_parameters[drug]['Vhalf']
    Kmax = param_lib.binding_parameters[drug]['Kmax']
    Ku = param_lib.binding_parameters[drug]['Ku']
    Hill_n = param_lib.binding_parameters[drug]['N']
    half_effect_conc = param_lib.binding_parameters[drug]['EC50']

    all_params = [Vhalf, Kmax, Ku, Hill_n, half_effect_conc]

    # Define parameter values to the system
    orig_param_values = pd.DataFrame(all_params, index=param_names)
    orig_param_values = orig_param_values.T
    ComparisonController = modelling.ModelComparison(orig_param_values)

    # Define the range of parameter values to explore
    param_values = orig_param_values
    param_range = SA_model.param_explore(parameter_interest, res_points=5)
    param_fullrange = SA_model.param_explore_gaps(param_range, 3, 'N')

    # Check for completed simulations to prevent repetition
    filename = 'SA_' + drug + '_' + parameter_interest + '.csv'
    if os.path.exists(data_dir + filename):
        saved_results_df = pd.read_csv(data_dir + filename,
                                       header=[0, 1], index_col=[0],
                                       skipinitialspace=True)
        ran_values = saved_results_df['param_values'][
            parameter_interest].values
    else:
        ran_values = []

    param_fullrange = [i for i in param_fullrange if i not in ran_values]

    # Evaluate the RMSD and MD between APD90s of a synthetic drug with
    # changing Hill coefficient from the ORd-SD model and the ORd-CS model
    n_workers = 8
    evaluator = pints.ParallelEvaluator(param_evaluation,
                                        n_workers=n_workers,
                                        args=[param_values])
    for i in range(int(np.ceil(len(param_fullrange) / n_workers))):
        print('Running samples ', n_workers * i, 'to',
              n_workers * (i + 1) - 1)
        big_df = evaluator.evaluate(
            param_fullrange[i * n_workers: (i + 1) * n_workers])

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
