import matplotlib.pyplot as plt
import myokit
import numpy as np
import pandas as pd

import modelling

saved_data_filepath = '../../simulation_data/protocol_test/'

# Model directory
model_dir = '../../model/'
current_model_filename = 'ohara-cipa-v1-2017-IKr.mmt'
current_model_filepath = model_dir + current_model_filename
AP_model_filename = 'ohara-cipa-v1-2017.mmt'
AP_model_filepath = model_dir + AP_model_filename

# Load hERG model
model, _, x = myokit.load(current_model_filepath)

protocol_params = modelling.ProtocolParameters()
protocol = protocol_params.protocol_parameters['Milnes']['function']

drug_model = modelling.BindingKinetics(model)
drug_model.protocol = protocol

# Load AP model
APmodel, _, x = myokit.load(AP_model_filepath)
AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')
pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)

# Parameters used in simulations
offset = 50
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

    prot_id = param_values['prot_id'][0]
    param_values = param_values.drop(columns=['prot_id'])

    ComparisonController.drug_param_values = param_values

    Hill_curve_coefs, drug_conc_Hill, peaks_norm = \
        ComparisonController.compute_Hill(drug_model,
                                          parallel=True)
    # parameters of Hill's curve are based on normalised drug concentration

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
        # Simulate action potentials
        # try:
        APD_trapping, APD_conductance, drug_conc_AP = \
            ComparisonController.APD_sim(
                AP_model, Hill_curve_coefs, drug_conc=drug_conc_AP,
                EAD=True, abs_tol=1e-8, rel_tol=1e-10)

        RMSError = ComparisonController.compute_RMSE(APD_trapping,
                                                     APD_conductance)
        MAError = ComparisonController.compute_MAE(APD_trapping,
                                                   APD_conductance)
        # except myokit.SimulationError:
        #     APD_trapping = [float("Nan")] * APD_points
        #     APD_conductance = [float("Nan")] * APD_points
        #     RMSError = float("Nan")
        #     MAError = float("Nan")

    # Create dataframe to save results
    conc_Hill_ind = ['conc_' + str(i) for i, _ in
                     enumerate(drug_conc_Hill)]
    conc_AP_ind = ['conc_' + str(i) for i, _ in enumerate(drug_conc_AP)]
    index_dict = {'prot_id': ['prot_id'],
                  'drug_conc_Hill': conc_Hill_ind,
                  'peak_current': conc_Hill_ind,
                  'Hill_curve': ['Hill_coef', 'IC50'],
                  'param_values': param_names, 'drug_conc_AP': conc_AP_ind,
                  'APD_trapping': conc_AP_ind,
                  'APD_conductance': conc_AP_ind, 'RMSE': ['RMSE'],
                  'MAE': ['MAE']}
    all_index = [(i, j) for i in index_dict.keys() for j in index_dict[i]]
    index = pd.MultiIndex.from_tuples(all_index)

    big_df = pd.DataFrame(
        [prot_id] + drug_conc_Hill + list(peaks_norm) +
        list(Hill_curve_coefs) + list(param_values.values[0]) +
        list(drug_conc_AP) + APD_trapping + APD_conductance +
        [RMSError] + [MAError], index=index)

    return big_df


drug = 'verapamil'
Vhalf = param_lib.binding_parameters[drug]['Vhalf']
Kmax = param_lib.binding_parameters[drug]['Kmax']
Ku = param_lib.binding_parameters[drug]['Ku']
Hill_n = param_lib.binding_parameters[drug]['N']
half_effect_conc = param_lib.binding_parameters[drug]['EC50']

all_params = [Vhalf, Kmax, Ku, Hill_n, half_effect_conc]

param_values = pd.DataFrame(all_params, index=param_names)
param_values = param_values.T

t_max = 25e3

protocol1 = myokit.Protocol()
protocol1.schedule(-80, 0, 800, period=t_max)
protocol1.schedule(-90, 800, 100, period=t_max)
protocol1.schedule(-80, 900, 100, period=t_max)
protocol1.schedule(-80, 11000, 14000, period=t_max)
print(protocol1.characteristic_time())

protocol2 = myokit.Protocol()
protocol2.schedule(-80, 0, 800)
protocol2.schedule(-90, 800, 100)
protocol2.schedule(-80, 900, 100)
protocol2.schedule(-80, 11000, 14000)
print(protocol2.characteristic_time())

protocol3 = myokit.Protocol()
protocol3.schedule(-80, 0, 800, period=t_max)
protocol3.schedule(-90, 800, 100, period=t_max)
protocol3.schedule(-80, 900, 100, period=t_max)
protocol3.schedule(-80, 11000, 14000 - 1, period=t_max)
print(protocol3.characteristic_time())

protocols = [protocol1, protocol2, protocol3]
for i in range(3):
    drug_model.protocol = protocols[i]

    param_values['prot_id'] = i
    df = param_evaluation(param_values)
    if i == 0:
        combined_df = df.T
    else:
        combined_df = pd.concat([combined_df, df.T])

combined_df.to_csv(saved_data_filepath + 'APD_verapamil.csv')

for i in range(len(combined_df.index)):

    drug_conc_AP = combined_df.iloc[[i]]['drug_conc_AP'].values[0]
    APD_trapping = combined_df.iloc[[i]]['APD_trapping'].values[0]
    APD_conductance = combined_df.iloc[[i]]['APD_conductance'].values[0]

    plt.figure()
    plt.plot(drug_conc_AP, APD_trapping, 'o')
    plt.plot(drug_conc_AP, APD_conductance, '^')
    plt.xscale('log')
    plt.savefig(saved_data_filepath + 'verapamil_APD_prot' + str(i) + '.pdf')
    plt.close()

plt.figure()
for i in range(len(combined_df.index)):

    drug_conc_Hill = combined_df.iloc[[i]]['drug_conc_Hill'].values[0]
    peak_current = combined_df.iloc[[i]]['peak_current'].values[0]

    plt.plot(drug_conc_Hill, peak_current)
plt.xscale('log')
plt.savefig(saved_data_filepath + 'verapamil_Hill_prot' + str(i) + '.pdf')
plt.close()