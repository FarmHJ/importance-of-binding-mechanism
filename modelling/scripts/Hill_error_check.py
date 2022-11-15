# Happens only in global sensitivity analysis
# Not important for now

import matplotlib.pyplot as plt
import myokit
import numpy as np
import os
import pandas as pd

import modelling

saved_data_dir = '../../'
param_names = ['Vhalf', 'Kmax', 'Ku', 'N', 'EC50']
param_values = np.loadtxt(saved_data_dir + "param_value_samples.txt")

model = '../../model/ohara-cipa-v1-2017-IKr.mmt'
model, _, x = myokit.load(model)

protocol_params = modelling.ProtocolParameters()
protocol = protocol_params.protocol_parameters['Milnes']['function']

BKmodel = modelling.BindingKinetics(model)
BKmodel.protocol = protocol

AP_model_filepath = '../../model/ohara-cipa-v1-2017.mmt'
APmodel, _, x = myokit.load(AP_model_filepath)
AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')

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

    ComparisonController.drug_param_values = orig_param_values
    Hill_curve_coefs, drug_conc_Hill, _ = \
        ComparisonController.compute_Hill(BKmodel, parallel=False)

    if isinstance(Hill_curve_coefs, str):
        return float("nan"), float("nan"), np.inf, 0

    drug_conc_AP = 10**np.linspace(np.log10(drug_conc_Hill[1]),
                                   np.log10(max(drug_conc_Hill)),
                                   20)
    try:
        APD_trapping, APD_conductance, _ = ComparisonController.APD_sim(
            AP_model, Hill_curve_coefs, drug_conc=drug_conc_AP)
        RMSError = ComparisonController.compute_RMSE(APD_trapping,
                                                     APD_conductance)
        MAError = ComparisonController.compute_MAE(APD_trapping,
                                                   APD_conductance)

    except myokit.SimulationError:
        RMSError = float("Nan")
        MAError = float("nan")

    return RMSError, MAError

print(param_values)
