#
# Compute I_net and qNet for each example drug with the AP-SD model and the
# AP-CS model.
#

import myokit
import numpy as np
import os
import pandas as pd

import modelling

# I_net and qNet currents
current_list = ['inal.INaL', 'ical.ICaL', 'ikr.IKr', 'iks.IKs', 'ik1.IK1',
                'ito.Ito']

# Set up AP model
APmodel = '../math_model/ohara-cipa-v1-2017-opt.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')

# Define current protocol
pulse_time = 2000
AP_model.protocol = myokit.pacing.blocktrain(pulse_time, 0.5)
base_conductance = APmodel.get('ikr.gKr').value()
prepace = 1000

# Define constants for simulations
save_signal = 1
abs_tol = 1e-7
rel_tol = 1e-8

# Simulate AP-SD model at control conditions
control_log = AP_model.conductance_simulation(
    base_conductance, prepace, timestep=0.01, abs_tol=abs_tol,
    rel_tol=rel_tol)

# Define drugs
param_lib = modelling.BindingParameters()
drug_list = param_lib.drug_compounds[:-1]

for drug in drug_list:
    # Define the range of drug concentration for a given drug
    drug_conc_lib = modelling.DrugConcentrations()
    drug_conc = drug_conc_lib.drug_concentrations[drug]['fine']

    # Get Hill curve coefficient
    data_dir = '../simulation_data/model_comparison/' + drug + '/'
    result_filename = 'Hill_curves.csv'
    Hill_curves_df = pd.read_csv(data_dir + result_filename, index_col=0)
    Hill_coefs = Hill_curves_df[Hill_curves_df['protocol'] == 'Milnes']
    Hill_coefs = np.array(Hill_coefs.values.tolist()[0][:-1])
    Hill_model = modelling.HillsModel()

    # Define directories to save simulated data
    data_dir = data_dir + '/Milnes/qNet/'
    if not os.path.isdir(data_dir):
        os.makedirs(data_dir)

    # Simulate AP after addition of example drugs for the AP-SD model and
    # AP-CS model
    # Compute I_net and qNet for each condition
    qNet_SD_arr = []
    qNet_CS_arr = []
    for i in range(len(drug_conc)):
        print('simulating concentration: ' + str(drug_conc[i]))
        log = AP_model.drug_simulation(
            drug, drug_conc[i], prepace + save_signal,
            log_var=['engine.time', 'membrane.V'] + current_list,
            timestep=0.01, abs_tol=abs_tol, rel_tol=rel_tol,
            set_state=control_log)

        inet = 0
        for c in current_list:
            inet += log[c]  # pA/pF

        qNet = np.trapz(inet, x=log.time()) * 1e-3  # pA/pF*s
        qNet_SD_arr.append(qNet)

        reduction_scale = Hill_model.simulate(Hill_coefs[:2], drug_conc[i])
        d2 = AP_model.conductance_simulation(
            base_conductance * reduction_scale, prepace + save_signal,
            save_signal=save_signal, timestep=0.01,
            log_var=['engine.time', 'membrane.V'] + current_list,
            abs_tol=abs_tol, rel_tol=rel_tol, set_state=control_log)

        inet = 0
        for c in current_list:
            inet += d2[c]  # pA/pF

        qNet = np.trapz(inet, x=d2.time()) * 1e-3  # pA/pF*s
        qNet_CS_arr.append(qNet)

        print('done concentration: ' + str(drug_conc[i]))

    # Save qNet data
    qNet_data = {'drug_conc': drug_conc, 'SD': qNet_SD_arr, 'CS': qNet_CS_arr}
    qNet_df = pd.DataFrame(data=qNet_data)
    qNet_df.to_csv(data_dir + 'qNets.csv')
