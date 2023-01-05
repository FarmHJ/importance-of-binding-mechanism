import myokit
import numpy as np
import os
import pandas as pd
import sys

import modelling

# I_net and qNet currents
current_list = ['inal.INaL', 'ical.ICaL', 'ikr.IKr', 'iks.IKs', 'ik1.IK1',
                'ito.Ito']

# Define drug
drug = sys.argv[1]

# Define the range of drug concentration for a given drug
drug_conc_lib = modelling.DrugConcentrations()
drug_conc = drug_conc_lib.drug_concentrations[drug]['fine']
prepace = 1000

# Get Hill curve coefficient
data_dir = '../simulation_data/model_comparison/' + \
    drug + '/Milnes/'
result_filename = 'Hill_curve.txt'
Hill_coefs = np.loadtxt(data_dir + result_filename, unpack=True)
Hill_coefs = np.array(Hill_coefs)
Hill_model = modelling.HillsModel()

# Define directories to save simulated data
data_dir = data_dir + 'qNet/'
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)

# Set AP model
APmodel = '../math_model/ohara-cipa-v1-2017-opt.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')

# Define current protocol
pulse_time = 2000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)
base_conductance = APmodel.get('ikr.gKr').value()

# Simulate AP of the ORd-SD model and the ORd-CS model
# Compute APD90
save_signal = 1
abs_tol = 1e-7
rel_tol = 1e-8
qNet_SD_arr = []
qNet_CS_arr = []
for i in range(len(drug_conc)):
    print('simulating concentration: ' + str(drug_conc[i]))
    log = AP_model.drug_simulation(
        drug, drug_conc[i], prepace + save_signal, save_signal=save_signal,
        log_var=['engine.time', 'membrane.V'] + current_list, timestep=0.01,
        abs_tol=abs_tol, rel_tol=rel_tol)
    log.save_csv(data_dir + 'SD_AP_inetcurrents_' + str(i) + '.csv')

    inet = 0
    for c in current_list:
        inet += log[c]  # pA/pF
    inet_data = {'time': log.time(), 'inet': inet}
    inet_df = pd.DataFrame(inet_data, columns=['time', 'inet'])
    inet_df.to_csv(data_dir + 'SD_inet_' + str(i) + '.csv', index=False)

    qNet = np.trapz(inet, x=log.time()) * 1e-3  # pA/pF*s
    qNet_SD_arr.append(qNet)

    reduction_scale = Hill_model.simulate(Hill_coefs[:2], drug_conc[i])
    d2 = AP_model.conductance_simulation(
        base_conductance * reduction_scale, prepace + save_signal,
        save_signal=save_signal, timestep=0.01,
        log_var=['engine.time', 'membrane.V'] + current_list,
        abs_tol=abs_tol, rel_tol=rel_tol)
    d2.save_csv(data_dir + 'CS_AP_inetcurrents_' + str(i) + '.csv')

    inet = 0
    for c in current_list:
        inet += d2[c]  # pA/pF
    inet_data = {'time': d2.time(), 'inet': inet}
    inet_df = pd.DataFrame(inet_data, columns=['time', 'inet'])
    inet_df.to_csv(data_dir + 'CS_inet_' + str(i) + '.csv', index=False)

    qNet = np.trapz(inet, x=d2.time()) * 1e-3  # pA/pF*s
    qNet_CS_arr.append(qNet)

    print('done concentration: ' + str(drug_conc[i]))

qNet_data = {'drug_conc': drug_conc, 'SD': qNet_SD_arr, 'CS': qNet_CS_arr}
qNet_df = pd.DataFrame(data=qNet_data)
qNet_df.to_csv(data_dir + 'qNets.csv')
