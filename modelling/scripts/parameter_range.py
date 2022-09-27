# To plot the range of parameter values
# Drug binding-related parameters

import itertools
import myokit
import numpy as np
import pandas as pd

import modelling

testing_fig_dir = '../../figures/testing/'
final_fig_dir = '../../figures/binding_kinetics_comparison/' + \
    'OHaraCiPA_model/sensitivity_analysis'

saved_fig_dir = testing_fig_dir

param_lib = modelling.BindingParameters()
drug_list = param_lib.drug_compounds

Vhalf = []
Kmax = []
Ku = []
N = []
EC50 = []

for drug in drug_list:
    Vhalf.append(param_lib.binding_parameters[drug]['Vhalf'])
    Kmax.append(param_lib.binding_parameters[drug]['Kmax'])
    Ku.append(param_lib.binding_parameters[drug]['Ku'])
    N.append(param_lib.binding_parameters[drug]['N'])
    EC50.append(param_lib.binding_parameters[drug]['EC50'])

all_params = [Vhalf, Kmax, Ku, N, EC50]
labels = ['Vhalf', 'Kmax', 'Ku', 'N', 'EC50']
param_ranges = [[-150, -50], [1e2, 1e6], [5e-3, 7e-2], [0.8, 1.2], [1e3, 1e8]]

# Plot of categorisation of parameters
# nrow = 2
# fig = modelling.figures.FigureStructure(
#     figsize=(10, 2.5 * nrow),
#     gridspec=(nrow, 3), hspace=0.3,
#     wspace=0.2,
#     height_ratios=[1] * nrow)

# for i, label in enumerate(labels):
#     fig.axs[int(i / 3)][i % 3].axhline(
#         param_ranges[i][0], xmin=0, xmax=len(all_params[i]),
#         c='red')
#     fig.axs[int(i / 3)][i % 3].axhline(
#         param_ranges[i][1], xmin=0, xmax=len(all_params[i]),
#         c='red')
#     fig.axs[int(i / 3)][i % 3].scatter(
#         np.arange(len(all_params[i])), all_params[i])
#     fig.axs[int(i / 3)][i % 3].set_title(label)

# fig.axs[0][1].set_yscale("log", nonpositive='clip')
# fig.axs[0][2].set_yscale("log", nonpositive='clip')
# fig.axs[1][1].set_yscale("log", nonpositive='clip')

# fig.savefig(saved_fig_dir + 'parameter_categorisation.pdf')

# Set AP model
APmodel = '../../model/ohara-cipa-v1-2017.mmt'

# Load model
APmodel, _, x = myokit.load(APmodel)

AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')
pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)
base_conductance = APmodel.get('ikr.gKr').value()

offset = 50
repeats = 800
save_signal = 2

# Simple SA - taking mean value of each category
category_means = []
for i, param in enumerate(all_params):
    low_cat = [k for k in param if k < param_ranges[i][0]]
    med_cat = [k for k in param if k >= param_ranges[i][0]
               and k <= param_ranges[i][1]]
    high_cat = [k for k in param if k > param_ranges[i][1]]
    category_means.append([np.mean(low_cat), np.mean(med_cat),
                           np.mean(high_cat)])

drug_conc = 10**np.linspace(-1, 1, 3)
# drug_conc = 10**np.linspace(-1, 7, 20)

for x in itertools.product(*category_means):
    param_values = pd.DataFrame(x, index=labels)
    param_values = param_values.T

    # Run simulations
    APD_trapping = []
    for i in range(len(drug_conc)):
        log = AP_model.CiPA_simulation(
            param_values, drug_conc[i], repeats, timestep=0.1,
            save_signal=save_signal,
            log_var=['engine.time', 'membrane.V'])

        # Compute APD90
        APD_trapping_pulse = []
        for pulse in range(save_signal):
            apd90 = AP_model.APD90(log['membrane.V', pulse], offset, 0.1)
            APD_trapping_pulse.append(apd90)
        APD_trapping.append(APD_trapping_pulse)
    APD_trapping = [max(i) for i in APD_trapping]
    print(APD_trapping)
