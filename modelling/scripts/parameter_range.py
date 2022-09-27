# To plot the range of parameter values
# Drug binding-related parameters

import matplotlib.pyplot as plt
import numpy as np

import modelling

testing_fig_dir = '../../figures/testing/'
# final_fig_dir = '../../figures/binding_kinetics_comparison/' + drug + '/'

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

fig, axs = plt.subplots(2, 3)

# rectangular box plot
for i, label in enumerate(labels):
    axs[int(i / 3), i % 3].scatter(np.arange(len(all_params[i])), all_params[i])
    axs[int(i / 3), i % 3].set_title(label)
plt.savefig(saved_fig_dir + 'parameter_range.pdf')
