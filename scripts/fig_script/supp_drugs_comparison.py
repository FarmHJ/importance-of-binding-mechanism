#
# Figure S1-S10
# Compares the (A) APD90 and (B) qNet between the AP-SD model and the
# AP-CS model for given drugs.
# (C) Plots the Hill curves of given drugs when stimulated with Milnes',
# Pneg80, P0 and P40 protocols.
#

import numpy as np
import os
import pandas as pd
import sys

import modelling

# Define drug and protocol
drug = sys.argv[1]
protocol_name = 'Milnes'

# Define directory to save figure and load results
fig_dir = '../../figures/supp_mat/model_comparison/'
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

data_dir = '../../simulation_data/model_comparison/' + drug + '/' + \
    protocol_name + '/'

# Set up structure of figure
fig = modelling.figures.FigureStructure(figsize=(10, 2.5),
                                        gridspec=(1, 3),
                                        wspace=0.35,
                                        plot_in_subgrid=True)
plot = modelling.figures.FigurePlot()

subgridspecs = [(1, 1)] * 3
subgs = []
for i in range(3):
    subgs.append(fig.gs[i].subgridspec(*subgridspecs[i]))
axs = [[[fig.fig.add_subplot(subgs[k][i, j]) for j in range(
    subgridspecs[k][1])] for i in range(subgridspecs[k][0])] for
    k in range(len(subgs))]

# Left panel
# Plot the APD90 of the AP-SD model and the AP-CS model for a given drug
panel5 = axs[0]

# Load the APD90 computed by the two models
APD_trapping = pd.read_csv(data_dir + 'SD_APD_fine.csv')
APD_conductance = pd.read_csv(data_dir + 'CS_APD_fine.csv')

# Organise data for plotting and identify EAD-like behaviour
drug_conc = APD_trapping['drug concentration'].values.tolist()
APD_trapping = [max(APD_trapping.loc[i].values.tolist()[1:-1]) for i in
                range(APD_trapping.shape[0])]
APD_conductance = [max(APD_conductance.loc[i].values.tolist()[1:-1])
                   for i in range(APD_conductance.shape[0])]
EAD_marker = [1050 if (i >= 1000 or j >= 1000) else None for (i, j)
              in zip(APD_trapping, APD_conductance)]

# Plot APD90 of both models
panel5[0][0].plot(drug_conc, APD_trapping, 'o', color='orange',
                  label='AP-SD model')
panel5[0][0].plot(drug_conc, APD_conductance, '^', color='blue',
                  label='AP-CS model', alpha=0.8)
panel5[0][0].scatter(drug_conc, EAD_marker, marker=(5, 2), color='k',
                     label='EAD-like AP')

# Adjust figure details
panel5[0][0].set_xscale("log", nonpositive='clip')
panel5[0][0].set_xlabel('Drug concentration (nM)')
panel5[0][0].set_ylabel(r'APD$_{90}$ (ms)')
panel5[0][0].legend(handlelength=1)

l_lim, r_lim = panel5[0][0].get_xlim()

# Middle panel
# Plot the qNet of the AP-SD model and the AP-CS model for a given drug
panel1 = axs[1]

# Load qNets data
qNets = pd.read_csv(data_dir + 'qNets.csv')

param_lib = modelling.BindingParameters()
Cmax = param_lib.binding_parameters[drug]['Cmax']

drug_conc = np.array(qNets['drug_conc'])
SD_qNet = np.array(qNets['SD'])
CS_qNet = np.array(qNets['CS'])

SD_qNet = [None if APD_trapping[i] >= 1000 else SD_qNet[i] for i in
           range(len(SD_qNet))]
CS_qNet = [None if APD_conductance[i] >= 1000 else CS_qNet[i] for i in
           range(len(CS_qNet))]

# Plot qNet of both models
panel1[0][0].plot(drug_conc, SD_qNet, 'o', color='orange',
                  label='AP-SD model')
panel1[0][0].plot(drug_conc, CS_qNet, '^', color='blue',
                  label='AP-CS model', alpha=0.8)
panel1[0][0].set_xlabel("Drug concentration (nM)")
panel1[0][0].set_xscale("log", nonpositive='clip')
panel1[0][0].set_xlim(left=l_lim, right=r_lim)
panel1[0][0].set_ylabel('qNet (C/F)')
panel1[0][0].legend(handlelength=1)

# Right panel
# Hill curves of the drug under different protocol stimulation of the SD model
panel6 = axs[2]

# Set up Hill equation and protocols
Hill_model = modelling.HillsModel()
protocol_list = modelling.ProtocolParameters().protocols
line_pattern = ['solid', 'dashed', 'dotted', 'dashdot']
color = ['orange', 'blue', 'red', 'green']

# Load Hill curve coefficients
Hill_coef_df = pd.read_csv(data_dir + '../Hill_curves.csv')

# Define range of drug concentration
drug_conc_lib = modelling.DrugConcentrations()
conc_grid = drug_conc_lib.drug_concentrations[drug]['fine']

# Simulate and plot the Hill curve
for p, prot in enumerate(protocol_list):
    Hill_eq = Hill_coef_df.loc[
        Hill_coef_df['protocol'] == prot]
    Hill_eq = Hill_eq.values.tolist()[0][1:-1]

    panel6[0][0].plot(conc_grid, Hill_model.simulate(
        Hill_eq, conc_grid), linestyle=line_pattern[p], color=color[p],
        label=prot)

# Adjust figure details
panel6[0][0].set_xscale("log", nonpositive='clip')
panel6[0][0].set_xlabel('Drug concentration (nM)')
panel6[0][0].set_ylabel('Normalised\npeak current')
panel6[0][0].legend(handlelength=3)

# Add panel labels
fig.fig.text(0.1, 0.9, '(A)', fontsize=11)
fig.fig.text(0.38, 0.9, '(B)', fontsize=11)
fig.fig.text(0.66, 0.9, '(C)', fontsize=11)

# Save figure
fig.savefig(fig_dir + drug + '_APD_Hill.pdf')
