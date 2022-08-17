# Compares the hERG current and action potential of state dependen drug block
# and conductance scaling
import pandas as pd

import modelling

# Set up directory for figures
testing_fig_dir = '../../figures/testing/'
final_fig_dir = '../../figures/binding_kinetics_comparison/'

saved_fig_dir = final_fig_dir

# Load data
drug = 'dofetilide'
protocol_name = 'Milnes'
saved_data_dir = '../../simulation_data/binding_kinetics_comparison/' + \
    drug + '/' + protocol_name + '/'

APD_trapping = pd.read_csv(saved_data_dir + 'CiPA_APD_large.csv')
APD_conductance = pd.read_csv(saved_data_dir + 'conductance_APD_large.csv')

drug_conc = APD_trapping['drug concentration'].values.tolist()
APD_trapping = [max(APD_trapping.loc[i].values.tolist()[1:-1]) for i in
                range(APD_trapping.shape[0])]
APD_conductance = [max(APD_conductance.loc[i].values.tolist()[1:-1])
                   for i in range(APD_conductance.shape[0])]
EAD_marker = [1050 if (i >= 1000 or j >= 1000) else None for (i, j)
              in zip(APD_trapping[1:], APD_conductance[1:])]

drug_conc = drug_conc[:-5]
APD_trapping = APD_trapping[:-5]
APD_conductance = APD_conductance[:-5]
EAD_marker = EAD_marker[:-5]

# Set up figure
fig = modelling.figures.FigureStructure(figsize=(8, 4), gridspec=(1, 2),
                                        wspace=0.3)
plot = modelling.figures.FigurePlot()

# Left panel
fig.axs[0][0].plot(drug_conc[1:], APD_trapping[1:],
                   'o', color='orange', label='state dependent drug block')
fig.axs[0][0].plot(drug_conc[1:], APD_conductance[1:],
                   '^', color='blue', label='conductance scaling', alpha=0.8)
fig.axs[0][0].plot(drug_conc[1:], EAD_marker, 'o', color='k',
                   marker=(5, 2), label='EAD-like AP')
fig.axs[0][0].set_xscale("log", nonpositive='clip')
fig.axs[0][0].set_xlabel('Drug concentration (nM)')
fig.axs[0][0].set_ylabel(r'APD$_{90}$ (ms)')
fig.axs[0][0].set_title('dofetilide-like drug')
fig.axs[0][0].legend(handlelength=1)

drug = 'verapamil'
saved_data_dir = '../../simulation_data/binding_kinetics_comparison/' + \
    drug + '/' + protocol_name + '/'

APD_trapping = pd.read_csv(saved_data_dir + 'CiPA_APD_large.csv')
APD_conductance = pd.read_csv(saved_data_dir + 'conductance_APD_large.csv')

drug_conc = APD_trapping['drug concentration'].values.tolist()
APD_trapping = [max(APD_trapping.loc[i].values.tolist()[1:-1]) for i in
                range(APD_trapping.shape[0])]
APD_conductance = [max(APD_conductance.loc[i].values.tolist()[1:-1])
                   for i in range(APD_conductance.shape[0])]
EAD_marker = [1050 if (i >= 1000 or j >= 1000) else None for (i, j)
              in zip(APD_trapping[1:], APD_conductance[1:])]

# Right panel
fig.axs[0][1].plot(drug_conc[1:], APD_trapping[1:], 'o', color='orange',
                   label='state dependent drug block')
fig.axs[0][1].plot(drug_conc[1:], APD_conductance[1:], '^', color='blue',
                   label='conductance scaling', alpha=0.8)
fig.axs[0][1].plot(drug_conc[1:], EAD_marker, 'o', color='k',
                   marker=(5, 2), label='EAD-like AP')
fig.axs[0][1].set_xscale("log", nonpositive='clip')
fig.axs[0][1].set_xlabel('Drug concentration (nM)')
fig.axs[0][1].set_ylabel(r'APD$_{90}$ (ms)')
fig.axs[0][1].set_title('verapamil-like drug')
fig.axs[0][1].legend(handlelength=1)

fig.fig.set_size_inches(8, 3)
fig.fig.text(0.075, 0.925, '(a)', fontsize=11)
fig.fig.text(0.5, 0.925, '(b)', fontsize=11)

fig.savefig(saved_fig_dir + 'APD_compare.svg', format='svg')
