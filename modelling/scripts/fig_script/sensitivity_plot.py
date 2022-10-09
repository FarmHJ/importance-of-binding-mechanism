# To plot the range of parameter values
# Drug binding-related parameters

import matplotlib.pyplot as plt
import pandas as pd

import modelling


testing_fig_dir = '../../figures/testing/'
final_fig_dir = '../../figures/binding_kinetics_comparison/' + \
    'OHaraCiPA_model/sensitivity_analysis/'

saved_fig_dir = final_fig_dir

saved_data_dir = '../../simulation_data/sensitivity_analysis/'

param_interest = 'N'

# for drug in ['dofetilide', 'terfenadine', 'cisapride', 'bepridil']:
for drug in ['dofetilide']:

    filename = 'SA_' + drug + '_' + param_interest + '.csv'
    df = pd.read_csv(saved_data_dir + filename,
                     header=[0, 1], index_col=[0],
                     skipinitialspace=True)
    # data included: drug_conc_Hill, peak_current, Hill_curve, param_values,
    # drug_conc_AP, APD_trapping, APD_conductance and MSE

    param_lib = modelling.BindingParameters()
    param_true = param_lib.binding_parameters[drug][param_interest]

    # Plot Hill curve
    # Can plot the change in APD90 in jupyter notebook
    for row_ind in range(len(df.index)):
        # fig = modelling.figures.FigureStructure(
        #     figsize=(10, 3),
        #     gridspec=(1, 2),
        #     wspace=0.2)

        # figname = drug + '_' + param_interest + '_' + str(row_ind) + '_APD.pdf'

        # APD_trapping = df.iloc[[row_ind]]['APD_trapping'].values[0]
        # APD_conductance = df.iloc[[row_ind]]['APD_conductance'].values[0]
        # fig.axs[0][0].plot(drug_conc_AP, APD_trapping, 'o-', color='red')
        # fig.axs[0][0].plot(drug_conc_AP, APD_conductance, '^-', color='red')

        # fig.axs[0][1].set_xscale('log')
        # fig.axs[0][0].set_xscale('log')
        # plt.savefig(saved_fig_dir + figname)
        # plt.close()

    # runs = len(df.index)
    # interest_param_values = df['param_values'][param_interest].values
    # RMSEs = df['RMSE']['RMSE'].values
    # MAEs = df['MAE']['MAE'].values

    # runs_check = len(df_check.index)
    # interest_param_values2 = df_check['param_values'][param_interest].values
    # RMSEs2 = df_check['RMSE']['RMSE'].values
    # MAEs2 = df_check['MAE']['MAE'].values

    # figname = drug + '_' + param_interest + '_RMSE.pdf'
    # plt.figure()
    # plt.plot(interest_param_values, RMSEs, 'o-')
    # plt.plot(interest_param_values, RMSEs2, '^-')
    # # plt.vlines(param_true, min(RMSEs), max(RMSEs), colors='red')
    # plt.xscale('log')
    # plt.savefig(saved_fig_dir + figname)

    # figname = drug + '_' + param_interest + '_MAE.pdf'
    # plt.figure()
    # plt.plot(interest_param_values, MAEs, 'o-')
    # plt.plot(interest_param_values, MAEs2, '^-')
    # # plt.vlines(param_true, min(MAEs), max(MAEs), colors='red')
    # plt.xscale('log')
    # plt.savefig(saved_fig_dir + figname)
