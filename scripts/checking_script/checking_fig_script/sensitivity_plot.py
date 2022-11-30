# To plot the range of parameter values
# Drug binding-related parameters

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import modelling


testing_fig_dir = '../../figures/testing/'
final_fig_dir = '../../figures/binding_kinetics_comparison/' + \
    'OHaraCiPA_model/sensitivity_analysis/'

saved_fig_dir = final_fig_dir

saved_data_dir = '../../simulation_data/sensitivity_analysis/'

SA_model = modelling.SensitivityAnalysis()
param_interest = 'N'

# for drug in ['dofetilide', 'terfenadine', 'cisapride', 'bepridil']:
# for drug in ['dofetilide']:

#     filename = 'SA_' + drug + '_' + param_interest + '.csv'
#     df = pd.read_csv(saved_data_dir + filename,
#                      header=[0, 1], index_col=[0],
#                      skipinitialspace=True) 
    # data included: drug_conc_Hill, peak_current, Hill_curve, param_values,
    # drug_conc_AP, APD_trapping, APD_conductance and MSE

    # ran_values = df['param_values'][param_interest].values
    # print(ran_values)
    # param_range = SA_model.param_explore_drug(drug, param_interest)
    # param_range = [i for i in param_range if i not in ran_values]
    # print(param_range)
    # param_lib = modelling.BindingParameters()
    # param_true = param_lib.binding_parameters[drug][param_interest]

    # Plot Hill curve
    # Can plot the change in APD90 in jupyter notebook
    # for row_ind in range(len(df.index)):
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

# 3D plot in parameter space
# Plot for known drugs
param_lib = modelling.BindingParameters()
drug_list = param_lib.drug_compounds

SA_model = modelling.SensitivityAnalysis()
param_names = SA_model.param_names

starting_param_df = pd.DataFrame([1] * 5, index=param_names).T
ComparisonController = modelling.ModelComparison(starting_param_df)

discrete_colors = ['red', 'blue', 'black']
APD_diff_label = ['similar', 'SD higher', 'CS higher']

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

for drug in drug_list:
    saved_data_dir = '../../simulation_data/binding_kinetics_comparison/' + \
        drug + '/Milnes/'
    APD_trapping = pd.read_csv(saved_data_dir + 'CiPA_APD_pulses1000.csv')
    APD_conductance = pd.read_csv(
        saved_data_dir + 'conductance_APD_pulses1000.csv')

    drug_conc = APD_trapping['drug concentration'].values.tolist()
    APD_trapping = [max(APD_trapping.loc[i].values.tolist()[1 + 998:-1])
                    for i in range(APD_trapping.shape[0])]
    APD_conductance = [max(APD_conductance.loc[i].values.tolist()[1 + 998:-1])
                       for i in range(APD_conductance.shape[0])]
    RMSError = ComparisonController.compute_RMSE(APD_trapping,
                                                 APD_conductance)
    MAError = ComparisonController.compute_MAE(APD_trapping,
                                               APD_conductance)

    if RMSError < 100:
        color_code = 0
    elif MAError > 0:
        color_code = 1
    else:
        color_code = 2

    Vhalf = param_lib.binding_parameters[drug]['Vhalf']
    Kmax = param_lib.binding_parameters[drug]['Kmax']
    Ku = param_lib.binding_parameters[drug]['Ku']

    ax.scatter(Vhalf, np.log(Kmax), np.log(Ku),
               c=discrete_colors[color_code], s=100,
               label=APD_diff_label[color_code], marker='^',
               zorder=-10)

# saved_data_dir = '../../simulation_data/sensitivity_analysis/'
# filename = 'SA_allparam_0_copy.csv'
# df = pd.read_csv(saved_data_dir + filename,
#                  header=[0, 1], index_col=[0],
#                  skipinitialspace=True)

# # Exploring the space
# Vhalf_range = df['param_values']['Vhalf'].values
# Kmax_range = df['param_values']['Kmax'].values
# Ku_range = df['param_values']['Ku'].values

# RMSError = df['RMSE']['RMSE'].values
# MAError = df['MAE']['MAE'].values

# param_id = df['param_id']['param_id'].values
# color_code_list = []
# for i in range(len(param_id)):
#     if RMSError[i] < 100:
#         color_code = 0
#     elif MAError[i] > 0:
#         color_code = 1  # SD higher
#     else:
#         color_code = 2  # CS higher
#     color_code_list.append(color_code)


# for j in range(len(param_id)):
#     ax.scatter(Vhalf_range[j], np.log(Kmax_range[j]), np.log(Ku_range[j]),
#                c=discrete_colors[color_code_list[j]], s=50,
#                label=APD_diff_label[color_code_list[j]], alpha=0.5,
#                zorder=-10)

handles, labels = ax.get_legend_handles_labels()
unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if
          l not in labels[:i]]
ax.legend(*zip(*unique), loc='upper left', bbox_to_anchor=(1.0, 1.0))
# ax.set_facecolor('silver')
ax.set_xlabel('Vhalf')
ax.set_ylabel('Kmax')
ax.set_zlabel('Ku')
ax.set_rasterization_zorder(0)
plt.show()
