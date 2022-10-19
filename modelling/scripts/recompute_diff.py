import os
import pandas as pd

import modelling

# saved_data_dir = '../../simulation_data/sensitivity_analysis/N/'
saved_data_dir = '../../simulation_results/'

param_names = modelling.SensitivityAnalysis().param_names

starting_param_df = pd.DataFrame([1] * 5, index=param_names).T
ComparisonController = modelling.ModelComparison(starting_param_df)

# Get all data
file_prefix = 'SA_allparam_gaps'
result_files = [f for f in os.listdir(saved_data_dir) if
                f.startswith(file_prefix)]
for file in result_files:
    print(file)
    # Load data
    saved_results_df = pd.read_csv(saved_data_dir + file,
                                   header=[0, 1], index_col=[0],
                                   skipinitialspace=True)
    saved_results_df = saved_results_df.reset_index(drop=True)

    for r in range(len(saved_results_df.index)):

        # Extract APD value
        APD_trapping = saved_results_df.iloc[[r]]['APD_trapping'].values[0]
        APD_conductance = saved_results_df.iloc[[r]][
            'APD_conductance'].values[0]

        # Compute new RMSE and MAE
        RMSError = ComparisonController.compute_RMSE(APD_trapping,
                                                     APD_conductance)
        MAError = ComparisonController.compute_MAE(APD_trapping,
                                                   APD_conductance)

        # Update RMSE and MAE
        saved_results_df.loc[r, ('RMSE', 'RMSE')] = RMSError
        saved_results_df.loc[r, ('MAE', 'MAE')] = MAError

    # Save dataframe
    saved_results_df.to_csv(saved_data_dir + file)
