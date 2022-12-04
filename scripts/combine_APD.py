# Combine all simulated data of the parameter space with essential information
# for easy loading when plotting figures
import os
import pandas as pd

# Read saved data from parameter space exploration
data_dir = '../../simulation_data/SA_space/'
file_prefix = 'SA_allparam_'
result_files = [data_dir + f for f in os.listdir(data_dir) if
                f.startswith(file_prefix)]

# Combine all files with selected columns
first_iter = True
for file in result_files:
    df = pd.read_csv(file, header=[0, 1], index_col=[0],
                     skipinitialspace=True)
    small_df = df[['param_id', 'param_values', 'APD_trapping',
                   'APD_conductance', 'RMSE', 'ME']]
    small_df = small_df.drop([('param_values', "EC50"), ('param_values', "N")],
                             axis=1)

    if first_iter:
        combined_df = small_df
        first_iter = False
    else:
        combined_df = pd.concat([combined_df, small_df])

# Save combined data
combined_df.to_csv(data_dir + 'SA_APD.csv')

# Repeat for a different file name prefix
file_prefix = 'SA_allparam_gaps_'
result_files = [data_dir + f for f in os.listdir(data_dir) if
                f.startswith(file_prefix)]

first_iter = True
for file in result_files:
    df = pd.read_csv(file,
                     header=[0, 1], index_col=[0],
                     skipinitialspace=True)
    small_df = df[['param_id', 'param_values', 'APD_trapping',
                   'APD_conductance', 'RMSE', 'ME']]
    small_df = small_df.drop([('param_values', "EC50"), ('param_values', "N")],
                             axis=1)

    if first_iter:
        combined_df = small_df
        first_iter = False
    else:
        combined_df = pd.concat([combined_df, small_df])

combined_df.to_csv(data_dir + 'SA_APD_gaps.csv')
