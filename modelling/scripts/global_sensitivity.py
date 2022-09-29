from SALib.sample import saltelli
from SALib.analyze import sobol
import numpy as np
import os

import modelling

# param_names = ['Vhalf', 'Kmax', 'Ku', 'N', 'EC50']

# # Define the model inputs
# problem = {
#     'num_vars': 5,
#     'names': param_names,
#     'bounds': [[-210, 0],
#                [1, 1e9],
#                [1e-5, 1],
#                [0.65, 1.304],
#                [4.206e1, 1e9]]
# }

# # Generate samples
# param_values = saltelli.sample(problem, 1024)

# # Run model (example)
# Y = Ishigami.evaluate(param_values)

# # Perform analysis
# Si = sobol.analyze(problem, Y, print_to_console=True)

# # Print the first-order sensitivity indices
# print(Si['S1'])