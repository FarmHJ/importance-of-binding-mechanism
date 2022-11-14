# Replicate simulation error when tolerance value is small

import myokit
import pandas as pd

saved_data_dir = '../../simulation_results/'

param_values = pd.read_csv(saved_data_dir + 'paramvalue_problem.csv',
                           header=[0], index_col=[0])

# Set up AP model
APmodel = '../../model/ohara-cipa-v1-2017.mmt'
APmodel, _, x = myokit.load(APmodel)
protocol = myokit.pacing.blocktrain(1000, 1, offset=50)

# Input drug concentration
APmodel.get('ikr.D').set_state_value(param_values['drug_conc'].values[0])

# Input drug parameters
sim = myokit.Simulation(APmodel, protocol)
sim.set_constant('ikr.Vhalf', param_values['Vhalf'].values[0])
sim.set_constant('ikr.Kmax', param_values['Kmax'].values[0])
sim.set_constant('ikr.Ku', param_values['Ku'].values[0])
sim.set_constant('ikr.n', param_values['N'].values[0])
sim.set_constant('ikr.halfmax', param_values['EC50'].values[0])
sim.set_constant('ikr.Kt', 3.5e-5)

# Set tolerance and run simulation
sim.set_tolerance(abs_tol=1e-8, rel_tol=1e-4)
log = sim.run(1000 * 1000)

# Set smaller tolerance value and run simulation
sim.reset()
sim.set_tolerance(abs_tol=1e-8, rel_tol=1e-10)
log = sim.run(1000 * 1000)
