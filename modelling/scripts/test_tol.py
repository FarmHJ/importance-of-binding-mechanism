# Replicate simulation error when tolerance value is small
import myokit
import pandas as pd

# Load problematic parameter combination
problem_file = '../../simulation_results/paramvalue_problem.csv'
param_values = pd.read_csv(problem_file, header=[0], index_col=[0])

# Set up AP model
APmodel = '../../model/ohara-cipa-v1-2017.mmt'
APmodel, _, x = myokit.load(APmodel)
protocol = myokit.pacing.blocktrain(1000, 1, offset=50)
tmax = protocol.characteristic_time()

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

# Set looser tolerance and run simulation
# This runs ok
sim.set_tolerance(abs_tol=1e-8, rel_tol=1e-4)
log = sim.run(tmax * 1000)

# Set smaller tolerance value and run simulation
# One pulse at a time
# This runs ok too
sim.reset()
sim.set_tolerance(abs_tol=1e-9, rel_tol=1e-10)
for pace in range(999):
    empty_log = sim.run(tmax, log=myokit.LOG_NONE)
    sim.set_time(0)
log = sim.run(tmax)

# Set smaller tolerance value and run simulation
# 1000 pulses in one run
# But not this one
sim.reset()
sim.set_tolerance(abs_tol=1e-9, rel_tol=1e-10)
log = sim.run(tmax * 5)
