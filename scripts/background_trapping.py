# Introduces the idea of trapping and justifies the use of the Milnes protocol
import myokit
import os

import modelling

drugs = ['dofetilide', 'verapamil']
drug_concs = [30, 1000]  # nM
short_label = ['drug_free', 'dofetilide', 'verapamil']

data_dir = '../simulation_data/background/'
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)

# Simulating the SD model with Milnes' protocol for 10 pulses after addition
# of a dofetilide-like drug and verapamil-like drug

# Load hERG model
model = '../math_model/ohara-cipa-v1-2017-IKr-opt.mmt'
model, _, x = myokit.load(model)
current_model = modelling.BindingKinetics(model)

# Define Milnes' protocol
pulse_time = 25e3
protocol = myokit.Protocol()
protocol.schedule(-80, 0, 800, period=pulse_time)
protocol.schedule(-90, 800, 100, period=pulse_time)
protocol.schedule(-80, 900, 100, period=pulse_time)
protocol.schedule(-80, 11000, 14000 - 1, period=pulse_time)
current_model.protocol = protocol
repeats = 10

abs_tol = 1e-7
rel_tol = 1e-8

# Simulate control condition
control_log = current_model.drug_simulation(drugs[0], 0, 1000, save_signal=10,
                                            abs_tol=abs_tol, rel_tol=rel_tol)
control_log.save_csv(data_dir + 'control_Milnes_current_pulses10.csv')

# Simulate 10 pulses after drug addition from steady state of control condition
control_log_single = current_model.drug_simulation(drugs[0], 0, 1000)
for i, conc_i in enumerate(drug_concs):

    log = current_model.drug_simulation(
        drugs[i], conc_i, repeats, save_signal=repeats,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'],
        set_state=control_log_single, abs_tol=abs_tol, rel_tol=rel_tol)
    log.save_csv(data_dir + drugs[i] + '_Milnes_current_pulses10.csv')

# Simulating the SD model with Milnes' protocol for 1000 pulses till steady
# state under drug free, addition of dofetilide-like drug and verapamil-like
# drug conditions

# Define Milnes' protocol
Milnes_protocol = modelling.ProtocolLibrary().Milnes(pulse_time)
current_model.protocol = Milnes_protocol

# Simulate SD model under drug-free condition
repeats = 1000
log_control = current_model.drug_simulation(drugs[0], 0, repeats,
                                            abs_tol=abs_tol, rel_tol=rel_tol)
log_control.save_csv(data_dir + short_label[0] + '_Milnes_current.csv')

# Simulate the SD model under addition of drugs
for d in range(len(drugs)):
    log = current_model.drug_simulation(drugs[d], drug_concs[d], repeats,
                                        abs_tol=abs_tol, rel_tol=rel_tol)
    log.save_csv(data_dir + short_label[d + 1] + '_Milnes_current.csv')

# Simulating the SD model with AP clamp protocol for 1000 pulses till steady
# state under drug free, addition of dofetilide-like drug and verapamil-like
# drug conditions

# Load AP model
APmodel = '../math_model/ohara-cipa-v1-2017-opt.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.BindingKinetics(APmodel, current_head='ikr')

# Define current pulse
pulse_time = 1000
protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)
AP_model.protocol = protocol

# Simulate AP for AP clamp protocol
APclamp = AP_model.drug_simulation(drugs[1], drug_concs[1], repeats,
                                   abs_tol=abs_tol, rel_tol=rel_tol)
APclamp.save_csv(data_dir + 'APclamp.csv')

# Set up AP clamp protocol
times = APclamp['engine.time']
voltages = APclamp['membrane.V']
tmax = times[-1] + 1

# Simulate SD model with AP clamp protocol under drug free condition
log_control = current_model.drug_APclamp(drugs[0], 0, times, voltages,
                                         tmax, repeats, abs_tol=abs_tol,
                                         rel_tol=rel_tol)
log_control.save_csv(data_dir + short_label[0] + '_APclamp_current.csv')

# Simulate the SD model under addition of drugs
for d in range(len(drugs)):
    log = current_model.drug_APclamp(drugs[d], drug_concs[d], times,
                                     voltages, tmax, repeats, abs_tol=abs_tol,
                                     rel_tol=rel_tol)
    log.save_csv(data_dir + short_label[d + 1] + '_APclamp_current.csv')
