## Not complete
# Simple checking if simulated current matches literature
# In particular, check if the input parameters are correct

import myokit
import sys

import modelling


plot_fig = True

drug = sys.argv[1]
protocol_name = 'Milnes'
protocol_params = modelling.ProtocolParameters()
pulse_time = protocol_params.protocol_parameters[protocol_name]['pulse_time']
protocol = protocol_params.protocol_parameters[protocol_name]['function']
repeats = 10
drug_conc = modelling.DrugConcentrations().drug_concentrations[
    drug]['lit_default']
drug_labels = [str(i) + ' nM' for i in drug_conc]

testing_fig_dir = '../../figures/testing/'
final_fig_dir = '../../figures/input_check/' + drug + '/'

saved_fig_dir = testing_fig_dir

# Load IKr model
model = '../../model/ohara-cipa-v1-2017-IKr.mmt'
model, _, x = myokit.load(model)

drug_model = modelling.BindingKinetics(model)
drug_model.protocol = protocol

# Run simulation
total_log = []
peaks = []
for i in range(len(drug_conc)):
    log = drug_model.drug_simulation(
        drug, drug_conc[i], repeats, save_signal=repeats,
        log_var=['engine.time', 'membrane.V', 'ikr.IKr'])
    peak, _ = drug_model.extract_peak(log, 'ikr.IKr')
    print(log['ikr.IKr', 9])
    for j in range(repeats):
        log['membrane.V', j]
    # peaks.append(peak[-1])
    # total_log.append(log)

# Plot drug response against drug concentration curve
# peaks = (peaks - min(peaks)) / (max(peaks) - min(peaks))

# plt.rcParams.update({'font.size': 9})

# plt.figure(figsize=(4, 3))
# plt.plot(np.log(drug_conc), peaks, 'o')
# plt.xlabel('Drug concentration (log)')
# plt.ylabel('Normalised peak current')
# plt.tight_layout()
# plt.savefig(saved_fig_dir + "peak_hERG_trapping_" + drug + "_concs.pdf")

# # peak reduction
# hERG_reduc_trap = []
#     hERG_reduc_conduct = []
# hERG_peak_conc_trap = hERG_peak_trapping[i]
# hERG_reduc_trap.append((hERG_peak_conc_trap[0]
#                         - hERG_peak_conc_trap[-1])
#                         / hERG_peak_conc_trap[0])

# hERG_peak_conc_conduct = hERG_peak_conductance[i]
# hERG_reduc_conduct.append((hERG_peak_conc_conduct[0]
#                            - hERG_peak_conc_conduct[-1])
#                            / hERG_peak_conc_conduct[0])

#     plt.figure()
#     plt.plot(hERG_reduc_trap)
#     plt.plot(hERG_reduc_conduct)
#     plt.savefig(saved_fig_dir + 'peak_reduction.pdf')