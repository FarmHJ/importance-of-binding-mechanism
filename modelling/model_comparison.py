import numpy as np

import modelling


class ModelComparison(object):
    """
    To create a class to run the model comparison.
    """

    def __init__(self, drug_param_values):
        super(ModelComparison, self).__init__()

        # Need to create or check that there is a similar ForwardModel for
        # BindingKinetics
        # if not issubclass(BKmodel, modelling.ForwardModel):
        #     raise TypeError(
        #         'Model has to be a subclass of modelling.ForwardModel.')

        self.drug_param_values = drug_param_values
        param_names = ['Vhalf', 'Kmax', 'Ku', 'N', 'EC50']

        # Check if the parameters are correct
        if sorted(list(self.drug_param_values.columns)) != sorted(param_names):
            ValueError("Parameters are incorrect. The parameters are " +
                       "Vhalf, Kmax, Ku, N and EC50")

        self.Hill_model = modelling.HillsModel()
        self.optimiser = modelling.HillsModelOpt(self.Hill_model)

    def compute_Hill(self, BKmodel, drug_conc=None, steady_state_pulse=1000,
                     Hill_upper_thres=0.9, Hill_lower_thres=0.05):

        if drug_conc is None:
            drug_conc = list(np.append(0, 10**np.linspace(-1, 1, 5)))

        peaks = []
        for i in range(len(drug_conc)):
            log = BKmodel.custom_simulation(
                self.drug_param_values, drug_conc[i], steady_state_pulse,
                log_var=['engine.time', 'ikr.IKr'])
            peak, _ = BKmodel.extract_peak(log, 'ikr.IKr')
            peaks.append(peak[-1])

        peaks_norm = (peaks - min(peaks)) / (max(peaks) - min(peaks))

        # Make sure there are enough data points for the head of Hill curve
        data_pt_checker = [True if i > Hill_upper_thres else False
                           for i in peaks_norm]
        counter = 0
        while sum(data_pt_checker) < 3 and counter < 20:
            drug_conc.insert(1, drug_conc[1] / np.sqrt(10))
            log = BKmodel.custom_simulation(
                self.drug_param_values, drug_conc[1], steady_state_pulse,
                log_var=['engine.time', 'ikr.IKr'])
            peak, _ = BKmodel.extract_peak(log, 'ikr.IKr')
            peaks.insert(1, peak[-1])
            peaks_norm = (peaks - min(peaks)) / (max(peaks) - min(peaks))
            data_pt_checker = [True if i > Hill_upper_thres else False
                               for i in peaks_norm]
            counter += 1

        if counter == 20:
            return 'Hill curve did not form.', drug_conc, peaks_norm

        # Make sure there are enough data points for the tail of Hill curve
        data_pt_checker = [True if i < Hill_lower_thres else False
                           for i in peaks_norm]
        counter = 0
        while sum(data_pt_checker) < 3 and counter < 20:
            drug_conc = drug_conc + [max(drug_conc) * np.sqrt(10)]
            log = BKmodel.custom_simulation(
                self.drug_param_values, drug_conc[-1], steady_state_pulse,
                log_var=['engine.time', 'ikr.IKr'])
            peak, _ = BKmodel.extract_peak(log, 'ikr.IKr')
            peaks.append(peak[-1])
            peaks_norm = (peaks - min(peaks)) / (max(peaks) - min(peaks))
            data_pt_checker = [True if i < Hill_lower_thres else False
                               for i in peaks_norm]
            counter += 1

        if counter == 20:
            return 'Hill curve did not form.', drug_conc, peaks_norm

        # return 0, drug_conc, peaks_norm
        # Fit Hill curve
        Hill_curve, _ = self.optimiser.optimise(drug_conc, peaks_norm)

        return Hill_curve[:2], drug_conc, peaks_norm

    def APD_sim(self, AP_model, Hill_curve_coefs, drug_conc=None,
                steady_state_pulse=1000, save_signal=2, offset=50,
                data_points=20):

        base_conductance = AP_model.original_constants['gKr']
        APD_trapping = []
        APD_conductance = []
        if drug_conc is None:
            drug_conc = 10**np.linspace(-1, 5, data_points)

        for i in range(len(drug_conc)):
            # Run simulation for trapping model
            log = AP_model.custom_simulation(
                self.drug_param_values, drug_conc[i], steady_state_pulse,
                timestep=0.1, save_signal=save_signal, abs_tol=1e-7,
                rel_tol=1e-8, log_var=['engine.time', 'membrane.V'])

            # Compute APD90
            APD_trapping_pulse = []
            for pulse in range(save_signal):
                apd90 = AP_model.APD90(log['membrane.V', pulse], offset, 0.1)
                APD_trapping_pulse.append(apd90)
            APD_trapping.append(APD_trapping_pulse)

            # Run simulation for conductance model
            reduction_scale = self.Hill_model.simulate(
                Hill_curve_coefs, drug_conc[i])
            d2 = AP_model.conductance_simulation(
                base_conductance * reduction_scale, steady_state_pulse,
                timestep=0.1, save_signal=save_signal, abs_tol=1e-7,
                rel_tol=1e-8, log_var=['engine.time', 'membrane.V'])

            # Compute APD90
            APD_conductance_pulse = []
            for pulse in range(save_signal):
                apd90 = AP_model.APD90(d2['membrane.V', pulse], offset, 0.1)
                APD_conductance_pulse.append(apd90)
            APD_conductance.append(APD_conductance_pulse)

        APD_trapping = [max(i) for i in APD_trapping]
        APD_conductance = [max(i) for i in APD_conductance]

        return APD_trapping, APD_conductance

    def compute_RMSE(self, APD_trapping, APD_conductance):

        RMSError = np.sum((np.array(APD_trapping) -
                          np.array(APD_conductance))**2) / len(APD_trapping)
        RMSError = np.sqrt(RMSError)

        return RMSError

    def compute_MAE(self, APD_trapping, APD_conductance):

        MAError = np.sum(np.array(APD_trapping) - np.array(APD_conductance)) \
            / len(APD_trapping)

        return MAError
