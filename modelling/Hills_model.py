import numpy as np
import pints


class HillsModel(pints.ForwardModel):
    """
    To infer Hills' coefficient and IC50 from simulated data
    """

    def __init__(self):
        super(HillsModel, self).__init__()

    def n_outputs(self):
        # Returns number of model outputs, i.e. response (metric)
        return 1

    def n_parameters(self):
        # Returns number of parameters to be inferred
        # i.e. Hill's coefficient and IC50
        return 2

    def simulate(self, params, drug_conc):

        # if not any([i == 0 for i in drug_conc]):
        #     raise ValueError("Must not have zero drug concentration.")

        Hills_coef = params[0]
        IC50 = params[1]

        response = IC50 / (np.power(drug_conc, Hills_coef) + IC50)

        return response


class HillsModelOpt(object):

    def __init__(self, model):
        super(HillsModelOpt, self).__init__()

        self.model = model

    def optimise(self, drug_conc, inhibit_metric):

        if not any([i == 0 for i in drug_conc]):
            raise ValueError("Must have drug concentration = 0")

        # normalise metric
        metric_min = min(inhibit_metric)
        metric_max = max(inhibit_metric)
        inhibit_metric = (
            inhibit_metric - metric_min) / (metric_max - metric_min)

        zero_drug_conc = drug_conc.index(0)
        print(zero_drug_conc)
        drug_conc = np.delete(drug_conc, zero_drug_conc)
        inhibit_metric = np.delete(inhibit_metric, zero_drug_conc)
        # del drug_conc[zero_drug_conc], inhibit_metric[zero_drug_conc]

        # optimisation
        problem = pints.SingleOutputProblem(self.model, drug_conc,
                                            inhibit_metric)
        log_likelihood = pints.GaussianLogLikelihood(problem)

#         log_prior = pints.ComposedLogPrior(
#                 pints.UniformLogPrior(0.5, 1),
#                 pints.UniformLogPrior(5, 15),
#                 pints.TruncatedGaussianLogPrior(0, 0.25, 0, np.inf))

        transform = pints.LogTransformation(n_parameters=3)
#         initial_parameters = log_prior.sample()
        initial_parameters = [0.9, 10, 0.1]
        optimiser = pints.OptimisationController(
            function=log_likelihood,
            x0=initial_parameters,
            method=pints.CMAES,
            transformation=transform)
        optimiser.set_parallel(True)

        optimiser.set_max_iterations(1000)
        param_best, score_best = optimiser.run()

        return param_best, score_best

# class Metrics - create a class to extract all possible metrics,
# e.g. avg normalised current, peak, max normalised peak(?)
