# Reference:
# Li Z, Dutta S, Sheng J, Tran PN, Wu W, Chang K, Mdluli T, Strauss DG,
# Colatsky T. Improving the In Silico Assessment of Proarrhythmia Risk by
# Combining hERG (Human Ether-Ã -go-go-Related Gene) Channel-Drug Binding
# Kinetics and Multichannel Pharmacology. Circ Arrhythm Electrophysiol.
# 2017 Feb;10(2):e004628. doi: 10.1161/CIRCEP.116.004628.


import myokit
import numpy as np

import modelling


class BindingKinetics(object):
    """
    To create a library of all the dynamic hERG-binding parameters for
    different drug compounds.
    (ref. Table 2)
    """

    def __init__(self, model, protocol=None, current_head=None):
        super(BindingKinetics, self).__init__()

        self.model = model
        self.protocol = protocol
        self.sim = myokit.Simulation(self.model, self.protocol)
        # self.sim.set_default_state(self.sim.state())
        self.initial_state = self.sim.state()

        if current_head is None:
            self.current_head = next(iter(self.model.states())).parent()
        else:
            self.current_head = self.model.get(current_head)
        # Save model's original constants
        self.original_constants = {
            "Vhalf": self.model.get(self.current_head.var('Vhalf')).eval(),
            "Kmax": self.model.get(self.current_head.var('Kmax')).eval(),
            "Ku": self.model.get(self.current_head.var('Ku')).eval(),
            "n": self.model.get(self.current_head.var('n')).eval(),
            "EC50": self.model.get(self.current_head.var('halfmax')).eval(),
            "Kt": self.model.get(self.current_head.var('Kt')).eval(),
            "gKr": self.model.get(self.current_head.var('gKr')).eval(), }

    def drug_simulation(self, drug, drug_conc, repeats,
                        timestep=0.1, save_signal=1):
        param_lib = modelling.BindingParameters()

        Vhalf = param_lib.binding_parameters[drug]['Vhalf']
        Kmax = param_lib.binding_parameters[drug]['Kmax']
        Ku = param_lib.binding_parameters[drug]['Ku']
        N = param_lib.binding_parameters[drug]['N']
        EC50 = param_lib.binding_parameters[drug]['EC50']

        t_max = self.protocol.characteristic_time()

        concentration = self.model.get('ikr.D')
        concentration.set_state_value(drug_conc)

        self.sim = myokit.Simulation(self.model, self.protocol)
        self.sim.reset()
        # self.sim.set_state(self.initial_state)

        self.sim.set_constant(self.current_head.var('Vhalf'), Vhalf)
        self.sim.set_constant(self.current_head.var('Kmax'), Kmax)
        self.sim.set_constant(self.current_head.var('Ku'), Ku)
        self.sim.set_constant(self.current_head.var('n'), N)
        self.sim.set_constant(self.current_head.var('halfmax'), EC50)
        self.sim.set_constant(self.current_head.var('Kt'), 3.5e-5)
        self.sim.set_constant(self.current_head.var('gKr'),
                              self.original_constants["gKr"])

        self.sim.pre(t_max * (repeats - save_signal))
        log = self.sim.run(t_max * save_signal, log_interval=timestep)
        d2 = log.npview()
        if save_signal > 1:
            d2 = d2.fold(t_max)

        # self.sim.reset()

        return d2

    def conductance_simulation(self, conductance, repeats,
                               timestep=1, save_signal=1):
        self.sim = myokit.Simulation(self.model, self.protocol)
        self.sim.reset()

        self.sim.set_constant(self.current_head.var('Vhalf'),
                              self.original_constants["Vhalf"])
        self.sim.set_constant(self.current_head.var('Kmax'),
                              self.original_constants["Kmax"])
        self.sim.set_constant(self.current_head.var('Ku'),
                              self.original_constants["Ku"])
        self.sim.set_constant(self.current_head.var('n'),
                              self.original_constants["n"])
        self.sim.set_constant(self.current_head.var('halfmax'),
                              self.original_constants["EC50"])
        self.sim.set_constant(self.current_head.var('Kt'),
                              self.original_constants["Kt"])
        self.sim.set_constant(self.current_head.var('gKr'), conductance)
        t_max = self.protocol.characteristic_time()

        self.sim.pre(t_max * (repeats - save_signal))
        log = self.sim.run(t_max * save_signal, log_interval=timestep)
        d2 = log.npview()
        if save_signal > 1:
            d2 = d2.fold(t_max)

        # self.sim.reset()

        return d2

    def state_occupancy_plot(self, ax, signal_log, pulse=None, legend=True):

        if pulse is None:
            ax.stackplot(signal_log.time(),
                         *[signal_log[s] for s in self.model.states()
                         if str(s.parent()) == 'ikr' and s.name() != 'D'],
                         labels=[s.name() for s in self.model.states()
                         if str(s.parent()) == 'ikr' and s.name() != 'D'],
                         zorder=-10)
        else:
            ax.stackplot(signal_log.time(),
                         *[signal_log[s, pulse] for s in self.model.states()
                         if str(s.parent()) == 'ikr' and s.name() != 'D'],
                         labels=[s.name() for s in self.model.states()
                         if str(s.parent()) == 'ikr' and s.name() != 'D'],
                         zorder=-10)

        ax.set_xlabel('Time (ms)')

        label_list = []
        for t in ax.get_legend_handles_labels():
            label_list.append(t)
        new_list = label_list[1]
        new_list = ['O*' if x == 'Obound' else x for x in new_list]
        new_list = ['C*' if x == 'Cbound' else x for x in new_list]
        new_list = ['IO*' if x == 'IObound' else x for x in new_list]

        if legend:
            ax.legend(ncol=2, handles=label_list[0], labels=new_list,
                      loc="lower right", handlelength=1, columnspacing=1,
                      labelspacing=0.3)
        ax.set_rasterization_zorder(0)

        return ax

    def extract_peak(self, signal_log, current_name):
        peaks = []
        pulses = len(signal_log.keys_like(current_name))

        if pulses == 0:
            peaks.append(np.max(signal_log[current_name]))

        for i in range(pulses):
            peaks.append(np.max(signal_log[current_name, i]))

        if peaks[0] == 0:
            peak_reduction = 0
        else:
            peak_reduction = (peaks[0] - peaks[-1]) / peaks[0]

        return peaks, peak_reduction

    def APD90(self, signal, offset, timestep):
        APA = max(signal) - min(signal)
        APD90_v = min(signal) + 0.1 * APA
        index = np.abs(np.array(signal) - APD90_v).argmin()
        APD90 = index * timestep - offset
        if APD90 < 1:
            APD90 = len(signal) * timestep

        return APD90
