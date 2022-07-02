# Reference: 
# Li Z, Dutta S, Sheng J, Tran PN, Wu W, Chang K, Mdluli T, Strauss DG, Colatsky T. 
# Improving the In Silico Assessment of Proarrhythmia Risk by Combining hERG (Human Ether-à-go-go-Related Gene) 
# Channel-Drug Binding Kinetics and Multichannel Pharmacology. Circ Arrhythm Electrophysiol. 
# 2017 Feb;10(2):e004628. doi: 10.1161/CIRCEP.116.004628. 


import matplotlib.pyplot as plt
import myokit
import numpy as np

import modelling


class BindingKinetics(object):
    """
    To create a library of all the dynamic hERG-binding parameters for
    different drug compounds.
    (ref. Table 2)
    """

    def __init__(self, model):
        super(BindingKinetics, self).__init__()

        self.model = model

    def drug_simulation(self, drug, drug_conc, repeats, timestep=0.1):
        param_lib = modelling.BindingParameters()

        Vhalf = param_lib.binding_parameters[drug]['Vhalf']
        Kmax = param_lib.binding_parameters[drug]['Kmax']
        Ku = param_lib.binding_parameters[drug]['Ku']
        N = param_lib.binding_parameters[drug]['N']
        EC50 = param_lib.binding_parameters[drug]['EC50']

        t_max = 25e3
        protocol = self.Milnes_protocol(t_max)

        concentration = self.model.get('IKr.D')
        concentration.set_state_value(drug_conc)
        
        sim = myokit.Simulation(self.model, protocol)
        sim.reset()
        current_head = next(iter(self.model.states())).parent()
        sim.set_constant(current_head.var('Vhalf'), Vhalf)
        sim.set_constant(current_head.var('Kmax'), Kmax)
        sim.set_constant(current_head.var('Ku'), Ku)
        sim.set_constant(current_head.var('n'), N)
        sim.set_constant(current_head.var('halfmax'), EC50)
        sim.set_constant(current_head.var('Kt'), 3.5e-5)
        
        log = sim.run(t_max * repeats, log_interval=timestep)
        sim.reset()

        return log

    def conductance_simulation(self, conductance, repeats):
        t_max = 25e3
        protocol = self.Milnes_protocol(t_max)
        
        sim = myokit.Simulation(self.model, protocol)
        sim.reset()
        current_head = next(iter(self.model.states())).parent()
        sim.set_constant(current_head.var('GKr_b'), conductance)

        log = sim.run(t_max * repeats, log_interval=1)
        sim.reset()

        return log

    def Milnes_protocol(self, t_max):
        protocol = myokit.Protocol()
        protocol.schedule(-80, 0, 800, t_max)
        protocol.schedule(-90, 800, 100, t_max)
        protocol.schedule(-80, 900, 100, t_max)
        protocol.schedule(-80, 11000, 14000 - 2, t_max)

        return protocol


    def state_occupancy_plot(self, ax, signal_log, pulse):
        ax.stackplot(signal_log.time(),
                *[signal_log[s, pulse] for s in self.model.states() if s.name() != 'D'],
                labels = [s.name() for s in self.model.states() if s.name() != 'D'],
                zorder=-10)
        ax.set_xlabel('Time (ms)')
        ax.set_ylabel('Occupancy')
        ax.legend()
        ax.set_rasterization_zorder(0)
        
        return ax


    def extract_peak(self, signal_log, current_name):
        peaks = []
        pulses = len(signal_log.keys_like(current_name))
        for i in range(pulses):
            peaks.append(np.max(signal_log[current_name, i]))
       
        if peaks[0] == 0:
            peak_reduction = 0
        else:
            peak_reduction = (peaks[0] - peaks[-1]) / peaks[0]
        
        return peaks, peak_reduction

