# Reference:
# Li Z, Dutta S, Sheng J, Tran PN, Wu W, Chang K, Mdluli T, Strauss DG,
# Colatsky T. Improving the In Silico Assessment of Proarrhythmia Risk by
# Combining hERG (Human Ether-Ã -go-go-Related Gene) Channel-Drug Binding
# Kinetics and Multichannel Pharmacology. Circ Arrhythm Electrophysiol.
# 2017 Feb;10(2):e004628. doi: 10.1161/CIRCEP.116.004628.

import numpy as np

import modelling


class BindingParameters(object):
    """
    To create a library of all the dynamic hERG-binding parameters for
    different drug compounds.
    (ref. Table 2)
    """

    def __init__(self):
        super(BindingParameters, self).__init__()

        self.drug_compounds = ['dofetilide', 'bepridil', 'terfenadine',
                               'cisapride', 'verapamil', 'ranolazine',
                               'quinidine', 'sotalol', 'chlorpromazine',
                               'ondansetron', 'diltiazem', 'mexiletine']
        self.binding_parameters = {
            'dofetilide': {
                'Kmax': 1e8,
                'Ku': 1.79e-5,
                'EC50': 5.483e8,
                'N': 0.9999,
                'Vhalf': -1.147,
                'Cmax': 2
            },
            'bepridil': {
                'Kmax': 3.735e7,
                'Ku': 1.765e-4,
                'EC50': 1e9,
                'N': 0.9365,
                'Vhalf': -54.93,
                'Cmax': 33
            },
            'terfenadine': {
                'Kmax': 9884,
                'Ku': 8.18e-5,
                'EC50': 41380,
                'N': 0.65,
                'Vhalf': -77.49,
                'Cmax': 4
            },
            'cisapride': {
                'Kmax': 9.997,
                'Ku': 4.161e-4,
                'EC50': 42.06,
                'N': 0.9728,
                'Vhalf': -199.5,
                'Cmax': 2.6
            },
            'cisapride_kmax_verapamil': {
                'Kmax': 4.646e4,
                'Ku': 4.161e-4,
                'EC50': 42.06,
                'N': 0.9728,
                'Vhalf': -199.5,
                'Cmax': 2.6
            },
            'cisapride_kmax_bepridil': {
                'Kmax': 3.735e7,
                'Ku': 4.161e-4,
                'EC50': 42.06,
                'N': 0.9728,
                'Vhalf': -199.5,
                'Cmax': 2.6
            },
            'cisapride_EC50_ranolazine': {
                'Kmax': 9.997,
                'Ku': 4.161e-4,
                'EC50': 1.472e5,
                'N': 0.9728,
                'Vhalf': -199.5,
                'Cmax': 2.6
            },
            'cisapride_EC50_verapamil': {
                'Kmax': 9.997,
                'Ku': 4.161e-4,
                'EC50': 9.184e6,
                'N': 0.9728,
                'Vhalf': -199.5,
                'Cmax': 2.6
            },
            'cisapride_kmax_EC50_verapamil': {
                'Kmax': 4.646e4,
                'Ku': 4.161e-4,
                'EC50': 9.184e6,
                'N': 0.9728,
                'Vhalf': -199.5,
                'Cmax': 2.6
            },

            'verapamil': {
                'Kmax': 4.646e4,
                'Ku': 7.927e-4,
                'EC50': 9.184e6,
                'N': 1.043,
                'Vhalf': -100,
                'Cmax': 81
            },
            'ranolazine': {
                'Kmax': 55.84,
                'Ku': 1.929e-2,
                'EC50': 1.472e5,
                'N': 0.95,
                'Vhalf': -94.87,
                'Cmax': 1948.2
            },
            'mexiletine': {
                'Kmax': 9.996,
                'Ku': 9.967e-2,
                'EC50': 2.308e6,
                'N': 1.304,
                'Vhalf': -86.26,
                'Cmax': 4129
            },
            'quinidine': {
                'Kmax': 5770,
                'Ku': 1e-2,
                'EC50': 1e6,
                'N': 0.8311,
                'Vhalf': -64.87,
                'Cmax': 3237
            },
            'sotalol': {
                'Kmax': 2403,
                'Ku': 1.985e-2,
                'EC50': 9.619e6,
                'N': 0.7516,
                'Vhalf': -55,
                'Cmax': 14690
            },
            'chlorpromazine': {
                'Kmax': 206000,
                'Ku': 3.866e-2,
                'EC50': 5.677e7,
                'N': 0.8871,
                'Vhalf': -14.57,
                'Cmax': 38
            },
            'ondansetron': {
                'Kmax': 33540,
                'Ku': 2.325e-2,
                'EC50': 9.95e6,
                'N': 0.8874,
                'Vhalf': -82.11,
                'Cmax': 139
            },
            'diltiazem': {
                'Kmax': 251,
                'Ku': 2.816e-1,
                'EC50': 1e6,
                'N': 0.9485,
                'Vhalf': -90.89,
                'Cmax': 122
            },
        }
        self.Hill_curve = {
            'dofetilide': {
                'Hill_coef': 1.002,
                'IC50': 6.187}}


class ProtocolParameters(object):
    """
    To create a library of all the parameters for different protocols,
    especially default pulse time and function name.
    """

    def __init__(self):
        super(ProtocolParameters, self).__init__()

        self.protocols = ['Milnes', 'Pneg80', 'P0', 'P40']
        self.protocol_parameters = {
            'Milnes': {
                'pulse_time': 25e3,
                'function': modelling.ProtocolLibrary().Milnes(25e3),
                'voltage_points': [-80, 0],
            },
            'Pneg80': {
                'pulse_time': 5400,
                'function': modelling.ProtocolLibrary().Pneg80(5400),
                'voltage_points': [-80, -50, 20],
            },
            'P0': {
                'pulse_time': 5400,
                'function': modelling.ProtocolLibrary().P0(5400),
                'voltage_points': [-80, -60, 0],
            },
            'P40': {
                'pulse_time': 5400,
                'function': modelling.ProtocolLibrary().P40(5400),
                'voltage_points': [-80, -60, 40],
            },
        }


class DrugConcentrations(object):
    """
    Create a library for default list of drug concentrations
    """
    def __init__(self):
        super(DrugConcentrations, self).__init__()

        self.drug_compounds = ['dofetilide', 'bepridil', 'terfenadine',
                               'cisapride', 'verapamil', 'ranolazine',
                               'quinidine', 'sotalol', 'chlorpromazine',
                               'ondansetron', 'diltiazem', 'mexiletine']
        self.drug_concentrations = {
            'dofetilide': {
                'coarse': [0, 0.1, 1, 10, 30, 100, 300, 500, 1000],
                'fine': 10.0**np.linspace(-1, 2.5, 20),
                'lit_default': [1, 3, 10, 30]
            },
            'verapamil': {
                'coarse': [0, 0.1, 1, 30, 300, 500, 1000, 10000, 1e5],
                'fine': 10.0**np.linspace(-1, 5, 20),
                'lit_default': [30, 100, 300, 1000]
            },
            'bepridil': {
                'coarse': [0, 0.1, 1, 10, 30, 100, 300, 500, 1000, 10000],
                'fine': 10.0**np.linspace(-1, 3.5, 20),
                'lit_default': [10, 30, 100, 300]
            },
            'terfenadine': {
                'coarse': [0, 0.1, 1, 10, 30, 100, 300, 500, 1000, 10000],
                'fine': 10.0**np.linspace(-1, 3.5, 20),
                'lit_default': [3, 10, 30, 100]
            },
            'cisapride': {
                'coarse': [0, 0.1, 1, 10, 30, 100, 300, 500, 1000, 10000],
                'fine': 10.0**np.linspace(-1, 3, 20),
                'lit_default': [1, 10, 100, 300]
            },
            'ranolazine': {
                'coarse': [0, 1, 30, 300, 500, 1000, 10000, 1e5, 1e6],
                'fine': 10.0**np.linspace(1, 5.5, 20),
                'lit_default': [1000, 1e4, 3e4, 1e5]
            },
            'quinidine': {
                'coarse': [0, 1, 30, 300, 500, 1000, 10000, 1e5, 1e6],
                'fine': 10.0**np.linspace(-1, 5.5, 20),
                'lit_default': [100, 300, 1000, 10000]
            },
            'sotalol': {
                'coarse': [0, 1, 30, 300, 500, 1000, 10000, 1e5, 1e6],
                'fine': 10.0**np.linspace(-1, 5.5, 20),
                'lit_default': [1e4, 3e4, 1e5, 3e5]
            },
            'chlorpromazine': {
                'coarse': [0, 1, 30, 300, 500, 1000, 3000, 10000, 1e5],
                'fine': 10.0**np.linspace(-1, 4.5, 20),
                'lit_default': [100, 300, 1000, 3000]
            },
            'ondansetron': {
                'coarse': [0, 1, 30, 300, 500, 1000, 10000, 1e5, 1e6],
                'fine': 10.0**np.linspace(-1, 5.5, 20),
                'lit_default': [300, 1000, 3000, 1e4]
            },
            'diltiazem': {
                'coarse': [0, 1, 30, 300, 500, 1000, 10000, 1e5, 1e6],
                'fine': 10.0**np.linspace(-1, 5.5, 20),
                'lit_default': [3000, 1e4, 3e4, 1e5]
            },
            'mexiletine': {
                'coarse': [0, 1, 30, 300, 500, 1000, 10000, 1e5, 1e6],
                'fine': 10.0**np.linspace(-1, 5.5, 20),
                'lit_default': [1e4, 3e4, 1e5, 3e5]
            },
        }
