import myokit
import os
import pandas as pd


class DatasetLibrary(object):
    """
    A data library class that reads the experimental data
    for different protocols and different drugs
    """
    def __init__(self):
        super(DatasetLibrary, self).__init__()

        self._directory = os.path.join(
            os.path.dirname(os.path.dirname(
                os.path.dirname(__file__))),
            "220122_exp_data")

        self.protocol_choice = ["CIPA", "Pharm"]
        self.drug_choice = ["Cisapride", "Dofetilide", "Verapamil"]
        self.protocol_title = {"CIPA": "CiPA protocol",
                               "Pharm": "Roche's protocol"}
        self.compound_name = {
            "Cisapride": ["19"],
            "Dofetilide": ["110", "RO0319253-000-001"],
            "Verapamil": ["13"]}

    def exp_data_list(self, protocol, drug):
        """
        Returns file directories for different cells of given
        protocol and drug.
        """

        # check protocol choice
        if protocol not in self.protocol_choice:
            raise ValueError(
                "Choice of protocol must be either 'CIPA' or 'Pharm'")

        if drug not in self.drug_choice:
            raise ValueError(
                "Choice of drug must be one of 'Cisapride', 'Dofetilide' \
                    and 'Verapamil'")

        filepath = os.path.join(self._directory, protocol, drug)
        file_cells = os.listdir(filepath)

        # Find name of cells and their file path
        files = []
        cells = []
        for filename in file_cells:
            if filename.endswith(".csv") and not filename.startswith("OA"):
                cell_name = filename[-8:-5]
                cells.append(cell_name)
                files.append(os.path.join(filepath, filename))
            if filename.startswith("DMSO"):
                conc_file = os.path.join(filepath, filename)
                conc_info = pd.read_excel(conc_file)

        cell_name_path = pd.DataFrame(data={
            "cells": cells, "file_path": files,
            "drug_concentration": [0] * len(cells)})

        # Get drug concentration used for each cell
        for cell_num in cells:
            conc = conc_info.loc[conc_info["WELL ID"] == cell_num,
                                 "Concentration (Mol/L)"].iloc[0]
            cell_name_path.loc[cell_name_path["cells"] == cell_num,
                               "drug_concentration"] = conc

        return cell_name_path

    def exp_data_read(self, filepath):
        """
        Returns the dataframe of experimental data
        of given file path
        """
        data = pd.read_csv(filepath, header=2, sep=';')

        return data

    def detail_read(self, protocol, drug):
        """
        Returns parameters of the experiment.
        Including added compound name, its concentration,
        capacitance, seal resistance and series resistance.
        """
        # check protocol choice
        if protocol not in self.protocol_choice:
            raise ValueError(
                "Choice of protocol must be either 'CIPA' or 'Pharm'")

        if drug not in self.drug_choice:
            raise ValueError(
                "Choice of drug must be one of 'Cisapride', 'Dofetilide' \
                    and 'Verapamil'")

        filepath = os.path.join(self._directory, protocol, drug)
        file_cells = os.listdir(filepath)
        for filename in file_cells:
            if filename.endswith("_processed.csv"):
                file_dir = os.path.join(filepath, filename)
                data = pd.read_csv(file_dir, index_col=0)

        return data


class ProtocolLibrary(object):
    """
    A library class with known protocols
    """
    def __init__(self):
        super(ProtocolLibrary, self).__init__()

    def Milnes(self, t_max):
        protocol = myokit.Protocol()
        protocol.schedule(-80, 0, 800, period=t_max)
        protocol.schedule(-90, 800, 100, period=t_max)
        protocol.schedule(-80, 900, 100, period=t_max)
        protocol.schedule(-80, 11000, 13999, period=t_max)

        return protocol

    def current_impulse(self, t_max):
        return myokit.pacing.blocktrain(t_max, 1, offset=50)

    def validation3(self, t_max):
        protocol = myokit.Protocol()
        protocol.schedule(-80, 0, 100)  # , period=t_max)
        protocol.schedule(-40, 100, 50)  # , period=t_max)
        protocol.schedule(20, 150, 500)  # , period=t_max)
        protocol.schedule(-40, 650, 500)  # , period=t_max)
        protocol.schedule(-80, 1150, 200)  # , period=t_max)

        return protocol

    def hERG_validation(self, t_max):
        protocol = myokit.Protocol()
        protocol.schedule(-80, 0, 100)
        protocol.schedule(20, 100, 900)
        protocol.schedule(-80, 1000, 1000)

        return protocol

    def Pneg80(self, t_max):
        protocol = myokit.Protocol()
        protocol.schedule(-80, 0, 200, period=t_max)
        protocol.schedule(20, 200, 500, period=t_max)
        protocol.schedule(-50, 700, 200, period=t_max)
        protocol.schedule(-80, 900, 4500 - 1, period=t_max)

        return protocol

    def P0(self, t_max):
        protocol = myokit.Protocol()
        protocol.schedule(-80, 0, 100, period=t_max)
        protocol.schedule(-60, 5100, 200, period=t_max)
        protocol.schedule(-80, 5300, 100 - 1, period=t_max)

        return protocol

    def P40(self, t_max):
        protocol = myokit.Protocol()
        protocol.schedule(-80, 0, 100, period=t_max)
        protocol.schedule(40, 100, 5000, period=t_max)
        protocol.schedule(-60, 5100, 200, period=t_max)
        protocol.schedule(-80, 5300, 100 - 1, period=t_max)

        return protocol
