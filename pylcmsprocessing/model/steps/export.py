import os
import common.tools as ct


class ExportMzTab:

    def __init__(self, path_db, row, mztab_format):
        self.db = path_db
        self.mztab_name = os.path.join(os.environ['OUTPUT'], "data_" + row[3] + ".mzTab")
        self.mztab_format = mztab_format

    def command_line(self):
        pscript = ct.find_rscript()
        command_line = os.path.join(pscript, "export_mzTab.R")
        ####We give all the name of the grouping parameters implicated in a single file
        return "".join(
            [os.environ["RscriptString"], ' "', command_line,'" "', self.db,'" "', self.mztab_name,'" "', self.mztab_format,'"'])


class ExportIsoToMGF:
    def __init__(self, dm, mgf):
       #  self.hash = row[3]
        self.dm = dm
        self.mgf = mgf

    def command_line(self):
        pscript = ct.find_rscript()
        command_line = os.path.join(pscript, "export_iso_to_mgf.R")
        ####We give all the name of the grouping parameters implicated in a single file
        return "".join(
            [os.environ["RscriptString"], ' "', command_line, '" "', self.dm, '" "', self.mgf,'"'])