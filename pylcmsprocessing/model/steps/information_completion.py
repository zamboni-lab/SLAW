import os
import pandas as pd
import common.tools as ct


# args <- commandArgs(trailingOnly = TRUE)
# PATH_DB <- args[1]
# PATH_DM <- args[2]
# PATH_FILLED <- args[3]
# PATH_MODEL <- args[4]
# PATH_ISOTOPES <- args[5]
# MAX_ISO <- as.integer(args[6])
# QUANT <- args[7]
# MARGIN_MZ <- as.numeric(args[8])
# MAX_CHARGE <- as.integer(args[9])
# PPM <- as.numeric(args[10])
# DMZ <- as.numeric(args[11])

class InformationExpander:
    def __init__(self,db,dm,temp_dm,rt_model,isotopes,max_iso,max_charge,quant,margin_mz,ppm,dmz,num_workers):
        self.max_charge = 2
        self.db = db
        self.dm = dm
        self.temp_dm = temp_dm
        self.max_iso = max_iso
        self.max_charge = max_charge
        self.model = rt_model
        self.isotopes = isotopes
        self.quant = quant
        self.margin_mz = margin_mz
        self.ppm = ppm
        self.dmz = dmz
        self.num_workers = num_workers

    def need_computing(self):
        if not os.path.isfile(self.dm):
            return False
        else:
            tab = pd.read_csv(self.dm,nrows=10)
            if "raw_isotopic_pattern" in tab.columns:
                return False
            else:
                return True

    def get_output(self):
        return self.output

    def command_line(self):
        pscript = ct.find_rscript()
        command_line = os.path.join(pscript,"gap_filling.R")
        return " ".join(["Rscript",command_line,self.db,self.dm,
                         self.temp_dm,self.model,self.isotopes,
                         str(self.max_iso),str(self.max_charge),self.quant,str(self.margin_mz),
                         str(self.ppm),str(self.dmz),str(self.num_workers)])



