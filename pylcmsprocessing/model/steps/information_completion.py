import os
import pandas as pd
import common.tools as ct

# class InformationExpander:
#     def __init__(self,db,dm,temp_dm,rt_model,isotopes,max_iso,max_charge,quant,margin_mz,ppm,dmz,num_files,num_workers):
#         self.max_charge = 2
#         self.db = db
#         self.dm = dm
#         self.temp_dm = temp_dm
#         self.max_iso = max_iso
#         self.max_charge = max_charge
#         self.model = rt_model
#         self.isotopes = isotopes
#         self.quant = quant
#         self.margin_mz = margin_mz
#         self.ppm = ppm
#         self.dmz = dmz
#         self.num_files = num_files
#         self.num_workers = num_workers
#
#     def need_computing(self):
#         if not os.path.isfile(self.dm):
#             return False
#         else:
#             tab = pd.read_csv(self.dm,nrows=10)
#             if "raw_isotopic_pattern" in tab.columns:
#                 return False
#             else:
#                 return True
#
#     def get_output(self):
#         return self.output
#
#     def command_line(self):
#         pscript = ct.find_rscript()
#         command_line = os.path.join(pscript,"gap_filling.R")
#         return " ".join([os.environ["RscriptString"]," ",command_line,'"'+self.db+'"','"'+self.dm+'"',
#                          '"'+self.temp_dm+'"','"'+self.model+'"','"'+self.isotopes+'"',
#                          str(self.max_iso),str(self.max_charge),self.quant,str(self.margin_mz),
#                          str(self.ppm),str(self.dmz),str(self.num_files),str(self.num_workers)])

class InformationExpanderRefactored:
    def __init__(self,db,dm,temp,rt_model,hdf5,isotopes,max_iso,max_charge,quant,ppm,dmz,num_files,num_workers):
        self.max_charge = 2
        self.db = db
        self.dm = dm
        self.temp = temp
        self.hdf5 = hdf5
        self.max_iso = max_iso
        self.max_charge = max_charge
        self.model = rt_model
        self.isotopes = isotopes
        self.quant = quant
        self.num_files = num_files
        self.num_workers = num_workers
        self.ppm = ppm
        self.dmz = dmz

    def need_computing(self):
        if not os.path.isfile(self.dm):
            return True
        else:
            tab = pd.read_csv(self.dm,nrows=10)
            if "isotopic_pattern_abs" in tab.columns:
                return False
            else:
                return True

    def get_output(self):
        return self.output

    def command_line(self):
        pscript = ct.find_rscript()
        command_line = os.path.join(pscript,"gap_filling_refactored.R")
        return " ".join([os.environ["RscriptString"]," ",command_line,'"'+self.db+'"','"'+self.dm+'"',
                         '"'+self.temp+'"','"'+self.model+'"','"'+self.hdf5+'"','"'+self.isotopes+'"',self.quant,
                         str(self.max_iso),str(self.max_charge),str(self.ppm),str(self.dmz),str(self.num_files),"0",str(self.num_workers)])


