from model.helper.parameters_handler import ParametersFileHandler
import common.references as cr
import common.tools as ct
import model.helper.parallel_runner as pr
import os
from shutil import copyfile
import logging


#       - gap-filled data matrix
#       - data matrix
#       - mztab
#     type: str
# ms2:
#     priority: LOW
#     description: String. What MS2 information should be output, 'fused mgf with isotopes' refer to MGF file containing MS-MS under the MS2 tag and isotopic patterns under the MS1 tag. 'fused mgf' output a single mgf with one MS-MS spectrum by feature.
#     value: fused mgf
#     multiset:
#       - fused mgf with isotopes
#       - fused mgf
# fused mgf is always generated


def make_mztab_name(hash):
    return "slaw_mztab_"+hash+".mzTab"

class OutputFormatter:
    """This class format an output given a ParametersHandler object"""
    def __init__(self,ph):
        if not isinstance(ph,ParametersFileHandler):
            raise TypeError("'ph' should be a ParametersFileHandler object.")

        self.ms1_format = ph["output_format__ms1"]["value"]
        self.ms2_format = ph["output_format__ms2"]["value"]

    #To ensure correctness, all these steps shold derive information from the correct time object.
    def output_grouping(self,path_dm):
        """If necessary add a copy of the non filled data table before going further"""
        if "gap-filled data matrix" in self.ms1_format:
            #We copy the non filled data matrix
            dir_dm = os.path.dirname(path_dm)
            base_dm = os.path.basename(path_dm)
            gap_filled_path = os.path.join(dir_dm,"non_filled_"+base_dm)
            _ = copyfile(path_dm,gap_filled_path)


    def output_msms(self,path_full,fused_mgf):
        """If necessary, this add the MS1 isotopic pattern to the MGF file"""
        if "fused mgf with isotopes" in self.ms2_format:
            #In this case we have to remap the MGF to the MS2 datasets.


    def output_mztab(self,outdir,hash):
        """If necessary convert the table to an mzTab format"""
        if "mztab" in self.ms1_format:
            logging.info("Converting SLAW output to mzTab")
            omztab = os.path.join(outdir,make_mztab_name(hash))
            #We call the Rscript for conversion
            sconv = os.path.join(ct.find_rscript(),"convert_output_SLAW_to_mzTab.R")+" "+outdir+" "+omztab
            runner = pr.ParallelRunner(1)
            runner.run(sconv, silent=True)





