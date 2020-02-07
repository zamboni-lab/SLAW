import common.references as cr
import common.tools as ct
import os
import common
import model.output_handler as oh


class IonAnnotater:
    ####Row is a row of the processing database and peak
    def __init__(self,hash,path_db,path_datamatrix,out,polarity,model,ncores,nfiles,
                 ppm,dmz,min_filter=2,adducts=None,main_adducts=None):
        self.hash = hash
        self.path_db = path_db
        self.path_datamatrix = path_datamatrix
        self.output_dm_reduced = os.path.join(out,cr.OUT["ANNOTATION"]+self.hash+"_reduced.csv")
        self.output_dm_full = os.path.join(out,cr.OUT["ANNOTATION"]+self.hash+"_full.csv")
        self.adducts = adducts
        self.main_adducts = main_adducts
        if not polarity in cr.DATA["IONMODE"]:
            raise Exception("Invalid ion mode argument. ",polarity)
        self.polarity=polarity
        self.model=model
        self.num_cores=ncores
        self.num_files= nfiles
        self.ppm=ppm
        self.dmz=dmz
        self.min_filter = min_filter

    ###Wrti both type of adduct to a temporry directory
    def write_adducts(self,ohv):
        path_add_full = ohv.getFile(cr.TEMP["IONANNOTATION"]["FULL"])
        path_add_main = ohv.getFile(cr.TEMP["IONANNOTATION"]["MAIN"])
        summary = open(path_add_full, "w+")
        tfiles = [self.adducts[i] + '\n' for i in range(len(self.adducts))]
        summary.writelines(tfiles)
        summary.close()

        summary = open(path_add_main, "w+")
        tfiles = [self.main_adducts[i]+'\n' for i in range(len(self.main_adducts))]
        summary.writelines(tfiles)
        summary.close()

    def need_computing(self):
        return not os.path.exists(self.output_dm_reduced)

    def get_output_datamatrices(self):
        return [self.output_dm_reduced,self.output_dm_full]

    def command_line(self,ohv):
        pscript = ct.find_rscript()
        pmatching = os.path.join(pscript,"cliques_matching.cpp")
        command_line = os.path.join(pscript,"annotating_mixed_method.R")
        ####We give all the name of the grouping parameters implicated in a single file
        cline = " ".join(["Rscript",command_line,self.path_datamatrix,self.path_db,self.output_dm_full,self.output_dm_reduced,
                 str(self.num_cores),common.references.DATA["IONANNOTATION"]["XCMS_MODEL"],ohv.getFile(cr.TEMP["IONANNOTATION"]["FULL"]),
                 ohv.getFile(cr.TEMP["IONANNOTATION"]["MAIN"]),self.polarity,str(self.ppm),str(self.dmz),
                          str(self.num_files),str(self.min_filter),pmatching])
        return cline
