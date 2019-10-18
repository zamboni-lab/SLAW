import common.references as cr
import common.tools as ct
import os
import common

class Grouper:
    ####Row is a row of the processing database and peak
    def __init__(self,row,temp,out,files,names,intensity,mztol,rttol):
        self.hash = row[3]
        self.temp = temp
        self.output_data=os.path.join(out,"datamatrix_"+self.hash+".csv")
        self.output_idx=os.path.join(out,"index_"+self.hash+".txt")
        self.files=files
        self.names=names
        self.intensity=intensity
        self.mztol=mztol
        self.rttol=rttol
        self.path_files =os.path.join(temp,"files.txt")


    def need_computing(self):
        return not os.path.exists(self.output_data)

    def make_temp_file(self):
        summary = open(self.path_files, "w+")
        tfiles = [self.files[i]+','+self.names[i]+'\n' for i in range(len(self.files))]
        summary.writelines(tfiles)

    def get_output_datamatrix(self):
        return self.output_data

    def get_output_index(self):
        return self.output_idx

    def get_temp(self):
        return self.temp

    def command_line(self):
        pscript = ct.find_rscript()
        command_line = os.path.join(pscript,"construct_peaktable.R")
        ####We give all the name of the grouping parameters implicated in a single file
        return " ".join(["Rscript",command_line,self.path_files,self.output_data,
                         self.output_idx,self.intensity,str(self.rttol),str(self.mztol)])
