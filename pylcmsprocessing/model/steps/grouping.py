import common.references as cr
import common.tools as ct
import os
import common

class Grouper:
    ####Row is a row of the processing database and peak
    def __init__(self,row,temp,out,files,names,intensity,mztol,rttol,filter_qc,fold_blank):
        self.hash = row[3]
        self.temp = temp
        self.output_data=os.path.join(out,"datamatrix_"+intensity+"_"+self.hash+".csv")
        self.output_idx=os.path.join(out,"index_"+self.hash+".txt")
        self.files=files
        self.names=names
        self.intensity=intensity
        self.mztol=mztol
        self.rttol=rttol
        self.filter_qc = filter_qc
        self.fold_blank = fold_blank
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

class OnlineGrouper:
    def __init__(self,row,path_db,
    blocks,alignment,out,intensity,mztol,ppm,rttol,
    num_ref,alpha,ms2_mz_tol,ms2_rt_tol,fused_mgf,temp_dm_1,
                 temp_dm_2,num_workers,outfig,filter_qc=0,fold_blank=3):
        self.hash = row[3]
        self.db=path_db
        self.output_data=os.path.join(out,"datamatrix_"+row[3]+".csv")
        #self.output_idx=os.path.join(out,"index_"+self.hash+".txt")
        self.blocks=blocks
        self.alignment=alignment
        self.intensity=intensity
        self.mztol=mztol
        self.ppm = ppm
        self.alpha=alpha
        self.ms2_mz_tol=ms2_mz_tol
        self.ms2_rt_tol=ms2_rt_tol
        self.rttol=rttol
        self.num_ref=num_ref
        self.fused_mgf=fused_mgf+row[3]+".mgf"
        self.dm_1=temp_dm_1
        self.dm_2=temp_dm_2
        self.num_workers=num_workers
        self.figure=outfig+row[3]+".pdf"
        self.filter_qc = filter_qc
        self.fold_blank = fold_blank

    def need_computing(self):
        return not os.path.exists(self.output_data)

    def get_output_datamatrix(self):
        return self.output_data

    def get_fused_mgf(self):
        return self.fused_mgf

    def command_line_aligning(self):
        pscript = ct.find_rscript()
        command_line = os.path.join(pscript,"construct_peaktable_online.R")
        ####We give all the name of the grouping parameters implicated in a single file
        return " ".join(["Rscript",command_line,self.db,self.blocks,self.alignment,
                        self.output_data,self.intensity,str(self.rttol),str(self.mztol),str(self.ppm),
                        str(self.num_ref),str(self.alpha),str(self.num_workers),self.figure])

    def command_line_filtering(self):
        pscript = ct.find_rscript()
        command_line = os.path.join(pscript,"filter_datamatrix.R")
        # args < - commandArgs(trailingOnly=TRUE)
        # PATH_DB < - args[1]
        # PATH_DM < - args[2]
        # TEMP_DM < - args[3]
        # FOLD_BLANK < - as.numeric(args[4])
        # QC_FRACTION < - as.numeric(args[5])
        ####We give all the name of the grouping parameters implicated in a single file
        return " ".join(["Rscript",command_line,self.db,self.output_data,
                         self.dm_1,str(self.fold_blank),str(self.filter_qc)])

    def command_line_fusing_msms(self):
        pscript = ct.find_rscript()
        command_line = os.path.join(pscript,"fusing_msms_spectra.R")
        return " ".join(["Rscript",command_line,self.db,str(self.ms2_mz_tol),str(self.ms2_rt_tol),str(self.num_workers),
                self.fused_mgf,self.dm_1,self.dm_2])
