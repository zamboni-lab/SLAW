import common.tools as ct
import os

class OnlineGrouper:
    def __init__(self,row,path_db,
    blocks,alignment,out,intensity,mztol,ppm,rttol,
    num_ref,alpha,ms2_mz_tol,ms2_rt_tol,fused_mgf,temp_dm_1,
                 temp_dm_2,path_hdf5,num_workers,outfig,filter_qc=0,fold_blank=3):
        self.hash = row[3]
        self.db=path_db
        self.output_data=os.path.join(out,"data_filled_"+row[3]+".csv")
        self.output_data_non_filled = os.path.join(out,"data_prefill_"+row[3]+".csv")
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
        self.hdf5 = path_hdf5
        self.num_workers=max(int(num_workers/2),1)
        self.figure=outfig+row[3]+".pdf"
        self.filter_qc = filter_qc
        self.fold_blank = fold_blank

    def need_computing(self):
        return not os.path.exists(self.output_data)

    def get_output_datamatrix(self):
        return self.output_data

    def get_non_filleddatamatrix(self):
        """This function is only used to output the non filled data matrix"""
        return self.output_data_non_filled

    def get_fused_mgf(self):
        return self.fused_mgf

    def command_line_aligning(self):
        pscript = ct.find_rscript()
        command_line = os.path.join(pscript,"construct_peaktable_online.R")
        ####We give all the name of the grouping parameters implicated in a single file
        return " ".join([os.environ["RscriptString"]," ",command_line,'"'+self.db+'" "'+self.blocks+'" "'+self.alignment+'"',
                        '"'+self.output_data_non_filled+'"',self.intensity,str(self.rttol),str(self.mztol),str(self.ppm),
                        str(self.num_ref),str(self.alpha),str(self.num_workers),'"'+self.figure+'"'])

    def command_line_filtering(self):
        pscript = ct.find_rscript()
        command_line = os.path.join(pscript,"filter_datamatrix.R")
        # args < - commandArgs(trailingOnly=TRUE)

        ####We give all the name of the grouping parameters implicated in a single file
        return " ".join([os.environ["RscriptString"]," ",command_line,'"'+self.db+'"','"'+self.output_data_non_filled+'"',
                         '"'+self.dm_1+'"',str(self.fold_blank),str(self.filter_qc)])

    def command_line_fusing_msms(self):
        pscript = ct.find_rscript()
        command_line = os.path.join(pscript,"fusing_msms_spectra.R")
        return " ".join([os.environ["RscriptString"],command_line,'"'+self.db+'"',str(self.ms2_mz_tol),str(self.ms2_rt_tol),str(self.num_workers),
                '"'+self.fused_mgf+'"','"'+self.dm_1+'"','"'+self.dm_2+'"','"'+self.hdf5+'"'])
