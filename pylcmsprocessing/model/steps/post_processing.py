import os
import common.tools as ct



class PostProcessing:

    def __init__(self,path_db,path_targets,path_fig_target,path_fig_summary,
    path_tab_rt,path_tab_int,path_output_hdf5,num_workers,raw_files,mztol,rttol):
        self.db = path_db
        self.target=path_targets
        self.path_hdf5 = path_output_hdf5
        self.path_out_summary=path_fig_summary
        self.path_out_target=path_fig_target
        self.path_rt_table=path_tab_rt
        self.path_int_table=path_tab_int
        self.raw_files=raw_files
        self.mztol=mztol
        self.rttol=rttol
        self.num_workers=num_workers

    def command_line(self):
        pscript = ct.find_rscript()
        command_line = os.path.join(pscript,"post_processing.R")
        ####We give all the name of the grouping parameters implicated in a single file
        return " ".join([os.environ["RscriptString"]," ",command_line,self.target,self.raw_files,self.db,
                        self.path_hdf5,self.path_out_target,self.path_out_summary,
                        self.path_rt_table, self.path_int_table,
                        str(self.mztol),str(self.rttol),str(self.num_workers)])



###Launche the various figures evetually.
