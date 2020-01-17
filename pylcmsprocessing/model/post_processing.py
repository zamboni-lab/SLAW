import sqlite3
import subprocess
import os
import glob
import common.references as cr
import common.tools as ct
import model.output_handler as oh
import model.parameters as params
import random
import model.inputbuilder as ib
import model.parallel_runner as pr
import model.peakpicking as mp
import model.grouping as mg
import model.evaluating as me
import model.comparing_evaluation as mce
import model.annotating_adducts_isotopes as mai
import pandas as pd
import shutil


class PostProcessing:

    def __init__(self,path_db,path_targets,path_fig_target,path_fig_summary,
    path_output_hdf5,num_workers,raw_files,mztol,rttol):
        self.db = path_db
        self.target=path_targets
        self.path_hdf5 = path_output_hdf5
        self.path_out_summary=path_fig_summary
        self.path_out_target=path_fig_target
        self.raw_files=raw_files
        self.mztol=mztol
        self.rttol=rttol
        self.num_workers=num_workers

    def command_line(self):
        pscript = ct.find_rscript()
        command_line = os.path.join(pscript,"post_processing.R")
        ####We give all the name of the grouping parameters implicated in a single file
        return " ".join(["Rscript",command_line,self.target,self.raw_files,self.db,
                        self.path_hdf5,self.path_out_target,self.path_out_summary,
                        str(self.mztol),str(self.rttol),str(self.num_workers)])



###Launche the various figures evetually.
