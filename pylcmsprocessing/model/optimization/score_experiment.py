import pandas as pd
import numpy as np
import scipy.spatial as sp
from math import log10,floor
import common.references as cr
import common.tools as ct
import model.steps.grouping as mg
import os
import subprocess


def get_scorer(name):
    if name=="ipopeak":
        return PeakpickingScorerIPO
    if name=="ipoalign":
        return AlignmentScorerIPO
    if name=="groupcv":
        return reproducibleCVscorer
    return reproducibleCVscorer


class ExperimentScorer:
    def __init__(self,exp):
        ####We compute the retnetion time
            self.exp = exp

    def score(self,output):
        print("Need to be implemented")


class PeakpickingScorerIPO(ExperimentScorer):
    def __init__(self,exp):
        self.exp = exp
    def score(self,output,ppm=10,tol_rt=0.015):
        all_peaktables = self.exp.get_query("SELECT output_ms FROM processing WHERE valid=1")

        ##We adapt the output
        all_peaktables=[p[0] for p in all_peaktables]
        # all_peaktables=[p[0].replace("/output",output,1) for p in all_peaktables]
        ##We change all the value with their actual path in the sampling folder
        val = 0
        for path in all_peaktables:
            if not os.path.exists(path):
                return -1.0
            vdata = pd.read_csv(path,header=0,sep=",")
            MASS_C13 = 1.003355
            MASS_CH3 = 15.023475
            MASS_CH2 = 14.01565
            RATIO_C13 = 1.109

            ##Counter of term used by IPO
            RP = 0
            LIP = 0
            all_peaks = vdata.shape[0]

            ###Lok for all the possible C13 and sec
            cnames = "intensity"
            int_tab = vdata[cnames]
            int_vec = int_tab.apply(func=lambda x: np.nanmean(x))
            ##We estimate the noise level as iin the IPO paper.
            noise_value = np.quantile(int_vec,0.03)

            order = sorted(range(len(int_vec)), key=lambda k: int_vec[k])
            used = [False]*len(order)

            ###Get the order of pthe peak by intensity
            spatial_mz = vdata[["mz"]]
            spatial_rt = vdata[["rt"]]
            mz_tree = sp.KDTree(spatial_mz)
            rt_tree = sp.KDTree(spatial_rt)
            for io in order:
                # print(io)
                mz = spatial_mz.mz[order[io]]
                int = int_vec[order[io]]

                ##We first check fi th epeak is a LIP
                maxC13 = floor((mz-2*MASS_CH3)/MASS_CH2)+2
                maxInt = int*maxC13*RATIO_C13
                if maxInt<noise_value:
                    LIP+=1
                    continue
                mz_C13 = mz+MASS_C13
                mz_epsilon = spatial_mz.mz[order[io]]*ppm/1e6
                idx_mz = mz_tree.query_ball_point((mz_C13,),mz_epsilon,p=1)
                smz = len(idx_mz)
                idx_mz = [f for f in idx_mz if not used[f] and order[f]>io]
                if len(idx_mz)==0:
                    continue
                rt = spatial_rt.rt[order[io]]
                idx_rt = rt_tree.query_ball_point((rt,),tol_rt,p=1)
                idx_rt = [f for f in idx_rt if not used[f] and order[f]>io and f in idx_mz]
                srt = len(idx_rt)
                if len(idx_rt)==0:
                    continue
                ###We check if there is an intersection

                idx_int = [f for f in idx_rt if int_vec[f]<maxInt]
                if len(idx_int)>0:
                    RP += 1
                    used[idx_int[0]]=True
                    used[order[io]]=True
            val += RP**2/(all_peaks-LIP)
        return val

class AlignmentScorerIPO(ExperimentScorer):
    def score(self,output):
        ###We build a data matrix with all tthe retention time
        fake_pp = ("","","","temphash")
        dir_blocks = self.exp.output.getDir(cr.TEMP["GROUPING"]["BLOCKS"])
        # dir_blocks = dir_blocks.replace("/output",output,1)
        dir_alignment = self.exp.output.getFile(cr.TEMP["GROUPING"]["ALIGNMENT"])
        # dir_alignment = dir_alignment.replace("/output", output, 1)
        dir_datamatrix = self.exp.output.getDir(cr.OUT["DATAMATRIX"])
        # dir_datamatrix = dir_datamatrix.replace("/output", output, 1)
        path_fig = self.exp.output.getFile(cr.OUT["FIGURES"]["RT_DEV"])
        # path_fig = path_fig.replace("/output", output, 1)
        hdat = "rt_cor"
        ppg = mg.OnlineGrouper(fake_pp, self.exp.db, dir_blocks, dir_alignment,
                               dir_datamatrix, hdat, 0.01, 15, 0.01, 150, 0.01,
                              0.05, 0.05, "NONE",
                               "NONE", "NONE", 1, path_fig)
        cli = ppg.command_line_aligning()
        dm_path = ppg.get_output_datamatrix()

        ###The data matrix  is creconstructed as the alignment is done already.
        subprocess.call(cli, shell=True, timeout=900)
        if not os.path.exists(dm_path):
            return -1.0,-1.0

        tdata = pd.read_csv(dm_path, header=0)
        print(tdata)
        cnames = [cc for cc in tdata.columns if cc.startswith(hdat)]
        num_sample = len(cnames)

        rt_tab = tdata[cnames]
        val = rt_tab.apply(lambda x: np.absolute(np.nanmean(x-np.nanmedian(x))), axis=1)
        ARTS = np.nanmean(val)
        RCS = 1/ARTS

        thresh_val = floor(len(cnames)*0.8)
        rps = np.sum((tdata.num_detection > thresh_val) & (tdata.total_detection == tdata.num_detection))
        total_peaks = tdata.shape[0]
        GS = (rps**2)/total_peaks
        ###We just have to comnpare the column

        return RCS,GS


class reproducibleCVscorer(ExperimentScorer):

    def score(self,output,num_skipped=1):
        path_dm = self.exp.get_datamatrix()[0][0]
        # path_dm = path_dm.replace("/output", output, 1)
        if not os.path.exists(path_dm):
            return -1.0
        data = pd.read_csv(path_dm, header=0)

        cnames = [cc for cc in data.columns if cc.startswith("int")]
        num_sample = len(cnames)

        int_tab = data[cnames]
        num_detect = num_sample-int_tab.isnull().sum(axis=1)
        CV = int_tab.apply(lambda x: np.nanstd(x)/np.nanmean(x), axis=1)

        weighted_CV = (num_detect*CV).sum()/num_detect.sum()
        ###We return the sum
        n_reproducible = sum([nn for nn in num_detect if nn>=(num_sample-num_skipped)])
        score = (0.5-weighted_CV)*log10(n_reproducible)
        return score
