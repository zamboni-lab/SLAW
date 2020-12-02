import pandas as pd
import numpy as np
import scipy.spatial as sp
from math import log10,floor
import common.references as cr
import model.steps.grouping as mg
import os
import subprocess
import logging


def get_scorer(name):
    if name=="ipopeak":
        return PeakpickingScorerIPO
    if name=="slawpeak":
        return PeakpickingScorerSLAW
    if name=="combined":
        return CombinedPeakScorer
    if name=="ipoalign":
        return AlignmentScorerIPO
    if name=="cvalign":
        return reproducibleCVscorer
    if name=="expalign":
        return exponentialScorer


class ExperimentScorer:
    def __init__(self,exp):
        ####We compute the retnetion time
            self.exp = exp

    def score(self,output):
        print("Need to be implemented")

    def require_alignement(self):
        return False

class PeakpickingScorer(ExperimentScorer):
    def require_alignement(self):
        return False

    def get_weight(self):
        temp = [1.0]
        return temp

class AlignmentScorer(ExperimentScorer):
    def require_alignement(self):
        return True

    @staticmethod
    def get_weight():
        temp = [1.0,2.0]
        return temp

class PeakpickingScorerIPO(PeakpickingScorer):
    def __init__(self,exp):
        self.exp = exp
    def score(self,ppm=10,tol_rt=0.015):
        all_peaktables = self.exp.get_query("SELECT output_ms FROM processing WHERE valid=1")
        all_peaktables=[p[0] for p in all_peaktables]
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

class PeakpickingScorerSLAW(PeakpickingScorer):
    def __init__(self, exp):
        self.exp = exp

    def score(self, threshold=0.2):
        all_peaktables = self.exp.get_query("SELECT output_ms FROM processing WHERE valid=1")
        all_peaktables = [p[0] for p in all_peaktables]
        ##We change all the value with their actual path in the sampling folder
        all_peaktables = [pd.read_csv(path, header=0, sep=",") for path in all_peaktables]
        for i in range(len(all_peaktables)):
            pp = all_peaktables[i]
            pp["sample"] = [i] * pp.shape[0]

        ###We concatenate the data

        # 5ppm = 0.02 peaks peakwidth =
        btab = pd.concat(all_peaktables)
        if len(btab.rt)==0:
            return -1.0
        norm_rt = max(btab.rt) / 20
        btab.mz = btab.mz / 0.005
        btab.rt = btab.rt / norm_rt

        ####We get the n
        sum_btab = btab[["mz", "rt"]]
        ktree = sp.KDTree(sum_btab)
        neighbours = ktree.query(sum_btab, k=len(all_peaktables))[1]
        ###We check if all the sample are different
        sel_idx = np.apply_along_axis(lambda x, idx: len(np.unique(idx.iloc[x])), 1, neighbours,
                                      idx=btab['sample'])

        ###To avoid 0 cases just select the one with
        vmax = max(sel_idx)
        sel_idx = sel_idx==vmax
        neighbours = neighbours[sel_idx]
        if np.sum(neighbours)==0:
            print("No feature found in every sample.")
            return -1

        ####THe thresold is the number of feature with very few

        peakwidths = np.apply_along_axis(lambda x, idx: idx.iloc[x], 0, neighbours, idx=btab["peakwidth"])
        cv_peakwidth = np.apply_along_axis(lambda x: np.std(x) / np.mean(x), 1, peakwidths)
        num_correct = np.sum(cv_peakwidth < threshold)
        if num_correct==0:
            # print("No feature with reproducible peakwidth found.")
            return -1
        return ((num_correct ** 2) / len(btab))

class CombinedPeakScorer(PeakpickingScorer):
    def __init__(self,exp):
        self.sc1 = PeakpickingScorerSLAW(exp)
        self.sc2 = PeakpickingScorerIPO(exp)

    def get_weight(self):
        return (1,1)

    def score(self,output):
        v1 = self.sc1.score()
        v2 = self.sc2.score()
        if v1==-1 or v2==-1:
            return (-1,-1)
        logging.debug("PEAKPICKING SCORE, INTEGRATION:"+str(v1)+" ISOTOPES:"+str(v2))
        return (v1,v2)

class AlignmentScorerIPO(AlignmentScorer):
    def score(self,output):
        ###We build a data matrix with all tthe retention time
        fake_pp = ("","","","temphash")
        dir_blocks = self.exp.output.getDir(cr.TEMP["GROUPING"]["BLOCKS"])
        dir_alignment = self.exp.output.getFile(cr.TEMP["GROUPING"]["ALIGNMENT"])
        dir_datamatrix = self.exp.output.getDir(cr.OUT["DATAMATRIX"])
        path_fig = self.exp.output.getFile(cr.OUT["FIGURES"]["RT_DEV"])
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

        tdata = pd.read_csv(dm_path,sep="\t",header=0)
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



def modif_exp(x,alpha=6):
    return (np.exp(alpha*x)-1)/(np.exp(alpha)-1)

class exponentialScorer(AlignmentScorer):
    def score(self, output, alpha = 10):
        ###We build a data matrix with all tthe retention time
        fake_pp = ("", "", "", "temphash")
        dir_blocks = self.exp.output.getDir(cr.TEMP["GROUPING"]["BLOCKS"])
        # dir_blocks = dir_blocks.replace("/output",output,1)
        dir_alignment = self.exp.output.getFile(cr.TEMP["GROUPING"]["ALIGNMENT"])
        # dir_alignment = dir_alignment.replace("/output", output, 1)
        dir_datamatrix = self.exp.output.getDir(cr.OUT["DATAMATRIX"])
        # dir_datamatrix = dir_datamatrix.replace("/output", output, 1)
        path_fig = self.exp.output.getFile(cr.OUT["FIGURES"]["RT_DEV"])
        # path_fig = path_fig.replace("/output", output, 1)
        hdat = "rt_cor"
        ###We copy the database for examintation

        ppg = mg.OnlineGrouper(fake_pp, self.exp.db, dir_blocks, dir_alignment,
                               dir_datamatrix, hdat, 0.01, 15, 0.01, 150, 0.01,
                               0.05, 0.05, "NONE",
                               "NONE", "NONE", 1, path_fig)
        cli = ppg.command_line_aligning()
        dm_path = ppg.get_output_datamatrix()

        ###The data matrix  is creconstructed as the alignment is done already.
        subprocess.call(cli, shell=True, timeout=900)
        if not os.path.exists(dm_path):
            return -1.0, -1.0
        tdata = pd.read_csv(dm_path, header=0,sep="\t")
        cnames = [cc for cc in tdata.columns if cc.startswith(hdat)]
        rt_tab = tdata[cnames]
        ###In vevery case we consider that a deviation superior to 0.3% of the maximum
        val = rt_tab.apply(lambda x: np.nanmean(np.absolute(x - np.nanmedian(x))), axis=1)
        # rt_tab.apply(lambda x: x, axis=1)
        ###rt is equl to htis one
        ARTS = np.nanmean(val)
        RCS = np.log10(1 / ARTS)
        summed_probas = np.sum(modif_exp(tdata.num_detection / len(cnames),alpha=alpha))
        logging.debug("ALIGNMENT SCORE, DRT:"+str(RCS)+" REPRODUCIBILITY:"+str(summed_probas))
        ###We just have to comnpare the column
        return RCS, summed_probas

class reproducibleCVscorer(AlignmentScorer):

    def score(self,output,num_skipped=1):
        path_dm = self.exp.get_datamatrix()[0][0]
        # path_dm = path_dm.replace("/output", output, 1)
        if not os.path.exists(path_dm):
            return -1.0
        data = pd.read_csv(path_dm,sep="\t",header=0)

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
