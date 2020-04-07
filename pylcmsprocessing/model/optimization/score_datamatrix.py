import pandas as pd
import numpy as np
import scipy.spatial as sp
from math import log10,floor
import os


def get_scorer(name):
    if name=="cv":
        return reproducibleCVscorer
    if name=="ipo":
        return IPOscorer
    return IPOscorer


class scorerDataMatrix:
    def __init__(self,path):
        self.path = path
        if os.path.exists(path):
            self.data = pd.read_csv(path,header=0)
        else:
            self.data=None


    ####Score a data mnatrix
    def score_datamatrix(self):
        print("Method need to be implmented")

class reproducibleCVscorer(scorerDataMatrix):

    def score_datamatrix(self,num_skipped=1):
        if self.data is None:
            return -1.0

        cnames = [cc for cc in self.data.columns if cc.startswith("int")]
        num_sample = len(cnames)

        int_tab = self.data[cnames]
        num_detect = num_sample-int_tab.isnull().sum(axis=1)
        CV = int_tab.apply(lambda x: np.nanstd(x)/np.nanmean(x), axis=1)

        weighted_CV = (num_detect*CV).sum()/num_detect.sum()
        ###We return the sum
        n_reproducible = sum([nn for nn in num_detect if nn>=(num_sample-num_skipped)])
        score = (0.5-weighted_CV)*log10(n_reproducible)
        return score

class IPOscorer(scorerDataMatrix):
    def score_datamatrix(self,ppm=10,tol_rt=0.015):
        MASS_C13 = 1.003355
        MASS_CH3 = 15.023475
        MASS_CH2 = 14.01565
        RATIO_C13 = 1.109

        ##Counter of term used by IPO
        RP = 0
        LIP = 0
        all_peaks = self.data.shape[0]

        ###Lok for all the possible C13 and sec
        cnames = [cc for cc in self.data.columns if cc.startswith("int")]
        int_tab = self.data[cnames]
        int_vec = int_tab.apply(lambda x: np.nanmean(x), axis=1)
        ##We estimate the noise level as iin the IPO paper.
        noise_value = np.quantile(int_vec,0.03)

        order = sorted(range(len(int_vec)), key=lambda k: int_vec[k])
        used = [False]*len(order)

        ###Get the order of pthe peak by intensity
        spatial_mz = self.data[["mz"]]
        spatial_rt = self.data[["rt"]]
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

        return RP**2/(all_peaks-LIP)

if __name__=="__main__":
    PATH_TEST = "C:/Users/dalexis/Documents/dev/docker_sing_output/datamatrix_2e6d0ad4448475fcc09d1644dcc36740.csv"
    scorer = get_scorer("cv")
    scm = scorer(PATH_TEST)
    scm.score_datamatrix()

