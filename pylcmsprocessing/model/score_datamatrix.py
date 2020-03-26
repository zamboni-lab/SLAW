import pandas as pd
import numpy as np
from math import log10
import os


class scorerDataMatrix:
    def __init__(self,path):
        self.path = path
        if os.path.exists(path):
            self.data = pd.read_csv(path,header=0)
        else:
            self.data=None

    ####Score a data mnatrix
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

