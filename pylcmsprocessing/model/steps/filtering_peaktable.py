
import pandas as pd

# Reference term
ABSOLUTE_PREFIX = "absolute"
RELATIVE_PREFIX = "relative"
METRICS = [ABSOLUTE_PREFIX+"_height",
    RELATIVE_PREFIX+"_height",
    ABSOLUTE_PREFIX+"_intensity",
    RELATIVE_PREFIX+"_intensity",
    "SN",
    "peakwidth",
    "mz",
    "rt"]

UNORDERED_OPERATORS = [
    "greater",
    "lower",
    "equal",
]
ORDERED_OPERATORS = [
    "top",
    "bottom"
]
OPERATORS = UNORDERED_OPERATORS+ORDERED_OPERATORS

#We define all the filter as function (Could be better but th epython name are not clear engouh)
def filter_greater(serie,value):
    return serie>=value

def filter_lower(serie,value):
    return serie<=value

def filter_equal(serie,value):
    return serie==value

def filter_top(serie,value):
    filter_idx = pd.Series([False]*len(serie))
    idx = np.argsort(serie)
    sel_idx = idx[-int(value):]
    filter_idx[sel_idx] = True
    return filter_idx

def filter_bottom(serie,value):
    filter_idx = pd.Series([False]*len(serie))
    idx = np.argsort(serie)
    sel_idx = idx[:int(value)]
    filter_idx[sel_idx] = True
    return filter_idx

#These functions parse the filter string with metric,operator,value field
def parse_one_filter(fstring):
    sargs = fstring.split(" ")
    if len(sargs)>3:
        raise ValueError("A 'filter' string should be composed of three elements, a Metric, an Operator and a Value")
    metric,operator,value = sargs
    if metric not in METRICS:
        raise ValueError("Invalid metric {}, authorized Metrics are {}".format(metric,",".join([x for x in METRICS])))
    if operator not in OPERATORS :
        raise ValueError("Invalid operator {}, authorized Operators are {}".format(operator,",".join([x for x in OPERATORS])))
    try:
        value = float(value)
    except ValueError:
        raise(ValueError("Invalid threshold value {}, it should be convertible to a float.".format(value)))
    return metric,operator,value

def parse_filters(fstring):
    """This functin take a filtering step and remov ethe rest of the data"""
    sargs = fstring.split("&")
    return [parse_one_filter(ff.strip()) for ff in sargs]

def apply_filter(peaktable,metric,operator,value):
    prefix = "absolute"
    if metric.startswith(ABSOLUTE_PREFIX) or metric.startswith(RELATIVE_PREFIX):
        prefix,metric = metric.split("_")
    raw_values = peaktable[metric]
    if prefix=="relative":
        raw_values = raw_values/max(raw_values)

    func_filter = "filter_"+operator
    sel_elems = eval(func_filter+"(raw_values,value)")
    return peaktable[sel_elems.tolist()]



class PeaktableFilter:
    def __init__(self,fstring):
        self.filters = parse_filters(fstring)

    def filter_peaktable(self,peaktable):
        for metric,operator,value in self.filters:
            peaktable = apply_filter(peaktable,metric,operator,value)
        return peaktable





if __name__=="__main__":
    import pandas as pd
    import numpy as np
    PATH_DT = "U:/processing/out/lcms_red2_redo2/CENTWAVE/peaktables/KO1_BN051189-91_MN-agar_7_DCM_Pos_R1.csv"
    peaktable = pd.read_csv(PATH_DT)
    TEST_STRING = "absolute_height greater 10000 & relative_height bottom 100"
    pf = PeaktableFilter(TEST_STRING)
    res_peaktable = pf.filter_peaktable(peaktable)
    res_peaktable

