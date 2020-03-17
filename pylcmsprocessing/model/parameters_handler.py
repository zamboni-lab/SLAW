import yaml
import numpy as np
import common.references as cr

def recur_node(x,path=None,all_args=None):
    if all_args is None:
        all_args = []
    if path is None:
        path = []
    if isinstance(x,dict) and "value" in x:
        return all_args.append(tuple(path))
    elif isinstance(x,dict):
        for k in x:
            recur_node(x[k],path+[k],all_args)
    else:
        pass

def get_all_parameters(ex):
    all_args=[]
    recur_node(ex,all_args=all_args)
    return all_args

def get_val(raw,path):
    if len(path)==0:
        return raw
    k = path[0]
    return get_val(raw[k],path[1:])

def set_val(raw,path,value,field="value"):
    if len(path)==0:
        raw[field]=value
        return
    k = path[0]
    set_val(raw[k],path[1:],value)

class ParametersFileHandler:
    SEP="__"
    def __init__(self,path=None):
        self.path=path
        with open(self.path, 'r') as stream:
            self.dic = yaml.safe_load(stream)
        self.param_path = get_all_parameters(self.dic)

    def __getitem__(self, item):
        if isinstance(item,str):
            item = item.split(ParametersFileHandler.SEP)
        return get_val(self.dic,item)

    def __setitem__(self, key, value):
        if isinstance(key,str):
            key = key.split(ParametersFileHandler.SEP)
        set_val(self.dic,key,value)

    def get_optimizable_parameters(self,string=True):
        to_optimize = {}
        for path in self.param_path:
            val = self[path]
            if "range" in val:
                ppath = path
                if string:
                    ppath = ParametersFileHandler.SEP.join(ppath)
                else :
                    ppath = tuple(path)
                to_optimize[ppath]=val["range"]
        return to_optimize

    def get_parameters(self):
        return self.param_path

    def get_parameters(self,string=True):
        if string:
            return [ParametersFileHandler.SEP.join(pp) for pp in self.param_path]
        return self.param_path

    def write_parameters(self,path):
        with open(path, 'w') as outfile:
            yaml.dump(path, outfile, default_flow_style=False)

if __name__=="__main__":
    PATH_PARAMS = "C:/Users/dalexis/Documents/dev/lcmsprocessing_docker/pylcmsprocessing/data/parameters_set.yaml"
    pfh = ParametersFileHandler(PATH_PARAMS)
    pfh.get_optimizable_parameters()
    all_par = pfh.get_parameters(string=False)
    all_par_str = pfh.get_parameters()

    print("Old value  of ",all_par[2],"is",pfh[all_par[2]]["value"])
    pfh[all_par[2]] = 10
    print("New value  of ",all_par[2],"is",pfh[all_par[2]]["value"])
    pfh[all_par_str[2]] = 30
    print("New value with str  of ",all_par[2],"is",pfh[all_par[2]]["value"])

