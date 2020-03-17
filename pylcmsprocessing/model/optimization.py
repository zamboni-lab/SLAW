###This take as input the paramtters ot ptimize and return the rest of the data
import hashlib
import os
import shutil
import yaml
from random import sample
import subprocess
import numpy as np
import pandas as pd

import model.experiment as me
import common.references as cr
import model.score_datamatrix as ms
import model.optimizer as mlp
import model.parameters_handler as ph
import common.tools as ct
from model.UI import UI

###Function used to optimize the paramters
def create_temp_directory(path_exp,params_archive,*argv):
    '''
    :param argv: The arguments ot be hashed
    :return: A triplet containing the hash of the paramters, the experiemnt directory and the paramters name
    '''
    ###We create a hash based on the joined paramters
    hash_val = "_".join([str(pp) for pp in argv])
    hash_val = int(hashlib.sha1(hash_val.encode()).hexdigest(), 16) % (10 ** 8)
    temp_dir = os.path.join(path_exp,str(hash_val))
    os.makedirs(temp_dir, exist_ok=True)
    stored_param = os.path.join(params_archive,"param_"+str(hash_val)+".yaml")
    temp_db = os.path.join(temp_dir,"processing_db.sqlite")
    temp_save_db = os.path.join(temp_dir,"save_processing_db.sqlite")
    stored_xml = os.path.join(params_archive,"xml_"+str(hash_val)+".xml")
    return hash_val,temp_dir,temp_db,temp_save_db,stored_param,stored_xml


def convert_val(x):
    if type(x) is np.ndarray:
        if isinstance(x[0], np.int):
            return int(x.item())
        if isinstance(x[0], np.float64):
            return float(x.item())
        return x.item()
    return x



# def peak_picking_alignment_scoring(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,fixed_params):
def peak_picking_alignment_scoring(peakpicking__noise_level_ms1,peakpicking__noise_level_ms2,peakpicking__traces_construction__ppm,
                                   peakpicking__traces_construction__dmz,peakpicking__traces_construction__min_scan,
                                   peakpicking__peaks_deconvolution__SN,peakpicking__peaks_deconvolution__peak_width_min,
                                   peakpicking__peaks_deconvolution__rt_wavelet_min,peakpicking__peaks_deconvolution__rt_wavelet_max,
                                   peakpicking__peaks_deconvolution__coefficient_area_threshold,grouping__ppm,grouping__drt,grouping__dmz,
                                   grouping__alpha,grouping__num_references,fixed_params):
    args_refs = locals()
    del args_refs["fixed_params"]

    def is_processed(stored_param, summary_table):
        if not os.path.isfile(summary_table):
            return False, 0
        rr = pd.read_csv(summary_table)
        sub_row = rr[[stored_param in xx for xx in rr.path_parameters]]
        if sub_row.shape[0]==0:
            return False, 0
        return True, sub_row.iloc[0,1]

    def convert_val(x):
        if type(x) is np.ndarray:
            if isinstance(x[0], np.int):
                return int(x.item())
            if isinstance(x[0], np.float64):
                return float(x.item())
            return x.item()
        return x

    # p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14=params
    DIR_TEMP,STORAGE_YAML,SUMMARY_YAML,num_cpus,polarity,path_samples,initial_yaml=fixed_params
    call_list = [DIR_TEMP,STORAGE_YAML]+list(args_refs.values())
    ###we create the directory andget the path which will be used in the data.
    hash_val,OUTPUT_DIR,PATH_DB,PATH_SAVE_DB,stored_param,stored_xml = create_temp_directory(*call_list)

    processed = is_processed(stored_param,SUMMARY_YAML)
    if processed[0]:
        return float(processed[1])

    ###We create a temporary directory to put the data in
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    LOG_PATH = os.path.join(OUTPUT_DIR,"log.txt")

    ####We load the orginial yaml file
    parameters = ph.ParametersFileHandler(initial_yaml)
    for k in args_refs.keys():
        parameters[k] = args_refs[k]
    ##We update the parameters to reflect the rest of the data
    # raw_yaml["peakpicking"]["noise_level_ms1"]["value"]=float(convert_val(peakpicking_noise_level_ms1))
    # raw_yaml["peakpicking"]["noise_level_ms2"]["value"]=float(convert_val(peakpicking_noise_level_ms2))
    # raw_yaml["peakpicking"]["traces_construction"]["ppm"]["value"]=float(convert_val(peakpicking_traces_construction_ppm))
    # raw_yaml["peakpicking"]["traces_construction"]["dmz"]["value"]=float(convert_val(peakpicking_traces_construction_dmz))
    # raw_yaml["peakpicking"]["traces_construction"]["min_scan"]["value"]=int(convert_val(peakpicking_traces_construction_min_scan))
    # raw_yaml["peakpicking"]["peaks_deconvolution"]["SN"]["value"]=float(convert_val(peakpicking_peaks_deconvolution_SN))
    # raw_yaml["peakpicking"]["peaks_deconvolution"]["peak_width_min"]["value"]=float(convert_val(peakpicking_peaks_deconvolution_peak_width_min))
    # raw_yaml["peakpicking"]["peaks_deconvolution"]["rt_wavelet_min"]["value"]=float(convert_val(peakpicking_peaks_deconvolution_rt_wavelet_min))
    # raw_yaml["peakpicking"]["peaks_deconvolution"]["rt_wavelet_max"]["value"]=float(convert_val(peakpicking_peaks_deconvolution_rt_wavelet_max))
    # raw_yaml["peakpicking"]["peaks_deconvolution"]["coefficient_area_threshold"]["value"]=float(convert_val(peakpicking_peaks_deconvolution_coefficient_area_threshold))
    ###We dump the yaml
    parameters.write_parameters(stored_param)

    PATH_DB = os.path.join(OUTPUT_DIR,"temp_processing"+str(hash_val)+"_db.sqlite")
    PATH_SAVE_DB = os.path.join(OUTPUT_DIR,"processing_db.sqlite")
    PATH_XML = os.path.join(OUTPUT_DIR,"temp_par"+str(hash_val)+".xml")

    exp = me.Experiment(PATH_DB,save_db = PATH_SAVE_DB,reset=False)
    ###We create the UI
    vui = UI(OUTPUT_DIR, path_samples, polarity=polarity, mass_spec="Exactive", num_workers=num_cpus, path_yaml=stored_param)
    vui.generate_yaml_files(cr.DATA["YAML"])
    ###The output database is generated.
    vui.generate_MZmine_XML(path_xml=PATH_XML)
    exp.initialise_database(num_cpus, OUTPUT_DIR, polarity, path_samples, ["ADAP"], 1)
    exp.building_inputs_single_processing(PATH_XML)
    try:
        exp.run("/MZmine-2.52-Linux", int(num_cpus), log=LOG_PATH)
        exp.correct_conversion()
        exp.post_processing_peakpicking()
    except Exception:
        subprocess.call("rm -r " + OUTPUT_DIR)
        # shutil.rmtree(OUTPUT_DIR)
        return -1

    exp.group_online(intensity="intensity",
                     ppm=float(convert_val(grouping__ppm)),
                     mztol=float(convert_val(grouping__dmz)),
                     rttol=float(convert_val(grouping__drt)),
                     n_ref=int(convert_val(grouping__num_references)),
                     alpha=float(convert_val(grouping__alpha)),
                     log=LOG_PATH)
    ###We load all the datamtrices.
    datamatrices = exp.get_datamatrix()
    path_datamatrix = datamatrices[0][0]
    try:
        msd = ms.scorerDataMatrix(path_datamatrix)
        ###We authorize 1 jump
        vscore = msd.score_datamatrix(1)
    except ValueError:
        subprocess.call("rm -r " + OUTPUT_DIR)
        return -1

    ###We write the score in the table
    to_join = list(args_refs.values())
    to_join = [str(pp) for pp in to_join]
    to_write = stored_param+","+str(vscore)+",".join(to_join)+"\n"

    ###We wrtie the header if it is the first one.
    if not os.path.isfile(SUMMARY_YAML):
        ##We write the path name.
        to_join=["path_parameters","score"]+list(args_refs.keys())
        joined_summary = ",".join(to_join)
        with open(SUMMARY_YAML, "w") as ff:
            ff.write(joined_summary)

    with open(SUMMARY_YAML,"a") as ff:
        ff.write(to_write)
    ##After each evaluation we remove the the directory.
    shutil.rmtree(OUTPUT_DIR)
    return vscore

class ParametersOptimizer:
    def __init__(self,exp,output,nrounds=10,input_par=None):
        if input_par is None:
            input_par = cr.DATA["YAML"]
        self.input_par = input_par
        self.output = output
        self.nrounds=nrounds
        # At the initialisation a temp directory for the output is created.
        self.path_exp=os.path.join(self.output,"temp_exp")
        os.makedirs(self.path_exp, exist_ok=True)
        self.params_archive=os.path.join(self.output,"parameters")
        os.makedirs(self.params_archive, exist_ok=True)
        self.path_samples = os.path.join(self.output,"mzML")
        os.makedirs(self.path_samples, exist_ok=True)
        self.path_summary = os.path.join(self.output,"summary_par.csv")
        self.temp_polarity = os.path.join(self.output,"polarity.csv")
        self.polarity = exp.get_polarity()
        self.num_workers = exp.get_workers()
        self.samples = exp.get_query("SELECT path FROM samples")
        self.path_db = exp.db
        self.temp_yaml = os.path.join(self.output,"temp_yaml.yaml")

    # Select the smaples ot optimize evetually.
    def select_samples(self,num_files=None):
        if num_files is None:
            # As many as worker to speed up the peakpicking process
            num_files = max(self.num_workers-1,1)
        if len(self.samples) < num_files:
            num_files = len(self.samples)
        sel_samples = sample(self.samples,num_files)
        for ss in sel_samples:
            shutil.copy(ss[0],self.path_samples)

    def determine_initial_parameters(self):
        pscript = os.path.join(ct.find_rscript(),"initial_parameters.R")

        # we tun the optimization ins a signle thread.
        path_db = self.path_db
        output_par = self.temp_yaml
        initial_par = self.input_par
        num_cores = self.num_workers
        cli = " ".join(["Rscript",pscript,path_db,initial_par,output_par,str(num_cores)])
        self.input_par=output_par
        ###We read the polarity directly
        subprocess.call(cli, shell=True, env=os.environ.copy())


    def do_optimize_parameters(self,optimizer="LIPO",**kwargs):
        ###We read the standar doptimization interval eventually
        optim_func = mlp.get_optimizers(optimizer)

        pfh = ph.ParametersFileHandler(self.temp_yaml)
        to_optimize = pfh.get_optimizable_parameters(string=True)
        all_parameters = pfh.get_parameters()
        lb = [x[0] for x in to_optimize.values()]
        ub = [x[1] for x in to_optimize.values()]

        ###if the yaml has not been optimize we create the first value
        if not os.path.isfile(self.temp_yaml):
            self.temp_yaml = self.input_par

        # Paramters which are always fixed.
        fixed_params = self.path_exp, self.params_archive, self.path_summary, self.num_workers, self.polarity, self.path_samples, self.temp_yaml
        dic_fixed = {"fixed_params":fixed_params}

        ###We fix the paramter which are not yet fixed
        for ll in all_parameters:
            if ll not in to_optimize:
                dic_fixed[ll]=pfh[ll]

        # we optimize the peakpicking
        voptim = optim_func(lb,ub,peak_picking_alignment_scoring,fixed_arguments=dic_fixed,**kwargs)
                         #LIPO algoirhtm  ,max_call=self.nrounds,initial_points=3,)
        final_dic = dic_fixed
        for io in range(len(voptim)):
            final_dic[to_optimize[io]] = convert_val(voptim[io])
        # The best paramters is stroed
        self.best_par = final_dic

    def export_best_parameters(self,outpath):
        pfh = ph.ParametersFileHandler(self.input_par)
        for k in self.best_par:
            pfh[k] = convert_val(self.best_par[k])
        pfh.write_parameters(outpath)


    def optimize_parameters(self,output_par,optimizer="LIPO",**kwargs):
        ##Supplementary arugment for LIPO
        # max_call=self.nrounds,initial_points=3
        ##Supplementary for random sampling
        # num_points=1000,num_cores = 1)
        self.determine_initial_parameters()
        self.select_samples()
        print("Finished initial parameters estimation")
        print("Optimizing remaining parameters")
        self.do_optimize_parameters(optimizer=optimizer,**kwargs)
        print("Finished  optimization")
        self.export_best_parameters(output_par)
