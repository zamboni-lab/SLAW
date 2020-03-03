###This take as input the paramtters ot ptimize and return the rest of the data
import hashlib
import os
import shutil
import yaml
from random import sample
import subprocess

import model.experiment as me
import common.references as cr
import model.score_datamatrix as ms
import model.LIPO as mlp
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

def peak_picking_alignment_scoring(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,fixed_params):
    ###All the params to be optimized and computed
    # p1:peakpicking/noise_level_ms1
    # p2:peakpicking/traces_construction/ppm
    # p3:peakpicking/traces_construction/dmz
    # p4:peakpicking/traces_construction/min_scan
    # p5:peakpicking/peaks_deconvolution/SN
    # p6:peakpicking/peaks_deconvolution/peak_width_min
    # p7:peakpicking/peaks_deconvolution/rt_wavelet_min
    # p8:peakpicking/peaks_deconvolution/rt_wavelet_max
    # p9:peakpicking/peaks_deconvolution/coefficient_area_threshold
    # p10:grouping/ppm
    # p11:grouping/drt
    # p12: grouping/dmz
    # p13:grouping/alpha
    # p14:grouping/num_references

    # p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14=params
    DIR_TEMP,STORAGE_YAML,SUMMARY_YAML,num_cpus,polarity,path_samples,initial_yaml=fixed_params

    ###we create the directory andget the path which will be used in the data.
    hash_val,OUTPUT_DIR,PATH_DB,PATH_SAVE_DB,stored_param,stored_xml = create_temp_directory(DIR_TEMP,STORAGE_YAML,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14)


    ###We create a temporary directory to put the data in
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    LOG_PATH = os.path.join(OUTPUT_DIR,"log.txt")

    def convert_val(x):

        if isinstance(x,float) or isinstance(x,int):
            return x
        return x.item()

    ####We load the orginial yaml file
    raw_yaml=None
    with open(initial_yaml, 'r') as stream:
        raw_yaml = yaml.safe_load(stream)
    ##We update the parameters to reflect the rest of the data
    raw_yaml["peakpicking"]["noise_level_ms1"]["value"]=float(convert_val(p1))
    raw_yaml["peakpicking"]["traces_construction"]["ppm"]["value"]=float(convert_val(p2))
    raw_yaml["peakpicking"]["traces_construction"]["dmz"]["value"]=float(convert_val(p3))
    raw_yaml["peakpicking"]["traces_construction"]["min_scan"]["value"]=int(convert_val(p4))
    raw_yaml["peakpicking"]["peaks_deconvolution"]["SN"]["value"]=float(convert_val(p5))
    raw_yaml["peakpicking"]["peaks_deconvolution"]["peak_width_min"]["value"]=float(convert_val(p6))
    raw_yaml["peakpicking"]["peaks_deconvolution"]["rt_wavelet_min"]["value"]=float(convert_val(p7))
    raw_yaml["peakpicking"]["peaks_deconvolution"]["rt_wavelet_max"]["value"]=float(convert_val(p8))
    raw_yaml["peakpicking"]["peaks_deconvolution"]["coefficient_area_threshold"]["value"]=float(convert_val(p9))
    ###We dump the file in a filde evnetually
    with open(stored_param, 'w') as outfile:
        yaml.dump(raw_yaml, outfile, default_flow_style=False)

    PATH_DB = "/temp_processing"+str(hash_val)+"_db.sqlite"
    if "CLUSTER" in os.environ:
        PATH_DB = os.path.join(OUTPUT_DIR,"/temp_processing"+str(hash_val)+"_db.sqlite")
    # LOG = "/log.txt"
    #Procesinf of the pipeline eventually.
    PATH_SAVE_DB = os.path.join(OUTPUT_DIR,"processing_db.sqlite")
    PATH_XML = os.path.join(OUTPUT_DIR,"/temp_par"+str(hash_val)+".xml")

    exp = me.Experiment(PATH_DB,save_db = PATH_SAVE_DB,reset=False)
    ###We create the UI
    vui = UI(OUTPUT_DIR, path_samples, polarity=polarity, mass_spec="Exactive", num_workers=num_cpus, path_yaml=stored_param)
    vui.generate_yaml_files(cr.DATA["YAML"])
    ###The output database is generated.
    # to_join = [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14]
    # print(to_join)
    # tto_join = [convert_val(pp) for pp in to_join]
    # print(tto_join)
    vui.generate_MZmine_XML(path_xml=PATH_XML)
    exp.initialise_database(num_cpus, OUTPUT_DIR, polarity, path_samples, ["ADAP"], 1)
    exp.building_inputs_single_processing(PATH_XML)
    try:
        exp.run("/MZmine-2.52-Linux", int(num_cpus), log=LOG_PATH)
        exp.correct_conversion()
        exp.post_processing_peakpicking()
    except Exception:
        shutil.rmtree(OUTPUT_DIR)
        return -1
    exp.group_online(intensity="intensity",
                     ppm=float(convert_val(p10)),
                     mztol=float(convert_val(p12)),
                     rttol=float(convert_val(p11)),
                     n_ref=int(convert_val(p14)),
                     alpha=float(convert_val(p13)),
                     log=LOG_PATH)
    ###We load all the datamtrices.
    datamatrices = exp.get_datamatrix()
    path_datamatrix = datamatrices[0][0]
    try:
        msd = ms.scorerDataMatrix(path_datamatrix)
        ###We authorize 1 jump
        vscore = msd.score_datamatrix(1)
    except ValueError:
        shutil.rmtree(OUTPUT_DIR,ignore_errors=True)
        ###Case where no peaktable has been found.
        return -1

    ###We write the score in the table
    to_join = [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14]
    print(to_join)
    to_join = [str(pp) for pp in to_join]
    to_write = stored_param+","+str(vscore)+",".join(to_join)+"\n"
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
            num_files = self.num_workers
        if len(self.samples)<num_files:
            num_files = len(self.samples)

        # In a first sketch we just select 5 files randomly.
        sel_samples = sample(self.samples,num_files)
        for ss in sel_samples:
            shutil.copy(ss[0],self.path_samples)

    def optimize_initial_parameters(self):
        pscript = os.path.join(ct.find_rscript(),"initial_parameters.R")

        # we tun the optimization ins a signle thread.
        path_db = self.path_db
        output_par = self.temp_yaml
        initial_par = self.input_par
        num_cores = self.num_workers
        cli = " ".join(["Rscript",pscript,path_db,initial_par,output_par,str(num_cores)])
        ###We read the polarity directly
        subprocess.call(cli, shell=True, env=os.environ.copy())

    def optimize_tricky_parameters(self,limits=None):
        if limits is None:
            ###We read the standar doptimization interval eventually
            limits = cr.ALGORITHMS_TABLE["ADAP"][3]
            with open(limits, 'r') as stream:
                limits = yaml.safe_load(stream)
        # We build the bound object
        lb = [limits[ll][0] for ll in limits]
        ub = [limits[ll][1] for ll in limits]

        # We update the limits vector
        ###if the yaml has not been optimize we create the first value
        if not os.path.isfile(self.temp_yaml):
            self.temp_yaml = self.input_par

        # We create the same dataset eventually
        fixed_params = self.path_exp, self.params_archive, self.path_summary, self.num_workers, self.polarity, self.path_samples, self.temp_yaml

        # we optimize the point
        voptim = mlp.LIPO(lb,ub,peak_picking_alignment_scoring,self.nrounds,3,fixed_arguments={"fixed_params":fixed_params})

        # The best paramters is passed along the data
        self.best_par = voptim


    def export_best_parameters(self,outpath):
        p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14 = self.best_par
        with open(self.input_par, 'r') as stream:
            raw_yaml = yaml.safe_load(stream)
        ##We update the parameters to reflect the rest of the data
        raw_yaml["peakpicking"]["noise_level_ms1"]["value"] = p1
        raw_yaml["peakpicking"]["traces_construction"]["ppm"]["value"] = p2
        raw_yaml["peakpicking"]["traces_construction"]["dmz"]["value"] = p3
        raw_yaml["peakpicking"]["traces_construction"]["min_scan"]["value"] = p4
        raw_yaml["peakpicking"]["peaks_deconvolution"]["SN"]["value"] = p5
        raw_yaml["peakpicking"]["peaks_deconvolution"]["peak_width_min"]["value"] = p6
        raw_yaml["peakpicking"]["peaks_deconvolution"]["rt_wavelet_min"]["value"] = p7
        raw_yaml["peakpicking"]["peaks_deconvolution"]["rt_wavelet_max"]["value"] = p8
        raw_yaml["peakpicking"]["peaks_deconvolution"]["coefficient_area_threshold"]["value"] = p9
        raw_yaml["grouping"]["ppm"]["value"] = p10
        raw_yaml["grouping"]["drt"]["value"] = p11
        raw_yaml["grouping"]["dmz"]["value"] = p12
        raw_yaml["grouping"]["alpha"]["value"] = p13
        raw_yaml["grouping"]["num_references"]["value"] = p14
        ###We dump the file in a filde evnetually
        with open(outpath, 'w') as outfile:
            yaml.dump(raw_yaml, outfile, default_flow_style=False)

    def optimize_parameters(self,output_par):
        self.optimize_initial_parameters()
        self.select_samples()
        print("Finished initial parameters estimation.")
        print("Optimizing remaining parameters.")
        self.optimize_tricky_parameters()
        print("Finished  optimization")
        self.export_best_parameters(output_par)


        ###We strat by doing the initial optimization of parameters
