###This take as input the paramtters ot ptimize and return the rest of the data
import hashlib
import os
import shutil
import yaml
from random import sample,seed
import subprocess
import numpy as np
import pandas as pd
import time
import math

import model.experiment as me
import common.references as cr
import model.optimization.score_experiment as ms
import model.optimization.sampling as mos
import model.helper.parameters_handler as ph
import common.tools as ct
from model.helper.UI import UI

###Function used to optimize the paramters
def create_temp_directory(path_exp,params_archive,dir_db,*argv):
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
    temp_db = os.path.join(dir_db,"processing_db_"+str(hash_val)+".sqlite")
    temp_save_db = os.path.join(temp_dir,"save_processing_db_"+str(hash_val)+".sqlite")
    stored_xml = os.path.join(params_archive,"xml_"+str(hash_val)+".xml")
    return hash_val,temp_dir,temp_db,temp_save_db,stored_param,stored_xml


def convert_val(x):
    if type(x) is np.ndarray:
        if isinstance(x[0], np.int):
            return int(x.item())
        if isinstance(x[0], np.float64):
            return float(x.item())
        return x.item()
    elif type(x) is np.float64:
        return float(x.item())
    return x



# def peak_picking_alignment_scoring(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,fixed_params):
def peak_picking_alignment_scoring(peakpicking__noise_level_ms1,peakpicking__noise_level_ms2,peakpicking__traces_construction__ppm,
                                   peakpicking__traces_construction__dmz,peakpicking__traces_construction__min_scan,
                                   peakpicking__peaks_deconvolution__SN,peakpicking__peaks_deconvolution__noise_level,
                                   peakpicking__peaks_deconvolution__peak_width__const,
                                   peakpicking__peaks_deconvolution__peak_width__add,
                                   peakpicking__peaks_deconvolution__rt_wavelet__const,
                                   peakpicking__peaks_deconvolution__rt_wavelet__add,
                                   peakpicking__peaks_deconvolution__coefficient_area_threshold,grouping__ppm,
                                   grouping__drt,grouping__dmz,grouping__alpha,grouping__num_references,fixed_params,
                                   **kwargs):
    args_refs = locals()
    del args_refs["fixed_params"]
    del args_refs["kwargs"]
    def is_processed(stored_param, summary_table):
        if not os.path.isfile(summary_table):
            return False, 0
        try:
            rr = pd.read_csv(summary_table)
        except pd.errors.EmptyDataError:
            return False, 0
        sub_row = rr[[stored_param in xx for xx in rr.path_parameters]]
        if sub_row.shape[0]==0:
            return False, 0
        return True, sub_row.iloc[0,4]

    def convert_val(x):
        if type(x) is np.ndarray:
            if isinstance(x[0], np.int):
                return int(x.item())
            if isinstance(x[0], np.float64):
                return float(x.item())
            return x.item()
        if type(x) is np.float64:
            return float(x)
        return x

    # p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14=params
    DIR_TEMP,STORAGE_YAML,STORAGE_DB,SUMMARY_YAML,num_cpus,polarity,path_samples,initial_yaml,parallel,scorer,pdb,reset=fixed_params
    if parallel:
        num_cpus=1
    call_list = [DIR_TEMP,STORAGE_YAML,STORAGE_DB]+list(args_refs.values())


    ###we create the directory andget the path which will be used in the data.
    hash_val,OUTPUT_DIR,PATH_DB,PATH_SAVE_DB,stored_param,stored_xml = create_temp_directory(*call_list)
    if not reset:
        processed = is_processed(stored_param,SUMMARY_YAML)
        if processed[0]:
            return float(processed[1])

    ###We write the header
    if not os.path.isfile(SUMMARY_YAML):
        ##We write the path name.
        to_join = ["path_parameters", "db", "db_save","output_dir", "score"] + list(args_refs.keys())
        joined_summary = ",".join(to_join)
        with open(SUMMARY_YAML, "w") as ff:
            ff.write(joined_summary + "\n")
    ###We create a temporary directory to put the data in
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    LOG_PATH = os.path.join(OUTPUT_DIR,"log.txt")

    ####We load the orginial yaml file
    parameters = ph.ParametersFileHandler(initial_yaml)
    l_values = [convert_val(ll) for ll in args_refs.values()]

    parameters.set_all_parameters(list(args_refs.keys()),list(l_values))

    ###The numbver of scans is handled separately as it is an integer.
    parameters["peakpicking__traces_construction__min_scan"] =int(convert_val(args_refs["peakpicking__traces_construction__min_scan"]))

    ###We dump the yaml
    parameters.write_parameters(stored_param)
    PATH_XML = os.path.join(OUTPUT_DIR,"temp_par"+str(hash_val)+".xml")

    exp = me.Experiment(PATH_DB,save_db = PATH_SAVE_DB,reset=False)
    ###We create the UI
    vui = UI(OUTPUT_DIR, path_samples, polarity=polarity, mass_spec="Exactive", num_workers=num_cpus,
             path_yaml=stored_param)
    # vui.generate_yaml_files(cr.DATA["YAML"])
    ###The output database is generated.
    vui.generate_MZmine_XML(path_xml=PATH_XML)

    ###If needed we tune the cups
    if num_cpus!=1:
        lf = os.listdir(path_samples)
        if num_cpus>len(lf):
            num_cpus=len(lf)

    exp.initialise_database(num_cpus, OUTPUT_DIR, polarity, path_samples, ["ADAP"], 1)
    exp.building_inputs_single_processing(PATH_XML)

    to_join = list(args_refs.values())
    to_join = [str(pp) for pp in to_join]


    if os.path.isfile(pdb):
        exp.rebase_experiment(pdb)
    else:
        try:
            exp.run("/MZmine-2.52-Linux", int(num_cpus), log=LOG_PATH)
            exp.correct_conversion()
            exp.post_processing_peakpicking()
        except Exception as e:
            print(e)
            try:
                if not reset:
                    os.remove(PATH_DB)
                    subprocess.call("rm -r " + OUTPUT_DIR)
            except Exception:
                pass
            # shutil.rmtree(OUTPUT_DIR)
            to_write = stored_param + "," + PATH_DB + "," + PATH_SAVE_DB + "," + OUTPUT_DIR + "," + str(
                -1) + "," + ",".join(to_join) + "\n"
            with open(SUMMARY_YAML, "a") as ff:
                ff.write(to_write)
            return -1

    ##We check te scorer if it is a peak sorer we don t need the data matrix
    if scorer!="ipopeak":
        exp.group_online(intensity="int",
                         ppm=float(convert_val(grouping__ppm)),
                         mztol=float(convert_val(grouping__dmz)),
                         rttol=float(convert_val(grouping__drt)),
                         n_ref=int(convert_val(grouping__num_references)),
                         alpha=float(convert_val(grouping__alpha)),
                         log=LOG_PATH)
    tvscore = -1.0
    try:
        msd = ms.get_scorer(scorer)
        msd = msd(exp)
        ###We authorize 1 jump
        vscore = msd.score(output=OUTPUT_DIR)
        tvscore = vscore
        if not isinstance(tvscore,float):
            tvscore = "|".join(["{0:.10g}".format(sc) for sc in tvscore])
    except Exception as e:
            ###We write the score in the table
        to_write = stored_param + "," +PATH_DB+","+PATH_SAVE_DB+","+OUTPUT_DIR+","+str(tvscore) +","+",".join(to_join) + "\n"
        with open(SUMMARY_YAML, "a") as ff:
            ff.write(to_write)
        if parallel:
            try:
                if not reset:
                    os.remove(PATH_DB)
                    shutil.rmtree(OUTPUT_DIR)
            except Exception:
                pass
        ###In all case we write the parameter in
        return -1
    ###We write the score in the table
    to_write = stored_param + "," +PATH_DB+","+PATH_SAVE_DB+","+OUTPUT_DIR+","+str(tvscore) +","+",".join(to_join) + "\n"


    with open(SUMMARY_YAML, "a") as ff:
        ff.write(to_write)
    if parallel:
        try:
            if not reset:
                os.remove(PATH_DB)
                shutil.rmtree(OUTPUT_DIR)
        except Exception:
            pass
    return vscore


def parse_optim_option(string):
    psampler,poptimizer,pscorer,gsampler,goptimizer,gscorer = string.split("_")
    return mos.getSampler(psampler),mos.getOptimizer(poptimizer),pscorer,mos.getSampler(gsampler),mos.getOptimizer(goptimizer),gscorer

class ParametersOptimizer:
    def __init__(self,exp,output,db_storage,num_workers=1,input_par=None):
        if input_par is None:
            input_par = cr.DATA["YAML"]
        self.input_par = input_par
        self.output = output
        # At the initialisation a temp directory for the output is created.
        self.path_exp=os.path.join(self.output,"temp_exp")
        self.dir_db=db_storage
        os.makedirs(self.path_exp, exist_ok=True)
        self.params_archive=os.path.join(self.output,"parameters")
        os.makedirs(self.params_archive, exist_ok=True)
        self.path_samples = os.path.join(self.output,"mzML")
        os.makedirs(self.path_samples, exist_ok=True)
        self.path_summary = os.path.join(self.output,"summary_par")
        self.temp_polarity = os.path.join(self.output,"polarity.csv")
        self.polarity = exp.get_polarity()
        self.num_workers = num_workers
        self.samples = exp.get_query("SELECT path FROM samples WHERE types='QC'")
        if len(self.samples)==0:
            self.samples = exp.get_query("SELECT path FROM samples WHERE types='sample'")
        self.path_db = exp.db
        self.temp_yaml = os.path.join(self.output,"temp_yaml.yaml")
        shutil.copy(input_par,self.temp_yaml)


    # Select the smaples ot optimize evetually.
    def select_samples(self,num_files=None):
        if num_files is None:
            # As many as worker to speed up the peakpicking process
            num_files = 10
        if len(self.samples) < num_files:
            num_files = len(self.samples)
        seed(num_files)
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


    ###Thoe optimization is always two steps
    def do_optimize_parameters(self,optim_string="bbd_rsm",max_its=10,num_points=50,**kwargs):

        ### We get the optimization algorithm.
        psampler,poptimizer,pscorer,gsampler,goptimizer,gscorer = parse_optim_option(optim_string)
        psampler = psampler()
        poptimizer = poptimizer()
        gsampler = gsampler()
        goptimizer = goptimizer()
        pdb = "NONE"

        ## Reading parameters to optimize and ranges
        pfh = ph.ParametersFileHandler(self.temp_yaml)
        to_optimize = pfh.get_optimizable_parameters(string=True)
        ###Parameterd are split between peakpciking and grouping
        to_optimize_peakpicking = {ll:to_optimize[ll] for ll in to_optimize if ll.startswith("peakpicking")}
        to_optimize_grouping = {ll:to_optimize[ll] for ll in to_optimize if ll.startswith("grouping")}
        all_parameters = pfh.get_parameters_values()

        final_dic = all_parameters.copy()

        ###if the yaml has not been optimized we create the first value
        if not os.path.isfile(self.temp_yaml):
            self.temp_yaml = self.input_par

        ##We get the bounds of the optimizable parameters.
        summary_peakpicking = self.path_summary + "_peakpicking.csv"

        if len(to_optimize_peakpicking) > 0 and pscorer!="none":
            lb_peakpicking = [x[0] for x in to_optimize_peakpicking.values()]
            ub_peakpicking = [x[1] for x in to_optimize_peakpicking.values()]
            bounds_peakpicking = mos.bounds(lb_peakpicking, ub_peakpicking, list(to_optimize_peakpicking.keys()))

            # Paramters which are always fixed.
            fixed_params_peakpicking = (self.path_exp, self.params_archive, self.dir_db, summary_peakpicking,\
                           self.num_workers, self.polarity, self.path_samples, self.temp_yaml, psampler.is_parallel(), pscorer, pdb, False)

            dic_fixed_peakpicking = {"fixed_params":fixed_params_peakpicking}
            ###We fix the paramters which will not be optimized
            for ll in all_parameters:
                if ll not in to_optimize_peakpicking:
                    dic_fixed_peakpicking[ll] = all_parameters[ll]
            try:
                psoptim = mos.samplingOptimizer(psampler, poptimizer, bounds_peakpicking, fixed_arguments=dic_fixed_peakpicking)
                pvoptim = psoptim.optimize(peak_picking_alignment_scoring,max_its=max_its,num_points=num_points, num_cores=self.num_workers)
            except Exception as e:
                print("Exception occured:",e)
                pass
            print("Peakpicking optimization finished.")
        else:
            print("No peakpicking optimization required.")

        if len(to_optimize_grouping) > 0 and gscorer!="none":
            ###We compute the peak picking a single time
            dic_pp = all_parameters.copy()
            fixed_pp_single = (self.path_exp, self.params_archive, self.dir_db, summary_peakpicking,\
                           self.num_workers, self.polarity, self.path_samples, self.temp_yaml, False, "ipopeak", pdb, True)

            if ("pvoptim" in locals()):
                for ik in pvoptim:
                    dic_pp[ik] = convert_val(pvoptim[ik])
            ####We add the correct fixed parameters
            dic_pp["fixed_params"] = fixed_pp_single
            dic_pp["reset"]=True
            vss = peak_picking_alignment_scoring(**dic_pp)

            lb_grouping = [x[0] for x in to_optimize_grouping.values()]
            ub_grouping = [x[1] for x in to_optimize_grouping.values()]
            bounds_grouping = mos.bounds(lb_grouping, ub_grouping, list(to_optimize_grouping.keys()))
            summary_alignment = self.path_summary+"_alignement.csv"
            fixed_params_grouping = (self.path_exp, self.params_archive, self.dir_db, summary_alignment,\
                           self.num_workers, self.polarity, self.path_samples, self.temp_yaml, gsampler.is_parallel(), gscorer, pdb, False)
            dic_fixed_grouping = {"fixed_params":fixed_params_grouping}

            ##We read the correct parameter to copy the experiment
            if os.path.isfile(summary_peakpicking):
                pda = pd.read_csv(summary_peakpicking)
                idmax = pda.score.idxmax()
                pdb = pda.db[idmax]

                ###In any case we recompute a single parameters with enough cores
                ###In every case we recompute a single paramters with the best peakpicking
                dic_call = dic_fixed_grouping

                fixed_params_grouping = (self.path_exp, self.params_archive, self.dir_db, summary_alignment, \
                                         self.num_workers, self.polarity, self.path_samples, self.temp_yaml,gsampler.is_parallel(), gscorer,pdb, False)
                dic_fixed_grouping = {"fixed_params":fixed_params_grouping}
            ###We fix the paramter which are not yet fixed
            for ll in all_parameters:
                if ll not in to_optimize_grouping:
                    dic_fixed_grouping[ll] = all_parameters[ll]

            # THe best peakpicking parameters are inject into the grouping parameters.
            if ("pvoptim" in locals()):
                for ik in pvoptim:
                    dic_fixed_grouping[ik] = convert_val(pvoptim[ik])

            gsoptim = mos.samplingOptimizer(gsampler, goptimizer, bounds_grouping, fixed_arguments=dic_fixed_grouping)
            try:
                gvoptim = gsoptim.optimize(peak_picking_alignment_scoring,max_its=max_its, num_points=num_points, num_cores=self.num_workers)
                print("Grouping optimization finished.")
                final_dic = dic_fixed_grouping
                for ik in gvoptim:
                    final_dic[ik] = convert_val(gvoptim[ik])
            except Exception:
                pass
        # The best paramters is stroed
        self.best_par = final_dic
        print(final_dic)

    def export_best_parameters(self,outpath):
        pfh = ph.ParametersFileHandler(self.input_par)
        ll_values = [convert_val(ll) for ll in self.best_par.values()]
        pfh.set_all_parameters(list(self.best_par.keys()), ll_values)
        # for k in self.best_par:
        #     if k == "fixed_params":
        #         continue
        #     pfh[k] = convert_val(self.best_par[k])
        pfh.dic["optimized"] = True
        pfh.write_parameters(outpath)


    def optimize_parameters(self,output_par,optimizer="lipo_rsm",max_its=10,num_files= 10,num_points=100,**kwargs):

        ###The memory of JAVA is always more limited for optimization, 700 less Mbs.
        memory_by_core = int(os.environ["MEMORY"])*0.7

        ###We set the JAVA option for the peak picking evnetually
        os.environ["JAVA_OPTS"] = "-Xms" + str(math.floor(memory_by_core / 2)) + "m -Xmx" + str(
            math.floor(memory_by_core) - 200) + "m"

        self.determine_initial_parameters()
        self.select_samples(num_files = num_files)
        print("Finished initial parameters estimation")
        print("Optimizing remaining parameters")
        self.do_optimize_parameters(optim_string=optimizer,max_its=max_its,num_points=num_points)
        print("Finished  optimization")
        self.export_best_parameters(output_par)

        ###We reset the jAVA memoery limit
        memory_by_core = int(os.environ["MEMORY"])

        ###We set the JAVA option for the peak picking evnetually
        os.environ["JAVA_OPTS"] = "-Xms" + str(math.floor(memory_by_core / 2)) + "m -Xmx" + str(
            math.floor(memory_by_core) - 200) + "m"