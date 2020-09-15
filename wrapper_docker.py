
import os
import multiprocessing
import yaml
import sys
import psutil
import math
import shutil
import logging

##We import the module of python processing.
sys.path.append('/pylcmsprocessing')
from pylcmsprocessing.model.helper.UI import UI
from pylcmsprocessing.model.experiment import Experiment,check_peakpicking
from pylcmsprocessing.model.steps.optimize_parameters import ParametersOptimizer
from pylcmsprocessing.model.helper.parameters_handler import ParametersFileHandler
from pylcmsprocessing.common.time_evaluation import Timer
import pylcmsprocessing.common.references as pcr

if __name__=="__main__":
##Two thing to check the number of CPUs and the ocnsumed meory eventually.
    ###We determine the amount of memory allocated to each process
    logging.basicConfig(format = "%(asctime)s|%(levelname)s: %(message)s",datefmt='%Y-%m-%d|%H:%M:%S',level=logging.INFO)
    timer = Timer()
    timer.store_point("wstart")
    avail_memory = (psutil.virtual_memory()[1] >> 20)
    ###We allocate the memory to each process
    num_cpus = int(multiprocessing.cpu_count()-1)
    ###Two limits to check, the number of CPUs and the memory consumption eventually.
    #1.5 Go
    memory_by_core = 1048*1.2

    if "MEMORY" in os.environ:
        memory_by_core = int(math.floor(float(os.environ["MEMORY"])))
    else:
        ##We save it ofr optimization
        os.environ["MEMORY"] = str(math.floor(memory_by_core))


    ##This is the number of thread
    ncores = avail_memory//memory_by_core

    ###Now we check if this number is bigger than the number of thread
    if ncores <= num_cpus:
        num_cpus = ncores
        if not "MEMORY" in os.environ:
            memory_by_core = avail_memory/num_cpus

    ###Specific case on EULER cluster
    if "LSB_MAX_NUM_PROCESSORS" in os.environ and int(os.environ["LSB_MAX_NUM_PROCESSORS"])<num_cpus:
        num_cpus = int(os.environ["LSB_MAX_NUM_PROCESSORS"])

    ##Command line arguments
    if "NCORES" in os.environ and int(os.environ["NCORES"])<num_cpus:
        num_cpus = int(os.environ["NCORES"])

    if "MEMORY" in os.environ:
        memory_by_core = int(os.environ["MEMORY"])

    ###We set the JAVA option for the peak picking evnetually
    os.environ["JAVA_OPTS"] = "-Xms"+str(math.floor(memory_by_core/2))+"m -Xmx"+str(math.floor(memory_by_core)-200)+"m"
    ##We output System information
    logging.info("Total memory available: "+str(avail_memory)+" and "+str( multiprocessing.cpu_count())+" cores. The workflow will use "+str(math.floor(memory_by_core))+ " Mb by cores on "+str(num_cpus)+" cores.")

    MANDATORY_ARGS = ["INPUT", "OUTPUT"]
    if os.environ['OUTPUT'].startswith('/sauer1') or os.environ['INPUT'].startswith('/sauer1'):
        MANDATORY_ARGS.append("USERNAME")

    if not all(env in os.environ for env in MANDATORY_ARGS):
        raise Exception(",".join(MANDATORY_ARGS)+' and  are mandatory arguments')
    #The path of the mounted directory is always output
    OUTPUT_DIR = os.environ['OUTPUT']
    if OUTPUT_DIR.startswith("/sauer1"):
        if not os.path.isdir(OUTPUT_DIR):
            logging.warning("Output directory "+OUTPUT_DIR+"does not exist, it will be created.")

    LOG = os.path.join(OUTPUT_DIR,"log.txt")

    # subprocess.call("java "+os.environ["JAVA_OPTS"]+" -XshowSettings:vm -version  >> "+LOG+" 2>&1",shell=True)
    #The raw files are always mounted into raw_files
    INPUT = os.environ['INPUT']
    if INPUT.startswith("/sauer1"):
        if not os.path.isdir(INPUT):
            raise Exception('Directory '+INPUT+' does not exist' )

    ##The yaml file is always putt in the paramters.text
    PATH_YAML = os.path.join(OUTPUT_DIR,"parameters.txt")
    PATH_XML = os.path.join(OUTPUT_DIR,"batch_xml_adap.xml")
    PATH_TEMP_XML = os.path.join(OUTPUT_DIR,"batch_xml_adap_temp.xml")

    PATH_TARGET = os.path.join(INPUT,"target.csv")
    # if os.path.isfile(PATH_TARGET):
    #     logging.warning("Detected target list")

    setup_params = False
    #If the yaml parameter file already exist we just read it, else. we don t read it
    #We put a message if the output is not here.
# THE sample database is always calculated before doing any processing
    PATH_DB = os.path.join(OUTPUT_DIR, "temp_processing_db.sqlite")
    DB_STORAGE = os.path.join(OUTPUT_DIR, "temp_optim", "db_storage")

    if os.path.isdir("/sauer1") or "CLUSTER" in os.environ:
        PATH_DB = "/temp_processing_db.sqlite"
        DB_STORAGE = "/db_storage"

    # LOG = "/log.txt"
    #Procesinf of the pipeline eventually.
    path_save_db = os.path.join(OUTPUT_DIR,"processing_db.sqlite")
    if os.path.isfile(path_save_db):
        shutil.copyfile(path_save_db,PATH_DB)
    exp = Experiment(PATH_DB,save_db = path_save_db,reset=False,temp_outdir=OUTPUT_DIR)
    ###The polarity computed at this step does not need to mahke any sense.
    path_ms2 = None
    if "MS2" in os.environ:
        path_ms2 = os.environ["MS2"]
    ###We try to guess the polarity form the middle file.
    pol = exp.guess_polarity(INPUT)
    logging.info("Polarity detected: " + exp.polarity)
    exp.initialise_database(num_cpus, OUTPUT_DIR, pol, INPUT, ["ADAP"], 1, path_ms2=path_ms2)
    timer.store_point("initialisation")
    timer.print_point("initialisation")
    vui = UI(OUTPUT_DIR, INPUT, polarity=os.environ["POLARITY"], mass_spec="Exactive", num_workers=num_cpus,
         path_yaml=PATH_YAML)

    ###In all case the first table is generated.
    if not os.path.exists(vui.path_yaml):
        ###We just create a an empty yaml file
        vui.generate_yaml_files()
        vui.initialize_yaml_polarity(PATH_YAML, pol)
        # self.determine_initial_parameters()
        logging.info("Empty directory detected, performing initial guess of SLAW parameters.")
        temp_opt = ParametersOptimizer(exp, "FAKE", DB_STORAGE, num_workers=num_cpus, input_par=PATH_YAML)
        temp_opt.determine_initial_parameters(output_par=PATH_YAML)
        logging.info("Parameters file generated please check the parameters values and/or parameters range before optimization.")
        exit(0)
    else:
        ph = ParametersFileHandler(vui.path_yaml)
        ##We check what is the given peakpicker.
        peakpicking = ph.get_peakpicking()
        peakpicking = check_peakpicking(peakpicking)

        if not ph.is_optimized():
            vui.generate_yaml_files()
            vui.initialize_yaml_polarity(PATH_YAML, pol)
            PATH_INITIAL = os.path.join(OUTPUT_DIR, pcr.OUT["INITIAL_PARAMETERS"])
            dummy = shutil.copyfile(PATH_YAML,PATH_INITIAL)

            with open(vui.path_yaml, 'r') as stream:
                raw_yaml = yaml.safe_load(stream)

            num_points = int(raw_yaml["optimization"]["number_of_points"]["value"])
            max_its = int(raw_yaml["optimization"]["num_iterations"]["value"])
            optim_files = int(raw_yaml["optimization"]["files_used"]["value"])
            initial_estimation = bool(raw_yaml["optimization"]["initial_estimation"]["value"])

            PATH_OPTIM = os.path.join(OUTPUT_DIR, "temp_optim")
            if not os.path.isdir(DB_STORAGE):
                os.makedirs(DB_STORAGE)
            ###We optimize the parameters
            num_cpus = int(num_cpus)
            par_opt = ParametersOptimizer(exp, PATH_OPTIM,DB_STORAGE,num_workers=num_cpus, input_par=PATH_YAML)
            optim_string = "balanced_rsm_combined_balanced_rsm_expalign"
            if "SAMPLER" in os.environ:
                optim_string = os.environ["SAMPLER"]+"_rsm_combined_"+os.environ["SAMPLER"]+"_rsm_expalign"
                logging.info("The optimisation string is:"+optim_string)
            par_opt.optimize_parameters(output_par=vui.path_yaml, optimizer=optim_string, max_its=max_its,
                                        num_points=num_points,num_files =optim_files,num_cores=num_cpus, initial_estimation=initial_estimation)
            timer.store_point("optimization")
            timer.print_point("optimization")
            ###If there was optimiz\ation we have ot reset the otpoimization proces
            exp.reset_processing()

    if not os.path.isfile(PATH_XML):
        if peakpicking=="ADAP":
            vui.generate_MZmine_XML(path_xml=PATH_XML)
            logging.info("An ADAP batch file has been generated in the "+OUTPUT_DIR+" directory, you ca use it to refine peakpicking parameters using MZmine.")
    exp.initialise_database(num_cpus,OUTPUT_DIR,vui.polarity,INPUT,["ADAP"], 1)
    # exp.building_inputs_single_processing(PATH_XML)
    ###We always read the yaml paramters file.
    with open(vui.path_yaml, 'r') as stream:
        raw_yaml = yaml.safe_load(stream)

    if peakpicking=="ADAP":
        exp.run_mzmine("/MZmine-2.52-Linux",PATH_XML,int(num_cpus*3),log = LOG)
        exp.correct_conversion()
        exp.post_processing_peakpicking_mzmine()
        ###If there is an MS2 folder we process it
        if "MS2" in os.environ:
            exp.extract_ms2(noise_level=float(raw_yaml["peakpicking"]["noise_level_ms2"]["value"]), output=pcr.OUT["ADAP"]["MSMS"])

    if peakpicking=="OPENMS":
        ###In this case arugment are read directly
        min_fwhm = float(raw_yaml["peakpicking"]["peaks_deconvolution"]["peak_width"]["value"][0])*60
        max_fwhm = float(raw_yaml["peakpicking"]["peaks_deconvolution"]["peak_width"]["value"][1])*60
        sn = float(raw_yaml["peakpicking"]["peaks_deconvolution"]["SN"]["value"])
        min_int = float(raw_yaml["peakpicking"]["noise_level_ms1"]["value"])
        min_scan = math.floor(float(raw_yaml["peakpicking"]["traces_construction"]["min_scan"]["value"]))
        max_outlier = math.floor(float(raw_yaml["peakpicking"]["traces_construction"]["num_outliers"]["value"]))
        quant = raw_yaml["grouping"]["extracted_quantity"]["value"]
        fwhm_fac = 1.2
        if "peak_width_fac" in raw_yaml["peakpicking"]["peaks_deconvolution"]:
            fwhm_fac = float(raw_yaml["peakpicking"]["peaks_deconvolution"]["peak_width_fac"]["value"])
        ppm = float(raw_yaml["peakpicking"]["traces_construction"]["ppm"]["value"])
        exp.run_openms(min_fwhm, max_fwhm, fwhm_fac, sn, ppm, min_int, max_outlier, min_scan, quant,log = LOG)
        exp.extract_ms2(noise_level=float(raw_yaml["peakpicking"]["noise_level_ms2"]["value"]),output=pcr.OUT["OPENMS"]["MSMS"],all=True)
        exp.post_processing_peakpicking_openms()

    if peakpicking=="CENTWAVE":
        min_peakwidth = float(raw_yaml["peakpicking"]["peaks_deconvolution"]["peak_width"]["value"][0]) * 60
        max_peakwidth = float(raw_yaml["peakpicking"]["peaks_deconvolution"]["peak_width"]["value"][1]) * 60
        sn = float(raw_yaml["peakpicking"]["peaks_deconvolution"]["SN"]["value"])
        min_int = float(raw_yaml["peakpicking"]["noise_level_ms1"]["value"])
        min_scan = math.floor(float(raw_yaml["peakpicking"]["traces_construction"]["min_scan"]["value"]))
        if "peak_width_fac" in raw_yaml["peakpicking"]["peaks_deconvolution"]:
            fwhm_fac = float(raw_yaml["peakpicking"]["peaks_deconvolution"]["peak_width_fac"]["value"])
        ppm = float(raw_yaml["peakpicking"]["traces_construction"]["ppm"]["value"])
        exp.run_xcms(min_peakwidth, max_peakwidth, sn, ppm, min_int, min_scan, log=LOG)
        exp.extract_ms2(noise_level=float(raw_yaml["peakpicking"]["noise_level_ms2"]["value"]),
                        output=pcr.OUT["CENTWAVE"]["MSMS"], all=True)
        exp.post_processing_peakpicking_xcms()

    ###As openMS does not give nativelyt MS-MS
    timer.store_point("peakpicking")
    timer.print_point("peakpicking")

    intensity = str(raw_yaml["grouping"]["extracted_quantity"]["value"])
    exp.group_online(intensity=intensity,
        ppm = float(raw_yaml["grouping"]["ppm"]["value"]),
        mztol=float(raw_yaml["grouping"]["dmz"]["value"]),
        rttol=float(raw_yaml["grouping"]["drt"]["value"]),
        n_ref = int(raw_yaml["grouping"]["num_references"]["value"]),
        alpha=float(raw_yaml["grouping"]["alpha"]["value"]),
        ms2_mz_tol = float(raw_yaml["peakpicking"]['peaks_deconvolution']["ms2_mz_tol"]["value"]),
        ms2_rt_tol = float(raw_yaml["peakpicking"]['peaks_deconvolution']["ms2_rt_tol"]["value"]),
        filter_qc= float(raw_yaml["filtering"]['frac_qc']["value"]),
        fold_blank= float(raw_yaml["filtering"]['fold_blank']["value"]),
        log=LOG)
    timer.store_point("alignment")
    timer.print_point("alignment")

    ###Gap filling and isotopic pattern extraction
    exp.add_missing_informations(max_iso=int(raw_yaml["ion_annotation"]['max_isotopes']["value"]),
                                 max_charge=int(raw_yaml["ion_annotation"]['max_charge']["value"]),
                                 quant=intensity,
                                 ppm=float(raw_yaml["ion_annotation"]["ppm"]["value"]),
                                 dmz=float(raw_yaml["ion_annotation"]["dmz"]["value"]))

    timer.store_point("gap-filling")
    timer.print_point("gap-filling")

    main_adducts_str=raw_yaml["ion_annotation"]["main_adducts_"+exp.polarity]["value"]
    adducts_str = raw_yaml["ion_annotation"]["adducts_"+exp.polarity]["value"]
    successfully_processed = exp.annotate_ions(int(raw_yaml["ion_annotation"]["num_files"]["value"]),float(raw_yaml["ion_annotation"]["ppm"]["value"]),
        float(raw_yaml["ion_annotation"]["dmz"]["value"]),min_filter=raw_yaml["ion_annotation"]["min_filter"]["value"],
                adducts=adducts_str,main_adducts=main_adducts_str, max_workers=num_cpus)
    timer.store_point("annotation")
    timer.print_point("annotation")
    if successfully_processed:
        # exp.post_processing(PATH_TARGET)
        exp.clean()
