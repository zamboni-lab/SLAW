
import os
import multiprocessing
import yaml
import sys
import psutil
import math
import shutil
import time


##We import the module of python processing.
sys.path.append('/pylcmsprocessing')
from pylcmsprocessing.model.UI import UI
from pylcmsprocessing.model.experiment import Experiment
from pylcmsprocessing.model.optimization import ParametersOptimizer

if __name__=="__main__":
##Two thing to check the number of CPUs and the ocnsumed meory eventually.
    ###We determine the amount of memory allocated to each process
    time_start = time.clock()
    avail_memory = (psutil.virtual_memory()[1] >> 20)
    ###We allocate the memory to each process
    num_cpus = multiprocessing.cpu_count()-1
    ###Two limits to check, the number of CPUs and the memory consumption eventually.
    #1.5 Go
    memory_by_core = 1048*1.2

    if "MEMORY" in os.environ:
        memory_by_core = int(os.environ["MEMORY"])

    ##This is the number of thread
    ncores = avail_memory//memory_by_core

    ###Now we check if this number is bigger than the number of thread
    if ncores <= num_cpus:
        num_cpus = ncores
        if not "MEMORY" in os.environ:
            memory_by_core = avail_memory/num_cpus

    if "LSB_MAX_NUM_PROCESSORS" in os.environ:
        num_cpus = int(os.environ["LSB_MAX_NUM_PROCESSORS"])

    if "NCORES" in os.environ:
        num_cpus = int(os.environ["NCORES"])

    if "MEMORY" in os.environ:
        memory_by_core = int(os.environ["MEMORY"])

    ###We set the JAVA option for the peak picking evnetually
    os.environ["JAVA_OPTS"] = "-Xms"+str(math.floor(memory_by_core/2))+"m -Xmx"+str(math.floor(memory_by_core))+"m"
    ##We output System information
    print("Total memory available: "+str(avail_memory)+" and "+str( multiprocessing.cpu_count())+" cores. The workflow will use "+str(math.floor(memory_by_core))+ " Mb by cores on "+str(num_cpus)+" cores.")

    MANDATORY_ARGS = ["INPUT", "OUTPUT"]
    if os.environ['OUTPUT'].startswith('/sauer1') or os.environ['INPUT'].startswith('/sauer1'):
        MANDATORY_ARGS.append("USERNAME")

    if not all(env in os.environ for env in MANDATORY_ARGS):
        raise Exception(",".join(MANDATORY_ARGS)+' and  are mandatory arguments.')
    #The path of the mounted directory is always output
    OUTPUT_DIR = os.environ['OUTPUT']
    if OUTPUT_DIR.startswith("/sauer1"):
        if not os.path.isdir(OUTPUT_DIR):
            print("Output directory "+OUTPUT_DIR+"does not exist.")

    LOG = os.path.join(OUTPUT_DIR,"log.txt")

    # subprocess.call("java "+os.environ["JAVA_OPTS"]+" -XshowSettings:vm -version  >> "+LOG+" 2>&1",shell=True)
    #The raw files are always mounted into raw_files
    INPUT = os.environ['INPUT']
    if INPUT.startswith("/sauer1"):
        if not os.path.isdir(INPUT):
            raise Exception('Directory '+INPUT+' does not exist.' )

    ##The yaml file is always putt in the paramters.text
    PATH_YAML = os.path.join(OUTPUT_DIR,"parameters.txt")
    PATH_XML = os.path.join(OUTPUT_DIR,"batch_xml_adap.xml")
    PATH_TEMP_XML = os.path.join(OUTPUT_DIR,"batch_xml_adap_temp.xml")

    PATH_TARGET = os.path.join(INPUT,"target.csv")
    if os.path.isfile(PATH_TARGET):
        print("Detected target list.")

    vui = UI(OUTPUT_DIR, INPUT, polarity=os.environ["POLARITY"], mass_spec="Exactive", num_workers=num_cpus, path_yaml = PATH_YAML)
    setup_params = False
    #If the yaml parameter file already exist we just read it, else. we don t read it
    #We put a message if the output is not here.
# THE sample database is always calculated before doing any processing
    PATH_DB = "/temp_processing_db.sqlite"
    if "CLUSTER" in os.environ:
        PATH_DB = os.path.join(OUTPUT_DIR,"temp_processing_db.sqlite")
    # LOG = "/log.txt"
    #Procesinf of the pipeline eventually.
    path_save_db = os.path.join(OUTPUT_DIR,"processing_db.sqlite")
    if os.path.isfile(path_save_db):
        shutil.copyfile(path_save_db,PATH_DB)
    exp = Experiment(PATH_DB,save_db = path_save_db,reset=False)
    exp.initialise_database(num_cpus, OUTPUT_DIR, vui.polarity, INPUT, ["ADAP"], 1)
    ###In all case the first table is generated.
    if not os.path.exists(vui.path_yaml):
        vui.generate_yaml_files()
        num_iter = 10
        if "NOPTIM" in os.environ:
            num_iter = int(num_iter)
        PATH_OPTIM = os.path.join(OUTPUT_DIR, "temp_optim")
        os.makedirs(PATH_OPTIM)
        par_opt = ParametersOptimizer(exp, PATH_OPTIM, nrounds=num_iter, input_par=None)
        par_opt.optimize_parameters(vui.path_yaml)


        ###In this case we optimize the parameter
        ##We first check fi there is anumber of iteration defined

    if not os.path.isfile(PATH_XML):
        vui.generate_MZmine_XML(path_xml=PATH_XML)
        print("An ADAP batch file has been generated in the "+OUTPUT_DIR+" directory, you ca use it to refine peakpicking parameters using MZmine.")
        # print("A parameters.txt file has been generated in the "+OUTPUT_DIR+" directory, please complete it and rerun the docker.")

    exp.initialise_database(num_cpus,OUTPUT_DIR,vui.polarity,INPUT,["ADAP"], 1)
    exp.building_inputs_single_processing(PATH_XML)
    exp.run("/MZmine-2.52-Linux",int(num_cpus),log = LOG)
    exp.correct_conversion()
    exp.post_processing_peakpicking()
    time_peakpicking = time.clock()

    exp.group_online(intensity=str(raw_yaml["grouping"]["extracted_quantity"]["value"]),
        ppm = float(raw_yaml["grouping"]["ppm"]["value"]),
        mztol=float(raw_yaml["grouping"]["dmz"]["value"]),
        rttol=float(raw_yaml["grouping"]["drt"]["value"]),
        n_ref = int(raw_yaml["grouping"]["num_references"]["value"]),
        alpha=float(raw_yaml["grouping"]["alpha"]["value"]),
        log=LOG)
    time_grouping = ()
    polarity = raw_yaml["ion_annotation"]["polarity"]["value"]
    main_adducts_str=raw_yaml["ion_annotation"]["main_adducts_"+polarity]["value"]
    adducts_str = raw_yaml["ion_annotation"]["adducts_"+polarity]["value"]
    successfully_processed = exp.annotate_ions(int(raw_yaml["ion_annotation"]["num_files"]["value"]),float(raw_yaml["ion_annotation"]["ppm"]["value"]),
        float(raw_yaml["ion_annotation"]["dmz"]["value"]),min_filter=raw_yaml["ion_annotation"]["min_filter"]["value"],
                adducts=adducts_str,main_adducts=main_adducts_str, max_workers=num_cpus)
    time_annotation = time.clock()
    if successfully_processed:
        exp.post_processing(PATH_TARGET)
        exp.clean()
