
import os
import multiprocessing
import yaml
import sys
import psutil
import math

##We import the module of python processing.
sys.path.append('/pylcmsprocessing')
from pylcmsprocessing.model.UI import UI
from pylcmsprocessing.model.experiment import Experiment

if __name__=="__main__":
##Two thing to check the number of CPUs and the ocnsumed meory eventually.

    ###We determine the amount of memory allocated to each process
    avail_memory = (psutil.virtual_memory()[1] >> 20)

    ###We allocate the memory to each process
    num_cpus = multiprocessing.cpu_count()-1

    ###Two limits to check, the number of CPUs and the memory consumption eventually.
    #1.5 Go

    ###1.5 Go by default

    memory_by_core = 1048*1.2

    if "MEMORY" in os.environ:
        memory_by_core = int(os.environ["MEMORY"])

    ##This is the number of thread
    ncores = avail_memory//memory_by_core

    ###Now we check if this number is bigger than the number of thread
    if ncores <= num_cpus:
        num_cpus = ncores
        memory_by_core = avail_memory/num_cpus
    percent_mem = math.floor(memory_by_core*100/avail_memory)
    ###We set the JAVA option for the peak picking evnetually
    os.environ["JAVA_OPTS"] = "-XX:InitialRAMPercentage="+str(percent_mem)+" -XX:MinRAMPercentage="+str(percent_mem)+" -XX:MaxRAMPercentage="+str(percent_mem)

    ##We output System information
    print("Memory available: "+str(avail_memory)+" with "+str(num_cpus)+" cores used allocated percent_mem:"+str(percent_mem))
    MANDATORY_ARGS = ["INPUT","OUTPUT","USERNAME"]

    if not all(env in os.environ for env in MANDATORY_ARGS):
        raise Exception('INPUT, OUTPUT and USERNAME are all mandatory arguments.')
    #The path of the mounted directory is always output
    OUTPUT_DIR = os.environ['OUTPUT']
    if OUTPUT_DIR.startswith("/sauer1"):
        if not os.path.isdir(OUTPUT_DIR):
            print("Output directory "+OUTPUT_DIR+"does not exist.")
    #The raw files are always mounted into raw_files
    INPUT = os.environ['INPUT']
    if INPUT.startswith("/sauer1"):
        if not os.path.isdir(INPUT):
            raise Exception('Directory '+INPUT+' does not exist.' )

    ##The yaml file is always putt in the paramters.text
    PATH_YAML = os.path.join(OUTPUT_DIR,"parameters.txt")
    PATH_XML = os.path.join(OUTPUT_DIR,"batch_xml_adap.xml")

    #Path of the outpu pytho file resuming the processing eventually.
    PATH_PYTHON = "python_script.py"
    PATH_DB = os.path.join(OUTPUT_DIR,"database_lcms_processing.sqlite")

    vui = UI(OUTPUT_DIR, INPUT, polarity="positive", mass_spec="Exactive", num_workers=num_cpus, path_yaml = PATH_YAML)

    setup_params = False
    #If the yaml parameter file already exist we just read it, else. we don t read it
    #We put a message if the output is not here.
    if not os.path.exists(vui.path_yaml):
        setup_params = True
        vui.generate_yaml_files()
        print("A parameters.txt file has been generated in the specified output directory, please complete it and rerun the docker")
    else:
        #In very case we generate an adate MZmine XML file.
        vui.generate_MZmine_XML(path_xml=PATH_XML)

        ##We read the yaml file
        with open(vui.path_yaml, 'r') as stream:
            raw_yaml = yaml.safe_load(stream)

        PATH_DB = "/processing_db.sqlite"
        #Procesinf of the pipeline eventually.
        exp = Experiment(PATH_DB,reset=True)
        ##We check how many file there is
        INPUT

        exp.initialise_database(num_cpus,OUTPUT_DIR,vui.polarity,INPUT,["ADAP"], 1)
        exp.building_inputs_single_processing(PATH_XML)
        exp.run("/MZmine-2.52-Linux",int(num_cpus),log = "/log.txt")
        exp.correct_conversion()
        exp.post_processing_peakpicking()
        exp.group(max_workers=1,mztol=float(raw_yaml["grouping"]["dmz"]["value"]),
            rttol=float(raw_yaml["grouping"]["drt"]["value"]),
            intensity=str(raw_yaml["grouping"]["extracted_quantity"]["value"]))
        polarity = raw_yaml["ion_annotation"]["polarity"]["value"]
        main_adducts_str=raw_yaml["ion_annotation"]["main_adducts_"+polarity]["value"]
        adducts_str = raw_yaml["ion_annotation"]["adducts_"+polarity]["value"]
        exp.annotate_ions(int(raw_yaml["ion_annotation"]["num_files"]["value"]),float(raw_yaml["ion_annotation"]["ppm"]["value"]),
            float(raw_yaml["ion_annotation"]["dmz"]["value"]),min_filter=raw_yaml["ion_annotation"]["min_filter"]["value"],
                    adducts=adducts_str,main_adducts=main_adducts_str, max_workers=num_cpus)
        exp.clean()
