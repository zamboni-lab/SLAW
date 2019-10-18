
import os
import multiprocessing
import yaml
import sys

##We import the module of python processing.
sys.path.append('/pylcmsprocessing')
from pylcmsprocessing.model.UI import UI
from pylcmsprocessing.model.experiment import Experiment

if __name__=="__main__":

    #The path of the mounted directory is always output
    OUTPUT_DIR = "/output"

    #The raw files are always mounted into raw_files
    RAW_FILES = "/rawfiles"

    ##The yaml file is always putt in the paramters.text
    PATH_YAML = os.path.join(OUTPUT_DIR,"parameters.txt")
    PATH_XML = os.path.join(OUTPUT_DIR,"batch_xml_adap.xml")

    #Path of the outpu pytho file resuming the processing eventually.
    PATH_PYTHON = "python_script.py"
    PATH_DB = os.path.join(OUTPUT_DIR,"database_lcms_processing.sqlite")
    num_cpus = multiprocessing.cpu_count()-1
    if "NCORES" in os.environ:
        num_cpus = int(os.environ['NCORES'])
    else:
        print(f'No NCORES environment variables, the number of cores used for parallel processing has been automatically set to {num_cpus}')
    vui = UI(OUTPUT_DIR, RAW_FILES, polarity="positive", mass_spec="Exactive", num_workers=num_cpus, path_yaml = PATH_YAML)

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

        #Procesinf of the pipeline eventually.
        exp = Experiment(PATH_DB,reset=True)
        exp.initialise_database(num_cpus,OUTPUT_DIR,vui.polarity,RAW_FILES,["ADAP"], 1)
        exp.building_inputs_single_processing(PATH_XML)
        exp.run("/MZmine-2.51-Linux",num_cpus)
        exp.correct_conversion()
        exp.group(max_workers=1,mztol=float(raw_yaml["grouping"]["dmz"]["value"]),
            rttol=float(raw_yaml["grouping"]["drt"]["value"]),
            intensity=str(raw_yaml["grouping"]["extracted_quantity"]["value"]))
        polarity = raw_yaml["ion_annotation"]["polarity"]["value"]
        main_adducts_str=raw_yaml["ion_annotation"]["main_adducts_"+polarity]["value"]
        adducts_str = raw_yaml["ion_annotation"]["adducts_"+polarity]["value"]
        exp.annotate_ions(int(raw_yaml["ion_annotation"]["num_files"]["value"]),float(raw_yaml["ion_annotation"]["ppm"]["value"]),
            float(raw_yaml["ion_annotation"]["dmz"]["value"]),min_filter=raw_yaml["ion_annotation"]["min_filter"]["value"],
                    adducts=adducts_str,main_adducts=main_adducts_str, max_workers=num_cpus)
