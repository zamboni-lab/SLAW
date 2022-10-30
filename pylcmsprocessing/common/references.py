import os

if not('SLAWROOT' in os.environ):
    os.environ['SLAWROOT'] = '/' # this is the default for the docker container

ALGORITHMS_TABLE = {
    "ADAP": ["MZmine", "adap.json", os.environ['SLAWROOT']+ "pylcmsprocessing/data/batch_adap_2_52.xml"],
    "BASELINE": ["MZmine", "baseline.json", os.environ['SLAWROOT']+ "pylcmsprocessing/data/batch_baseline_2_52.xml"],
    "FeatureFinderMetabo": ["openMS", "openms.json", None, None],
    "centWave": ["centWave",None,None,None]
}

MZMINE_ALGORITHM = ["BASELINE","ADAP"]

TEMPLATE_TYPE = os.environ['SLAWROOT']+ "pylcmsprocessing/data/templates/template_types.yaml"
TEMPLATES = {
    "baseline": os.environ['SLAWROOT']+ "pylcmsprocessing/data/templates/template_adap.yaml",
    "adap":os.environ['SLAWROOT']+ "pylcmsprocessing/data/templates/template_adap.yaml",
    "openms":os.environ['SLAWROOT']+ "pylcmsprocessing/data/templates/template_openms.yaml",
    "centwave":os.environ['SLAWROOT']+ "pylcmsprocessing/data/templates/template_centwave.yaml"
}

SUMMARY_DIR = "summary_parameters"
SUMMARY_COMBINATIONS = "summary_combinations"

OUT={"DB":"processing.sqlite",
     "LOG":"log.txt",
     "POLARITY":"polarity",
     "DONE":"done",
     "BASELINE":
         {
             "SUMMARY":"BASELINE/BASELINE_parameters.csv",
             "SUMMARY_TEMPLATES": "BASELINE/summary_templates.csv",
             "JSON":"BASELINE/temp_json.json",
             "XML": "BASELINE/xml",
             "XML_TEMPLATES":"BASELINE/xml_templates",
             "XML_MODEL": "xml_params_mzmine.xml",
             "CANDIDATES": "BASELINE/candidates.csv",
             "PEAKTABLES": "BASELINE/peaktables",
             "MSMS": "BASELINE/msms"
         },"ADAP":
         {
             "SUMMARY":"ADAP/ADAP_parameters.csv",
             "SUMMARY_TEMPLATES": "ADAP/summary_templates.csv",
             "JSON":"ADAP/temp_json.json",
             "XML": "ADAP/xml",
             "XML_TEMPLATES":"ADAP/xml_templates",
             "XML_MODEL": "xml_params_mzmine.xml",
             "CANDIDATES": "ADAP/candidates.csv",
             "PEAKTABLES": "ADAP/peaktables",
             "MSMS": "ADAP/msms"
         },"OPENMS":{
        "PEAKTABLES": "OPENMS/peaktables",
        "MSMS": "OPENMS/msms"
    },"CENTWAVE":{
        "PEAKTABLES": "CENTWAVE/peaktables",
        "MSMS": "CENTWAVE/msms"
    },
     "FIGURES":{"RT_DEV":"figures/rt_dev_",
                "PEAKS":"figures/peaks",
                "DIAGNOSIS":"figures/diagnosis"},
     "DATAMATRIX": "temp",
     "FUSED_MSMS": "spectra_",
     "ANNOTATION":"data_",
     "EVALUATION": "evaluations",
     "TARGET":{"RT":"targeted_rt","INT":"targeted_int"},
     "RES_EVALUATION":{
         "FILE":"results_evaluations.csv",
         "PARAM": "best_param.xml",
         "PEAKTABLES" : "best_data",
         "FIGURE":"figure_evaluation.pdf"
     },
     "EIC":"EICs/eics.hdf5",
     "INITIAL_PARAMETERS":"initial_parameters.txt",
     "PARAMETERS":"parameters.yaml",
     "PYTHON_SCRIPT":"processing.py",
     "INCORRECT_FILES":"invalidMZML.txt"
     }

TEMP={"DIR":"temp",
      "CONVERSION":"temp/temp_names.csv",
      "GROUPING":{
          "TEMP":"temp/temp_grouping",
          "BLOCKS":"temp/blocks",
          "ALIGNMENT":"temp/alignment.rds"},
      "MISSING":{"HDF5":"temp/gap_filling.hdf5"},
      "FUSING":{"TEMP1":"temp/temp_dm_1",
                "TEMP2":"temp/temp_dm_2",
                "HDF5":"temp/fusing_msms.hdf5"},
      "REPLICATES":"temp/replicates",
      "IONANNOTATION":{"FULL":"temp/adducts.csv",
                       "MAIN":"temp/main_adducts.csv"},
      "POSTPROCESSING":"temp/raw_files.txt",
      "ISOTOPES":"temp/isotopes.txt"
      }

TO_CLEAN = ["ADAP","OPENMS","CENTWAVE",
            os.path.dirname(OUT["FIGURES"]["PEAKS"]),OUT["DATAMATRIX"],
            TEMP["DIR"],OUT["DONE"]]


##This is the order of the optimization variables, they are optimized by batch.
ORDER_VARIABLES_PEAKPICKING = ["peakpicking__traces_construction__ppm",
                               "peakpicking__peaks_deconvolution__peak_width__const",
                               "peakpicking__peaks_deconvolution__peak_width__add",
                               "peakpicking__peaks_deconvolution__peak_width__fac",
                               "peakpicking__peaks_deconvolution__SN",
                               "peakpicking__traces_construction_min_scan",
                               "peakpicking__traces_construction__dmz",
                               "peakpicking__peaks_deconvolution__rt_wavelet__const",
                               "peakpicking__peaks_deconvolution__rt_wavelet__add",
                               "peakpicking__peaks_deconvolution__coefficient_area_threshold",
                               "peakpicking__traces_construction__num_outliers",
                               "peakpicking__traces_construction__num_outliers",
                               "peakpicking__peaks_deconvolution__noise_level",
                               "peakpicking__noise_level__ms1","peakpicking__noise_level__ms2"]

ORDER_VARIABLES_GROUPING = ["grouping__ppm","grouping__drt","grouping__dmz",
                            "grouping__alpha","grouping__num_references"]

DATA = {
    "OPTIMIZATION":{"BALANCED_POINTS":os.environ['SLAWROOT']+ "pylcmsprocessing/data/balanced_samples.pickle"},
    "IONANNOTATION":{
        "XCMS_MODEL" : os.environ['SLAWROOT']+ "pylcmsprocessing/data/xcms_raw_model.RData",
        "positive":{"FULL":os.environ['SLAWROOT']+ "pylcmsprocessing/data/adducts_pos.txt","MAIN":os.environ['SLAWROOT']+ "pylcmsprocessing/data/adducts_main_pos.txt"},
        "negative":{"FULL":os.environ['SLAWROOT']+ "pylcmsprocessing/data/adducts_neg.txt","MAIN":os.environ['SLAWROOT']+ "pylcmsprocessing/data/adducts_main_neg.txt"}
    },
    "ISOTOPES":os.environ['SLAWROOT']+ "pylcmsprocessing/data/isotopes.tsv",
    "IONMODE":["positive","negative"],
    "YAML": os.environ['SLAWROOT']+ "pylcmsprocessing/data/parameters_set.yaml"
}

def default_adducts_positive():
    lines = [line.rstrip('\n') for line in open(DATA["IONANNOTATION"]["positive"]["FULL"])]
    return(lines)

def default_adducts_main_positive():
    lines = [line.rstrip('\n') for line in open(DATA["IONANNOTATION"]["positive"]["MAIN"])]
    return(lines)

def default_adducts_negative():
    lines = [line.rstrip('\n') for line in open(DATA["IONANNOTATION"]["negative"]["FULL"])]
    return(lines)

def default_adducts_main_negative():
    lines = [line.rstrip('\n') for line in open(DATA["IONANNOTATION"]["negative"]["MAIN"])]
    return(lines)

MASS_SPEC={"Exactive":"parameters_set_exactive",
           "TOF":"parameters_set_tof"}

CONSTANT={"PEAKPICKING_TIMOUT":900}
