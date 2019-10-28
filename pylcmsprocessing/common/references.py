ALGORITHMS_TABLE = {
    "ADAP": ["MZmine", "adap.json", "/pylcmsprocessing/data/batch_adap_2_52.xml"],
    "SAVGOL": ["MZmine", "savgol.json", "/pylcmsprocessing/data/batch_savgol.xml"],
    "FeatureFinderMetabo": ["openMS", "openms.json", None]
}

SUMMARY_DIR = "summary_parameters"
SUMMARY_COMBINATIONS = "summary_combinations"

OUT={"DB":"processing.sqlite",
    "ADAP":
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
     },"CENTWAVE":
          {"SUMMARY":"ADAP/ADAP_parameters.csv",
    "SUMMARY_TEMPLATES": "ADAP/summary_templates.csv",
    "JSON":"ADAP/temp_json.json",
    "XML": "ADAP/xml",
    "XML_TEMPLATES":"ADAP/xml_templates",
    "CANDIDATES": "ADAP/candidates.csv",
    "PEAKTABLES": "ADAP/peaktables",
     },
     "DATAMATRIX": "datamatrices",
     "ANNOTATION":"annotated_peaktable_",
     "EVALUATION": "evaluations",
     "RES_EVALUATION":{
         "FILE":"results_evaluations.csv",
         "PARAM": "best_param.xml",
         "PEAKTABLES" : "best_data",
        "FIGURE":"figure_evaluation.pdf"
     },
     "PARAMETERS":"parameters.yaml",
     "PYTHON_SCRIPT":"processing.py"
}


DATA = {
    "IONANNOTATION":{
        "XCMS_MODEL" : "/pylcmsprocessing/data/xcms_raw_model.RData",
        "positive":{"FULL":"/pylcmsprocessing/data/adducts_pos.txt","MAIN":"/pylcmsprocessing/data/adducts_main_pos.txt"},
        "negative":{"FULL":"/pylcmsprocessing/data/adducts_neg.txt","MAIN":"/pylcmsprocessing/data/adducts_main_neg.txt"}
    },
    "IONMODE":["positive","negative"],
    "YAML": "/pylcmsprocessing/data/parameters_set.yaml"
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
    lines = [line.rstrip('\n') for line in open(DATA["IONANNOTATION"]["positive"]["MAIN"])]
    return(lines)


TEMP={"CONVERSION":
      "temp/temp_names.csv",
      "GROUPING":
      "temp/temp_grouping",
      "REPLICATES":"temp/replicates",
      "IONANNOTATION":{"FULL":"temp/adducts.csv",
                       "MAIN":"temp/main_adducts.csv"}
      }


MASS_SPEC={"Exactive":"parameters_set_exactive",
"TOF":"parameters_set_tof"}
