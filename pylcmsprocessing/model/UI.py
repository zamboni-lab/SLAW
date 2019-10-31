import common.references as cr
import yaml
import os
import xml.etree.ElementTree as ET
import model.output_handler as oh
import multiprocessing

####This is the set of functions whici are furnished to the Users
class UI:
    def __init__(self,output_dir,raw_files,polarity="positive",mass_spec = "Exactive",path_yaml = None, num_workers=None):
        output_dir = os.path.normpath(output_dir)
        self.path_rawfiles = raw_files
        self.output = oh.OutputDirectory(output_dir)
        if path_yaml is None:
            self.path_yaml = self.output.getFile(cr.OUT["PARAMETERS"])
        else:
            self.path_yaml = path_yaml
        if num_workers is None:
            self.num_workers = self.availableCpus()
        elif num_workers>self.availableCpus():
            self.num_workers = self.availableCpus()
        else:
            self.num_workers=num_workers
        self.polarity=polarity
        if not mass_spec in cr.MASS_SPEC:
            raise Exception("Unknown mass spectrometer.")


    def generate_yaml_files(self,force=False):
        #We rewrite the yaml file to be completed by the user.
        if not force:
            raw_yaml = None
            with open(cr.DATA["YAML"], 'r') as stream:
                raw_yaml = yaml.safe_load(stream)
            raw_yaml["ion_annotation"]["polarity"]['value'] = self.polarity
            if self.polarity=="positive":
                raw_yaml["ion_annotation"]["adducts_positive"]["value"] = cr.default_adducts_positive()
                raw_yaml["ion_annotation"]["main_adducts_positive"]["value"] = cr.default_adducts_main_positive()
            else:
                raw_yaml["ion_annotation"]["adducts_negative"]["value"] = cr.default_adducts_negative()
                raw_yaml["ion_annotation"]["main_adducts_negative"]["value"] = cr.default_adducts_main_negative()
            #We write the yaml in the output directory
            with open(self.path_yaml, 'w') as outfile:
                yaml.dump(raw_yaml,outfile, default_flow_style=False)

    def generate_MZmine_XML(self,path_xml=None):
        if path_xml is None:
            path_xml = self.output.getFile(cr.OUT["ADAP"]["XML_MODEL"])
        raw_yaml = None
        with open(self.path_yaml, 'r') as stream:
            raw_yaml = yaml.safe_load(stream)
        PATH_XML = os.path.join("data",cr.ALGORITHMS_TABLE["ADAP"][2])
        tree = ET.parse(PATH_XML)
        root = tree.getroot()

        ###Centroidization noise level MS1
        root[1][2][0][0].text=str(raw_yaml['peakpicking']['centroidization']['noise_level_ms1']['value'])

        ###Centroidization noise level MS2
        root[2][2][0][0].text=str(raw_yaml['peakpicking']['centroidization']['noise_level_ms1']['value'])

        ###Mass traces construction
        root[3][3].text=str(raw_yaml['peakpicking']['traces_construction']['min_scan']['value'])
        root[3][4].text=str(raw_yaml['peakpicking']['centroidization']['noise_level_ms1']['value'])
        root[3][5].text=str(raw_yaml['peakpicking']['centroidization']['noise_level_ms1']['value'])
        root[3][6][1].text=str(raw_yaml['peakpicking']['traces_construction']['ppm']['value'])
        root[3][6][0].text=str(raw_yaml['peakpicking']['traces_construction']['dmz']['value'])

        ###Peak deconvolution
        root[4][2][5][0].text=str(raw_yaml['peakpicking']['peaks_deconvolution']['SN']['value'])
        root[4][2][5][2].text=str(raw_yaml['peakpicking']['peaks_deconvolution']['noise_level']['value'])
        root[4][2][5][4][0].text=str(raw_yaml['peakpicking']['peaks_deconvolution']['peak_width_min']['value'])
        root[4][2][5][4][1].text=str(raw_yaml['peakpicking']['peaks_deconvolution']['peak_width_max']['value'])
        root[4][2][5][5][0].text=str(raw_yaml['peakpicking']['peaks_deconvolution']['rt_wavelet_min']['value'])
        root[4][2][5][5][0].text=str(raw_yaml['peakpicking']['peaks_deconvolution']['rt_wavelet_max']['value'])
        root[4][4].text=str(raw_yaml['peakpicking']['peaks_deconvolution']['ms2_mz_tol']['value'])
        root[4][5].text=str(raw_yaml['peakpicking']['peaks_deconvolution']['ms2_rt_tol']['value'])

        ###We write the XML file somewhere
        tree.write(path_xml)


    def openYamlParameters(self):
        raw_yaml = None
        with open(self.path_yaml, 'r') as stream:
            raw_yaml = yaml.safe_load(stream)
        return raw_yaml

    def availableCpus(self):
        return multiprocessing.cpu_count()-1

    ####Generate the python file used to process the data eventually
    def generate_python_file(self,path_python = None):
        #We get the dayabase storage stuff eventually
        path_db = self.output.getFile(cr.OUT["DB"])
        raw_yaml = self.openYamlParameters()
        if path_python is None:
            path_python = self.output.getFile(cr.OUT["PYTHON_SCRIPT"])
        if self.num_workers is None:
            self.num_workers=self.availableCpus()
        vfile = ('from model.experiment import Experiment\n'
                 '\n'
                 'if __name__ == "__main__":\n'
                 '    exp = Experiment("') + path_db + '''",reset=True)
    exp.initialise_database(''' + str(multiprocessing.cpu_count()-1) +''',
    "''' + os.path.normpath(self.output.getRoot())+'''",
    "''' + self.polarity+'''",
    "''' + os.path.normpath(self.path_rawfiles) +'''",
    ["ADAP"], 1)
    exp.building_inputs_single_processing(
        "''' + self.path_xml +'''")
    ####The path to MZmine is always input in full
    exp.run("''' + str(raw_yaml["constant"]["path_mzmine"]["value"]) + '''", ''' + str(self.availableCpus()) + ''')
    exp.correct_conversion()
    exp.group(max_workers=1,mztol='''+str(raw_yaml["grouping"]["dmz"]["value"])+\
    ',rttol='+str(raw_yaml["grouping"]["drt"]["value"])+',intensity="'+\
    str(raw_yaml["grouping"]["extracted_quantity"]["value"])+'''")
    exp.annotate_ions('''+str(raw_yaml["ion_annotation"]["num_files"]["value"])+''','''+str(raw_yaml["ion_annotation"]["ppm"]["value"]) + ','+ \
        str(raw_yaml["ion_annotation"]["dmz"]["value"]) +',min_filter='+str(raw_yaml["ion_annotation"]["min_filter"]["value"])+\
                ',adducts=["' + '","'.join(raw_yaml["ion_annotation"]["adducts"]["value"]) + '"],main_adducts=["'+ \
        '","'.join(raw_yaml["ion_annotation"]["main_adducts"]["value"])+'"], max_workers='+str(self.num_workers)+")"
        with open(path_python,"w") as py_file:
            py_file.write(vfile)
