import common.references as cr
import yaml
import os
import xml.etree.ElementTree as ET
import model.helper.output_handler as oh
import model.helper.parameters_handler as ph
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

    def parameters_exist(self):
        return os.path.isfile(self.path_yaml)

    def initialize_yaml_polarity(self, polarity, path_yaml=None):
        if path_yaml is None:
            path_yaml = self.path_yaml
        with open(path_yaml, 'r') as stream:
            raw_yaml = yaml.safe_load(stream)
        #if polarity == "positive": NEWLY IT FILLS ALL ADDUCTS FOR BOTH MODES. SPLITTING IS CONFUSING.
        if raw_yaml["ion_annotation"]["adducts_positive"]["value"] == "NONE":
            raw_yaml["ion_annotation"]["adducts_positive"]["value"] = cr.default_adducts_positive()
            raw_yaml["ion_annotation"]["main_adducts_positive"]["value"] = cr.default_adducts_main_positive()
        #if polarity == "positive": NEWLY IT FILLS ALL ADDUCTS FOR BOTH MODES. SPLITTING IS CONFUSING.else:
        if raw_yaml["ion_annotation"]["adducts_negative"]["value"] == "NONE":
            raw_yaml["ion_annotation"]["adducts_negative"]["value"] = cr.default_adducts_negative()
            raw_yaml["ion_annotation"]["main_adducts_negative"]["value"] = cr.default_adducts_main_negative()
        # We write the yaml in the output directory
        with open(path_yaml, 'w') as outfile:
            yaml.safe_dump(raw_yaml, outfile, default_flow_style=False, sort_keys=False)


    def generate_yaml_files(self,force=False):
        #We rewrite the yaml file to be completed by the user.
        with open(cr.DATA["YAML"], 'r') as stream:
            raw_yaml = yaml.safe_load(stream)
        #We write the yaml in the output directory
        if not os.path.isfile(self.path_yaml) or force:
            with open(self.path_yaml, 'w') as outfile:
                yaml.safe_dump(raw_yaml,outfile, default_flow_style=False, sort_keys=False)
            self.initialize_yaml_polarity(self.polarity,self.path_yaml)


    def generate_MZmine_XML(self,path_xml=None,algorithm="ADAP"):
        if algorithm not in cr.MZMINE_ALGORITHM:
            algorithm="ADAP"
        if path_xml is None:
            path_xml = self.output.getFile(cr.OUT[algorithm]["XML_MODEL"])
        pfh = ph.ParametersFileHandler(self.path_yaml)
        PATH_XML = os.path.join("data",cr.ALGORITHMS_TABLE[algorithm][2])
        tree = ET.parse(PATH_XML)
        root = tree.getroot()


        ###Centroidization noise level MS1
        root[1][2][0][0].text=str(pfh[('peakpicking','noise_level_ms1')]['value'])
        ###Centroidization noise level MS2
        root[2][2][0][0].text=str(pfh[('peakpicking','noise_level_ms2')]['value'])

        ###Mass traces construction
        root[3][3].text=str(pfh[('peakpicking','traces_construction','min_scan')]['value'])
        root[3][4].text=str(pfh[('peakpicking','noise_level_ms1')]['value'])
        root[3][5].text=str(pfh[('peakpicking','noise_level_ms1')]['value'])
        root[3][6][1].text=str(pfh[('peakpicking','traces_construction','ppm')]['value'])
        root[3][6][0].text=str(pfh[('peakpicking','traces_construction','dmz')]['value'])

        ###This is different depending of the algorithm
        ###Peak deconvolution
        if algorithm=="ADAP":
            root[4][2][5][0].text=str(pfh[('peakpicking','peaks_deconvolution','SN')]['value'])
            root[4][2][5][2].text=str(pfh[('peakpicking','peaks_deconvolution','noise_level')]['value'])
            root[4][2][5][3].text=str(pfh[('peakpicking','peaks_deconvolution','coefficient_area_threshold')]['value'])
            root[4][2][5][4][0].text=str(pfh[('peakpicking','peaks_deconvolution','peak_width')]['value'][0])
            root[4][2][5][4][1].text=str(pfh[('peakpicking','peaks_deconvolution','peak_width')]['value'][1])
            root[4][2][5][5][0].text=str(pfh[('peakpicking','peaks_deconvolution','rt_wavelet')]['value'][0])
            root[4][2][5][5][1].text=str(pfh[('peakpicking','peaks_deconvolution','rt_wavelet')]['value'][1])
            root[4][4].text=str(pfh[('peakpicking','peaks_deconvolution','ms2_mz_tol')]['value'])
            root[4][5].text=str(pfh[('peakpicking','peaks_deconvolution','ms2_rt_tol')]['value'])
        if algorithm=="BASELINE":
            root[4][2][0][0].text=str(pfh[('peakpicking','peaks_deconvolution','noise_level')]['value'])
            root[4][2][0][1][0].text=str(pfh[('peakpicking','peaks_deconvolution','peak_width')]['value'][0])
            root[4][2][0][1][1].text=str(pfh[('peakpicking','peaks_deconvolution','peak_width')]['value'][1])
            root[4][2][0][2].text=str(pfh[('peakpicking','noise_level_ms1')]['value'])

        ###We write the XML file somewhere
        tree.write(path_xml)

    def openYamlParameters(self):
        raw_yaml = None
        with open(self.path_yaml, 'r') as stream:
            raw_yaml = yaml.safe_load(stream)
        return raw_yaml

    def availableCpus(self):
        return multiprocessing.cpu_count()-1

