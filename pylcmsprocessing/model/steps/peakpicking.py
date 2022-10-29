from abc import ABC, abstractmethod
import os
import model.helper.inputbuilder as hi
import hashlib
import logging
from common.tools import get_platform,find_rscript

class PeakPicking(ABC):
    @abstractmethod
    def need_computing(self):
        """Check if the processing need to be done (eg the peaktable does not exist already)"""
        pass
    @abstractmethod
    def command_line_processing(self):
        """Return the CLI used by SLAW"""
        pass

    @abstractmethod
    def get_output(self):
        """Get the filename of the resulting peaktable."""
        pass

    def filtering_peaktable(self,filtering):
        """Filter the peaktable if it exists"""
        pass

###We eill test this one a billion time because we don t really need the rest
class PeakPickingMZmine(PeakPicking):
    ###At the class level we always consider the maximum number of parameters
    ###input is the batch and output is not useful for MZmine
    def __init__(self,dbrow,path_software):
        self.computing= not os.path.exists(dbrow[5])
        self.output = dbrow[5]
        if self.computing:
            self.software = path_software
            self.input = dbrow[3]
            self.step = dbrow[6]

    def need_computing(self):
        return self.computing

    def need_conversion(self):
        return self.step == 1

    def command_line_processing(self):
        platform = get_platform()
        script_launch = ""
        redirect = ""
        software_part = self.software
        if platform=="Windows":
            script_launch="startMZmine_Windows.bat"
        elif platform=="Linux":
            script_launch="startMZmine-Linux"
        elif platform=="Mac":
            script_launch="startMZmine-MacOSX.command"
        else:
            script_launch="startMZmine-Linux"
        software_part =os.path.join(software_part,script_launch)
        return software_part+" "+'"'+self.input+'"'

    def get_output(self):
        return self.output

def hash_list(in_list):
    e = locals()
    del e["a"]
    str_val = "|".join([str(ee) for ee in in_list])
    str_val = str_val.encode()
    m = hashlib.md5(str_val)
    return m.hexdigest()

class PeakPickingOpenMS(PeakPicking):
    def __init__(self,row,min_fwhm,max_fwhm,fwhm_fac,snt,ppm,min_int,max_outlier,min_points,quant):
        self.input = row[3]
        self.min_fwhm = min_fwhm
        self.max_fwhm = max_fwhm
        self.fwhm= min_fwhm+(max_fwhm-min_fwhm)*fwhm_fac
        self.snt = snt
        self.ppm = ppm
        self.min_int = min_int
        self.max_outlier = max_outlier
        self.min_points = min_points
        self.quant = quant
        self.output = row[5]

    def need_computing(self):
        return not os.path.isfile(self.output)

    def get_output(self):
        return self.output

    def command_line_processing(self):
        snt_str = "true"
        if self.snt<=0:
            snt_str = "false"
            self.snt=0
        cli =  " ".join(["FeatureFinderMetabo",'-in "'+self.input+'" -out "'+self.output+
                         '" -algorithm:common:chrom_peak_snr',str(self.snt),
                         "-algorithm:epd:masstrace_snr_filtering",snt_str,"-algorithm:epd:max_fwhm",
                         str(self.max_fwhm)," -algorithm:epd:min_fwhm",str(self.min_fwhm),
                         "-algorithm:common:chrom_fwhm",str(self.fwhm),
                         "-algorithm:ffm:use_smoothed_intensities true -algorithm:epd:width_filtering fixed -algorithm:mtd:mass_error_ppm",
                         str(self.ppm),"-algorithm:mtd:quant_method",self.quant,"-algorithm:mtd:trace_termination_criterion outlier -algorithm:mtd:trace_termination_outliers",
                         str(self.max_outlier),"-algorithm:mtd:min_trace_length",str(self.min_fwhm)])
        return cli

class PeakPickingXCMS(PeakPicking):
    def __init__(self,row,min_peakwidth,max_peakwidth,snt,ppm,min_int,min_points):
        self.input = row[3]
        self.min_peakwidth = min_peakwidth
        self.max_peakwidth = max_peakwidth
        self.snt = snt
        self.ppm = ppm
        self.min_int = min_int
        self.min_points = min_points
        self.point_prefilter = max(int(min_points/2),3)
        self.int_prefilter = min_int/2
        self.output = row[5]

    def need_computing(self):
        return not os.path.isfile(self.output)

    def get_output(self):
        return self.output

    def command_line_processing(self):
        pjoin = os.path.join(find_rscript(), "wrapper_xcms_peak_picking.R")
        cli_args = [os.environ["RscriptString"]," ",pjoin,'"'+self.input+'"','"'+self.output+'"',self.ppm,self.min_peakwidth,
                    self.max_peakwidth,self.snt,
                    self.point_prefilter,self.int_prefilter,self.min_int]
        cli_args = [str(arg) for arg in cli_args]
        cli =  " ".join(cli_args)
        return cli