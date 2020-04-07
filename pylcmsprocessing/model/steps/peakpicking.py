from abc import ABC, abstractmethod
import os
from common.tools import get_platform,find_rscript


class PeakPicking(ABC):

    @abstractmethod
    def need_computing(self):
        pass

    @abstractmethod
    def command_line_processing(self):
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

    def command_line_processing(self,hide=True):
        platform = get_platform()
        script_launch = ""
        redirect = ""
        software_part = self.software
        if platform=="Windows":
            script_launch="startMZmine_Windows.bat"
            if hide:
                redirect=" > NUL 2>&1"
        elif platform=="Linux":
            script_launch="startMZmine-Linux"
            if hide:
                redirect=" &>/dev/null"
        elif platform=="Mac":
            script_launch="startMZmine-MacOSX.command"
            if hide:
                redirect=" > NUL 2>&1"
        else:
            script_launch="startMZmine-Linux"
            if hide:
                redirect=" &>/dev/null"


        software_part =os.path.join(software_part,script_launch)

        return software_part+" "+self.input
        return software_part+" "+self.input+redirect

    def get_output(self):
        return self.output

# class PeakPickingXCMS(PeakPicking):
