import common.tools as ct
import os


class evaluationComparator():
    def __init__(self,path_db,best_param,figure,peaktables):
        self.db = path_db
        self.param = best_param
        self.figure = figure
        self.peaktable = peaktables

    def command_line(self):
        pscript = ct.find_rscript()
        command_line = os.path.join(pscript,"evaluationSummary.R")
        ####We give all the name of the grouping parameters implicated in a single file
        return "Rscript "+command_line+" "+self.db+" "+self.figure+" "+self.param+" "+self.peaktable
