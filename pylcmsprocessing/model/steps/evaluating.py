import common.tools as ct
import common.references as cr
import os

class Evaluate:

    ####Row is a row of the processing database and peak
    def __init__(self,datamatrix,output,replicates):
        self.input = datamatrix
        self.output = output
        self.replicates = replicates

    def need_computing(self):
        return not os.path.exists(self.output)


    def get_output(self):
        return self.output

    def get_input(self):
        return self.input

    def command_line(self):
        pscript = ct.find_rscript()
        command_line = os.path.join(pscript,"evaluationFromPeakTable.R")
        ####We give all the name of the grouping parameters implicated in a single file
        return "Rscript "+command_line+" "+self.input+" "+self.replicates+" "+self.output
