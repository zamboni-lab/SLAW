import itertools
import random
import csv
import json

class ParametersGenerator:
    def __init__(self,json, max_parameters= 50, seed = 512):
        self.json = json
        self.seed = seed
        self.max_parameters=max_parameters
        self.combinations = []

    def read_json(self):
        with open(self.json) as json_file:
            return json.load(json_file)


    def get_list(self):
        fjson = self.read_json()
        return list(fjson.values())

    def get_names(self):
        fjson = self.read_json()
        return list(fjson.keys())

    def make_combination(self):
        all_args = self.get_list()
        self.combinations = list(itertools.product(*all_args))


    def sample_combinations(self):
        if len(self.combinations)==0:
            raise Exception("Calculate the combination before sampling it")
        elif len(self.combinations) > self.max_parameters:
            random.seed(self.seed)
            idx = random.sample(range(len(self.combinations)),self.max_parameters)
            print("sampling ",self.max_parameters)
            self.combinations=[self.combinations[x] for x in idx]


    def export_combinations(self,outfile):
        with open(outfile, 'w', newline="") as csvFile:
            writer = csv.writer(csvFile,delimiter=";")
            writer.writerow(self.get_names())
            writer.writerows(self.combinations)
        csvFile.close()
        print("The ",len(self.combinations)," parameters has been written in "+outfile)

