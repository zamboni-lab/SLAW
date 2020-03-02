import common.tools
import os
import common.references
import subprocess
import pandas as pd
import sqlite3


class inputBuilder:
    def __init__(self):
        pass

    def build_inputs(self):
        pass

    def filter_inputs(self):
        pass

    def clean(self):
        pass


###build the nput and return a
class MZMineBuilder(inputBuilder):
    def __init__(self, tuple, db, output, id):
        self.xml = tuple[4]
        self.json = tuple[3]
        self.summary = tuple[5]
        self.db = db
        self.output = output
        self.id = id
        super(inputBuilder, self).__init__()

    def build_inputs(self):
        # OUT = {"ADAP":
        #            {"SUMMARY": "ADAP/ADAP_parameters.csv",
        #             "SUMMARY_TEMPLATES": "ADAP/candidates.csv",
        #             "XML": "ADAP/xml",
        #             "XML_TEMPLATES": "ADAP/xml_templates",
        #             "CANDIDATES": "ADAP/candidates.csv"}
        #        }

        ###We create the necessary directories
        xml_templates = self.output.getDir(common.references.OUT["ADAP"]["XML_TEMPLATES"])
        summary_templates = self.output.getFile(common.references.OUT["ADAP"]["SUMMARY_TEMPLATES"])
        xml = self.output.getDir(common.references.OUT["ADAP"]["XML"])
        candidates = self.output.getFile(common.references.OUT["ADAP"]["CANDIDATES"])
        peaktables = self.output.getDir(common.references.OUT["ADAP"]["PEAKTABLES"])
        msms = self.output.getDir(common.references.OUT["ADAP"]["MSMS"])
        prscript = common.tools.find_rscript()

        scriptMZmine = os.path.join(prscript, "wrapper_MZmine_peak_picking.R")

        commandline = " ".join(
            [scriptMZmine, self.db, self.json, self.xml, xml, self.summary, summary_templates, xml_templates,
             candidates, peaktables, msms, str(self.id)])

        print("Computing MZmine parameters")
        subprocess.call("Rscript " + commandline, shell=True)

    def build_inputs_single_parameter_set(self,xml_file):

        ###We create the necessary directories
        pxml = self.output.getDir(common.references.OUT["ADAP"]["XML"])
        candidates = self.output.getFile(common.references.OUT["ADAP"]["CANDIDATES"])
        peaktables = self.output.getDir(common.references.OUT["ADAP"]["PEAKTABLES"])
        msms = self.output.getDir(common.references.OUT["ADAP"]["MSMS"])
        prscript = common.tools.find_rscript()
        tjson = self.output.getFile(common.references.OUT["ADAP"]["JSON"])
        scriptMZmine = os.path.join(prscript, "wrapper_MZmine_peak_picking_xml.R")
        commandline = " ".join([scriptMZmine, self.db, xml_file, pxml, candidates, peaktables, msms, tjson])

        print("Computing MZmine parameters")
        subprocess.call("Rscript " + commandline, shell=True)
        hash_val = common.tools.md5sum(xml_file)

        conn = sqlite3.connect(self.db)
        c = conn.cursor()

        try:
            ctuple = (1, 1, xml_file, hash_val, "", "", "", "", "", "", "", "")
            c.execute("""INSERT INTO peakpicking VALUES (?,?,?,?,?,?,?,?,?,?,?, ?)""", ctuple)
        except sqlite3.IntegrityError as e:
            pass

        ###We read the candidates table output by the R script.
        path_candidates = self.output.getFile(common.references.OUT["ADAP"]["CANDIDATES"])
        candidates = pd.read_csv(path_candidates,sep=";")

        pid = 1
        counter_processed = 0
        counter_to_process = 0

        for row in candidates.itertuples(index=False):
            ###If the algorithm does not exist we skip
            ###We try to insert it
            try:
                ctuple = (pid, 1, int(row[1]), row[2], row[3], row[4], row[5], 1)
                pid += 1
                c.execute("""INSERT INTO processing VALUES (?,?,?,?,?,?,?,?)""", ctuple)
                counter_to_process += 1
            except sqlite3.IntegrityError as e:
                counter_processed += 1

        print(counter_processed, " existing peak picking ", counter_to_process, " added.")
        conn.commit()
        conn.close()



    def get_first_id(self, c, table):
        c.execute("SELECT MAX(id) FROM " + table)
        vid = c.fetchone()
        if vid[0] is None:
            vid = 1
        else:
            vid = vid + 1
        return vid

    def filter_inputs(self):
        conn = sqlite3.connect(self.db)
        c = conn.cursor()

        cid = self.get_first_id(c, "peakpicking")

        ###We check if the output already exists. If they do we romve the files
        summary_templates = self.output.getFile(common.references.OUT["ADAP"]["SUMMARY_TEMPLATES"])
        print("opening ",summary_templates)
        peakpickings = pd.read_csv(summary_templates,sep=";")

        ###We first check if th has of the algorithm exist, if they do
        existing_processing = set([])

        # c.execute('''CREATE TABLE IF NOT EXISTS peakpicking
        # (id INTEGER PRIMARY KEY,
        # algorithm INTEGER,
        # parameter TEXT,
        # hash TEXT,
        # peaktable TEXT,
        # hash_peaktable TEXT,
        # evaluation TEXT,
        # hash_evaluation TEXT,
        # CONSTRAINT known_software,
        # FOREIGN KEY (algorithm)
        # REFERENCES algorithms(id)),
        # CONSTRAINT hashing UNIQUE (hash)''')

        # We create a vector storing the new index to process
        new_index = [0] * 1000
        counter_index = 0
        for row in peakpickings.itertuples(index=False):
            ###We try to insert if if it is not exisiting
            try:
                ctuple = (cid, self.id, row[0], row[2], "", "", "", "")
                new_index[counter_index] = cid
                cid = cid + 1
                c.execute("""INSERT INTO peakpicking VALUES (?,?,?,?,?,?,?,?)""", ctuple)
            except sqlite3.IntegrityError as e:
                existing_processing.add(row[0])
                new_index[counter_index] = None
            counter_index+=1

        new_index = new_index[0:counter_index]

        ###Adding the required elements to the processing table

        # c.execute('''CREATE TABLE IF NOT EXISTS processing
        #   (id INTEGER PRIMARY KEY,
        #   peakpicking INTEGER,
        #   sample INTEGER,
        #   input TEXT NOT NULL,
        #   hash_input TEXT NOT NULL,
        #   output TEXT NOT NULL,
        #   hash_output TEXT NOT NULL,
        #   step INTEGER NOT NULL,
        #   CONSTRAINT fk_peakpicking
        #   FOREIGN KEY (peakpicking)
        #   REFERENCES peakpicking(id),
        #   CONSTRAINT fk_sample
        #   FOREIGN KEY (sample)
        #   REFERENCES samples(id),
        #   CONSTRAINT hashing_input UNIQUE (hash_input),
        #   CONSTRAINT hashing_output UNIQUE (hash_output)
        #   )''')

        path_candidates = self.output.getFile(common.references.OUT["ADAP"]["CANDIDATES"])
        candidates = pd.read_csv(path_candidates,sep=";")

        pid = self.get_first_id(c, "processing")
        counter_processed = 0
        counter_to_process = 0

        for row in candidates.itertuples(index=False):
            ###If the algorithm does not exist we skip
            print(row)
            if row[0] in existing_processing:
                continue
            ###We try to insert it
            try:
                # print("_",row[0],"_",len(row))
                ctuple = (pid, new_index[row[0] - 1], row[1], row[2], row[3], row[4], 1)
                pid += 1
                c.execute("""INSERT INTO processing VALUES (?,?,?,?,?,?,?)""", ctuple)
                counter_to_process += 1
            except sqlite3.IntegrityError as e:
                counter_processed += 1

        print(counter_processed, " existing peak picking ", counter_to_process, " added.")

        conn.commit()
        conn.close()

        # iter_csv = pd.read_csv(candidates, iterator=True, chunksize=500)
        # to_remove = pd.concat([chunk[self.output.exists(chunk['output'])] for chunk in iter_csv])

    def clean(self):
        summary_templates = self.output.getFile(common.references.OUT["ADAP"]["SUMMARY_TEMPLATES"])
        os.remove(summary_templates)
        path_candidates = self.output.getFile(common.references.OUT["ADAP"]["CANDIDATES"])
        os.remove(path_candidates)
