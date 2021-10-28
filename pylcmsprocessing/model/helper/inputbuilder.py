import common.tools
import os
import common.references as cr
import subprocess
import pandas as pd
import hashlib
import sqlite3
import logging
import shutil

# The builders files are just here to do


###Given a hash conrresponding to an experiemnet, handle the removal of the data if needed.
class inputChecker:
    def __init__(self,output_handler,hash,path_db):
        ###Given an
        self.hash = hash
        self.output = output_handler
        self.db = path_db

    ## We chekc fi the database need reset
    def need_reset(self):
        conn = sqlite3.connect(self.db)
        queries = conn.execute("SELECT * FROM peakpicking").fetchall()
        conn.close()
        if len(queries)>=2:
            return True
        hashes = [q[3] for q in queries]
        return not(self.hash in hashes)

    ## If needed we rester the processing table and the hash table
    def prepare_db(self):
        if self.need_reset():
            conn = sqlite3.connect(self.db)
            ###We empty the processing and peakpicking table
            conn.execute("DELETE FROM processing")
            conn.execute("DELETE FROM peakpicking")
            conn.commit()
            conn.close()

    ###If it is over we clean the directory
    def prepare_directories(self):
        if self.need_reset():
            for cl in cr.TO_CLEAN:
                pcl = self.output.getPath(cl)
                if os.path.isdir(pcl):
                    shutil.rmtree(pcl)
                elif os.path.isfile(pcl):
                    os.remove(pcl)

    def prepare_outputs(self):
        self.prepare_db()
        self.prepare_directories()


class inputBuilder:
    def __init__(self):
        pass

    def build_inputs(self):
        pass

    def filter_inputs(self):
        pass

    def clean(self):
        pass

    def get_first_id(self, c, table):
        c.execute("SELECT MAX(id) FROM " + table)
        vid = c.fetchone()
        if vid[0] is None:
            vid = 1
        else:
            vid = vid[0] + 1
        return vid


###build the nput and return a
class MZMineBuilder(inputBuilder):
    def __init__(self, tuple, db, output, id, algorithm = "ADAP"):
        self.xml = tuple[4]
        self.json = tuple[3]
        self.summary = tuple[5]
        self.algorithm = algorithm
        self.db = db
        self.output = output
        self.id = id
        super(inputBuilder, self).__init__()

    def build_inputs(self):
        ###We create the necessary directories
        xml_templates = self.output.getDir(cr.OUT["ADAP"]["XML_TEMPLATES"])
        summary_templates = self.output.getFile(cr.OUT["ADAP"]["SUMMARY_TEMPLATES"])
        xml = self.output.getDir(cr.OUT["ADAP"]["XML"])
        candidates = self.output.getFile(cr.OUT["ADAP"]["CANDIDATES"])
        peaktables = self.output.getDir(cr.OUT["ADAP"]["PEAKTABLES"])
        msms = self.output.getDir(cr.OUT["ADAP"]["MSMS"])
        prscript = common.tools.find_rscript()

        scriptMZmine = os.path.join(prscript, "wrapper_MZmine_peak_picking.R")

        commandline = " ".join(
            [scriptMZmine, self.db, self.json, self.xml, xml, self.summary, summary_templates, xml_templates,
             candidates, peaktables, msms, str(self.id)])

        logging.info("Computing MZmine parameters")
        subprocess.call("Rscript " + commandline, shell=True)

    def build_inputs_single_parameter_set_optimization(self,xml_file):
        hash_val = common.tools.md5sum(xml_file)
        ###Checking if input folder exists.
        icheck = inputChecker(self.output, hash_val, self.db)
        icheck.prepare_directories()
        icheck.prepare_db()

        ###In this case we ust have to remap the parameters
        ###We create the necessary directories
        pxml = self.output.getDir(cr.OUT[self.algorithm]["XML"])
        candidates = self.output.getFile(cr.OUT[self.algorithm]["CANDIDATES"])
        peaktables = self.output.getDir(cr.OUT[self.algorithm]["PEAKTABLES"])
        msms = self.output.getDir(cr.OUT[self.algorithm]["MSMS"])
        prscript = common.tools.find_rscript()
        tjson = self.output.getFile(cr.OUT[self.algorithm]["JSON"])
        scriptMZmine = os.path.join(prscript, "wrapper_MZmine_peak_picking_xml.R")
        commandline = " ".join([scriptMZmine, self.db, xml_file, pxml, candidates, peaktables, msms, tjson])
        # logging.debug(commandline)
        logging.info("Computing MZmine parameters")
        subprocess.call("Rscript " + commandline, shell=True)
        conn = sqlite3.connect(self.db)
        c = conn.cursor()
        try:
            ctuple = (1, 1, xml_file, hash_val, "", "", "", "")
            c.execute("""INSERT INTO peakpicking VALUES (?,?,?,?,?,?,?,?)""", ctuple)
        except sqlite3.IntegrityError as e:
            pass

        ###We read the candidates table output by the R script.
        path_candidates = self.output.getFile(cr.OUT[self.algorithm]["CANDIDATES"])
        candidates = pd.read_csv(path_candidates,sep=";")

        pid = 1
        counter_processed = 0
        counter_to_process = 0

        for row in candidates.itertuples(index=False):
            ###If the algorithm does not exist we skip
            ###We try to insert it
            try:
                ctuple = (pid, 1, int(row[1]), row[2], row[3], row[4], row[5], 1, 1)
                pid += 1
                c.execute("""INSERT INTO processing VALUES (?,?,?,?,?,?,?,?,?)""", ctuple)
                counter_to_process += 1
            except sqlite3.IntegrityError as e:
                counter_processed += 1

        logging.info(str(counter_processed)+" existing peak picking "+ str(counter_to_process)+ " added.")
        conn.commit()
        conn.close()

    def build_inputs_single_parameter_set(self,xml_file):
        hash_val = common.tools.md5sum(xml_file)
        ###Checking if input folder exists.
        icheck = inputChecker(self.output, hash_val, self.db)
        icheck.prepare_directories()
        icheck.prepare_db()

        ###We create the necessary directories
        pxml = self.output.getDir(cr.OUT[self.algorithm]["XML"])
        candidates = self.output.getFile(cr.OUT[self.algorithm]["CANDIDATES"])
        peaktables = self.output.getDir(cr.OUT[self.algorithm]["PEAKTABLES"])
        msms = self.output.getDir(cr.OUT[self.algorithm]["MSMS"])
        prscript = common.tools.find_rscript()
        tjson = self.output.getFile(cr.OUT[self.algorithm]["JSON"])
        scriptMZmine = os.path.join(prscript, "wrapper_MZmine_peak_picking_xml.R")
        commandline = " ".join([scriptMZmine, self.db, xml_file, pxml, candidates, peaktables, msms, tjson])
        logging.info("Computing MZmine parameters")
        subprocess.call("Rscript " + commandline, shell=True)
        conn = sqlite3.connect(self.db)
        c = conn.cursor()
        try:
            ctuple = (1, 1, xml_file, hash_val, "", "", "", "")
            c.execute("""INSERT INTO peakpicking VALUES (?,?,?,?,?,?,?,?)""", ctuple)
        except sqlite3.IntegrityError as e:
            pass

        ###We read the candidates table output by the R script.
        path_candidates = self.output.getFile(cr.OUT[self.algorithm]["CANDIDATES"])
        candidates = pd.read_csv(path_candidates,sep=";")

        pid = 1
        counter_processed = 0
        counter_to_process = 0

        for row in candidates.itertuples(index=False):
            ###If the algorithm does not exist we skip
            ###We try to insert it
            try:
                ctuple = (pid, 1, int(row[1]), row[2], row[3], row[4], row[5], 1, 1)
                pid += 1
                c.execute("""INSERT INTO processing VALUES (?,?,?,?,?,?,?,?,?)""", ctuple)
                counter_to_process += 1
            except sqlite3.IntegrityError as e:
                counter_processed += 1

        logging.info(str(counter_processed)+" existing peak picking "+ str(counter_to_process)+ " added.")
        conn.commit()
        conn.close()


    def filter_inputs(self):
        conn = sqlite3.connect(self.db)
        c = conn.cursor()

        cid = self.get_first_id(c, "peakpicking")

        ###We check if the output already exists. If they do we romve the files
        summary_templates = self.output.getFile(cr.OUT[self.algorithm]["SUMMARY_TEMPLATES"])
        logging.info("Opening ",summary_templates)
        peakpickings = pd.read_csv(summary_templates,sep=";")

        ###We first check if th has of the algorithm exist, if they do
        existing_processing = set([])
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
        path_candidates = self.output.getFile(common.references.OUT[self.algorithm]["CANDIDATES"])
        candidates = pd.read_csv(path_candidates,sep=";")

        pid = self.get_first_id(c, "processing")
        counter_processed = 0
        counter_to_process = 0

        for row in candidates.itertuples(index=False):
            ###If the algorithm does not exist we skip
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
        logging.info(str(counter_processed)+" existing peak picking "+str(counter_to_process)+" added.")
        conn.commit()
        conn.close()

        # iter_csv = pd.read_csv(candidates, iterator=True, chunksize=500)
        # to_remove = pd.concat([chunk[self.output.exists(chunk['output'])] for chunk in iter_csv])

    def clean(self):
        summary_templates = self.output.getFile(cr.OUT[self.algorithm]["SUMMARY_TEMPLATES"])
        os.remove(summary_templates)
        path_candidates = self.output.getFile(cr.OUT[self.algorithm]["CANDIDATES"])
        os.remove(path_candidates)

def hash_list(in_list):
    str_val = "|".join([str(ee) for ee in in_list])
    str_val = str_val.encode()
    m = hashlib.md5(str_val)
    return m.hexdigest()

class openMSBuilder(inputBuilder):
    def __init__(self, db, output):
        self.db = db
        self.output = output

    def build_inputs_single_parameter_set(self,*args):

        conn = sqlite3.connect(self.db)
        c = conn.cursor()
        cid = self.get_first_id(c, "peakpicking")

        ###We create the necessary directories
        conn = sqlite3.connect(self.db)
        c = conn.cursor()
        hash_val = hash_list(args)
        icheck = inputChecker(self.output, hash_val, self.db)
        icheck.prepare_directories()
        icheck.prepare_db()
        path_peaktables = self.output.getDir(cr.OUT["OPENMS"]["PEAKTABLES"])
        path_msms = self.output.getDir(cr.OUT["OPENMS"]["MSMS"])
        try:
            ctuple = (cid, 2, "",hash_val, "", "", "", "")
            c.execute("""INSERT INTO peakpicking VALUES (?,?,?,?,?,?,?,?)""", ctuple)
        except sqlite3.IntegrityError as e:
            pass

        ###We just insert the sample name
        samples = c.execute("SELECT id,path,level FROM samples").fetchall()
        pid = 0
        counter_to_process = 0
        counter_processed = 0
        ###We create the processing table and add the input if they don t exists.
        for vid,path,level in samples:
            temp_l = list(args)+[path]
            hash_sample = hash_list(temp_l)
            path_ms = os.path.join(path_peaktables,"peaktable_"+str(vid)+"_"+hash_sample+".featureXML")
            if level=="MS2":
                path_ms = "NOT PROCESSED"
            else :
                pid += 1
            path_msms_s = os.path.join(path_msms,"msms_"+str(vid)+"_"+hash_sample+".mgf")
            try:
                ctuple = (pid, 1, vid, path, hash_sample, path_ms, path_msms_s, 1, 1)
                c.execute("""INSERT INTO processing VALUES (?,?,?,?,?,?,?,?,?)""", ctuple)
                counter_to_process += 1
            except sqlite3.IntegrityError as e:
                counter_processed += 1

        logging.info(str(counter_to_process)+" peakpicking added")
        conn.commit()
        conn.close()

class xcmsBuilder(inputBuilder):
    def __init__(self, db, output):
        self.db = db
        self.output = output

    def build_inputs_single_parameter_set(self,*args):
        conn = sqlite3.connect(self.db)
        c = conn.cursor()
        cid = self.get_first_id(c, "peakpicking")
        ###We create the necessary directories
        conn = sqlite3.connect(self.db)
        c = conn.cursor()
        hash_val = hash_list(args)
        icheck = inputChecker(self.output, hash_val, self.db)
        icheck.prepare_directories()
        icheck.prepare_db()
        path_peaktables = self.output.getDir(cr.OUT["CENTWAVE"]["PEAKTABLES"])
        path_msms = self.output.getDir(cr.OUT["CENTWAVE"]["MSMS"])
        try:
            ctuple = (cid, 3, "", hash_val, "", "", "", "")
            c.execute("""INSERT INTO peakpicking VALUES (?,?,?,?,?,?,?,?)""", ctuple)
        except sqlite3.IntegrityError as e:
            pass

        ###We just insert the sample name
        samples = c.execute("SELECT id,path,level FROM samples").fetchall()
        pid = 0
        counter_to_process = 0
        counter_processed = 0
        ###We create the processing table and add the input if they don t exists.
        for vid,path,level in samples:
            temp_l = list(args)+[path]
            hash_sample = hash_list(temp_l)
            # / output / ADAP / msms / msms_1_3146b653da32a8d0806ec9c635c629e5.mgf
            # / output / ADAP / peaktables / peaktable_1_3146b653da32a8d0806ec9c635c629e5.csv
            path_ms = os.path.join(path_peaktables,"peaktable_"+str(vid)+"_"+hash_sample+".csv")
            if level=="MS2":
                path_ms = "NOT PROCESSED"
            else :
                pid += 1
            path_msms_s = os.path.join(path_msms,"msms_"+str(vid)+"_"+hash_sample+".mgf")
            try:
                ctuple = (pid, 1, vid, path, hash_sample, path_ms, path_msms_s, 1, 1)
                c.execute("""INSERT INTO processing VALUES (?,?,?,?,?,?,?,?,?)""", ctuple)
                counter_to_process += 1
            except sqlite3.IntegrityError as e:
                counter_processed += 1
        logging.info(str(counter_to_process)+" peakpicking to do.")
        conn.commit()
        conn.close()
