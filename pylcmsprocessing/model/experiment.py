import sqlite3
import multiprocessing
import os
import logging
import time
from functools import partial
import common.references as cr
import common.tools as ct
import model.helper.output_handler as oh
import model.helper.parameters as params
import random
import model.helper.inputbuilder as ib
import model.helper.parallel_runner as pr
import model.steps.peakpicking as mp
import model.steps.grouping as mg
import model.steps.annotating_adducts_isotopes as mai
import model.steps.post_processing as pp
import model.steps.information_completion as ic
import model.steps.filtering_peaktable as fp
import common.slaw_exception as cs
import pandas as pd
import math
import shutil

def is_converted(csvfile):
    try:
        tpd = pd.read_csv(csvfile, nrows=0, sep=",", header=0)
    except pd.io.common.CParserError:
        tpd = pd.read_csv(csvfile, nrows=1, sep=",", header=0)
    if "SN" in tpd.columns:
        return True
    return False


####We ensure that ther is a single databse connection at the same timnb
class Experiment:

    #####Databases peaks
    def __init__(self, db,save_db=None, reset=False, temp_outdir = None):

        path_db = os.path.abspath(db)
        if os.path.isfile(path_db):
            if reset:
                os.remove(path_db)
        self.db = path_db
        self.path_save_db = None
        if save_db is not None:
            self.path_save_db = os.path.abspath(save_db)
        self.parameters = True
        self.peakpicking = []
        self.output = None
        if temp_outdir is not None :
            self.output = oh.OutputDirectory(temp_outdir)
        else:
            if self.db_exists():
                self.output = oh.OutputDirectory(self.get_outdir())
            elif self.path_save_db is not None and os.path.isfile(self.path_save_db):
                self.load_db()
                self.output = oh.OutputDirectory(self.get_outdir())


    def db_exists(self):
        return os.path.exists(self.db)

    ##Saving the DB outside the docker eventually
    ##Write DB outside the docker
    def save_db(self):
        if self.path_save_db is not None:
            shutil.copy(self.db,self.path_save_db)

    ##Copy the DB outside the docker.
    def load_db(self):
        if not os.path.isfile(self.db):
            if os.path.isfile(self.path_save_db):
                shutil.copy(self.path_save_db,self.db)

    def get_workers(self, open=True):
        if open:
            self.open_db()
        c = self.conn.cursor()
        c.execute('''SELECT max_jobs FROM common''')
        val = c.fetchall()[0][0]
        if open:
            self.close_db()
        return val

    def get_outdir(self, open=True):
        if open:
            self.open_db()
        c = self.conn.cursor()
        c.execute('''SELECT outdir FROM common''')
        val = c.fetchall()[0][0]
        if open:
            self.close_db()
        return val

    def guess_polarity(self,input_path=None):
        path_raw = None
        output = self.output.getFile(cr.OUT["POLARITY"])
        pscript = os.path.join(ct.find_rscript(), "get_polarity.R")

        if input_path is None:
            raw_files = self.get_query("SELECT path FROM samples WHERE level='MS1'")
            sel_samp = len(raw_files)//2
            path_raw = raw_files[sel_samp][0]
        else:
            all_path = [os.path.join(input_path,pp) for pp in os.listdir(input_path) if
                        pp.upper().endswith("MZML") or pp.upper().endswith("MZXML")]
            path_raw = all_path[0]

        logging.info("Guessing polarity from file:"+os.path.basename(path_raw))
        args = ["Rscript",pscript,'"'+path_raw+'"', '"'+output+'"']
        cli = " ".join(args)
        ##We call the script eventually.
        par_run = pr.run_cl_solo(cli)
        # subprocess.call(cli,shell=True)
        ##We now read the output file
        with open(output, "r") as f:
            polarity = f.readline().rstrip()
        ###We remoe the path.
        os.remove(output)
        self.polarity = polarity
        os.environ["POLARITY"]=polarity
        return polarity


    def get_polarity(self, open=True):
        if open:
            self.open_db()
        c = self.conn.cursor()
        c.execute('''SELECT polarity FROM common''')
        val = c.fetchall()[0][0]
        if open:
            self.close_db()
        return val

    def add_initial_infos(self, max_jobs=None, out_dir=None, polarity=None):
        os.makedirs(os.path.dirname(self.db), exist_ok=True)
        os.makedirs(out_dir, exist_ok=True)

        self.open_db()
        c = self.conn.cursor()
        if not polarity in cr.DATA["IONMODE"]:
            raise Exception("Unknown polarity, it should be one of " + ",".join(cr.DATA["IONMODE"]))

        ###We creat ethe common table
        c.execute('''CREATE TABLE IF NOT EXISTS common(max_jobs INTEGER,outdir TEXT,polarity TEXT)''')

        ###If the value exists
        c.execute('''SELECT * FROM common''')
        common = c.fetchall()
        if len(common)==0:
            if max_jobs is None or out_dir is None:
                raise Exception("Missing 'max_jobs' or 'out_dir' arguments please furnish them")
            values = (int(max_jobs), out_dir, polarity)
            c.execute('INSERT INTO common VALUES (?,?,?)', values)
        if len(common)!=0:
            out_dir = common[0][1]

        ###Creating the output directory
        out_dir = os.path.abspath(out_dir)
        self.output = oh.OutputDirectory(out_dir)
        self.close_db()

    ####Opening db plus foreign key constraints
    def open_db(self):
        self.conn = sqlite3.connect(self.db)
        self.conn.execute("PRAGMA foreign_keys = 1")

    def close_db(self, commit=True):
        if commit:
            self.conn.commit()
        self.conn.close()

    def get_query(self, query):
        self.open_db()
        c = self.conn.cursor()
        queries = self.conn.execute(query).fetchall()
        self.close_db()
        return queries

    def update_query_construction(self, table, id, idvals, field, value):
        if type(value) is str:
            value = '"' + value + '"'
        if type(idvals) is str:
            idvals = '"' + idvals + '"'
        query = "UPDATE " + table + " SET " + field + "=" + str(value) + " WHERE " + id + "=" + str(idvals)
        return query

    ##We consider that the database i open
    def table_is_empty(self, table):
        res = self.conn.execute("SELECT COUNT(*) FROM " + table + "  LIMIT 1").fetchone()
        return res[0] == 0

    def add_processing_method(self, algorithms, num_parameters=50):
        ####Processing method will all be furnishe dbut they are now handled at the beginning
        softwares = ""
        try:
            softwares = [cr.ALGORITHMS_TABLE[x] for x in algorithms]
        except KeyError:
            logging.warning("Unknown algorithm "+ ",".join(algorithms)+ " known algorithms are "+
                  ",".join(cr.ALGORITHMS_TABLE))
        ####We check that all the parameters exists
        dir_data = ct.find_data()
        for idx in range(len(softwares)):
            s = softwares[idx]
            s[1] = os.path.join(dir_data, s[1])
            s[2] = os.path.join(dir_data, s[2])
            if not os.path.isfile(s[1]):
                raise OSError("File", s[1], "does not exists.")
            softwares[idx] = s
        ###We evetually add them to the databse
        self.open_db()
        c = self.conn.cursor()
        ###We creat ethe common table
        c.execute('''CREATE TABLE IF NOT EXISTS algorithms
                             (id INTEGER PRIMARY KEY,
                             name TEXT,
                             software TEXT,
                              json TEXT,
                              xml TEXT,
                              summary TEXT,
                              num_parameters INTEGER,
                              seed INTEGER)''')

        ###We verifye if there is already sample in the databse
        if self.table_is_empty("algorithms"):
            ###Completion of the input tuples
            if type(num_parameters) is int:
                num_parameters = [num_parameters] * len(algorithms)

            values = [(idx + 1, algorithms[idx], alg[0], alg[1], alg[2], "", num_parameters[idx],
                       random.randint(0, 10000)) for idx, alg in enumerate(softwares)]

            ###We check all the values of algoirthms
            c.execute("SELECT * FROM algorithms")
            res = c.fetchall()

            for v in values:
                c.execute('INSERT INTO algorithms VALUES (?,?,?,?,?,?,?,?)', v)
        self.close_db()

    def get_samples(self):
        return self.get_query("SELECT path FROM samples WHERE level='MS1'")

    ###Return a list of the data matrix eventually
    def get_datamatrix(self):
        return self.get_query("SELECT peaktable FROM peakpicking")


    def build_samples(self, path_samples, path_ms2=None, is_optim=False):

        ###We build the majority
        pscript = os.path.join(ct.find_rscript(), "createSQLiteexperiment.R")
        cline = "Rscript " + pscript + " -d '" + path_samples + "' -b '" + self.db + "'"
        if path_ms2 is not None:
            cline = cline + " -o '"+path_ms2+"'"
        if "SLAWSUMMARY" in os.environ and not is_optim:
            cline = cline + " -s "+ os.environ["SLAWSUMMARY"]
        logging.info(cline)
        pr.run_cl_solo(cline,error=False,output=False)
        # subprocess.call(cline, shell=True)
        ##We check if the database has been created
        self.open_db()

        cursor = self.conn.cursor()
        # get the count of tables with the name
        cursor.execute(''' SELECT count(name) FROM sqlite_master WHERE type='table' AND name='samples' ''')
        # if the count is 1, then table exists
        if cursor.fetchone()[0] != 1:
            raise FileNotFoundError("Table ''samples' was not created correctly. Check your input files format and your summary.csv file.")
        cursor.close()
        self.close_db()


    def build_peakpicking(self):
        ###For each kind of algorithm we add the
        self.open_db()
        c = self.conn.cursor()

        ###Creation of the table
        c.execute('''CREATE TABLE IF NOT EXISTS peakpicking
        (id INTEGER PRIMARY KEY,
        algorithm INTEGER,
        parameter TEXT,
        hash TEXT,
        peaktable TEXT,
        fused_msms TEXT,
        annotated_peaktable_full TEXT,
        annotated_peaktable_reduced TEXT,
        CONSTRAINT known_software
        FOREIGN KEY (algorithm)
        REFERENCES algorithms(id),
        CONSTRAINT hashing UNIQUE (hash)
        )''')

        self.close_db()


    def build_processing(self):
        self.open_db()
        c = self.conn.cursor()
        ###We just create the table
        c.execute('''CREATE TABLE IF NOT EXISTS processing
          (id INTEGER PRIMARY KEY,
          peakpicking INTEGER,
          sample INTEGER,
          input TEXT NOT NULL,
          hash TEXT NOT NULL,
          output_ms TEXT NOT NULL,
          output_ms2 TEXT NOT NULL,
          valid INTEGER,
          step INTEGER NOT NULL,
          CONSTRAINT fk_peakpicking
          FOREIGN KEY (peakpicking)
          REFERENCES peakpicking(id),
          CONSTRAINT fk_sample
          FOREIGN KEY (sample)
          REFERENCES samples(id),
          CONSTRAINT hashing_input UNIQUE (hash)
          )''')
        self.close_db()

    def initialise_database(self, max_jobs, outdir, polarity, path_samples, algorithms, num_parameters=50, path_ms2=None, is_optim=False):
        ###We construct the initia software in the method
        self.add_initial_infos(max_jobs, outdir, polarity)
        self.build_samples(path_samples,path_ms2=path_ms2,is_optim=is_optim)
        ###If the output directory is not existing we create it
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        self.add_processing_method(algorithms, num_parameters)
        self.build_peakpicking()
        self.build_processing()

    #Parameters handling part

    #Sample the required number of parameters form the
    def sample_parameters(self, seed=512):
        #For each processing method we sample the parameter
        algs = self.get_query("SELECT * FROM algorithms")

        self.open_db()
        c = self.conn.cursor()

        for alg in algs:
            #We first check if there is not already a tabl

            dir_summary = cr.OUT[alg[1]]["SUMMARY"]
            path_summary = self.output.getFile(dir_summary)

            #We create sampler object
            pgenerator = params.ParametersGenerator(alg[3], max_parameters=alg[6], seed=seed)
            pgenerator.make_combination()
            pgenerator.sample_combinations()
            pgenerator.export_combinations(path_summary)

            #We update the dtaabase every time
            up_query = self.update_query_construction("algorithms", "id", alg[0], "summary", path_summary)
            c.execute(up_query)

        self.close_db()

    ###Creation of the input part for MZmine

    def building_inputs_single_processing(self, xml_file, algorithm="ADAP"):
        algs = self.get_query("SELECT * FROM algorithms")
        for alg in algs:
            if alg[2] == "MZmine":
                builder = ib.MZMineBuilder(alg, self.db, self.output, alg[0], algorithm=algorithm)
                builder.build_inputs_single_parameter_set(xml_file)

    def reset_processing(self):
        query = "DELETE FROM processing"
        self.open_db()
        c = self.conn.cursor()
        self.conn.execute(query).fetchall()
        self.close_db()


    def building_inputs_evaluations(self):
        algs = self.get_query("SELECT * FROM algorithms")
        for alg in algs:
            if alg[2] == "MZmine":
                builder = ib.MZMineBuilder(alg, self.db, self.output, alg[0])
                builder.build_inputs()
                builder.filter_inputs()
                builder.clean()

    ###Running the peak picking

    def find_software_from_peakpicking(self):
        self.open_db()
        c = self.conn.cursor()
        c.execute("SELECT algorithm FROM peakpicking")
        algos = c.fetchall()
        c.execute("SELECT software FROM algorithms")
        softwares = c.fetchall()

        ###We build a list giving the name of the software for each peakpicking method.
        softwares = [softwares[algos[i][0] - 1][0] for i in range(0, len(algos))]

        self.close_db()
        return softwares

    def run_mzmine(self, pmzmine, xml_file, algorithm = "ADAP", batch_size=2000, silent=True, log = None, input_only=False):
        logging.info("Building input")
        self.building_inputs_single_processing(xml_file,algorithm=algorithm)
        if input_only:
            return None
        num_workers = self.get_workers()
        runner = pr.ParallelRunner(num_workers)
        logging.info("Building softwares")
        softwares = self.find_software_from_peakpicking()

        ####we generate the command line for all the inputs
        self.open_db()
        c = self.conn.cursor()
        c.execute("SELECT * FROM processing WHERE output_ms!='NOT PROCESSED'")

        ###We modify it to
        num_incorrect_files = 0
        tot_files = 0
        while True:
            rows = c.fetchmany(batch_size)
            if not rows: break
            logging.info("Building peakpicking")
            peakpickings = [mp.PeakPickingMZmine(row, pmzmine) for row in rows]
            need_processing = [x.need_computing() for x in peakpickings]
            processing = 0
            while any(need_processing):
                ####Peak picking
                clis = [x.command_line_processing() for x in peakpickings if x.need_computing()]
                ####We run the jobs actually
                if len(clis) > 0:
                    runner.run(clis, silent=silent, log = log, timeout = cr.CONSTANT["PEAKPICKING_TIMOUT"])
                names_output = [x.get_output() + "\n" for x in peakpickings]
                name_temp = cr.TEMP["CONVERSION"]
                path_temp = self.output.getFile(name_temp)
                summary = open(path_temp, "w+")
                summary.writelines(names_output)
                #
                pjoin = os.path.join(ct.find_rscript(), "wrapper_MZmine_peak_table_conversion.R ")
                # ###Calling the script on all the processed path
                cline = "Rscript " + pjoin + " " + path_temp + " " + str(self.get_workers(open=False))
                pr.run_cl_solo(cline)
                # subprocess.call(cline, shell=True)
                need_processing = [not os.path.exists(row[5]) for row in rows]
                processing += 1
                if processing >= 2:
                    break
                tot_files += len(rows)
            if processing >=2:
                num_incorrect_files += sum(need_processing)
                ###Now we jus thave to update the peak table
                for irow in range(len(need_processing)):
                    row = rows[irow]
                    if need_processing[irow]:
                        pid = row[0]
                        update_query = self.update_query_construction("processing","id",str(pid),"valid","0")
                        c.execute(update_query)
        logging.info("Peak picking finished "+str(num_incorrect_files)+" files not processed on a total of "+str(tot_files))
        self.close_db()

    def run_openms(self,min_fwhm,max_fwhm,fwhm_fac,snt,ppm,min_int,max_outlier,min_points,quant,silent=True, log = None,batch_size=10000):
        ###We collect the output path
        builder = ib.openMSBuilder(self.db,output=self.output)
        builder.build_inputs_single_parameter_set(min_fwhm,max_fwhm,fwhm_fac,snt,ppm,min_int,max_outlier,min_points,quant)
        num_workers = self.get_workers()
        runner = pr.ParallelRunner(num_workers)
        ####we generate the command line for all the inputs
        self.open_db()
        c = self.conn.cursor()
        c.execute("SELECT * FROM processing WHERE output_ms!='NOT PROCESSED'")

        while True:
            rows = c.fetchmany(batch_size)
            if not rows: break
            if quant=="intensity":
                quant = "area"
            peakpickings = [mp.PeakPickingOpenMS(row, min_fwhm,max_fwhm,fwhm_fac,snt,ppm,min_int,max_outlier,min_points,quant) for row in rows]
            ids = [row[0] for row in rows]
            need_processing = [x.need_computing() for x in peakpickings]
            ids_to_convert = [idv for idv,pp in zip(ids,need_processing) if pp]
            ####Peak picking
            temp = [(x.get_output(),x.command_line_processing()) for x in peakpickings if x.need_computing()]
            clis = [x[1] for x in temp]
            names_output = [x[0] for x in temp]

            ####We run the jobs actually
            if len(clis) > 0:
                runner.run(clis, silent=silent, log = log, timeout = cr.CONSTANT["PEAKPICKING_TIMOUT"])
            #We convert the peaktables back to the data.
            pjoin = os.path.join(ct.find_rscript(), "FromFeaturesMLToDf.R")
            # ###CWe create the command line to allow the peakpicking
            if len(names_output) > 0:
                names_converted = [x.split(".")[0] + ".csv" for x in names_output]
                clis_conversion = ["Rscript " + pjoin + ' "'+old+ '" "' +new+ '"' for old,new in zip(names_output,names_converted)]
                runner.run(clis_conversion, silent=silent, log = log, timeout = cr.CONSTANT["PEAKPICKING_TIMOUT"])
                to_retry = []
                for pid,nn,o in zip(ids_to_convert,names_converted,names_output):
                    if os.path.isfile(nn):
                        update_query = self.update_query_construction("processing", "id", str(pid), "output_ms", nn)
                        c.execute(update_query)
                    else:
                        to_retry.append((pid,nn,o))
                ###If they are not converted  we wait 5s and do it.
                if len(to_retry)>0:
                    ###5 second waiting time to unlock file if needed.
                    time.sleep(2)
                    clis_conversion_retry = ["Rscript " + pjoin +' "' + old +'" "'+ new+'"' for id,new,old in
                                       to_retry]
                    runner.run(clis_conversion_retry, silent=silent, log=log, timeout=cr.CONSTANT["PEAKPICKING_TIMOUT"])
                    for pid,nn,o in to_retry:
                        if os.path.isfile(nn):
                            update_query = self.update_query_construction("processing", "id", str(pid), "output_ms", nn)
                            c.execute(update_query)
                        else:
                            update_query = self.update_query_construction("processing", "id", str(pid), "valid", "0")
                            c.execute(update_query)
        self.close_db()

    def run_xcms(self,min_peakwidth,max_peakwidth,snt,ppm,min_int,min_points,silent=True, log = None,batch_size=10000):
        ###We collect the output path
        builder = ib.xcmsBuilder(self.db,output=self.output)
        builder.build_inputs_single_parameter_set(min_peakwidth,max_peakwidth,snt,ppm,min_int,min_points)
        num_workers = self.get_workers()
        runner = pr.ParallelRunner(num_workers)
        ####we generate the command line for all the inputs
        self.open_db()
        c = self.conn.cursor()
        c.execute("SELECT * FROM processing WHERE output_ms!='NOT PROCESSED'")
        while True:
            rows = c.fetchmany(batch_size)
            if not rows: break
            #row, min_fwhm, max_fwhm, snt, ppm, min_int, min_points
            peakpickings = [mp.PeakPickingXCMS(row, min_peakwidth,max_peakwidth,snt,ppm,min_int,min_points) for row in rows]
            ids = [row[0] for row in rows]
            need_processing = [x.need_computing() for x in peakpickings]
            ids_to_process = [idv for idv,pp in zip(ids,need_processing) if pp]
            ####Peak picking
            temp = [(x.get_output(),x.command_line_processing()) for x in peakpickings if x.need_computing()]
            clis = [x[1] for x in temp]
            names_output = [x[0] for x in temp]
            ####We run the jobs actually
            if len(clis) > 0:
                runner.run(clis, silent=silent, log=log, timeout = cr.CONSTANT["PEAKPICKING_TIMOUT"])
            if len(names_output) > 0:
                for nn,vid in zip(names_output,ids_to_process):
                    if not os.path.isfile(nn):
                        update_query = self.update_query_construction("processing", "id", str(vid), "valid", "0")
                        c.execute(update_query)
        self.close_db()

    ###Try to correct he MZmine ocrrection in case processing was interrupted
    def correct_conversion(self):
        path_temp = self.output.getFile(cr.TEMP["CONVERSION"])
        self.open_db()
        c = self.conn.cursor()
        c.execute("SELECT output_ms FROM processing WHERE output_ms!='NOT PROCESSED' AND valid=1")
        processings = c.fetchall()
        p_conversion = os.path.join(ct.find_rscript(), "wrapper_MZmine_peak_table_conversion.R ")
        to_convert = [x[0] + "\n" for x in processings if not is_converted(x[0])]
        if len(to_convert) > 0:
            # logging.warning("correcting " + str(len(to_convert)) + " files:","|".join(to_convert))
            with open(path_temp, "w+") as summary:
                summary.writelines(to_convert)
            cline = "Rscript " + p_conversion + " " + path_temp + " " + str(self.get_workers(open=False))
            pr.run_cl_solo(cline)
            # subprocess.call(cline, shell=True)
        else:
            logging.info("No corrections to do")
        self.close_db()
        self.save_db()

    ###This step is always done before grouping
    def post_processing_peakpicking_mzmine(self,algorithm="ADAP"):
        ###We rename the peaktable to get more meaningful names.
        self.open_db()
        c = self.conn.cursor()
        # c.execute("SELECT path,id,output_ms,output_ms2 FROM samples INNER JOIN processing on samples.id=processing.sample WHERE output_ms!='NOT PROCESSED' AND valid=1")
        # all_infos = c.fetchall()
        c.execute("SELECT id,output_ms,output_ms2 FROM processing WHERE output_ms!='NOT PROCESSED' AND valid=1")
        all_peaktable = c.fetchall()

        c.execute("SELECT path FROM samples WHERE level='MS1'")
        all_samples = c.fetchall()
        dirname = self.output.getDir(cr.OUT[algorithm]["PEAKTABLES"])
        dirname_msms = self.output.getDir(cr.OUT[algorithm]["MSMS"])
        ###Now we jsut rename all the files
        for ip in range(0,len(all_peaktable)):
            ##Correcting if needed convert he ms files
            new_name = os.path.join(dirname,os.path.splitext(os.path.basename(all_samples[ip][0]))[0]+".csv")
            if os.path.isfile(new_name):
                continue
            old_name = all_peaktable[ip][1]
            os.rename(old_name,new_name)
            query = self.update_query_construction("processing", "id", str(all_peaktable[ip][0]), "output_ms", new_name)
            c.execute(query)
            full_name = all_peaktable[ip][2]
            try:
                if os.path.getsize(full_name)==0:
                    os.remove(full_name)
                else:
                    new_name = os.path.join(dirname_msms,os.path.splitext(os.path.basename(all_samples[ip][0]))[0]+".mgf")
                    os.rename(full_name,new_name)
                    query = self.update_query_construction("processing", "id", str(all_peaktable[ip][0]), "output_ms2", new_name)
                    c.execute(query)
            except FileNotFoundError: #File has already been converted.
                continue
            name_quant = os.path.splitext(full_name)[0]+"_quant.csv"
            if os.path.isfile(name_quant):
                os.remove(name_quant)
        self.close_db()
        self.save_db()

    ###This step is always done before grouping
    def post_processing_peakpicking_openms(self):
        ###We rename the peaktable to get more meaningful names.
        self.open_db()
        c = self.conn.cursor()
        c.execute("SELECT processing.id,path,output_ms,output_ms2 FROM samples INNER JOIN processing on samples.id=processing.sample WHERE level='MS1' AND output_ms!='NOT PROCESSED' AND valid=1")
        res = c.fetchall()
        dirname = self.output.getDir(cr.OUT["OPENMS"]["PEAKTABLES"])
        dirname_msms = self.output.getDir(cr.OUT["OPENMS"]["MSMS"])
        ###Now we jsut rename all the files
        for id,sample,peaktable,ms2 in res:
            ##Correcting if needed convert he ms files
            out_dir = os.path.dirname(peaktable)
            raw_samp = os.path.basename(sample).split(".")[0]
            new_name = os.path.join(out_dir,raw_samp+".csv")
            if os.path.isfile(new_name):
                continue
            os.rename(peaktable,new_name)
            query = self.update_query_construction("processing", "id", id, "output_ms", new_name)
            c.execute(query)

        self.close_db()
        self.save_db()

    ###This step is always done before grouping
    def post_processing_peakpicking_xcms(self):
        ###We rename the peaktable to get more meaningful names.
        self.open_db()
        c = self.conn.cursor()
        c.execute("SELECT processing.id,path,output_ms,output_ms2 FROM samples INNER JOIN processing on samples.id=processing.sample WHERE level='MS1' AND output_ms!='NOT PROCESSED' AND valid=1")
        res = c.fetchall()
        ###Now we jsut rename all the files
        for id,sample,peaktable,ms2 in res:
            ##Correcting if needed convert he ms files
            out_dir = os.path.dirname(peaktable)
            raw_samp = os.path.basename(sample).split(".")[0]
            new_name = os.path.join(out_dir,raw_samp+".csv")
            if os.path.isfile(new_name):
                continue
            os.rename(peaktable,new_name)
            query = self.update_query_construction("processing", "id", id, "output_ms", new_name)
            c.execute(query)
        self.close_db()
        self.save_db()

    def filter_peaktables(self,filtering_string):
        """Filter all the peaktables according to the filtering string"""
        peak_filter = fp.PeaktableFilter(filtering_string)
        self.open_db()
        c = self.conn.cursor()
        c.execute("SELECT output_ms FROM samples INNER JOIN processing on samples.id=processing.sample WHERE level='MS1' AND output_ms!='NOT PROCESSED' AND valid=1")
        all_peaktables = c.fetchall()
        self.close_db()
        self.save_db()
        all_peaktables = [x[0] for x in all_peaktables]
        num_workers = self.get_workers()
        logging.info("Starting peaktable filtration")
        with multiprocessing.Pool(min(num_workers, len(all_peaktables))) as executor:
            executor.map(partial(fp.par_peaktable_filtering,peak_filter=peak_filter),all_peaktables)
        logging.info("Done peaktables filtration")


    def extract_ms2(self,output,noise_level=0,all=False):
        num_workers = self.get_workers()
        dir_out= self.output.getDir(output)
        pscript = os.path.join(ct.find_rscript(),"extractingMGFfromMZML.R")
        cli = " ".join(["Rscript",pscript,'"'+self.db+'"','"'+dir_out+'"',str(noise_level),str(num_workers),str(all)])
        pr.run_cl_solo(cli)
        # subprocess.call(cli,shell=True)
        logging.info("MS2 extraction finished")

    ###At the moment there is a single grouping method
    def group(self, max_workers=2, silent=False, intensity="height",mztol=0.007,rttol=0.02, log=None):
        num_workers = self.get_workers()
        runner = pr.ParallelRunner(min(num_workers, max_workers))
        ####We create all the grouper eventually
        self.open_db()
        c = self.conn.cursor()
        c.execute("SELECT * FROM peakpicking")
        all_peakpicking = c.fetchall()
        groupers = [0] * len(all_peakpicking)
        countgroup = 0
        dir_temp = self.output.getDir(cr.TEMP["GROUPING"])
        dir_datamatrix = self.output.getDir(cr.OUT["DATAMATRIX"])
        for pp in all_peakpicking:
            ###name of file
            flists = c.execute("SELECT output_ms FROM processing WHERE peakpicking = " + str(pp[0]))
            flists = [f[0] for f in flists]

            ###getting samples
            nlists = c.execute("SELECT path FROM samples WHERE level='MS1'")
            nlists = [os.path.basename(nn[0]) for nn in nlists]

            ppg = mg.Grouper(pp, dir_temp, dir_datamatrix, flists, nlists, intensity,mztol,rttol)
            ###We update the peaktable direction of the file
            poutput_dm = ppg.get_output_datamatrix()
            poutput_idx = ppg.get_output_index()
            ppg.make_temp_file()
            query = self.update_query_construction("peakpicking", "id", str(pp[0]), "peaktable", poutput_dm)
            c.execute(query)
            query = self.update_query_construction("peakpicking", "id", str(pp[0]), "index_file", poutput_idx)
            c.execute(query)
            groupers[countgroup] = ppg
            countgroup += 1
        if countgroup != 0:
            groupers = groupers[0:countgroup]
            clis = [g.command_line() for g in groupers if g.need_computing()]
            if len(clis) > 0:
                runner.run(clis, silent=silent, log=log)
        self.close_db()
        self.save_db()
        logging.info("Grouping finished")


    def group_online(self,intensity="int",
    ppm = 15,mztol=0.007,rttol=0.02,n_ref = 150,
                     ms2_mz_tol=0.01,ms2_rt_tol=0.05,
                     alpha = 0.1,filter_qc=0.5,fold_blank=3,log=None,post_processing=True):
        num_workers = self.get_workers()
        runner = pr.ParallelRunner(num_workers)
        ####We create all the grouper eventually
        self.open_db()
        c = self.conn.cursor()
        c.execute("SELECT * FROM peakpicking")
        all_peakpicking = c.fetchall()
        self.close_db()
        groupers = [0] * len(all_peakpicking)
        countgroup = 0

        #Alignement output
        dir_blocks = self.output.getDir(cr.TEMP["GROUPING"]["BLOCKS"])
        dir_alignment = self.output.getFile(cr.TEMP["GROUPING"]["ALIGNMENT"])
        dir_datamatrix = self.output.getDir(cr.OUT["DATAMATRIX"])
        path_fig = "/output/figure"

        ###Spectra fusing output
        path_temp_1 = self.output.getFile(cr.TEMP["FUSING"]["TEMP1"])
        path_temp_2 = self.output.getFile(cr.TEMP["FUSING"]["TEMP2"])
        path_fused_msms = self.output.getFile(cr.OUT["FUSED_MSMS"])

        for pp in all_peakpicking:
            ###name of file
            self.open_db()
            c = self.conn.cursor()
            ###getting samples
            ppg = mg.OnlineGrouper(pp,self.db,dir_blocks, dir_alignment,
            dir_datamatrix, intensity, mztol, ppm, rttol, n_ref, alpha,
                                   ms2_mz_tol,ms2_rt_tol,path_fused_msms,
                                   path_temp_1,path_temp_2,num_workers,path_fig,
                                   filter_qc,fold_blank)
            ###We update the peaktable path
            poutput_dm = ppg.get_output_datamatrix()
            poutput_mgf = ppg.get_fused_mgf()
            query_align = self.update_query_construction("peakpicking", "id", str(pp[0]), "peaktable", poutput_dm)
            c.execute(query_align)
            query_fusing = self.update_query_construction("peakpicking", "id", str(pp[0]), "fused_msms", poutput_mgf)
            c.execute(query_fusing)
            self.close_db()
            groupers[countgroup] = ppg
            countgroup += 1
        if countgroup != 0:
            groupers = groupers[0:countgroup]
            clis_align = [g.command_line_aligning() for g in groupers if g.need_computing()]
            clis_filtering = [g.command_line_filtering() for g in groupers if g.need_computing()]
            clis_fusing = [g.command_line_fusing_msms() for g in groupers if g.need_computing()]
            if len(clis_align) > 0:
                logging.info("Aligning")
                runner.run(clis_align, log=log)
                if post_processing:
                    logging.info("Filtering")
                    runner.run(clis_filtering, log=log)
                    logging.info("Extracting consensus MS-MS spectra")
                    runner.run(clis_fusing, log=log)
        self.save_db()
        logging.info("Alignment finished")
    #Parallelism is handled by R always a single trhead in this case
    def annotate_ions(self, nfiles, ppm, dmz, adducts=None, main_adducts=None, max_workers=2, min_filter = 2):
        num_workers = self.get_workers()
        polarity = self.get_polarity()
        #We create all the grouper eventually
        self.open_db()
        c = self.conn.cursor()
        c.execute("SELECT * FROM peakpicking")
        all_peakpicking = c.fetchall()
        #We create the ouput of the adducts files.
        path_temp_adducts = self.output.getFile(cr.TEMP["IONANNOTATION"]["FULL"])
        path_temp_adducts_main = self.output.getFile(cr.TEMP["IONANNOTATION"]["MAIN"])
        runner = pr.ParallelRunner(min(num_workers, max_workers))
        successfull_processing = True
        #We get the polarity
        if adducts is None:
            if polarity == "positive":
                adducts = cr.default_adducts_positive()
                main_adducts = cr.default_adducts_main_positive()
            elif polarity == "negative":
                adducts = cr.default_adducts_negative()
                main_adducts = cr.default_adducts_main_negative()

        annotaters = [0] * len(all_peakpicking)
        count_annot = 0

        #If the adducts are not specified we load them
        for pp in all_peakpicking:
            path_datamatrix = self.output.getDir(cr.OUT["DATAMATRIX"])
            ### If the data matrix does not exist we skip to the next iteration
            if not os.path.isfile(pp[4]):
                successfull_processing=False
                continue

            ppg = mai.IonAnnotater(pp[3], self.db, pp[4], path_datamatrix, polarity, cr.DATA["IONANNOTATION"]["XCMS_MODEL"], num_workers, nfiles,
                                   ppm, dmz, min_filter, adducts, main_adducts)

            #We update the output of the processing into the database.
            ppg.write_adducts(self.output)
            path_data_matrices = ppg.get_output_datamatrices()
            query = self.update_query_construction("peakpicking", "id", str(pp[0]), "annotated_peaktable_full",
                                                   path_data_matrices[1])
            c.execute(query)
            query = self.update_query_construction("peakpicking", "id", str(pp[0]), "annotated_peaktable_reduced",
                                                   path_data_matrices[0])
            c.execute(query)
            annotaters[count_annot] = ppg
            count_annot += 1
        if count_annot != 0:
            annotaters = annotaters[0:count_annot]
            clis = [ann.command_line(self.output) for ann in annotaters if ann.need_computing()]
            if len(clis) > 0:
                runner.run(clis, silent=True)
        self.close_db()
        self.save_db()
        logging.info("Annotation finished")
        return successfull_processing

    def add_missing_informations(self,max_iso, max_charge, quant, ppm, dmz):
        num_workers = self.get_workers()
        runner = pr.ParallelRunner(num_workers)
        ###Data of previous steps
        path_isotopes = cr.DATA["ISOTOPES"]
        path_rt_model = self.output.getFile(cr.TEMP["GROUPING"]["ALIGNMENT"])
        path_temp_1 = self.output.getFile(cr.TEMP["FUSING"]["TEMP1"])
        self.open_db()
        c = self.conn.cursor()
        c.execute("SELECT * FROM peakpicking")
        all_peakpicking = c.fetchall()
        margin_mz = 0.005

        expanders = [0] * len(all_peakpicking)
        count_expand = 0

        #The number of files is determined based on the available memory
        num_files = 3
        if "TOTAL_SLAW_MEMORY" in os.environ:
            num_files = max(math.floor(int(os.environ["TOTAL_SLAW_MEMORY"])/2000),3)
        if num_files > 10:
            num_files = 10
        
        for pp in all_peakpicking:
            ###pp 4,10,11
            ie = ic.InformationExpander(self.db, pp[4], path_temp_1, path_rt_model, path_isotopes,
                                         max_iso, max_charge, quant,
                                         margin_mz, ppm, dmz,num_files, num_workers)
            expanders[count_expand] = ie
            count_expand += 1
        if count_expand != 0:
            expanders = expanders[0:count_expand]
            clis = [iexp.command_line() for iexp in expanders if iexp.need_computing()]
            if len(clis) > 0:
                runner.run(clis, silent=True)
        logging.info("Gap filling and isotopic pattern extraction finished.")

    def post_processing(self,targets,path_raw_files=None,mztol=0.05,rttol=0.03):
        if path_raw_files is None:
            self.open_db()
            c = self.conn.cursor()
            raw_files=c.execute("SELECT path FROM samples WHERE level='MS1'")
            raw_files=[rr[0] for rr in raw_files]
            self.close_db()
            self.save_db()
            path_raw_files = self.output.getFile(cr.TEMP["POSTPROCESSING"])

        runner = pr.ParallelRunner(1)
        ###We check if a targetted table exist
        num_workers = self.get_workers()
        self.open_db()
        c = self.conn.cursor()
        c.execute("SELECT * FROM peakpicking")
        all_peakpicking = c.fetchall()
        self.close_db()

        pprocessors = [0] * len(all_peakpicking)
        count_pprocessors = 0

        for app in all_peakpicking:
            path_dm = app[4]
            if not os.path.isfile(path_dm):
                continue
            path_peaks = self.output.getFile(cr.OUT["FIGURES"]["PEAKS"]+app[3]+".pdf")
            path_diagnosis = self.output.getFile(cr.OUT["FIGURES"]["DIAGNOSIS"]+app[3]+".pdf")
            path_tab_rt = os.path.join(self.output.getDir(cr.OUT["DATAMATRIX"]),cr.OUT["TARGET"]["RT"]+app[3]+".csv")
            path_tab_int = os.path.join(self.output.getDir(cr.OUT["DATAMATRIX"]),cr.OUT["TARGET"]["INT"]+app[3]+".csv")
            path_hdf5 = self.output.getFile(cr.OUT["EIC"]+app[3]+".pdf")
            ppv = pp.PostProcessing(self.db,path_targets=targets,path_fig_target=path_peaks,
            path_fig_summary=path_diagnosis,path_tab_rt=path_tab_rt,path_tab_int=path_tab_int,
            path_output_hdf5=path_hdf5,num_workers=num_workers,raw_files=path_raw_files,mztol=mztol,rttol=rttol)
            pprocessors[count_pprocessors] = ppv
            count_pprocessors += 1
        if count_pprocessors != 0:
            pprocessors = pprocessors[0:count_pprocessors]
            clis = [ppv.command_line() for ppv in pprocessors]
            if len(clis) > 0:
                runner.run(clis, silent=True)
        logging.info("Diagnosis figures printed")

    ###change the path on an experiment to
    def rebase_experiment(self,path_db):
        vconn = sqlite3.connect(path_db)
        c = vconn.cursor()

        self.open_db()
        cnew = self.conn.cursor()
        res = vconn.execute("SELECT id,output_ms,valid FROM processing").fetchall()
        vconn.close()
        for idx, path, valid in res:
            query = ''' UPDATE processing SET output_ms = ? , valid = ? WHERE id = ?'''
            cnew.execute(query, (path, valid, idx))
        self.close_db()

    def clean(self):
        ###We clena the files only if the ions have been correctly annotated.
        self.open_db()
        c = self.conn.cursor()
        c.execute("SELECT annotated_peaktable_reduced FROM peakpicking")
        path_annotation = c.fetchall()[0][0]
        self.close_db()
        if os.path.isfile(path_annotation):
            to_rm= [cr.TEMP["IONANNOTATION"]["FULL"],cr.TEMP["IONANNOTATION"]["MAIN"],
            cr.TEMP["CONVERSION"],cr.OUT["ADAP"]["JSON"],cr.OUT["ADAP"]["CANDIDATES"]]
            # to_rm= [cr.TEMP["GROUPING"]["TEMP"],cr.TEMP["IONANNOTATION"]["FULL"],cr.TEMP["IONANNOTATION"]["MAIN"],
            # cr.TEMP["CONVERSION"],cr.OUT["ADAP"]["JSON"],cr.OUT["ADAP"]["CANDIDATES"],cr.TEMP["DIR"]]
            for waste in to_rm:
                pwaste = self.output.getPath(waste)
                if os.path.isdir(pwaste):
                    shutil.rmtree(pwaste, ignore_errors=True)
                elif os.path.exists(pwaste):
                    os.remove(pwaste)
