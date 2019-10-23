import sqlite3
import subprocess
import os
import common.references as cr
import common.tools as ct
import model.output_handler as oh
import model.parameters as params
import random
import model.inputbuilder as ib
import model.parallel_runner as pr
import model.peakpicking as mp
import model.grouping as mg
import model.evaluating as me
import model.comparing_evaluation as mce
import model.annotating_adducts_isotopes as mai
import pandas as pd


def is_converted(csvfile):
    tpd = pd.read_csv(csvfile, nrows=2, sep=",", header=1)
    if len(tpd.columns) >= 6:
        return False
    return True


####We ensure that ther is a single databse connection at the same timnb
class Experiment:

    ###THe path of the data is useless
    # The maximum number of workoers ot cal l durieng paralle processing
    ###parameters JSON for open MS and MZmine ex "X:/Documents/dev/script/tools/lcmsprocessingtools/tests/openMS/json_parameters.json"
    ###OUTDIR an output directory it needs ot be empty, it is created
    ###Number of tested parsmeters set.

    #####Databases peaks

    def __init__(self, db, reset=False):
        path_db = os.path.abspath(db)
        if os.path.isfile(path_db):
            if reset:
                os.remove(path_db)
        self.db = path_db
        self.parameters = True
        self.peakpicking = []
        self.output = None
        if self.db_exists():
            self.output = oh.OutputDirectory(self.get_outdir())

    def db_exists(self):
        return os.path.exists(self.db)

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
        if len(common) is 0:
            if max_jobs is None or out_dir is None:
                raise Exception("Missing 'max_jobs' or 'out_dir' arguments please furnish them")
            values = (int(max_jobs), out_dir, polarity)
            c.execute('INSERT INTO common VALUES (?,?,?)', values)
        if len(common) is not 0:
            out_dir = common[0][1]

            print("common: ", common)

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
            print("Unknown algorithm ", ",".join(algorithms), " known algorithms are ",
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
            for r in res:
                print(r)

            for v in values:
                c.execute('INSERT INTO algorithms VALUES (?,?,?,?,?,?,?,?)', v)
        self.close_db()

    def get_samples(self):
        return self.get_query("SELECT path FROM sample")

    def build_samples(self, path_samples):

        ###We build the majority
        pscript = os.path.join(ct.find_rscript(), "createSQLiteexperiment.R")
        cline = "Rscript " + pscript + " -d " + path_samples + " -b " + self.db
        subprocess.call(cline, shell=True)

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
        hash_peaktable TEXT,
        index_file TEXT,
        evaluation TEXT,
        hash_evaluation TEXT,
        annotated_peaktable_full TEXT,
        annotated_peaktable_reduced TEXT,
        CONSTRAINT known_software
        FOREIGN KEY (algorithm)
        REFERENCES algorithms(id),
        CONSTRAINT hashing UNIQUE (hash)
        )''')

        self.close_db()

    ###We jsut initalize the peaktable
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
          output TEXT NOT NULL,
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

    def initialise_database(self, max_jobs, outdir, polarity, path_samples, algorithms, num_parameters=50):
        ###We construct the initia software in the method
        self.add_initial_infos(max_jobs, outdir, polarity)
        self.build_samples(path_samples)
        ###If the output directory is not existing we create it
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        self.add_processing_method(algorithms, num_parameters)
        self.build_peakpicking()
        self.build_processing()
        print("Database intialized")

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

    def building_inputs_single_processing(self, xml_file):
        algs = self.get_query("SELECT * FROM algorithms")
        for alg in algs:
            if alg[2] == "MZmine":
                builder = ib.MZMineBuilder(alg, self.db, self.output, alg[0])
                builder.build_inputs_single_parameter_set(xml_file)

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

    def run(self, pmzmine, batch_size=40, silent=True):
        num_workers = self.get_workers()
        runner = pr.ParallelRunner(num_workers)
        softwares = self.find_software_from_peakpicking()

        ####we generate the command line for all the inputs
        self.open_db()
        c = self.conn.cursor()
        c.execute("SELECT * FROM processing")

        ###We modify it to

        while True:
            rows = c.fetchmany(batch_size)
            if not rows: break

            peakpickings = [mp.PeakPickingMZmine(row, pmzmine) for row in rows]
            need_processing = [x.need_computing() for x in peakpickings]

            while any(need_processing):
                ####Peak picking
                clis = [x.command_line_processing(hide=False) for x in peakpickings if x.need_computing()]
                ####We run the jobs actually
                if len(clis) > 0:
                    runner.run(clis, silent=silent)

                names_output = [x.get_output() + "\n" for x in peakpickings]

                name_temp = cr.TEMP["CONVERSION"]
                path_temp = self.output.getFile(name_temp)

                summary = open(path_temp, "w+")

                summary.writelines(names_output)
                #
                pjoin = os.path.join(ct.find_rscript(), "wrapper_MZmine_peak_table_conversion.R ")
                #
                # ###Calling the script on all the processed path
                cline = "Rscript " + pjoin + " " + path_temp + " " + str(self.get_workers(open=False))
                subprocess.call(cline, shell=True)
                need_processing = [not os.path.exists(row[5]) for row in rows]
            else:
                print("no processing")
        print("Peak picking finished")
        self.close_db()

    ###Try to correct he MZmine ocrrection in case processing was interrupted
    def correct_conversion(self):
        path_temp = self.output.getFile(cr.TEMP["CONVERSION"])
        self.open_db()
        c = self.conn.cursor()
        c.execute("SELECT output FROM processing")
        processings = c.fetchall()
        p_conversion = os.path.join(ct.find_rscript(), "wrapper_MZmine_peak_table_conversion.R ")
        to_convert = [x[0] + "\n" for x in processings if not is_converted(x[0])]
        if len(to_convert) > 0:
            print("correcting " + str(len(to_convert)) + " files.")
            with open(path_temp, "w+") as summary:
                summary.writelines(to_convert)
            cline = "Rscript " + p_conversion + " " + path_temp + " " + str(self.get_workers(open=False))
            subprocess.call(cline, shell=True)
        else:
            print("nothing to correct")
        self.close_db()

        ###At the moment there is a single grouping method

    def group(self, max_workers=2, silent=False, intensity="height",mztol=0.007,rttol=0.02):
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
            flists = c.execute("SELECT output FROM processing WHERE peakpicking = " + str(pp[0]))
            flists = [f[0] for f in flists]

            ###getting samples
            nlists = c.execute("SELECT path FROM samples")
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
            clis = [g.command_line() for g in groupers]
            if len(clis) > 0:
                runner.run(clis, silent=silent)
        self.close_db()
        print("Grouping finished")

    def evaluate(self, max_workers=2, silent=False):
        num_workers = self.get_workers()
        runner = pr.ParallelRunner(min(num_workers, max_workers))

        ####We create all the grouper eventually
        self.open_db()
        c = self.conn.cursor()
        c.execute("SELECT id,hash,peaktable FROM peakpicking")
        all_eval = c.fetchall()
        to_evaluate = [0] * len(all_eval)
        count_eval = 0
        file_replicates = self.output.getFile(cr.TEMP["REPLICATES"])
        ###Wrting replicates
        c.execute("SELECT replicate FROM samples")
        all_replicates = [str(a[0]) + "\n" for a in c.fetchall()]
        with open(file_replicates, "w+") as summary:
            summary.writelines(all_replicates)
        dir_evaluation = self.output.getDir(cr.OUT["EVALUATION"])
        for eval in all_eval:

            output = os.path.join(dir_evaluation, "eval_" + eval[1] + ".csv")
            oev = me.Evaluate(eval[2], output, file_replicates)

            if oev.need_computing():
                query = self.update_query_construction("peakpicking", "id", str(eval[0]), "evaluation",
                                                       oev.get_output())
                c.execute(query)
                to_evaluate[count_eval] = oev
                count_eval += 1

        if count_eval != 0:
            evals = to_evaluate[0:count_eval]
            clis = [g.command_line() for g in evals]
            if len(clis) > 0:
                runner.run(clis, silent=silent)
        self.close_db()

        ####We check all the smaple to evaluate in a loop first

    ####This part is not done in paralllel
    def compare_evaluation(self):

        cfigure = self.output.getFile(cr.OUT["RES_EVALUATION"]["FIGURE"])
        cparam = self.output.getFile(cr.OUT["RES_EVALUATION"]["PARAM"])
        cpeaktables = self.output.getDir(cr.OUT["RES_EVALUATION"]["PEAKTABLES"])

        ###We perform the comparison

        ce = mce.evaluationComparator(self.db, cparam, cfigure, cpeaktables)

        cli = ce.command_line()

        ###The comparison is always done on single ocore
        subprocess.call(cli, shell=True)

    #Annotations
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
            clis = [ann.command_line(self.output) for ann in annotaters]
            if len(clis) > 0:
                runner.run(clis, silent=True)
        self.close_db()
        print("Annotation finished")
