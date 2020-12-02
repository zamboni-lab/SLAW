import multiprocessing as mp
import concurrent.futures
import subprocess
import os
import logging

def run_cl_solo(cl,timeout=None,error=True,output=True):
    ###We always include the environement
    my_env = os.environ.copy()
    if timeout is None:
        process = subprocess.Popen(cl,shell=True, env=my_env, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
    else:
        process = subprocess.Popen(cl, shell=True, env=my_env, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    try:
        str_out = stdout.decode('UTF-8').strip()
        str_err = stderr.decode('UTF-8').strip()
        if len(str_out)>1 and output:
            logging.debug(str_out)
        if len(str_err) > 1 and error:
            logging.debug(str_err)
    except UnicodeDecodeError:
        pass

def run_cl(cl,timeout=None):
    ###We always include the environement
    my_env = os.environ.copy()
    if timeout is None:
        process = subprocess.Popen(cl[0], shell=True, env=my_env, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
    else:
        process = subprocess.Popen(cl[0], shell=True, env=my_env, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,timeout=timeout)
    stdout, stderr = process.communicate()
    str_out = stdout.decode('UTF-8').strip()
    str_err = stderr.decode('UTF-8').strip()
    if len(str_out)>1:
        logging.debug(str_out)
    if len(str_err) > 1:
        logging.debug(str_err)

class ParallelRunner:
    def __init__(self,max_jobs,chunksize=1,timeout=None):
        ###The maximum number of jbos can never go above th enumber of corese
        if max_jobs>mp.cpu_count():
            max_jobs = mp.cpu_count()-1
        self.max_jobs=int(max_jobs)
        self.chunksize=1

    def run(self,command_lines,silent=True,log=None,timeout=None):
        if len(command_lines)==1:
            run_cl_solo(command_lines[0],timeout=timeout)
        else:
            supp_str=""
            # if log is not None:
            #     supp_str = " >> "+log+" 2>&1"
            largs = [(cl+supp_str,timeout) for cl in command_lines]
            with mp.Pool(min(self.max_jobs,len(command_lines))) as executor:
                executor.map(run_cl, largs)
