import multiprocessing as mp
import concurrent.futures
import subprocess
import os

def run_cl(cl,timeout=None):
    ###We always include the environement
    my_env = os.environ.copy()
    # print(cl[0])
    if timeout is None:
        subprocess.call(cl[0],shell=True,env=my_env)
    else:
        subprocess.call(cl[0],shell=True,timeout=timeout,env=my_env)

class ParallelRunner:
    def __init__(self,max_jobs,chunksize=1):
        ###The maximum number of jbos can never go above th enumber of corese
        if max_jobs>mp.cpu_count():
            max_jobs = mp.cpu_count()-1
        self.max_jobs=max_jobs
        self.chunksize=1


    def run(self,command_lines,silent=True,log=None,timeout=None):
        try:
            from subprocess import DEVNULL  # py3k
        except ImportError:
            import os
            DEVNULL = open(os.devnull, 'wb')
        supp_str=""
        if log is not None:
            supp_str = " >> "+log+" 2>&1"
        ###We recompute the chunksize to give as many as possibl to each size
        # chunksize = min([len(command_lines)//self.max_jobs])
        # if chunksize != len(command_lines)/self.max_jobs:
        #     chunksize <- chunksize+1
        largs = [(cl+supp_str,timeout) for cl in command_lines]
        with concurrent.futures.ProcessPoolExecutor(min(self.max_jobs,len(command_lines))) as executor:
            executor.map(run_cl, largs)#,chunksize=chunksize)
