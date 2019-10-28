import multiprocessing as mp
import subprocess

def run_cl(cl):
    print(cl)
    subprocess.call(cl,shell=True)



class ParallelRunner:
    def __init__(self,max_jobs,chunksize=20):
        ###The maximum number of jbos can never go above th enumber of corese
        if max_jobs>mp.cpu_count():
            max_jobs = mp.cpu_count()-1
        self.max_jobs=max_jobs
        self.chunksize=chunksize


    def run(self,command_lines,silent=True,log=None):
        try:
            from subprocess import DEVNULL  # py3k
        except ImportError:
            import os
            DEVNULL = open(os.devnull, 'wb')
        supp_str=""
        if log is not None:
            supp_str = " > "+log

        with mp.Pool(self.max_jobs) as pool:
            results = pool.starmap(run_cl, [(cl+supp_str,) for cl in command_lines])
        return results


# if __name__=="__main__":
#     pl = parallel_runner(100,100)
#
#     clines = ["echo "+str(i)+"t" for i in range(0,10)]
#
#     tres = pl.run(clines)
#
#     print(tres)
