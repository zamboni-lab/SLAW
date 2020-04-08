import os
import csv

 ###Create command line
#path_input	cores	optimization	noptim	repet	xcms	comment


PATH_OUTPUT = "/cluster/scratch/dalexis/res_evaluation_paper"

PATH_CLI = os.path.join(PATH_OUTPUT,"clis.txt")
PATH_EXPERIMENTS = os.path.join(PATH_OUTPUT,"experiments.csv")
PATH_OUTPUT = os.path.join(PATH_OUTPUT,"coded_experiments")
PATH_LOGS = os.path.join(PATH_OUTPUT,"logs")
PATH_EXPERIMENTS = "C:/Users/dalexis/Documents/dev/lcmsprocessing_docker/external/experiments.csv"
###We create
with open(PATH_EXPERIMENTS, newline='') as csvfile:
   spamreader = csv.reader(csvfile, delimiter=';')
   next(spamreader,None)
   for row in spamreader:
       input,num_cores,optim,noptim,memory,xcms,ignore2 = row
       name_dir = os.path.basename(os.path.dirname(input))
       id = name_dir+"n"+str(num_cores)+"_o"+optim+"_no"+noptim+"_x"+xcms
       log_err = os.path.join(PATH_LOGS,id+".err")
       log_out = os.path.join(PATH_LOGS,id+".out")
        cline = 'bsub -n '+str(num_cores)+' -N -R singularity -R "rusage[mem=',str(int(memory)+500),']" -o '+
        log_out+' -e '+log_err+' -W 24:00 "SINGULARITYENV_OPTIM='+optim+' SINGULARITYENV_MEMORY=3500'

'''bsub -n 51 -N -R singularity -R "rusage[mem=''',str(memory),
''']" -o ~/logs/random_rsm_ipo.log -e ~/logs/random_rsm_ipo.err -W 24:00 "SINGULARITYENV_OPTIM=random_rsm_ipo SINGULARITYENV_MEMORY=3500 SINGULARITYENV_NCORES=50 SINGULARITYENV_NOPTIM=100 singularity run -C -W /cluster/scratch/dalexis/data_philipp_vincenth -B /cluster/scratch/dalexis/data_philipp_vincenth/res_random:/output  -B /cluster/scratch/dalexis/data_philipp_vincenth/mzML:/input lcms_workflow_zamboni.simg"
'''





