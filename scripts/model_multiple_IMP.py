import os, sys, re
#SHIVA
#import configparser as ConfigParser
#HYDRA
import ConfigParser
import optparse
import shutil
import subprocess
import time
import numpy as np

# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.dirname(__file__))

# Append scripts path to python path #
sys.path.append(scripts_path)

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Get python path #
python = os.path.join(config.get("Paths", "python_path"), "python")


# Import my functions #
import functions


#-------------#
# Functions   #
#-------------#


def get_time():

    from datetime import datetime
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y  %H:%M:%S")
    return dt_string




def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python model_multiple_IMP.py -i input_folder  -o output [--info log_file --dummy dummy_dir  --reuse  --parallel --verbose ]")

    parser.add_option("-i", action="store", type="string", dest="input_folder", default=None, help="Input folder, contains the inputs orftopology and restraints to model macro-complexes", metavar="{directory}")
    parser.add_option("-o", action="store", type="string", dest="output", default=None, help="Output folder rootname  ", metavar="{directory}")
    parser.add_option("--info",default=None,action="store", type="string", dest="info",help="Information LOG file of MODELS that have failed and have been completed (default modelling_IMP_execution.log)")
    parser.add_option("--reuse",default=False, action="store_true", dest="reuse", help="Reuse the information files. If the flag is used then models that had failed will remain as FAILED, otherwise it tries to redo them (default=False)")
    parser.add_option("--parallel", default=False, action="store_true", dest="parallel", help="Run in paralleli models of the directory (default = False)")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)", metavar="{boolean}")
    parser.add_option("--max_models",default=None,action="store", dest="max_models",help="Maximum number of models to submit or execute (default is None to model all). WARNING if no limit is used the number of models and the computational requirements may be large")
  

    (options, args) = parser.parse_args()

    if options.input_folder is None or options.output is None :
        parser.error("missing arguments: type option \"-h\" for help")

    return options


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options   = parse_options()
    parallel  = options.parallel
    info_file = options.info
    verbose   = options.verbose
    reuse     = options.reuse
    max_models= options.max_models
    input_folder = options.input_folder
    output    = options.output

    #Collect data from configuration
    if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
    else: cluster_queue=config.get("Cluster", "cluster_queue")
    max_jobs_in_queue = int(config.get("Cluster","max_jobs_in_queue"))
    cluster_submit    = config.get("Cluster","cluster_submit")
    cluster_qstat     = config.get("Cluster","cluster_qstat")
    command_queue     = os.path.join(scripts_path,config.get("Cluster","command_queue3"))


    if not input_folder.startswith("/"):  input_folder = os.path.abspath(input_folder)
    if not output.startswith("/"):   output = os.path.abspath(output)
    if info_file is None: info_file=os.path.join(input_folder,"modelling_IMP_execution.log")



    # Create list of models to do
    premodels=[]
    for input_file in os.listdir(input_folder):
        if input_file.endswith(".topology.txt"): 
            premodels.append(input_file.rstrip(".topology.txt"))

    # Confirm set of models
    set_of_models=set()
    for m in premodels:
        if os.path.exists(os.path.join(input_folder,m+".topology.fasta")) and os.path.exists(os.path.join(input_folder,m+".restraints.csv")): 
            set_of_models.add(m)


    # Iterate untill all model sequences are modelled
    submitted=set()
    if options.verbose: print("Start iteration to model with IMP")
    info_file=options.info
    if info_file is None: info_file=os.path.join(input_folder,"modelling_IMP_execution.log")
    if reuse and functions.fileExist(info_file):
         if options.verbose: print ("Reuse previous information of runs from file %s"%info_file)
    else:
         if options.verbose: print ("Open to write %s"%info_file)
         log_file = open(info_file,"w")
         log_file.write("#List of MODELS\n")
         log_file.close()


    done    = functions.done_jobs(info_file)
    iterate = functions.check_done(done,set_of_models)
    n_done  = 0

    if max_models is None: list_of_models=[x for x in set_of_models]
    else:                  list_of_models=[x for x in set_of_models][:int(max_models)]
    start_time = float(time.time())
    while(iterate):
        for model in list_of_models:
            if model in submitted:continue
            submitted.add(model)
            if verbose: print("Modelling %s ..."%model)
            if reuse:
              log_file = open(info_file,"r")
              skip_adding=False
              for line in log_file:
                       if model in line.split(): skip_adding=True
              if skip_adding and verbose: print("\t-- Already in the information file ")
              if skip_adding: submitted.add(model)
              if skip_adding: continue
              log_file.close()
            #paremeter values
            parameters =   " -i %s "%input_folder
            parameters =   parameters + " -m %s "%model
            parameters =   parameters + " -o %s "%output
            parameters =   parameters + " --info %s "%info_file
            logfile    =   output+"_"+model+".log"
            program=os.path.join(scripts_path,"model_IMP.py")
            if parallel:
              if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
              else: cluster_queue=config.get("Cluster", "cluster_queue")
              #python = os.path.join(config.get("Paths", "python_path"), "python")
              if verbose: print("\t-- Submit  %s %s" % (program,parameters))
              functions.submit_command_to_queue("%s %s %s >& %s " % ("python",program,parameters,logfile), cluster_queue, max_jobs_in_queue,command_queue,options.dummy_dir,cluster_submit,cluster_qstat)
              submitted.add(model)      
            else:
              if verbose: print("\t-- Execute  %s %s" % (program,parameters))
              os.system("%s %s %s >& %s " % ("python",program,parameters,logfile))
              submitted.add(model)      
        #Check next iteration, profiles submitted and profiles done
        done    = functions.done_jobs(info_file)
        iterate = functions.check_done(done,set_of_models) 
        current_time = float(time.time())
        # if more than 1.5 days running stop iterating
        if current_time - start_time > 125000: 
            iterate=False
            sys.stdout.write("\t-- Exceeding time (limited to 1.5 days)\n")
        
        if len(done) > n_done or np.mod(current_time,7200)==0:
           n_done=len(done)
           if options.verbose: 
                sys.stdout.write("TIME %s\n"%get_time())
                sys.stdout.write("Number of IMP models already done %d\n"%n_done)
                sys.stdout.write("\t-- Check files done ...\n")
                log_file = open(info_file,"r")
                for line in log_file:
                   print("\t\t-- %s"%line.strip())
                log_file.flush()
                log_file.close()
                sys.stdout.write("\t-- Still running model IMP  %s ...\n"%functions.check_done(done,set_of_models))
                sys.stdout.write("\t-- Continue iteration %s ...\n"%iterate)
                sys.stdout.flush()


    #Done
    if verbose:sys.stdout.write("Done\n")



