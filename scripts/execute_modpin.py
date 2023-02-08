import os, sys, re
import ConfigParser
import argparse
import shutil
import subprocess
import time


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


      
     
def parse_user_arguments(*args, **kwds):
    parser = argparse.ArgumentParser("execute_modpin -i PPI_LIST -o OUTPUT_DIR -seq SEQUENCES")
    parser.add_argument('-i', '--query_list', dest = 'query_list', action = 'store',
                        help = 'Input file with a list of pairs of proteins to group and test')
    parser.add_argument('-l', '--label_name', dest = 'label', action = 'store', default = None,
                        help = 'Label to store the analyses (default is the name of the input query)')
    parser.add_argument('-seq', '--sequences', dest = 'sequence_file', action = 'store',
                        help = 'Input file with a list of sequences in FASTA format')
    parser.add_argument('-o', '--output_directory', dest = 'output_dir', action = 'store', default = 'ModPPI',
                        help = 'Output directory (default is ModPPI)')
    parser.add_argument('-n', '--number_of_models', dest = 'nmodels', action = 'store', default = 1, type=int,
                        help = 'Number of models for each template (default is 1)')
    parser.add_argument('-d', '--dummy_dir', dest = 'dummy_dir', action = 'store', default = '/tmp/modppi_dummy',
                        help = 'Specifies the location of the dummy folder (default is /tmp/modppi_dummy)')
    parser.add_argument('-opt', '--optimize', dest = 'optimize', action = 'store_true',
                        help = 'Flag to allow model optimization (default is False)')
    parser.add_argument('-3did', '--use_domain_interactions', dest = 'did', action = 'store_true',
                        help = 'Flag to include domain-domain interactions from 3DiD (default is False)')
    parser.add_argument('-v', '--verbose', dest = 'verbose', action = 'store_true',
                        help = 'Flag for verbose mode (default is False)')
    parser.add_argument("-j","--parallel", default=False, action="store_true", dest="parallel", 
                        help="Submit JOBS to Queues in parallel (default = False)")
    parser.add_argument('-hydro','--hydrogens', dest = 'hbplus', action = 'store_true',
                        help = 'Flag to include hydrogens in modelling (default is False, and always True for analyses)')
    parser.add_argument('-r','--renumerate', dest = 'renumerate' , action = 'store_true',
                        help = 'Flag to renumber the sequences as in the original FastA (default is False, and always True for analyses)')
    parser.add_argument("-info", "--information_file", default=None, action="store",  dest="info",
                        help="Information LOG file of MODELS that have failed and have been completed")
    parser.add_argument("-reuse","--reusing_models", default=False, action="store_true", dest="reuse", 
                        help="Reuse the information files. If the flag is used then models that had failed will remain as FAILED, otherwise it tries to redo them (default=False)")
    parser.add_argument("--complete",default=0.95, action="store", dest="complete", 
            help="Ratio of completness over the total number of interactions to be done(default= 0.95). This is useful in the server to stop MODPIN when the time of execution exceeds more than 48 hours ")
 
 
    options = parser.parse_args()
    return options


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options   = parse_user_arguments()
    parallel  = options.parallel
    info_file = options.info
    verbose   = options.verbose
    reuse     = options.reuse
    complete  = float(options.complete)
    dummy_dir = options.dummy_dir
    if not dummy_dir.startswith("/"): dummy_dir = os.path.abspath(options.dummy_dir)
    if not os.path.exists(dummy_dir): os.makedirs(dummy_dir)
    modppi    = config.get("Paths","modppi")
    
    # Iterate untill all protein sequences are modelled
    ppi_list = set()
    query    = open(options.query_list,"r")
    for line in query:
        if line.startswith("#"):continue
        ppi = line.strip().split()
        ppi_list.add(ppi[0]+"::"+ppi[1])
    query.close()
    submitted=set()
    n_done=0
    if options.verbose: print("Start iteration to model")
    info_file=options.info
    if info_file is None: info_file=os.path.join(options.output_dir,"modelling_ppi_list_execution.log")
    if reuse and functions.fileExist(info_file):
         if options.verbose: print ("Reuse previous information of runs from file %s"%info_file)
    else:
         if options.verbose: print ("Open to write %s"%info_file)
         log_file = open(info_file,"w")
         log_file.write("#List of MODELS\n")
         if os.path.exists(os.path.join(options.output_dir,"interactions_done.list")):
            interactions_done = open(os.path.join(options.output_dir,"interactions_done.list"),"r")
            for line in interactions_done:
                if line.startswith("#"):continue
                words=line.split()
                if len(words)<2: continue
                log_file.write("%s::%s\t%s\n"%(words[0],words[1],words[2]))
            interactions_done.close()
         log_file.close()


    done    = functions.done_jobs(info_file)
    iterate = functions.check_done(done,ppi_list)
    n_done  = 0
    maxtime =  3600 * 4
    start_time = d_time = 0
    while(iterate):
        for ppi in ppi_list:
            if ppi in submitted:continue
            submitted.add(ppi)
            if verbose: print("Modelling %s ..."%ppi)
            if reuse:
              log_file = open(info_file,"r")
              skip_adding=False
              for line in log_file:
                       if ppi in line.split(): skip_adding=True
              if skip_adding and verbose: print("\t-- Already in the information file ")
              if skip_adding: submitted.add(ppi)
              if skip_adding: continue
              log_file.close()
            #paremeter values
            ppi_file   =  os.path.join(dummy_dir,ppi+".ppi")
            if not os.path.exists(ppi_file):
               ppi_input  =  open(ppi_file,"w")
               ppi_input.write("%s\t%s\n"%(ppi.split("::")[0], ppi.split("::")[1]))
               ppi_input.close()
            
            parameters= " -ppi " + ppi_file
            parameters= parameters + " -skip "
            parameters= parameters + " -seq " + options.sequence_file
            parameters= parameters + " -o " + options.output_dir
            parameters= parameters + " -n " + str(options.nmodels)
            parameters= parameters + " -d " + dummy_dir
            if options.label is not None: parameters = parameters + " -l " + options.label
            if options.did :              parameters = parameters + " -3did "
            if options.verbose:           parameters = parameters + " --verbose "
            if options.optimize:          parameters = parameters + " --optimize "
            if options.renumerate:        parameters = parameters + " --renumerate "
            if options.hbplus:            parameters = parameters + " --hydrogens "
            if parallel:
              if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
              else: cluster_queue=config.get("Cluster", "cluster_queue")
              program=modppi
              python = os.path.join(config.get("Paths", "python_path"), "python")
              if verbose: print("\t-- Submit  %s %s" % (program,parameters))
              functions.submit_command_to_queue("%s %s %s" % (python,program,parameters), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
              submitted.add(ppi)      
            else:
              os.system("%s %s %s" % (python,program,parameters))
              submitted.add(ppi)      
            
        #Check next iteration, profiles submitted and profiles done
        log_file = open(info_file,"w")
        log_file.write("#List of MODELS\n")
        if os.path.exists(os.path.join(options.output_dir,"interactions_done.list")):
            interactions_done = open(os.path.join(options.output_dir,"interactions_done.list"),"r")
            for line in interactions_done:
                if line.startswith("#"):continue
                words=line.split()
                if len(words)<2: continue
                log_file.write("%s::%s\t%s\n"%(words[0],words[1],words[2]))
            interactions_done.close()
        log_file.close()
        done    = functions.done_jobs(info_file)
        iterate = functions.check_done(done,ppi_list) 
        if len(done) > n_done:
           n_done=len(done)
           if n_done> 1 and start_time==0: start_time = float(time.time())
           if options.verbose: 
                sys.stdout.write("Number of PPIs already done %d\n"%n_done)
                if d_time>0: sys.stdout.write("Time: %f\n"%d_time)
                sys.stdout.write("\t-- Check files done ...\n")
                log_file = open(info_file,"r")
                for line in log_file:
                   print("\t\t-- %s"%line.strip())
                log_file.flush()
                log_file.close()
                sys.stdout.write("\t-- Still running PPI models  %s ...\n"%functions.check_done(done,ppi_list))
                sys.stdout.write("\t-- Continue iteration %s ...\n"%iterate)
                sys.stdout.flush()
        #Check next iteration, if exceeding time and enough profiles stop iteration
        if start_time > 0:
           current_time = float(time.time())
           d_time = current_time - start_time
           if float(  len(done) ) / len(ppi_list) > complete  and d_time > maxtime: iterate=False
           if d_time > maxtime and iterate and complete < 1.0: 
               if options.verbose: sys.stdout.write("Time: %f Done: %d (%d) Ratio to end: %f\n"%(d_time,len(done),len(ppi_list),complete))
               complete = complete - 0.001
               #Limit the time to 1.5 days maximum
               if complete>0 and maxtime<125000: maxtime  = maxtime  + 60
               sys.stdout.flush()



    #Done
    if verbose:sys.stdout.write("Done\n")





