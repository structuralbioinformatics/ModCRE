import os, sys, re
import ConfigParser
import numpy
import optparse
import pandas
import shutil
import socket
import subprocess
import pandas as pd
import pickle
import cPickle

# Get module path #
scripts_path = os.path.abspath(os.path.dirname(__file__))
sys.path.append(scripts_path)

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Import from module #
import pbm, contacts, functions, nr, spotentials, dssp, x3dna, tomtom, triads
import pwm_pbm as PWM

# Imports jbonet's module #
from SBI.data import aminoacids1to3, aminoacids_polarity_boolean, nitrogenous_bases
from SBI.structure import PDB

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python2.7 %prog --pbm=<pbm_dir> --pdb=<pdb_dir> --pwm=<pwm_dir> [--dummy=<dummy_dir> -o <output_dir> --start=<start_step> --stop=<stop_step> -v]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="<dummy_dir>")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="<output_dir>")
    parser.add_option("--pbm", action="store", type="string", dest="pbm_dir", help="PBM directory (output directory from pbm.py)", metavar="<pbm_dir>")
    parser.add_option("--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (output directory from pdb.py)", metavar="<pdb_dir>")
    parser.add_option("--pwm", action="store", type="string", dest="pwm_dir", help="PWM directory (from CIS-BP)", metavar="<pwm_dir>")
    parser.add_option("--start", default=1, action="store", type="int", dest="start_step", help="Start at a given step (default = 1; first)", metavar="<start_step>")
    parser.add_option("--stop", default=7, action="store", type="int", dest="stop_step", help="Stop at a given step (default = 9; last)", metavar="<stop_step>")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")
    parser.add_option("-p", "--paralelize", default=False, action="store_true", dest="paralelize", help="Run paralel jobs in cluster")
    parser.add_option("-a", "--all_templates", default=False, action="store_true", dest="all_templates", help="Use as many models as templates per TF, otherwise only the best model is built per TF (default is False)", metavar="boolean")


    (options, args) = parser.parse_args()

    if options.pbm_dir is None or options.pdb_dir is None or options.pwm_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options


if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()


    # Create output "main" subdirs #
    for subdir in ["clusters", "models", "dssp", "x3dna", "nr", "potentials", "msa", "pwms", "tomtom", "results", "contacts", "triads", "cis-bp_pwms"]:
        if not os.path.exists(os.path.join(os.path.abspath(options.output_dir), subdir)):
            os.makedirs(os.path.join(os.path.abspath(options.output_dir), subdir))

    ###############################
    # 1. Create grid search files #
    ###############################
    if options.verbose: sys.stdout.write("\nCreate files for grid search...\n")
    # Skip if starts later #
    if options.start_step <= 1:
        # Verbose mode #
        if options.verbose: sys.stdout.write("\n")
        ##############################
        # 1.1 Extract TF-PWM motifs  #
        ##############################
        if options.verbose: sys.stdout.write("\t\t-- extracting TF motifs...\n")
        # Initialize #
        families = {}
        # For each TF... #
        for file_name in os.listdir(os.path.join(os.path.abspath(options.pbm_dir), "tfs")):
            # Get TF object #
            tf_obj = pbm.TF(file_name=os.path.join(os.path.abspath(options.pbm_dir), "tfs", file_name))
            # For each sequence... #
            if options.verbose: sys.stdout.write("\t\t\t--read %s\n"%file_name) 
            # Get formatted family name #
            family = None
            family_list = tf_obj.get_families()
            for family in family_list:
                if options.verbose: sys.stdout.write("\t\t\t--family %s\n"%str(family)) 
            if len(family_list) > len(tf_obj._sequences):
                print "Error: wrong size between families",len(family_list),"sequences",len(tf_obj._sequences)
                continue
            for i in range(len(tf_obj._sequences)):
                # Skip if k-mers file does not exist #
                kmers_file = os.path.join(os.path.abspath(options.pbm_dir), "kmers", tf_obj.get_id() + "." + str(i) + ".txt")
                if not os.path.exists(kmers_file): continue
                try:
                    family = family_list[i]
                except:
                    family = family_list[0]
                families.setdefault(family, set()).add(("%s.%s" % (tf_obj.get_id(), str(i)), tf_obj._motifs[i]))        
        ##############################
        # 1.3 Build TF-DNA models    #
        ############################## 
        if options.verbose: sys.stdout.write("\t\t-- build TF-DNA models...\n")
        # Initialize #
        # For each family... #
        for family in families.keys():
            # For each TF sequence... #
            for tf_id, motif in sorted(families[family]):
                # Initialize #
                sequence_file = os.path.join(os.path.abspath(options.pbm_dir), "sequences", "%s.fa" % tf_id)
                label = tf_id.split(".")[0] + "_mL_" + motif + "_mL_" + family.replace(" ","_").replace("/","=") + "_mL_"
                # Verbose mode #
                if options.verbose: sys.stdout.write("\t\t\t-- modelling %s\n" % label)
                if options.all_templates:
                   flag="--n-model=1"
                else:
                   flag="--best"
                if options.paralelize:
                    # Submit to queue #
                    functions.submit_command_to_queue("python " + scripts_path + "/model_protein.py -i " + sequence_file + " -o " + os.path.join(options.output_dir, "models") + " -l " + label + " --p=" + os.path.abspath(options.pdb_dir) + " --dummy=" + options.dummy_dir + " --renumerate -v -s " + flag, None, int(config.get("Cluster", "max_jobs_in_queue")), os.path.join(scripts_path, config.get("Cluster", "command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                # Else... #
                else:
                    # Submit to queue #
                    os.system("python " + scripts_path + "/model_protein.py -i " + sequence_file + " -o " + os.path.join(options.output_dir, "models") + " -l " + label + " --p=" + os.path.abspath(options.pdb_dir) + " --dummy=" + options.dummy_dir + " --renumerate -v  -s " + flag)
                    
        # Exit if stops here #
        if options.paralelize:
           sys.stdout.write("Wait until all models are done and continue with start=2\n")
           sys.stdout.write("Exiting...\n\n")
           exit(0)
    if options.stop_step == 1 :
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ###########################################
    # 2.  Get DSSP, 3DNA, contacts and triads #
    ###########################################
    if options.start_step <= 2:
        if options.verbose: sys.stdout.write("\t\t-- get DSSP, 3DNA, contacts and triads information...\n")
        # For PDB file... #
        for pdb_file in os.listdir(os.path.join(os.path.abspath(options.output_dir), "models")):
            # If PDB file... #
            if pdb_file.endswith(".pdb"):
                # Skip if DSSP file already exists #
                pdb_obj = PDB(os.path.join(os.path.abspath(options.output_dir), "models", pdb_file))
                label = pdb_file.split(":")[0]
                dssp_file = os.path.join(os.path.abspath(options.output_dir), "dssp", label + ".txt")
                if not os.path.exists(dssp_file):
                    dssp_obj = dssp.get_dssp_obj(os.path.join(os.path.abspath(options.output_dir), "models", pdb_file))
                    dssp_obj.write(dssp_file)
                else:
                    dssp_obj = dssp.DSSP(file_name=dssp_file)
                # Skip if X3DNA file already exists #
                x3dna_file = os.path.join(os.path.abspath(options.output_dir), "x3dna", label + ".txt")
                if not os.path.exists(x3dna_file):
                    x3dna_obj = x3dna.get_x3dna_obj(os.path.join(os.path.abspath(options.output_dir), "models", pdb_file))
                    x3dna_obj.write(x3dna_file)
                else:
                    x3dna_obj = x3dna.X3DNA(file_name=x3dna_file)
                # Skip if contacts file already exists #
                contacts_file = os.path.join(os.path.abspath(options.output_dir), "contacts", label + ".txt")
                if not os.path.exists(contacts_file):
                    contacts_obj = contacts.get_contacts_obj(pdb_obj, x3dna_obj)
                    contacts_obj.write(contacts_file)
                else:
                    contacts_obj = contacts.Contacts(file_name=contacts_file)
                if len(contacts_obj.get_contacts())<1:
                    print("Missing Protein-DNA contacts in file: " + str(os.path.join(os.path.abspath(options.output_dir), "models", pdb_file)) + " ...\n")
                    continue
                # Skip if triads file already exists #
                triads_file = os.path.join(os.path.abspath(options.output_dir), "triads", label + ".txt")
                if not os.path.exists(triads_file):
                    triads_obj = triads.get_triads_obj(pdb_obj, dssp_obj, x3dna_obj, contacts_obj)
                    triads_obj.write(triads_file)

    # Exit if stops here #
    if options.stop_step == 2:
        sys.stdout.write("Exiting...\n\n")
        exit(0)


    #########################
    # 3.  Generate PWMs     #
    #########################
    if options.start_step <= 3:
        if options.verbose: sys.stdout.write("\t\t-- get all PWMs...\n")
        # For each PDB file #
        for pdb_file in os.listdir(os.path.join(os.path.abspath(options.output_dir), "models")):
            if pdb_file.endswith(".pdb"):
                # Get all complementary files # 
                label = pdb_file.split(":")[0]
                pdb_path = os.path.join(os.path.abspath(options.output_dir), "models", pdb_file)
                dssp_file = os.path.join(os.path.abspath(options.output_dir), "dssp", label + ".txt")
                x3dna_file = os.path.join(os.path.abspath(options.output_dir), "x3dna", label + ".txt")
                contacts_file = os.path.join(os.path.abspath(options.output_dir), "contacts", label + ".txt")
                triads_file = os.path.join(os.path.abspath(options.output_dir), "triads", label + ".txt")
                # Send job to queue to create pwms #
                if options.verbose: sys.stdout.write("\t\t\t-- run 'grid_search_get_pwms' to get PWM of %s ...\n"%label)
                if options.paralelize:
                    #print("%s %s -p %s -d %s -x %s -c %s -t %s --pdb=%s --pbm=%s -o %s" % ("python", os.path.join(scripts_path, "grid_search_get_pwms.py"), pdb_path, dssp_file, x3dna_file, contacts_file, triads_file, os.path.abspath(options.pdb_dir), os.path.abspath(options.pbm_dir), os.path.abspath(options.output_dir)))
                    functions.submit_command_to_queue("%s %s -v -p %s -d %s -x %s -c %s -t %s --pdb=%s --pbm=%s -o %s" % ("python", os.path.join(scripts_path, "grid_search_get_pwms.py"), pdb_path, dssp_file, x3dna_file, contacts_file, triads_file, os.path.abspath(options.pdb_dir), os.path.abspath(options.pbm_dir), os.path.abspath(options.output_dir)), None, int(config.get("Cluster", "max_jobs_in_queue")), os.path.join(scripts_path, config.get("Cluster", "command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                else:
                    os.system("%s %s -v -p %s -d %s -x %s -c %s -t %s --pdb=%s --pbm=%s -o %s" % ("python", os.path.join(scripts_path, "grid_search_get_pwms.py"), pdb_path, dssp_file, x3dna_file, contacts_file, triads_file, os.path.abspath(options.pdb_dir), os.path.abspath(options.pbm_dir), options.output_dir))

        # Exit if stops here #
        if options.paralelize:
           sys.stdout.write("Wait until all models are done and continue with start=4\n")
           sys.stdout.write("Exiting...\n\n")
           exit(0)
    if options.stop_step == 3 :
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    
    ###################################
    # 4.  Get meme formated database  #
    ###################################
    if options.start_step <= 4:
        if options.verbose:sys.stdout.write("\t--  Open %s\n"%os.path.abspath(options.pwm_dir))
        for pwm_file in os.listdir(os.path.abspath(options.pwm_dir)):
            if options.verbose:sys.stdout.write("\t\t-- Get MEME of %s\n"%pwm_file) 
            input_pwm = os.path.join(os.path.abspath(options.pwm_dir), pwm_file)
            output_pwm = os.path.join(os.path.abspath(options.output_dir), "cis-bp_pwms", pwm_file[:-4] + ".meme.s")
            try:
             pwm_obj=PWM.nMSA(input_pwm,pwm_file[:-4],"txt")
             pwm_obj.write(output_pwm,"meme",overwrite=True)
            except:
             if options.verbose:sys.stdout.write("\t\t\t-- Skip due to wrong format: %s\n"%pwm_file) 


    # Exit if stops here #
    if options.stop_step == 4:
        sys.stdout.write("Exiting...\n\n")
        exit(0)


    ###################################
    # 5.    Execute tomtom            #
    ###################################
    if options.start_step <= 5:
        # Build the cis-bp_pwms database
        jobs=[]
        cisbp_pwm_dir = os.path.join(os.path.abspath(options.output_dir), "cis-bp_pwms")
        database=os.path.join(os.path.abspath(options.output_dir), "cis-bp_pwms","tf_motifs.txt")
        if  os.path.exists(database): shutil.move(database,database+".old")
        for motif in os.listdir(cisbp_pwm_dir):
                if motif.endswith(".meme.s"):
                   os.system("cat %s >> %s" % (os.path.join(os.path.abspath(options.output_dir), "cis-bp_pwms",motif), database))
        for modcre_pwm_file in os.listdir(os.path.join(os.path.abspath(options.output_dir), "pwms")):
            if options.verbose:sys.stdout.write("\t\t-- parsing file %s\n"%modcre_pwm_file)
            if modcre_pwm_file.endswith(".meme.s"):
                modcre_pwm=modcre_pwm_file.rstrip(".meme.s")
                tf_id     = modcre_pwm.split("_mL_")[0] 
                motif     = modcre_pwm.split("_mL_")[1]
                family_tf = modcre_pwm.split("_mL_")[2]
                template  = modcre_pwm.split("_mL_")[3]
                potential = modcre_pwm.split("_mL_")[4]
                data      = modcre_pwm.split("_mL_")[5]
                dist      = modcre_pwm.split("_mL_")[6]
                threshold = modcre_pwm.split("_mL_")[7]
                pdb_chain = modcre_pwm.split("_mL_")[8]
                if options.verbose:sys.stdout.write("\t\t\t-- TF %s %s \n"%(tf_id,motif))
                if options.verbose:sys.stdout.write("\t\t\t-- data %s \n"%(data))
                if options.verbose:sys.stdout.write("\t\t\t-- dist %s \n"%(dist))
                if options.verbose:sys.stdout.write("\t\t\t-- threshold %s \n"%(threshold))
                if options.verbose:sys.stdout.write("\t\t\t-- chain %s \n"%(pdb_chain))
                if int(threshold) not in range(70, 101, 2):continue
                if "general" in potential:                
                    fam = "general" 
                    family = "general"
                else:
                    fam = "family"
                    family = template
                if options.verbose:sys.stdout.write("\t\t\t-- Family %s %s \n"%(fam,family))
                if "taylor" in potential:
                    taylor = "taylor."
                else:
                    taylor = ""
                if ".acc"  in potential:
                    comp ="acc"
                elif ".bins" in potential:
                    comp ="bins"
                else:
                    comp = ""
                if ".pmf."  in potential:
                    pmf ="pmf."
                else:
                    pmf = ""
                computation= pmf+taylor+comp
                if options.verbose:sys.stdout.write("\t\t\t-- Computation %s \n"%(computation))
                cisbp_pwm  = database
                # the correct thing is to search against all cis-bp and select later the pval of the motif
                tomtom_file = os.path.join(os.path.abspath(options.output_dir), "tomtom", tf_id + "_mL_" + motif + "_mL_" + family_tf + "_mL_" + template + "_mL_" + family + "." + computation + "_mL_" + data + "_mL_" + str(dist) + "_mL_" + str(threshold) + "_mL_" + pdb_chain + ".txt")
                if options.verbose:sys.stdout.write("\t\t\t-- tomtom of %s \n"%(tomtom_file))
                # Execute tomtom #
                if options.paralelize:
                  python=os.path.join(config.get("Paths","python_path"),"python")
                  if not os.path.isfile(tomtom_file):
                    if modcre_pwm_file != None and tomtom_file != None:
                       pwm_file=os.path.join(os.path.abspath(options.output_dir), "pwms", modcre_pwm_file)
                       if options.verbose:sys.stdout.write("\t\t\t-- Add job %s \n"%(pwm_file))
                       jobs.append('%s %s  -d %s -i %s -o %s --dummy=%s'%(python,os.path.join(scripts_path,"tomtom.py"),cisbp_pwm,pwm_file,tomtom_file,os.path.abspath(options.dummy_dir)))
                else:
                  if not os.path.isfile(tomtom_file):
                    if modcre_pwm_file != None and tomtom_file != None:
                        try:
                            tomtom_obj = tomtom.get_tomtom_obj(database_file=cisbp_pwm, pwm_file=os.path.join(os.path.abspath(options.output_dir), "pwms", modcre_pwm_file), dummy_dir=os.path.abspath(options.dummy_dir))
                            tomtom_obj.write(file_name=tomtom_file)
                        except:
                            print("ERROR executing tomtom for files: " + modcre_pwm_file + "  " + cisbp_pwm + "\n")
                    else:
                        print("Running tomtom repeatedly on the same file: " + tomtom_file) 

        if options.paralelize and len(jobs)>0:
                  # Create tmp directory #
                  tmp_dir = os.path.join(options.dummy_dir, "sh")
                  if os.path.exists(tmp_dir) == False:
                     os.mkdir(tmp_dir)
                  else:
                     shutil.rmtree(tmp_dir)
                     os.mkdir(tmp_dir)

                  # Initialize #
                  min_jobs_in_queue=int(config.get("Cluster","min_jobs_in_queue"))
                  if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
                  else: cluster_queue=config.get("Cluster", "cluster_queue")
                  command_queue=config.get("Cluster","command_queue")
                  submit=config.get("Cluster","cluster_submit")
                  n_run=0
                  n_jobs=0
                  for job in jobs:
                     if n_jobs == min_jobs_in_queue:
                        n_jobs=0
                        n_run=n_run+1
                        out.close()
                     if n_jobs==0:
                        process="p_%d.sh"%(n_run)
                        dummy_file = os.path.join(tmp_dir, process)
                        shutil.copy(os.path.join(scripts_path,config.get("Cluster","command_queue")),dummy_file)
                        out = open(dummy_file , "a")
                     out.write("%s\n"%(job))
                     n_jobs=n_jobs+1
                  out.close()
                  # Get current working directory #
                  cwd = os.getcwd()
                  # Change directory #
                  os.chdir(tmp_dir)
                  # Submit all jobs
                  for n in xrange(0,n_run):
                   sh_file="p_%d.sh"%(n)
                   if cluster_queue is None:
                    os.system("%s  %s" % (submit,sh_file))
                   else:
                    os.system("%s -q %s %s" % (submit,cluster_queue,sh_file))
                  # Return to original directory #
                  os.chdir(cwd)
    # Exit if stops here #
    if options.stop_step == 5:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ###################################
    # 6.  Parse tomtom files          #
    ###################################
    if options.start_step <= 6:
        # Load data into a dataframe #
        series_dict = {}
        if options.verbose:sys.stdout.write("\t-- Parsing TOMTOM folder\n")
        tomtom_lists=[]
        batch=0
        nbatch=0
        for tfile in os.listdir(os.path.join(os.path.abspath(options.output_dir), "tomtom")):
          if tfile.endswith("txt"):
           try:
            tomtom_obj = tomtom.Tomtom(open(os.path.join(os.path.abspath(options.output_dir), "tomtom", tfile), "r").readlines())
            if options.verbose:sys.stdout.write("\t\t-- File %s\n"%(tfile))
            tf_id      = tfile.split("_mL_")[0] 
            motif      = tfile.split("_mL_")[1]
            family_tf  = tfile.split("_mL_")[2]
            template   = tfile.split("_mL_")[3]
            potential  = tfile.split("_mL_")[4]
            family     = potential.split(".")[0]
            comput     = ".".join(potential.split(".")[1:])
            data       = tfile.split("_mL_")[5]
            distance   = tfile.split("_mL_")[6]
            threshold  = tfile.split("_mL_")[7]
            pdb_chain  = tfile.split("_mL_")[8].rstrip(".txt")
            if int(threshold) not in range(70, 101, 2):continue
            if "general" in family:
                fam_potentials = "general"
            else:
                fam_potentials = "family"
            if "taylor" in comput:
                taylor = True
            else:
                taylor = False
            computation=comput.split(".")[-1]
            pval=1.0
            rank=len(tomtom_obj.get_hits())
            for tomhit in tomtom_obj.get_hits():
                if tomhit.get_hit().endswith(".txt"):
                   motif_query = tomhit.get_hit().strip(".txt")
                else:
                   motif_query = tomhit.get_hit()
                if motif == motif_query:
                   pval = tomhit.get_p_value()
                   rank = tomhit.get_rank()
                   break
            tomtom_dict = {
                "tf_id": tf_id,
                "motif": motif,
                "family": family_tf,
                "fam_potentials": fam_potentials,
                "data": data,
                "computation": computation,
                "taylor": taylor,
                "distance": distance,
                "threshold": threshold,
                "p_value": pval,
                "rank":rank
            }
            tomtom_lists.append(tomtom_dict)
            #accumulate results for general parameters only when using general potential
            if "general" in family:
                fam_potentials = "general"
                family_tf      = "general"
                tomtom_dict = {
                 "tf_id": tf_id,
                 "motif": motif,
                 "family": family_tf,
                 "fam_potentials": fam_potentials,
                 "data": data,
                 "computation": computation,
                 "taylor": taylor,
                 "distance": distance,
                 "threshold": threshold,
                 "p_value": pval,
                 "rank":rank
                }
                tomtom_lists.append(tomtom_dict)
            nbatch=nbatch+1
            #if options.verbose:sys.stdout.write("\t\t\t-- batch results %d\n"%nbatch)
            if nbatch >= 10000:
                dump_file_batch=os.path.join(os.path.abspath(options.output_dir), "results","grid_search_df."+str(batch)+".pickle")
                if options.verbose:sys.stdout.write("\t\t-- dump batch results on %s\n"%dump_file_batch)
                fdump_file_batch=open(dump_file_batch,"wb")
                try:
                   cPickle.dump(tomtom_lists,fdump_file_batch)
                   fdump_file_batch.close()
                except:
                   if options.verbose:sys.stdout.write("\t\t-- Failed to DUMP results\n")
                nbatch=0
                batch=batch+1
                tomtom_lists=[]
           except:
            if options.verbose:sys.stdout.write("\t\t-- Skip file %s\n"%(tfile))
                 
        #Save dump
        if nbatch > 0:
            try:
                dump_file_batch=os.path.join(os.path.abspath(options.output_dir), "results","grid_search_df."+str(batch)+".pickle")
                if options.verbose:sys.stdout.write("\t\t-- dump batch results on %s\n"%dump_file_batch)
                fdump_file_batch=open(dump_file_batch,"wb")
                cPickle.dump(tomtom_lists,fdump_file_batch)
                fdump_file_batch.close()
            except:
                if options.verbose:sys.stdout.write("\t\t-- Failed to DUMP results\n")
    # Exit if stops here #
    if options.stop_step == 6:
        sys.stdout.write("Exiting...\n\n")
        exit(0)


    ###################################
    # 7.  Create outputs              #
    ###################################
    if options.start_step <= 7:
        # Read the csv data frame #
        #Read CPickle
        tomtom_lists=[]
        for tfile in os.listdir(os.path.join(os.path.abspath(options.output_dir), "results")):
            if tfile.startswith("grid_search_df"):
               if options.verbose:sys.stdout.write("\t-- reading %s\n"%tfile)
               tomtom_batch=cPickle.load( open(os.path.join(os.path.abspath(options.output_dir), "results", tfile),"rb") )
               tomtom_lists.extend(tomtom_batch)
        if options.verbose:sys.stdout.write("\t-- create PANDAS data frame\n")
        dataframe = pd.DataFrame(tomtom_lists)
        # Save dataframe as a textfile #
        if options.verbose:sys.stdout.write("\t-- writing data frame CSV file grid_search_df.csv\n")
        dataframe.to_csv(path_or_buf=os.path.join(os.path.abspath(options.output_dir), "results", "grid_search_df.csv"))
        """
        ### Columns of the transposed dataframe
        0 -> computation 
        1 -> data
        2 -> distance
        3 -> fam_potentials
        4 -> family
        5 -> motif
        6 -> p_value
        7 -> taylor
        8 -> tf_id
        9 -> threshold
       10 -> rank
        """
        families = list(set(dataframe.transpose().loc['family']))
        datas = list(set(dataframe.transpose().loc['data']))
        fam_potentials = list(set(dataframe.transpose().loc['fam_potentials']))
        taylors = list(set(dataframe.transpose().loc['taylor']))
        computations = list(set(dataframe.transpose().loc['computation']))
        distances = list(set(dataframe.transpose().loc['distance']))
        thresholds = sorted(list(set(dataframe.transpose().loc['threshold'])))  
        results_dict = {}
        min_results_dict = {}
        all_results_dict = {}
        for family in families:
            if family == "family": continue
            if options.verbose:sys.stdout.write("\t\t-- check family %s\n"%(family))
            max_positives = None
            max_average = None
            max_positives_rank = None
            min_average_rank = None
            max_p_list=set()
            max_pr_list=set()
            max_a_list=set()
            min_ar_list=set()
            for data in datas:
                if data == "data": continue
                for fam_potential in fam_potentials:
                    if fam_potential == "fam_potentials": continue
                    for taylor in taylors:
                        if taylor == "taylor": continue
                        for computation in computations:
                            if computation == "computation": continue
                            for distance in distances:
                                if distance == "distance": continue
                                for threshold in thresholds:
                                    if threshold == "threshold": continue
                                    # Get positives: number of tomtom outputs with p-value under 0.05 #
                                    if options.verbose:sys.stdout.write("\t\t\t-- condition %s\n"%(str(data)+" + Taylor:"+str(taylor)+" + Counting:"+str(computation)+" + Distance:"+str(distance)+" + PWM:"+str(threshold)))
				    p_value_set=set()
                                    rank_set=set()
                                    for tomtom_dict in tomtom_lists:
                                      if str(tomtom_dict['family'])!=str(family): continue
                                      if str(tomtom_dict['data'])!=str(data): continue
                                      if str(tomtom_dict['computation'])!=str(computation): continue
                                      if str(tomtom_dict['fam_potentials'])!=str(fam_potential): continue
                                      if str(tomtom_dict['taylor'])!=str(taylor): continue
                                      if abs(float(tomtom_dict['distance'])-float(distance))>0.1: continue
                                      if abs(float(tomtom_dict['threshold'])-float(threshold))>0.1: continue
                                      rank_set.add(int(tomtom_dict['rank']))
                                      pval=tomtom_dict['p_value']
                                      if type(pval) == type("a"): continue  
                                      if str(pval) == "nan": continue
                                      if not isinstance(pval,float) and not isinstance(pval,int):continue
                                      p_value_set.add(float(tomtom_dict['p_value']))
				    p_value_list=list(p_value_set)
                                    if len(p_value_list)<=0:continue
                                    ranking_list=list(rank_set)
                                    if len(ranking_list)<=0:continue
                                    # get positives by p_value
                                    positives = 0
                                    pval_list = []
                                    for pval in p_value_list:
                                        if pval<0: continue
                                        if pval<1.0e-10: pval=1.0e-10
                                        pval_list.append(numpy.log10(pval)*-1)
                                        if float(pval) < 0.05:
                                            positives += 1
                                    # Get the average p_value #
                                    if len(pval_list)>0:
                                       average = numpy.mean(pval_list)
                                    else:
                                       average=None

                                    # get positives by ranking motifs
                                    positives_rank = 0
                                    rank_list = []
                                    max_rank=max(ranking_list)
                                    for rank in ranking_list:
                                        rank_list.append(float(rank))
                                        if rank < 1+int(max_rank/100):
                                           positives_rank += 1
                                    # Get the average rank #
                                    if len(rank_list)>0:
                                       average_rank = numpy.mean(rank_list)
                                    else:
                                       average_rank=None

                                    # select best values and thresholds
                                    if (max_positives == None) and (max_average == None) and (max_positives_rank == None) and (min_average_rank == None):
                                        max_positives = positives
                                        max_average = average
                                        max_positives_rank = positives_rank
                                        min_average_rank = average_rank
                                        max_p_list.add(max_positives)
                                        max_pr_list.add(max_positives_rank)
                                        max_a_list.add(max_average)
                                        min_ar_list.add(min_average_rank)
                                        best_parameters = [str(data), str(fam_potential), str(taylor), str(computation), str(distance), str(threshold), str(max_positives), str(max_average), str(max_positives_rank), str(min_average_rank)]
                                        min_threshold = float(threshold)
                                        min_parameters = [str(data), str(fam_potential), str(taylor), str(computation), str(distance), str(threshold), str(max_positives), str(max_average), str(max_positives_rank), str(min_average_rank)]
                                    else:
                                        max_positives_o=max([x for x in max_p_list])
                                        max_positives_rank_o=max([x for x in max_pr_list])
                                        max_average_o=max([x for x in max_a_list])
                                        min_average_rank_o=min([x for x in min_ar_list])
                                        if positives > max_positives_o and positives_rank > max_positives_rank_o and average_rank < min_average_rank_o:
                                            max_positives = positives
                                            max_average = average
                                            max_positives_rank = positives_rank
                                            min_average_rank = average_rank
                                            best_parameters = [str(data), str(fam_potential), str(taylor), str(computation), str(distance), str(threshold), str(positives), str(average), str(positives_rank), str(average_rank)]
                                            min_threshold = float(threshold)
                                            min_parameters = [str(data), str(fam_potential), str(taylor), str(computation), str(distance), str(threshold), str(positives), str(average), str(positives_rank), str(average_rank)]
                                        if ((positives > max_positives_o * 0.95 ) and (average > max_average_o)):
                                            max_average = average
                                            max_positives = positives
                                            best_parameters = [str(data), str(fam_potential), str(taylor), str(computation), str(distance), str(float(threshold)), str(positives), str(average), str(positives_rank), str(average_rank)] 
                                        if ((positives_rank > max_positives_rank_o * 0.95 ) and (positives > max_positives_o * 0.85 ) and (average_rank < min_average_rank_o)):
                                            min_average_rank = average_rank
                                            max_positives_rank = positives_rank
                                            best_parameters = [str(data), str(fam_potential), str(taylor), str(computation), str(distance), str(float(threshold)), str(positives), str(average), str(positives_rank), str(average_rank)] 
                                        if ((positives > max_positives_o * 0.95 ) and (average > max_average_o * 0.95) and (float(threshold) < min_threshold)):
                                            min_threshold = float(threshold)
                                            max_positives = positives
                                            min_parameters = [str(data), str(fam_potential), str(taylor), str(computation), str(distance), str(threshold), str(positives), str(average), str(positives_rank), str(average_rank)]
                                        if ((positives_rank > max_positives_rank_o * 0.95) and (positives > max_positives_o * 0.85 ) and (average > max_average_o * 0.85) and ( average_rank < min_average_rank_o ) and (float(threshold) < min_threshold)):
                                            min_threshold = float(threshold)
                                            max_positives_rank = positives_rank
                                            min_parameters = [str(data), str(fam_potential), str(taylor), str(computation), str(distance), str(threshold), str(positives), str(average), str(positives_rank), str(average_rank)]
                                        max_p_list.add(max_positives)
                                        max_pr_list.add(max_positives_rank)
                                        max_a_list.add(max_average)
                                        min_ar_list.add(min_average_rank)
                                    # Get all the data into an output file, so we can make further analysis on it #
                                    all_results_dict.setdefault(family, [])
                                    all_results_dict[family].append([str(data), str(fam_potential), str(taylor), str(computation), str(distance), str(threshold), str(positives), str(average), str(positives_rank), str(average_rank)])
            results_dict[family] = best_parameters
            min_results_dict[family] = min_parameters
        # Store the results dict as a pickled dictionary #
        pickle.dump(results_dict, open(os.path.join(os.path.abspath(options.output_dir), "results", "results_dict.p"), "wb")) 
        pickle.dump(min_results_dict, open(os.path.join(os.path.abspath(options.output_dir), "results", "min_results_dict.p"), "wb")) 
        pickle.dump(all_results_dict, open(os.path.join(os.path.abspath(options.output_dir), "results", "all_data_dict.p"), "wb"))
        # Write a text file with the results of the grid search #
        out_file = os.path.join(os.path.abspath(options.output_dir), "results", "grid_search_results.txt")
        o_file = open(out_file, "w")
        o_file.write("Family,Data,Fam_potentials,Taylor,Computation,Distance,Threshold,Positives,Average,Positive_ranks,Average_rank\n")
        if results_dict.has_key("general"):
            o_file.write("general," + ",".join(results_dict["general"]) + "\n")
        for family in sorted(results_dict.keys()):
            if family != "general":
                o_file.write(family + "," + ",".join(results_dict[family]) + "\n")
        o_file.close()
        # Get results with highest number of positives and the lowest threshold #
        min_file = os.path.join(os.path.abspath(options.output_dir), "results", "grid_search_results_min_threshold.txt")
        m_file = open(min_file, "w")
        m_file.write("Family,Data,Fam_potentials,Taylor,Computation,Distance,Threshold,Positives,Average,Positive_ranks,Average_rank\n")
        if min_results_dict.has_key("general"):
            m_file.write("general," + ",".join(min_results_dict["general"]) + "\n")
        for family in sorted(min_results_dict.keys()):
            if family != "general":
                m_file.write(family + "," + ",".join(min_results_dict[family]) + "\n")
        m_file.close()
        # Get results with highest number of positives and the lowest threshold second version#
        min_file = os.path.join(os.path.abspath(options.output_dir), "results", "grid_search_results_best_min_threshold.txt")
        m_file = open(min_file, "w")
        m_file.write("Family,Data,Fam_potentials,Taylor,Computation,Distance,Threshold,Positives,Average,Positive_ranks,Average_rank\n")
        for family in sorted(all_results_dict.keys()):
            top_50 = sorted(all_results_dict[family], key=lambda x: float(x[6]), reverse=True)[:50]
            top_25 = sorted(top_50, key=lambda x: float(x[7]), reverse=True)[:25]
            top_10 = sorted(top_25, key=lambda x: float(x[8]), reverse=True)[:10]
            top_5 = sorted(top_10, key=lambda x: float(x[9]))[:5]
            top_1 = sorted(top_5, key=lambda x: float(x[5]))[0]
            m_file.write(family + "," + ",".join(top_1) + "\n")
        m_file.close()
        # Get results with highest number of positives and the lowest threshold third version#
        min_file = os.path.join(os.path.abspath(options.output_dir), "results", "grid_search_results_rank_min_threshold.txt")
        m_file = open(min_file, "w")
        m_file.write("Family,Data,Fam_potentials,Taylor,Computation,Distance,Threshold,Positives,Average,Positive_ranks,Average_rank\n")
        for family in sorted(all_results_dict.keys()):
            top_50 = sorted(all_results_dict[family], key=lambda x: float(x[8]), reverse=True)[:50]
            top_25 = sorted(top_50, key=lambda x: float(x[6]), reverse=True)[:25]
            top_10 = sorted(top_25, key=lambda x: float(x[8]), reverse=True)[:10]
            top_5 = sorted(top_10, key=lambda x: float(x[9]))[:5]
            top_1 = sorted(top_5, key=lambda x: float(x[5]))[0]
            m_file.write(family + "," + ",".join(top_1) + "\n")
        m_file.close()
        # Write a text file with all the results of the analysis #
        data_file = os.path.join(os.path.abspath(options.output_dir), "results", "grid_search_data.txt")
        d_file = open(data_file, "w")
        d_file.write("Family,Data,Fam_potentials,Taylor,Computation,Distance,Threshold,Positives,Average,Positive_ranks,Average_rank\n")
        for family in sorted(all_results_dict.keys()):
            for param in all_results_dict[family]:
                d_file.write(family + "," + ",".join(param) + "\n")
        d_file.close()
        # Write a top 10 dictionary #
        for family in sorted(all_results_dict.keys()):
            top_file = os.path.join(os.path.abspath(options.output_dir), "results", "grid_search_top10_"+family+".txt")
            t_file = open(top_file, "w")
            t_file.write("Family,Data,Fam_potentials,Taylor,Computation,Distance,Threshold,Positives,Average,Positive_ranks,Average_rank\n")
            top_10 = sorted(all_results_dict[family], key=lambda x: float(x[6]), reverse=True)[:10]
            for t in top_10:
                t_file.write(family + "," + ",".join(t) + "\n")
            t_file.close()
        # Write a top 10 ranking dictionary #
        for family in sorted(all_results_dict.keys()):
            top_file = os.path.join(os.path.abspath(options.output_dir), "results", "grid_search_top10_ranking_"+family+".txt")
            t_file = open(top_file, "w")
            t_file.write("Family,Data,Fam_potentials,Taylor,Computation,Distance,Threshold,Positives,Average,Positive_ranks,Average_rank\n")
            top_10 = sorted(all_results_dict[family], key=lambda x: float(x[8]), reverse=True)[:10]
            for t in top_10:
                t_file.write(family + "," + ",".join(t) + "\n")
            t_file.close()

    # Exit if stops here #
    if options.stop_step == 7:
        sys.stdout.write("Exiting...\n\n")
        exit(0)
    


    ###################################
    # 8.  Create configuration file   #
    ###################################
    if options.start_step <= 8:
        # Write a text file with the results of the grid search #
        if options.verbose:sys.stdout.write("Select parameters...\n")
        selected={}
        input_file = os.path.join(os.path.abspath(options.output_dir), "results", "grid_search_results.txt")
        i_file = open(input_file, "r")
        best_parameters={}
        for line in i_file:
            if line.startswith("Family"):continue
            family=line.split(",")[0]
            best_parameters.setdefault(family,[]).append(line.strip().split(",")[1:])
        i_file.close()
        bmin_parameters={}
        min_file = os.path.join(os.path.abspath(options.output_dir), "results", "grid_search_results_best_min_threshold.txt")
        m_file = open(min_file, "r")
        for line in m_file:
            if line.startswith("Family"):continue
            family=line.split(",")[0]
            bmin_parameters.setdefault(family,[]).append(line.strip().split(",")[1:])
        m_file.close()
        rmin_parameters={}
        min_file = os.path.join(os.path.abspath(options.output_dir), "results", "grid_search_results_rank_min_threshold.txt")
        m_file = open(min_file, "r")
        for line in m_file:
            if line.startswith("Family"):continue
            family=line.split(",")[0]
            rmin_parameters.setdefault(family,[]).append(line.strip().split(",")[1:])
        m_file.close()
        mmin_parameters={}
        min_file = os.path.join(os.path.abspath(options.output_dir), "results", "grid_search_results_min_threshold.txt")
        m_file = open(min_file, "r")
        for line in m_file:
            if line.startswith("Family"):continue
            family=line.split(",")[0]
            mmin_parameters.setdefault(family,[]).append(line.strip().split(",")[1:])
        m_file.close()

        topp_parameters={}
        topr_parameters={}
        top_parameters={}
        for family in best_parameters.iterkeys():
            if options.verbose:sys.stdout.write("\t--best parameters for %s...\n"%family)
            best_file=os.path.join(os.path.abspath(options.output_dir), "results","grid_search_top10_"+family+".txt")
            b_file=open(best_file,"r")
            for line in b_file:
              if line.startswith("Family"):continue
              family=line.split(",")[0]
              top_parameters.setdefault(family,[]).append(line.strip().split(",")[1:])
              topp_parameters.setdefault(family,[]).append(line.strip().split(",")[1:])
            b_file.close()
            best_file=os.path.join(os.path.abspath(options.output_dir), "results","grid_search_top10_ranking_"+family+".txt")
            b_file=open(best_file,"r")
            for line in b_file:
              if line.startswith("Family"):continue
              family=line.split(",")[0]
              top_parameters.setdefault(family,[]).append(line.strip().split(",")[1:])
              topr_parameters.setdefault(family,[]).append(line.strip().split(",")[1:])
            b_file.close()
        for family in best_parameters.iterkeys():
            if top_parameters.has_key(family):
               max_pp=max([float(x[6]) for x in topp_parameters.get(family)])
               max_pr=max([float(x[8]) for x in topr_parameters.get(family)])
               min_pp=min([float(x[6]) for x in topp_parameters.get(family)])
               min_pr=min([float(x[8]) for x in topr_parameters.get(family)])
            if best_parameters.has_key(family):
               b_pp=max([float(x[6]) for x in best_parameters.get(family)])
               b_pr=max([float(x[8]) for x in best_parameters.get(family)])
            if mmin_parameters.has_key(family):
               mm_pp=max([float(x[6]) for x in mmin_parameters.get(family)])
               mm_pr=max([float(x[8]) for x in mmin_parameters.get(family)])
            if bmin_parameters.has_key(family):
               bm_pp=max([float(x[6]) for x in bmin_parameters.get(family)])
               bm_pr=max([float(x[8]) for x in bmin_parameters.get(family)])
            if rmin_parameters.has_key(family):
               rm_pp=max([float(x[6]) for x in rmin_parameters.get(family)])
               rm_pr=max([float(x[8]) for x in rmin_parameters.get(family)])

            accept_best=False
            if max_pp>0 and min_pp>0:
               d_pp=max_pp-min_pp
               d_pr=max_pr-min_pr
               if (b_pp > min_pp - d_pp/2) and (b_pr >= min_pr-d_pr): accept_best=True
            else:
               accept_best=True

            accept_min =False
            accept_mmin =False
            accept_bmin =False
            accept_rmin =False
            min_parameters=[]
            if (max_pp>0 and min_pp>0) or (max_pr>0 and min_pr>0):
               d_pp=max_pp-min_pp
               d_pr=max_pr-min_pr
               if (mm_pp > min_pp - d_pp/2) and (mm_pr > min_pr - d_pr/2): accept_mmin =True
               if (bm_pp >= min_pp - d_pp/2) and (bm_pr >= min_pr - d_pr/2): accept_bmin =True
               if (rm_pp >= min_pp - d_pp/2) and (rm_pr >= min_pr - d_pr/2): accept_rmin =True
               if accept_mmin: min_parameters.append(mmin_parameters[family][0])
               if accept_bmin: min_parameters.append(bmin_parameters[family][0])
               if accept_rmin: min_parameters.append(rmin_parameters[family][0])
               dummy_min=[x for x in sorted(min_parameters, key=lambda x: float(x[5])) if x[0]=="pbm"]
               if len(dummy_min)>0:
                  min_parameters=dummy_min
               else:
                  min_parameters=sorted(min_parameters, key=lambda x: float(x[5]))
               if len( min_parameters )>0: accept_min =True
            else:
               min_parameters.append(mmin_parameters[family][0])
               min_parameters.append(bmin_parameters[family][0])
               min_parameters.append(rmin_parameters[family][0])
               dummy_min=[x for x in sorted(min_parameters, key=lambda x: float(x[5])) if x[0]=="pbm"]
               if len(dummy_min)>0:
                  min_parameters=dummy_min
               else:
                  min_parameters=sorted(min_parameters, key=lambda x: float(x[5]))
               if len( min_parameters )>0: accept_min =True
               
            if accept_min and accept_best:
                if min_parameters[0][0] == "pdb" and best_parameters[family][0][0] == "pbm": accept_min = False
            

            #select min if acceptable, otherwise accept best, otherwise select the minimum threhold among top
            if accept_min:
               if options.verbose:sys.stdout.write("\t\t-- Minimum threshold parameters for %s...\n"%family)
               selected.setdefault(family,min_parameters[0])
            elif accept_best:
               if options.verbose:sys.stdout.write("\t\t-- Best P-value parameters for %s...\n"%family)
               selected.setdefault(family,best_parameters[family][0])
            else:
               t=[x for x in sorted(top_parameters[family], key=lambda x: float(x[5])) if x[0]=="pbm"]
               b=[x for x in sorted(topp_parameters[family], key=lambda x: float(x[5])) if x[0]=="pbm"]
               r=[x for x in sorted(topr_parameters[family], key=lambda x: float(x[5])) if x[0]=="pbm"]
               if len(b)>0:
                  if options.verbose:sys.stdout.write("\t\t-- Best PBM P-value with lowest threshold parameters for %s...\n"%family)
                  selected.setdefault(family,b[0])
               elif len(r)>0:
                  if options.verbose:sys.stdout.write("\t\t-- Best PBM ranking with lowest threshold parameters for %s...\n"%family)
                  selected.setdefault(family,r[0])
               elif len(t)>0:
                  if options.verbose:sys.stdout.write("\t\t-- Lowest PBM threshold parameters among best rank and P-value for %s...\n"%family)
                  selected.setdefault(family,t[0])
               else:
                  if options.verbose:sys.stdout.write("\t\t-- Lowest threshold parameters among best rank and P-value for %s...\n"%family)
                  selected.setdefault(family,sorted(topp_parameters[family], key=lambda x: float(x[5]))[0])
        #Write in parameters file
        parameter_file=os.path.join(os.path.abspath(options.output_dir), "results","grid_search_parameters.txt")
        o_file=open(parameter_file,"w")
        for family,parameters in selected.iteritems():
            data          = parameters[0]
            fam_potential = parameters[1]
            if parameters[2] == "True": 
                taylor= "taylor"
            else:
                taylor = " "
            computation   = parameters[3]
            distance      = "%3d"%(int(float(parameters[4])))
            threshold     = "%3.2f"%(float(parameters[5])/100)
            o_file.write("%s = %s,%s,%s,%s,%s,%s\n"%(family.replace("=","/").replace("_"," "),data,fam_potential,taylor,computation,threshold,distance))
        o_file.close()
         


    # Exit if stops here #
    if options.stop_step == 8:
        sys.stdout.write("Exiting...\n\n")
        exit(0)
    

