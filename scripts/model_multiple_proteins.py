import os, sys, re
import ConfigParser
import optparse
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

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python model_protein.py -i input_file -p pdb_dir [--dummy=dummy_dir --n-model=n_model --n-total=n_total -o output_dir] [-a -d -f -m -e -r resolution_file -s --dimer --monomer --unbound_fragments --unrestrictive] --info LOG_FILE")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file with several sequences in FASTA format (sequence) or THREADING format (i.e. from threader.py)", metavar="{filename}")
    parser.add_option("--n-model", default=1, action="store", type="int", dest="n_model", help="Number of models per template (default = 1)", metavar="{int}")
    parser.add_option("--n-total", action="store", type="int", dest="n_total", default=None, help="Total number of models per execution of this program (if not \"None\", automatically sets the \"--n-model\" parameter; supersedes option \"--n-model\"; default = None)", metavar="{int}")
    parser.add_option("--best", action="store_true", dest="best", default=False, help="It only uses the first hit of blast (supersedes n-model and n-total to 1; default = False)", metavar="{boolean}")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="{directory}")
    parser.add_option("-t", "--threading", default=False, action="store_true", dest="threading", help="Use a list of threading files (default=False)", metavar="{boolean}")
    parser.add_option("-p", "--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (i.e. output dir from pdb.py)", metavar="{directory}")
    parser.add_option("-l", action="store", type="string", dest="label", help="Label to include in the output models name", metavar="{str}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode. If not selected the dummy directory will be removed (default = False)", metavar="{boolean}")
    parser.add_option("--force", default=False, action="store_true", dest="force", help="Force to do the modelling even if the file already exists (default = False)", metavar="{boolean}")
    parser.add_option("--opt", default=False, action="store_true", dest="optimization", help="Run modelling with optimization (default = False). Note: use it with care, this is not advisable as structures can take missconfomations", metavar="{boolean}")
    parser.add_option("--parallel", default=False, action="store_true", dest="parallel", help="Run in parallel if the input is a directory (default = False)")
    parser.add_option("--info",default=None,action="store", type="string", dest="info",help="Information LOG file of MODELS that have failed and have been completed")
    parser.add_option("--reuse",default=False, action="store_true", dest="reuse", help="Reuse the information files. If the flag is used then models that had failed will remain as FAILED, otherwise it tries to redo them (default=False)")


    group = optparse.OptionGroup(parser, "FASTA options", "Provide a FASTA file as input.")
    group.add_option("-a", "--all", default=False, action="store_true", dest="all_models", help="Use ALL possible templates (default = False)", metavar="{boolean}")
    group.add_option("-d", "--dna", default=False, action="store_true", dest="dna", help="DNA modelling (if \"True\", remodels DNA using 3DNA; default = False)", metavar="{boolean}")
    group.add_option("-f","--full", default=False, action="store_true", dest="full_mode", help="Model all protein fragments with template that bind DNA (default = False)", metavar="{boolean}")
    group.add_option("--unbound_fragments", default=False, action="store_true", dest="full_fragments_mode", help="Model all protein fragments with template even unbound to DNA (default = False)", metavar="{boolean}")
    group.add_option("-e", "--twilight", default=False, action="store_true", dest="filter_twilight_zone_hits", help="Select template-hits over the twilight zone unless there are no models (default = False)", metavar="{boolean}")
    group.add_option("--restrictive", default=False, action="store_true", dest="filter_restrictive", help="Restrict modelling to only template-hits over the twilight zone when option --twilight is selected (default = False)", metavar="{boolean}")
    group.add_option("--unrestrictive", default=False, action="store_true", dest="no_restrict", help="Allow models with wrong alignment in the DNA-binding interface (default = False)", metavar="{boolean}")
    group.add_option("-m", action="store", dest="mutation", default=None, help="Mutate input sequence (e.g. \"R88C\" mutates the \"R\" at position \"88\" for a \"C\")", metavar="{str}")
    group.add_option("-r", action="store", default=None, dest="resolution_file", help="Resolution file (if provided, sorts similar alignments by RMSD of the template PDB; default = None)", metavar="{filename}")
    group.add_option("-s","--select", default=False, action="store_true", dest="select_mode", help="Select mode (if True, selects either monomers or dimers, whichever has more hits; default = False)", metavar="{boolean}")
    group.add_option("--dimers", default=False, action="store_true", dest="dimer_mode", help="If True forces to work only with dimers; default = False)", metavar="{boolean}")
    group.add_option("--monomers", default=False, action="store_true", dest="monomer_mode", help="If True forces to work only with monomers; default = False)", metavar="{boolean}")
    group.add_option("--renumerate", default=False, dest = 'renumerate' , action = 'store_true', help = 'Flag to renumber the sequences as in the original FastA input (default is False)')
    group.add_option("--chains_fixed",default=False, action="store_true", dest="chains_fixed", help="Replace the names of the protein chains to A-B-C-D... and DNA chains to a-b-c-d....")


    parser.add_option_group(group)

    (options, args) = parser.parse_args()

    if options.input_file is None  or options.pdb_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")
    if options.mutation is not None:
        m = re.search("^([a-zA-Z]\d+[a-zA-Z])$", options.mutation)
        if not m:
            parser.error("incorrect provided mutation \"%s\": e.g.\"R88C\" mutates the \"R\" at position 88 for a \"C\"" % options.mutation)

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
    threading = options.threading

    #collect data
    dummy_dir = options.dummy_dir
    if not dummy_dir.startswith("/"): dummy_dir = os.path.abspath(options.dummy_dir)
    if options.input_file  is not None:
     input_file = options.input_file
     if not input_file.startswith("/"): input_file = os.path.abspath(options.input_file)
    output_dir = options.output_dir
    if not output_dir.startswith("/"): output_dir = os.path.abspath(options.output_dir)
    pdb_dir = options.pdb_dir
    if not pdb_dir.startswith("/"): pdb_dir = os.path.abspath(options.pdb_dir)

    #Check input requirements
    if not functions.fileExist(input_file) or not os.path.exists(pdb_dir):
       sys.stdout.write("Missing input files, please check help with option -h\n")
       exit(0)

    # Create output and dummy directory #
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(dummy_dir):
        os.makedirs(dummy_dir)

    # Iterate untill all protein sequences are modelled
    submitted=set()
    n_done=0
    if options.verbose: print("Start iteration to model")
    info_file=options.info
    if info_file is None: info_file=os.path.join(options.output_dir,"modelling_list_execution.log")
    if reuse and functions.fileExist(info_file):
         if options.verbose: print ("Reuse previous information of runs from file %s"%info_file)
    else:
         if options.verbose: print ("Open to write %s"%info_file)
         log_file = open(info_file,"w")
         log_file.write("#List of MODELS\n")
         log_file.close()

    if threading:
      set_of_proteins=set()
      thread_input = open(input_file,"r")
      for thread in thread_input:
          if thread.startswith("#"): continue
          set_of_proteins.add(thread.strip())
    else:
      set_of_proteins=set()
      for name,sequence in functions.parse_fasta_file(input_file):
        protein = options.label+"_"+name.replace("|","_").split()[0]
        input_fasta=open(os.path.join(dummy_dir,protein+".fa"),"w")
        input_fasta.write(">%s\n%s\n"%(protein,sequence))
        input_fasta.close()
        set_of_proteins.add(protein)


    done    = functions.done_jobs(info_file)
    iterate = functions.check_done(done,set_of_proteins)
    n_done  = 0
    start_time = float(time.time())

    while(iterate):
        for protein in set_of_proteins:
            if protein in submitted:continue
            submitted.add(protein)
            if verbose: print("Modelling %s ..."%protein)
            if reuse:
              log_file = open(info_file,"r")
              skip_adding=False
              for line in log_file:
                       if protein in line.split(): skip_adding=True
              if skip_adding and verbose: print("\t-- Already in the information file ")
              if skip_adding: submitted.add(protein)
              if skip_adding: continue
              log_file.close()
            #paremeter values
            if threading:
              thread_file= protein
              label = os.path.basename(protein).rstrip(".txt")
              parameters =   " -i %s "%thread_file
            else:
              fasta_file =  os.path.join(dummy_dir,protein+".fa")
              label      =  protein
              parameters =   " -i %s "%fasta_file
            n_model    =  str(options.n_model)
            parameters =  parameters + " -o %s "%output_dir
            parameters =  parameters + " -p %s "%pdb_dir
            parameters =  parameters + " -l %s "%label
            parameters =  parameters + " --n-model %s "%n_model
            parameters =  parameters + " --info %s "%info_file
            parameters =  parameters + " --dummy %s "%dummy_dir
            if threading: parameters = parameters + " --threading "
            #None default
            if options.n_total is not None:          parameters = parameters + " --n-total %s "%options.n_total
            if options.mutation  is not None:        parameters = parameters + " -m %s "%options.mutation
            if options.resolution_file  is not None: parameters = parameters + " -r %s "%options.resolution_file
            #boolean
            if options.best :                     parameters = parameters + " --best "
            if options.verbose:                   parameters = parameters + " --verbose "
            if options.force:                     parameters = parameters + " --force "
            if options.optimization:              parameters = parameters + " --opt "
            if options.all_models:                parameters = parameters + " --all "
            if options.dna:                       parameters = parameters + " --dna "
            if options.full_mode:                 parameters = parameters + " --full "
            if options.full_fragments_mode:       parameters = parameters + " --unbound_fragments "
            if options.filter_twilight_zone_hits: parameters = parameters + " --twilight "
            if options.filter_restrictive:        parameters = parameters + " --restrictive "
            if options.no_restrict:               parameters = parameters + " --unrestrictive "
            if options.select_mode:               parameters = parameters + " --select "
            if options.dimer_mode:                parameters = parameters + " --dimer "
            if options.monomer_mode:              parameters = parameters + " --momnomer "
            if options.renumerate:                parameters = parameters + " --renumerate "
            if options.chains_fixed:              parameters = parameters + " --chains_fixed "
            program=os.path.join(scripts_path,"model_protein.py")
            if parallel:
              if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
              else: cluster_queue=config.get("Cluster", "cluster_queue")
              python = os.path.join(config.get("Paths", "python_path"), "python")
              if verbose: print("\t-- Submit  %s %s" % (program,parameters))
              functions.submit_command_to_queue("%s %s %s" % (python,program,parameters), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
              submitted.add(protein)      
            else:
              if verbose: print("\t-- Execute  %s %s" % (program,parameters))
              os.system("%s %s %s" % (python,program,parameters))
              submitted.add(protein)      
        #Check next iteration, profiles submitted and profiles done
        done    = functions.done_jobs(info_file)
        iterate = functions.check_done(done,set_of_proteins) 
        current_time = float(time.time())
        # if more than 1.5 days running stop iterating
        if current_time - start_time > 125000: 
            iterate=False
            sys.stdout.write("\t-- Exceeding time (limited to 1.5 days)\n")

        if len(done) > n_done:
           n_done=len(done)
           if options.verbose: 
                sys.stdout.write("Number of models already done %d\n"%n_done)
                sys.stdout.write("\t-- Check files done ...\n")
                log_file = open(info_file,"r")
                for line in log_file:
                   print("\t\t-- %s"%line.strip())
                log_file.flush()
                log_file.close()
                sys.stdout.write("\t-- Still running protein profiles  %s ...\n"%functions.check_done(done,set_of_proteins))
                sys.stdout.write("\t-- Continue iteration %s ...\n"%iterate)
                sys.stdout.flush()


    #Done
    shutil.rmtree(dummy_dir)
    if verbose:sys.stdout.write("Done\n")



    




