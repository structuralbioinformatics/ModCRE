import os, sys, re
import ConfigParser
import json
import numpy as np
import optparse
import shutil
import subprocess
from time import time
import pickle
import cPickle



#-------------#
# Options     #
#-------------#

def parse_options():
    '''
    This function parses the command line arguments and returns an optparse object.
    '''

    parser = optparse.OptionParser("scan.py [--dummy=DUMMY_DIR] -i INPUT_FILE [-l LABEL -o OUTPUT_DIR] --pbm=PBM_dir --pdb=PDB_DIR -s SPECIE [-v]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="DUMMY_DIR")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input FASTA file", metavar="INPUT_FILE")
    parser.add_option("-l", action="store", type="string", dest="label", help="Label ID (this identifies both output directory and files, if necessary; default = None)", metavar="LABEL")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="OUTPUT_DIR")
    parser.add_option("--pbm", action="store", type="string", dest="pbm_dir", help="PBM directory (i.e. output dir from pbm.py)", metavar="PBM_DIR")
    parser.add_option("--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (i.e. output dir from pdb.py)", metavar="PDB_DIR")
    parser.add_option("-s","--specie", action="store", type="string", dest="specie",  default=None, help="Specie to obtain specific orthologs (i.e. taxon/code/common_name as 9606/HUMAN/'Homo sapiens',default = 9606)", metavar="SPECIE")
    parser.add_option("--scan_family", action="store", type="string", dest="scan_family",  default=None, help="Scan with specific familes of TFs, separated by ',' (i.e. from pbm/families.txt, default uses all)", metavar="SCAN_FAMILY")
    parser.add_option("--ft", action="store", type="float", dest="fimo_pvalue_threshold", default=None, help="P-value threhold for fimo matches", metavar="PVALUE")
    parser.add_option("--max", action="store", type="int",dest="max_stored_matches",  default=None, help="Maximum number of matches stored", metavar="{integer}")
    parser.add_option("--db", default=None, action="store", type="string", dest="database", help="Database of PWMs that will be used to scan the DNA", metavar="DATABASE")
    parser.add_option("--external", default=None, action="store", type="string", dest="external", help="Database folder and accumulated MEME database of PWMs, its associated PDB and a directory with homologs (files are separated by coma, e.g.: '--external database_pwm_folder, database_of_pwm.txt, association_file.txt, homologs_folder' the first three are mandatory, but if homologs_folder is skipped or None, then the folder 'homologs' from pbm is used)", metavar="EXTERNAL_DB")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")
    parser.add_option("-c", "--complexes", default=False, action="store_true", dest="cluster_complexes", help="Cluster complexes into connected binding sites (default = False)")
    parser.add_option("--parallel", default=False, action="store_true", dest="parallel", help="Run in parallel (default = False)")
    parser.add_option("--info",default=None,action="store", type="string", dest="info",help="Information LOG file of SCANS that have failed and have been completed")
    parser.add_option("--reuse",default=False, action="store_true", dest="reuse", help="Reuse the information files. If the flag is used then scans that had failed will remain as FAILED, otherwise it tries to redo them (default=False)")
    parser.add_option("--rank",default=False, action="store_true", dest="rank", help="Rank by binding scoring-energy. If the flag is not used then all orthologs get anull score and the programs runs faster, otherwise the sacn is slow but the selection of orthologs is improved (default=False)")



    group = optparse.OptionGroup(parser, "Statistical potentials", "Select your statistical potentials of choice. In \"--auto\" mode, the program uses S3DCdd family potentials derived from both PDB and PBM data and approached by Taylor.\"-a\" option overrides options \"-f\", \"-p\", \"-m\", \"-r\" and \"-t\".")
    group.add_option("-a", "--auto", default=False, action="store_true", dest="auto_mode", help="Automate the selection of statistical potentials (default = False)")
    group.add_option("-f", "--family", default=False, action="store_true", dest="family_potentials", help="Use family potentials (default = False)")
    group.add_option("-p", default=False, action="store_true", dest="pbm_potentials", help="Use potentials derived from both PBM + PDB data (default = False)")
    group.add_option("-m", "--pmf", default=False, action="store_true", dest="pmf_potentials", help="Use of raw mean-force potentials with no Z-scoring (default = False)")
    group.add_option("--pot", default=None, action="store", type="string", dest="split_potential", help="Split-potential to be used (3d, 3dc, s3dc, s3dc_dd, s3dc_di, pair; default from configuration file)", metavar="{string}")
    group.add_option("--taylor", default=False, action="store_true", dest="taylor_approach", help="Approach PMF by Taylor (default = False)")
    group.add_option("-t", action="store", type="float", default=None, dest="score_threshold", help="Threshold on the scaled score to consider positive k-mers (default = 0.95)", metavar="{float}")
    group.add_option("-b", "--bins", default=False, action="store_true",  dest="bins", help="Computate the potentials: by bins (if selected) or accumulative (default).", metavar="{boolean}")
    group.add_option("--file",default=None, action="store", type="string", dest="potential_file", help="Use potentials from specific file (default = None)", metavar="{string}")
    group.add_option("--radius",default=0, action="store", type="float", dest="radius", help="Maximum contact distance to calculate interactions (default=0 implies the use of 'max_contact_distance' or family-specific radius from configuration", metavar="{string}")

   
  

    (options, args) = parser.parse_args()

    if options.input_file is None or options.pbm_dir is None or options.pdb_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

#-------------#
# Main        #
#-------------#

# Add "." to sys.path #
scripts_path =  os.path.abspath(os.path.dirname(__file__))
sys.path.append(scripts_path)


# Imports jbonet's module #
from SBI.structure import PDB
from SBI.data      import aminoacids1to3

# Imports my functions #
import functions,fimo, tomtom, interface, x3dna, contacts, triads,  model_protein, nr, threader,dssp
import pwm_pbm as PWM
import homologs as HOMO


# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Define which python to be used #
python = os.path.join(config.get("Paths", "python_path"), "python")


if __name__ == "__main__":

    start_time = time()
    # Arguments & Options #
    options = parse_options()
    dummy_dir             = os.path.abspath(options.dummy_dir)
    input_file            = os.path.abspath(options.input_file)
    label                 = options.label
    output_dir            = os.path.abspath(options.output_dir)
    pbm_dir               = os.path.abspath(options.pbm_dir)
    pdb_dir               = os.path.abspath(options.pdb_dir)
    specie                = options.specie
    verbose               = options.verbose
    scan_family           = options.scan_family
    fimo_pvalue_threshold = options.fimo_pvalue_threshold
    max_stored_matches    = options.max_stored_matches
    database              = options.database
    external              = options.external
    cluster_complexes     = options.cluster_complexes
    reuse                 = options.reuse
    rank                  = options.rank  
    info_file             = options.info
    auto_mode             = options.auto_mode
    family_potentials     = options.family_potentials  
    pbm_potentials        = options.pbm_potentials
    pmf_potentials        = options.pmf_potentials
    split_potential       = options.split_potential
    taylor_approach       = options.taylor_approach
    score_threshold       = options.score_threshold
    bins                  = options.bins
    potential_file        = options.potential_file
    radius                = options.radius
    parallel              = options.parallel

    if not os.path.exists(dummy_dir):  os.makedirs(dummy_dir)
    if not os.path.exists(output_dir): os.mkdir(output_dir)
    if info_file is None: info_file=os.path.join(output_dir,"scanner.log")
    window  = int(config.get("Parameters","window_scanner"))
    overlap = int(config.get("Parameters","overlap_scanner"))
    step    = max( (window - overlap),2)
    

    #Split the DNA sequence in windows
    if verbose: print("SCANNER: Parse original DNA fasta file...")
    # For each header, sequence... #
    for fasta_sequence in functions.parse_fasta_file(input_file):
        dna_header   = fasta_sequence[0]
        dna_sequence = fasta_sequence[1]
    
    number_of_fragments = int(len(dna_sequence)/step)
    dna_files = []
    dna_fragments = {}
    end   = 0
    extension=".txt"
    if input_file.endswith(".txt"):   extension = ".txt"
    if input_file.endswith(".fa"):    extension = ".fa"
    if input_file.endswith(".fasta"): extension = ".fasta"

    for i in range(number_of_fragments):
        start = max((end - overlap),0)
        if start >= len(dna_sequence):continue
        end   = min((start + window),len(dna_sequence))
        words_header=dna_header.split()
        words_header[0] = words_header[0]+":"+str(i+1)+":"+str(start+1)+"-"+str(end)
        f_dna_header   = " ".join(words_header)
        f_dna_sequence = dna_sequence[start:end]
        dna_file = input_file.rstrip(extension)+":"+str(i+1)+":"+str(start+1)+"-"+str(end)+extension
        fo = open(dna_file,"w")
        fo.write(">%s\n%s\n"%(f_dna_header,f_dna_sequence))
        fo.close()
        dna_files.append(dna_file)
        dna_fragments.setdefault(start+1,dna_file)

    dna_set=set(dna_files)
    submitted=set()
    n_done=0
    if options.verbose: print("Start iteration to scan")
    if reuse and functions.fileExist(info_file):
         if options.verbose: print ("Reuse previous information of runs from file %s"%info_file)
    else:
         if options.verbose: print ("Open to write %s"%info_file)
         log_file = open(info_file,"w")
         log_file.write("#List of DNA framents\n")
         log_file.close()


    done    = functions.done_jobs(info_file)
    iterate = functions.check_done(done,dna_set)
    n_done  = 0


    while(iterate):
        for start_index in dna_fragments.keys():
            dna_file=dna_fragments[start_index]
            dna=os.path.basename(dna_file)
            if dna in submitted:continue
            submitted.add(dna)
            if verbose: print("Scanning %s ..."%dna)
            if reuse:
              log_file = open(info_file,"r")
              skip_adding=False
              for line in log_file:
                       if dna in line.split(): skip_adding=True
              if skip_adding and verbose: print("\t-- Already in the information file ")
              if skip_adding: submitted.add(dna)
              if skip_adding: continue
              log_file.close()
            #paremeter values
            parameters= " -i " + dna_file
            parameters= parameters + " --index " + str(start_index)
            parameters= parameters + " --pbm " + pbm_dir
            parameters= parameters + " --pdb " + pdb_dir
            parameters= parameters + " -o " + output_dir+"_"+str(start_index)
            parameters= parameters + " --dummy " + dummy_dir
            parameters= parameters + " --info " + info_file
            if rank:                       parameters = parameters + " --rank "
            if bins :                      parameters = parameters + " --bins "
            if verbose:                    parameters = parameters + " --verbose "
            if auto_mode:                  parameters = parameters + " --auto "
            if pbm_potentials:             parameters = parameters + " -p "
            if pmf_potentials:             parameters = parameters + " --pmf "
            if taylor_approach:            parameters = parameters + " --taylor "
            if family_potentials: parameters =    parameters + " --family "
            if score_threshold is not None:       parameters = parameters + " -t " + str(score_threshold)
            if potential_file is not None:        parameters = parameters + " --file " + potential_file
            if radius is not None:                parameters = parameters + " --radius " + str(radius)
            if split_potential is not None:       parameters = parameters + " --pot " + split_potential
            if max_stored_matches is not None:    parameters = parameters + " --max " + str(max_stored_matches)
            if fimo_pvalue_threshold is not None: parameters = parameters + " --ft " + str(fimo_pvalue_threshold)
            if database is not None:              parameters = parameters + " --db " + database
            if external is not None:              parameters = parameters + " --external " + external
            if scan_family is not None:           parameters = parameters + " --scan_family " + scan_family
            if specie is not None:                parameters = parameters + " --specie "+ specie
            if label is not None:                 parameters = parameters + " -l " + label

            if parallel:
              if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
              else: cluster_queue=config.get("Cluster", "cluster_queue")
              program=os.path.join(scripts_path,"scan.py")
              python = os.path.join(config.get("Paths", "python_path"), "python")
              if verbose: print("\t-- Submit  %s %s" % (program,parameters))
              functions.submit_command_to_queue("%s %s %s" % (python,program,parameters), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
              submitted.add(dna)      
            else:
              os.system("%s %s %s" % (python,program,parameters))
              submitted.add(dna)      
            
        #Check next iteration, profiles submitted and profiles done
        done    = functions.done_jobs(info_file)
        iterate = functions.check_done(done,dna_set) 
        if len(done) > n_done:
           n_done=len(done)
           if options.verbose: 
                sys.stdout.write("Number of DNA fragments already done %d\n"%n_done)
                sys.stdout.write("\t-- Check files done ...\n")
                log_file = open(info_file,"r")
                for line in log_file:
                   print("\t\t-- %s"%line.strip())
                log_file.flush()
                log_file.close()
                sys.stdout.write("\t-- Still running scans  %s ...\n"%functions.check_done(done,dna_set))
                sys.stdout.write("\t-- Continue iteration %s ...\n"%iterate)
                sys.stdout.flush()

    # Collect output files
    # 1) threading files
    thread2index={}
    if not os.path.exists(os.path.join(output_dir,"aux_files")): os.mkdir(os.path.join(output_dir,"aux_files"))
    for start_index in dna_fragments.keys():
           folder = os.path.join(output_dir+"_"+str(start_index),"aux_files")
           for thread_file in os.listdir(folder):
               thread2index.setdefault(os.path.basename(thread_file),start_index)
               shutil.copy(os.path.join(output_dir+"_"+str(start_index),"aux_files",thread_file),os.path.join(output_dir,"aux_files",thread_file))
    # 2) orthologs
    if os.path.exists(os.path.join(output_dir,  "orthologs.json")) and verbose: sys.stdout.write("\t--Use existing 'orthologs.json' file\n")
    if not os.path.exists(os.path.join(output_dir,  "orthologs.json")):
       orthologs=[]
       for start_index in dna_fragments.keys():
           folder      = output_dir+"_"+str(start_index)
           f_orthologs = json.loads(''.join([line for line in functions.parse_file(os.path.join(folder, "orthologs.json"))]))
           orthologs.extend(f_orthologs)
       out = open(os.path.join(output_dir, "orthologs.json"), "wt")
       out.write(json.dumps(orthologs, separators=(',', ':'), indent=2))
       out.close()
    if os.path.exists(os.path.join(output_dir,  "orthologs.txt")) and verbose: sys.stdout.write("\t--Use existing 'orthologs.txt' file\n")
    if not os.path.exists(os.path.join(output_dir,  "orthologs.txt")):
       out = open(os.path.join(output_dir, "orthologs.txt"), "w")
       n=0
       for start_index in dna_fragments.keys():
           folder      = output_dir+"_"+str(start_index)
           t_orthologs = open(os.path.join(folder, "orthologs.txt"))
           for line in t_orthologs:
               if line.startswith(">"):
                   n = n+1
                   data = line.rstrip().split("|")
                   info = "|".join(data[1:])
                   out.write(">%s|%s\n"%(str(n),info))
                   continue
               if line.startswith("/"):
                   out.write("//\n")
                   continue
               thread = line.strip().split("\t")
               if str(thread[0]) == "None" and str(thread[1])=="None": 
                   continue
               else:
                   out.write("%s\t%s\n"%(str(thread[0]),str(thread[1])))
       out.close()
    # 3) orthologs with best templates
    if os.path.exists(os.path.join(output_dir,  "orthologs_with_best_templates.json")) and verbose: sys.stdout.write("\t--Use existing 'orthologs_with_best_templates.json' file\n")
    if not os.path.exists(os.path.join(output_dir,  "orthologs_with_best_templates.json")):
       orthologs=[]
       for start_index in dna_fragments.keys():
           folder      = output_dir+"_"+str(start_index)
           f_orthologs = json.loads(''.join([line for line in functions.parse_file(os.path.join(folder, "orthologs_with_best_templates.json"))]))
           orthologs.extend(f_orthologs)
       out = open(os.path.join(output_dir, "orthologs_with_best_templates.json"), "wt")
       out.write(json.dumps(orthologs, separators=(',', ':'), indent=2))
       out.close()
    if os.path.exists(os.path.join(output_dir,  "orthologs_with_best_templates.txt")) and verbose: sys.stdout.write("\t--Use existing 'orthologs_with_best_templates.txt' file\n")
    if not os.path.exists(os.path.join(output_dir,  "orthologs_with_best_templates.txt")):
       out = open(os.path.join(output_dir, "orthologs_with_best_templates.txt"), "w")
       n=0
       for start_index in dna_fragments.keys():
           folder      = output_dir+"_"+str(start_index)
           t_orthologs = open(os.path.join(folder, "orthologs_with_best_templates.txt"))
           for line in t_orthologs:
               if line.startswith(">"):
                   n = n+1
                   data = line.rstrip().split("|")
                   info = "|".join(data[1:])
                   out.write(">%s|%s\n"%(str(n),info))
                   continue
               if line.startswith("/"):
                   out.write("//\n")
                   continue
               thread = line.strip().split("\t")
               if str(thread[0]) == "None" and str(thread[1])=="None": 
                   continue
               else:
                   out.write("%s\t%s\n"%(str(thread[0]),str(thread[1])))
       out.close()

    #######################################################################
    # Check to continue or finish
    #######################################################################

    if not cluster_complexes:
       print("No clustering required")
       print("Done")
       print("total process time: %s seconds" % str(time() - start_time))
       exit(0)
    else:
       print("\n######################### CONTINUE WITH THE CLUSTERING OF CLOSE BINDING COMPLEXES #################################\n")


    #Initialize parameters to connect complexes and scan

    max_dimer_distance=int(config.get("Parameters", "max_dimer_distance")) 
    max_position_scan =int(config.get("Parameters", "max_position_scan"))
    max_complex_connect = int(config.get("Parameters", "max_complex_connect"))
    meme_path = config.get("Paths", "meme_path")
    max_ratio_loss_scan = config.get("Parameters", "max_ratio_loss_scan")

    ########################################################################
    # 1.      Select EXTENDED database and scan DNA sequence               #
    ########################################################################
    if verbose: sys.stdout.write("--Scan DNA sequence with EXTENDED database of PWMs...\n")
    # Execute fimo #
    if options.database is None and external is None:
        # Use the default databases #
        database_extended_file = os.path.join(pbm_dir, "pwms", "database_extended.txt")
        homologs_dir  = os.path.join(pbm_dir,"homologs") 
        if not functions.fileExist(database_extended_file):
            if verbose: sys.stdout.write("\t--DBD database extended not found, we generate it with pwms\n")  
            for pwm_file in os.listdir(os.path.join(pbm_dir, "pwms")):
                if pwm_file.endswith("meme"): 
                   m = re.search("(\S+)_(\d+)(\S+).meme", pwm_file)
                   if m:
                      os.system("cat %s >> %s" % (os.path.join(pbm_dir, "pwms",pwm_file), database_extended_file))
        else:
            if verbose: sys.stdout.write("\t--Database extended for all complexes is already available as 'database_extended.txt'...\n")
    elif external is None:
        # Use the database choosen by the user #
        database_extended_file = os.path.abspath(options.database)
        # By default homologs directory is that from pbm directory
        homologs_dir  = os.path.join(pbm_dir,"homologs") 
    else:
        # External data is originally None
        database_file = association_file = None
        # By default homologs directory is that from pbm directory
        homologs_dir  = os.path.join(pbm_dir,"homologs") 
        # read external data files
        data_external = external.split(",")
        if len(data_external) == 4: database_dir, database_file, association_file, homologs_dir = data_external
        if len(data_external) == 3: database_dir, database_file, association_file = data_external
        if len(data_external) < 3 : 
            print("Error: missing external data, at least 1 folder and 2 files are required")
            info=open(info_file,"a")
            info.write("%s\tFAIL\n"%(os.path.basename(input_file)))
            info.close()
            exit(0)
        if not database_file.startswith("/"):      database_file    = os.path.abspath(database_file)
        if not database_dir.startswith("/"):       database_dir     = os.path.abspath(database_dir)
        if not association_file.startswith("/"):   association_file = os.path.abspath(association_file)
        if not homologs_dir.startswith("/"):       homologs_dir     = os.path.abspath(homologs_dir)
        fold2pwm={}
        for pwm_file in os.listdir(os.path.join(pbm_dir, "pwms")):
            if not pwm_file.endswith(".meme"): continue
            if not pwm_file.endswith(".meme.s"): continue
            hit = pwm_file.rstrip(".meme").rstrip(".meme.s")
            # Get chain PDB obj #
            m = re.search("(\S+)_(\d+)(\S+)",hit)
            if m:
                  pdb_name  = m.group(1)
                  helix_hit = m.group(2)
                  pdb_chain = m.group(3)
            else:
                  pdb_name  = hit[:4]
                  helix_hit = hit[5:6]
                  pdb_chain = hit[6:]
            fold_hit  = pdb_name+"_"+pdb_chain
            fold2pwm.setdefault(fold_hit,set()).add(hit)
        associated_pdb={}
        associated_pdb_inv={}
        associated_helix={}
        assoc=open(association_file,"r")
        for line in assoc:
            data=line.split()
            associated_pdb.setdefault(data[0],data[1])
            associated_pdb_inv.setdefault(data[1],set()).add(data[0])
            if fold2pwm.has_key(data[1]):
               associated_helix.setdefault(data[0],str(sorted([int(pdb_helix_chain[5:6]) for pdb_helix_chain in fold2pwm.get(data[1])])[0]))
            else:
               associated_helix.setdefault(data[0],"1")
        assoc.close()
        database_extended_file = database_file

    try:
        fimo_ext = fimo.get_fimo_obj(database_extended_file, input_file,fimo_pvalue_threshold,max_stored_matches,os.path.abspath(options.dummy_dir))
    except:
        if verbose == True: sys.stdout.write("\n")
        sys.stdout.write("\"fimo\" execution failed!")
        if verbose == True: sys.stdout.write("\n\n")
        exit(0)
    print("total process time: %s seconds" % str(time() - start_time))


    ########################################################################
    #  2.     Get FIMO EXTENDED OBJECT and group the clusters              #
    ########################################################################
    if verbose: sys.stdout.write("--TF ensembles described by PWMs of complexes...\n")
    if os.path.exists(os.path.join(output_dir,  "orthologs_complexes.json")) and verbose: sys.stdout.write("\t--Use existing 'orthologs_complexes.json' file\n")
    if not os.path.exists(os.path.join(output_dir,  "orthologs_complexes.json")):
      complex_pdb_hits={}
      complex_pdb_chains={}
      # For each hit... #
      for hit in fimo_ext.get_hits(sort=True):
        # Get folds name  and chains#
        if verbose: sys.stdout.write("\t-- Check hit %s ...\n"%hit.get_hit())
        if external is not None:
          if not associated_pdb.has_key(hit.get_hit()): continue
          fold_hit = associated_pdb.get(hit.get_hit())
          helix_hit= associated_helix.get(hit.get_hit())
          pdb_hit  = fold_hit.split("_")[0]
          chain_hit= fold_hit.split("_")[1]
        else:
          m = re.search("(\S+)_(\d+)(\S+)",hit.get_hit())
          if m:
             pdb_hit   = m.group(1)
             helix_hit = m.group(2)
             chain_hit = m.group(3)
          else:
             pdb_hit   = hit.get_hit()[0:4] 
             helix_hit = hit.get_hit()[5:6]
             chain_hit = hit.get_hit()[6:]
          if scan_family is not None:
             if pdb_hit not in select_tf_codes: continue
          fold_hit  = pdb_hit + "_" +  chain_hit
        # Skip if next hit does not belong to selected families #
        if scan_family is not None:
           if pdb_hit not in select_tf_codes: 
              if verbose: sys.stdout.write("\t\t-- Skip %s not in the family selection\n"%(pdb_hit))
              continue
        if verbose: sys.stdout.write("\t\t-- Use TF %s in range: %d - %d with chains %s \n"%(pdb_hit,int(hit.get_start()),int(hit.get_end()),chain_hit))
        #Group the hits around an interval of DNA for each fold
        complex_tf=chain_hit+"_"+helix_hit+"_"+str(hit.get_start())+"_"+str(hit.get_end())
        complex_pdb_hits.setdefault(pdb_hit,[]).append(complex_tf)
        complex_pdb_chains.setdefault(pdb_hit,set()).add(chain_hit)

      complex_cluster_hits={}
      for pdb_hit,complex_list in complex_pdb_hits.iteritems():

        if verbose: sys.stdout.write("\t-- Check hits for fold %s ( %d hits )...\n"%(pdb_hit,len(complex_list)))
        largest_chain="".join(set([x for x in complex_pdb_chains[pdb_hit]]))

        # Obtain the adjacency matrix of hits with the same part of a TF complex
        pairing = [[False for x in range(len(complex_list))] for y in range(len(complex_list))]
        for i in range(len(complex_list)): pairing[i][i]=True
        for i in range(len(complex_list)):
           complex_hit= complex_list[i]
           start_hit  = int(complex_hit.split("_")[2])
           end_hit    = int(complex_hit.split("_")[3])
           helix_hit  = complex_hit.split("_")[1]
           chains_hit = set([x for x in complex_hit.split("_")[0]])
           binding_interval = end_hit - start_hit
           #max_complex_connect      =  max(max_position_scan, max_dimer_distance)
           range_hit  =set(range(start_hit,end_hit+1))
           for j in range(i+1,len(complex_list)):
               complex_check= complex_list[j]
               start_check  = int(complex_check.split("_")[2])
               end_check    = int(complex_check.split("_")[3])
               helix_check  = complex_check.split("_")[1]
               chains_check = set([x for x in complex_check.split("_")[0]])
               range_check  =set(range(start_check ,end_check+1 ))
               range_test   =set(range(start_check - max_complex_connect ,end_check  + max_complex_connect + 1 ))
               if len(chains_check.intersection(chains_hit))>0 and len(range_test.intersection(range_hit))>0:
                  pairing[i][j]=True
                  pairing[j][i]=True

        #Form the clusters of TFs
        cluster_complex={}
        checked_chain=[False for x in range(len(complex_list))]
        for n in range(len( largest_chain ),0,-1):
         for i in range(len(complex_list)):
           complex_hit= complex_list[i]
           start_hit  = int(complex_hit.split("_")[2])
           end_hit    = int(complex_hit.split("_")[3])
           helix_hit  = complex_hit.split("_")[1]
           chains_hit = set([x for x in complex_hit.split("_")[0]])
           if len(chains_hit) == n and not checked_chain[i]:
              cluster_complex.setdefault(complex_hit,set()).add(i)
              if verbose: sys.stdout.write("\t\t-- Starting cluster around %d - %d  with chains %s \n"%(start_hit,end_hit,chains_hit))
              checked_chain[i]=True
              for j in  range(len(complex_list)):
                  if pairing[i][j] and not checked_chain[j]:
                     cluster_complex.setdefault(complex_hit,set()).add(j)
                     if verbose: sys.stdout.write("\t\t\t-- Add hit %s \n"%(complex_list[j]))
                     checked_chain[j]=True

        #Rename the main node of each cluster
        cluster_hits={}
        for key, complex_number in cluster_complex.iteritems():
           chains_in_cluster=set()
           range_in_cluster =set()
           helix_in_cluster=set()
           for i in complex_number:
               complex_hit =  complex_list[i]
               start_hit   = int(complex_hit.split("_")[2])
               end_hit     = int(complex_hit.split("_")[3])
               helix_hit   = complex_hit.split("_")[1]
               chains_hit  = set([x for x in complex_hit.split("_")[0]])
               range_in_cluster.update(set(range(start_hit,end_hit+1)))
               chains_in_cluster.update(chains_hit)
               helix_in_cluster.update(helix_hit)
           chains  = "".join([x for x in chains_in_cluster])
           helices = "".join([x for x in helix_in_cluster])
           start   = min(range_in_cluster)
           end     = max(range_in_cluster)
           name_cluster =   chains + "_" + helices + "_" + str(start) + "_" + str(end)
           cluster_hits.setdefault(name_cluster,[])
           if verbose: sys.stdout.write("\t\t-- Rename cluster %s  as %s \n"%(key,name_cluster))
           for i in complex_number:
               complex_hit =  complex_list[i]
               cluster_hits[name_cluster].append(complex_hit)
               if verbose: sys.stdout.write("\t\t\t-- Add hit %s  to %s  \n"%(complex_hit,name_cluster))

        #Fill the complex_cluster_hits dictionary
        complex_cluster_hits.setdefault(pdb_hit,cluster_hits)

      # Generate list of ortholog complexes  
      orthologs_complexes=[]
      complex_threaded={}
      threaded=set()
      for file_thread in os.listdir(os.path.join(output_dir,"aux_files")):
          threaded.add(file_thread)
      for pdb_hit,cluster_hits in complex_cluster_hits.iteritems():
       for name_cluster, complex_list in cluster_hits.iteritems():
         tf_chain = {}
         pdb_name_cluster = str(pdb_hit)+"_"+str(name_cluster)
         for file_thread in threaded:
           mm = re.search("^%s/(\S+).(\S{4})_(\d+)(\S+).(\d+)-(\d+).txt$" % os.path.join(output_dir, "aux_files"), file_thread )
           try:
             if mm:
                 homolog_seq  = mm.group(1)
                 pdb_thread   = mm.group(2)
                 helix_thread = mm.group(3)
                 chains_thread= set([x for x in mm.group(4)])
                 thread_chain = mm.group(4)
                 start_thread = int(mm.group(5))
                 end_thread   = int(mm.group(6))
             else:
                 homolog_seq  = str(os.path.basename(file_thread).split(".")[0])
                 pdb_thread   = str(os.path.basename(file_thread).split(".")[1].split("_")[0])
                 chains_thread= set([x for x in os.path.basename(file_thread).split(".")[1].split("_")[1]])
                 helix_thread = str(os.path.basename(file_thread).split(".")[1].split("_")[1][0])
                 thread_chain = str(os.path.basename(file_thread).split(".")[1].split("_")[1][-1])
                 start_thread = int(os.path.basename(file_thread).split(".")[2].split("-")[0])
                 end_thread   = int(os.path.basename(file_thread).split(".")[2].split("-")[1])
           except:
               raise ValueError("Error handling threading file: " + str(file_thread) + "\n")
               continue
           range_thread=set(range(start_thread,end_thread))
           if pdb_hit==pdb_thread:
              if verbose: sys.stdout.write("\t--Results for PDB %s (threading %s)...\n"%(pdb_hit,os.path.basename(file_thread)))
              word     = name_cluster.split("_")
              start    = int(word[2])
              end      = int(word[3])
              helix    = set([x for x in word[1]])
              chains   = set([x for x in word[0]])
              for complex_hit in complex_list:
                  start_hit   = int(complex_hit.split("_")[2])
                  end_hit     = int(complex_hit.split("_")[3])
                  chains_hit  = set([x for x in complex_hit.split("_")[0]])
                  range_complex=set(range((start_hit-max_position_scan),(end_hit+max_position_scan)))
                  if verbose: sys.stdout.write("\t\t-- Cluster members check. Chains: %s (data %s) Binding site: %d (%d) - %d (%d) \n"%(list(chains_hit),"".join([x for x in chains_thread]),start_hit,start_thread,end_hit,end_thread))
                  if len(chains_thread.intersection(chains_hit))>0 and len(range_thread.intersection(range_complex))>0:
                    if verbose: sys.stdout.write("\t\t\t--Add threading file %s\n"%os.path.basename(file_thread))
                    tf_chain.setdefault(thread_chain,set()).add(file_thread)
              range_complex=set(range((start-max_position_scan),(end+max_position_scan)))
              if verbose: sys.stdout.write("\t\t-- Full length check. Chains: %s (data %s) Binding site: %d (%d) - %d (%d) \n"%(list(chains),"".join([x for x in chains_thread]),start,start_thread,end,end_thread))
              if  len(chains_thread.intersection(chains))>0 and len(range_thread.intersection(range_complex))>0:
                 #Assign chains with a potential threading file (i.e at least one sequence in the species can be modelled)
                 if verbose: sys.stdout.write("\t\t\t--Add threading file %s\n"%os.path.basename(file_thread))
                 tf_chain.setdefault(thread_chain,set()).add(file_thread)
         if verbose: sys.stdout.write("\t\t Complex with cluster named %s \n"%(pdb_name_cluster))
         complex_threaded.setdefault(pdb_name_cluster,tf_chain)

      for pdb_name_cluster,tf_chain in complex_threaded.iteritems():
          chain_thread_files={}
          for chain,thread_files in tf_chain.iteritems():
              chain_thread_files.setdefault(chain,[]).extend(list(thread_files))
          pdb_fold      = pdb_name_cluster.split("_")[0]
          chains        = pdb_name_cluster.split("_")[1] 
          helix_cluster = pdb_name_cluster.split("_")[2]
          start_cluster = pdb_name_cluster.split("_")[3]
          end_cluster   = pdb_name_cluster.split("_")[4]
          hits =[]
          for h in helix_cluster:
              hits.append(pdb_fold + "_" + h + chains)
          orthologs_complexes.append({'pdb_cluster':pdb_name_cluster,'pdb':pdb_fold, 'helices':helix_cluster, 'start':start_cluster,'end':end_cluster, 'hits': hits, 'chains':chain_thread_files})

      # Create json file #
      out = open(os.path.join(output_dir, "orthologs_complexes.json"), "wt")
      out.write(json.dumps(orthologs_complexes, separators=(',', ':'), indent=2))
      out.close()    
    print("total process time: %s seconds" % str(time() - start_time))

    ########################################################################
    # 3.      Create TF clusters with the output for complexes             #
    ########################################################################
    if verbose: sys.stdout.write("--TF clusters around binding sites...\n")
    # Initialize #
    orthologs_complexes = json.loads(''.join([line for line in functions.parse_file(os.path.join(output_dir, "orthologs_complexes.json"))]))
    orthologs = json.loads(''.join([line for line in functions.parse_file(os.path.join(output_dir, "orthologs.json"))]))
    if os.path.exists(os.path.join(output_dir,  "TF_clusters.json")) and verbose: sys.stdout.write("\t--Use existing 'TF_clusters.json' file\n")
    if not os.path.exists(os.path.join(output_dir,  "TF_clusters.json")):
     # Create TF_clusters data #
     group_cluster={}
     group_threads={}
     for ortholog_complex in orthologs_complexes:
         range_complex   = set(range (int(ortholog_complex['start']),int(ortholog_complex['end'])+1))
         fragment= ortholog_complex['start'] + "_" + ortholog_complex['end']
         if verbose: sys.stdout.write("\t--Check complexes in region %s - %s involving hits of %s ...\n"%(ortholog_complex['start'],ortholog_complex['end'],ortholog_complex['pdb_cluster']))
         complex_threads = set()
         complex_hits    = set()
         for chain, thread_complex_files in ortholog_complex['chains'].iteritems():
          for thread_file in thread_complex_files:
              complex_seq,complex_hit,complex_fragment,dummy_txt = thread_file.split(".")
              complex_hits.add(complex_hit)
              complex_threads.add(os.path.basename(thread_file))
         if verbose: sys.stdout.write("\t\t-- Use complex hits: %s\n"%list(complex_hits))
         if verbose: sys.stdout.write("\t\t-- Use complex threads: %s\n"%list(complex_threads))
         if verbose: sys.stdout.write("\t\t-- Use complex range:  %s\n"%list(range_complex))
         # For each cluster... #
         for i in range(len(orthologs)):
          ortholog        = orthologs[i]
          cluster_threads = set()
          for j in range(2):
              for thread_file_data in ortholog['orthologs'][j]:
                  thread_file=thread_file_data[0]
                  if thread_file is not None:
                     start_index = int(thread2index.get(os.path.basename(thread_file)))
                     cluster_threads.add(os.path.basename(thread_file))
          start_index   = start_index - 1
          range_cluster = set( range(int(start_index+int(ortholog["start"])),int(start_index+int(ortholog["end"]))+1))
          #if verbose: sys.stdout.write("\t--Check orthologs in region  (INDEX %s ) %s - %s for ID %s\n"%(str(start_index),ortholog["start"],ortholog["end"],ortholog['hit_name']))
          if len( range_complex.intersection(range_cluster) ) >0 and  len( complex_threads.intersection(cluster_threads) ) >0 :
             if verbose: sys.stdout.write("\t\t--Add cluster %d (by hit %s in [%d-%d])...\n"%(i,ortholog['hit_name'],int(start_index+int(ortholog["start"])),int(start_index+int(ortholog["end"]))))
             group_cluster.setdefault(fragment,set()).add(i+1)
             group_threads.setdefault(fragment,set()).update(complex_threads)
             group_threads.setdefault(fragment,set()).update(cluster_threads)

     if verbose: sys.stdout.write("\n\t--Merging complexes into fragments of DNA....\n")
     merge_fragments={}
     number_of_fragments=0
     for fragment in sorted(group_cluster.keys(), key=lambda x: int(x.split("_")[0])):
        binding        = fragment.split("_")
        range_fragment = set(range( int(binding[0]), int(binding[1])+1))
        if group_threads.has_key(fragment):
            thread_files_fragments = group_threads[fragment]
        skip = False
        if verbose: sys.stdout.write("\t--Check merging fragment %s with previous...\n"%fragment)
        for i in range(number_of_fragments+1):
            if merge_fragments.has_key(i):
               add_to_fragment_done = False
               for fragment_done in merge_fragments[i]:
                   binding_done = fragment_done.split("_")
                   range_done   = set(range( int(binding_done[0])-max_complex_connect, int(binding_done[1])+max_complex_connect+1))
                   if group_threads.has_key(fragment_done):
                       thread_files_fragments_done = group_threads[fragment_done]
                   if len(range_fragment.intersection(range_done))>0 and len(thread_files_fragments.intersection(thread_files_fragments_done))>0:
                       add_to_fragment_done = True
               if  add_to_fragment_done:
                   if verbose: sys.stdout.write("\t\t--Add fragment %s to segment %d\n"%(fragment,i))
                   merge_fragments[i].append( fragment )
                   skip = True
                   continue    
        if not skip:
           number_of_fragments = number_of_fragments + 1
           if verbose: sys.stdout.write("\t\t--Start segment %d\n"%(number_of_fragments))
           merge_fragments.setdefault(number_of_fragments,[]).append(fragment)

     merged_group_cluster={}
     merged_group_threads={}
     if verbose: sys.stdout.write("\n\t--Check DNA segments and group clusters of TFs....\n")
     for i in range(1,number_of_fragments):
         segment_binding=set()
         if merge_fragments.has_key(i):
            if verbose: sys.stdout.write("\t--Check segment %d\n"%i)
            for fragment in merge_fragments[i]:
                binding         = fragment.split("_")
                segment_binding.update( set(range( int(binding[0]), int(binding[1])+1)))
                if group_threads.has_key(fragment):
                   thread_files_fragments = group_threads[fragment]
                   for  thread_file in thread_files_fragments:
                      try:
                        thread_seq,thread_hit,thread_fragment,dummy_txt = thread_file.split(".")
                        binding = thread_fragment.split("-")
                        segment_binding.update( set(range( int(binding[0]), int(binding[1])+1)))
                      except:
                        if verbose: sys.stdout.write("\t\t--Skip threading file %s \n"%(thread_file))
            new_fragment=str(min(segment_binding))+"_"+str(max(segment_binding))
            for fragment in merge_fragments[i]:
                if group_cluster.has_key(fragment):
                   if verbose: sys.stdout.write("\t\t--Add clusters %s to largest fragment %d - %d\n"%(list(group_cluster[fragment]),min(segment_binding),max(segment_binding)))
                   merged_group_cluster.setdefault(new_fragment,set()).update(group_cluster[fragment])
                if group_threads.has_key(fragment):
                   merged_group_threads.setdefault(new_fragment,set()).update(group_threads[fragment])

     tf_clusters=[]
     if verbose: sys.stdout.write("\n\t--Check DNA segments and group threaded TFs....\n")
     for fragment in sorted(merged_group_cluster.iterkeys(), key=lambda x: int(x.split("_")[0])):
      binding= fragment.split("_")
      group_threads_by_pdb = {}
      pdb_chains={}
      tf_group_threads_by_pdb = {}
      tf_pdb_chains={}
      list_of_hits=set()
      if verbose: sys.stdout.write("\t--Check cluster around region %s - %s \n"%(binding[0],binding[1]))
      if merged_group_threads.has_key(fragment):
        for thread_file in merged_group_threads[fragment]:
            data             = thread_file.split(".")
            sequence         = data[0]
            fimo_hit         = data[1]
            pdb,helix_chain  = fimo_hit.split("_")
            helix            = helix_chain[0]
            chain            = helix_chain[1]
            if verbose: sys.stdout.write("\t\t--Add threaded file %s ...\n"%(thread_file))
            group_threads_by_pdb.setdefault((pdb,chain),set()).add(thread_file)
            pdb_chains.setdefault(pdb,set()).add(chain)
        for pdb,chain in group_threads_by_pdb.iterkeys():
            tf_group_threads_by_pdb.setdefault(pdb+"_"+chain,list(group_threads_by_pdb[(pdb,chain)]))
        for pdb,chains in  pdb_chains.iteritems():
            tf_pdb_chains.setdefault(pdb,list(chains))   
        for nclust in merged_group_cluster[fragment]:
            cluster = orthologs[nclust-1]
            list_of_hits.add(cluster['hit_name'])
        tf_clusters.append({'start':binding[0],'end':binding[1],'clusters': list(merged_group_cluster[fragment]), 'node_hits':list(list_of_hits), 'threads': tf_group_threads_by_pdb, 'pdb_chains': tf_pdb_chains})

     # Create json file #
     out = open(os.path.join(output_dir, "TF_clusters.json"), "wt")
     out.write(json.dumps(tf_clusters, separators=(',', ':'), indent=2))
     out.close()    

    #Create output file
    if verbose: sys.stdout.write("--Create output of merged TF clusters around binding sites...\n")
    output_file = os.path.join(output_dir, "TF_clusters.txt")
    out = open(output_file, "wt")
    tf_clusters = json.loads(''.join([line for line in functions.parse_file(os.path.join(output_dir, "TF_clusters.json"))]))
    for tf_cluster in tf_clusters:
     out.write(">%s-%s | %s\n"%(tf_cluster['start'],tf_cluster['end'],tf_cluster['clusters']))
     tf_group_threads_by_pdb = tf_cluster['threads']
     tf_pdb_chains=tf_cluster['pdb_chains']
     for pdb,chains in tf_pdb_chains.iteritems():
        out.write("\tPDB %s:\n"%pdb)
        for chain in chains:
            out.write("\tChain  %s:\n"%chain)
            pdb_chain=pdb+"_"+chain
            if tf_group_threads_by_pdb.has_key(pdb_chain):
               for thread_file in tf_group_threads_by_pdb[pdb_chain]:
                   out.write("\t\t%s\n"%thread_file)
     out.write("//\n")
    out.close()

    
    ########################################################################
    # 4.    DONE                                                           #
    ########################################################################
    if verbose:sys.stdout.write("Done\n")
        

