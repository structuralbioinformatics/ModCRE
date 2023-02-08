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
# Functions   #
#-------------#

def fileExist(file):
    '''
    Check existing files
    '''
    if file is not None:
        return os.path.exists(file) and os.path.isfile(file)
    else:
        return False

def make_subdirs(main, subdirs):
    '''
    This function makes all subdirs listed in "subdirs".
    '''

    for subdir in subdirs:
        if not os.path.exists(os.path.join(main, subdir)):
            os.makedirs(os.path.join(main, subdir))

def remove_files(files):
    '''
    This function removes all files listed in "files".
    '''

    for each_file in files:
        if os.path.exists(each_file):
            os.remove(each_file)


def parse_best_orthologs(input_file,rank):

    orthologs_list = []

    # Get parameters from configuration
    max_orthologs = int(config.get("Parameters","max_orthologs"))
    #store families per pdb_chain
    families={}
    for line in functions.parse_file(os.path.join(pdb_dir, "families.txt")):
            if line.startswith("#"): continue
            pdb_chain, family = line.split(";")
            families[pdb_chain] = family


    inp = open(input_file,"r")
    for line in inp:
        if line.startswith(">"): 
            #initialize
            ortholog_dict={}
            fragment_dict={}
            monomers     =[]
            dimers       =[]
            #store data
            data=line.strip().split("|")
            interval = data[1]
            start    = int(interval.split("-")[0])
            fragment_dict.setdefault("start",start)
            end      = int(interval.split("-")[1])
            fragment_dict.setdefault("end",end)
            pval     = float(data[2])
            fragment_dict.setdefault("pval",pval)
            hit_name = data[3]
            fragment_dict.setdefault("hit_name",hit_name)
            fragment_dict.setdefault("proteins",set())
            fragment_dict.setdefault("monomer",set())
            fragment_dict.setdefault("dimer",set())
            fragment_dict.setdefault("families",set())
            fragment_dict.setdefault("pdb_chains",set())
            continue
        if line.startswith("//"):
            #print "NEW ORTHOLOG LIST"
            if len(monomers)>0:
               #print monomers
               monomer_sorted = [monomer[1] for monomer in sorted(monomers,key=lambda x: x[0])]
               if not rank: max_orthologs=len(monomer_sorted)
               for thread in monomer_sorted[:min(max_orthologs,len(monomer_sorted))]:
                  uid,gene,thread_file,score,d_score = thread.split(";")
                  fragment_dict["proteins"].add((uid,gene))
                  fragment_dict["monomer"].add((thread_file,score,d_score))
                  pdb_Hchain = thread_file.split(".")[1]
                  pdb_chain  = pdb_Hchain[0:4] + "_" + pdb_Hchain[6:7]
                  fragment_dict["pdb_chains"].add(pdb_chain)
                  if families.has_key(pdb_chain): fragment_dict["families"].add(families[pdb_chain])
            if len(dimers)>0:
               #print dimers
               dimer_sorted   = [dimer[1] for dimer in sorted(dimers,key=lambda x: x[0])]
               if not rank: max_orthologs=len(dimer_sorted)
               for thread in dimer_sorted[:min(max_orthologs,len(dimer_sorted))]:
                  #read thread files 
                  monomer_A,monomer_B = thread
                  uid_A,gene_A,thread_A,score_A,d_score_A = monomer_A.split(";")
                  uid_B,gene_B,thread_B,score_B,d_score_B = monomer_B.split(";")
                  #define fragment_dict for uidA
                  fragment_dict["proteins"].add((uid_A,gene_A))
                  fragment_dict["monomer"].add((thread_A,score_A,d_score_A))
                  pdb_Hchain = thread_A.split(".")[1]
                  pdb_chain  = pdb_Hchain[0:4] + "_" + pdb_Hchain[6:7]
                  fragment_dict["pdb_chains"].add(pdb_chain)
                  if families.has_key(pdb_chain): fragment_dict["families"].add(families[pdb_chain])
                  #define fragment_dict for uidB
                  fragment_dict["proteins"].add((uid_B,gene_B))
                  fragment_dict["monomer"].add((thread_B,score_B,d_score_B))
                  pdb_Hchain = thread_B.split(".")[1]
                  pdb_chain  = pdb_Hchain[0:4] + "_" + pdb_Hchain[6:7]
                  fragment_dict["pdb_chains"].add(pdb_chain)
                  if families.has_key(pdb_chain): fragment_dict["families"].add(families[pdb_chain])
                  #define fragment_dict for dimer
                  dimer=((thread_A,score_A,d_score_A),(thread_B,score_B,d_score_B))
                  idimer=((thread_B,score_B,d_score_B),(thread_A,score_A,d_score_A))
                  if dimer in fragment_dict["dimer"]: continue
                  if idimer in fragment_dict["dimer"]: continue
                  fragment_dict["dimer"].add(dimer)
            if len(dimers)>0 or len(monomers)>0:
                  fragment_dict["dimer"]      = [x for x in fragment_dict["dimer"] ]
                  fragment_dict["monomer"]    = [x for x in fragment_dict["monomer"] ]  
                  fragment_dict["proteins"]   = [x for x in fragment_dict["proteins"] ]
                  fragment_dict["pdb_chains"] = [x for x in fragment_dict["pdb_chains"] ]
                  fragment_dict["families"]   = [x for x in fragment_dict["families"] ]
                  ortholog_dict = fragment_dict
            if ortholog_dict != {}: orthologs_list.append(ortholog_dict)
            continue
        if not line.startswith("//") and not line.startswith(">"):
            monomer_A,monomer_B = line.strip().split()
            if monomer_A == "None" or len(monomer_A.split(";"))<5: monomer_A=None
            if monomer_B == "None" or len(monomer_B.split(";"))<5: monomer_B=None
            if monomer_A is not None:
                uid_A,gene_A,thread_A,score_A,d_score_A = monomer_A.split(";")
                if score_A != 0 and d_score_A != 0 : monomers.append((float(score_A),monomer_A))
            if monomer_B is not None:
                uid_B,gene_B,thread_B,score_B,d_score_B = monomer_B.split(";")
                if score_B != 0 and d_score_B != 0 : monomers.append((float(score_B),monomer_B))
            if monomer_B is not None and monomer_A is not None:
                if score_B != 0 and d_score_B != 0 and score_A != 0 and d_score_A != 0 :
                   dimers.append((min(float(score_A),float(score_B)),(monomer_A,monomer_B)))

    inp.close()

    return orthologs_list



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
    parser.add_option("--info",default=None,action="store", type="string", dest="info",help="Information LOG file of SCANS that have failed and have been completed")
    parser.add_option("--reuse",default=False, action="store_true", dest="reuse", help="Reuse the information files. If the flag is used then scans that had failed will remain as FAILED, otherwise it tries to redo them (default=False)")
    parser.add_option("--rank",default=False, action="store_true", dest="rank", help="Rank by binding scoring-energy. If the flag is not used then all orthologs get a  null score and the programs runs faster, otherwise the sacn is slow but the selection of orthologs is improved (default=False)")
    parser.add_option("--index", action="store", type="int",dest="index",  default=1, help="Index coordinates of the DNA sequence to locate the starting position of the sequence", metavar="{integer}")


    group = optparse.OptionGroup(parser, "Statistical potentials", "Select your statistical potentials of choice. In \"--auto\" mode, the program uses S3DCdd family potentials derived from both PDB and PBM data and approached by Taylor.\"-a\" option overrides options \"-f\", \"-p\", \"-m\", \"-r\" and \"-t\".")
    group.add_option("-a", "--auto", default=False, action="store_true", dest="auto_mode", help="Automate the selection of statistical potentials (default = False)")
    group.add_option("-f", "--family", default=False, action="store_true", dest="family_potentials", help="Use family potentials (default = False)")
    group.add_option("-p", default=False, action="store_true", dest="pbm_potentials", help="Use potentials derived from both PBM + PDB data (default = False)")
    group.add_option("-m", "--pmf", default=False, action="store_true", dest="pmf_potentials", help="Use of raw mean-force potentials with no Z-scoring (default = False)")
    group.add_option("--pot", default=None, action="store", type="string", dest="split_potential", help="Split-potential to be used (3d, 3dc, s3dc, s3dc_dd, s3dc_di, pair; default from configuration file)", metavar="{string}")
    group.add_option("--taylor", default=False, action="store_true", dest="taylor_approach", help="Approach PMF by Taylor (default = False)")
    group.add_option("-t", action="store", type="float", default=None, dest="score_threshold", help="Threshold on the scaled score to consider positive k-mers (default = 0.95)", metavar="{float}")
    #group.add_option("-k","--known", default=False, action="store_true", dest="known", help="The name is of a known PDB file, with format 'code_chain' (default = False)")
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

start_time = time()

# Add "." to sys.path #
src_path =  os.path.abspath(os.path.dirname(__file__))
sys.path.append(src_path)

# Imports jbonet's module #
from SBI.structure import PDB
from SBI.data      import aminoacids1to3

# Imports my functions #
import functions, fimo, tomtom, interface, x3dna, contacts, triads,  model_protein, nr, threader,dssp
import pwm_pbm as PWM
import homologs as HOMO

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(src_path, "config.ini")
config.read(config_file)

# Define which python to be used #
python = os.path.join(config.get("Paths", "python_path"), "python")

# Arguments & Options #
options = parse_options()
dummy_dir = os.path.abspath(options.dummy_dir)
if not os.path.exists(dummy_dir): os.makedirs(dummy_dir)
input_file = os.path.abspath(options.input_file)
label = options.label
output_dir = os.path.abspath(options.output_dir)
pbm_dir = os.path.abspath(options.pbm_dir)
pdb_dir = os.path.abspath(options.pdb_dir)
specie = options.specie
verbose = options.verbose
scan_family = options.scan_family
external = options.external
cluster_complexes =options.cluster_complexes
reuse = options.reuse
rank  = options.rank
info_file = options.info
if info_file is None: info_file=os.path.join(output_dir,"scanner.log")
if not fileExist(info_file):
        info=open(info_file,"w")
        info.write("#Scanning  results\n")
        info.close()

start_index = options.index - 1

#Always try as known and let the algorithm choose otherwise
#known = options.knwon
known = True
auto_mode = options.auto_mode

# Create output directory #
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
# Create output aux_files subdirs #
make_subdirs(output_dir, subdirs=["aux_files"])

########################
# Start run            #
########################

if verbose: sys.stdout.write("\nSCAN of DNA sequence with a database of PWMs\n  Structural Bioinformatics Lab 2017\n")

########################
# 1. Parse dimers file #
########################
if verbose: sys.stdout.write("\n--Parse dimers file...\n")
# Initialize #
dimers = {}
# For each line... #
for line in functions.parse_file(os.path.join(pdb_dir, "dimers.txt")):
    chain_A, chain_B, contacts_dimer, overlap = line.split(";")
    dimers[chain_A] = chain_B
    dimers[chain_B] = chain_A

#########################
# 2. Parse DNA sequence #
#########################
if verbose: sys.stdout.write("--Parse DNA fasta file...\n")
# For each header, sequence... #
for fasta_sequence in functions.parse_fasta_file(input_file):
    dna_header = fasta_sequence[0]
    dna_sequence = fasta_sequence[1]

########################
# 3. Scan DNA sequence #
########################
if verbose: sys.stdout.write("--Scan DNA sequence to database of PWMs...\n")
# Execute fimo #
meme_path = os.path.join(src_path, config.get("Paths", "meme_path"))
if options.database is None and external is None:
    # Use the default databases #
    database_file = os.path.join(pbm_dir, "pwms", "database.txt")
    database_extended_file = os.path.join(pbm_dir, "pwms", "database_extended.txt")
    homologs_dir  = os.path.join(pbm_dir,"homologs") 
    if not fileExist(database_extended_file):
        if verbose: sys.stdout.write("\t--DBD database extended not found, we generate it with pwms\n")  
        for pwm_file in os.listdir(os.path.join(pbm_dir, "pwms")):
            if pwm_file.endswith("meme"): 
               m = re.search("(\S+)_(\d+)(\S+).meme", pwm_file)
               if m:
                  os.system("cat %s >> %s" % (os.path.join(pbm_dir, "pwms",pwm_file), database_extended_file))
    else:
        if verbose: sys.stdout.write("\t--Database extended for all complexes is already available as 'database_extended.txt'...\n")
    if not fileExist(database_file):
        if verbose: sys.stdout.write("\t--DBD database not found, we generate it with pwms\n")  
        # Check TF families
        tf_codes=set()
        if os.path.exists(os.path.join(pdb_dir, "families.txt")):
          for line in functions.parse_file(os.path.join(pdb_dir, "families.txt")):
              if line.startswith("#"):continue
              pdb_name, family = line.split(";")
              if family != "Unknown" and family != "Undefined": tf_codes.add(pdb_name[0:4].lower())
        for pwm_file in os.listdir(os.path.join(pbm_dir, "pwms")):  
            if pwm_file.endswith("meme"): 
               m = re.search("(\S+)_(\d+)(\S+).meme", pwm_file)
               if m: 
                 chains =  m.group(3)
                 pdb_code = m.group(1).lower()
                 if len(chains) == 1 and pdb_code in tf_codes: os.system("cat %s >> %s" % (os.path.join(pbm_dir, "pwms",pwm_file), database_file))

    else:
        if verbose: sys.stdout.write("\t--Database for DBDs is already available as 'database.txt'...\n")
elif external is None:
    # Use the database choosen by the user #
    database_file = os.path.abspath(options.database)
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
    print("SCAN ASSOCIATION OPEN  "+association_file)
    for line in assoc:
        data=line.split()
        m = re.search("(\S+)_(\d+)(\S+)",data[1])
        if m:
            pdb_name  = m.group(1)
            helix_hit = m.group(2)
            pdb_chain = m.group(3)
            fold_hit  = pdb_name+"_"+pdb_chain
        else:
            fold_hit  = data[1]
        associated_pdb.setdefault(data[0],fold_hit)
        associated_pdb_inv.setdefault(fold_hit,set()).add(data[0])
        if fold2pwm.has_key(data[1]):
           associated_helix.setdefault(data[0],str(sorted([int(pdb_helix_chain[5:6]) for pdb_helix_chain in fold2pwm.get(data[1])])[0]))
        else:
           associated_helix.setdefault(data[0],"1")
        print("SCAN ASSOCIATION PDB "+data[0]+" <==> "+data[1]+ " FOLD_HIT " + fold_hit + " HELIX "+associated_helix.get(data[0]))
    assoc.close()
    database_extended_file = database_file
try:
        fimo_obj = fimo.get_fimo_obj(database_file, input_file,options.fimo_pvalue_threshold,options.max_stored_matches,os.path.abspath(options.dummy_dir))
        fimo_ext = fimo.get_fimo_obj(database_extended_file, input_file,options.fimo_pvalue_threshold,options.max_stored_matches,os.path.abspath(options.dummy_dir))
except:
        if verbose == True: sys.stdout.write("\n")
        sys.stdout.write("\"fimo\" execution failed!")
        if verbose == True: sys.stdout.write("\n\n")
        info=open(info_file,"a")
        info.write("%s\tFAIL\n"%(os.path.basename(input_file)))
        info.close()
        exit(0)
print("total process time: %s seconds" % str(time() - start_time))

#########################################
# 4. Get FIMO OBJECT and built clusters #
#########################################
if verbose: sys.stdout.write("--Use PWMs of folds...\n")
# Initialize #
members = []   
clusters = [] 
max_dimer_distance=int(config.get("Parameters", "max_dimer_distance")) 
max_position_scan =int(config.get("Parameters", "max_position_scan"))
#Select special family
if scan_family is not None:
   select_families = set(scan_family.split(","))
   select_tf_codes=set()
   if os.path.exists(os.path.join(pdb_dir, "families.txt")):
      for line in functions.parse_file(os.path.join(pbm_dir, "families.txt")):
          if line.startswith("#"):continue
          pdb_name, family = line.split(";")
          if family in select_families: select_tf_codes.add(pdb_name[0:4].lower())

#Show total number of found  
if verbose: sys.stdout.write("\t--Total number of hits: %d\n"%(len(fimo_obj.get_hits())))
members_storage=os.path.join(options.output_dir,"members.pickle")
cluster_storage=os.path.join(options.output_dir,"clusters.pickle")
if not fileExist(members_storage) or not fileExist(cluster_storage):
 # For each hit... #
 for hit in fimo_obj.get_hits(sort=True):
    # Skip if already done #
    if (hit.get_hit(), hit.get_start(), hit.get_end()) in members: continue
    print("SCAN FIMO HIT "+hit.get_hit())
    # Initialize #
    # Get folds file #
    fold = []
    if external is not None:
      if not associated_pdb.has_key(hit.get_hit()): continue
      fold_hit   = associated_pdb.get(hit.get_hit())
      helix_hit  = associated_helix.get(hit.get_hit())
      pdb_hit    = fold_hit.split("_")[0]
      chain_hit  = fold_hit.split("_")[1]
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
    folds_file = os.path.join(pdb_dir, "folds", fold_hit + ".txt")
    if not os.path.exists(folds_file): continue
    fold.append(fold_hit)
    if os.path.exists(folds_file):
      for line in functions.parse_file(folds_file):
        pdb_chain, tm_score = line.split(";")
        pdb_chain_file=os.path.join(pdb_dir, "split",pdb_chain+".pdb")
        if not os.path.exists(pdb_chain_file):continue
        fold.append(pdb_chain)
    print("SCAN CLUSTER ADD "+pdb_hit+" FOLD HIT "+ fold_hit+ " FIMO HIT "+hit.get_hit()+ " START "+str( hit.get_start()))
    clusters.append({'pdb' : pdb_hit, 'helix' : helix_hit, 'chain' :  chain_hit, 'fold' : fold_hit, 'hit_name' : hit.get_hit(), 'from' : hit.get_start(), 'to' : hit.get_end(), 'members' : [[hit, None]]})
    range_cluster=set(range( hit.get_start(),  hit.get_end()+1))
    if verbose:sys.stdout.write("\t--Use hit %s\n"%hit.get_hit())
    # For each hit... #
    for next_hit in fimo_obj.get_hits(sort=True):
        # Skip if already done #
        if (next_hit.get_hit(), next_hit.get_start(), next_hit.get_end()) in members: continue
        # Initialize #
        if external is not None:
           if not associated_pdb.has_key(next_hit.get_hit()): continue
           next_fold_hit   = associated_pdb.get(next_hit.get_hit())
           next_helix_hit  = associated_helix.get(next_hit.get_hit())
           next_hit_file   = os.path.join(database_dir,next_hit.get_hit() + ".meme")
           next_pdb_hit    = next_fold_hit.split("_")[0]
           next_chain_hit  = next_fold_hit.split("_")[1]
        else:
           mm = re.search("(\S+)_(\d+)(\S+)",next_hit.get_hit())
           if mm:
             next_pdb_hit   = mm.group(1)
             next_helix_hit = mm.group(2)
             next_chain_hit = mm.group(3)
           else:
             next_pdb_hit   = next_hit.get_hit()[0:4] 
             next_helix_hit = next_hit.get_hit()[5:6]
             next_chain_hit = next_hit.get_hit()[6:]
           next_fold_hit =  next_pdb_hit + "_" + next_chain_hit 
           next_hit_file = os.path.join(pbm_dir, "pwms", next_hit.get_hit() + ".meme")
        # Skip if next hit does not belong to hit fold #
        if next_fold_hit not in fold: continue
        # Skip if next hit does not belong to selected families #
        if scan_family is not None:
           if next_pdb_hit not in select_tf_codes: continue
        if verbose:sys.stdout.write("\t\t--Check to include hit %s\n"%next_hit.get_hit())
        range_next_hit = set(range (next_hit.get_start()-max_position_scan,next_hit.get_end()+max_position_scan+1))
        if len(range_next_hit.intersection(range_cluster)) <= 0 :
           if verbose:sys.stdout.write("\t\t\t--Skip %s unconnected regions %d-%d vs %d-%d\n"%(next_hit.get_hit(),min(list(range_next_hit)),max(list(range_next_hit)),min(list(range_cluster)),max(list(range_cluster))))
           continue
        # For each pair of members... #
        for member_pair in clusters[-1]["members"]:
            # For each member... #
            for member in member_pair:
                # Skip if undefined member #
                if member is None: continue
                # Skip if member is next hit #
                if member == next_hit: continue
                # Skip if already done
                if (next_hit.get_hit(), next_hit.get_start(), next_hit.get_end()) in members: continue
                if verbose:sys.stdout.write("\t\t\t--Compare with %s \n"%member.get_hit())
                # Execute tomtom #
                meme_path = os.path.join(src_path, config.get("Paths", "meme_path"))
                if external is not None:
                  member_file = os.path.join(database_dir, member.get_hit() + ".meme")
                else:
                  member_file = os.path.join(pbm_dir, "pwms", member.get_hit() + ".meme")
                try:
                    tomtom_obj= tomtom.get_tomtom_obj(member_file, next_hit_file,os.path.abspath(options.dummy_dir))
                   # process = subprocess.check_output([meme_path + "tomtom", "-text", member_file, next_hit_file], stderr=subprocess.PIPE)
                except:
                    if verbose == True: sys.stdout.write("\n")
                    sys.stdout.write('"tomtom" execution failed! Skip file')
                    if verbose == True: sys.stdout.write("\n")
                    #exit(0)
                    continue
                for tomtom_hit in tomtom_obj.get_hits(sort=True):
                    # Skip if tomtom overlap does not equal fimo overlap #
                    if hit.get_strand() == next_hit.get_strand() and tomtom_hit.get_strand() != "+": continue
                    if hit.get_strand() != next_hit.get_strand() and tomtom_hit.get_strand() != "-": continue
                    if abs( hit.get_start()  - (next_hit.get_start() + (int(tomtom_hit.get_offset()) * int(hit.get_strand() + "1"))) )> max_position_scan: continue
                    if verbose:sys.stdout.write("\t\t\t\t--Checking strand-orientation %s vs %s (TOMTOM %s) and start position %d vs %d \n"%(hit.get_strand(),next_hit.get_strand(),tomtom_hit.get_strand(),hit.get_start(),next_hit.get_start() + (tomtom_hit.get_offset() * int(hit.get_strand() + "1"))))
                    # Expand cluster if necessary... #
                    if next_hit.get_start() < clusters[-1]['from']: clusters[-1]['from'] = next_hit.get_start()
                    if next_hit.get_end()   > clusters[-1]['to']  : clusters[-1]['to']   = next_hit.get_end()
                    if verbose:sys.stdout.write("\t\t\t\t--Accepted to modify cluster limits %d - %d\n"%( clusters[-1]['from'],clusters[-1]['to']))
                    # Add to cluster #
                    if member == member_pair[0]:
                        if verbose:sys.stdout.write("\t\t\t\t--Include hit %s\n"%next_hit.get_hit())
                        clusters[-1]["members"].append([next_hit, None])
                    if member == member_pair[1]:
                        if verbose:sys.stdout.write("\t\t\t\t--Include hit %s\n"%next_hit.get_hit())
                        clusters[-1]["members"].append([None, next_hit])
                    # Add to members #
                    members.append((next_hit.get_hit(), next_hit.get_start(), next_hit.get_end()))
        # If next hit dimerizes... #
        if dimers.has_key(next_fold_hit):
            # Get interfaces #
            interface_A = interface.Interface(os.path.join(pdb_dir, "interfaces", next_fold_hit + ".txt"))
            interface_B = interface.Interface(os.path.join(pdb_dir, "interfaces", dimers[next_fold_hit] + ".txt"))
            # Get from and to #
            from_B = next_hit.get_start() - (interface_A.get_start() - interface_B.get_start())
            to_B = next_hit.get_end() - (interface_A.get_end() - interface_B.get_end())
            sense = next_hit.get_strand()
            sequence = dna_sequence[from_B - 1:to_B]
            # For each hit... #
            if verbose:sys.stdout.write("\t\t--Dimer fold %s\n"%next_fold_hit)
            for next_next_hit in fimo_obj.get_hits(sort=True):
                # Skip if already done #
                if (next_next_hit.get_hit(), next_next_hit.get_start(), next_next_hit.get_end()) in members: continue
                # If the 2nd monomer is among hits... #
                if external is not None:
                  if not associated_pdb.has_key(next_next_hit.get_hit()): continue
                  next_next_fold_hit  = associated_pdb.get(next_next_hit.get_hit())
                  next_next_helix_hit = associated_helix.get(next_next_hit.get_hit())
                  next_next_hit_file  = os.path.join(database_dir,next_next_hit.get_hit() + ".meme")
                  next_next_pdb_hit   = next_next_fold_hit.split("_")[0]
                  next_next_chain_hit = next_next_fold_hit.split("_")[1]
                else:
                  mmm = re.search("(\S+)_(\d+)(\S+)",next_next_hit.get_hit())
                  if mmm:
                    next_next_pdb_hit   = mmm.group(1)
                    next_next_helix_hit = mmm.group(2)
                    next_next_chain_hit = mmm.group(3)
                  else:
                    next_next_pdb_hit   = next_next_hit.get_hit()[0:4] 
                    next_next_helix_hit = next_next_hit.get_hit()[5:6]
                    next_next_chain_hit = next_next_hit.get_hit()[6:]
                  next_next_fold_hit =  next_next_pdb_hit + "_" + next_next_chain_hit
                binding_interval = next_next_hit.get_end()- next_next_hit.get_start()
                if next_next_fold_hit == dimers[next_fold_hit] and abs(next_next_hit.get_start() - from_B) < max_dimer_distance  and abs(next_next_hit.get_end() - to_B) < max_dimer_distance  and next_next_hit.get_strand() == sense and next_next_helix_hit==next_helix_hit:
                    # Expand cluster if necessary... #
                    if next_next_hit.get_start() < clusters[-1]['from']: clusters[-1]['from'] = next_next_hit.get_start()
                    if next_next_hit.get_end() > clusters[-1]['to']: clusters[-1]['to'] = next_next_hit.get_end()
                    # Add to cluster #
                    if clusters[-1]["members"][-1][0] is None:
                       if verbose:sys.stdout.write("\t\t\t--Include hit %s on dimer side  [ %d - %d ]\n"%(next_next_hit.get_hit(),next_next_hit.get_start(),next_next_hit.get_end()))
                       clusters[-1]["members"][-1][0] = next_next_hit
                    if clusters[-1]["members"][-1][1] is None:
                       if verbose:sys.stdout.write("\t\t\t--Include hit %s on dimer side  [ %d - %d ]\n"%(next_next_hit.get_hit(),next_next_hit.get_start(),next_next_hit.get_end()))
                       clusters[-1]["members"][-1][1] = next_next_hit
                    # Add to members #
                    members.append((next_next_hit.get_hit(), next_next_hit.get_start(), next_next_hit.get_end()))
                    continue
            # If not done create the dimer hit artificially... #
            dimer_hits=set()
            if external is not None: 
                dimer_fold_hit = dimers[next_fold_hit][0:4]+"_"+dimers[next_fold_hit][-1]
                if associated_pdb_inv.has_key(dimer_fold_hit):
                    for dimer_pwm in associated_pdb_inv[dimer_fold_hit]:
                        dimer_hits.add(dimer_pwm)
            else: dimer_hits.add(dimers[next_fold_hit][0:4]+"_"+next_helix_hit+dimers[next_fold_hit][-1])
            for dimer_hit in dimer_hits:
              if (dimer_hit, from_B, to_B) not in members:
                if clusters[-1]["members"][-1][0] is None:
                   if verbose:sys.stdout.write("\t\t\t--Include hit %s on dimer side [ %d - %d ]\n"%(dimer_hit,int(from_B),int(to_B)))
                   clusters[-1]["members"][-1][0] = fimo.FimoHit(dimer_hit,int(from_B),int(to_B),sense,None,None,sequence)
                if clusters[-1]["members"][-1][1] is None:
                   if verbose:sys.stdout.write("\t\t\t--Include hit %s on dimer side [ %d - %d ]\n"%(dimer_hit,int(from_B),int(to_B)))
                   clusters[-1]["members"][-1][1] = fimo.FimoHit(dimer_hit,int(from_B),int(to_B),sense,None,None,sequence)
                members.append((dimer_hit, from_B, to_B))
    # Add to members #
    members.append((hit.get_hit(), hit.get_start(), hit.get_end()))
 # Save members and clusters
 dump_members=open(members_storage,"wb")
 if verbose:sys.stdout.write("\t-- Save MEMBERS %s\n"%members_storage)
 cPickle.dump(members,dump_members)
 dump_members.close()
 dump_cluster=open(cluster_storage,"wb")
 if verbose:sys.stdout.write("\t-- Save CLUSTERS %s\n"%cluster_storage)
 cPickle.dump(clusters,dump_cluster)
 dump_cluster.close()
else:
 # Load  members and clusters
 if verbose:sys.stdout.write("\t-- Load MEMBERS %s\n"%members_storage)
 dump_members=open(members_storage,"rb")
 members=cPickle.load(dump_members)
 dump_members.close()
 if verbose:sys.stdout.write("\t-- Load CLUSTER %s\n"%cluster_storage)
 dump_cluster=open(cluster_storage,"rb")
 clusters=cPickle.load(dump_cluster)
 dump_cluster.close()


print("total process time: %s seconds" % str(time() - start_time))

#####################################################
# 5. Thread orthologs and built homologs dictionary #
#####################################################
if verbose: sys.stdout.write("--Thread orthologous transcription factors...\n")
# Initialize #
homologs = {}
orthologs = []
emboss_path = os.path.join(src_path, config.get("Paths", "emboss_path"))
# For each cluster... #
homologs_storage=os.path.join(options.output_dir,"homologs.pickle")
orthologs_storage=os.path.join(options.output_dir,"orthologs.pickle")
if not fileExist(homologs_storage) or not fileExist(orthologs_storage):
 for cluster in clusters:
    # Initialize #
    orthologs.append({'pdb_chain' : cluster['fold'], 'hit_name' : cluster['hit_name'],'start' : cluster['from'], 'end' : cluster['to'], 'orthologs' : [[], [], []]})
    # For each pair of members... #
    if verbose: sys.stdout.write("\t--Check hit %s\n"%cluster['hit_name'])
    for members in cluster['members']:
        # For each member... #
        for i in range(2):
            # Skip if monomer is None #
            if members[i] is None: continue
            if external is not None:
               if not associated_pdb.has_key(members[i].get_hit()): continue
               fold_hit = associated_pdb.get(members[i].get_hit())
               helix_hit= associated_helix.get(members[i].get_hit())
               pdb_name = fold_hit.split("_")[0]
               pdb_chain= fold_hit.split("_")[1]
            else:
               # Get chain PDB obj #
               m = re.search("(\S+)_(\d+)(\S+)",members[i].get_hit())
               if m:
                  pdb_name  = m.group(1)
                  helix_hit = m.group(2)
                  pdb_chain = m.group(3)
               else:
                  pdb_name  = members[i].get_hit()[:4]
                  helix_hit = members[i].get_hit()[5:6]
                  pdb_chain = members[i].get_hit()[6:]
               fold_hit  = pdb_name+"_"+pdb_chain
            if verbose: sys.stdout.write("\t\t--Using fold %s for threading...\n"%fold_hit)
            if not os.path.exists(os.path.join(pdb_dir, "split",  fold_hit + ".pdb")):continue
            pdb_obj = PDB(os.path.join(pdb_dir, "split",  fold_hit + ".pdb"))
            # Get sequence/PDB correlation #
            sequence_to_crystal, crystal_to_sequence = model_protein.get_sequence_to_crystal_correlations(pdb_obj, pdb_chain)
            # Get interface object #
            if not os.path.exists(os.path.join(pdb_dir, "interfaces", fold_hit + ".txt")): continue
            interface_obj = interface.Interface(os.path.join(pdb_dir, "interfaces", fold_hit + ".txt"))
            # Get helix 
            dna_helix=helix_hit
            # Get binding site #
            binding_site = range(interface_obj.get_start(), interface_obj.get_end() + 1)
            # Initialize #
            if not os.path.exists( os.path.join(homologs_dir, fold_hit + ".txt")): continue
            homologs_file = os.path.join(homologs_dir, fold_hit + ".txt")
            try:
               homologs_obj  = HOMO.Homologs(homologs_file)
               if verbose: sys.stdout.write("\t\t--homologs object created properly for %s\n"%fold_hit)
            except:
                if verbose: sys.stdout.write("\t\t--skip homologs of %s\n"%fold_hit)
                continue
            # For each ortholog... #
            ortholog_hits=homologs_obj.get_orthologs(specie)
            # Get all hits according to the lowest restruction of Rost curve (parameters 0,None)
            for ortholog in ortholog_hits.get_hits((0,None)):
                # Add homolog #
                mm = re.search("^(tr|sp|TR|SP)\|(\S+)\|\S+\_(\S+)", ortholog.sequenceID)
                if mm:
                   ortholog_name = mm.group(2)
                else:
                   ortholog_name=ortholog.sequenceID.strip().split()[0].split("|")[0]
                # If threaded file does not exist... #
                if external is not None:
                    if associated_pdb.has_key(members[i].get_hit()): 
                        fold_hit = associated_pdb.get(members[i].get_hit())
                        threading_file = os.path.join(output_dir, "aux_files", ortholog_name + "." + fold_hit + "." + str(start_index+members[i].get_start()) + "-" + str(start_index+members[i].get_end()) + ".txt")
                    else:
                        threading_file = os.path.join(output_dir, "aux_files", ortholog_name + "." + members[i].get_hit() + "." + str(start_index+members[i].get_start()) + "-" + str(start_index+members[i].get_end()) + ".txt")
                else:
                   threading_file = os.path.join(output_dir, "aux_files", ortholog_name + "." + members[i].get_hit() + "." + str(start_index+members[i].get_start()) + "-" + str(start_index+members[i].get_end()) + ".txt")
                if os.path.exists(threading_file):
                    homologs.setdefault(ortholog_name, ortholog)
                    threaded_ortholog = threader.Threaded(threading_file)
                    # Add to ortholog file #
                    if verbose: sys.stdout.write("\t\t\t--Add threading file %s ...\n"%(os.path.basename(threading_file)))
                    orthologs[-1]['orthologs'][i].append([threading_file, 0, 0])
                else:
                    threaded_ortholog = threader.Threaded(None)
                    # Initialize #
                    try:
                     template_alignment = ortholog.sequences[0].sequence.upper()
                     ortholog_alignment = ortholog.sequences[1].sequence.upper()
                     coverage = 100* float (ortholog.get_min_coverage_of_full_sequence_segment())
                     identity = ortholog.identities
                    except:
                     if verbose: sys.stdout.write("\t\t\t--Skip %s ...\n"%ortholog.sequenceID)
                     continue
                    if verbose: sys.stdout.write("\t\t\t--Threading %s ...\n"%(ortholog.sequenceID))
                    homologs.setdefault(ortholog_name, ortholog)
                    # Get the correct aminoacid number
                    first = pdb_obj.get_chain_by_id(pdb_chain).gapped_protein_sequence.find(template_alignment.replace("-", ""))
                    if first<0:first=0
                    positions = range(first, len(template_alignment.replace("-", "")) + first)
                    # Thread ortholog sequence
                    protein = {}
                    for position in range(len(template_alignment)):
                      # Skip if hit alignment position is gapped #
                      if template_alignment[position] != "-":
                        sequence_position = positions.pop(0)
                        if sequence_position in sequence_to_crystal:
                           aminoacid_number = sequence_to_crystal[sequence_position]
                           protein.setdefault((pdb_chain,aminoacid_number), ortholog_alignment[position])
                    # Thread DNA
                    dna  = {}
                    dnae = {}
                    dna_seq_dummy=members[i].get_sequence()
                    if len(dna_seq_dummy) != interface_obj.get_interface_length():
                       if  len(dna_seq_dummy)  < interface_obj.get_interface_length():
                             dna_seq_interface               = dna_seq_dummy + "N"*(interface_obj.get_interface_length()-len(dna_seq_dummy))
                             dna_seq_interface_complementary = triads.get_complementary_dna_sequence(dna_seq_interface) + "N"*(interface_obj.get_interface_length()-len(dna_seq_dummy))
                             if (members[i].get_start()+interface_obj.get_interface_length()) <= len(dna_sequence):
                                   dna_seq_extended          = dna_sequence[members[i].get_start():members[i].get_start()+interface_obj.get_interface_length()]
                             else:
                                   dna_seq_extended          = dna_sequence[members[i].get_start():] 
                             dna_seq_extended_complementary  = triads.get_complementary_dna_sequence(dna_seq_extended)
                             
                       if  len(dna_seq_dummy)  > interface_obj.get_interface_length():
                             dna_seq_interface               = dna_seq_dummy
                             dna_seq_interface_complementary = triads.get_complementary_dna_sequence(dna_seq_dummy)
                             dna_seq_extended                = dna_seq_dummy[0:interface_obj.get_interface_length()]
                             dna_seq_extended_complementary  = triads.get_complementary_dna_sequence(dna_seq_extended)
                    else:
                        dna_seq_interface               = dna_seq_dummy
                        dna_seq_interface_complementary = triads.get_complementary_dna_sequence(dna_seq_dummy)
                        dna_seq_extended                = dna_seq_dummy
                        dna_seq_extended_complementary  = triads.get_complementary_dna_sequence(dna_seq_extended)
                    if members[i].get_strand() == "+":
                        dna.setdefault(dna_seq_interface,interface_obj.get_start())
                        dnae.setdefault(dna_seq_extended,interface_obj.get_start())
                    else:
                        dna.setdefault(dna_seq_interface_complementary,interface_obj.get_start())
                        dnae.setdefault(dna_seq_extended_complementary,interface_obj.get_start())
                    # Create threaded object #
                    threaded_ortholog.set_pdb_name(pdb_name)
                    threaded_ortholog.set_pdb_chain(pdb_chain)
                    threaded_ortholog.set_dna_helix(dna_helix)
                    threaded_ortholog.set_identity(identity)
                    threaded_ortholog.set_coverage(coverage)
                    threaded_ortholog.set_protein(protein)
                    threaded_ortholog.set_dna(dna)
                    threaded_ortholog.set_dna_fixed(dnae)
                    threaded_ortholog.set_query_ali(ortholog_alignment)
                    threaded_ortholog.set_pdb_ali(template_alignment)

                    # Create threaded file #
                    threaded_ortholog.write(threading_file)
                
                    # Add to ortholog file #
                    if verbose: sys.stdout.write("\t\t\t--Add threading file %s ...\n"%(os.path.basename(threading_file)))
                    orthologs[-1]['orthologs'][i].append([threading_file, 0, 0])
 # Save members and clusters
 dump_homologs=open(homologs_storage,"wb")
 if verbose:sys.stdout.write("\t-- Save HOMOLOGS %s\n"%homologs_storage)
 cPickle.dump(homologs,dump_homologs)
 dump_homologs.close()
 dump_orthologs=open(orthologs_storage,"wb")
 if verbose:sys.stdout.write("\t-- Save ORTHOLOGS %s\n"%orthologs_storage)
 cPickle.dump(orthologs,dump_orthologs)
 dump_orthologs.close()
else:
 # Load  members and clusters
 if verbose:sys.stdout.write("\t-- Load HOMOLOGS %s\n"%homologs_storage)
 dump_homologs=open(homologs_storage,"rb")
 homologs=cPickle.load(dump_homologs)
 dump_homologs.close()
 if verbose:sys.stdout.write("\t-- Load ORTHOLOGS %s\n"%orthologs_storage)
 dump_orthologs=open(orthologs_storage,"rb")
 orthologs=cPickle.load(dump_orthologs)
 dump_orthologs.close()

print("total process time: %s seconds" % str(time() - start_time))

#######################
# 6. Score orthologs  #
#######################
if verbose: sys.stdout.write("--Score orthologous transcription factors...\n")
# Initialize #
families = {}
# For each line... #
for line in functions.parse_file(os.path.join(pdb_dir, "families.txt")):
    if line.startswith("#"):continue
    pdb_chain, family = line.split(";")
    families[pdb_chain] = family
if options.split_potential is None:
   split_potential = config.get("Parameters", "split_potential")
if os.path.exists(os.path.join(output_dir,  "orthologs.json")) and verbose: sys.stdout.write("\t--Use existing 'orthologs.json' file\n")
if not os.path.exists(os.path.join(output_dir,  "orthologs.json")):
    # For each cluster... #
    for i in range(len(clusters)):
        # Get cluster #
        cluster = clusters[i]
        if verbose: sys.stdout.write("\t--Score on fragment %d-%d centered on %s\n"%(cluster['from'],cluster['to'],cluster['hit_name']))
        # For each pair of members... #
        for members in cluster['members']:
            # For each monomer... #
            for j in range(2):
                # Skip if monomer is None #
                if members[j] is None: continue
                # Get fold/PDB/chain name
                m = re.search("(\S+)_(\d+)(\S+)",members[j].get_hit())
                if verbose: sys.stdout.write("\t-- Evaluate hit %s\n"%members[j].get_hit())
                if external is not None:
                   if not associated_pdb.has_key(members[j].get_hit()): continue
                   fold_hit = associated_pdb.get(members[j].get_hit())
                   helix_hit= associated_helix.get(members[j].get_hit())
                   pdb_name = fold_hit.split("_")[0]
                   pdb_chain= fold_hit.split("_")[1]
                else:
                   # Get chain PDB obj #
                   m = re.search("(\S+)_(\d+)(\S+)",members[j].get_hit())
                   if m:
                      pdb_name  = m.group(1)
                      helix_hit = m.group(2)
                      pdb_chain = m.group(3)
                   else:
                      pdb_name  = members[j].get_hit()[:4]
                      helix_hit = members[j].get_hit()[5:6]
                      pdb_chain = members[j].get_hit()[6:]
                   fold_hit  = pdb_name+"_"+pdb_chain
                # Get x3dna info #
                if verbose: sys.stdout.write("\t\t--Get DNA object %s...\n"%pdb_name)
                if not os.path.exists(os.path.join(pdb_dir, "x3dna", pdb_name + ".txt")):
                   if verbose: sys.stdout.write("\t\t\t--Skip (not found)...\n")
                   continue
                x3dna_obj = x3dna.X3DNA(os.path.join(pdb_dir, "x3dna", pdb_name + ".txt"))
                # Get DSSP object #
                if options.verbose: sys.stdout.write("\t\t--Get DSSP %s...\n"%pdb_name)
                if not os.path.exists(os.path.join(pdb_dir, "dssp",  pdb_name+ ".txt")):
                   if verbose: sys.stdout.write("\t\t\t--Skip (not found)...\n")
                   continue
                dssp_obj = dssp.DSSP(os.path.join(pdb_dir, "dssp",  pdb_name+ ".txt"))
                # Get contacts #
                if options.verbose: sys.stdout.write("\t\t--Get contacts %s...\n"%pdb_name)
                if not os.path.exists(os.path.join(pdb_dir, "contacts", pdb_name + ".txt" )):
                   if verbose: sys.stdout.write("\t\t\t--Skip (not found)...\n")
                   continue
                contacts_obj = contacts.Contacts(os.path.join(pdb_dir, "contacts", pdb_name + ".txt"))
                # Get PDB object
                if options.verbose: sys.stdout.write("\t\t--Get PDB %s...\n"%fold_hit)
                if not os.path.exists(os.path.join(pdb_dir, "split", fold_hit + ".pdb")):
                   if verbose: sys.stdout.write("\t\t\t--Skip (not found)...\n")
                   continue
                pdb_obj = PDB(os.path.join(pdb_dir, "split", fold_hit + ".pdb"))
                # Add the DNA chain of the corresponding helix
                dna_pdb_file = pdb_name+".dna."+helix_hit+".pdb"
                if options.verbose: sys.stdout.write("\t\t--Get PDB %s...\n"%dna_pdb_file)
                if not os.path.exists(os.path.join(pdb_dir, "split",dna_pdb_file)):
                   if verbose: sys.stdout.write("\t\t\t--Skip (not found)...\n")
                   continue
                dna_pdb_obj  = PDB(os.path.join(pdb_dir, "split",dna_pdb_file))
                for dna_chain_obj in dna_pdb_obj.nucleotides:
                   pdb_obj.add_chain(dna_chain_obj)
                # Get triads object #
                triads_file = fold_hit + ".txt"
                if options.verbose: sys.stdout.write("\t\t--Get triads %s...\n"%triads_file)
                if not os.path.exists(os.path.join(pdb_dir, "triads",triads_file)):
                   if verbose: sys.stdout.write("\t\t\t--Skip (not found)...\n")
                   continue
                triads_obj = triads.Triads(os.path.join(pdb_dir, "triads",triads_file))
                #triads_obj = triads.get_triads_obj(pdb_obj, dssp_obj, x3dna_obj, contacts_obj, complementary=True)
                # Get potential #
                if options.verbose:sys.stdout.write("\t\t--Load potentials...\n")
                try:
                    #potentials, thresholds, radii, structural_homologs_by_chain = PWM.load_statistical_potentials(pdb_obj, pdb_dir, pbm_dir, families, options.radius, None, split_potential, True, False, False, None ,False, False, False, True, None, dummy_dir, verbose)
                    potentials, thresholds, radii, structural_homologs_by_chain = PWM.load_statistical_potentials(pdb_obj, pdb_dir, pbm_dir, families, options.radius, options.potential_file, split_potential, auto_mode, options.family_potentials, options.pbm_potentials, options.score_threshold, options.taylor_approach, options.pmf_potentials , options.bins, known, None, dummy_dir,options.verbose)
                except:
                    if options.verbose:sys.stdout.write("\t\t--WARNING: Failed to use the potential in automated mode\n")
                    continue
                potential=potentials[pdb_chain]
                radius=0
                if radii.has_key(pdb_chain): radius=radii.get(pdb_chain)
                if radius <= 0: radius=float(config.get("Parameters", "max_contact_distance"))
                # For each ortholog... 
                for ortholog in orthologs[i]['orthologs'][j]:
                    # Skip energy ranking if option rank was not selected
                    if not rank: continue
                    # Skip if incorrect ortholog #
                    if external is not None:
                        if not associated_pdb.has_key(members[j].get_hit()): continue
                        fold_hit = associated_pdb.get(members[j].get_hit())
                        helix_hit= associated_helix.get(members[j].get_hit())
                        pdb_name = fold_hit.split("_")[0]
                        pdb_chain= fold_hit.split("_")[1]
                        if pdb_name not in ortholog[0]: 
                            if options.verbose:sys.stdout.write("\t\t--Skip hit %s (PDB: %s) out of %s\n"%(members[j].get_hit(), pdb_name, os.path.basename(ortholog[0])))
                            continue
                    else:
                        if members[j].get_hit() not in ortholog[0]: 
                           if options.verbose:sys.stdout.write("\t\t--Skip hit %s out of %s\n"%(members[j].get_hit(), os.path.basename(ortholog[0])))
                           continue   
                    # Get threaded object #
                    threaded_obj = threader.Threaded(os.path.join(output_dir, "aux_files", ortholog[0]))
                    # Create new contact files #
                    for kmer in threaded_obj.get_kmers_fixed():
                        # Get interface #
                        interface = threaded_obj.get_kmer_fixed_interface(kmer)
                        if options.verbose:sys.stdout.write("\t\t--Check %s with DNA %s (binding %s) using hit %s\n"%( os.path.basename(ortholog[0]),kmer,interface,orthologs[i]['hit_name']))
                        # Skip if k-mer has no interface #
                        if interface is None: continue
                        # For each contact... #
                        for triad_obj in triads_obj.get_triads():
                            # Initialize #
                            a_oa, b_ob, distance, residue_A, residue_B = triad_obj.return_as_string().split(";")
                            protein = a_oa
                            dna = b_ob
                            chain, residue_pos = residue_A.split("-")
                            if chain != pdb_chain: continue
                            # Skip if far contact... #
                            dab = np.floor(float(distance))
                            if dab > radius: continue
                            residue_num=int(residue_pos)
                            bp_1F,bp_1R,bp_2F,bp_2R=residue_B.split(",")
                            # Skip if environment could not be obtained for amino acid #
                            if "None" in a_oa: continue
                            # Skip if environment could not be obtained for dinucleotide #
                            if "None" in b_ob: continue
                            # The following arguments arguments are defined as in the papers #
                            dab = np.floor(float(distance))
                            a, hydrophobicity, degree_of_exposure, secondary_structure = a_oa.split("-")
                            aminoacid_three_letters = a
                            oa = "%s-%s-%s" % (hydrophobicity, degree_of_exposure, secondary_structure)
                            b, nitrogenous_bases, dna_strand, dna_groove, dna_chemical_group = b_ob.split("-")
                            dinucleotide_sequence = b
                            ob = "%s-%s-%s-%s" % (nitrogenous_bases, dna_strand, dna_groove, dna_chemical_group)
                            a_b = "%s;%s" % (a, b)
                            a_b_oa_ob = "%s;%s" % (a_oa, b_ob)
                            oa_ob = "%s;%s" % (oa, ob)
                            ### Compute potentials for template PDB structure (protein and DNA) ###
                            try:
                             template_score = potential.get_score(split_potential, key=a_b_oa_ob, distance=dab)
                            except:
                             template_score=None
                            ### Compute potentials for ortholog ### 
                            # Skip if threaded residue does not exist #
                            if not threaded_obj.has_aminoacid(pdb_chain, residue_num): continue
                            # Skip if threaded residue is gapped or weird #
                            if threaded_obj.get_threaded_aminoacid(pdb_chain, residue_num) not in aminoacids1to3 or threaded_obj.get_threaded_aminoacid(pdb_chain, residue_num) == "X": continue
                            #Change amino-acid
                            aminoacid_single_letter = threaded_obj.get_threaded_aminoacid(pdb_chain, residue_num)
                            aminoacid_three_letters = aminoacids1to3[aminoacid_single_letter]
                            aminoacid_environment   = triads.AminoAcidEnvironment(aminoacid_single_letter,degree_of_exposure,secondary_structure)
                            hydrophobicity          = aminoacid_environment.get_hydrophobicity()
                            # Get new protein part of triad #
                            protein = aminoacid_three_letters + "-" + hydrophobicity + "-" + degree_of_exposure + "-" + secondary_structure
                            # Get dinucleotide of binding site #
                            chain_dna,base_number_1=bp_1F.split("-")
                            chain_dna,base_number_2=bp_2F.split("-")
                            basepair_1 = x3dna_obj.get_residue_basepair(chain_dna,int(base_number_1))
                            basepair_2 = x3dna_obj.get_residue_basepair(chain_dna,int(base_number_2))
                            dinucleotide = x3dna_obj.get_inverse_dinucleotide(basepair_1,basepair_2)
                            # Change dinucleotide
                            dinucleotide_sequence = ""
                            if x3dna_obj.get_dinucleotide(dinucleotide) == None: continue
                            for basepair in x3dna_obj.get_dinucleotide(dinucleotide):
                                # Skip if threaded basepair does not exist #
                                if threaded_obj.has_basepair_fixed(kmer, basepair):
                                    dinucleotide_sequence += threaded_obj.get_threaded_nucleotide_fixed(kmer, basepair)
                            # Skip if k-mer does not cover the whole dinucleotide #
                            if len(dinucleotide_sequence) != 2: continue
                            nucleotide_environment = triads.DinucleotideEnvironment(dinucleotide_sequence, dna_strand, dna_groove, dna_chemical_group)
                            nitrogenous_bases      = nucleotide_environment.get_nitrogenous_bases()
                            # Get new dna part of triad #
                            dna = dinucleotide_sequence + "-" + nitrogenous_bases + "-" + dna_strand + "-" + dna_groove + "-" + dna_chemical_group
                            a_b_oa_ob = "%s;%s" % (protein, dna)
                            # Get score #
                            try:
                             score = potential.get_score(split_potential, key=a_b_oa_ob, distance=dab)
                             if template_score is not None: dscore= score - template_score
                            except:
                             score=None
                            if score is not None:
                                ortholog[1] += score
                                ortholog[2] += dscore
                        if options.verbose:sys.stdout.write("\t\t\t--Score energy %10.5f (differential vs template %10.5f)\n"%(ortholog[1],ortholog[2]))
                        #break
    # Create json file #
    out = open(os.path.join(output_dir, "orthologs.json"), "wt")
    out.write(json.dumps(orthologs, separators=(',', ':'), indent=2))
    out.close()
print("total process time: %s seconds" % str(time() - start_time))


#######################
# 7. Create outputs   #
#######################
if verbose: sys.stdout.write("--Create output by folds and clean unnecessary files...\n")
# Initialize #
max_ratio_loss_scan = config.get("Parameters", "max_ratio_loss_scan")
threaded = set()
output_file = os.path.join(output_dir, "orthologs_with_best_templates.txt")
orthologs_best_json = os.path.join(output_dir, "orthologs_with_best_templates.json")
orthologs = json.loads(''.join([line for line in functions.parse_file(os.path.join(output_dir, "orthologs.json"))]))
# Store a pickled dictionary #
#out_pickle = os.path.join(output_dir, "orthologs.p")
#pickle.dump(orthologs, open(out_pickle, "wb"))

# Create output file #
out = open(output_file, "wt")
# For each cluster... #
for i in range(len(clusters)):
    # Get cluster #
    cluster = clusters[i]
    # Get cluster main hit #
    hit = fimo_obj.get_hit(cluster['hit_name'])
    # Print header #   
    try:
       out.write(">%s|%s-%s|%s|%s\n" % (str(i + 1), str(start_index+cluster['from']), str(start_index+cluster['to']), str(hit.get_p_value()),cluster['hit_name']))
    except:
       if verbose: sys.stdout.write("\t--Skip hit %s\n"%cluster['hit_name'])
       continue
    # Initialize #
    unique_orthologs = [{},{}]
    # For each monomer... #
    if verbose: sys.stdout.write("\t--Cluster %d on  %s \n"%(i,cluster['hit_name']))
    try:
      if verbose: sys.stdout.write("\t--Select output on hit %s (%d orthologs)\n"%(cluster['hit_name'],len(orthologs[i]['orthologs'][0])+len(orthologs[i]['orthologs'][1])))
    except:
      sys.stderr.write("ERROR: Inconsistency of files. Delete file 'orthologs.json' and run again\n")
      info=open(info_file,"a")
      info.write("%s\tFAIL\n"%(os.path.basename(input_file)))
      info.close()
      exit(0)
    #for j in range(len(orthologs[i]['orthologs'])):
    for j in range(2):
        # For each ortholog... #
        for ortholog in orthologs[i]['orthologs'][j]:
            m = re.search("^%s/(\S+).\S{4}_(\S+).\d+-\d+.txt$" % os.path.join(output_dir, "aux_files"), ortholog[0] )
            if m:
                if float(ortholog[2]) > 0:
                  if abs(float(ortholog[1])-float(ortholog[2]))> 0.01:
                     ratio=float(ortholog[2])/abs(float(ortholog[1])-float(ortholog[2]))
                  else:
                     continue
                else:
                   ratio=0
                if ratio > max_ratio_loss_scan: continue
                if unique_orthologs[j].has_key(m.group(1)):
                  # If better score #
                  if verbose: sys.stdout.write("\t\t--Check %s score %10.5f (total number of templates %d)\n"%(os.path.basename(ortholog[0]),float(ortholog[1]),len(unique_orthologs[j][m.group(1)])))
                  if float(ortholog[1]) < float(unique_orthologs[j][m.group(1)][0][1]):
                    unique_orthologs[j][m.group(1)]=[]
                    unique_orthologs[j][m.group(1)].append(ortholog)
                  elif abs(float(ortholog[1]) - float(unique_orthologs[j][m.group(1)][0][1])) < 0.1:
                    unique_orthologs[j][m.group(1)].append(ortholog)
                  else:
                    if verbose: sys.stdout.write("\t\t\t Skip %s from the best\n"%(os.path.basename(ortholog[0])))
                else:
                   if verbose: sys.stdout.write("\t\t--Check %s score %10.5f (first template)\n"%(os.path.basename(ortholog[0]),float(ortholog[1])))
                   unique_orthologs[j].setdefault(m.group(1),[]).append(ortholog)
            else:
                if verbose: sys.stdout.write("ortholog %s does not match regular expression!\n"%ortholog[0])
    # Initialize #
    matrix = [[[None], [None]] for j in range(max([len([unique_orthologs[k].iterkeys()]) for k in range(2)]))]
    # For each monomer... #
    #for j in range(len(orthologs[i]['orthologs'])):
    for j in range(2):
        scores = set([unique_orthologs[j][ortholog][0][1] for ortholog in sorted(unique_orthologs[j], key=lambda x: unique_orthologs[j][x][0][1])])
        # For each score... #
        for score in sorted(scores):
            for database in ["sp", "tr","SP","TR"]:
                # For each ortholog... #
                for ortholog in unique_orthologs[j].iterkeys():
                    # Skip if ortholog does not belong to database or has a different score #
                    try:
                       header = homologs[ortholog].sequenceID
                    except:
                       continue
                    m = re.search("^(tr|sp|TR|SP)\|(\S+)\|(\S+)\_(\S+)", header)
                    if m:
                       ortholog_database = m.group(1)
                       ortholog_name     = m.group(2)
                       ortholog_genename = m.group(3)
                       ortholog_specie   = m.group(4)
                    if ortholog_database != database or unique_orthologs[j][ortholog][0][1] != score or ortholog_name != ortholog : continue
                    for unique in unique_orthologs[j][ortholog]:
                      # For each element in matrix... #
                      for k in matrix:
                         data="%s;%s;%s;%s;%s" % (ortholog, ortholog_genename , unique[0][len(output_dir) + 1  :], str(unique[1]), str(unique[2]))
                         k[j].append(data)
                         threaded.add(unique[0][len(output_dir) + 1:])
    for k in matrix:
      for kn in k[0]:
        for km in k[1]:
          if kn is None and km is None and len(k[0])>1 and len(k[1])>1: continue
          out.write("%s\t%s\n" % (kn,km))
    out.write("//\n")
out.close()

orthologs_with_best_templates_list = parse_best_orthologs(output_file,rank)
# Create json file #
out = open(orthologs_best_json, "wt")
out.write(json.dumps(orthologs_with_best_templates_list, separators=(',', ':'), indent=2))
out.close()


output_file = os.path.join(output_dir, "orthologs.txt")
# Create output file #
out = open(output_file, "wt")
# For each cluster... #
for i in range(len(clusters)):
    # Get cluster #
    cluster = clusters[i]
    # Get cluster main hit #
    hit = fimo_obj.get_hit(cluster['hit_name'])
    # Print header #   
    out.write(">%s|%s-%s|%s|%s\n" % (str(i + 1), str(start_index+cluster['from']), str(start_index+cluster['to']), str(hit.get_p_value()),str(cluster['hit_name'])))
    # Initialize #
    all_orthologs = [{},{}]
    # For each monomer... #
    if verbose: sys.stdout.write("\t--Add orthologs for hit %s (%d orthologs)\n"%(cluster['hit_name'],len(orthologs[i]['orthologs'][0])+len(orthologs[i]['orthologs'][1])))
    #for j in range(len(orthologs[i]['orthologs'])):
    for j in range(2):
        # For each ortholog... #
        for ortholog in orthologs[i]['orthologs'][j]:
            if float(ortholog[1]) >= 0 and ortholog[0][len(output_dir) + 1:] not in  threaded: continue
            m = re.search("^%s/(\S+).\S{4}_(\S+).\d+-\d+.txt$" % os.path.join(output_dir, "aux_files"), ortholog[0] )
            if m:
                if verbose: sys.stdout.write("\t\t--  %s score %10.5f \n"%(os.path.basename(ortholog[0]),float(ortholog[1])))
                all_orthologs[j].setdefault(m.group(1),[]).append(ortholog)
            else:
                if verbose: sys.stdout.write("ortholog %s does not match regular expression!\n"%ortholog[0])
    # Initialize #
    matrix = [[[None], [None]] for j in range(max([len([all_orthologs[k].iterkeys()]) for k in range(2)]))]
    # For each monomer... #
    for j in range(2):
        scores=set( [ ortho[1] for ortholog in sorted(all_orthologs[j], key=lambda x: all_orthologs[j][x][0][1]) for ortho in all_orthologs[j][ortholog] ] )
        # For each score... #
        for score in sorted(scores):
            for database in ["sp", "tr","SP","TR"]:
                # For each ortholog... #
                for ortholog in all_orthologs[j].iterkeys():
                    # Skip if ortholog does not belong to database or has a different score #
                    header = homologs[ortholog].sequenceID
                    m = re.search("^(tr|sp|TR|SP)\|(\S+)\|(\S+)\_(\S+)", header)
                    if m:
                       ortholog_database = m.group(1)
                       ortholog_name     = m.group(2)
                       ortholog_genename = m.group(3)
                       ortholog_specie   = m.group(4)
                    for ortho in all_orthologs[j][ortholog]:
                      if ortholog_database != database or ortho[1] != score or ortholog_name != ortholog : continue
                      # For each element in matrix... #
                      for k in matrix:
                          data="%s;%s;%s;%s;%s" % (ortholog, ortholog_genename , ortho[0][len(output_dir) + 1  :], str(ortho[1]), str(ortho[2]))
                          k[j].append(data)
                          threaded.add(ortho[0][len(output_dir) + 1:])
    for k in matrix:
      for kn in k[0]:
        for km in k[1]:
          if (kn is None  and len(k[0])>1) and (km is None and len(k[1])>1): continue
          out.write("%s\t%s\n" % (kn,km))
    out.write("//\n")
out.close()


# For each file check and clean... #
cleanset=set()
for threading_file in os.listdir(os.path.join(output_dir, "aux_files")):
    m = re.search("^(\S+.\S{4}_\S+.\d+-\d+).txt$", os.path.basename(threading_file))
    # If threading file is not required among selected skip cleaning... #
    if m:
       save=False 
       for file_thread in threaded:
         if threading_file in cleanset: continue
         if os.path.basename(threading_file) in file_thread: save=True 
       if not save: 
           if verbose: sys.stdout.write("\t--Clean unused threaded file %s %s\n"%(threading_file,file_thread))
           cleanset.add(threading_file)
           os.remove(os.path.join(output_dir, "aux_files",threading_file))
print("total process time: %s seconds" % str(time() - start_time))

if not options.cluster_complexes:
    print("No clustering required")
    print("Done")
    print("total process time: %s seconds" % str(time() - start_time))
    info=open(info_file,"a")
    info.write("%s\tDONE\n"%(os.path.basename(input_file)))
    info.close()
    exit(0)
else:
    print("\n######################### CONTINUE WITH THE CLUSTERING OF CLOSE BINDING COMPLEXES #################################\n")


#########################################################
# 8. Get FIMO EXTENDED OBJECT and group the clusters    #
#########################################################
max_complex_connect = int(config.get("Parameters", "max_complex_connect"))
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
    complex_tf=chain_hit+"_"+helix_hit+"_"+str(start_index+hit.get_start())+"_"+str(start_index+hit.get_end())
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
       name_cluster =   chains + "_" + helices + "_" + str(start_index+start) + "_" + str(start_index+end)
       cluster_hits.setdefault(name_cluster,[])
       if verbose: sys.stdout.write("\t\t-- Rename cluster %s  as %s \n"%(key,name_cluster))
       for i in complex_number:
           complex_hit =  complex_list[i]
           cluster_hits[name_cluster].append(complex_hit)
           if verbose: sys.stdout.write("\t\t\t-- Add hit %s  to %s  \n"%(complex_hit,name_cluster))

    #Fill the complex_cluster_hits dictionary
    complex_cluster_hits.setdefault(pdb_hit,cluster_hits)
          


  orthologs_complexes=[]
  complex_threaded={}
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
####################################
# 9. Create outputs for complexes  #
####################################
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
     for i in range(len(clusters)):
          cluster         = clusters[i]
          range_cluster   = set( range(int(start_index+cluster['from']),int(start_index+cluster['to'])+1))
          ortholog        = orthologs[i]
          cluster_hits    = set()
          cluster_threads = set()
          for j in range(2):
              for members in cluster['members']:
                  if members[j] is not None:
                     cluster_hits.add(members[j].get_hit())
              for thread_file_data in ortholog['orthologs'][j]:
                  thread_file=thread_file_data[0]
                  if thread_file is not None:
                     cluster_threads.add(os.path.basename(thread_file))
          #if verbose: sys.stdout.write("\t\t-- Compare with cluster [%d] hits: %s\n"%(i+1,list(cluster_hits)))
          #if verbose: sys.stdout.write("\t\t-- Compare with cluster [%d] threads: %s\n"%(i+1,list(cluster_threads)))
          #if verbose: sys.stdout.write("\t\t-- Compare with cluster [%d] range: %s\n"%(i+1,list(range_cluster)))
          #if verbose: sys.stdout.write("\t\t-- Common range: %s Common hits: %s Common threads: %s \n"%(range_complex.intersection(range_cluster),complex_hits.intersection(cluster_hits),complex_threads.intersection(cluster_threads)))
          if len( range_complex.intersection(range_cluster) ) >0 and  len( complex_threads.intersection(cluster_threads) ) >0 :
             if verbose: sys.stdout.write("\t\t--Add cluster %d (by hit %s in [%d-%d])...\n"%(i,cluster['hit_name'],int(start_index+cluster['from']),int(start_index+cluster['to'])))
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
            cluster = clusters[nclust-1]
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
print("total process time: %s seconds" % str(time() - start_time))
info=open(info_file,"a")
info.write("%s\tDONE\n"%(os.path.basename(input_file)))
info.close()
print("Done")
exit(0)

