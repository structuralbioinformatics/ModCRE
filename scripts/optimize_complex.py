import os, sys, re
import ConfigParser
import optparse
import shutil
import subprocess
import numpy as np

scripts_path = os.path.abspath(os.path.dirname(__file__))
sys.path.append(scripts_path)

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Get python and SRC paths in 2.7 #
python     = os.path.join(config.get("Paths", "python_path"), "python")
modeller   = config.get("Paths", "modeller_path")
src_path   = config.get("Paths", "src_path")
modpy_path = config.get("Paths", "modpy_path")
 
if not os.path.exists(modpy_path):
       modpy_path    =  os.path.join(config.get("Paths","src_path"),config.get("Paths", "modpy_path"))


# Import my functions #
from functions import submit_command_to_queue
from functions import done_jobs
from functions import check_done
from functions import parse_fasta_file

#Initialize
from SBI.structure.chain import Chain
from SBI.structure.chain import ChainOfProtein
from SBI.structure.chain import ChainOfNucleotide
from SBI.sequence import Sequence
from SBI.structure import PDB
from model_protein import remove_files
from model_protein import renumber_pdb
from model_protein import create_mapping



def parse_arguments():  # This function passes the arguments given by the user.

    parser = optparse.OptionParser()
    parser.add_option('-d', action="store",  dest="output_folder", default=None, help="Path to the folder with the Complexes.")
    parser.add_option("--parallel", default=False, action="store_true", dest="parallel", help="Run in parallel (default = False)")
    parser.add_option("--reuse", default=False, action="store_true", dest="reuse", help="Use previous optimizations (default = False)")

    options,args = parser.parse_args()

    if options.output_folder is None  :
        parser.error("missing arguments: type option \"-h\" for help")


    return options


def clean(options):

 print("Clean folder "+options.output_folder)
 for output_fragment in os.listdir(options.output_folder):
   if output_fragment.startswith("fragment") and os.path.isdir(os.path.join(options.output_folder,output_fragment)):

      print("\t-- Fragment "+output_fragment)
      output_dir = os.path.join(options.output_folder,output_fragment)
      dummy_dir  = os.path.join(options.output_folder,output_fragment,"dummy_optimization")
      shutil.rmtree(dummy_dir)
      for pdb_file in os.listdir(os.path.join(options.output_folder,output_fragment)):
        if pdb_file.endswith(".pdb") and pdb_file.startswith("dna"):
            if pdb_file.endswith(".non_optimized.pdb"):continue
            code=pdb_file.rstrip("pdb")
            if os.path.exists(os.path.join(output_dir,code+"ini")):
                  os.remove(os.path.join(output_dir,code+"ini"))
            if os.path.exists(os.path.join(output_dir,code+"rsr")):
                  os.remove(os.path.join(output_dir,code+"rsr"))

def renumber(options):

 print("Renumber structures  "+options.output_folder)
 job_dir = os.path.join(os.path.dirname(options.output_folder),"../..")
 sequences_file = os.path.join(job_dir,"pwm_database", "sequences","sequences.fasta")
 uniprot_file   = os.path.join(config.get("Paths", "files_path"),config.get("Paths", "uniprot"))
 fasta_file     = os.path.join(options.output_folder,"sequences.fa")
 if os.path.exists(fasta_file): os.remove(fasta_file)
 for output_fragment in os.listdir(options.output_folder):
   if output_fragment.startswith("fragment") and os.path.isdir(os.path.join(options.output_folder,output_fragment)):

      print("\t-- Fragment "+output_fragment)
      output_dir = os.path.join(options.output_folder,output_fragment)
      dummy_dir  = os.path.join(options.output_folder,output_fragment,"dummy_optimization")
      set_of_models=set()
      for pdb_file in os.listdir(os.path.join(options.output_folder,output_fragment)):
        if pdb_file.endswith(".pdb") and pdb_file.startswith("dna"):
            if pdb_file.endswith(".non_optimized.pdb"):continue
            set_of_models.add(pdb_file)

      #Renumber the chains according to FASTA sequence
      chain2protein      = {}
      proteincodes       = set()
      sequences_input    = {}
      sequences_complex  = {}
      for pdb_file in set_of_models:
          original =  os.path.join(options.output_folder,output_fragment,pdb_file.rstrip("pdb")+"non_optimized.pdb") 
          opt_file =  os.path.join(options.output_folder,output_fragment,pdb_file)
          if not os.path.exists(original): continue
          if not os.path.exists(opt_file): continue
          fo = open(original,"r")
          for line in fo:
              if line.startswith("REMARK"):
                 chain   = line.split("=")[0].split()[-1]
                 protein = ".".join(line.split("=")[1].split()[0].split(":")[0].split(".")[:-2])
                 binding = line.split("=")[1].split()[0].split(":")[0].split(".")[-1]
                 if protein == "" or protein == " ": continue
                 print("\t\t-- Protein "+protein+" in chain "+chain)
                 chain2protein.setdefault(chain,protein)
                 proteincodes.add(protein)
          fo.close()
      n_count=1
      if os.path.exists(sequences_file) and os.path.isfile(sequences_file): 
         print("\t-- Read sequence file "+sequences_file)
         for header,sequence in parse_fasta_file(sequences_file,gz=False,clean=True):
             for protein in proteincodes:
                if protein in set(sequences_input.keys()):continue
                if protein in header: 
                   sequences_input.setdefault(protein,sequence)
         renumerate=True
      elif os.path.exists(uniprot_file) and os.path.isfile(uniprot_file):
         print("\t-- Try to read sequence file "+uniprot_file)
         try:
           for header,sequence in parse_fasta_file(uniprot_file,gz=False,clean=True):
             if np.mod(n_count,1000)==0 : sys.stdout.write(".")
             if np.mod(n_count,50000)==0: sys.stdout.write("\n")
             n_count=n_count+1
             if proteincodes==set(sequences_input.keys()): break 
             for protein in proteincodes:
                if protein in set(sequences_input.keys()):continue
                if protein in header: 
                   sys.stdout.write("\nFound %s\n"%(protein))
                   n_count=1
                   sequences_input.setdefault(protein,sequence)
           renumerate=True
         except:
           print("\t-- Try as GZIP compressed file")
           try:
             for header,sequence in parse_fasta_file(uniprot_file,gz=True,clean=True):
               if np.mod(n_count,1000)==0 : sys.stdout.write(".")
               if np.mod(n_count,50000)==0: sys.stdout.write("\n")
               n_count=n_count+1
               if proteincodes==set(sequences_input.keys()): break 
               for protein in proteincodes:
                 if protein in set(sequences_input.keys()):continue
                 if protein in header: 
                   sys.stdout.write("\nFound %s\n"%(protein))
                   n_count=1
                   sequences_input.setdefault(protein,sequence)
             renumerate=True   
           except:
             renumerate=False 
      elif os.path.exists(uniprot_file+".gz")  and os.path.isfile(uniprot_file+".gz"):
         print("\t-- Try to read sequence file "+uniprot_file+".gz")
         for header,sequence in parse_fasta_file(uniprot_file+".gz",gz=True,clean=True):
             if np.mod(n_count,1000)==0 : sys.stdout.write(".")
             if np.mod(n_count,50000)==0: sys.stdout.write("\n")
             n_count=n_count+1
             if proteincodes==set(sequences_input.keys()): break 
             for protein in proteincodes:
               if protein in set(sequences_input.keys()):continue
               if protein in header:
                  sys.stdout.write("\nFound %s in %s\n"%(protein,header))
                  n_count=1
                  sequences_input.setdefault(protein,sequence)
         renumerate=True    
      else:
         renumerate=False
      fa=open(fasta_file,"a")
      for p,s in sequences_input.iteritems():
        fa.write(">%s\n%s\n"%(p,s))
      fa.close()
      if renumerate:
         pdb_obj=PDB(opt_file)
         sequences_complex  = {}
         for chain_id in pdb_obj.chain_identifiers:
             chain=pdb_obj.get_chain_by_id(chain_id)
             if chain.chaintype != "P": 
                sequences_complex.setdefault(chain_id,(chain_id,None))
                continue
             if chain_id in chain2protein.keys():
                header  = chain2protein[chain_id]
                sequences_complex.setdefault(chain_id,(header,sequences_input.get(header)))
             else:
                sequences_complex.setdefault(chain_id,(chain_id,chain.protein_sequence))
         pdb_renumber=PDB()
         pdb_renumber=renumber_pdb(opt_file,sequences_complex,dummy_dir) 
         if len(pdb_renumber.chain_identifiers)>1:
            pdb_renumber.write(opt_file,force=True)
         else:
            print("Failed to renumber "+pdb_file)

 

def optimize(options):

 print("Open folder "+options.output_folder)
 for output_fragment in os.listdir(options.output_folder):
   if output_fragment.startswith("fragment") and os.path.isdir(os.path.join(options.output_folder,output_fragment)):

      print("\t-- Fragment "+output_fragment)
      output_dir = os.path.join(options.output_folder,output_fragment)
      dummy_dir  = os.path.join(options.output_folder,output_fragment,"dummy_optimization")
      info_file  = os.path.join(output_dir,"optimized.log")

      if not os.path.exists(dummy_dir):
          os.makedirs(dummy_dir)

      reuse=options.reuse
      if os.path.exists(info_file) and reuse:
          print("\t\t-- Reuse previous information of runs from file %s"%info_file)
      else:
          log_file = open(info_file,"w")
          log_file.write("#List of MODELS\n")
          log_file.close()

      submitted=set()
      set_of_models=set()

      for pdb_file in os.listdir(os.path.join(options.output_folder,output_fragment)):
        if pdb_file.endswith(".pdb") and pdb_file.startswith("dna"):
            if pdb_file.endswith(".non_optimized.pdb"):continue
            set_of_models.add(pdb_file)

      done    = done_jobs(info_file)
      iterate = check_done(done,set_of_models)
      n_done  = 0

      while(iterate):
        for pdb_file in set_of_models:
            if pdb_file in submitted: continue
            submitted.add(pdb_file)
            if reuse:
              print("\t\t-- Reusing optimizations")
              log_file = open(info_file,"r")
              skip_adding=False
              for line in log_file:
                       if pdb_file in line.split(): skip_adding=True
              if skip_adding: print("\t\t-- Already in the information file ")
              if skip_adding: submitted.add(pdb_file)
              if skip_adding: continue
              log_file.close()
            output_opt_file = pdb_file.rstrip("pdb")+"opt.pdb"
            parameters      =  " --pdb  "    + os.path.join(output_dir,pdb_file)
            parameters      =  parameters + " --output  " + os.path.join(output_dir,output_opt_file)
            parameters      =  parameters + " --log %s "%info_file
            modpysh         =  os.path.join(modeller,"modpy.sh")
            program         =  os.path.join(modpy_path,"optimize.py")

            if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
            else: cluster_queue=config.get("Cluster", "cluster_queue")

            if options.parallel:
              print("\t\t-- Submit  %s %s" % (program,parameters))
              submit_command_to_queue("%s %s %s %s" % (modpysh,python,program,parameters), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
            else:
              print("\t\t-- Execute  %s %s %s %s" % (modpysh,python,program,parameters))
              os.system("%s %s %s %s" % (modpysh,python,program,parameters))
            submitted.add(pdb_file)

        done    = done_jobs(info_file)
        iterate = check_done(done,set_of_models) 
        if len(done) > n_done:
           n_done=len(done)
           sys.stdout.write("Number of models already done %d\n"%n_done)
           sys.stdout.write("\t\t-- Check files done ...\n")
           log_file = open(info_file,"r")
           for line in log_file:
                   print("\t\t\t-- %s"%line.strip())
           log_file.flush()
           log_file.close()
           sys.stdout.write("\t\t-- Still running protein profiles  %s ...\n"%check_done(done,set_of_models))
           sys.stdout.write("\t\t-- Continue iteration %s ...\n"%iterate)
           sys.stdout.flush()

      for pdb_file in set_of_models:
            output_opt_file = pdb_file.rstrip("pdb")+"opt.pdb"
            print("\t\t-- Rename file "+os.path.join(options.output_folder,output_fragment,output_opt_file))
            if os.path.exists(os.path.join(options.output_folder,output_fragment,output_opt_file)):
              try:
                if not os.path.exists(os.path.join(options.output_folder,output_fragment,pdb_file.rstrip("pdb")+"non_optimized.pdb")):
                   shutil.move(os.path.join(options.output_folder,output_fragment,pdb_file),os.path.join(options.output_folder,output_fragment,pdb_file.rstrip("pdb")+"non_optimized.pdb"))  
                else:
                   os.remove(os.path.join(options.output_folder,output_fragment,pdb_file))
                shutil.move(os.path.join(options.output_folder,output_fragment,output_opt_file),os.path.join(options.output_folder,output_fragment,pdb_file))
              except:
                print("\t\t-- Skip "+os.path.join(options.output_folder,output_fragment,output_opt_file))

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

     options=parse_arguments()
     print("\n*******************\n Optimization\n*******************\n")
     optimize(options)
     print("\n*******************\n Renumbering\n*******************\n")
     renumber(options)
     print("\n*******************\n Cleaning folders\n*******************\n")
     clean(options)
     print("Done")

