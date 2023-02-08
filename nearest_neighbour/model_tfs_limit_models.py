import os, sys, re
from Bio import motifs as mm
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import ConfigParser
import gzip
import numpy
import optparse
import shutil
import socket
import subprocess

# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.join(os.path.dirname(__file__),"../scripts"))
# Append scripts path to python path #
sys.path.append(scripts_path)

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Imports my functions #
import functions


home=os.path.join(scripts_path,"..")
tf_dir=os.path.join(home,"nearest_neighbour","sequences")
command=os.path.join(home,"scripts","model_protein_nn.py")
python=os.path.join(config.get("Paths", "python_path"), "python")
dummy=os.path.join(home,"dummy_model_tf")
parameters=" --n-total=100  --renumerate --full --all -v "
pdb_dir=os.path.join(home,"pdb")
out_dir_global=os.path.join(home,"nearest_neighbour","TF_MODELS")

for ratio in [0.9,0.8]:
#for ratio in [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]:
  out_dir = out_dir_global + "_" + str(ratio)
  if not os.path.exists(out_dir):
    os.makedirs(out_dir)
  for tf_seq in os.listdir(tf_dir):
    if tf_seq.endswith("fa"):
        input_file=os.path.join(tf_dir,tf_seq)
        if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
        else: cluster_queue=config.get("Cluster", "cluster_queue")
        print "Submit %s -i %s --pdb=%s --id_threshold %f --dummy=%s -o %s %s"%(command,input_file,pdb_dir,ratio,dummy,out_dir,parameters)
        functions.submit_command_to_queue("%s %s -i %s --pdb=%s --id_threshold %f --dummy=%s -o %s %s"%(python,command,input_file,pdb_dir,ratio, dummy,out_dir,parameters),cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),dummy,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))


