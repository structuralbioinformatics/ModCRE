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
tf_dir=os.path.join(home,"nearest_neighbour","TF_MODELS")
python=os.path.join(config.get("Paths", "python_path"), "python")
dummy=os.path.join(home,"dummy_pwm_tf")
pdb_dir=os.path.join(home,"pdb")
out_dir=os.path.join(home,"nearest_neighbour","PWM_MODELS")
if not os.path.exists(out_dir):
   os.makedirs(out_dir)
pbm_dir=os.path.join(home,"pbm_2.0")
radius=30.0

for tf_seq in os.listdir(tf_dir):
    if tf_seq.endswith("pdb"):
      input_file=os.path.join(tf_dir,tf_seq)
      output_file=os.path.join(out_dir,tf_seq.rstrip(".pdb"))
      output_pwm  = output_file+".pwm"
      output_meme = output_file+".meme"
      output_msa  = output_file+".msa"
      output_logo = output_file+".logo"
      output_fwd_logo = output_file+".logo.fwd.png"
      output_rev_logo = output_file+".logo.rev.png"
      skip=False
      if os.path.exists(output_pwm) and os.path.exists(output_meme) and os.path.exists(output_msa) and os.path.exists(output_fwd_logo) and os.path.exists(output_rev_logo):
            skip=True
      if not skip:
        if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
        else: cluster_queue=config.get("Cluster", "cluster_queue")
        print "%s %s --radius=%f -i %s -o %s --pdb=%s --pbm=%s --auto -v --dummy=%s " % ( os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "pwm_pbm.py"), radius, input_file, output_file, os.path.abspath(pdb_dir), os.path.abspath(pbm_dir) , os.path.abspath(dummy))
        functions.submit_command_to_queue("%s %s --radius=%f -i %s -o %s --pdb=%s --pbm=%s --auto -v --dummy=%s " % ( os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "pwm_pbm.py"), radius, input_file, output_file, os.path.abspath(pdb_dir), os.path.abspath(pbm_dir) , os.path.abspath(dummy) ) , cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),dummy,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))


