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
tf_dir_global=os.path.join(home,"nearest_neighbour","TF_MODELS")
python=os.path.join(config.get("Paths", "python_path"), "python")
dummy=os.path.join(home,"dummy_pwm_tf")
pdb_dir=os.path.join(home,"pdb")
out_dir_global=os.path.join(home,"nearest_neighbour","PWM_MODELS")
pbm_dir=os.path.join(home,"pbm_2.0")
radius=30.0

list_ratio=[1.0,0.8,0.6,0.4,0.2,0.1]
#list_ratio=[1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]:
for ratio in list_ratio:
  if ratio==1.0:
     out_dir = out_dir_global
     tf_dir  = tf_dir_global 
  else:
     out_dir = out_dir_global + "_" + str(ratio)
     tf_dir  = tf_dir_global + "_" + str(ratio)
  if not os.path.exists(out_dir):
     os.makedirs(out_dir)
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
      for pre in list_ratio:
          if ratio >= pre : continue
          if pre==1.0: 
             pre_dir  = out_dir_global 
          else:
             pre_dir  = out_dir_global + "_" + str(pre)
          pre_file = os.path.join(pre_dir,tf_seq.rstrip(".pdb"))
          pre_pwm  = pre_file+".pwm"
          pre_meme = pre_file+".meme"
          pre_msa  = pre_file+".msa"
          pre_logo = pre_file+".logo"
          pre_fwd_logo = pre_logo+".fwd.png"
          pre_rev_logo = pre_logo+".rev.png"
          if os.path.exists(pre_pwm): shutil.copy(pre_pwm,output_pwm)
          if os.path.exists(pre_meme ): shutil.copy(pre_meme ,output_meme )  
          if os.path.exists(pre_msa  ): shutil.copy(pre_msa ,output_msa )  
          if os.path.exists(pre_fwd_logo ): shutil.copy(pre_fwd_logo ,output_fwd_logo )  
          if os.path.exists(pre_rev_logo ): shutil.copy(pre_rev_logo ,output_rev_logo )  

      if os.path.exists(output_pwm) and os.path.exists(output_meme):
            skip=True
      if not skip:
        if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
        else: cluster_queue=config.get("Cluster", "cluster_queue")
        print "%s %s --radius=%f -i %s -o %s --pdb=%s --pbm=%s --auto -v --dummy=%s " % ( os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "pwm_pbm.py"), radius, input_file, output_file, os.path.abspath(pdb_dir), os.path.abspath(pbm_dir) , os.path.abspath(dummy))
        functions.submit_command_to_queue("%s %s --radius=%f -i %s -o %s --pdb=%s --pbm=%s --auto -v --dummy=%s " % ( os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "pwm_pbm.py"), radius, input_file, output_file, os.path.abspath(pdb_dir), os.path.abspath(pbm_dir) , os.path.abspath(dummy) ) , cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),dummy,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))


