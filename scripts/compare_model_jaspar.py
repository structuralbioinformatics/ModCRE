import os, sys, re
import ConfigParser
import optparse
import shutil
import subprocess
import json
import numpy as np
import pandas as pd
import string
import argparse
from scipy import stats
import cPickle
from collections import Counter



# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.dirname(__file__))

# Append scripts path to python path #
sys.path.append(scripts_path)

# Append scripts path
#scripts_path="/home/boliva/sit_sbi/ModCRE/scripts"
#sys.path.append(scripts_path)

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Get python path #
python = os.path.join(config.get("Paths", "python_path"), "python")

# Import my functions #
import functions

# Import my modules #
import pwm_pbm as PWM
import tomtom as TOMTOM


# Imports jbonet's module #
from SBI.structure import PDB
from SBI.structure.chain import ChainOfProtein, ChainOfNucleotide
from SBI.structure.residue import ResidueOfNucleotide
from SBI.structure.atom import AtomOfNucleotide

# Imports my functions #
import SeqIO as SEQ




def jaspar_tf_mdls( omdltmt, tf_file, pwm_models_dir, jaspar_db, tmp):

      tf_mdl_tomtom={}
      data_info=set()
      fl=open(tf_file,"r")
      for line in fl:
        if line.startswith("#"):continue
        data  = line.split()
        tfid  = data[0]
        tfmdl = data[1]
        motif_cisbp = data[2]
        motif_mdl   = data[3]
        data_info.add((tfid,tfmdl,motif_mdl))
      fl.close()
      for (tfid,tfmdl,motif_mdl) in data_info:
        pwm_model   = os.path.join(pwm_models_dir,motif_mdl)
        try:
               tomtom_obj  = TOMTOM.get_tomtom_obj(jaspar_db, pwm_model, tmp)
        except Exception as e:
               sys.stderr.write("TOMTOM HIT %s %s Error %s\n"%(jaspar_db,pwm_model,e))
               continue    
        tf_mdl_tomtom.setdefault(tfmdl,tomtom_obj)
      outmt=open(omdltmt,"wb")
      cPickle.dump(tf_mdl_tomtom,outmt)
      outmt.close()

      return tf_mdl_tomtom




#-------------#
# Options     #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python compare_model_jaspar.py  --tfid=tfid --nn_motifs=nn_dir  --jaspar_pwms=jaspar_pwms_dir  --models_pwms=models_pwms_dir [--dummy -o output_dir]",
                                   epilog      = '@Oliva\'s lab 2018',
                                   description = "The program uses previous results from nearest neighbours and compares the PWMs of all models of a TFID with JASPAR DB")

    parser.add_option("--dummy", default="./tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("--nn_motifs", action="store", type="string", dest="nn_dir", help="Folder with the output of the nearest neighbour analyses of CisBP", metavar="{directory}")
    parser.add_option("--jaspar_pwms", action="store", type="string", dest="jaspar_pwms_dir", help="Folder with PWMS in meme format of JASPAR", metavar="{directory}")
    parser.add_option("--models_pwms", action="store", type="string", dest="pwm_models_dir", help="Folder with PWMS in meme format of CisBP modelled with ModCRE", metavar="{directory}")
    parser.add_option("--tfid", action="store", type="string", dest="tfid", help="Name of the TFID with models", metavar="{name}")
    parser.add_option("-o", action="store", type="string", dest="output_dir", help="Folder name for outputing the comparisons", metavar="{directory}")
    parser.add_option("-v","--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")

    (options, args) = parser.parse_args()

    if options.pwm_models_dir is None:
        print("Missing models_pwms")
        exit()
    if options.jaspar_pwms_dir is None:
        print("Missing jaspar_pwms")
        exit()
    if options.nn_dir is None:
        print("Missing nn_dir")
        exit()
    if options.output_dir is None:
        print("Missing output_dir")
        exit()

    return options




def fileExist(file):
    '''
    Check existing files
    '''
    if file is not None:
        return os.path.exists(file) and os.path.isfile(file)
    else:
        return False


#-------------#
# Main        #
#-------------#

def main():

    # Arguments & Options #
    options = parse_options()

    try:
     print options
     PWM_jaspar  = options.jaspar_pwms_dir
     if not PWM_jaspar.startswith("/"):
        PWM_jaspar = os.path.abspath(PWM_jaspar)
     print PWM_jaspar
     nn_dir  = options.nn_dir
     if not nn_dir.startswith("/"):
        nn_dir=os.path.abspath(nn_dir)
     print nn_dir
     tfmotifs_dir       = os.path.join(nn_dir,"models") 
     print tfmotifs_dir
     tf_tomtom_dir      = os.path.join(nn_dir,"models_tomtom")
     print tf_tomtom_dir 
     output_dir=options.output_dir
     print output_dir
     tfid=options.tfid
     print tfid
     tfid_jaspar_pickle = os.path.join(output_dir,tfid+".pickle")
     print tfid_jaspar_pickle
     pwm_models_dir     = options.pwm_models_dir
     if not pwm_models_dir.startswith("/"):
         pwm_models_dir=os.path.abspath(pwm_models_dir)
     print pwm_models_dir
    except:
     print("Input error: check missing data")
     exit()

    tmp          = options.dummy_dir
    if not tmp.startswith("/"): tmp= os.path.abspath(options.dummy_dir)
    verbose      = options.verbose


    jaspar_db=os.path.join(PWM_jaspar,"jaspar_db.txt")
    if not fileExist(jaspar_db):
       for motif in os.listdir(PWM_jaspar):
           os.system("cat %s >> %s\n"%(motif,jaspar_db))

    tf_file=os.path.join(tfmotifs_dir,tfid+".dat")

    if options.verbose: print("Run TOMTOM comparison of %s in %s"%(tf_file,jaspar_db))

    tf_mdl_jaspar_tomtom = jaspar_tf_mdls(tfid_jaspar_pickle,tf_file, pwm_models_dir, jaspar_db, tmp)

    if options.verbose: print("Done %s"%tfid_jaspar_pickle)


if __name__ == "__main__":
    main()


