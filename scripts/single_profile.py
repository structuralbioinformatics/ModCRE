import os, sys, re
from collections import Counter
import ConfigParser
import itertools
import numpy as np
import optparse
import shutil
import random
import threader
import cPickle
import pandas as pd

import matplotlib as mpl 
# agg backend is used to create plot as a .png file
mpl.use('Agg')
import matplotlib.pyplot as plt 

# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.dirname(__file__))

# Append scripts path to python path #
sys.path.append(scripts_path)

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Import my functions #
import functions

# Import my modules #
import profile as PROFILE


#-------------#
# Options     #
#-------------#

def parse_options():
    '''
    This function parses the command line arguments and returns an optparse object.
    '''

    parser = optparse.OptionParser("Usage: single_profile.py [--dummy=DUMMY_DIR] -i INPUT_FILE -d DNA_FASTA [-l LABEL -o OUTPUT_NAME --output_dir OUTPUT_DIR --info INFO_FILE] --pbm=PBM_dir --pdb=PDB_DIR [-v --save --plot --meme --reset] [--model_accuracy]  [-a -f -p -s SPLIT_POTENTIAL -e ENERGY_PROFILE -t THRESHOLD -k -b --taylor --file POTENTIAL --radius RADIUS --fragment FRAGMENT]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="DUMMY_DIR")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input PDB/Thread file. Mandatory", metavar="INPUT_FILE")
    parser.add_option("-l", "--label", action="store", type="string", default=None, dest="label", help="Label to organize the output files as 'label.energies.txt' ", metavar="LABEL")
    parser.add_option("--meme", default=False, action="store_true", dest="meme", help="Use 'uniprobe2meme' to calculate the PWM matrix for 'FIMO' (default = False)")
    parser.add_option("-r","--reset", default=False, action="store_true", dest="reset", help="Clean the sequences of the original MSA and reset them by a random selection in accordance with the PWM (default = False)")
    parser.add_option("-o", "--output", default=None, action="store", type="string", dest="output_name", help="Output name for tables and compressed PICKLE format (default = as PDB input)", metavar="OUTPUT_NAME")
    parser.add_option("--output_dir", default=None, action="store", type="string", dest="output_dir", help="Output directory (default = as PDB input)", metavar="OUTPUT_DIR")
    parser.add_option("--pbm", action="store", type="string", default=None, dest="pbm_dir", help="PBM directory (i.e. output dir from pbm.py). This is ,mandatory unless using option --file on potentials", metavar="PBM_DIR")
    parser.add_option("--pdb", action="store", type="string", default=None, dest="pdb_dir", help="PDB directory (i.e. output dir from pdb.py). This is ,mandatory unless using option --template", metavar="PDB_DIR")
    parser.add_option("-d","--dna",default=None, action="store", type="string", dest="dna_file", help="File of a DNA sequence in FASTA format to profile", metavar="FASTA")
    parser.add_option("--fragment", default=None, action="store", type="string", dest="fragment", help="Fragment of protein to apply the potential. Format is 'a-b;c-d': two regions between residues a-b and c-d. (Default is None it applies to all amino-acids)")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")
    parser.add_option("--save", default=False, action="store_true", dest="save", help="Save PDB models and scores while scanning the DNA (default = False)")
    parser.add_option("--plot", default=False, action="store_true", dest="plot", help="Plot profiles (default = False)")
    parser.add_option("--info", default=None, action="store", type="string", dest="info_file", help="File to store information of SUCCESS/FAILURE of the run (default = standard output)")
    parser.add_option("--model_accuracy", default=False, action="store_true", dest="model_accuracy", help="Calculate the scores with 3D-models to obtain the highest accuracy on distances (default = False, it threads the sequence of each DNA fragment without recalculating distances and reduces the time in half)")
    parser.add_option("--threading", default=False, action="store_true", dest="threading", help="Use a threading file as input (default is False)")

    group = optparse.OptionGroup(parser, "Statistical potentials", "Select your statistical potentials of choice. By default it uses S3DCdd general potential derived from PDB (the simplest one). In \"--auto\" mode, the program uses S3DCdd family potentials derived from both PDB and/or PBM data and/or approached by Taylor as selected in Potentials configuration file. In case family potentials cannot be applied, the program uses general potentials derived from both PDB and PBM data and approached by Taylor. \"-a\" option overrides options \"-f\", \"-p\" and \"-t\".")
    group.add_option("-a", "--auto", default=False, action="store_true", dest="auto_mode", help="Automate the selection of statistical potentials (default = False)")
    group.add_option("-f", "--family", default=False, action="store_true", dest="family_potentials", help="Use family potentials (default = False)")
    group.add_option("-p", default=False, action="store_true", dest="pbm_potentials", help="Use potentials derived from both PBM + PDB data (default = False)")
    group.add_option("-s","--potential", default="s3dc_dd", action="store", type="string", dest="split_potential", help="Split-potential to be used for PWM and FIMO (3d, 3dc, s3dc, s3dc_dd, s3dc_di, pair; default = s3dc_dd)", metavar="{string}")
    group.add_option("-e","--energy", default="all", action="store", type="string", dest="energy_profile", help="Select a specific Split-potential to be used for the profile (all, 3d, 3dc, s3dc, s3dc_dd, s3dc_di, pair; default = all)", metavar="{string}")
    group.add_option("-t", action="store", type="float", default=None, dest="score_threshold", help="Threshold on the scaled score to consider positive k-mers (default = 0.95)", metavar="{float}")
    group.add_option("-k","--known", default=False, action="store_true", dest="known", help="The name is of a known PDB file, with format 'code_chain' (default = False)")
    group.add_option("--taylor", default=False, action="store_true", dest="taylor_approach", help="Approach PMF by Taylor (default = False)")
    group.add_option("-b", "--bins", default=False, action="store_true",  dest="bins", help="Computate the potentials: by bins (if selected) or accumulative (default).", metavar="{boolean}")
    group.add_option("--file",default=None, action="store", type="string", dest="potential_file", help="Use potentials from specific file (default = None)", metavar="{string}")
    group.add_option("--radius",default=0, action="store", type="float", dest="radius", help="Maximum contact distance to calculate interactions (default=0 implies the use of 'max_contact_distance' from configuration", metavar="{string}")
    group.add_option("-m", "--pmf", default=False, action="store_true", dest="pmf", help="Use of raw mean-force potentials with no Z-scoring (default = False)")
    group.add_option("--methylation", default=False, action="store_true", dest="methylation", help="Use methylated cytosine specificities of binding/non-binding (default = False)")


    parser.add_option_group(group)

    
    (options, args) = parser.parse_args()

    if options.input_file is None or (options.pbm_dir is None and options.potential_file is None) or options.dna_file is None or options.pdb_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")
    if not functions.fileExist(options.dna_file):
        parser.error("Missing DNA fasta file")
    if not re.search("^3d$|^3dc$|^local$|^pair$|^s3dc$|^s3dc_dd$|^s3dc_di$|^all$", options.split_potential):
        parser.error("incorrect value for -s argument: type option \"-h\" for help")
        raise Exception("SINGLE_PROFILE: Wrong input")
    return options


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

# Arguments & Options #
    options = parse_options()
    verbose = options.verbose
    threading = options.threading
    save=options.save
    plot=options.plot
    dummy_dir = options.dummy_dir
    if not dummy_dir.startswith("/"): dummy_dir = os.path.abspath(options.dummy_dir)
    if not os.path.exists(dummy_dir): os.makedirs(dummy_dir)
    pdb_dir = options.pdb_dir
    if pdb_dir is not None:
       if not pdb_dir.startswith("/"): pdb_dir = os.path.abspath(options.pdb_dir)
    pbm_dir = options.pbm_dir
    if pbm_dir is not None:
       if not pbm_dir.startswith("/"): pbm_dir = os.path.abspath(options.pbm_dir)
    input_file = options.input_file
    if not input_file.startswith("/"): input_file = os.path.abspath(input_file)
    output_dir=os.path.dirname(input_file)
    if options.output_dir is not None: output_dir = options.output_dir
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    if threading:
       thread_file = input_file
       output_name=os.path.basename(thread_file)
    else:
       pdb_file = input_file
       output_name=os.path.basename(pdb_file)
    if output_name.endswith(".txt"): 
       output_name=output_name.rstrip(".txt")
    elif output_name.endswith(".thread"): 
       output_name=output_name.rstrip(".thread")
    elif output_name.endswith(".pdb"): 
       output_name=output_name.rstrip(".pdb")
    if options.output_name is not None: output_name = options.output_name
    if options.label is not None: output_name = output_name +"."+ options.label
    output_file = os.path.join(output_dir,output_name)
    if options.info_file is not None:
       info  = open(os.path.join(output_dir,options.info_file),"a")
       if verbose: print("Information of the run is stored in file %s"%options.info_file)
    else:
       if verbose: print("Information of the run is stored in standard output")
       info  = sys.stdout
    dna_file=options.dna_file
    if dna_file is not None:
      if not dna_file.startswith("/"): dna_file=os.path.abspath(dna_file)
    families = {}
    if options.verbose:sys.stdout.write("Check families...\n")
    if not os.path.exists(os.path.join(pdb_dir, "families.txt")):
       sys.stdout.write("Families file %s is not used\n"%(os.path.join(pdb_dir, "families.txt")))
    else:
      for line in functions.parse_file(os.path.join(options.pdb_dir, "families.txt")):
        if line.startswith("#"): continue
        pdb_chain, family = line.split(";")
        families[pdb_chain] = family
    # Potential input parameters
    potential_file    =options.potential_file
    radius            =options.radius
    split_potential   =options.split_potential
    energy_profile    =options.energy_profile
    auto_mode         =options.auto_mode
    family_potentials =options.family_potentials
    pbm_potentials    =options.pbm_potentials
    score_threshold   =options.score_threshold
    taylor_approach   =options.taylor_approach
    pmf               =options.pmf
    bins              =options.bins
    known             =options.known
    meme              =options.meme
    reset             =options.reset
    methylation       =options.methylation

    # Restrictions on TF fragment
    fragment_restrict=None
    if options.fragment is not None:
       fragment_restrict={}
       segment_fragment=options.fragment.split(";")
       for interval in segment_fragment:
           ac,bc = interval.split("-")
           if len(ac.split("_"))>0: a,chain1=ac.split("_")
           if len(bc.split("_"))>0: b,chain2=bc.split("_")
           if chain1==chain2:
              fragment_restrict.setdefault(chain1,[]).append((int(a),int(b)))
    # No restrictions on binding are allowed
    binding_restrict=None
    # Use accurated models
    model_accuracy=options.model_accuracy
    # Define fimo_thresholds
    fimo_thresholds=config.get("Parameters", "fimo_profile_thresholds").rstrip().split(",")

    try:

      if threading:
            profile=PROFILE.calculate_single_profile_of_thread(dna_file,fimo_thresholds,energy_profile,thread_file,output_file,pbm_dir,pdb_dir,families,potential_file, radius,fragment_restrict, binding_restrict, split_potential,auto_mode,family_potentials,pbm_potentials,score_threshold,taylor_approach,pmf,bins,known,meme,reset, dummy_dir,verbose,save,methylation)
      else:
         if model_accuracy:
            profile=PROFILE.calculate_single_profile_by_models(dna_file,fimo_thresholds,energy_profile,pdb_file,output_file,pbm_dir,pdb_dir,families,potential_file, radius,fragment_restrict, binding_restrict, split_potential,auto_mode,family_potentials,pbm_potentials,score_threshold,taylor_approach,pmf,bins,known,meme,reset, dummy_dir,verbose,save,methylation)
         else:
            profile=PROFILE.calculate_single_profile_by_thread(dna_file,fimo_thresholds,energy_profile,pdb_file,output_file,pbm_dir,pdb_dir,families,potential_file, radius,fragment_restrict, binding_restrict, split_potential,auto_mode,family_potentials,pbm_potentials,score_threshold,taylor_approach,pmf,bins,known,meme,reset, dummy_dir,verbose,save,methylation)

# Write outputs
      # Name the files
      output = output_file
      # Write profile in pickle format
      if functions.fileExist(output+".pickle"): os.remove(output+".pickle")
      out=open(output+".pickle","wb")
      cPickle.dump(profile,out)
      out.close()
      # Write table in CSV format
      if functions.fileExist(output+".csv"): os.remove(output+".csv")
      profile.write_table(output+".csv")
      # Plot profiles
      if plot: profile.plot(output)
      # Write information of success
      info.write("%s\tDONE\n"%output_name)
      info.flush()
   
    except Exception as e:

      info.write("%s\tFAILED\n"%output_name)
      info.flush()
      sys.stderr.write("ERROR %s"%e)




      



