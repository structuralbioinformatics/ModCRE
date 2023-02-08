import os, sys, re
import ConfigParser
import optparse
import string
import argparse
import cPickle
import numpy as np
import subprocess
from scipy import stats
from collections import Counter

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
import pwm_pbm as PWM
import tomtom as TOMTOM

#-------------#
# Options     #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser(
		epilog      = '@Oliva\'s lab 2018', 
		description = "The program compares a set of PWMs (models) with a database of PWMs (motifs). All options (except for --dummy and -v) are mandatory.")

    parser.add_option("--dummy", default="./tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="File of correspondences between FastA sequences of TFs and their PWMs in CisBP", metavar="{filename}")
    parser.add_option("-m", action="store", type="string", dest="input_models", help="File of correspondences between TF names and their modelled PWMs ", metavar="{filename}")
    parser.add_option("--query", action="store", type="string", dest="query", help="Name of the TF to analyse", metavar="{string}")
    parser.add_option("--motif_database", action="store", type="string", dest="pwm_db", help="File of the database of PWM's motifs ", metavar="{filename}")
    parser.add_option("--table", action="store", type="string", dest="table_output", help="Name of the output file table with the summary of TOTOM comparisons", metavar="{filename}")
    parser.add_option("--tom", action="store", type="string", dest="totom_output", help="Name of the cPickle-compressed file of all TOTOM comparisons", metavar="{filename}")
    parser.add_option("--pbm", action="store", type="string", dest="pbm_dir", help="PBM directory (i.e. output dir from pbm.py)", metavar="{directory}")
    parser.add_option("--pwm_models", action="store", type="string", dest="models_dir", default=None, help="PWM directory of modelled motifs", metavar="{directory}")
    parser.add_option("-t", action="store", type="string", dest="tfs_file", help="TFs file (from CIS-BP; i.e. cisbp_1.02.tfs.sql)", metavar="{filename}")
    parser.add_option("-f", action="store", type="string", dest="family_file", help="TFs Family file (from CIS-BP; i.e. cisbp_1.02.tf_families.sql)", metavar="{filename}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")

   
    (options, args) = parser.parse_args()
    if options.input_file is None  or options.pbm_dir is None  or options.models_dir is None or options.family_file is None or options.tfs_file is None or options.input_models is None or options.pwm_db is None or options.query is None or options.totom_output is None or options.table_output is None:
        parser.error("\t Missing arguments, type option \"-h\" for help")

    return options


def create_tf_models(omdltab,omdltmt,pwm_db,tf_families,tf_families_cisbp,tf_name_query,tf_nameseq,tf_nameseq_rev,tf_has_motif,tf_has_model,models_dir,verbose,tmp):
      tf_models={}
      tf_mdl_tomtom={}
      if not tf_has_model.has_key(tf_name_query):
       if  tf_nameseq.has_key(tf_name_query):
         tf_dummy=tf_nameseq.get(tf_name_query)
         if tf_has_model.has_key(tf_nameseq_rev.get(tf_dummy)[1]):
           tf_name_query_use=tf_nameseq_rev.get(tf_nameseq.get(tf_name_query))[1]
           if verbose: sys.stdout.write("Use model of %s\n"%tf_name_query_use)
           model_to_use=tf_has_model.get(tf_name_query_use)
         else:
           if verbose: sys.stdout.write("No model is defined for %s\n"%tf_name_query)
           return  (tf_models,tf_mdl_tomtom)
       else:
         if verbose: sys.stdout.write("No model is defined for %s\n"%tf_name_query)
         return  (tf_models,tf_mdl_tomtom)
      else:
         model_to_use=tf_has_model.get(tf_name_query)
      
      outab=open(omdltab,"w")
      outab.write("#%14s\t%45s\t%15s\t%45s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n"%("tf_query", "tf_model","motif_query","motif_model","p_value","e_value","q_value","p_score","e_score","q_score","Ranking","Family","CisBP_Fam"))
      for model_pair in model_to_use:
        if not tf_nameseq.has_key( tf_name_query ): continue
        tf_query=tf_nameseq.get(tf_name_query)
        tf_name_family=tf_nameseq_rev.get(tf_query)[1].upper()
        if not tf_has_motif.has_key(tf_query):continue
        motif_query = tf_has_motif.get(tf_query)
        tf_model    = model_pair[0]
        motif_model = model_pair[1]
        pwm_model   = os.path.join(models_dir,motif_model)
        try:
           tomtom_obj  = TOMTOM.get_tomtom_obj(pwm_db, pwm_model, tmp)
        except Exception as e:
           sys.stderr.write("TOMTOM HIT %s %s Error %s\n"%(pwm_db,pwm_model,e))
           continue
        tf_mdl_tomtom.setdefault(tf_model,tomtom_obj)
        if tomtom_obj.has_hit(motif_query):
           tomhit      = tomtom_obj.get_hit(motif_query)
           p_score = 20.0
           e_score = 20.0
           q_score = 20.0
           try:
              if float(tomhit.get_p_value()) > 1.0e-20 : p_score = -np.log(float(tomhit.get_p_value()))
              if float(tomhit.get_e_value()) > 1.0e-20 : e_score = -np.log(float(tomhit.get_e_value()))
              if float(tomhit.get_q_value()) > 1.0e-20 : q_score = -np.log(float(tomhit.get_q_value()))
           except Exception as e:
              sys.stderr.write("Skip TOMTOM HIT %s %s Error %s\n"%(tf_query,motif_query,e))
           rank = tomhit.get_rank()
           ranking = float(tomtom_obj.get_size()-rank +1)
           tf_models.setdefault((tf_name_query,tf_model),(p_score,e_score,q_score,ranking,0.0,0.0))
           tf_query_family="Unknown"
           tf_query_family_cisbp="Unknown"
           if tf_families.has_key(tf_name_family): tf_query_family="+".join(["_".join(str(x.strip()).split()) for x in tf_families.get(tf_name_family)])
           if tf_families_cisbp.has_key(tf_name_family): tf_query_family_cisbp="+".join(["_".join(str(x.strip()).split()) for x in tf_families_cisbp.get(tf_name_family)])
           if verbose: sys.stdout.write("TF_query: %s (motif %s) TF_MODEL: %s  (motif %s)  p_value: %e e_value: %e q_value %e p_score: %f e_score: %f q_score: %f  Ranking: %f Family: %s Family-code in CisBP: %s \n"%(tf_name_query,motif_query, tf_model,motif_model,tomhit.get_p_value(),tomhit.get_e_value(),tomhit.get_q_value(),p_score,e_score,q_score,ranking,tf_query_family,tf_query_family_cisbp))
           outab.write("%15s\t%45s\t%15s\t%45s\t%10.1e\t%10.1e\t%10.1e\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%10s\t%10s\n"%(tf_name_query, tf_model,motif_query,motif_model,tomhit.get_p_value(),tomhit.get_e_value(),tomhit.get_q_value(),p_score,e_score,q_score,ranking,tf_query_family,tf_query_family_cisbp))
      outab.close()
      outmt=open(omdltmt,"wb")
      cPickle.dump(tf_mdl_tomtom,outmt)
      outmt.close()
      return (tf_models,tf_mdl_tomtom)


def parse_tfs_motifs(input_file,pbm_dir):
    tf_has_motif={}
    tf_nameseq={}
    tf_nameseq_rev={}
    fd=open(input_file,"r")
    for line in fd:
        if line.startswith("#"): continue
        tf,motif=line.split()
        tf_has_motif.setdefault(tf,motif)
        tf_file  = pbm_dir+"/sequences/"+tf
        if not fileExist(tf_file): continue
        fo=open(tf_file,"r")
        for line in fo:
            if line.startswith(">"): 
                #tf_nameseq.setdefault(line.lstrip(">").split()[0],tf)
                #tf_nameseq_rev.setdefault(tf,line.lstrip(">").split()[0])
                tf_nameseq.setdefault(tf.rstrip(".fa"),tf)
                tf_nameseq_rev.setdefault(tf,(tf.rstrip(".fa"),line.lstrip(">").split()[0]))
        fo.close()
    fd.close()
    return tf_has_motif,tf_nameseq,tf_nameseq_rev

def parse_tfs_models(input_models):

    tf_has_model={}
    fd=open(input_models,"r")
    for line in fd:
        if line.startswith("#"): continue
        #tf, pdb, motif_a, motif_b=line.split()
        tf, motif_a=line.split()
        if motif_a.endswith(".meme.s"): pdb=motif_a.rstrip(".meme.s")+".pdb"
        if motif_a.endswith(".meme"): pdb=motif_a.rstrip(".meme")+".pdb"
        tf_has_model.setdefault(tf.rstrip(".txt"),[]).append((pdb,motif_a))
        #if tf_nameseq.has_key(tf):
           #tf_has_model.setdefault(tf_nameseq.get(tf),(pdb,motif_a))
    fd.close()

    return tf_has_model
      


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
    verbose = options.verbose
    tf_name_query = options.query
    omdltab = options.table_output
    omdltmt = options.totom_output
    pwm_db  = options.pwm_db
    models_dir = options.models_dir
    # Initialize #
    if not os.path.exists(options.dummy_dir): os.makedirs(options.dummy_dir)
    tmp     = options.dummy_dir
    #Read CisBP families
    if verbose: sys.stdout.write("-- Read CisBP families\n")
    cisbp_families={}
    for line in functions.parse_file(os.path.abspath(options.family_file)):
        m = re.search("\('(.+)', '(.+)', '.+', .+, .+\)", line)
        if m:           
            cisbp_families.setdefault(m.group(1).upper(),set()).add(m.group(2))
    #Read CisBP TFs
    tf_species  = {}
    tf_families = {}
    tf_families_cisbp = {}
    for line in functions.parse_file(os.path.abspath(options.tfs_file)):
         m = re.search("\('(.+)', '(.+)', '.+', '.+', '.+', '(.+)', '[DIN]'\),*", line)
         if m:
             tf_species.setdefault(m.group(1), set()).add(m.group(3).replace("_", " ").upper())
             tf_families_cisbp.setdefault(m.group(1).upper(), set()).add(m.group(2).upper())
             tf_families.setdefault(m.group(1).upper(), set()).update(cisbp_families[m.group(2).upper()])

    #Create TF-Motif correspondence   
    if verbose: sys.stdout.write("-- Parsing input of TFs and their motifs\n")
    tf_has_motif,tf_nameseq,tf_nameseq_rev= parse_tfs_motifs(options.input_file,options.pbm_dir)
 
    #Create TF-model correspondence   
    if verbose: sys.stdout.write("-- Parsing input of TFs and their models\n")
    tf_has_model=parse_tfs_models(options.input_models)

#Run the TOTOM comparison of models with the database of motifs

    tf_models_tf,tf_mdl_tomtom=create_tf_models(omdltab,omdltmt,pwm_db,tf_families,tf_families_cisbp,tf_name_query,tf_nameseq,tf_nameseq_rev,tf_has_motif,tf_has_model,models_dir,verbose,tmp)

#Finish

    sys.stdout.write("Done. \n")
    exit(0)


        
if __name__ == "__main__":
    main()

