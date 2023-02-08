import os, sys, re
import ConfigParser
import optparse
import subprocess
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import string
import argparse
from scipy import stats
import cPickle
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


def fileExist(file):
    '''
    Check existing files
    '''
    if file is not None:
        return os.path.exists(file) and os.path.isfile(file)
    else:
        return False

def clean_name(name):
    new_name=""
    for x in name:
        if x == "/" or x == "," or x == ";" or x == " ": new_name= new_name + "_"
        else: new_name = new_name + x
    return new_name

#-------------#
# Options     #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python nearest_neighbours.py -i tf_motifs_file -m tf_models_file --pwm=pwm_dir --pdb=pbm_dir --pwm_models=models_dir [-o output_dir]",
                                   epilog      = '@Oliva\'s lab 2018',
                                   description = "The program calculates the closest similar sequences (nearest neighbour) of each TF and compares their PWMs using TOMTOM. It also compares the PWMs modelled for each TF  with the dataset of PWMs. Boxplots are built to compare both success using different conditions and scores derived from TOMTOM results")

    parser.add_option("--dummy", default="./tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="File of correspondences between FastA sequences of TFs and their PWMs in CisBP", metavar="{filename}")
    parser.add_option("-m", action="store", type="string", dest="input_models", help="File of correspondences between TF names and their modelled PWMs ", metavar="{filename}")
    parser.add_option("-o", action="store", default="output_neighbours", type="string", dest="output_file", help="Root name for output files and directories (default = 'output_neighbopurs')", metavar="{rootname}")
    parser.add_option("-d", action="store", default=10, type="int", dest="plot_step", help="Interval size of ID percentage to group nieghbours (default = 10 )", metavar="{integer}")
    parser.add_option("-s", "--selected_score",action="store", default=3, type="int", dest="selected_score", help="Select the TOMTOM score to do the plots ( 0 is -log(P value); 1 is -log(E value); 2 is -log(Q value); 3 is ranking;  4 is enrichment and forces enrichment flag; 5 is ranking of enrichment and forces enrichment flag; default is 3 )", metavar="{integer}")
    parser.add_option("--top", default=1, action="store", type="int", dest="top_threshold", help="Threshold on the number of top motifs to calculate the ratio of True-Positives (only applies when option -s 3 is set, default is 1)", metavar="{integer}")
    parser.add_option("--etop", default=1, action="store", type="int", dest="top_enrichment", help="Threshold on the number of TOP motifs to calculate the enrichment (default is 1)", metavar="{integer}")
    parser.add_option("--pbm", action="store", type="string", dest="pbm_dir", help="PBM directory (i.e. output dir from pbm.py)", metavar="{directory}")
    parser.add_option("--pwm", action="store", type="string", dest="pwm_dir", default=None, help="PWM directory with PWM motifs from CisBP", metavar="{directory}")
    parser.add_option("--pwm_models", action="store", type="string", dest="models_dir", default=None, help="PWM directory of modelled motifs", metavar="{directory}")
    parser.add_option("-t", action="store", type="string", dest="tfs_file", help="TFs file (from CIS-BP; i.e. cisbp_1.02.tfs.sql)", metavar="{filename}")
    parser.add_option("-f", action="store", type="string", dest="family_file", help="TFs Family file (from CIS-BP; i.e. cisbp_1.02.tf_families.sql)", metavar="{filename}")
    parser.add_option("--use_CisBP_family", default=False, action="store_true", dest="cisbp_family",help="Use family codes from CisBP instead of common family names of TFs (default = False)", metavar="{boolean}")
    parser.add_option("-n", "--normalize", default=False, action="store_true", dest="normalize", help="Normalize TOMTOM score from 0 to 100 (default = False) ")
    parser.add_option("-e", "--enrichment", default=False, action="store_true", dest="enrichment", help="Flag to calculate  and use enrichment (Options s=4 and s=5 on score force automatically this flag, default = False) ")
    parser.add_option("--average", default=False, action="store_true", dest="average", help="Use the average scores obtained with all models of a TF (default = False it uses all models) ")
    parser.add_option("--maximum", default=False, action="store_true", dest="maximum", help="Use the maximum score among all models of a TF (default = False it uses all models) ")
    parser.add_option("--minimum", default=False, action="store_true", dest="minimum", help="Use the minimum score among all models of a TF (default = False it uses all models) ")
    parser.add_option("--skip", default=False,action="store_true", dest="skip_compare",help="Skip the comparison between TFs and stop after storing the TOTOM comparisons with models (default = False)")
    parser.add_option("-p","--parallel",default=False, action="store_true", dest="parallel", help="Run in parallel the comparisons between the database of motifs and all models of each TF (this automatically forces 'skip' flag, default is False)")
    parser.add_option("-v","--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")

   
    (options, args) = parser.parse_args()
    if (options.input_file is None ) or options.pbm_dir is None or options.pwm_dir is None or options.models_dir is None or options.family_file is None or options.tfs_file is None or options.input_models is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def parse_mmseq_comparisons(input_file):
    compair={}
    if fileExist(input_file):
       fi=open(input_file,"r")
       for line in fi:
           if line.startswith("#"):continue
           data=line.split()
           query = data[0]
           hit   = data[1]
           perc  = float(data[2])
           blast = float(data[10])
           triad  = (hit,perc,blast)
           compair.setdefault(query,set()).add(triad)
       fi.close()
       return compair
    else:
       sys.stderr.write("Input file with MMSEQ comparisons %s not found\n"%input_file)
       return None

def parse_table(table):
    tf_compare={}
    fo=open(table,"r")
    for line in fo:
      if line.startswith("#"): continue
      tf_query, tf_target,motif_query,motif_target,tf_id,tf_blast,p_value,e_value,q_value,p_score,e_score,q_score,ranking,tf_query_family,tf_query_family_cisbp = line.split()
      tf_compare.setdefault((tf_query,tf_target),(p_score,e_score,q_score,ranking,tf_id,tf_blast))
    fo.close()
    return tf_compare

def parse_table_models(table):
    tf_models={}
    fo=open(table,"r")
    for line in fo:
      if line.startswith("#"): continue
      tf_query, tf_model,motif_query,motif_model,p_value,e_value,q_value,p_score,e_score,q_score,ranking,tf_query_family,tf_query_family_cisbp = line.split()
      tf_models.setdefault((tf_query,tf_model),(p_score,e_score,q_score,ranking,0.0,0.0))
    fo.close()
    return tf_models

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
      
def normalizing(x,mx,mn):
    y=x
    if mx>mn:
       y = 100* (x - mn) / (mx - mn)
    if y<0:
       y = 0
    return y


def extract_TP_ratio(filename,suffix,size,top,step):
       tpr={} 
       for perc in xrange(100,step,-step):
         middle=int(perc-step/2)
         tp = 0
         dim= 0
         if fileExist(filename+str(middle)+suffix):
           fo=open(filename+str(middle)+suffix,"r")
           for line in fo:
             dim = dim + 1
             dummy,value=line.split()
             rank = size - float(value) + 1
             if rank <= top: tp = tp + 1
           fo.close()
           tpr.setdefault(int(middle), float(tp) / float(dim))
       return tpr

def plotcurve(p,a,b,out,title):
  x = np.array(p)
  y = np.array(a)
  z = np.array(b)
  plt.plot(x,y,'-o', ms=5, lw=2, alpha=0.7, mfc='cyan')
  plt.plot(x,z,'-o', ms=5, lw=2, alpha=0.7, mfc='firebrick')
  plt.xlabel('percentage of identity')
  plt.ylabel('ratio of true positives')
  plt.title("%s"%title)
  plt.grid(True)
  plt.savefig(out)
  plt.close()

def printverbose(f,flag,message):
 """Define verbose to print in file 'f', if 'flag'=True, message given as 'message'"""
 if flag: f.write("%s"%(message))

def create_inputboxplot(input_name,label,plot_step,fileroot,selected_score,verbose):
    file_root=fileroot+"_"+label+"_"
    if verbose: sys.stdout.write("Global PLOT, input %s \n"%(input_name))
    fo=open(input_name,"w")
    for perc in xrange(100,plot_step,-plot_step):
        middle=int(perc-plot_step/2)
        fo.write("%s\t%s\n"%(str(perc),file_root+str(middle)+".dat"))
        fo.write("%s\t%s\n"%(str(perc)+"M",file_root+str(middle)+"M.dat"))
        if fileExist(file_root+str(middle)+".dat"): os.remove(file_root+str(middle)+".dat")
        if fileExist(file_root+str(middle)+"M.dat"): os.remove(file_root+str(middle)+"M.dat")
    fo.close()

def create_seqboxplot_global_enrichment(tf,plot_step,fileroot,selected_score,normalize,maximum,minimum,verbose):
    found_interval={}
    label = "global"
    for tf_enrich, scores in tf.iteritems():
        tf_name,tf_id=tf_enrich
        if verbose: sys.stdout.write("Write files for global PLOT %s\n"%(fileroot+"_"+label+"_"+str(tf_id)+".dat"))
        scr = float(scores[selected_score])
        if normalize:
           scr = normalizing(float(scores[selected_score]),maximum,minimum)
        fo=open(fileroot+"_"+label+"_"+str(tf_id)+".dat","a+")
        fo.write("%s\t%f\n"%(tf_name,scr))
        found_interval.setdefault(tf_id,set()).add(tf_name)
        fo.close()
    return found_interval

 
def create_seqboxplot_global(tf,plot_step,fileroot,selected_score,normalize,maximum,minimum,verbose):
    if selected_score > 3: selected_score=3
    found_interval={}
    label = "global"
    for tf_pair, scores in tf.iteritems():
        tf_id=100*float(scores[4])
        for perc in xrange(plot_step,100,plot_step):
            middle=int(perc+plot_step/2)
            if tf_id > perc and tf_id <= perc+plot_step:
                if verbose: sys.stdout.write("Write files for global PLOT %s\n"%(fileroot+"_"+label+"_"+str(middle)+".dat"))
                scr = float(scores[selected_score])
                if normalize:
                   scr = normalizing(float(scores[selected_score]),maximum,minimum)
                fo=open(fileroot+"_"+label+"_"+str(middle)+".dat","a+")
                fo.write("%s\t%f\n"%("+".join(tf_pair),scr))
                found_interval.setdefault(middle,set()).add(tf_pair[0])
                fo.close()
    return found_interval

def create_mdlboxplot_global(tf,found_interval,plot_step,fileroot,selected_score,normalize,maximum,minimum,verbose):
    for tf_pair, scores in tf.iteritems():
        for perc in xrange(plot_step,100,plot_step):
            middle=int(perc+plot_step/2)
            if not  found_interval.has_key(middle): continue
            if tf_pair[0] in found_interval.get(middle):
              if verbose: sys.stdout.write("Write files for global PLOT %s\n"%(fileroot+"_global_"+str(middle)+"M.dat"))
              scr = float(scores[selected_score])
              if normalize:
                   scr = normalizing(float(scores[selected_score]),maximum,minimum)
              fo=open(fileroot+"_global_"+str(middle)+"M.dat","a+")
              fo.write("%s\t%f\n"%("+".join(tf_pair),scr))
              fo.close()

def create_seqboxplot_family_enrichment(tf,tf_nameseq,tf_nameseq_rev,families,plot_step,fileroot,selected_score,normalize,maximum,minimum,verbose):
    for tf_enrich, scores in tf.iteritems():
        tf_name,tf_id=tf_enrich
        tf_query=tf_name.upper() 
        tf_query_use=get_tf_family(tf_query,tf_nameseq,tf_nameseq_rev)
        if families.has_key(tf_query_use):
          fam=families.get(tf_query_use)
          family=clean_name("+".join(["_".join(str(x).split()) for x in fam]))
          scr = float(scores[selected_score])
          if normalize:
                      scr = normalizing(float(scores[selected_score]),maximum.get(family),minimum.get(family))
          if verbose: sys.stdout.write("Write files for family %s  PLOT %s\n"%(family,fileroot+"_"+family+"_"+str(tf_id)+".dat"))
          fo=open(fileroot+"_"+family+"_"+str(tf_id)+".dat","a+")
          fo.write("%s\t%f\n"%(tf_name,scr))
          fo.close()

def create_seqboxplot_family(tf,tf_nameseq,tf_nameseq_rev,families,plot_step,fileroot,selected_score,normalize,maximum,minimum,verbose):
    if selected_score > 3: selected_score=3
    for tf_pair, scores in tf.iteritems():
        tf_query=tf_pair[0].upper()
        tf_query_use=get_tf_family(tf_query,tf_nameseq,tf_nameseq_rev)
        tf_id=100*float(scores[4])
        if families.has_key(tf_query_use):
          fam=families.get(tf_query_use)
          family=clean_name("+".join(["_".join(str(x).split()) for x in fam]))
          for perc in xrange(plot_step,100,plot_step):
            middle=int(perc+plot_step/2)
            if tf_id > perc and tf_id <= perc+plot_step:
                scr = float(scores[selected_score])
                if normalize:
                      scr = normalizing(float(scores[selected_score]),maximum.get(family),minimum.get(family))
                if verbose: sys.stdout.write("Write files for family %s  PLOT %s\n"%(family,fileroot+"_"+family+"_"+str(middle)+".dat"))
                fo=open(fileroot+"_"+family+"_"+str(middle)+".dat","a+")
                fo.write("%s\t%f\n"%("+".join(tf_pair),scr))
                fo.close()

def create_mdlboxplot_family(tf,tf_nameseq,tf_nameseq_rev,families,found_interval,plot_step,fileroot,selected_score,normalize,maximum,minimum,verbose):
    for tf_pair, scores in tf.iteritems():
        tf_query=tf_pair[0].upper()
        tf_query_use=get_tf_family(tf_query,tf_nameseq,tf_nameseq_rev)
        if families.has_key(tf_query_use):
          fam=families.get(tf_query_use)
          family=clean_name("+".join(["_".join(str(x).split()) for x in fam]))
          for perc in xrange(plot_step,100,plot_step):
            middle=int(perc+plot_step/2)
            if fileExist(fileroot+"_"+family+"_"+str(middle)+".dat") and tf_pair[0] in found_interval.get(middle):
               scr = float(scores[selected_score])
               if normalize:
                      scr = normalizing(float(scores[selected_score]),maximum.get(family),minimum.get(family))
               if verbose: sys.stdout.write("Write files for family %s  PLOT %s\n"%(family,fileroot+"_"+family+"_"+str(middle)+"M.dat"))
               fo=open(fileroot+"_"+family+"_"+str(middle)+"M.dat","a+")
               fo.write("%s\t%f\n"%("+".join(tf_pair),scr))
               fo.close()

def modelExist(tf_name_query,tf_nameseq,tf_nameseq_rev,tf_has_motif,tf_has_model,verbose):
      if not tf_has_model.has_key(tf_name_query):
        if  tf_nameseq.has_key(tf_name_query):
          tf_dummy=tf_nameseq.get(tf_name_query)
          if tf_has_model.has_key(tf_nameseq_rev.get(tf_dummy)[1]):
            if verbose: sys.stdout.write("Use model of %s\n"%tf_name_query_use)
            return True
          else:
            if verbose: sys.stdout.write("No model is defined for %s\n"%tf_name_query)
            return False
        else:
         if verbose: sys.stdout.write("No model is defined for %s\n"%tf_name_query)
         return False  
      else:
        return True

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
        tomtom_obj  = TOMTOM.get_tomtom_obj(pwm_db, pwm_model, tmp)
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

def create_average_models(tf_models):
    tf_list_mdl_scores={}
    tf_average_models={}
    tf_maxim_models={}
    tf_minim_models={}
    for tf_pair,scores in tf_models.iteritems():
            tf_list_mdl_scores.setdefault(tf_pair[0],[]).append(scores)
    for tf_name_query, list_scores in tf_list_mdl_scores.iteritems():
      if len(list_scores) > 0:
        p_score = sum([float(x[0]) for x in list_scores])/len(list_scores)
        e_score = sum([float(x[1]) for x in list_scores])/len(list_scores)
        q_score = sum([float(x[2]) for x in list_scores])/len(list_scores)
        ranking = sum([float(x[3]) for x in list_scores])/len(list_scores)
        enrich  = sum([float(x[4]) for x in list_scores])/len(list_scores)
        erank   = sum([float(x[5]) for x in list_scores])/len(list_scores)
        tf_average_models.setdefault((tf_name_query,"average"),(p_score,e_score,q_score,ranking,enrich,erank))
        p_score = max([float(x[0]) for x in list_scores])
        e_score = max([float(x[1]) for x in list_scores])
        q_score = max([float(x[2]) for x in list_scores])
        ranking = max([float(x[3]) for x in list_scores])
        enrich  = max([float(x[4]) for x in list_scores])
        erank   = max([float(x[5]) for x in list_scores])
        tf_maxim_models.setdefault((tf_name_query,"maximum"),(p_score,e_score,q_score,ranking,enrich,erank))
        p_score = min([float(x[0]) for x in list_scores])
        e_score = min([float(x[1]) for x in list_scores])
        q_score = min([float(x[2]) for x in list_scores])
        ranking = min([float(x[3]) for x in list_scores])
        enrich  = min([float(x[4]) for x in list_scores])
        erank   = min([float(x[5]) for x in list_scores])
        tf_minim_models.setdefault((tf_name_query,"minimum"),(p_score,e_score,q_score,ranking,enrich,erank))

    return (tf_average_models,tf_maxim_models,tf_minim_models)

def create_tf_comparisons(outable,dmdltmt,pwm_out,pwm_db,tf_neighbours,tf_nameseq,tf_nameseq_rev,tf_has_motif,tf_families,tf_families_cisbp,verbose,tmp):
     if verbose:sys.stdout.write("Create table of comparisons %s \n"%outable)
     outab=open(outable,"w")
     outab.write("#%14s\t%15s\t%15s\t%15s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n"%("tf_query", "tf_target","motif_query","motif_target","%id","BLAST_e_val","p_value","e_value","q_value","p_score","e_score","q_score","Ranking","Family","CisBP_Fam"))
     tf_tomtom={}
     #Store TOMTOM comparisons of ALL TFs
     for tf_name_query,tf_relatives in tf_neighbours.iteritems(): 
        if not tf_nameseq.has_key( tf_name_query ): continue
        tf_query=tf_nameseq.get(tf_name_query)
        if not tf_has_motif.has_key(tf_query):continue
        motif_query = tf_has_motif.get(tf_query)
        pwm_query   = os.path.abspath(os.path.join(pwm_out,motif_query.rstrip(".txt")+".meme"))
        omdltmt = dmdltmt + "/"+motif_query.rstrip(".txt")+".pickle"
        if fileExist(omdltmt):
          if verbose: sys.stdout.write("Use DUMP of TOTOM MOTIF %s\n"%(omdltmt))
        else:
          if verbose: sys.stdout.write("Runing TOMTOM for TF_query: %s (%s vs %s)\n"%(tf_name_query,pwm_db,pwm_query))
          tomtom_obj  = TOMTOM.get_tomtom_obj(os.path.abspath(pwm_db), pwm_query, os.path.abspath(tmp))
          if verbose: sys.stdout.write("  -- Query: %s Rank: %d\n"%(tomtom_obj.get_query(),tomtom_obj.get_size()))
          outmt=open(omdltmt,"wb")
          cPickle.dump(tomtom_obj,outmt)
          outmt.close()
        tf_tomtom.setdefault(motif_query,omdltmt)
     tf_compare={}
     if verbose: sys.stdout.write("--Construct tables of comparison\n")
     for tf_name_query,tf_relatives in tf_neighbours.iteritems():  
      if not tf_nameseq.has_key( tf_name_query ): continue
      tf_query=tf_nameseq.get(tf_name_query)
      tf_name_query_family=tf_nameseq_rev.get(tf_query)[1].upper()
      if not tf_has_motif.has_key(tf_query):continue
      motif_query = tf_has_motif.get(tf_query)
      pwm_query   = pwm_out+"/"+motif_query.rstrip(".txt")+".meme"
      for tf_name_target,tf_id,tf_blast in tf_relatives:
          if not tf_nameseq.has_key(tf_name_target):continue
          tf_target=tf_nameseq.get(tf_name_target)
          tf_name_target_family=tf_nameseq_rev.get(tf_target)[1].upper()
          if not tf_has_motif.has_key(tf_target):continue
          motif_target= tf_has_motif.get(tf_target)
          pwm_target  = pwm_out+"/"+motif_target.rstrip(".txt")+".meme"
          if not tf_tomtom.has_key(motif_target): continue
          tomtom_obj_pickle  = open(tf_tomtom.get(motif_target),"rb")
          tomtom_obj         = cPickle.load(tomtom_obj_pickle)
          tomtom_obj_pickle.close()
          if tomtom_obj.has_hit(motif_query):
           tomhit      = tomtom_obj.get_hit(motif_query)
           p_score = -np.log(1.0e-20)
           e_score = -np.log(1.0e-20)
           q_score = -np.log(1.0e-20)
           try:
              if float(tomhit.get_p_value()) > 1.0e-20 : p_score = -np.log(float(tomhit.get_p_value()))
              if float(tomhit.get_e_value()) > 1.0e-20 : e_score = -np.log(float(tomhit.get_e_value()))
              if float(tomhit.get_q_value()) > 1.0e-20 : q_score = -np.log(float(tomhit.get_q_value()))
           except Exception as e:
              sys.stderr.write("Skip TOMTOM HIT %s %s Error %s\n"%(tf_query,tf_target,e))
           rank = tomhit.get_rank()
           ranking = float(tomtom_obj.get_size()-rank +1)
           tf_compare.setdefault((tf_name_query,tf_name_target),(p_score,e_score,q_score,ranking,tf_id,tf_blast))
           tf_query_family="Unknown"
           tf_query_family_cisbp="Unknown"
           if tf_families.has_key(tf_name_query_family): tf_query_family="+".join(["_".join(str(x.strip()).split()) for x in tf_families.get(tf_name_query_family)])
           if tf_families_cisbp.has_key(tf_name_query_family): tf_query_family_cisbp="+".join(["_".join(str(x.strip()).split()) for x in tf_families_cisbp.get(tf_name_query_family)])
           if verbose: sys.stdout.write("TF_query: %s (motif %s) TF_hit: %s  (motif %s) ID: %f BLAST e_value: %e p_value: %e e_value: %e q_value %e p_score: %f e_score: %f q_score: %f  Ranking: %f Family: %s Family-code in CisBP: %s \n"%(tf_name_query,motif_query, tf_name_target,motif_target,tf_id,tf_blast,tomhit.get_p_value(),tomhit.get_e_value(),tomhit.get_q_value(),p_score,e_score,q_score,ranking,tf_query_family,tf_query_family_cisbp))
           outab.write("%15s\t%15s\t%15s\t%15s\t%10.5f\t%10.1e\t%10.1e\t%10.1e\t%10.1e\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%10s\t%10s\n"%(tf_name_query, tf_name_target,motif_query,motif_target,tf_id,tf_blast,tomhit.get_p_value(),tomhit.get_e_value(),tomhit.get_q_value(),p_score,e_score,q_score,ranking,tf_query_family,tf_query_family_cisbp))
     outab.close()

     return tf_compare

def get_tf_family(tf_query,tf_nameseq,tf_nameseq_rev):
      tf_query_use=tf_query
      if  tf_nameseq.has_key(tf_query):
              tf_dummy=tf_nameseq.get(tf_query)
              if tf_nameseq_rev.has_key(tf_dummy):
                 tf_query_use=tf_nameseq_rev.get(tf_nameseq.get(tf_query))[1]
      return tf_query_use 
       
def calculate_maxmin(tf_compare,tf_models,tf_nameseq,tf_nameseq_rev,tf_families,set_of_families,selected_score):

      #if selected_score >3: 
      #    selected_score_seq = 3
      #else:
      #    selected_score_seq = selected_score

      selected_score_seq = selected_score

      seq_scores_in_family={}
      for tf_pair, scores in tf_compare.iteritems():    
          tf_query=tf_pair[0].upper() 
          seq_scores_in_family.setdefault("global",[]).append(float(scores[selected_score_seq]))
          tf_query_use=get_tf_family(tf_query,tf_nameseq,tf_nameseq_rev)
          if tf_families.has_key(tf_query_use): 
             fam=tf_families.get(tf_query_use)
             family=clean_name("+".join(["_".join(str(x).split()) for x in fam]))
             seq_scores_in_family.setdefault(family,[]).append(float(scores[selected_score_seq]))

      mdl_scores_in_family={}
      for tf_pair, scores in tf_models.iteritems(): 
          tf_query=tf_pair[0].upper()
          mdl_scores_in_family.setdefault("global",[]).append(float(scores[selected_score]))
          tf_query_use=get_tf_family(tf_query,tf_nameseq,tf_nameseq_rev)
          if tf_families.has_key(tf_query_use): 
             fam=tf_families.get(tf_query_use)
             family=clean_name("+".join(["_".join(str(x).split()) for x in fam]))
             mdl_scores_in_family.setdefault(family,[]).append(float(scores[selected_score]))

      maximum_seq={}
      if seq_scores_in_family.has_key("global"): maximum_seq.setdefault("global",max([float(x) for x in seq_scores_in_family.get("global")]))
      for family in set_of_families:
        if seq_scores_in_family.has_key(family):
          maximum_seq.setdefault(family,max([float(x) for x in seq_scores_in_family.get(family)]))

      minimum_seq={}
      if seq_scores_in_family.has_key("global"): minimum_seq.setdefault("global",min([float(x) for x in seq_scores_in_family.get("global")]))
      for family in set_of_families:
        if seq_scores_in_family.has_key(family):
          minimum_seq.setdefault(family,min([float(x) for x in seq_scores_in_family.get(family)]))

      maximum_mdl={}
      if mdl_scores_in_family.has_key("global"): maximum_mdl.setdefault("global",max([float(x) for x in mdl_scores_in_family.get("global")]))
      for family in set_of_families:
        if mdl_scores_in_family.has_key(family):
          maximum_mdl.setdefault(family,max([float(x) for x in mdl_scores_in_family.get(family)]))

      minimum_mdl={}
      if mdl_scores_in_family.has_key("global"): minimum_mdl.setdefault("global",min([float(x) for x in mdl_scores_in_family.get("global")]))
      for family in set_of_families:
        if mdl_scores_in_family.has_key(family):
          minimum_mdl.setdefault(family,min([float(x) for x in mdl_scores_in_family.get(family)]))

      return maximum_seq,minimum_seq,maximum_mdl,minimum_mdl


def calculate_enrichment_neighbours(outable,top,normalize,plot_step,number_of_motifs,verbose):

    rank_target={}
    set_targets={}
    tf_scores={}
    toprank = int(top)
    fo=open(outable,"r")
    for line in fo:
      if line.startswith("#"): continue
      tf_query, tf_target,motif_query,motif_target,tf_id_ratio,tf_blast,p_value,e_value,q_value,p_score,e_score,q_score,ranking,tf_query_family,tf_query_family_cisbp = line.split()
      tf_id=100*float(tf_id_ratio)
      rank_target.setdefault(tf_target,[]).append((tf_query,float(ranking)))
      tf_scores.setdefault((tf_query,tf_target),(p_score,e_score,q_score,ranking,tf_id_ratio,tf_blast))
      for perc in xrange(plot_step,100,plot_step):
            middle=int(perc+plot_step/2)
            if tf_id > perc and tf_id <= perc+plot_step:
               set_targets.setdefault((tf_query,middle),set()).add(tf_target)
    fo.close()
    tf_number_of_neighbours={}
    tf_list_of_hits_by_neighbour={}
    tf_pval_list_of_hits={}
    tf_eval_list_of_hits={}
    tf_qval_list_of_hits={}
    for tf_perc,set_tf in  set_targets.iteritems():
        tf,percentage_id=tf_perc
        tf_number_of_neighbours.setdefault(tf_perc,len(set_tf))
        for target in set_tf:
          if rank_target.has_key(target):
            for result in rank_target.get(target):
                query,ranking=result
                if ranking > number_of_motifs-toprank:
                   tf_list_of_hits_by_neighbour.setdefault(tf_perc,[]).append(query)
                   if query==tf:
                      tf_pval_list_of_hits.setdefault(tf_perc,[]).append(float(tf_scores.get((query,target))[0]))
                      tf_eval_list_of_hits.setdefault(tf_perc,[]).append(float(tf_scores.get((query,target))[1]))
                      tf_qval_list_of_hits.setdefault(tf_perc,[]).append(float(tf_scores.get((query,target))[2]))
         
    
    tf_enrichment={}
    for tf_perc, list_of_hits in tf_list_of_hits_by_neighbour.iteritems():
        tf_query,percentage_id=tf_perc
        group=Counter(list_of_hits)
        if not tf_number_of_neighbours.has_key(tf_perc): continue
        if group.has_key(tf_query) and tf_number_of_neighbours.get(tf_perc)>0:
           e=100*float(group.get(tf_query))/tf_number_of_neighbours.get(tf_perc)
        else:
           e=0.0
        count= 0
        n    = 0
        #Find the rank (top_rank) of the query on the top selected results
        rank = number_of_motifs 
        for hit,k in sorted(group.iteritems(),key=lambda x:x[1],reverse=True):
            if count != int(k): 
                count=int(k)
                n    =n+1
            if hit == tf_query: 
                rank = float(n)
        top_rank=int(number_of_motifs-rank+1)
        #Find the minimum enrichment with p-value lower than the random selection (1/number_of_motifs): lower than this implies the enrichment is as good as random
        if e>0 and group.has_key(tf_query) and tf_number_of_neighbours.get(tf_perc)>0:
          min_enrichment,worst_enrichment=best_enrichment(tf_number_of_neighbours.get(tf_perc),number_of_motifs,toprank)
          if verbose: sys.stdout.write("\t-- %s enrichment: %f (lower limit %f) number of neighbours %d at ID %d\n"%(tf_query,e,min_enrichment,tf_number_of_neighbours.get(tf_perc),percentage_id)) 
        else:
          min_enrichment=0
          worst_enrichment=0
        #Get the maximum number of hits
        max_enrichment=max([x for x in group.itervalues()])
        #Get the average of scores (p_value, e_value,q_value) with all the selected solutions in the enrichment set
        p_score=0.0
        e_score=0.0
        q_score=0.0
        if tf_pval_list_of_hits.has_key(tf_perc): p_score = sum([float(x) for x in tf_pval_list_of_hits.get(tf_perc)])/len(tf_pval_list_of_hits.get(tf_perc))
        if tf_eval_list_of_hits.has_key(tf_perc): e_score = sum([float(x) for x in tf_eval_list_of_hits.get(tf_perc)])/len(tf_eval_list_of_hits.get(tf_perc))
        if tf_qval_list_of_hits.has_key(tf_perc): q_score = sum([float(x) for x in tf_qval_list_of_hits.get(tf_perc)])/len(tf_qval_list_of_hits.get(tf_perc))
        #store the data in a new dictionary named tf_enrichment
        if e<min_enrichment:
          tf_enrichment.setdefault(tf_perc,(p_score,e_score,q_score,0,e,0))
        else:   
          tf_enrichment.setdefault(tf_perc,(p_score,e_score,q_score,float(top_rank),e,float(top_rank)))

    return tf_enrichment


    
 
def calculate_enrichment(tf_models,tf_mdl_tomtom,tf_has_motif,tf_nameseq,top,normalize,number_of_motifs,verbose):

    toprank = int(top)

    tf_number_of_models={}
    tf_list_of_models={}
    tf_size_tomtom={}

    for tf_pair,scores in tf_models.iteritems():
        tf_query,tf_model = tf_pair
        tf_list_of_models.setdefault(tf_query,[]).append(tf_model)

    for tf_query,list_of_models in tf_list_of_models.iteritems():
        tf_number_of_models.setdefault(tf_query,len(list_of_models))

    tf_list_of_hits_by_models={}
    tf_pval_list_of_hits={}
    tf_eval_list_of_hits={}
    tf_qval_list_of_hits={}
    tf_rank_list_of_hits={}
    for tf_pair,scores in tf_models.iteritems():
        tf_name_query,tf_model = tf_pair
        tomtom_obj             = tf_mdl_tomtom.get(tf_model)
        if not tf_nameseq.has_key(tf_name_query): continue
        tf_query=tf_nameseq.get(tf_name_query)
        tf_size_tomtom.setdefault(tf_name_query,float(tomtom_obj.get_size()))
        if not tf_has_motif.has_key(tf_query): continue
        motif_query = tf_has_motif.get(tf_query)
        for tomhit in tomtom_obj.get_hits():
            if tomhit.get_rank() < toprank:
               tf_list_of_hits_by_models.setdefault(tf_name_query,[]).append(tomhit.get_hit())
               if tomhit.get_hit() == motif_query:
                  tf_pval_list_of_hits.setdefault(tf_name_query,[]).append(float(scores[0]))
                  tf_eval_list_of_hits.setdefault(tf_name_query,[]).append(float(scores[1]))
                  tf_qval_list_of_hits.setdefault(tf_name_query,[]).append(float(scores[2]))
                  tf_rank_list_of_hits.setdefault(tf_name_query,[]).append(float(scores[3]))
        #tf_list_of_hits_by_models.setdefault(tf_query,[]).extend([tomhit.get_hit() for tomhit in tomtom_obj.get_hits() if tomhit.get_rank() < toprank])

    tf_enrichment={}
    for tf_name_query, list_of_hits in tf_list_of_hits_by_models.iteritems():
        if not tf_nameseq.has_key(tf_name_query): continue
        tf_query=tf_nameseq.get(tf_name_query)
        if not tf_has_motif.has_key(tf_query): continue
        motif_query = tf_has_motif.get(tf_query)
        #Get the collection of hits group= {tf_motif:times_found}
        group=Counter(list_of_hits)
        #Get the enrichment of the motif_query on the top selected results
        if group.has_key(motif_query) and tf_number_of_models.get(tf_name_query)>0:
           e=100*float(group.get(motif_query))/tf_number_of_models.get(tf_name_query)
        else:
           e=0.0
        count= 0
        n    = 0
        #Find the rank (top_rank) of the motif_query on the top selected results
        rank = tf_size_tomtom.get(tf_name_query)
        for hit,k in sorted(group.iteritems(),key=lambda x:x[1],reverse=True):
            if count != int(k): 
                count=int(k)
                n    =n+1
            if hit == motif_query: 
                rank = float(n)
        top_rank=int(tf_size_tomtom.get(tf_name_query)-rank+1)
        #Find the minimum enrichment with p-value lower than the random selection (1/number_of_motifs): lower than this implies the enrichment is as good as random
        if e>0 and group.has_key(motif_query) and tf_number_of_models.get(tf_name_query)>0:
          min_enrichment,worst_enrichment=best_enrichment(tf_number_of_models.get(tf_name_query),number_of_motifs,toprank)
          if verbose: sys.stdout.write("\t-- %s enrichment: %f (lower limit %f) number of models %d\n"%(tf_name_query,e,min_enrichment,tf_number_of_models.get(tf_name_query))) 
        else:
          min_enrichment=0
          worst_enrichment=0
        #Get the maximum number of hits
        max_enrichment=max([x for x in group.itervalues()])
        #Get the average of scores (p_value, e_value,q_value) with all the selected solutions in the enrichment set
        p_score=0.0
        e_score=0.0
        q_score=0.0
        ranking=0.0
        if tf_pval_list_of_hits.has_key(tf_name_query): p_score = sum([float(x) for x in tf_pval_list_of_hits.get(tf_name_query)])/len(tf_pval_list_of_hits.get(tf_name_query))
        if tf_eval_list_of_hits.has_key(tf_name_query): e_score = sum([float(x) for x in tf_eval_list_of_hits.get(tf_name_query)])/len(tf_eval_list_of_hits.get(tf_name_query))
        if tf_qval_list_of_hits.has_key(tf_name_query): q_score = sum([float(x) for x in tf_qval_list_of_hits.get(tf_name_query)])/len(tf_qval_list_of_hits.get(tf_name_query))
        if tf_rank_list_of_hits.has_key(tf_name_query): ranking = sum([float(x) for x in tf_rank_list_of_hits.get(tf_name_query)])/len(tf_rank_list_of_hits.get(tf_name_query))
        #store the data in a new dictionary named tf_enrichment
        if e<min_enrichment:
          tf_enrichment.setdefault((tf_name_query,"enrichment"),(p_score,e_score,q_score,ranking,e,0))
        else:   
          tf_enrichment.setdefault((tf_name_query,"enrichment"),(p_score,e_score,q_score,ranking,e,float(top_rank)))
    return tf_enrichment

def log_p_value(m,e,n,x):
    y=None
    k=(int)(m*e)
    if x<n and k<=m:
       fx = np.log(float(x)/float(n))
       fz = np.log(float(n-x)/float(n))
       y  = log_combinatorial(m,k)  +  k*fx  + (m-k)*fz 
    return y

def log_factorial(n):
    x=0
    if n>0:
       for i in xrange(1,int(n+1)): 
           x= x + np.log(float(i)) 
    return x

def log_combinatorial(n,m):
    x=0
    if n>m:
       x =  log_factorial(n) - log_factorial(m) - log_factorial(n-m)
    return x


def best_enrichment(m,n,x):
    max_LogpValue= -np.log(float(n))
    y = 1
    z = 1
    max_pvalue=0
    skip=False
    for p in xrange(100,1,-1):
        e= (float)(p)/100
        log_p  = log_p_value(m,e,n,x)
        if log_p is None: continue
        if np.exp(log_p)>max_pvalue: 
           z=p
           max_pvalue=np.exp(log_p)
        if log_p >= max_LogpValue and not skip:
           y=p
           skip=True
    return (y,z)

def boxplot(input_name,boxplot_name,title,output_name,verbose):
 from matplotlib.patches import Polygon
 if fileExist(input_name):
  dataname=[]
  dataplot=[]
  fp=open(input_name,"r")
  for line in fp:
   (name,datafile)=line.strip().split()
   if fileExist(datafile):
     dataname.append(name)
     datavalue=[]
     fi=open(datafile,"r")
     if verbose: sys.stdout.write("Open %s\n"%(datafile))
     for datainfo in fi:
       (info,value)=datainfo.strip().split()
       datavalue.append(float(value))
     dataplot.append(np.array(datavalue))
     fi.close()
   else:
     dataname.append(name)
     datavalue=[0.0]
     dataplot.append(np.array(datavalue))
  fp.close()
 else:
  raise IOError("Input File %s not found"%input_name)

 if len(dataplot)<=0: 
  raise IOError("Input Data of %s is Empty"%input_name)
 
 fig, ax1 = plt.subplots(figsize=(10,6))
 fig.canvas.set_window_title(title)
 plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
 bp = plt.boxplot(dataplot, notch=0, sym='+', vert=1, whis=1.5)
 plt.setp(bp['boxes'], color='black')
 plt.setp(bp['whiskers'], color='black')
 plt.setp(bp['fliers'], color='red', marker='+')
 ax1.set_axisbelow(True)
 ax1.set_title(title)
 ax1.set_xlabel('percentage of identity')
 ax1.set_ylabel('score values')
 boxColors = ['cyan','firebrick']
 numBoxes = len(dataname)
 medians = range(numBoxes)
 k=0
#Boxplot
 for i in range(numBoxes):
   box = bp['boxes'][i]
   boxX = []
   boxY = []
   for j in range(5):
      boxX.append(box.get_xdata()[j])
      boxY.append(box.get_ydata()[j])
   boxCoords = zip(boxX,boxY)
  # Group by color and hashed lines
   k = i % 2
   boxPolygon = Polygon(boxCoords, facecolor=boxColors[k])
   ax1.add_patch(boxPolygon)
  # Now draw the median lines back over what we just filled in
   med = bp['medians'][i]
   medianX = []
   medianY = []
   for j in range(2):
       medianX.append(med.get_xdata()[j])
       medianY.append(med.get_ydata()[j])
       plt.plot(medianX, medianY, color='green',linewidth=1.0)
       medians[i] = medianY[0]
  # Finally, overplot the sample averages, with horizontal alignment
  # in the center of each box
   plt.plot([np.average(med.get_xdata())], [np.average(dataplot[i])],
           color='white', marker='*', markeredgecolor='black')
# Set the axes ranges and axes labelsa
 ax1.set_xlim(0.5, numBoxes+0.5)
 top = max([ max([x for x in dataplot[i]]) for i in range(numBoxes) ])
 bottom = min([ min([x for x in dataplot[i]]) for i in range(numBoxes) ])
 ax1.set_ylim(bottom, top)
 xtickNames = plt.setp(ax1, xticklabels=dataname)
 plt.setp(xtickNames, rotation=45, fontsize=12)
 plt.savefig(boxplot_name)
 plt.close()

# Statistics: compare all vs all data
 ks={}
 mw={} 
 for i in range(numBoxes):
  for j in range(numBoxes):
   try:
    ks.setdefault((i,j),stats.ks_2samp(dataplot[i],dataplot[j]))
   except:
    ks.setdefault((i,j),(0.0,1.0))
    e = sys.exc_info()[0]
    #sys.stdout.write("Fail Kolmogorov Smirnov Test %s\n"%repr(e))
   try:
    mw.setdefault((i,j),stats.mannwhitneyu(dataplot[i],dataplot[j]))
   except:
    mw.setdefault((i,j),(0.0,0.5))
    e = sys.exc_info()[0]
    #sys.stdout.write("Fail Mann-Whitney Test %s\n"%e)
 fo=open(output_name,"w")
 fo.write("Kolmogorov-Smirnov Statistic\n") 
 fo.write("+--------------------------+\n") 
 for i in range(numBoxes):
  for j in range(numBoxes):
    if j>=i :
     fo.write("%s versus %s : Statistic %f p-value %f Average(%s)= %f Average(%s)= %f\n"%(dataname[i],dataname[j],ks[(i,j)][0],ks[(i,j)][1],dataname[i],np.average(dataplot[i]),dataname[j],np.average(dataplot[j])))

 fo.write("\nMann-Whitney U Statistic\n") 
 fo.write("+----------------------+\n") 
 for i in range(numBoxes):
  for j in range(numBoxes):
    if j>=i :
     fo.write("%s versus %s : Statistic %f p-value %f Average(%s)= %f Average(%s)= %f\n"%(dataname[i],dataname[j],mw[(i,j)][0],2*mw[(i,j)][1],dataname[i],np.average(dataplot[i]),dataname[j],np.average(dataplot[j])))
 fo.close()
 

#-------------#
# Main        #
#-------------#

def main():

    # Arguments & Options #
    options = parse_options()
    cdhit=config.get("Paths", "cd-hit")
    mmseqs=config.get("Paths", "mmseqs")
    verbose=options.verbose
    enrichment=options.enrichment
    # Initialize #
    if not os.path.exists(options.dummy_dir): os.makedirs(options.dummy_dir)
    tmp     =options.dummy_dir
    pwm_dir =options.pwm_dir
    pwm_out =options.output_file + "/pwms"
    seq_dir =options.output_file + "/fasta"
    outable =options.output_file + ".dat"
    dmdltab =options.output_file + "/models"
    dmdltmt =options.output_file + "/models_tomtom"
    outgraph=options.output_file + "_graph"
    outplot =options.output_file + "_plot"
    oTPgraph=options.output_file + "_TP_graph"
    oTPplot =options.output_file + "_TP_plot"
    if not os.path.exists(options.output_file): 
       os.makedirs(options.output_file)
       os.makedirs(pwm_out)
       os.makedirs(seq_dir)
       os.makedirs(dmdltab)
       os.makedirs(dmdltmt)
    if not os.path.exists(pwm_out):os.makedirs(pwm_out)
    if not os.path.exists(seq_dir):os.makedirs(seq_dir)
    if not os.path.exists(dmdltab):os.makedirs(dmdltab)
    if not os.path.exists(dmdltmt):os.makedirs(dmdltmt)
    plot_step      = options.plot_step
    selected_score = options.selected_score
    if selected_score > 3: enrichment=True
 
    #Read CisBP families
    if verbose: sys.stdout.write("-- Read CisBP families\n")
    cisbp_families={}
    for line in functions.parse_file(os.path.abspath(options.family_file)):
        m = re.search("\('(.+)', '(.+)', '.+', .+, .+\)", line)
        if m:           
            cisbp_families.setdefault(m.group(1).upper(),set()).add(m.group(2))
    fo=open("cisbp_families.txt","w")
    for cis,fam in cisbp_families.iteritems():
        fo.write("%s\t%s\n"%(cis, ";".join([str(x) for x in fam])))
    fo.close()

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
    fo=open("families.txt","w")
    for tf,fam in tf_families.iteritems():
        fo.write("%s\t%s\n"%(tf, ";".join([str(x) for x in fam]) ))
    fo.close()
    fo=open("families_cisbp.txt","w")
    for tf,fam in tf_families_cisbp.iteritems():
        fo.write("%s\t%s\n"%(tf, ";".join([str(x) for x in fam])))
    fo.close()

    #Create TF-Motif correspondence   
    if verbose: sys.stdout.write("-- Parsing input of TFs and their motifs\n")
    tf_has_motif,tf_nameseq,tf_nameseq_rev= parse_tfs_motifs(options.input_file,options.pbm_dir)
 
    #Create TF-model correspondence   
    if verbose: sys.stdout.write("-- Parsing input of TFs and their models\n")
    tf_has_model=parse_tfs_models(options.input_models)

    #Create sequence database
    fasta_file=seq_dir+"/tf.fa"
    fasta_db  =seq_dir+"/tf.fa.db"
    search_db =seq_dir+"/tf.search.db"
    compared  =seq_dir+"/tf.compared"
    pwm_db    =pwm_out+"/tf_motifs.db"
    if verbose: sys.stdout.write("Check files %s and %s\n"%(pwm_db,fasta_file))
    if not fileExist(pwm_db) or not fileExist(fasta_file):
     fa=open(fasta_file,"w")
     for tf,motif in tf_has_motif.iteritems():
       try:
        pwm_file = options.pwm_dir+"/"+motif
        tf_file  = options.pbm_dir+"/sequences/"+tf
        if fileExist(pwm_file) and fileExist(tf_file) and tf_nameseq_rev.has_key(tf):
           title_seq=tf_nameseq_rev.get(tf)[0]
           pwm = PWM.nMSA(pwm_file,motif,"txt")
           name= pwm_out+"/"+os.path.basename(pwm_file).rstrip(".txt")+".pwm"
           if not fileExist(name): pwm.write(name,"pwm")
           name= pwm_out+"/"+os.path.basename(pwm_file).rstrip(".txt")+".meme"
           if not fileExist(name): pwm.write(name,"meme")
           os.system("cat %s >> %s\n"%(name,pwm_db)) 
           seq = open(tf_file,"r")
           for line in seq:
             if line.startswith(">"):
               fa.write(">%s\n"%title_seq)
             else:
               fa.write("%s\n"%line.strip())
           seq.close()
        else:
           if verbose: sys.stderr.write("Files TF %s MOTIF %s are not found\n"%(tf,motif))
       except :
        sys.stderr.write("Error to parse TF %s MOTIF %s\n"%(tf,motif))
     fa.close()
     if not fileExist(compared):
       if verbose: sys.stdout.write("%s createdb  %s %s >& createdb.log\n"%(mmseqs,fasta_file, fasta_db))
       os.system("%s createdb  %s %s >& createdb.log"%(mmseqs,fasta_file, fasta_db))
       if verbose: sys.stdout.write("%s createindex %s %s >& createindex.log \n"%( mmseqs,fasta_db,tmp))
       os.system("%s createindex %s %s >& createindex.log"%( mmseqs,fasta_db,tmp))
       if verbose: sys.stdout.write("%s search %s %s %s %s --threads 1 -s  7.5 --max-seq-id 1.0 --num-iterations 4 >& search.log\n"%(mmseqs,fasta_db,fasta_db,search_db,tmp))
       os.system("%s search %s %s %s %s --threads 1 -s 7.5 --max-seq-id 1.0 --num-iterations 4 >& search.log"%(mmseqs,fasta_db,fasta_db,search_db,tmp))
       process = subprocess.Popen( [mmseqs, "convertalis", fasta_db, fasta_db, search_db,compared ],  stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
       process.communicate()
     else:
       if verbose: sys.stdout.write("Using file %s\n"%(compared))
    else:
     if verbose: sys.stdout.write("Using files %s and %s \n"%(fasta_file,pwm_db))

#Parse MMSEQ file into  dictionary TF ->{pairs (tf,%id)}
    if verbose: sys.stdout.write("-- Parsing TF x TF sequence comparison by MMseqs in %s\n"%compared)
    tf_neighbours=parse_mmseq_comparisons(compared)
    if tf_neighbours is None:
        exit(0)

#Compare Models
#For each TF query in  tf_neighbours 
    tf_models={}
    if enrichment: tf_enrichment={}
    for tf_name_query in tf_neighbours.iterkeys():
      omdltab = dmdltab + "/"+tf_name_query+".dat"
      omdltmt = dmdltmt + "/"+tf_name_query+".pickle"
      if fileExist(omdltab) and fileExist(omdltmt):
        if verbose: sys.stdout.write("Use table of models %s\n"%(omdltab))
        tf_models_tf=parse_table_models(omdltab)
        if verbose: sys.stdout.write("Use DUMP of TOTOM comparisons of models %s\n"%(omdltmt))
        tf_mdl_tomtom=cPickle.load(open(omdltmt,"rb"))
      else:
        if options.parallel:
         if modelExist(tf_name_query,tf_nameseq,tf_nameseq_rev,tf_has_motif,tf_has_model,verbose):
            program=os.path.join(scripts_path,"compare_tf_models_with_motifs.py")
            cluster_queue=None
            if config.get("Cluster", "cluster_queue") != "None": cluster_queue=config.get("Cluster", "cluster_queue")
            functions.submit_command_to_queue("%s %s -i %s -m %s -t %s -f %s --pbm=%s --pwm_model=%s --motif_database=%s --query=%s --table=%s --tom=%s --dummy=%s -v "%(os.path.join(config.get("Paths", "python_path"), "python"),program,options.input_file,options.input_models,options.tfs_file, options.family_file,options.pbm_dir,options.models_dir,pwm_db,tf_name_query,omdltab,omdltmt,options.dummy_dir), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
            continue
         else:
            continue
        else:
         tf_models_tf,tf_mdl_tomtom=create_tf_models(omdltab,omdltmt,pwm_db,tf_families,tf_families_cisbp,tf_name_query,tf_nameseq,tf_nameseq_rev,tf_has_motif,tf_has_model,options.models_dir,verbose,tmp)
      #Accumulate tf_models in a single dictionary
      for key,value in tf_models_tf.iteritems():
          tf_models.setdefault(key,value)
      #Calculate enrichment
      if enrichment:
         number_of_motifs=len([x for x in tf_neighbours.iterkeys()])
         tf_enrichment_tf=calculate_enrichment(tf_models_tf,tf_mdl_tomtom,tf_has_motif,tf_nameseq,options.top_enrichment,options.normalize,number_of_motifs,verbose)
         #Accumulate tf_enrichment in a single dictionary
         for key,value in tf_enrichment_tf.iteritems():
            tf_enrichment.setdefault(key,value)
#Skip next analysis if flags parallel or skip are ON
    if options.parallel:
       sys.stdout.write("Wait until all queues are done and run again. \n")
       exit(0)
    if options.skip_compare: 
       sys.stdout.write("Done. \n")
       exit(0)
#Average scores of Models and also maximum and minimum scores so we characterize the limits
    if enrichment:
      tf_average_models,tf_maxim_models,tf_minim_models=create_average_models(tf_enrichment)
    else:
      tf_average_models,tf_maxim_models,tf_minim_models=create_average_models(tf_models)

#Compare neighbouring motifs
    if fileExist(outable):
     if verbose: sys.stdout.write("Use table of comparisons %s\n"%(outable))
     tf_compare=parse_table(outable)
    else:
     if verbose: sys.stdout.write("Number of motifs to create Table of comparisons %d\n"%len([x for x in tf_neighbours.iterkeys()]))
     tf_compare=create_tf_comparisons(outable,dmdltmt,pwm_out,pwm_db,tf_neighbours,tf_nameseq,tf_nameseq_rev,tf_has_motif,tf_families,tf_families_cisbp,verbose,tmp)

#Calculate enrichment by neighbouring motifs
    if enrichment:
     tf_percentage_enrichment = calculate_enrichment_neighbours(outable,options.top_enrichment,options.normalize,plot_step,number_of_motifs,verbose)

#Make the graph distribution of neighbour analyses
    boxplot_exe  = config.get("Paths", "boxplot")
    python       = os.path.join(config.get("Paths", "python_path"), "python")
    nc           =  2
    nh           =  2*(100/options.plot_step) +1

#variables to Normalize scores
    set_of_families=set([clean_name("+".join(["_".join(str(x).split()) for x in fam])) for fam in tf_families.itervalues()])
    if enrichment:
        maximum_seq,minimum_seq,maximum_mdl,minimum_mdl=calculate_maxmin(tf_percentage_enrichment,tf_enrichment,tf_nameseq,tf_nameseq_rev,tf_families,set_of_families,selected_score)
    else:
        maximum_seq,minimum_seq,maximum_mdl,minimum_mdl=calculate_maxmin(tf_compare,tf_models,tf_nameseq,tf_nameseq_rev,tf_families,set_of_families,selected_score)

#Make Input Global
    input_name   =  outplot+"_global.dat"
    create_inputboxplot(input_name,"global",plot_step,outgraph,selected_score,verbose)
    if enrichment:
      found_interval=create_seqboxplot_global_enrichment(tf_percentage_enrichment,plot_step,outgraph,selected_score,options.normalize,maximum_seq.get("global"),minimum_seq.get("global"),verbose)
    else:
      found_interval=create_seqboxplot_global(tf_compare,plot_step,outgraph,selected_score,options.normalize,maximum_seq.get("global"),minimum_seq.get("global"),verbose)
    if options.average:
      create_mdlboxplot_global(tf_average_models,found_interval,plot_step,outgraph,selected_score,options.normalize,maximum_mdl.get("global"),minimum_mdl.get("global"),verbose)
    elif options.maximum:
      create_mdlboxplot_global(tf_maxim_models,found_interval,plot_step,outgraph,selected_score,options.normalize,maximum_mdl.get("global"),minimum_mdl.get("global"),verbose)
    elif options.minimum:
      create_mdlboxplot_global(tf_minim_models,found_interval,plot_step,outgraph,selected_score,options.normalize,maximum_mdl.get("global"),minimum_mdl.get("global"),verbose)
    else:
      if enrichment:
        create_mdlboxplot_global(tf_enrichment,found_interval,plot_step,outgraph,selected_score,options.normalize,maximum_mdl.get("global"),minimum_mdl.get("global"),verbose)
      else:
        create_mdlboxplot_global(tf_models,found_interval,plot_step,outgraph,selected_score,options.normalize,maximum_mdl.get("global"),minimum_mdl.get("global"),verbose)

#Make Inputs families
    for fam in tf_families.itervalues():
      family=clean_name("+".join(["_".join(str(x).split()) for x in fam]))
      input_name   =  outplot+"_"+family+".dat"
      if fileExist(input_name): continue
      if verbose: sys.stdout.write("Families PLOT input %s\n"%( input_name))
      create_inputboxplot(input_name,family,plot_step,outgraph,selected_score,verbose)
    if enrichment:
      create_seqboxplot_family_enrichment(tf_percentage_enrichment,tf_nameseq,tf_nameseq_rev,tf_families,plot_step,outgraph,selected_score,options.normalize,maximum_seq,minimum_seq,verbose)
    else:
      create_seqboxplot_family(tf_compare,tf_nameseq,tf_nameseq_rev,tf_families,plot_step,outgraph,selected_score,options.normalize,maximum_seq,minimum_seq,verbose)
    if options.average:
      create_mdlboxplot_family(tf_average_models,tf_nameseq,tf_nameseq_rev,tf_families,found_interval,plot_step,outgraph,selected_score,options.normalize,maximum_mdl,minimum_mdl,verbose)
    elif options.maximum:
      create_mdlboxplot_family(tf_maxim_models,tf_nameseq,tf_nameseq_rev,tf_families,found_interval,plot_step,outgraph,selected_score,options.normalize,maximum_mdl,minimum_mdl,verbose)
    elif options.minimum:
      create_mdlboxplot_family(tf_minim_models,tf_nameseq,tf_nameseq_rev,tf_families,found_interval,plot_step,outgraph,selected_score,options.normalize,maximum_mdl,minimum_mdl,verbose)
    else:
      if enrichment:
        create_mdlboxplot_family(tf_enrichment,tf_nameseq,tf_nameseq_rev,tf_families,found_interval,plot_step,outgraph,selected_score,options.normalize,maximum_mdl,minimum_mdl,verbose)
      else:
        create_mdlboxplot_family(tf_models,tf_nameseq,tf_nameseq_rev,tf_families,found_interval,plot_step,outgraph,selected_score,options.normalize,maximum_mdl,minimum_mdl,verbose)


#Make the plots
# Global
    input_name   =  outplot+"_global.dat"
    boxplot_name =  outplot+"_global"
    output_name  =  outplot+"_global.log" 
    if verbose: sys.stdout.write('\t-- boxplot INPUT %s OUTPUT %s STATISTICS %s\n'%(input_name,boxplot_name,output_name))
    try:
      boxplot(input_name,boxplot_name,"global",output_name,verbose)
    except Exception as e:
      if verbose: sys.stderr.write("Failed global %s\n"%e)
      
    #os.system('%s %s -l %s -o %s -nc %d -nh %d -t "Distribution of %s neighbours" >& %s'%(python,boxplot_exe,input_name,boxplot_name,nc,nh,"global",output_name))
# Families
    for fam in tf_families.itervalues():
        family=clean_name("+".join(["_".join(str(x).split()) for x in fam]))
        input_name   =  outplot+"_"+family+".dat"
        boxplot_name =  outplot+"_"+family
        output_name  =  outplot+"_"+family+".log"
        if fileExist(output_name) or fileExist(boxplot_name+".png"): continue
        if verbose: sys.stdout.write('\t-- boxplot INPUT %s OUTPUT %s STATISTICS %s\n'%(input_name,boxplot_name,output_name))
        try:
          boxplot(input_name,boxplot_name,family,output_name,verbose)
        except Exception as e:
          if verbose: sys.stderr.write("Failed %s %s\n"%(family,e))
        #os.system('%s %s -l %s -o %s -nc %d -nh %d -t "Distribution of %s neighbours" >& %s'%(python,boxplot_exe,input_name,boxplot_name,nc,nh,family,output_name))

# TPR curves
    if selected_score == 3 or selected_score == 5: 
       size_db_seq = maximum_seq.get("global")
       size_db_mdl = maximum_mdl.get("global")
       if options.normalize:
          size_db_seq = 100.0
          size_db_mdl = 100.0
       tpr_mdl=extract_TP_ratio(outgraph+"_global_","M.dat",size_db_mdl,options.top_threshold,options.plot_step)
       tpr_seq=extract_TP_ratio(outgraph+"_global_",".dat",size_db_seq,options.top_threshold,options.plot_step)
       fo = open (oTPplot+"_global.dat","w")
       plotname=oTPplot+"_global"
       title="True Positive ratio with all models"
       x=[]
       y=[]
       z=[]
       for middle in sorted([int(key) for key in tpr_mdl.iterkeys()]):
           if tpr_seq.has_key(middle) and tpr_mdl.has_key(middle):
              x.append(float(middle))
              y.append(tpr_seq.get(middle))
              z.append(tpr_mdl.get(middle))
              fo.write("%5d\t%10.5f\t%10.5f\n"%(middle,tpr_seq.get(middle),tpr_mdl.get(middle)))
       fo.close()
       #os.system('%s %s -i %s -o %s -x %s -y %s -t %s'%(python,"plotcurve.py",oTPplot+"_global.dat",plotname,"percentage of identity","ratio of true positives",plotname))
       try:
         plotcurve(x,y,z,plotname,title)
       except Exception as e:
         if verbose: sys.stdout.write('Failed to plot global TP ratio due to %s\n'%(e))
       tpr_lst_mdl={}
       tpr_lst_seq={}
       for family in set_of_families:
           size_db_seq = maximum_seq.get(family)
           size_db_mdl = maximum_mdl.get(family)
           if options.normalize:
              size_db_seq = 100.0
              size_db_mdl = 100.0
           tpr_mdl=extract_TP_ratio(outgraph+"_"+family+"_","M.dat",size_db_mdl,options.top_threshold,options.plot_step)
           tpr_seq=extract_TP_ratio(outgraph+"_"+family+"_",".dat",size_db_seq,options.top_threshold,options.plot_step)
           if len([x for x in tpr_mdl.iterkeys()]) > 0 and len([x for x in tpr_seq.iterkeys()]) > 0:
              fo = open (oTPplot+"_"+family+".dat","w")
              plotname=oTPplot+"_"+family
              title="True Positive ratio of "+family+" set"
              x=[]
              y=[]
              z=[]
              for middle in sorted([int(key) for key in tpr_mdl.iterkeys()]):
                 if tpr_seq.has_key(middle) and tpr_mdl.has_key(middle):
                    x.append(float(middle))
                    y.append(tpr_seq.get(middle))
                    z.append(tpr_mdl.get(middle))
                    tpr_lst_mdl.setdefault(middle,[]).append(tpr_mdl.get(middle))
                    tpr_lst_seq.setdefault(middle,[]).append(tpr_seq.get(middle))
                    fo.write("%5d\t%10.5f\t%10.5f\n"%(middle,tpr_seq.get(middle),tpr_mdl.get(middle)))
              fo.close()
              #os.system('%s %s -i %s -o %s -x %s -y %s -t %s'%(python,"plotcurve.py",oTPplot+"_"+family+".dat",plotname,"percentage of identity","ratio of true positives",plotname))
              try:
                 plotcurve(x,y,z,plotname,title)
              except Exception as e:
                 if verbose: sys.stdout.write('Failed to plot family %s TP ratio due to %s\n'%(family,e))
       fo=open(oTPplot+"_average.dat","w")
       a=[]
       b=[]
       c=[]
       plotname=oTPplot+"_average"
       title="Average of True Positives with all families"
       for middle in sorted([int(x) for x in tpr_lst_seq.iterkeys()]):
            ave_seq=sum([x for x in tpr_lst_seq.get(middle)])/len([x for x in tpr_lst_seq.get(middle)])
            ave_mdl=sum([x for x in tpr_lst_mdl.get(middle)])/len([x for x in tpr_lst_mdl.get(middle)])
            rms_seq=0.0
            rms_mdl=0.0
            a.append(float(middle))
            b.append(ave_seq)
            c.append(ave_mdl)
            try:
             rms_seq=np.sqrt(sum([x*x for x in tpr_lst_seq.get(middle)])/len([x for x in tpr_lst_seq.get(middle)])-ave_seq*ave_seq)
            except Exception as e:
             if verbose: sys.stdout.write("Missing families RMS in %d by %s\n"%(middle,e))
            try:
             rms_mdl=np.sqrt(sum([x*x for x in tpr_lst_mdl.get(middle)])/len([x for x in tpr_lst_mdl.get(middle)])-ave_mdl*ave_mdl)
            except Exception as e:
             if verbose: sys.stdout.write("Missing families RMS in %d by %s\n"%(middle,e))
            fo.write("%5d\t%10.5f\t%10.5f\t%10.5f\t%10.5f\n"%(middle,ave_seq,ave_mdl,rms_seq,rms_mdl))
       fo.close()
       try:
          plotcurve(a,b,c,plotname,title)
       except Exception as e:
          if verbose: sys.stdout.write('Failed to plot average TP ratio of families due to %s\n'%(e))
       #os.system('%s %s -i %s -o %s -x %s -y %s -t %s'%(python,"plotcurve.py",oTPplot+"_average.dat",oTPplot+"_average","percentage of identity","ratio of true positives",oTPplot+"_average"))

#Remove dummy folders
    if not verbose: shutil.rmtree(dummy_dir)
 


        
if __name__ == "__main__":
    main()
