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
import pwm_pbm as PWM


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
import functions,tomtom



def fileExist(file):
    '''
    Check existing files
    '''
    if file is not None:
        return os.path.exists(file) and os.path.isfile(file)
    else:
        return False



def matches(a,b,offset,overlap):
    if offset<0:
        start=abs(offset)
        end=min(start+overlap+1,len(a))
        aa=a[start:end]
        end=min(overlap,len(b))
        bb=b[0:end]
    if offset>0:
        start=offset
        end=min(overlap,len(a))
        aa=a[0:end]
        end=min(start+overlap+1,len(b))
        bb=b[start:end]
    if offset==0:
        end=min(overlap,len(a))
        aa=a[0:end]
        end=min(overlap,len(b))
        bb=b[0:end]
    length=min(len(aa),len(bb))
    match=0
    for i in xrange(length):
        if aa[i]==bb[i]: 
            match=match+1
    return match

def get_data_tomtom(tf_mdl_tomtom,hit,query):

        tomtom_obj=tf_mdl_tomtom.get(query)
        pvalue=qvalue=evalue=1.0
        p_score=q_score=q_score=0.0
        overlap=offset=match=0
        rank=len(tomtom_obj.get_hits())
        seq_query=seq_hit=""

        if tomtom_obj.has_hit(hit):
           tomhit      = tomtom_obj.get_hit(hit)
           p_score = 20.0
           e_score = 20.0
           q_score = 20.0
           try:
              if float(tomhit.get_p_value()) > 1.0e-20 : p_score = -np.log(float(tomhit.get_p_value()))
              if float(tomhit.get_e_value()) > 1.0e-20 : e_score = -np.log(float(tomhit.get_e_value()))
              if float(tomhit.get_q_value()) > 1.0e-20 : q_score = -np.log(float(tomhit.get_q_value()))
           except Exception as e:
              sys.stderr.write("Skip TOMTOM HIT %s %s Error %s\n"%(query,hit,e))
           pvalue    = float(tomhit.get_p_value())
           qvalue    = float(tomhit.get_q_value())
           evalue    = float(tomhit.get_e_value())
           overlap   = tomhit.get_overlap()
           offset    = tomhit.get_offset()
           seq_query = tomhit.get_query_sequence()
           seq_hit   = tomhit.get_hit_sequence()
           rank      = tomhit.get_rank()
           match     = matches(seq_query,seq_hit,offset,overlap)

        
        x=np.array([pvalue,qvalue,evalue,p_score,q_score,q_score,seq_query,seq_hit,match,rank])

        return x

def logos_tf_models(outlogos, table_of_best_pv,table_of_best_match,PWM_modcre,PWM_cisbp,PWM_jaspar,dummy_dir):


    cisbp_families=set(table_of_best_pv["CisBP_Fam"].values.tolist())
    jaspar_motifs =set(table_of_best_pv["motif_jaspar"].values.tolist())
    for cisbp_fam in cisbp_families:

        output_dir=os.path.join(outlogos,cisbp_fam)
        if not os.path.exists(output_dir): os.makedirs(output_dir)

        tfid =table_of_best_pv[table_of_best_pv["CisBP_Fam"]==cisbp_fam]["TFID"]
        best_pv_jaspars    =set(table_of_best_pv[table_of_best_pv["CisBP_Fam"]==cisbp_fam]["motif_jaspar"].values.tolist())
        best_pv_cisbps     =set(table_of_best_pv[table_of_best_pv["CisBP_Fam"]==cisbp_fam]["motif_cisbp"].values.tolist())
        best_match_jaspars =set(table_of_best_match[table_of_best_match["CisBP_Fam"]==cisbp_fam]["motif_jaspar"].values.tolist())
        best_match_cisbps  =set(table_of_best_match[table_of_best_match["CisBP_Fam"]==cisbp_fam]["motif_cisbp"].values.tolist())
        #select pv_jaspar<1

        select=np.array(table_of_best_pv["pv_jaspar"].values.tolist(),dtype=float)<1
        best_pv_motifs     =[pair for pair in zip(table_of_best_pv[table_of_best_pv["CisBP_Fam"]==cisbp_fam]["motif_model"].values[select].tolist(),  table_of_best_pv[table_of_best_pv["CisBP_Fam"]==cisbp_fam]["motif_jaspar"].values[select].tolist())]

        #select match > 0

        select=np.array(table_of_best_pv["match_jaspar"].values.tolist(),dtype=float)>0
        best_match_motifs  =[pair for pair in zip(table_of_best_match[table_of_best_match["CisBP_Fam"]==cisbp_fam]["motif_model"].values[select].tolist(), table_of_best_match[table_of_best_match["CisBP_Fam"]==cisbp_fam]["motif_jaspar"].values[select].tolist())]
        
        for best_pv_motif,motif_jaspar in best_pv_motifs:
          if not fileExist(os.path.join(PWM_modcre,best_pv_motif)): continue
          msa_obj=PWM.nMSA(os.path.join(PWM_modcre,best_pv_motif),option="meme")
          PWM.write_logo(msa_obj,os.path.join(output_dir,best_pv_motif.rstrip(".meme")+".pv_"+motif_jaspar+".logo"),dummy_dir)

        for best_pv_jaspar in best_pv_jaspars:
          if not fileExist(os.path.join(PWM_jaspar,best_pv_jaspar+".meme")): continue
          msa_obj=PWM.nMSA(os.path.join(PWM_jaspar,best_pv_jaspar+".meme"),option="meme")
          PWM.write_logo(msa_obj,os.path.join(output_dir,best_pv_jaspar+".logo"),dummy_dir)

        for best_pv_cisbp in best_pv_cisbps:
          if not fileExist(os.path.join(PWM_cisbp,best_pv_cisbp)): continue
          msa_obj=PWM.nMSA(os.path.join(PWM_cisbp,best_pv_cisbp),option="txt")
          PWM.write_logo(msa_obj,os.path.join(output_dir,best_pv_cisbp.rstrip(".txt")+".logo"),dummy_dir)

        for best_match_motif,motif_jaspar in best_match_motifs:
          if not fileExist(os.path.join(PWM_modcre,best_match_motif)): continue
          msa_obj=PWM.nMSA(os.path.join(PWM_modcre,best_match_motif),option="meme")
          PWM.write_logo(msa_obj,os.path.join(output_dir,best_match_motif.rstrip(".meme")+".match_"+motif_jaspar+".logo"),dummy_dir)

        for best_match_jaspar in best_match_jaspars:
          if not fileExist(os.path.join(PWM_jaspar,best_match_jaspar+".meme")): continue
          msa_obj=PWM.nMSA(os.path.join(PWM_jaspar,best_match_jaspar+".meme"),option="meme")
          PWM.write_logo(msa_obj,os.path.join(output_dir,best_match_jaspar+".logo"),dummy_dir)

        for best_match_cisbp in  best_match_cisbps:
          if not fileExist(os.path.join(PWM_cisbp,best_match_cisbp)): continue
          msa_obj=PWM.nMSA(os.path.join(PWM_cisbp,best_match_cisbp),option="txt")
          PWM.write_logo(msa_obj,os.path.join(output_dir,best_match_cisbp.rstrip(".txt")+".logo"),dummy_dir)





def table_tf_models(omdltab,  tfid, tfinfo, tf_mdl_cisbp_tomtom, tf_mdl_jaspar_tomtom, cisbp_jaspar_tomtom, verbose,tmp):

    if verbose: print("Make tables for %s"%tfid)
    cisbp_families=set()
    for data in tfinfo:
        family,cisbp_fam,jaspar_id,tfmdl,motif_cisbp,motif_jaspar,motif_mdl=data
        cisbp_families.add(cisbp_fam)
    r=[]
    for data in tfinfo:
        family,cisbp_fam,jaspar_id,tfmdl,motif_cisbp,motif_jaspar,motif_mdl=data
        if not tf_mdl_cisbp_tomtom.has_key(tfmdl): continue
        if not tf_mdl_jaspar_tomtom.has_key(tfmdl): continue
        if not cisbp_jaspar_tomtom.has_key(motif_cisbp): continue
        x_init          = np.array([tfid])
        data_tmpl       = np.array(data)
        x = np.hstack((x_init,data_tmpl))
        if verbose: print("\t--Check model %s  vs CisBP Motif %s"%(tfmdl,motif_cisbp))
        mdl_vs_cisbp    = get_data_tomtom(tf_mdl_cisbp_tomtom,motif_cisbp,tfmdl)
        x = np.hstack((x,mdl_vs_cisbp))
        if verbose: print("\t--Check model %s  vs JASPAR Motif %s"%(tfmdl,motif_jaspar))
        mdl_vs_jaspar   = get_data_tomtom(tf_mdl_jaspar_tomtom,motif_jaspar,tfmdl)
        x = np.hstack((x,mdl_vs_jaspar))
        if verbose: print("\t--Check CisBP Motif %s  vs JASPAR Motif %s"%(motif_cisbp,motif_jaspar))
        cisbp_vs_jaspar = get_data_tomtom(cisbp_jaspar_tomtom,motif_jaspar,motif_cisbp)
        x = np.hstack((x,cisbp_vs_jaspar))
        r.append(x)
    result=np.array(r)
    cisbp_family_codes=result[:,2]
    jaspar_motif_codes=result[:,6]
    r_best_pv=[]
    r_best_match=[]
    r_data=[]
    column_data=["TFID","Family","CisBP_Fam", "JASPAR_ID","TF_model","motif_cisbp","motif_jaspar","motif_model","pv_cisbp","ev_cisbp","qv_cisbp","ps_cisbp","es_cisbp","qs_cisbp","seq_mdl_cisbp","seq_cisbp","match_cisbp","rk_cisbp","pv_jaspar","ev_jaspar","qv_jaspar","ps_jaspar","es_jaspar","qs_jaspar","seq_mdl_jaspar","seq_jaspar","match_jaspar","rk_jaspar","pv_jaspar_cisbp","ev_jaspar_cisbp","qv_jaspar_cisbp","ps_jaspar_cisbp","es_jaspar_cisbp","qs_jaspar_cisbp","seq_cisbp_jaspar","seq_jaspar_cisbp","match_jaspar_cisbp","rk_jaspar_cisbp"]
    column_parsed=["TFID","Family","CisBP_Fam", "JASPAR_ID","motif_cisbp","motif_jaspar","pv_mean","ev_mean","qv_mean","ps_mean","es_mean","qs_mean","rk_mean","match_mean","pv_std","ev_std","qv_std","ps_std","es_std","qs_std","rk_std","match_std","ratio_pv_jaspar","ratio_ev_jaspar","ratio_qv_jaspar","ratio_pv_cisbp","ratio_ev_cisbp","ratio_qv_cisbp","ps_jaspar_cisbp_mean","match_jaspar_cisbp_mean","ps_jaspar_cisbp_std","match_jaspar_cisbp_std"]
    done=set()
    for cisbp_fam in cisbp_families:
      result_cisbp=result[cisbp_family_codes==cisbp_fam]
      if len(result_cisbp) <=0: continue
      for jaspar_motif in jaspar_motif_codes:
        family_result=pd.DataFrame(result_cisbp[jaspar_motif_codes==jaspar_motif])
        family_result.columns=column_data
        if len(family_result.iloc[:,0].values) <=0: continue
        output_table_dir=os.path.join(omdltab,cisbp_fam)
        if not os.path.exists(output_table_dir): os.makedirs(output_table_dir)
        output=os.path.join(output_table_dir,tfid+"_"+jaspar_motif+".csv")
        if not fileExist(output):
           if verbose: print("\t-- Store table %s"%(output))
           family_result.to_csv(output)

        if (tfid,cisbp_fam,jaspar_motif) in done: continue

        best_pv_jaspar=family_result["pv_jaspar"].values.argmin()
        best_match_jaspar=family_result["match_jaspar"].values.argmax()
        
        data_best_pv    = family_result.iloc[best_pv_jaspar,:].values.tolist()
        data_best_match = family_result.iloc[best_match_jaspar,:].values.tolist()
        r_best_pv.append(data_best_pv)
        r_best_match.append(data_best_match)

        data_jaspar=family_result[["pv_jaspar","ev_jaspar","qv_jaspar","ps_jaspar","es_jaspar","qs_jaspar","rk_jaspar","match_jaspar"]].values
        data_cisbp =family_result[["pv_cisbp","ev_cisbp","qv_cisbp","ps_cisbp","es_cisbp","qs_cisbp","rk_cisbp","match_cisbp"]].values
        cmp_cisbp  =family_result[["ps_jaspar_cisbp","match_jaspar_cisbp"]].values
        data_info  =np.hstack((family_result.iloc[1,:4].values,family_result.iloc[1,5:7].values))
        data_scr   =np.array(data_jaspar,dtype=float)
        data_scrb  =np.array(data_cisbp,dtype=float)
        data_cmp   =np.array(cmp_cisbp,dtype=float)
        ratio_jaspar=np.array([float(sum(data_scr[:,0]<0.05))/len(data_scr[:,0]),float(sum(data_scr[:,1]<0.05))/len(data_scr[:,0]),float(sum(data_scr[:,2]<0.05))/len(data_scr[:,0])])
        ratio_cisbp =np.array([float(sum(data_scrb[:,0]<0.05))/len(data_scrb[:,0]),float(sum(data_scrb[:,1]<0.05))/len(data_scrb[:,0]),float(sum(data_scrb[:,2]<0.05))/len(data_scrb[:,0])])
        data_mean  =data_scr.mean(0)
        data_std   =data_scr.std(0)
        cmp_mean   =data_cmp.mean(0)
        cmp_std    =data_cmp.std(0)
        x =np.hstack((data_info,data_mean))
        x =np.hstack((x,data_std))
        x =np.hstack((x,ratio_jaspar))
        x =np.hstack((x,ratio_cisbp ))
        x =np.hstack((x,cmp_mean))
        x =np.hstack((x,cmp_std))
        r_data.append(x)
        done.add((tfid,cisbp_fam,jaspar_motif))


    result_best_pv=np.array(r_best_pv)
    result_best_match=np.array(r_best_match)
    result_data=np.array(r_data)
    table_of_best_pv=pd.DataFrame(columns=column_data)
    for i in range(len(result_best_pv)): 
       table_of_best_pv.loc[i]=result_best_pv[i]
    table_of_best_match=pd.DataFrame(columns=column_data)
    for i in range(len(result_best_match)):
       table_of_best_match.loc[i]=result_best_match[i]
    table_of_mean=pd.DataFrame(columns=column_parsed)
    for i in range(len(result_data)):
       table_of_mean.loc[i]=result_data[i]


#Use a filtered set to remove failed causes problems on those that fail match but not pv or fail pv but not match
#    result_best_pv_filtered    = result_best_pv[np.array(result_best_pv[:,18],dtype=float)<1]
#    result_best_match_filtered = result_best_match[np.array(result_best_match[:,26],dtype=float)>0]
#    result_data_filtered       = result_data[np.array(result_data[:,9],dtype=float)>0]
#    table_of_best_pv=pd.DataFrame(columns=column_data)
#    for i in range(len(result_best_pv_filtered)): 
#       table_of_best_pv.loc[i]=result_best_pv_filtered[i]
#    table_of_best_match=pd.DataFrame(columns=column_data)
#    for i in range(len(result_best_match_filtered)):
#       table_of_best_match.loc[i]=result_best_match_filtered[i]
#    table_of_mean=pd.DataFrame(columns=column_parsed)
#    for i in range(len(result_data_filtered)):
#       table_of_mean.loc[i]=result_data_filtered[i]


    return table_of_best_pv,table_of_best_match,table_of_mean




def dict_jaspar_cisbp(omdltmt, selection, jaspar_db, pwm_cisbp_dir, verbose,tmp):

      tf_mdl_tomtom={}
      
      for name in os.listdir(pwm_cisbp_dir):
          if name.endswith("meme"):
              cisbp_motif=name.rstrip("meme")+("txt")
              #if verbose: print("Check CisBP motif %s"%cisbp_motif)
              if cisbp_motif not in selection: continue
              pwm_cisbp = os.path.join(pwm_cisbp_dir,name)
              if verbose: print("Use CisBP motif meme %s"%name)
              try:
                  tomtom_obj  = TOMTOM.get_tomtom_obj(jaspar_db, pwm_cisbp , tmp)
              except Exception as e:
                  sys.stderr.write("TOMTOM HIT %s %s Error %s\n"%(jaspar_db,pwm_model,e))
                  continue
              tf_mdl_tomtom.setdefault(cisbp_motif,tomtom_obj)

      outmt=open(omdltmt,"wb")
      cPickle.dump(tf_mdl_tomtom,outmt)
      outmt.close()

      return tf_mdl_tomtom



def jaspar_tf_models( omdltmt, tfid,  tf_data, pwm_models_dir, jaspar_db, tmp):

      tf_mdl_tomtom={}

      for tfid_info in tf_data.get(tfid):
            family,cisbp_fam,jaspar_id,tfmdl,motif_cisbp,motif_jaspar,motif_mdl=tfid_info
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

    parser = optparse.OptionParser("python jaspar_comparison.py --jaspar_info=jaspar_info --cisbp_info=cisbp_info --nn_motifs=nn_dir --map=idmapping --seq_dir=tf_seq_dir --jaspar_pwms=jaspar_pwms_dir --cisbp_pwms=cisbp_pwms_dir --models_pwms=models_pwms_dir [--dummy -o output_dir]",
                                   epilog      = '@Oliva\'s lab 2018',
                                   description = "The program uses previous results from nearest neighbours and compares the PWM with JASPAR for the common sequences")

    parser.add_option("--dummy", default="./tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("--jaspar_info", action="store", type="string", dest="jaspar_info", help="JSON file of JASPAR with sequence, motif and UNiprot codes", metavar="{file}")
    parser.add_option("--cisbp_info", action="store", type="string", dest="cisbp_info", help="TF information from CisBP with family and TF ID codes", metavar="{file}")
    parser.add_option("--nn_motifs", action="store", type="string", dest="nn_dir", help="Folder with the output of the nearest neighbour analyses of CisBP", metavar="{directory}")
    parser.add_option("--map", action="store", type="string", dest="idmapping", help="File with cross-name sequence identifiers from UniProt", metavar="{file}")
    parser.add_option("--seq_dir", action="store", type="string", dest="tf_seq_dir", help="Folder with TF sequences in FastA format", metavar="{directory}")
    parser.add_option("--jaspar_pwms", action="store", type="string", dest="jaspar_pwms_dir", help="Folder with PWMS in meme format of JASPAR", metavar="{directory}")
    parser.add_option("--cisbp_pwms", action="store", type="string", dest="cisbp_pwms_dir", help="Folder with PWMS in meme format of CisBP", metavar="{directory}")
    parser.add_option("--models_pwms", action="store", type="string", dest="models_pwms_dir", help="Folder with PWMS in meme format of CisBP modelled with ModCRE", metavar="{directory}")
    parser.add_option("-o", action="store", type="string", dest="output_dir", help="Folder name for outputing the comparisons", metavar="{directory}")
    parser.add_option("-p","--parallel",default=False, action="store_true", dest="parallel", help="Run in parallel the comparisons between the database of motifs and all models of each TF (this automatically forces 'skip' flag, default is False)")
    parser.add_option("-v","--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")

    (options, args) = parser.parse_args()

    if options.models_pwms_dir is None:
        print("Missing models_pwms")
        exit()
    if options.cisbp_pwms_dir is None:
        print("Missing cisbp_pwms")
        exit()
    if options.jaspar_pwms_dir is None:
        print("Missing jaspar_pwms")
        exit()
    if options.nn_dir is None:
        print("Missing nn_dir")
        exit()
    if options.tf_seq_dir is None:
        print("Missing tf_seq_dir")
        exit()
    if options.idmapping is None:
        print("Missing idmapping")
        exit()
    if options.cisbp_info is None:
        print("Missing cisbp_info")
        exit()
    if options.jaspar_info is None:
        print("Missing jaspar_info")
        exit()
    if options.output_dir is None:
        print("Missing output_dir")
        exit()

    return options
#-------------#
# Main        #
#-------------#

def main():

    # Arguments & Options #
    options = parse_options()
    try:
     PWM_modcre   = options.models_pwms_dir
     if not PWM_modcre.startswith("/"):
         PWM_modcre=os.path.abspath(PWM_modcre)
     PWM_cisbp    = options.cisbp_pwms_dir
     if not PWM_cisbp.startswith("/"):
         PWM_cisbp=os.path.abspath(PWM_cisbp)
     PWM_jaspar   = options.jaspar_pwms_dir
     if not PWM_jaspar.startswith("/"):
         PWM_jaspar=os.path.abspath(PWM_jaspar)
     nn_dir       = options.nn_dir
     if not nn_dir.startswith("/"):
         nn_dir=os.path.abspath(nn_dir)
     tfmotifs_dir = os.path.join(nn_dir,"models") 
     tf_tomtom_dir= os.path.join(nn_dir,"models_tomtom")
     tf_seq_dir   = options.tf_seq_dir
     if not tf_seq_dir.startswith("/"):
         tf_seq_dir=os.path.abspath(tf_seq_dir)
     idmapping    = options.idmapping
     tfcisbp      = options.cisbp_info
     tfjaspar     = options.jaspar_info
     outdir       = options.output_dir
     if not outdir.startswith("/"):
         outdir=os.path.abspath(outdir)
     outjaspar    = os.path.join(outdir,"jaspar")
     outjs_cisbp  = os.path.join(outjaspar,"cisbp_tomtom")
     outjs_models = os.path.join(outjaspar,"models_tomtom")
     outmdltab    = os.path.join(outjaspar,"tables")
     outlogos     = os.path.join(outjaspar,"logos")
    except:
     print("Input error: check missing data")
     exit()
    tmp          = options.dummy_dir
    if not tmp.startswith("/"): tmp= os.path.abspath(options.dummy_dir)
    if not os.path.exists(tmp):  os.makedirs(tmp)
    verbose      = options.verbose

    # Output files
    if not os.path.exists(outjaspar):    os.makedirs(outjaspar)
    if not os.path.exists(outjs_cisbp):  os.makedirs(outjs_cisbp)
    if not os.path.exists(outjs_models): os.makedirs(outjs_models)
    if not os.path.exists(outmdltab):     os.makedirs(outmdltab)
    if not os.path.exists(outlogos):     os.makedirs(outlogos)


    #Addresses
    #PWM_modcre="/home/boliva/sit_sbi/ModCRE/nearest_neighbour/PWM_MODELS"
    #PWM_cisbp="/home/boliva/sit_sbi/ModCRE/pbm/CisBP_2019/pwms/"
    #PWM_jaspar="/home/boliva/sit_sbi/ModCRE/jaspar/Jaspar_2020/non-redundant/All_CORE/"
    #tfmotifs_dir="/home/boliva/sit_sbi/ModCRE/nearest_neighbour/TF_MOTIFS/models/"
    #tf_seq_dir="/home/boliva/sit_sbi/ModCRE/pbm_2.0/sequences/"

    #Files with IDs
    #idmapping="/home/boliva/sit_sbi/ModCRE/jaspar/Jaspar_2020/idmapping_selected.tab"
    #tfcisbp="/home/boliva/sit_sbi/ModCRE/pbm/CisBP_2019/TF_Information_PBM_Chip_SELEX.txt"
    #tfjaspar="/home/boliva/sit_sbi/ModCRE/jaspar/Jaspar_2020/JASPAR-profile-inference/files/all_core.uniprot.json"

    #Create data of JASPAR

    fl=open(tfjaspar,"r")
    jaspar_data=json.load(fl)
    fl.close()

    jaspar_ids=set([x.upper() for x in jaspar_data.iterkeys()])

    jaspar_db=os.path.join(PWM_jaspar,"jaspar_db.txt")
    if not fileExist(jaspar_db):
       for motif in os.listdir(PWM_jaspar):
           os.system("cat %s >> %s\n"%(motif,jaspar_db))

    # Data of CisBP

    fl=open(tfcisbp,"r")
    cisbp_ids={}
    for line in fl:
     if line.startswith("#"):continue
     data=line.split()
     tfid=data[0]
     tfname=data[1].upper()
     tfgene=data[4].upper()
     cisbp_ids.setdefault(tfname,tfid)
     cisbp_ids.setdefault(tfgene,tfid)
    fl.close()

    cisbp_db=os.path.join(PWM_cisbp,"cisbp_db.txt")
    if not fileExist(cisbp_db):
       for motif in os.listdir(PWM_cisbp): 
           pwm_file = os.path.join(PWM_cisbp,motif)
           if verbose: print("Add %s to DB %s"%(pwm_file,cisbp_db))
           try:
            pwm = PWM.nMSA(pwm_file,motif,"txt")
            name= pwm_file.rstrip("txt")+"meme"
            if not fileExist(name): pwm.write(name,"meme")
           except Exception as e:
            if verbose: print("Skip due to error %s"%e)
            continue
           os.system("cat %s >> %s\n"%(name,cisbp_db))

    

    #Cross-relation data cisbp-jaspar

    jaspar_codes={}
    cisbp_codes={}
    idmap_jaspar={}
    idmap_cisbp={}
    all_codes=set()
    jaspar_id_codes=set()
    cisbp_id_codes=set()
    jaspar_cisbp={}
    cisbp_jaspar={}

    fl=open(idmapping,"r")
    for line in fl:
     data_raw=line.split("\t")
     data = [x.upper() for x in data_raw]
     uac,uid,genid,refseq,gi,pdb,go,uref100,uref90,uref50 = data[:10]
     uac100=uref100.lstrip("UNIREF100_")
     uac90 =uref90.lstrip("UNIREF90_")
     uac50 =uref50.lstrip("UNIREF50_")
     unigene=data[14]
     embl,emblcds,ensembl,enstrs,enspro=data[15:20]
     found=False
     seq_codes=[uac,uid,genid,refseq,gi,uac100,uac90,uac50,unigene,embl,emblcds,ensembl,enstrs,enspro]
     for idmap_list in seq_codes:
      if found:continue
      for idmap in idmap_list.split(";"):
       if idmap.isalnum():
        if idmap in jaspar_ids and not found:
           jaspar_id=idmap
           found=True
     if found:
      jaspar_codes.setdefault(jaspar_id,seq_codes)
      for idmap_list in seq_codes:
       for idmap in idmap_list.split(";"):
        if idmap.isalnum():
         idmap_jaspar.setdefault(idmap,jaspar_id)
         all_codes.add(idmap)
         jaspar_id_codes.add(idmap)
     found=False
     for idmap_list in seq_codes:
       for idmap in idmap_list.split(";"):
         if idmap.isalnum():
           if cisbp_ids.has_key(idmap) and not found:
               cisbp_id=idmap
               found=True
     if found:
       cisbp_codes.setdefault(cisbp_id,seq_codes)
       for idmap_list in seq_codes:
          for idmap in idmap_list.split(";"):
            if idmap.isalnum():
              idmap_cisbp.setdefault(idmap,cisbp_id)
              all_codes.add(idmap)
              cisbp_id_codes.add(idmap)

    fl.close()


    for code in cisbp_id_codes.intersection(jaspar_id_codes):
     jaspar_id = idmap_jaspar.get(code)
     cisbp_id  = cisbp_ids.get(idmap_cisbp.get(code))
     jaspar_cisbp.setdefault( jaspar_id, set() ).add(cisbp_id)
     cisbp_jaspar.setdefault( cisbp_id, set() ).add(jaspar_id)
    


    # working data
    tfid_data={}
    if verbose:print("Reading data from nearest neighbours")
    list_of_tf_files = os.listdir(tfmotifs_dir)
    selected_cisbp_motifs=set()
    cisbp_to_classic_families={}
    for tf_file in list_of_tf_files:
      tf_short=".".join(tf_file.split(".")[:2])
      if not cisbp_jaspar.has_key(tf_short): continue 
      jaspar_id_set=cisbp_jaspar.get(tf_short)
      if verbose: print("Read models from %s"%tf_file)
      fl=open(os.path.join(tfmotifs_dir,tf_file),"r")
      for line in fl:
           if line.startswith("#"):continue
           data= line.split()
           tfid        = data[0]
           tfmdl       = data[1]
           motif_cisbp = data[2]
           motif_mdl   = data[3]
           family      = data[11]
           cisbp_fam   = data[12]
           seq_file    = os.path.join(tf_seq_dir,tfid+".fa")
           if verbose: print("Store families %s %s"%(cisbp_fam,family))
           cisbp_to_classic_families.setdefault(cisbp_fam,set()).add(family)
           tf_seq=""
           fseq=open(seq_file,"r")
           for line_seq in fseq:
             if line_seq.startswith(">"):continue
             tf_seq=tf_seq+line_seq.strip()
           fseq.close()
           for jaspar_id in jaspar_id_set:
             js_seq=jaspar_data[jaspar_id][1]
             if tf_seq not in js_seq or js_seq in tf_seq:
               for motif_jaspar in jaspar_data[jaspar_id][0]:
                 if not fileExist(os.path.join(PWM_jaspar,motif_jaspar+".meme")): continue
                 tfid_data.setdefault(tfid,set()).add((family,cisbp_fam,jaspar_id,tfmdl,motif_cisbp,motif_jaspar,motif_mdl))
                 selected_cisbp_motifs.add(motif_cisbp)
                 if verbose:
                    print("Add model %s to TFID data %s compared to motifs %s %s %s "%(motif_mdl,tfid,jaspar_id,motif_cisbp,motif_jaspar))
                    print("Add CisBP selected motif %s from family %s %s "%(motif_cisbp,cisbp_fam,family))
             else:
                if verbose: print("Different sequence JASPAR_ID %s  TFID %s SEQ %s %s"%(jaspar_id,tfid,js_seq,tf_seq)) 
                continue
      fl.close()
    nf=0
    table_families=open(os.path.join(outjaspar,"Table_families.txt"),"w")
    table_families.write("#   \t%10s\t%10s\n"%("CisBP name","Family"))
    for cisbf_fam,classic in cisbp_to_classic_families.iteritems():
         nf+=1
         table_families.write("%5d\t%10s"%(nf,cisbf_fam))
         for fam in classic:
             table_families.write("\t%10s"%fam)
         table_families.write("\n")
    table_families.close()

    submitted=set()

    cisbp_jaspar_pickle = os.path.join(outjs_cisbp,"cisbp_jaspar.pickle")
    if fileExist(cisbp_jaspar_pickle):
             if verbose: print("Use CisBP TOMTOMs on JASPAR %s"%cisbp_jaspar_pickle)
             cisbp_jaspar_tomtom = cPickle.load(open(cisbp_jaspar_pickle,"rb"))
    else:
             if verbose: print("Runing CisBP %s TOMTOMs on JASPAR %s"%(PWM_cisbp,jaspar_db))
             cisbp_jaspar_tomtom = dict_jaspar_cisbp(cisbp_jaspar_pickle, selected_cisbp_motifs , jaspar_db, PWM_cisbp, verbose,tmp)

    total_table_best_pv=[]
    total_table_best_match=[]
    total_table_mean=[]
    summary_table=[]
    for tfid,tfinfo in tfid_data.iteritems():
      #print("%s\t %s"%(tfid,str(tfinfo)))
           
      tfid_jaspar_pickle = os.path.join(outjs_models,tfid+".pickle")
      if fileExist(tfid_jaspar_pickle):
             if verbose: print("Use TOMTOM models %s on JASPAR DB"%tfid_jaspar_pickle)
             tf_mdl_jaspar_tomtom = cPickle.load(open(tfid_jaspar_pickle,"rb"))
      else:
             if options.parallel:
                if tfid_jaspar_pickle not in submitted:
                  program=os.path.join(scripts_path,"compare_model_jaspar.py")
                  cluster_queue=None
                  if config.get("Cluster", "cluster_queue") != "None": cluster_queue=config.get("Cluster", "cluster_queue")
                  options_tfid="--nn_motifs %s --jaspar_pwms %s --models_pwms %s --tfid %s --dummy %s -o %s "%(options.nn_dir,PWM_jaspar,PWM_modcre,tfid,tmp, outjs_models)
                  if verbose: print("Execute %s %s %s "%(os.path.join(config.get("Paths", "python_path"), "python"),program,options_tfid))
                  functions.submit_command_to_queue("%s %s %s "%(os.path.join(config.get("Paths", "python_path"), "python"),program,options_tfid), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),tmp,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                  submitted.add(tfid_jaspar_pickle)
                continue
             else:
                tf_mdl_jaspar_tomtom = jaspar_tf_models(tfid_jaspar_pickle , tfid,  tfid_data, pwm_models_dir, jaspar_db, tmp)
      tfid_cisbp_pickle=os.path.join(tf_tomtom_dir,tfid+".pickle")
      if not fileExist(tfid_cisbp_pickle): continue
      tf_mdl_cisbp_tomtom=cPickle.load(open(tfid_cisbp_pickle,"rb"))

      if verbose: print("Make Tables for %s"%tfid)
      table_of_best_pv,table_of_best_match,table_of_mean=table_tf_models(outmdltab, tfid, tfinfo, tf_mdl_cisbp_tomtom, tf_mdl_jaspar_tomtom, cisbp_jaspar_tomtom, verbose,tmp)
      if verbose: print("Make Logos for %s"%tfid)
      logos_tf_models(outlogos, table_of_best_pv,table_of_best_match,PWM_modcre,PWM_cisbp,PWM_jaspar,tmp)

      columns_best_pv=table_of_best_pv.columns
      columns_best_match=table_of_best_match.columns
      columns_mean=table_of_mean.columns

      
      if len( table_of_best_pv.iloc[:,:].values.tolist() )>0  and len( table_of_best_match.iloc[:,:].values.tolist() )>0:
         check=set([(x[0],x[5],x[6]) for x in total_table_best_pv])
         select=[ x not in check for x in zip(table_of_best_pv["TFID"].values.tolist(),table_of_best_pv["motif_cisbp"].values.tolist(),table_of_best_pv["motif_jaspar"].values.tolist())]
         if len(table_of_best_pv.iloc[:,:].values[select].tolist() )>0  and len( table_of_best_match.iloc[:,:].values[select].tolist() )>0:
             total_table_best_pv.extend(table_of_best_pv.iloc[:,:].values[select].tolist())
             total_table_best_match.extend(table_of_best_match.iloc[:,:].values[select].tolist())
             total_table_mean.extend(table_of_mean.iloc[:,:].values[select].tolist())
             short_table=table_of_mean[["TFID","Family","CisBP_Fam", "JASPAR_ID","motif_cisbp","motif_jaspar","ps_mean","ps_std","rk_mean","rk_std","match_mean","match_std","ps_jaspar_cisbp_mean","ps_jaspar_cisbp_std","match_jaspar_cisbp_mean","match_jaspar_cisbp_std","ratio_pv_jaspar","ratio_pv_cisbp"]]
             short_table.insert(6,"best_pv",table_of_best_pv["pv_jaspar"].values[select].tolist(),True)
             short_table.insert(6,"best_match",table_of_best_match["match_jaspar"].values[select].tolist(),True)
             short_table.assign(seq_best_match=table_of_best_match["seq_mdl_jaspar"].values[select].tolist())
             short_table.assign(seq_jaspar=table_of_best_match["seq_jaspar"].values[select].tolist())
             columns_summary=short_table.columns
             summary_table.extend(short_table.iloc[:,:].values.tolist())

      sys.stdout.flush()

    if len(submitted)>0:
          print("Wait until all runs are done")
          exit()

    if len(total_table_best_pv)>0:
       table_best_pv_total          =pd.DataFrame(total_table_best_pv)
       table_best_pv_total.columns  =columns_best_pv
    if len(total_table_best_match)>0:
       table_best_match_total         =pd.DataFrame(total_table_best_match)
       table_best_match_total.columns =columns_best_match
    if len(total_table_mean)>0:
       table_mean_total             =pd.DataFrame(total_table_mean)
       table_mean_total.columns     =columns_mean
    if len(summary_table)>0:
       table_summary                =pd.DataFrame(summary_table)
       table_summary.columns        =columns_summary   

    #Write Full tables

    out_best_pv=os.path.join(outmdltab,"table_best_pv.csv")
    out_best_match=os.path.join(outmdltab,"table_best_match.csv")
    out_mean=os.path.join(outmdltab,"table_averages.csv")
    out_summary=os.path.join(outmdltab,"table_summary.csv")

    if verbose: 
        print("Write Global Tables")
        print("\t-- %s"%out_best_pv)
        print("\t-- %s"%out_best_match)
        print("\t-- %s"%out_mean)
        print("\t-- %s"%out_summary)

    table_best_pv_total.to_csv(out_best_pv)
    table_best_match_total.to_csv(out_best_match)
    table_mean_total.to_csv(out_mean)
    table_summary.to_csv(out_summary)

    #Write Specific family tables

    summary_per_family={}
    summary_per_family_motif={}
    for cisbp_fam,cisbp_motif,jaspar_motif in zip(table_summary["CisBP_Fam"].values.tolist(),table_summary["motif_cisbp"].values.tolist(),table_summary["motif_jaspar"].values.tolist()):

      out_best_pv          =os.path.join(outmdltab,cisbp_fam,"table_best_pv_"+cisbp_motif.rstrip(".txt")+"_"+jaspar_motif+".csv")
      out_best_match       =os.path.join(outmdltab,cisbp_fam,"table_best_match_"+cisbp_motif.rstrip(".txt")+"_"+jaspar_motif+".csv")
      out_mean             =os.path.join(outmdltab,cisbp_fam,"table_averages_"+cisbp_motif.rstrip(".txt")+"_"+jaspar_motif+".csv")
      out_summary          =os.path.join(outmdltab,cisbp_fam,"table_summary_"+cisbp_motif.rstrip(".txt")+"_"+jaspar_motif+".csv")

      if verbose: 
        print("Write Tables per CisBP family")
        print("\t-- %s"%out_best_pv)
        print("\t-- %s"%out_best_match)
        print("\t-- %s"%out_mean)
        print("\t-- %s"%out_summary)

     
      best_pv_fam       =table_best_pv_total[table_best_pv_total.loc[:,"CisBP_Fam"]==cisbp_fam]
      best_pv_fam_cb    =best_pv_fam[best_pv_fam.loc[:,"motif_cisbp"]==cisbp_motif]
      table_best_pv_fam =best_pv_fam_cb[best_pv_fam.loc[:,"motif_jaspar"]==jaspar_motif]
     
      best_match_fam       =table_best_match_total[table_best_match_total.loc[:,"CisBP_Fam"]==cisbp_fam]
      best_match_fam_cb    =best_pv_fam[best_match_fam.loc[:,"motif_cisbp"]==cisbp_motif]
      table_best_match_fam =best_pv_fam_cb[best_match_fam.loc[:,"motif_jaspar"]==jaspar_motif]

      mean_fam       =table_mean_total[table_mean_total.loc[:,"CisBP_Fam"]==cisbp_fam]
      mean_fam_cb    =mean_fam[mean_fam.loc[:,"motif_cisbp"]==cisbp_motif]
      table_mean_fam =mean_fam_cb[mean_fam.loc[:,"motif_jaspar"]==jaspar_motif]

      summary_fam       =table_summary[table_summary.loc[:,"CisBP_Fam"]==cisbp_fam]
      summary_fam_cb    =summary_fam[summary_fam.loc[:,"motif_cisbp"]==cisbp_motif]
      table_summary_fam =summary_fam_cb[summary_fam.loc[:,"motif_jaspar"]==jaspar_motif]

      for data in table_summary_fam[["TFID","ps_mean","ps_jaspar_cisbp_mean","ratio_pv_jaspar","ratio_pv_cisbp"]].values.tolist():
          summary_per_family_motif.setdefault((cisbp_fam,cisbp_motif,jaspar_motif),[]).append(data)
          summary_per_family.setdefault(cisbp_fam,[]).append(data)


      table_best_pv_fam.to_csv(out_best_pv)
      table_best_match_fam.to_csv(out_best_match)
      table_mean_fam.to_csv(out_mean)
      table_summary_fam.to_csv(out_summary)


    #Summary per family

    table_summary_per_family_motif_list=[]
    for motifs,data_list in summary_per_family_motif.iteritems():
        data_all=np.array(data_list)
        data=np.array(data_all[:,1:],dtype=float)
        cisbp_fam,cisbp_motif,jaspar_motif=motifs
        tfids=data_all[:,0].tolist()
        number_of_tfs        =len(set(tfids))
        ps_mean,ps_jaspar_cisbp_mean,ratio_jaspar_mean,ratio_cisbb_mean=data.mean(0)
        success_ps           =sum(data[:,0]>2.99)
        success_r_cisbp      =sum(data[:,2]>0.5)
        success_r_jaspar     =sum(data[:,3]>0.5)
        success_min_cisbp    =sum(data[:,2]>0)
        success_min_jaspar   =sum(data[:,3]>0)
        table_summary_per_family_motif_list.append([cisbp_fam,cisbp_motif,jaspar_motif,number_of_tfs,ps_mean,ps_jaspar_cisbp_mean,ratio_jaspar_mean,ratio_cisbb_mean,success_ps,success_r_cisbp,success_r_jaspar,success_min_cisbp,success_min_jaspar])
    table_summary_per_family_motif=pd.DataFrame(table_summary_per_family_motif_list)
    table_summary_per_family_motif.columns=["cisbp_fam","cisbp_motif","jaspar_motif","#TFs","ps_mean","ps_jaspar_cisbp_mean","ratio_jaspar_mean","ratio_cisbb_mean","success_ps","success_ratio_cisbp","success_ratio_jaspar","minimum_success_cisbp","minimum_success_jaspar"]
    output_global=os.path.join(outmdltab,"table_summary_per_family_motif.csv")
    table_summary_per_family_motif.to_csv(output_global)

    table_summary_per_family_list=[]
    for cisbp_fam,data_list in summary_per_family.iteritems():
        data_all=np.array(data_list)
        data=np.array(data_all[:,1:],dtype=float)
        tfids=data_all[:,0].tolist()
        number_of_tfs        =len(set(tfids))
        number_of_instances  =len(data[:,0])
        ps_mean,ps_jaspar_cisbp_mean,ratio_jaspar_mean,ratio_cisbb_mean=data.mean(0)
        success_ps           =sum(data[:,0]>2.99)
        success_r_cisbp      =sum(data[:,2]>0.5)
        success_r_jaspar     =sum(data[:,3]>0.5)
        success_min_cisbp    =sum(data[:,2]>0)
        success_min_jaspar   =sum(data[:,3]>0)
        table_summary_per_family_list.append([cisbp_fam,number_of_tfs,number_of_instances,ps_mean,ps_jaspar_cisbp_mean,ratio_jaspar_mean,ratio_cisbb_mean,success_ps,success_r_cisbp,success_r_jaspar,success_min_cisbp,success_min_jaspar])
    table_summary_per_family=pd.DataFrame(table_summary_per_family_list)
    table_summary_per_family.columns=["cisbp_fam","#TFs","#Instances","ps_mean","ps_jaspar_cisbp_mean","ratio_jaspar_mean","ratio_cisbb_mean","success_ps","success_ratio_cisbp","success_ratio_jaspar","minimum_success_cisbp","minimum_success_jaspar"]
    output_global=os.path.join(outmdltab,"table_summary_per_family.csv")
    table_summary_per_family.to_csv(output_global)

    print("Done")

        
if __name__ == "__main__":
    main()







