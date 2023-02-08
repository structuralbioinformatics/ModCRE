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
import pandas as pd
import shutil

# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.join(os.path.dirname(__file__),"../scripts"))

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

def plotcurve(p,a,b,c,out,title):
  x = np.array(p)
  u = np.array(a)
  v = np.array(b)
  w = np.array(c)
  plt.plot(x,u,'-', ms=5, lw=2, alpha=0.7, mfc='black')
  plt.plot(x,v,'-', ms=5, lw=2, alpha=0.7, mfc='cyan')
  plt.plot(x,w,'-', ms=5, lw=2, alpha=0.7, mfc='red')
  plt.xlabel('-log(p-value)')
  plt.ylabel('score')
  plt.title("%s"%title)
  plt.grid(True)
  plt.savefig(out)
  plt.close()

#-------------#
# Options     #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python nn_similar_PWMs.py -i tf_motifs_file -m tf_models_file --pwm=pwm_dir [-o output_dir]",
                                   epilog      = '@Oliva\'s lab 2018',
                                   description = "The program calculates the plots of the nmber of neighbours of each TF and the acceptable rank enrichmnt as a function of the P-value between motifs")

    parser.add_option("--dummy", default="./tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", default=None, help="File of correspondences between FastA sequences of TFs and their PWMs in CisBP", metavar="{filename}")
    parser.add_option("-f", action="store", type="string", dest="family_file", default=None, help="TFs Family file (from CIS-BP; i.e. cisbp_1.02.tf_families.sql)", metavar="{filename}")
    parser.add_option("-t", action="store", type="string", dest="tfs_file", default=None, help="TFs file (from CIS-BP; i.e. cisbp_1.02.tfs.sql)", metavar="{filename}")
    parser.add_option("-o", action="store", default="output", type="string", dest="output", help="Root name for output files and directories (default = 'output')", metavar="{rootname}")
    parser.add_option("--pwm", action="store", type="string", dest="pwm_dir", default=None, help="PWM folder with PWM motifs from CisBP (pbm/CisBP_2019/pwms)", metavar="{directory}")
    parser.add_option("--use_CisBP_family", default=False, action="store_true", dest="cisbp_family",help="Use family codes from CisBP instead of common family names of TFs (default = False)", metavar="{boolean}")
    parser.add_option("-p","--parallel",default=False, action="store_true", dest="parallel", help="Run in parallel the comparisons between the database of motifs and all models of each TF (this automatically forces 'skip' flag, default is False)")
    parser.add_option("-v","--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")

   
    (options, args) = parser.parse_args()
    if (options.input_file is None ) or  options.pwm_dir is None or  options.family_file is None or options.tfs_file is None :
        parser.error("missing arguments: type option \"-h\" for help")

    return options


#-------------#
# Main        #
#-------------#

def main():

    # Arguments & Options #
    options = parse_options()
    verbose=options.verbose
    family_file=options.family_file
    tfs_file=options.tfs_file
    pwm_dir=options.pwm_dir
    input_file=options.input_file
    output=options.output
    tmp   =options.dummy_dir
    # Initialize #
    if not os.path.exists(options.dummy_dir): os.makedirs(options.dummy_dir)
    if not os.path.exists(output): os.makedirs(output)
    pwm_db = os.path.join(output,"pwms")
    if not os.path.exists(pwm_db): os.makedirs(pwm_db)
    tomtom_dir = os.path.join(output,"tomtom")
    if not os.path.exists(tomtom_dir): os.makedirs(tomtom_dir)
    table_dir = os.path.join(output,"tables")
    if not os.path.exists(table_dir): os.makedirs(table_dir)
    plot_dir = os.path.join(output,"graphs")
    if not os.path.exists(plot_dir): os.makedirs(plot_dir)


    #Read CisBP families
    if verbose: sys.stdout.write("-- Read CisBP families\n")
    cisbp_families={}
    for line in functions.parse_file(os.path.abspath(options.family_file)):
        m = re.search("\('(.+)', '(.+)', '.+', .+, .+\)", line)
        if m:           
            cisbp_families.setdefault(m.group(1).upper(),set()).add(m.group(2))
    
    #Read CisBP TFs
    tf_families = {}
    families = set()
    for line in functions.parse_file(os.path.abspath(tfs_file)):
         m = re.search("\('(.+)', '(.+)', '.+', '.+', '.+', '(.+)', '[DIN]'\),*", line)
         if m:
             if options.cisbp_family:
                tf_families.setdefault(m.group(1).upper(), set()).add(m.group(2).upper())
                families.add(m.group(2).upper())
             else:
                tf_families.setdefault(m.group(1).upper(), set()).update(cisbp_families[m.group(2).upper()])
                families.update(cisbp_families[m.group(2).upper()])

    
    #Read TF-Motifs
    tf_motifs = {}
    motifs_tf = {}
    for line in functions.parse_file(os.path.abspath(input_file)):
        if line.startswith("#"):continue
        data=line.strip().split()
        tf=data[0].strip(".fa")
        pwm=data[1]
        tf_motifs.setdefault(tf,pwm)
        motifs_tf.setdefault(pwm,tf)


    #Make DB of PWMs
    motifs_db = os.path.join(pwm_db,"motifs_db.txt")
    if not fileExist(motifs_db):
       for tf,pwm in tf_motifs.iteritems():
           pwm_file = os.path.join(pwm_dir,pwm+".meme")
           pwm_copy = os.path.join(pwm_db,pwm+".meme")
           shutil.copy(pwm_file,pwm_copy)
           os.system("cat %s >> %s"%(pwm_file,motifs_db))

    #TOMTOM comparison
    pwm_matrix={}
    tf_matrix={}
    size = 0
    for tf,pwm in tf_motifs.iteritems():
        size = size + 1
        #if size >10: continue
        tomtom_file = os.path.join(tomtom_dir,tf+".tomtom")
        if verbose: sys.stdout.write("-- Check TF %s on database %s => TOMTOM = %s\n"%(tf,motifs_db,tomtom_file))
        if not  fileExist(tomtom_file):
          pwm_file = os.path.join(pwm_db,pwm+".meme")
          try:
           tomtom_obj  = TOMTOM.get_tomtom_obj(motifs_db, pwm_file, tmp)
          except Exception as e:
           sys.stderr.write("TOMTOM HIT %s %s Error %s\n"%(motifs_db, pwm_file,e))
           continue
          tomtom_obj.write(tomtom_file)
        else:
          with open(tomtom_file,"r") as tomtom_data: 
               tomtom_obj  = TOMTOM.Tomtom(tomtom_data.readlines())
          tomtom_data.close()
        for tomhit in tomtom_obj.get_hits():
            hit    = tomhit.get_hit()
            score  = 50.0
            escore = 50.0
            qscore = 50.0
            if float(tomhit.get_p_value()) > 1.0e-50 : score = -np.log(float(tomhit.get_p_value()))
            if float(tomhit.get_e_value()) > 1.0e-50 : escore = -np.log(float(tomhit.get_e_value()))
            if float(tomhit.get_q_value()) > 1.0e-50 : qscore = -np.log(float(tomhit.get_q_value()))
            #print("Get "+pwm+" "+hit+" scores ",score,escore,qscore)
            pwm_matrix.setdefault((pwm,hit),(score,escore,qscore))
            if motifs_tf.has_key(hit):
               tf_hit = motifs_tf[hit]
               tf_a =".".join(tf.split(".")[:-1])
               tf_b =".".join(tf_hit.split(".")[:-1])
               tf_matrix.setdefault((tf_a,tf_b),(score,escore,qscore))
               #print("\t in  "+tf+" "+tf_hit+" scores ",score,escore,qscore)
            else:
               print("Fail "+hit+" not found TF")

    
    #Get graph (global and per family)
    families.add("global")
    for family in families:
      code=clean_name(family)
      skip=True
      graph_neighbor={}
      for score_hundred in range(0,100,1):
        score=float(score_hundred)/10.0
        pglobe={}
        eglobe={}
        qglobe={}
        n_data=0
        for pair, sc in tf_matrix.iteritems():
            n_data=n_data+1
            tf_a,tf_b = pair
            if not tf_families.has_key(tf_a):
               print("Family not found for "+tf_a)
               continue
            #print(family,tf_families[tf_a])
            if family in tf_families[tf_a] or family=="global":
               pglobe.setdefault(tf_a,set())
               eglobe.setdefault(tf_a,set())
               qglobe.setdefault(tf_a,set())
               if float(sc[0]) >= score: pglobe[tf_a].add(tf_b)
               if float(sc[1]) >= score: eglobe[tf_a].add(tf_b)
               if float(sc[2]) >= score: qglobe[tf_a].add(tf_b)
        if len([x for tf,x in pglobe.iteritems()])<=0:
           print("Missing data for family %s out of %d checks at threshold %d"%(family,n_data,score))
           continue
        else:
           skip=False
           print("Family "+code+" done.")
        pminimum = min([len(x) for tf,x in pglobe.iteritems()])
        eminimum = min([len(x) for tf,x in eglobe.iteritems()])
        qminimum = min([len(x) for tf,x in qglobe.iteritems()])
        paverage = sum([float(len(x)) for tf,x in pglobe.iteritems()])/len([len(x) for tf,x in pglobe.iteritems()])
        eaverage = sum([float(len(x)) for tf,x in eglobe.iteritems()])/len([len(x) for tf,x in eglobe.iteritems()])
        qaverage = sum([float(len(x)) for tf,x in qglobe.iteritems()])/len([len(x) for tf,x in qglobe.iteritems()])
        pnormal  = 100*float(size-pminimum+1)/float(size)
        enormal  = 100*float(size-eminimum+1)/float(size)
        qnormal  = 100*float(size-qminimum+1)/float(size)
        pnormalX = 100*float(size-paverage+1)/float(size)
        enormalX = 100*float(size-eaverage+1)/float(size)
        qnormalX = 100*float(size-qaverage+1)/float(size)
        perror   = 100-100*float(size-pminimum+1)/float(size)
        eerror   = 100-100*float(size-eminimum+1)/float(size)
        qerror   = 100-100*float(size-qminimum+1)/float(size)
        perrorX  = 100-100*float(size-paverage+1)/float(size)
        eerrorX  = 100-100*float(size-eaverage+1)/float(size)
        qerrorX  = 100-100*float(size-qaverage+1)/float(size)
        graph_neighbor.setdefault("threshold",[]).append(score)
        graph_neighbor.setdefault("number_pV_min",[]).append(pminimum)
        graph_neighbor.setdefault("number_pV_mean",[]).append(paverage)
        graph_neighbor.setdefault("normal_rank_pV_min",[]).append(pnormal)
        graph_neighbor.setdefault("normal_rank_pV_mean",[]).append(pnormalX)
        graph_neighbor.setdefault("error_rank_pV_min",[]).append(perror )
        graph_neighbor.setdefault("error_rank_pV_mean",[]).append(perrorX )
        graph_neighbor.setdefault("number_EV_min",[]).append(eminimum)
        graph_neighbor.setdefault("number_EV_mean",[]).append(eaverage)
        graph_neighbor.setdefault("normal_rank_EV_min",[]).append(enormal)
        graph_neighbor.setdefault("normal_rank_EV_mean",[]).append(enormalX)
        graph_neighbor.setdefault("error_rank_EV_min",[]).append(eerror )
        graph_neighbor.setdefault("error_rank_EV_mean",[]).append(eerrorX )
        graph_neighbor.setdefault("number_qV_min",[]).append(qminimum)
        graph_neighbor.setdefault("number_qV_mean",[]).append(qaverage)
        graph_neighbor.setdefault("normal_rank_qV_min",[]).append(qnormal)
        graph_neighbor.setdefault("normal_rank_qV_mean",[]).append(qnormalX)
        graph_neighbor.setdefault("error_rank_qV_min",[]).append(qerror )
        graph_neighbor.setdefault("error_rank_qV_mean",[]).append(qerrorX )
      if not skip:
        #print(graph_neighbor)
        table_neighbor=pd.DataFrame(graph_neighbor)
        table_file   = os.path.join(table_dir,code+".csv")
        table_error  = os.path.join(table_dir,code+".error.csv")
        table_neigh  = os.path.join(table_dir,code+".neighbors.csv")
        plot_number  = os.path.join(plot_dir,code+".nn.png")
        plot_numberX = os.path.join(plot_dir,code+".nn.average.png")
        plot_normal  = os.path.join(plot_dir,code+".normal.png")
        plot_normalX = os.path.join(plot_dir,code+".normal.average.png")
        plot_error   = os.path.join(plot_dir,code+".error.png")
        plot_errorX  = os.path.join(plot_dir,code+".error.average.png")
        table_neighbor.to_csv(table_file)
        select=[c for c in table_neighbor.columns.values.tolist() if "error" in c and "mean" in c]
        select.append("threshold")
        table_neighbor_errors = table_neighbor[select]
        table_neighbor_errors.to_csv(table_error)
        select=[c for c in table_neighbor.columns.values.tolist() if "number" in c and "mean" in c]
        select.append("threshold")
        table_neighbor_number = table_neighbor[select]
        table_neighbor_number.to_csv(table_neigh)
        plotcurve(table_neighbor["threshold"].tolist(),table_neighbor["number_pV_min"].tolist(),table_neighbor["number_EV_min"].tolist(),table_neighbor["number_qV_min"].tolist(),plot_number,"Minimum number of neighbors")
        plotcurve(table_neighbor["threshold"].tolist(),table_neighbor["normal_rank_pV_min"].tolist(),table_neighbor["normal_rank_EV_min"].tolist(),table_neighbor["normal_rank_qV_min"].tolist(),plot_normal,"Acceptable rank enrichment")
        plotcurve(table_neighbor["threshold"].tolist(),table_neighbor["error_rank_pV_min"].tolist(),table_neighbor["error_rank_EV_min"].tolist(),table_neighbor["error_rank_qV_min"].tolist(),plot_error,"Error of rank enrichment")
        plotcurve(table_neighbor["threshold"].tolist(),table_neighbor["number_pV_mean"].tolist(),table_neighbor["number_EV_mean"].tolist(),table_neighbor["number_qV_mean"].tolist(),plot_numberX,"Average number of neighbors")
        plotcurve(table_neighbor["threshold"].tolist(),table_neighbor["normal_rank_pV_mean"].tolist(),table_neighbor["normal_rank_EV_mean"].tolist(),table_neighbor["normal_rank_qV_mean"].tolist(),plot_normalX,"Acceptable rank enrichment")
        plotcurve(table_neighbor["threshold"].tolist(),table_neighbor["error_rank_pV_mean"].tolist(),table_neighbor["error_rank_EV_mean"].tolist(),table_neighbor["error_rank_qV_mean"].tolist(),plot_errorX,"Error of rank enrichment")
      
        

        
if __name__ == "__main__":
    main()
        
