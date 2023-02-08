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
import time

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

# Imports jbonet's module #
from SBI.data import aminoacids1to3, aminoacids_polarity_boolean, nitrogenous_bases
from SBI.structure import PDB

# Import my modules #
import contacts, dssp, interface,  spotentials, tmalign, triads, x3dna, scorer, model_protein, fimo, threader
import pwm_pbm as PWM

# Import my modules #
import profile as PROFILE


#-------------#
# Functions   #
#-------------#

def done_profiles(info_file):
    done=set()
    if not functions.fileExist(info_file): return done
    info=open(info_file,"r")
    for line in info:
        if line.startswith("#"):continue
        done.add(line.strip().split()[0])
    info.close()
    return done

def check_submitted(submitted,pdb_files):
    if len(submitted)>0: return (pdb_files!=pdb_files.intersection(submitted))
    else: return True

def check_done(done,name_of_profiles):
    if len(done)>0: return (name_of_profiles != name_of_profiles.intersection(done))
    else: return True

def classify_pdbs(pdb_files):

    class_monomer=[]
    class_dimer=[]
    for name_file in pdb_files:
        name      = os.path.basename(name_file)
        positions = name.split(":")
        if len(positions)<3: continue
        if len(positions)>3:class_dimer.append(name)
        else: class_monomer.append(name)

    number_of_clusters=0
    cluster={}
    original_classified=[]
    classified=[]
    iteration=0
    while(len(classified)<len(class_monomer) or len(original_classified)<len(classified)):
      original_classified=classified[:]
      iteration = iteration + 1
      for pdb_file in class_monomer:
        if pdb_file in classified: continue
        words=pdb_file.split(":")
        positions=[]
        for word in words:
            positions.append(word.split("_")[0])
        a  = int(positions[1])
        b  = int(positions[2])
        if len(cluster)<=0:
           number_of_clusters = number_of_clusters +1 
           cluster.setdefault(number_of_clusters,set()).add((pdb_file,int(a),int(b),0,0))
           classified.append(pdb_file)
           continue
        
        cluster_to_add=[]
        for number,values in cluster.iteritems():
            x=np.array([ int(p) for name,p,q,r,s in values ])
            y=np.array([ int(q) for name,p,q,r,s in values ])
            distance = max( abs(int(x.mean()) - a) , abs(int(y.mean()) - b) )
            set_c=set([z for z in range(int(x.mean()),int(y.mean()))])
            set_p=set([z for z in range(a,b)])
            ratio = 1.0
            if len(set_c)>0 and len(set_p)>0:
               ratio=1-float(len(set_c.intersection(set_p)))/min(len(set_c),len(set_p))
            if distance < 10 or ratio<0.2: 
                cluster_to_add.append((number,distance,ratio))
        if len(cluster_to_add)>0:
          mindi=min([dd for num,dd,rr in cluster_to_add])
          selected_clusters=[(num,rr) for num,dd,rr in cluster_to_add if dd==mindi]
          minra=min([rr for num,rr in selected_clusters])
          select=[num for num,rr in selected_clusters if rr==minra]
          number=select[0]
          cluster[number].add((pdb_file,int(a),int(b),0,0))
          classified.append(pdb_file)

      if original_classified==classified and iteration>3:
        skip=False
        for pdb_file in class_monomer:
           if pdb_file in classified: continue
           if skip :continue
           words=pdb_file.split(":")
           positions=[]
           for word in words:
               positions.append(word.split("_")[0])
           a  = int(positions[1])
           b  = int(positions[2])
           number_of_clusters = number_of_clusters +1 
           cluster.setdefault(number_of_clusters,set()).add((pdb_file,int(a),int(b),0,0))
           classified.append(pdb_file)
           skip=True
           iteration=0

    cluster2={}
    number_of_clusters=0
    original_classified=[]
    classified=[]
    iteration=0
    while(len(classified)<len(class_dimer) or len(original_classified)<len(classified)):
      original_classified=classified[:]
      iteration = iteration + 1
      for pdb_file in class_dimer:
        if pdb_file in classified: continue
        words=pdb_file.split(":")
        positions=[]
        for word in words:
            positions.append(word.split("_")[0])
        a  = int(positions[1])
        b  = int(positions[2])
        c  = d = 0
        if len(positions)>2: c  = int(positions[3])
        if len(positions)>3: d  = int(positions[4])
        if len(cluster2)<=0:
           number_of_clusters = number_of_clusters +1 
           cluster2.setdefault(number_of_clusters,set()).add((pdb_file,int(a),int(b),int(c),int(d)))
           classified.append(pdb_file)
           continue
        cluster_to_add=[]
        for number,values in cluster2.iteritems():
            x=np.array([ int(p) for name,p,q,r,s in values ])
            y=np.array([ int(q) for name,p,q,r,s in values ])
            u=np.array([ int(r) for name,p,q,r,s in values ])
            v=np.array([ int(s) for name,p,q,r,s in values ])
            distance_1 = max( abs(int(x.mean()) - a) , abs(int(y.mean()) - b) ) 
            distance_2 = max( abs(int(u.mean()) - c) , abs(int(v.mean()) - d) )
            set_c=set([z for z in range(int(x.mean()),int(y.mean()))])
            set_p=set([z for z in range(a,b)])
            ratio_1 = 1.0
            if len(set_c)>0 and len(set_p)>0:
               ratio_1 = 1-float(len(set_c.intersection(set_p)))/min(len(set_c),len(set_p))
            set_c=set([z for z in range(int(u.mean()),int(v.mean()))])
            set_p=set([z for z in range(c,d)])
            ratio_2 = 1.0
            if len(set_c)>0 and len(set_p)>0:
               ratio_2 = 1-float(len(set_c.intersection(set_p)))/min(len(set_c),len(set_p))
            if ( distance_1 < 10 and distance_2 < 10) or (ratio_1 < 0.2 and ratio_2 < 0.2):
                cluster_to_add.append((number,distance_1,ratio_1))
        if len(cluster_to_add)>0:
          mindi=min([dd for num,dd,rr in cluster_to_add])
          selected_clusters=[(num,rr) for num,dd,rr in cluster_to_add if dd==mindi]
          minra=min([rr for num,rr in selected_clusters])
          select=[num for num,rr in selected_clusters if rr==minra]
          number=select[0]
          cluster2[number].add((pdb_file,int(a),int(b),int(c),int(d)))
          classified.append(pdb_file)

      if original_classified==classified and iteration>3:
        skip=False
        for pdb_file in class_monomer:
           if pdb_file in classified: continue
           if skip :continue
           words=pdb_file.split(":")
           positions=[]
           for word in words:
               positions.append(word.split("_")[0])
           a  = int(positions[1])
           b  = int(positions[2])
           c  = d = 0
           if len(positions)>2: c  = int(positions[3])
           if len(positions)>3: d  = int(positions[4])
           number_of_clusters = number_of_clusters +1 
           cluster2.setdefault(number_of_clusters,set()).add((pdb_file,int(a),int(b),int(c),int(d)))
           classified.append(pdb_file)
           skip=True
           iteration=0

    cluster_all={}
    number_of_clusters = 0
    for key,values in cluster.iteritems():
        number_of_clusters = number_of_clusters + 1
        cluster_all.setdefault(number_of_clusters,values)
    for key,values in cluster2.iteritems():
        number_of_clusters = number_of_clusters + 1
        cluster_all.setdefault(number_of_clusters,values)

    cluster_domain={}
    for number_of_cluster in cluster_all.iterkeys():
      list_of_data = cluster_all.get(number_of_cluster)
      x=np.array([ int(a) for name,a,b,c,d in  list_of_data])
      y=np.array([ int(c) for name,a,b,c,d in  list_of_data])
      cluster_domain.setdefault(number_of_cluster,(x.mean(),y.mean()))

    order = [nc for nc,position in sorted(cluster_domain.items(),key=lambda x: x[1][0])]

    cluster_ordered={}

    number_of_cluster = 0
    for nc in order:
        number_of_cluster = number_of_cluster + 1
        cluster_ordered.setdefault(number_of_cluster,cluster_all.get(nc))

    return cluster_ordered

def calculate_protein_profiles(info_file,threading,pdb_files,name_of_profiles,label,output_dir,input_folder,dna_file,thresholds, energy_profile,pdb_file,pbm_dir,pdb_dir,families,potential_file, radius,fragment_restrict, binding_restrict, split_potential,auto_mode,family_potentials,pbm_potentials,score_threshold,taylor_approach,pmf,bins,known,meme,reset, dummy_dir,verbose,save,reuse,redofimo,complete,methylation=False): 

    # get dena sequences 
    dna_sequences=get_dna_sequences(dna_file)

    # Iterate untill all protein profiles are done
    submitted=set()
    n_done=0
    if verbose: print("Start iteration to check and run profiles")
    done=done_profiles(info_file)
    iterate = check_done(done,name_of_profiles)
    maxtime =  3600 * 3
    start_time = d_time = 0
    while( iterate ):
        for dna_name,dna_sequence,single_dna_file in dna_sequences:
          for pdb_file in pdb_files:
            # Define the output names for the protein profile
            if threading: output_pdb_name=pdb_file.rstrip(".txt")
            else:         output_pdb_name=pdb_file.rstrip(".pdb")
            if label is not None: output_pdb_name = output_pdb_name +"."+ label
            output_pdb_name = output_pdb_name + "." + dna_name
            output_file = os.path.join(output_dir,output_pdb_name)
            if output_pdb_name in submitted: continue
            # Skip if compressed results already exist
            if functions.fileExist(output_file+".pickle"):
                   submitted.add(output_pdb_name)
                   if verbose:print("\t-- Found protein profile %s"%output_pdb_name)
                   if redofimo and functions.fileExist(output_file+".meme"):
                      if verbose:print("\t\t-- Read profile %s"%output_pdb_name)
                      inp=open(output_file+".pickle","rb")
                      profile_protein=cPickle.load(inp)
                      inp.close()
                      new_pwm_file  = output_file+".meme"
                      new_logo_file = output_file+".logo"
                      if verbose:print("\t\t-- Read PWM %s"%(os.path.basename(new_pwm_file)))
                      new_msa_obj = PWM.nMSA(new_pwm_file,motif_name=output_pdb_name,option="meme")
                      if verbose:print("\t\t-- Set up MSA and logos of %s"%(os.path.basename(new_pwm_file)))
                      new_msa_obj.set_sequences()
                      new_msa_obj.write(file_name=output_file+".msa", option="msa", overwrite=True)
                      new_msa_obj.write(file_name=output_file+".pwm", option="pwm", overwrite=True)
                      PWM.write_logo(new_msa_obj,new_logo_file, dummy_dir)
                      if verbose:print("\t\t-- Run FIMO profiles")
                      #check the use of methylations
                      if meme: methylation_pwm = False
                      else:    methylation_pwm = methylation
                      #Modify the DNA sequence to standard only nucleotides 
                      if methylation_pwm:
                         dummy_dna_file = single_dna_file
                      else:
                         if not os.path.exists( os.path.join(dummy_dir,"DNA_FILES") ): os.makedirs(os.path.join(dummy_dir,"DNA_FILES"))
                         dummy_dna_file = os.path.join(dummy_dir,"DNA_FILES",os.path.basename(single_dna_file))
                         with open(dummy_dna_file,"w") as fo:
                           for header_dna, sequence_dna in functions.parse_fasta_file(single_dna_file):
                               sequence_dna_modified = sequence_dna.upper().upper().replace("X","C").replace("O","C").replace("J","G").replace("Q","G")
                               fo.write(">%s\n%s\n"%(header_dna,sequence_dna_modified))
                         fo.close()
                      for threshold in thresholds:
                          if verbose: sys.stdout.write("\t\t\t-- Set profile %s for threshold %s\n"%(output_pdb_name,threshold))
                          fimo_obj = fimo.get_fimo_obj(new_pwm_file, dummy_dna_file, float(threshold), None , dummy_dir)
                          profile_protein.set_fimo_scores_by_threshold(fimo_obj,threshold)
                      if verbose:print("\t\t-- Rewrite the new profile %s"%(os.path.basename(output_file)))
                      os.remove(output_file+".pickle")
                      out=open(output_file+".pickle","wb")
                      cPickle.dump(profile_protein,out)
                      out.close()
                   info = open(info_file,"r")
                   skip_adding=False
                   for line in info:
                       if output_pdb_name in line.split(): skip_adding=True
                   if skip_adding and verbose: print("\t\t-- Already in the information file as DONE")
                   if skip_adding: continue
                   if verbose: print("\t\t-- Add in the list of done")
                   info.close()
                   info = open(info_file,"a")
                   info.write("%s\tDONE\n"%output_pdb_name)
                   info.flush()
                   info.close()
                   continue
            else:
                   if redofimo and functions.fileExist(output_file+".meme"):
                      new_pwm_file  = output_file+".meme"
                      new_logo_file = output_file+".logo"
                      if verbose:print("\t\t-- Read PWM %s"%(os.path.basename(new_pwm_file)))
                      new_msa_obj = PWM.nMSA(new_pwm_file,motif_name=output_pdb_name,option="meme")
                      if verbose:print("\t\t-- Set up MSA and logos of %s"%(os.path.basename(new_pwm_file)))
                      new_msa_obj.set_sequences()
                      new_msa_obj.write(file_name=output_file+".msa", option="msa", overwrite=True)
                      new_msa_obj.write(file_name=output_file+".pwm", option="pwm", overwrite=True)
                      PWM.write_logo(new_msa_obj,new_logo_file, dummy_dir)
            #Skip if reusing the previous data
            if reuse:
                   if verbose:print("\t-- Not found protein profile %s"%output_pdb_name)
                   info = open(info_file,"r")
                   skip_adding=False
                   for line in info:
                       if output_pdb_name in line.split(): skip_adding=True
                   if skip_adding and verbose: print("\t-- Already in the information file but FAILED")
                   if skip_adding: submitted.add(output_pdb_name)
                   if skip_adding: continue
                   info.close()
            # Split calculation in parallel
            if parallel:
               # Skip if already submitted 
               if output_pdb_name in submitted:continue
               if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
               else: cluster_queue=config.get("Cluster", "cluster_queue")
               program=os.path.join(scripts_path,"single_profile.py")
               python=os.path.join(config.get("Paths", "python_path"), "python")
               parameters = " -i %s "%os.path.join(input_folder,pdb_file)
               parameters = parameters + " --dna %s "%single_dna_file
               parameters = parameters + " --pbm=%s "%pbm_dir
               parameters = parameters + " --pdb=%s "%pdb_dir
               parameters = parameters + " --dummy=%s "%dummy_dir
               parameters = parameters + " --output %s "%output_pdb_name
               parameters = parameters + " --output_dir %s"%output_dir
               parameters = parameters + " --info %s"%info_file
               parameters = parameters + " --potential %s "%split_potential
               parameters = parameters + " --energy %s "%energy_profile
               parameters = parameters + " --radius %f "%radius
               if fragment is not None        : parameters = parameters + " --fragment %s "%fragment
               if potential_file is not None  : parameters = parameters + " --file %s "%potential_file
               if score_threshold is not None : parameters = parameters + " --score_threshold %f "%score_threshold
               if threading        : parameters = parameters + " --threading "
               if pmf              : parameters = parameters + " --pmf "
               if auto_mode        : parameters = parameters + " --auto "
               if family_potentials: parameters = parameters + " --family "
               if pbm_potentials   : parameters = parameters + " -p "
               if known            : parameters = parameters + " --known "
               if taylor_approach  : parameters = parameters + " --taylor "
               if bins             : parameters = parameters + " --bins "
               if model_accuracy   : parameters = parameters + " --model_accuracy "
               if verbose          : parameters = parameters + " --verbose "
               if save             : parameters = parameters + " --save "
               if plot and save    : parameters = parameters + " --plot "
               if meme             : parameters = parameters + " --meme "
               if reset            : parameters = parameters + " --reset "
               if methylation      : parameters = parameters + " --methylation "
               if verbose:print("\t-- Submit protein profile calculation %s"%output_pdb_name)
               if verbose:print("\t\t-- %s %s %s" % (python,program,parameters))
               functions.submit_command_to_queue("%s %s %s" % (python,program,parameters), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
               submitted.add(output_pdb_name)
        
            # Linear calculation
            else:
               if verbose:print("\t-- Get protein profile %s"%output_pdb_name)
               input_file=os.path.join(input_folder,pdb_file)
               if threading:
                      profile=PROFILE.calculate_single_profile_of_thread(single_dna_file,thresholds,energy_profile,input_file,output_file,pbm_dir,pdb_dir,families,potential_file, radius,fragment_restrict, binding_restrict, split_potential,auto_mode,family_potentials,pbm_potentials,score_threshold,taylor_approach,pmf,bins,known,meme,reset, dummy_dir,verbose,save,methylation)
               else:
                   if model_accuracy:
                      profile_protein=PROFILE.calculate_single_profile_by_models(single_dna_file,thresholds,energy_profile,input_file,output_file,pbm_dir,pdb_dir,families,potential_file, radius,fragment_restrict, binding_restrict, split_potential,auto_mode,family_potentials,pbm_potentials,score_threshold,taylor_approach,pmf,bins,known,meme,reset, dummy_dir,verbose,save,methylation)
                   else:
                      profile_protein=PROFILE.calculate_single_profile_by_thread(single_dna_file,thresholds,energy_profile,input_file,output_file,pbm_dir,pdb_dir,families,potential_file, radius,fragment_restrict, binding_restrict, split_potential,auto_mode,family_potentials,pbm_potentials,score_threshold,taylor_approach,pmf,bins,known,meme,reset, dummy_dir,verbose,save,methylation)
               out=open(output_file+".pickle","wb")
               cPickle.dump(profile_protein,out)
               out.close()
               submitted.add(output_pdb_name)

        #Check next iteration, profiles submitted and profiles done
        done=done_profiles(info_file)
        iterate= check_done(done,name_of_profiles) 
        if len(done) > n_done:
           n_done=len(done)
           if n_done> 1 and start_time==0: start_time = float(time.time())
           if verbose: 
                sys.stdout.write("Number of profiles already done %d\n"%n_done)
                if d_time>0: sys.stdout.write("Time: %f\n"%d_time)
                sys.stdout.write("\t-- Check files done ...\n")
                info = open(info_file,"r")
                for line in info:
                   print("\t\t-- %s"%line.strip())
                info.flush()
                info.close()
                sys.stdout.write("\t-- Still running protein profiles  %s ...\n"%check_done(done,name_of_profiles))
                sys.stdout.write("\t-- Continue iteration %s ...\n"%iterate)
                sys.stdout.flush()
        #Check next iteration, if exceeding time and enough profiles stop iteration
        if start_time > 0:
           current_time = float(time.time())
           d_time = current_time - start_time
           #print("CHECK PROFILER TIME DONE %d TOTAL %d RATIO %f COMPLETE %f TIME %f MAXTIME %f"%(len(done),len(name_of_profiles),float(  len(done) ) / len(name_of_profiles),complete,d_time,maxtime))
           if float(  len(done) ) / len(name_of_profiles) > complete  and d_time > maxtime: iterate=False
           if d_time > maxtime and iterate and complete < 1.0: 
               if options.verbose: sys.stdout.write("Time: %f Done: %d (%d) Ratio to end: %f\n"%(d_time,len(done),len(name_of_profiles),complete))
               complete = complete - 0.001
               maxtime  = maxtime  + 300


    # Return checked files
    return done

def get_dna_sequences(dna_file):
    input_folder = os.path.dirname(os.path.abspath(dna_file))
    dna_sequences=[] 
    n_dna=0
    for header, sequence in functions.parse_fasta_file(dna_file):
        if verbose: sys.stdout.write("Read DNA named %s (lenght %d) \n"%(header,len(sequence)))
        dna_sequence = sequence
        dna_name = header.split()[0].replace("/","-")
        n_dna = n_dna + 1
        dummy_dna_file = os.path.join(input_folder,dna_name+"_"+str(n_dna)+".fa")
        dna_sequences.append((dna_name+"_"+str(n_dna),dna_sequence,dummy_dna_file))
    for dna_name,dna_sequence,dummy_dna_file in dna_sequences:
     if not functions.fileExist(dummy_dna_file):
        dna_dummy= open(dummy_dna_file,"w")
        dna_dummy.write(">%s\n%s\n"%(dna_name,dna_sequence))
        dna_dummy.close()
    return dna_sequences


#-------------#
# Options     #
#-------------#

def parse_options():
    '''
    This function parses the command line arguments and returns an optparse object.
    '''

    parser = optparse.OptionParser("Usage: xprofiler.py [--dummy=DUMMY_DIR] -i INPUT_FILE  -d DNA_FASTA [-l LABEL -o OUTPUT_NAME --info INFO_FILE --complete COMPLETE ] --pbm=PBM_dir --pdb=PDB_DIR [-v --save --plot --meme --reset] [--html --html_types HTML_SCORE_TYPES --html_energies HTML_ENERGIES]  [ --parallel --model_accuracy ] [-a -f -p -s SPLIT_POTENTIAL -e ENERGY_PROFILE -t THRESHOLD -k -b --taylor --file POTENTIAL --radius RADIUS --pmf --fragment FRAGMENT]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="DUMMY_DIR")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file with a list of folders containing PDB/Threading files. THIS INPUT IS MANDATORY", metavar="INPUT_FILE")
    parser.add_option("-l", action="store", type="string", default=None, dest="label", help="Label to organize the output files as 'label.energies.txt' ", metavar="LABEL")
    parser.add_option("--meme", default=False, action="store_true", dest="meme", help="Use 'uniprobe2meme' to calculate the PWM matrix for 'FIMO' , it removes non-standard nucleotides (default = False)")
    parser.add_option("-r","--reset", default=False, action="store_true", dest="reset", help="Clean the sequences of the original MSA and reset them by a random selection in accordance with the PWM (default = False)")
    parser.add_option("-o", "--output", default=None, action="store", type="string", dest="output_name", help="Output name for tables and plots (default is the name of the input FOLDER", metavar="OUTPUT_NAME")
    parser.add_option("--pbm", action="store", type="string", default=None, dest="pbm_dir", help="PBM directory (i.e. output dir from pbm.py). This is ,mandatory unless using option --file on potentials", metavar="PBM_DIR")
    parser.add_option("--pdb", action="store", type="string", default=None, dest="pdb_dir", help="PDB directory (i.e. output dir from pdb.py). This is ,mandatory unless using option --template", metavar="PDB_DIR")
    parser.add_option("-d","--dna",default=None, action="store", type="string", dest="dna_file", help="File of a DNA sequence in FASTA format to profile", metavar="FASTA")
    parser.add_option("--complete",default=1.00, action="store", type="float", dest="complete", help="Ratio of completness over the total number of profiles top be done(default= 0.95). This is useful in the server to stop the profiler when the time of execution exceeds more than 48 hours ", metavar="RATIO")
    parser.add_option("--fragment", default=None, action="store", type="string", dest="fragment", help="Fragment of protein to apply the potential. Format is 'a-b;c-d': two regions between residues a-b and c-d. (Default is None it applies to all amino-acids) !!! WARNING: all protein models MUST have the same CHAIN ID and numbering !!!")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")
    parser.add_option("--save", default=False, action="store_true", dest="save", help="Save PDB models and scores while scanning the DNA (default = False)")
    parser.add_option("--plot", default=False, action="store_true", dest="plot", help="Plot profiles (default = False)")
    parser.add_option("--html", default=False, action="store_true", dest="html", help="Plot profiles in HTML format (default = False)")
    parser.add_option("--html_header", default=False, action="store_true",  dest="html_header", help="Add a title in the plot header (default = False)")
    parser.add_option("--html_types", default="fimo_binding", action="store", dest="html_score_types", help="Plot in HTML the profiles of types selected (default = 'fimo_binding' ). Several profile-types are allowed if separated by coma: normal, energy, energy_best, energy_per_nucleotide, fimo_binding, fimo_score and fimo_log_score. If 'normal' is selected the statistical potentials selected will be normalized")
    parser.add_option("--html_energies", default="s3dc_dd", action="store", dest="html_energies", help="Plot in HTML  profiles of selected statistical potentials (default = 's3dc_dd' ). Several potentials are allowed if separated by coma: '3d','3dc','local','pair','s3dc','s3dc_di','s3dc_dd' ")
    parser.add_option("--info", default=None, action="store", type="string", dest="info_file", help="File to store information of SUCCESS/FAILURE of the run (default = standard output)")
    parser.add_option("--model_accuracy", default=False, action="store_true", dest="model_accuracy", help="Calculate the scores with 3D-models to obtain the highest accuracy on distances (default = False, it threads the sequence of each DNA fragment without recalculating distances and reduces the time in half)")
    parser.add_option("--parallel",default=False, action="store_true", dest="parallel", help="Run model profiles in parallel (default=False)")
    parser.add_option("--reuse",default=False, action="store_true", dest="reuse", help="Reuse the information files. If the flag is used then profiles that had failed will remain as FAILED, otherwise it tries to redo them (default=False)")
    parser.add_option("--redofimo",default=False, action="store_true", dest="redofimo", help="Redo FIMO results. If the flag is used then FIMO runs again, this is done when the PWM has been manually changed (default=False)")
    parser.add_option("--chains_fixed",default=False, action="store_true", dest="chains_fixed", help="Replace the names of the protein chains to A-B-C-D... and DNA chains to a-b-c-d....")
    parser.add_option("--threading",default=False, action="store_true", dest="threading", help="Use threading files instead of PDB files")

    group = optparse.OptionGroup(parser, "Statistical potentials", "Select your statistical potentials of choice. By default it uses S3DCdd general potential derived from PDB (the simplest one). In \"--auto\" mode, the program uses S3DCdd family potentials derived from both PDB and/or PBM data and/or approached by Taylor as selected in Potentials configuration file. In case family potentials cannot be applied, the program uses general potentials derived from both PDB and PBM data and approached by Taylor. \"-a\" option overrides options \"-f\", \"-p\" and \"-t\".")
    group.add_option("-a", "--auto", default=False, action="store_true", dest="auto_mode", help="Automate the selection of statistical potentials (default = False)")
    group.add_option("-f", "--family", default=False, action="store_true", dest="family_potentials", help="Use family potentials (default = False)")
    group.add_option("-p", default=False, action="store_true", dest="pbm_potentials", help="Use potentials derived from both PBM + PDB data (default = False)")
    group.add_option("-s", default="s3dc_dd", action="store", type="string", dest="split_potential", help="Split-potential to be used (all, 3d, 3dc, s3dc, s3dc_dd, s3dc_di, pair; default = s3dc_dd)", metavar="{string}")
    group.add_option("-e","--energy", default="all", action="store", type="string", dest="energy_profile", help="Select a specific Split-potential to be used for the profile (all, 3d, 3dc, s3dc, s3dc_dd, s3dc_di, pair; default = all)", metavar="{string}")
    group.add_option("-t", action="store", default=None, type="float", dest="score_threshold", help="Threshold on the scaled score to consider positive k-mers (default = 0.95)", metavar="{float}")
    group.add_option("-k","--known", default=False, action="store_true", dest="known", help="The name is of a known PDB file, with format 'code_chain' (default = False)")
    group.add_option("--taylor", default=False, action="store_true", dest="taylor_approach", help="Approach PMF by Taylor (default = False)")
    group.add_option("-b", "--bins", default=False, action="store_true",  dest="bins", help="Computate the potentials: by bins (if selected) or accumulative (default).", metavar="{boolean}")
    group.add_option("--file",default=None, action="store", type="string", dest="potential_file", help="Use potentials from specific file (default = None)", metavar="{string}")
    group.add_option("--radius",default=0, action="store", type="float", dest="radius", help="Maximum contact distance to calculate interactions (default=0 implies the use of 'max_contact_distance' from configuration", metavar="{string}")
    group.add_option("-m", "--pmf", default=False, action="store_true", dest="pmf", help="Use of raw mean-force potentials with no Z-scoring (default = False)")
    group.add_option("--methylation", default=False, action="store_true", dest="methylation", help="Use methylated cytosine specificities of binding/non-binding. Use option --meme to calculate the PWM with only standard nucleotides (default = False)")

    parser.add_option_group(group)

    
    (options, args) = parser.parse_args()

    if options.input_file is None or (options.pbm_dir is None and options.potential_file is None)  or options.pdb_dir is None:
        parser.error("missing mandatory arguments: type option \"-h\" for help")
    if not functions.fileExist(options.dna_file):
        parser.error("Missing DNA fasta file")
    if not re.search("^3d$|^3dc$|^local$|^pair$|^s3dc$|^s3dc_dd$|^s3dc_di$|^all$", options.split_potential):
        parser.error("incorrect value for -s argument: type option \"-h\" for help")
        exit(0)
    return options


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()
    verbose=options.verbose
    methylation = options.methylation
    threading=options.threading
    save=options.save
    plot=options.plot
    html=options.html
    chains_fixed=options.chains_fixed
    reuse=options.reuse
    redofimo=options.redofimo
    complete=float(options.complete)
    html_header=options.html_header
    html_score_types = [x for x in options.html_score_types.replace(" ","").split(",") if x!="normal"]
    html_energies    = options.html_energies.replace(" ","").split(",")
    html_normal      = ("normal" in options.html_score_types)
    dummy_dir = options.dummy_dir
    if not dummy_dir.startswith("/"): dummy_dir = os.path.abspath(options.dummy_dir)
    if not os.path.exists(dummy_dir): os.makedirs(dummy_dir)
    pdb_dir = options.pdb_dir
    if pdb_dir is not None:
       if not pdb_dir.startswith("/"): pdb_dir = os.path.abspath(options.pdb_dir)
    pbm_dir = options.pbm_dir
    if pbm_dir is not None:
       if not pbm_dir.startswith("/"): pbm_dir = os.path.abspath(options.pbm_dir)
    input_file=options.input_file
    if not input_file.startswith("/"): input_file = os.path.abspath(options.input_file)
    if options.output_name is not None: main_output_name = options.output_name
    else: main_output_name = os.path.basename(input_file)
    label=options.label
    if label is not None: main_output_name = main_output_name +"."+ label
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
    fragment          =options.fragment
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
    model_accuracy    =options.model_accuracy
    # Parallelize calculation
    parallel          =options.parallel
    # Define fimo_thresholds
    fimo_thresholds   =config.get("Parameters", "fimo_profile_thresholds").rstrip().split(",")

    #Get the DNA sequence (a single one in tf_profiler) 
    dna_sequences=[]
    if dna_file is not None: dna_sequences=get_dna_sequences( dna_file)
     
    #read input file of folders
    folders=[]
    fi=open(input_file,"r")
    for line in fi:
      if line.startswith("#"):continue
      folders.append(os.path.abspath(line.split()[0]))
    fi.close()

    #Make a dictionary with profiles of TF-domains and DNA targets =>  tf_profiles[(domain,dna)]
    tf_profiles={}
    tf_domains={}

    # For each folder
    for input_folder in folders:
            output_name = main_output_name +"_"+ os.path.basename(input_folder)
            # Check protein models of the folder
            pdb_all_files=set()
            if verbose: sys.stdout.write("Read PDB files from directory %s\n"%input_folder)
            for pdb_file in os.listdir(input_folder):
             try:
              if threading:
                 if verbose: sys.stdout.write("\t-- Add Threading file %s\n"%pdb_file)
              else:
                 if not pdb_file.endswith(".pdb"): continue
                 if chains_fixed:
                    pdb_obj=PDB(os.path.join(input_folder,pdb_file))
                    pdb_obj=model_protein.chains_fixed(pdb_obj)
                    pdb_obj.write(os.path.join(input_folder,pdb_file),force=True)
                 if verbose: sys.stdout.write("\t-- Add PDB file %s\n"%pdb_file)
              pdb_all_files.add(pdb_file)
             except:
              continue

            if threading:
                list_of_data=[]
                for threading_file in pdb_all_files:
                    thr_obj = threader.Threaded(threading_file=threading_file)
                    a= min([int(key[1]) for key in thr_obj.get_protein()])
                    b= max([int(key[1]) for key in thr_obj.get_protein()])
                    c=d=0
                    list_of_data.append((threading_file,a,b,c,d))
                cluster={0:list_of_data}
            else:
                cluster=classify_pdbs(pdb_all_files)

            for number_of_cluster in cluster.iterkeys():

                          list_of_data = cluster.get(number_of_cluster)
                          x=np.array([ int(a) for name,a,b,c,d in  list_of_data])
                          y=np.array([ int(b) for name,a,b,c,d in  list_of_data])
                          u=np.array([ int(c) for name,a,b,c,d in  list_of_data])
                          v=np.array([ int(d) for name,a,b,c,d in  list_of_data])

                          pdb_files=set() 
                          name_of_profiles=set()
                          name_of_profiles_per_dna={}
                          for pdb_file,a,b,c,d in list_of_data:
                            pdb_files.add(pdb_file)
                            for dna_name,dna_sequence,dummy_dna_file in dna_sequences:
                                     if threading: output_pdb_name=pdb_file.rstrip(".txt")
                                     else:         output_pdb_name=pdb_file.rstrip(".pdb")
                                     if label is not None: output_pdb_name = output_pdb_name +"."+ label
                                     output_pdb_name = output_pdb_name + "." + dna_name
                                     name_of_profiles.add(output_pdb_name)
                                     name_of_profiles_per_dna.setdefault(dna_name,set()).add(output_pdb_name)

                          profile_label_dir  = "."+str(int(x.mean()))+"_"+str(int(y.mean()))
                          if int(u.mean()) > 0 or int(v.mean()) > 0: profile_label_dir = profile_label_dir + "." + str(int(c)) + "_" + str(int(d))

            
                          output_dir=os.path.join(input_folder,output_name+profile_label_dir)
                          if not os.path.exists(output_dir): os.makedirs(output_dir)

                          #Start dictionary tf_profiles for cluster(domain)
                          tf_domain = os.path.basename(input_folder)+profile_label_dir
                          tf_domains.setdefault(os.path.basename(input_folder),[]).append(tf_domain)
                          

                          if options.info_file is not None: info_file  = os.path.join(output_dir,options.info_file)
                          else: info_file  = os.path.join(output_dir,"summary_"+output_name+profile_label_dir+".txt")

                          if reuse and functions.fileExist(info_file):
                            if verbose: print ("Reuse previous information of runs from file %s"%info_file)
                          else:
                            if verbose: print ("Open to write %s"%info_file)
                            info = open(info_file,"w")
                            info.write("#List of profiles\n")
                            info.close()

                          # Iterate untill all protein profiles are done
                          if len(dna_sequences)>0:
                            done=calculate_protein_profiles(info_file,threading,pdb_files,name_of_profiles,label,output_dir,input_folder,dna_file,fimo_thresholds, energy_profile,pdb_file,pbm_dir,pdb_dir,families,potential_file, radius,fragment_restrict, binding_restrict, split_potential,auto_mode,family_potentials,pbm_potentials,score_threshold,taylor_approach,pmf,bins,known,meme,reset, dummy_dir,verbose,save,reuse,redofimo,complete,methylation) 

                          #Profiles for each DNA sequence
                          profiles={}
                          for dna_name,dna_sequence,dummy_dna_file in dna_sequences:
                            if verbose:print("Making profiles with DNA %s\n\t-- %s"%(dna_name,dna_sequence))

                            # Define and Fill the profile
                            if not name_of_profiles_per_dna.has_key(dna_name):continue
                            profile = PROFILE.Profile(dna_sequence,thresholds=fimo_thresholds,potential=energy_profile)
                            for name in name_of_profiles_per_dna.get(dna_name):
                               name_profile=os.path.join(output_dir,name)
                               if verbose:print("\t\t-- Add protein profile %s "%name)
                               if not functions.fileExist(name_profile+".pickle"):continue
                               out=open(name_profile+".pickle","rb")
                               profile_protein = cPickle.load(out)
                               out.close()
                               if verbose: print("\t\t\t-- DNA in protein profile %s"%(profile_protein.get_dna()))
                               html_output = os.path.join(output_dir,name)
                               if html: 
                                                    if verbose: print("\t\t\t-- Plot HTML of protein profile %s"%name)
                                                    profile_protein.plot_html(html_score_types,html_normal,html_energies,html_output,html_header)
                                                    if verbose: print("\t\t\t-- Write CSV table of protein profile %s"%name)
                                                    profile_protein.write_table(html_output)
                               profile.add(profile_protein)

                            #Add the profile to the dictionaries
                            profiles.setdefault(dna_name,profile)
                            tf_profiles.setdefault((tf_domain,dna_name),profile)

                            #Analyse profiles with spepcific DNA
                            if len(profile.get_profiles())>0:
                              # Write table of MEAN in CSV format
                              output = os.path.join(output_dir,output_name+"."+dna_name+".mean")
                              if functions.fileExist(output+".csv"): os.remove(output+".csv")
                              if verbose:print("\t-- Write tables %s "%os.path.basename(output))
                              profile.write_mean_table(output+".csv")
                              # Write table of RMSD in CSV format
                              output = os.path.join(output_dir,output_name+"."+dna_name+".rmsd")
                              if functions.fileExist(output+".csv"): os.remove(output+".csv")
                              profile.write_rmsd_table(output+".csv")
                              # Plot profiles
                              output = os.path.join(output_dir,output_name+"."+dna_name+".average")
                              if verbose:print("\t-- Plot profiles of scores %s "%os.path.basename(output))
                              if plot: profile.plot(output)
                              if html: profile.plot_html(html_score_types,html_normal,html_energies,output,html_header)
                              output = os.path.join(output_dir,output_name+"."+dna_name+".accumulated")
                              if plot and save: profile.plot_all(output)
                            else:
                              if verbose: print("Check PDBs, all profiles have failed")

                          #Compare protein profiles plot in HTML
                          if html:
                            done_comparison=set()
                            for dna_name_one,dna_sequence_one,dummy_dna_file_one in dna_sequences:
                               for name_one in name_of_profiles_per_dna.get(dna_name_one):
                                  for dna_name_two,dna_sequence_two,dummy_dna_file_two in dna_sequences:
                                     if dna_name_one == dna_name_two: continue
                                     for name_two in name_of_profiles_per_dna.get(dna_name_two):
                                         protein_name_one = ".".join(name_one.split(".")[:-1])
                                         protein_name_two = ".".join(name_two.split(".")[:-1])
                                         if protein_name_one != protein_name_two: continue
                                         protein_name = protein_name_one
                                         if (protein_name,dna_name_one,dna_name_two) in done_comparison: continue
                                         name_profile_one=os.path.join(output_dir,name_one)
                                         name_profile_two=os.path.join(output_dir,name_two)
                                         if not functions.fileExist(name_profile_one+".pickle"):continue
                                         if not functions.fileExist(name_profile_two+".pickle"):continue
                                         if verbose: print("\t-- Plot HTML comparison of protein profiles %s (DNA: %s and %s)"%(protein_name,dna_name_one,dna_name_two))
                                         profile_one = PROFILE.Profile(dna_sequence_one,thresholds=fimo_thresholds,potential=energy_profile)
                                         out=open(name_profile_one+".pickle","rb")
                                         profile_protein_one = cPickle.load(out)
                                         out.close()
                                         profile_one.add(profile_protein_one)
                                         profile_two= PROFILE.Profile(dna_sequence_two,thresholds=fimo_thresholds,potential=energy_profile)
                                         out=open(name_profile_two+".pickle","rb")
                                         profile_protein_two = cPickle.load(out)
                                         out.close()
                                         profile_two.add(profile_protein_two)
                                         twoprofile=PROFILE.twoprofile(profile_one,profile_two,thresholds=fimo_thresholds,potential=energy_profile)
                                         html_output = os.path.join(output_dir,protein_name+"."+dna_name_one+"."+dna_name_two+".compare")
                                         twoprofile.plot_html(html_score_types,html_normal,html_energies,html_output,html_header)
                                         done_comparison.add((protein_name,dna_name_one,dna_name_two))
                                         done_comparison.add((protein_name,dna_name_two,dna_name_one))

                          #Compare profiles
                          done_comparison=set()
                          for dna_name_one,dna_sequence_one,dummy_dna_file_one in dna_sequences:
                             for dna_name_two,dna_sequence_two,dummy_dna_file_two in dna_sequences:
                                 profile_one = profiles.get(dna_name_one)
                                 profile_two = profiles.get(dna_name_two)
                                 if (dna_name_one,dna_name_two) in done_comparison: continue
                                 if dna_name_one == dna_name_two: continue
                                 if len(profile_one.get_profiles())<=0 or len(profile_one.get_profiles())<=0 : continue
                                 twoprofile=PROFILE.twoprofile(profile_one,profile_two,thresholds=fimo_thresholds,potential=energy_profile)
                                 # Write table of Statistics in CSV format
                                 output = os.path.join(output_dir,output_name+"."+dna_name_one+"."+dna_name_two+".stats")
                                 if verbose:print("\t-- Write tables %s "%os.path.basename(output))
                                 if functions.fileExist(output+".csv"): os.remove(output+".csv")
                                 try:
                                      twoprofile.write_stats_table(output+".csv")
                                 except:
                                      print("\t\tError: statistic comparison is not possible, check that the number of profiles is the same")
                                 # Write table of MEAN in CSV format
                                 output = os.path.join(output_dir,output_name+"."+dna_name_one+"."+dna_name_two+".mean")
                                 if functions.fileExist(output+".csv"): os.remove(output+".csv")
                                 try:
                                      twoprofile.write_mean_table(output+".csv")
                                 except:
                                      print("\t\tError: comparison of averaged profiles is not possible, check that all profiles are well defined")
                                 # Write table of RMSD in CSV format
                                 output = os.path.join(output_dir,output_name+"."+dna_name_one+"."+dna_name_two+".rmsd")
                                 if functions.fileExist(output+".csv"): os.remove(output+".csv")
                                 try:
                                      twoprofile.write_rmsd_table(output+".csv")
                                 except:
                                      print("\t\tError: comparison of standard deviation of profiles is not possible, check that all profiles are well defined")
                                 # Plot profiles
                                 output = os.path.join(output_dir,output_name+"."+dna_name_one+"."+dna_name_two+".statistical")
                                 if plot: 
                                                    try:
                                                            twoprofile.stplot(output)
                                                    except:
                                                            print("\t\tError: plot of statistic comparison is not possible, check that the number of profiles is the same")
                                 output = os.path.join(output_dir,output_name+"."+dna_name_one+"."+dna_name_two+".difference")
                                 if plot: 
                                                    try:
                                                            twoprofile.plot(output)
                                                    except:
                                                            print("\t\tError: plot average-comparison is not possible, check that all profiles are well defined")
                                 output = os.path.join(output_dir,output_name+"."+dna_name_one+"."+dna_name_two+".compare")
                                 if verbose:print("\t-- Plot the comparison of scored profiles %s "%os.path.basename(output))
                                 if plot: 
                                                    try:
                                                            twoprofile.twoplot(output)
                                                    except:
                                                            print("\t\tError: 2-plot comparison is not possible, check  that all profiles are well defined")
                                 if html: 
                                                    try:
                                                            twoprofile.plot_html(html_score_types,html_normal,html_energies,output,html_header)
                                                    except:
                                                            print("\t\tError: plot of HTML comparison is not possible, check  that all profiles are well defined")
                                 done_comparison.add((dna_name_one,dna_name_two))
                                 done_comparison.add((dna_name_two,dna_name_one))

    #Compare profiles between TFs
    print("Comparison of two or more TFs binding the same DNA")
    output_dir =  os.path.join(os.path.dirname(input_file), main_output_name+"_TF_comparison")
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    for dna_name,dna_sequence,dummy_dna_file in dna_sequences:
      done_comparison=set()
      print("-- DNA %s %s"%(dna_name,dna_sequence))
      for input_folder_one in folders:
         for input_folder_two in folders:
             tf_one = os.path.basename(input_folder_one)
             tf_two = os.path.basename(input_folder_two)
             if tf_one == tf_two: continue
             if not tf_domains.has_key(tf_one):continue
             if not tf_domains.has_key(tf_two):continue
             print("\t-- Compare %s vs %s"%(tf_one,tf_two))
             for domain_one in tf_domains.get(tf_one):
                 if not tf_profiles.has_key((domain_one,dna_name)): continue
                 profile_one = tf_profiles.get((domain_one,dna_name))
                 for domain_two in tf_domains.get(tf_two):
                     if not tf_profiles.has_key((domain_two,dna_name)): continue
                     if (domain_one,domain_two) in done_comparison: continue
                     profile_two = tf_profiles.get((domain_two,dna_name))
                     print("\t\t-- Domains %s vs %s"%(domain_one,domain_two))
                     if len(profile_one.get_profiles())<=0 or len(profile_one.get_profiles())<=0 : continue
                     twoprofile=PROFILE.twoprofile(profile_one,profile_two,thresholds=fimo_thresholds,potential=energy_profile)
                     # Write table of Statistics in CSV format
                     output = os.path.join(output_dir,main_output_name+"."+dna_name+"."+domain_one+"."+domain_two+".stats")
                     if verbose:print("\t\t\t-- Write tables %s "%os.path.basename(output))
                     if functions.fileExist(output+".csv"): os.remove(output+".csv")
                     try:
                                      twoprofile.write_stats_table(output+".csv")
                     except:
                                      print("\t\t\tError: statistic comparison is not possible, check that the number of profiles is the same")
                     # Write table of MEAN in CSV format
                     output = os.path.join(output_dir,main_output_name+"."+dna_name+"."+domain_one+"."+domain_two+".mean")
                     if functions.fileExist(output+".csv"): os.remove(output+".csv")
                     try:
                                      twoprofile.write_mean_table(output+".csv")
                     except:
                                      print("\t\t\tError: comparison of averaged profiles is not possible, check that all profiles are well defined")

                     # Write table of RMSD in CSV format
                     output = os.path.join(output_dir,main_output_name+"."+dna_name+"."+domain_one+"."+domain_two+".rmsd")
                     if functions.fileExist(output+".csv"): os.remove(output+".csv")
                     try:
                                      twoprofile.write_rmsd_table(output+".csv")
                     except:
                                      print("\t\t\tError: comparison of standard deviation of profiles is not possible, check that all profiles are well defined")
                     # Plot profiles
                     output = os.path.join(output_dir,main_output_name+"."+dna_name+"."+domain_one+"."+domain_two+".statistical")
                     if plot: 
                            try:
                                      twoprofile.stplot(output)
                            except:
                                      print("\t\t\tError: plot of statistic comparison is not possible, check that the number of profiles is the same")
                     output = os.path.join(output_dir,main_output_name+"."+dna_name+"."+domain_one+"."+domain_two+".difference")
                     if plot: 
                            try:
                                      twoprofile.plot(output)
                            except:
                                      print("\t\t\tError: plot average-comparison is not possible, check that all profiles are well defined")
                     output = os.path.join(output_dir,main_output_name+"."+dna_name+"."+domain_one+"."+domain_two+".compare")
                     if verbose:print("\t\t\t-- Plot the comparison of scored profiles %s "%os.path.basename(output))
                     if plot: 
                            try:
                                      twoprofile.twoplot(output)
                            except:
                                      print("\t\t\tError: 2-plot comparison is not possible, check  that all profiles are well defined")
                     if html: 
                            try:
                                      twoprofile.plot_html(html_score_types,html_normal,html_energies,output,html_header)
                            except:
                                      print("\t\t\tError: plot of HTML comparison is not possible, check  that all profiles are well defined")
                     done_comparison.add((domain_one,domain_two))
                     done_comparison.add((domain_two,domain_one))



    if verbose: print("Done")

