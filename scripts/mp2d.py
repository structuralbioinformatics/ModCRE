import os, sys, re
import ConfigParser
import json
import optparse
import shutil
import subprocess
import imp
import math
import numpy as np
import pickle
import time


# Add "scripts" to sys.path #
scripts_path = os.path.join(os.path.dirname(os.path.abspath(__file__)))

# Imports jbonet's module #
from SBI.structure          import PDB
from SBI.structure.chain    import Chain
from SBI.structure.chain    import ChainOfProtein
from SBI.data               import aminoacids1to3
from SBI.external.blast     import blast_parser


sys.path.append(scripts_path)
import blast, contacts, dssp, interface, model_protein, model_dna, x3dna, triads, functions, fimo, tomtom, JSONfile
import pwm_pbm as PWM
import profiler as PROFILER



#-------------#
#  Functions  #
#-------------#

def update_profiles(json_dict,output_dir,remote):

    python = os.path.join(config.get("Paths", "python_path"), "python")
    pbm_dir = config.get("Paths", "pbm_dir")
    pdb_dir = config.get("Paths", "pdb_dir")
    dummy_dir = os.path.join(output_dir, "dummy")

    if "input_pwms" in json_dict.keys(): input_pwms_data = handle_input_pwms(json_dict["input_pwms"])


    redofimo=False
    if len(input_pwms_data)>0: 
        print("MP2D UPDATE: Read USER PWMs")
        redofimo=True
    if len(json_dict["dna"]["seq"])>0:
       dna_fasta_file=os.path.join(output_dir, "DNA","input_dna.fa")
       print("Make input_dna.fa ")
       dna_fasta=open(dna_fasta_file,"w")
       for dna in json_dict["dna"]["seq"]:
           dna_id = dna["_id"].lstrip(">").split()[0].replace("|","_")
           dna_seq = dna["full_sequence"]
           dna_fasta.write(">%s\n%s\n"%(dna_id,dna_seq))
           print(">%s\n%s\n"%(dna_id,dna_seq))
       dna_fasta.close()
       parameters= " -i "+os.path.join(output_dir, "TF_modeling") + " --dummy=" + dummy_dir + " --pdb " + pdb_dir + " --pbm " + pbm_dir + " -d " + dna_fasta_file
       parameters= parameters + " -v  --auto --parallel "
       if redofimo: parameters= parameters + " --redofimo "
       else:        parameters= parameters + " --reuse "
       parameters= parameters + " --html --html_types normal,energy,fimo_log_score,fimo_binding --html_energies s3dc_dd "
       parameters= parameters + " --complete 0.95 "
       parameters= parameters + " -o profile "
       if json_dict['input']['type'] == 'PDB': parameters= parameters + " --model_accuracy "
       logfile = os.path.join(output_dir,"profiler.log")
       command = "profiler.py"
    else:
       parameters= " -i "+ os.path.join(output_dir, "TF_modeling") + " --dummy=" + dummy_dir + " --pdb " + pdb_dir + " --pbm " + pbm_dir + " -a -v " + " -o " + os.path.join(output_dir, "TF_modeling") + " --reuse --parallel "
       parameters= parameters + " --complete 0.95 "
       logfile = os.path.join(output_dir,"pwm.log")
       command="pwm_pbm.py"

    # Execute PWM in remote#
    if remote:
          print("RUN %s ON REMOTE for %s"%(command,os.path.basename(output_dir)))
          job = functions.execute_in_remote(command,parameters,config,output_dir,logfile,waiting=True)
          print("FINISHED JOB %s for %s"%(str(job),os.path.basename(output_dir)))
          if str(job).startswith("Error"):
               json_dict["error_msg"] = job
               raise ValueError(str(job)) 
    else:
    # Execute PWM in local#
          print("RUN %s ON LOCAL \n"%(command))
          print("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))
          os.system("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))


def get_binding_site_length(output_dir,pdb_file,remote=False):


    dummy_dir = os.path.join(output_dir,"dummy")
    interface_file = os.path.join(dummy_dir, os.path.basename(pdb_file).replace(".pdb","_interface.txt"))
    pdb_dir = config.get("Paths", "pdb_dir")

    #MAKE INTERFACE
    parameters = " -i " + pdb_file + " --dummy " + dummy_dir + " -d dinucleotides " + " -o " + interface_file
    command="interface.py"
    if not os.path.exists(interface_file):
    # Execute INTERFACE in remote#
      if remote:
          print("RUN %s ON REMOTE"%command)
          job = functions.execute_in_remote(command,parameters,config,output_dir,logfile,waiting=True)
          print("JOB: "+str(job))
      else:
    # Execute INTERFACE in local#
          print("RUN %s ON LOCAL \n"%(command))
          print("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))
          os.system("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))

    bs_length = 0
    binding_site = []
    interface = open(interface_file,"r")
    for line in interface:
        if line.startswith("#"): continue
        bs_length = bs_length + 1
        data=line.split(";")
        if int(data[1])>0: binding_site.append(int(data[0]))
    interface.close()

    print("INTERFACE of %s"%os.path.basename(pdb_file))
    print binding_site

    start = min(binding_site)
    end   = max(binding_site)
    binding_length = end - start + 1
    if binding_length <=0: binding_length = bs_length


    return binding_length





def generate_motif_list(pdb_files,families,output_dir,json_dict,dna_fa):

    motif_list=[]
    for i in range(0, len(pdb_files)):
        current_motif = {}
        pdb_path = pdb_files[i]
        pdb_name = pdb_path.split("/")[-1].replace(".pdb", "")
        current_motif["model"] = pdb_name
        current_motif["pdb_file"] = pdb_path.replace(root_path, "")
        current_motif["seq_identity"] = get_seq_identity(json_dict, pdb_path)
        
        if json_dict['input']['type'] == 'PDB':
            current_motif["template"] = "N/A"
            #current_motif["domain"] = "N/A"
            current_motif["domain"] = pdb_name.split(":")[1] + ":" + pdb_name.split(":")[2].split("_")[0]
        else:
            current_motif["template"] = pdb_name.split("_")[1] + "_" + pdb_name.split("_")[2]
            if pdb_name.count(":") == 2:
                current_motif["domain"] = pdb_name.split(":")[1] + ":" + pdb_name.split(":")[2].split("_")[0]
            elif pdb_name.count(":") == 4:
                current_motif["domain"] = pdb_name.split(":")[1] + ":" + pdb_name.split(":")[2] + "," + pdb_name.split(":")[3] + ":" + pdb_name.split(":")[4].split("_")[0]

        current_motif["ID"] = i+1
        current_motif["TF_name"] = json_dict["TF_name"]  
        # Get PDB object #
        pdb_obj = PDB(pdb_path)
        current_motif["mono_dim_hetero"] = get_mono_dim_hetero(json_dict, pdb_obj)
        try:
                current_motif["family"] = families[current_motif["template"]]
        except:
                current_motif["family"] = "Unknown"
        # Get a MSA object for this motif #  
        motif_dir = os.path.join(output_dir, pdb_name)
        if not os.path.exists(os.path.join(motif_dir, "PWM", pdb_name + ".meme")): continue
        msa_obj_s3dc_dd = PWM.nMSA(os.path.join(motif_dir, "PWM", pdb_name + ".meme"),pdb_name,"meme")
        msa_file=os.path.join(motif_dir, "PWM", pdb_name + ".msa")
        pwm_meme=os.path.join(motif_dir, "PWM", pdb_name + ".meme")
        pwm_file=os.path.join(motif_dir, "PWM", pdb_name + ".pwm")
        logo_file=os.path.join(motif_dir, "PWM", pdb_name + ".logo")
        current_motif["binding_site_length"] = get_binding_site_length(output_dir,pdb_path,remote=False)
        #current_motif["binding_site_length"] = msa_obj_s3dc_dd.get_binding_site_length()
        current_motif["complete_meme_pwm"] = pwm_meme
        current_motif["msa_file"] = msa_file.replace(root_path, "")
        current_motif["meme_pwm"] = pwm_meme.replace(root_path, "")
        current_motif["pwm_file"] = pwm_file.replace(root_path, "")
        current_motif["logos"] = [logo_file.replace(root_path, "") + ".fwd.png", logo_file.replace(root_path, "") + ".rev.png"] 
        try:
          if msa_obj_s3dc_dd.information_content()>0:
            current_motif["IC"] = msa_obj_s3dc_dd.information_content() 
            current_motif["information_content"] = round(msa_obj_s3dc_dd.information_content(), 3)
          else:
            current_motif["IC"] = 0
            current_motif["information_content"] =0
        except:
            current_motif["IC"] = 0
            current_motif["information_content"] =0
        current_motif["consensus_sequence"] = msa_obj_s3dc_dd.get_main_sequence()
        #Frequencies to view on website
        #current_motif["pwm_frequencies"]=(numpy.array(msa_obj_s3dc_dd.get_pwm())).astype(numpy.float).transpose()
        current_motif["pwm_frequencies"]=msa_obj_s3dc_dd.get_pwm()
        #Protein binding profile
        current_motif["binding_profiles"]={}
        for dna in json_dict["dna"]["seq"]:
          dna_id  = dna["_id"]
          current_motif["binding_profiles"].setdefault(dna_id,{})
          html_file = os.path.join(output_dir, pdb_name, "DNA", "binding_profiles", pdb_name + "." + dna_fa[dna_id] + ".html")
          csv_file  = os.path.join(output_dir, pdb_name, "DNA", "binding_profiles", pdb_name + "." + dna_fa[dna_id] + ".csv")
          current_motif["binding_profiles"][dna_id]={}
          current_motif["binding_profiles"][dna_id].setdefault("html",html_file.replace(root_path, ""))
          current_motif["binding_profiles"][dna_id].setdefault("csv",csv_file.replace(root_path, ""))
          current_motif["binding_profiles"][dna_id].setdefault("html_fullscreen",html_file.replace(root_path, ""))
        #Protein binding profile comparison
        done_comparison=set()
        for dna1 in json_dict["dna"]["seq"]:
         for dna2 in json_dict["dna"]["seq"]:
          dna1_id = dna1["_id"]
          dna2_id = dna2["_id"]
          if dna1_id == dna2_id: continue
          if (dna1_id,dna2_id) in done_comparison: continue
          html_file = os.path.join(output_dir, pdb_name, "DNA", "binding_profiles", pdb_name + "." + dna_fa[dna1_id] + "." + dna_fa[dna2_id] + ".compare.html")
          current_motif["binding_profiles"].setdefault("comparison", {})
          current_motif["binding_profiles"]["comparison"].setdefault(dna1_id + "_" + dna2_id,{})
          current_motif["binding_profiles"]["comparison"][dna1_id + "_" + dna2_id].setdefault("html",html_file.replace(root_path, ""))
          current_motif["binding_profiles"]["comparison"][dna1_id + "_" + dna2_id].setdefault("html_fullscreen",html_file.replace(root_path, ""))
          done_comparison.add((dna1_id,dna2_id))
        # Set the tf_model from this iteration for the next one, to avoid reloading potentials if possible #
        motif_list.append(current_motif)

    # Motifs order is defined by %ID with templates #
    motif_list=sorted(motif_list,key=lambda motif_list: motif_list["seq_identity"],reverse=True)


    return motif_list

def update_classified_pdbs(pdb_files):

    directories = set()
    for pdb_file in pdb_files:
          print("-- UPDATE_CLASS_PDB CHECK PDBs %s"%pdb_file)
          models_dir=os.path.dirname(pdb_file)
          directories.update(set([os.path.join(models_dir,d) for d in os.listdir(models_dir) if os.path.isdir(os.path.join(models_dir,d)) and d.startswith("profile")]))

    number_of_clusters=0
    cluster={}
    for folder in directories:
          print("-- UPDATE_CLASS_PDB CHECK FOLDER %s"%folder)
          profile_name=os.path.basename(folder)
          print("-- UPDATE_CLASS_PDB PROFILE %s"%profile_name)
          number_of_clusters= number_of_clusters + 1
          for meme_file in os.listdir(folder):
              if not meme_file.endswith(".meme"): continue
              for pdb_name in pdb_files:
                  if str(os.path.basename(pdb_name.rstrip(".pdb"))) in str(meme_file.rstrip(".meme")):
                          pdb_file = os.path.basename(pdb_name.rstrip(".pdb"))
                          print("-- UPDATE_CLASS_PDB CHECK FILE %s %s"%(meme_file,pdb_file))

              #pdb_file = [ pdb_name for pdb_name in pdb_files if str(pdb_name.rstrip(".pdb")) in str(meme_file.rstrip(".meme"))][0]
              words=meme_file.lstrip("model:").split(":")
              positions=[]
              for word in words:
                  print("-- UPDATE_CLASS_PDB CHECK POSITIONS %s"%word.split("_")[0])
                  positions.append(word.split("_")[0])
              a  = int(positions[0])
              b  = int(positions[1])
              c  = d = 0
              print("-- UPDATE_CLASS_PDB ADD TO CLUSTER %d %s %d %d %d %d"%(number_of_clusters,os.path.join(folder,pdb_file),int(a),int(b),int(c),int(d)))
              cluster.setdefault(number_of_clusters,set()).add((os.path.join(folder,pdb_file),int(a),int(b),int(c),int(d)))
    cluster_domain={}
    for number_of_cluster in cluster.keys():
      list_of_data = cluster.get(number_of_cluster)
      x=np.array([ int(a) for name,a,b,c,d in  list_of_data])
      y=np.array([ int(b) for name,a,b,c,d in  list_of_data])
      cluster_domain.setdefault(number_of_cluster,(x.mean(),y.mean()))
    order = [nc for nc,position in sorted(cluster_domain.items(),key=lambda x: x[1][0])]

    cluster_ordered={}

    number_of_cluster = 0
    for nc in order:
        number_of_cluster = number_of_cluster + 1
        cluster_ordered.setdefault(number_of_cluster,cluster.get(nc))

    return cluster_ordered

 


def get_cluster_fragments(clusters,fragment_list):

    groups_by_fragment={}
    if len(fragment_list)>0:
      for interval in fragment_list:
        if len(interval.split("."))>1:
              frg       =[int(z) for z in interval.split(".")[0].split("_")]
              set_frg   = set([z for z in range(frg[0],frg[1])])
              try:
                frg_dimer =[int(z) for z in interval.split(".")[1].split("_")]
                set_frg2  = set([z for z in range(frg_dimer[0],frg_dimer[1])])
              except Exception as e:
                print("INTERVAL: %s yields error reading %s"%(interval,e))
                set_frg2=None
              groups_by_fragment.setdefault(interval,(set_frg,set_frg2))
        else:
              frg       =[int(z) for z in interval.split("_")]
              set_frg   = set([z for z in range(frg[0],frg[1])])
              groups_by_fragment.setdefault(interval,(set_frg,None))

    for interval in fragment_list:
        print("MP2D: GET_CLUSTER_FRAGMENT  LIST "+interval)
        print groups_by_fragment[interval]

    groups_by_cluster={}
    for number_of_cluster in clusters.keys():
        list_of_data = clusters.get(number_of_cluster)
        x=np.array([ int(a) for name,a,b,c,d in  list_of_data])
        y=np.array([ int(b) for name,a,b,c,d in  list_of_data])
        u=np.array([ int(c) for name,a,b,c,d in  list_of_data])
        v=np.array([ int(d) for name,a,b,c,d in  list_of_data])
        set_frg = set([z for z in range(int(x.mean()),int(y.mean()))])
        if int(u.mean()) > 0 or int(v.mean()) > 0: 
           set_frg2 = set([z for z in range(int(u.mean()),int(v.mean()))])
           groups_by_cluster.setdefault(number_of_cluster,(set_frg,set_frg2))
        else:
           groups_by_cluster.setdefault(number_of_cluster,(set_frg,None))

    for clu in groups_by_cluster.keys():
        print("MP2D: GET_CLUSTER_FRAGMENT  CLUSTER "+str(clu))
        print groups_by_cluster[clu]

    fragment={}
    if len(fragment_list)>0:
      for number_of_cluster in clusters.keys():
        data_set=groups_by_cluster.get(number_of_cluster)
        if data_set[1] is not None:
           selected_fragment=None
           for interval in groups_by_fragment.keys():
             if selected_fragment is None: selected_fragment=interval
             else:
               if interval[1] is not None and selected_fragment[1] is not None: 
                  total     =len(data_set[0].intersection(groups_by_fragment[interval][0]))+len(data_set[1].intersection(groups_by_fragment[interval][1]))
                  previous  =len(data_set[0].intersection(groups_by_fragment[selected_fragment][0]))+len(data_set[1].intersection(groups_by_fragment[selected_fragment][1]))
               else: 
                  total     =len(data_set[0].intersection(groups_by_fragment[interval][0]))
                  previous  =len(data_set[0].intersection(groups_by_fragment[selected_fragment][0]))
               if total > previous : selected_fragment=interval
        else:
           selected_fragment=None
           for interval in groups_by_fragment.keys():
             if selected_fragment is None: selected_fragment=interval
             else:
               total     =len(data_set[0].intersection(groups_by_fragment[interval][0]))
               previous  =len(data_set[0].intersection(groups_by_fragment[selected_fragment][0]))
               if total > previous : selected_fragment=interval
        fragment.setdefault(number_of_cluster,selected_fragment)
    else:
      for number_of_cluster in clusters.keys():
        data_set=groups_by_cluster.get(number_of_cluster)
        a=min(data_set[0])
        b=1+max(data_set[0])
        profile_label_dir  = str(a)+"_"+str(b)
        if data_set[1] is not None:
          c=min(data_set[1])
          d=1+max(data_set[1])
          profile_label_dir = profile_label_dir + "." + str(c) + "_" + str(d)
        fragment.setdefault(number_of_cluster,profile_label_dir)

    for clu in fragment.keys():
        print("MP2D: GET_CLUSTER_FRAGMENT  FRAGMENT "+str(clu)+"  ....  "+fragment[clu])    

    return fragment

def group_motifs_by_clusters(clusters,motif_list):

    # Cluster motifs into clusters. 
    groups={}
    for number_of_cluster in clusters.iterkeys():
        list_of_data = clusters.get(number_of_cluster)
        x=np.array([ int(a) for name,a,b,c,d in  list_of_data])
        y=np.array([ int(b) for name,a,b,c,d in  list_of_data])
        u=np.array([ int(c) for name,a,b,c,d in  list_of_data])
        v=np.array([ int(d) for name,a,b,c,d in  list_of_data])
        profile_label_dir  = str(int(x.mean()))+"_"+str(int(y.mean()))
        if int(u.mean()) > 0 or int(v.mean()) > 0: profile_label_dir = profile_label_dir + "." + str(int(c)) + "_" + str(int(d))
        #print("Cluster %s Label %s"%(str(number_of_cluster),profile_label_dir))
        for mot in motif_list:
           #print("Motif %s"%str(mot["model"]))
           for name,a,b,c,d in clusters.get(number_of_cluster):
              name_pdb=os.path.basename(name).rstrip(".pdb")
              #print("NAME %s"%str(name_pdb))
              if  mot["model"]==name_pdb:
                mot['ID']=len(groups.setdefault(number_of_cluster,[]))+1
                groups.setdefault(number_of_cluster,[]).append(mot)
    return groups

def get_motif_info(motif_list):

    motif_info = {}
    for ID in range(1, len(motif_list)+1):
        clu = motif_list[ID-1]
        max_x = set()
        min_x = set()
        max_y = set()
        min_y = set()
        x_values = []
        y_values = []

        for mot in clu:
            # We are working with monomers #
            if len(mot["domain"].split(":")) == 2:
                domain = range(int(mot["domain"].split(":")[0]), int(mot["domain"].split(":")[1]))
                xvals = set()
                for d in domain:
                    xvals.add(d)
                x_values = list(xvals)
                y_values = [0]*len(x_values)
            # We are working with dimers #
            if len(mot["domain"].split(":")) == 3:
                dom1 = mot["domain"].split(",")[0]
                dom2 = mot["domain"].split(",")[1]
                max_x.add(max(range(int(dom1.split(":")[0]), int(dom1.split(":")[1]))))
                min_x.add(min(range(int(dom1.split(":")[0]), int(dom1.split(":")[1]))))
                max_y.add(max(range(int(dom2.split(":")[0]), int(dom2.split(":")[1]))))
                min_y.add(min(range(int(dom2.split(":")[0]), int(dom2.split(":")[1]))))
                x_values = range(min(range(int(dom1.split(":")[0]), int(dom1.split(":")[1]))), max(range(int(dom1.split(":")[0]), int(dom1.split(":")[1]))), 1)
                y_values = [((max(range(int(dom2.split(":")[0]), int(dom2.split(":")[1])))+min(range(int(dom2.split(":")[0]), int(dom2.split(":")[1]))))/2.0)]*len(x_values)

        motif_info[ID] = {}
        if max_x == set():
            # The cluster is only made of monomers #
            motif_info[ID]["x_values"] = x_values
            motif_info[ID]["y_values"] = y_values
        else:
            # The cluster contains dimers #
            motif_info[ID]["max_x"] = max(list(max_x))
            motif_info[ID]["min_x"] = min(list(min_x))
            motif_info[ID]["max_y"] = max(list(max_y))
            motif_info[ID]["min_y"] = min(list(min_y))
            motif_info[ID]["x_values"] = x_values
            motif_info[ID]["y_values"] = y_values

    return motif_info


def  generate_binding_profiles_json(clusters,fragment,output_dir,json_dict,dna_fa):

    jd={}
    for number_of_cluster in clusters.iterkeys():
      jd.setdefault(str(number_of_cluster),{})
      for dna in json_dict["dna"]["seq"]:
          dna_id  = dna["_id"]
          html_file = os.path.join(output_dir,"DNA","binding_profiles","profile."+dna_fa[dna_id]+".average."+fragment[number_of_cluster]+".html")
          csv_file  = os.path.join(output_dir,"DNA","binding_profiles","profile."+dna_fa[dna_id]+".average."+fragment[number_of_cluster]+".csv")
          jd[str(number_of_cluster)].setdefault(dna_id,{})
          jd[str(number_of_cluster)][dna_id].setdefault("html",html_file.replace(root_path, ""))
          jd[str(number_of_cluster)][dna_id].setdefault("html_fullscreen",html_file.replace(root_path, ""))
          jd[str(number_of_cluster)][dna_id].setdefault("csv",csv_file.replace(root_path, ""))
      done_comparison=set()
      for dna1 in json_dict["dna"]["seq"]:
        for dna2 in json_dict["dna"]["seq"]:
          dna1_id = dna1["_id"]
          dna2_id = dna2["_id"]
          if dna1_id == dna2_id: continue
          if (dna1_id,dna2_id) in done_comparison: continue
          html_file = os.path.join(output_dir,"DNA","binding_profiles","profile."+dna_fa[dna1_id]+"."+dna_fa[dna2_id]+".compare."+fragment[number_of_cluster]+".html")
          jd[str(number_of_cluster)].setdefault("comparison",{})
          jd[str(number_of_cluster)]["comparison"].setdefault(dna1_id + "_" + dna2_id,{})
          jd[str(number_of_cluster)]["comparison"][dna1_id + "_" + dna2_id].setdefault("html",html_file.replace(root_path, ""))
          jd[str(number_of_cluster)]["comparison"][dna1_id + "_" + dna2_id].setdefault("html_fullscreen",html_file.replace(root_path, ""))
          done_comparison.add((dna1_id,dna2_id))
    return jd


def generate_profiles(json_dict,output_dir,remote):

    python = os.path.join(config.get("Paths", "python_path"), "python")
    pbm_dir = config.get("Paths", "pbm_dir")
    pdb_dir = config.get("Paths", "pdb_dir")
    dummy_dir = os.path.join(output_dir, "dummy")

    if len(json_dict["dna"]["seq"])>0:
       dna_fasta_file=os.path.join(output_dir, "DNA","input_dna.fa")
       dna_fasta=open(dna_fasta_file,"w")
       for dna in json_dict["dna"]["seq"]:
           dna_id = dna["_id"].lstrip(">").split()[0]
           dna_seq = dna["full_sequence"]
           dna_fasta.write(">%s\n%s\n"%(dna_id,dna_seq))
           print(">%s\n%s\n"%(dna_id,dna_seq))
       dna_fasta.close()
       parameters= " -i "+os.path.join(output_dir, "TF_modeling") + " --dummy=" + dummy_dir + " --pdb " + pdb_dir + " --pbm " + pbm_dir + " -d " + dna_fasta_file
       parameters= parameters + " -v --chains_fixed --auto --parallel "
       parameters= parameters + " --html --html_types normal,energy,fimo_log_score,fimo_binding --html_energies s3dc_dd "
       parameters= parameters + " --complete 0.95 "
       parameters= parameters + " -o profile "
       if json_dict['input']['type'] == 'PDB': parameters= parameters + " --model_accuracy "
       logfile = os.path.join(output_dir,"profiler.log")
       command = "profiler.py"
    else:
       parameters= " -i "+ os.path.join(output_dir, "TF_modeling") + " --dummy=" + dummy_dir + " --pdb " + pdb_dir + " --pbm " + pbm_dir + " -a -v " + " -o " + os.path.join(output_dir, "TF_modeling") + " --reuse --parallel "
       parameters= parameters + " --complete 0.95 "
       logfile = os.path.join(output_dir,"pwm.log")
       command="pwm_pbm.py"

    # Execute PWM in remote#
    if remote:
          print("RUN %s ON REMOTE"%command)
          job = functions.execute_in_remote(command,parameters,config,output_dir,logfile,waiting=True)
          print("FINISHED REMOTE JOB: "+job)
          if str(job).startswith("Error"):
               json_dict["error_msg"] = job
               raise ValueError(str(job)) 
    else:
    # Execute PWM in local#
          print("RUN %s ON LOCAL \n"%(command))
          print("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))
          os.system("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))


def arrange_profiles_for_web(json_dict,output_dir,pdb_files):
  
    fragment_list=[]
    fragment_set =set()

    if len(json_dict["dna"]["seq"])>0:
      n_dna=0
      dna_names=[]
      for dna in json_dict["dna"]["seq"]:
          n_dna=n_dna+1
          dna_id = dna["_id"].lstrip(">").split()[0]
          dna_seq = dna["full_sequence"]
          dna_names.append(dna_id+"_"+str(n_dna))
      for pdb_file in pdb_files:
        models_dir=os.path.dirname(pdb_file)
        base_name=pdb_file.rstrip(".pdb")
        pdb_name=os.path.basename(base_name)
        pwm_file = pdb_name+".pwm"
        meme_file= pdb_name+".meme"
        msa_file = pdb_name+".msa"
        logo_file= pdb_name+".logo"
        logo_fwd = logo_file+".fwd.png"
        logo_rev = logo_file+".rev.png"
        print("\t-- Copy profiles of %s"%os.path.basename(base_name))
        # Create output directories #
        motif_dir = os.path.join(output_dir, pdb_name)
        if os.path.exists(motif_dir) == False:
            os.mkdir(motif_dir)
        if os.path.exists(os.path.join(motif_dir,"PWM")) == False:
            os.mkdir(os.path.join(motif_dir, "PWM"))
        if os.path.exists(os.path.join(motif_dir, "DNA")) == False:
            os.mkdir(os.path.join(motif_dir, "DNA"))    
        if os.path.exists(os.path.join(motif_dir, "DNA", "binding_profiles")) == False:
            os.mkdir(os.path.join(motif_dir, "DNA", "binding_profiles"))
        if os.path.exists(os.path.join(motif_dir, "DNA", "DNA_modeling")) == False:
            os.mkdir(os.path.join(motif_dir, "DNA", "DNA_modeling"))
        # Make accessible to evey user the folders
        os.system("chmod -R 777 %s"%motif_dir)
        # Hard-links to the original PWMs 
        directories = [d for d in os.listdir(models_dir) if os.path.isdir(os.path.join(models_dir,d))]
        for folder in directories:
         for dna_name in dna_names:
          pwm_dna_file = pdb_name+"."+dna_name+".pwm"
          meme_dna_file= pdb_name+"."+dna_name+".meme"
          msa_dna_file = pdb_name+"."+dna_name+".msa"
          logo_dna_file= pdb_name+"."+dna_name+".logo"
          logo_dna_fwd = logo_dna_file+".fwd.png"
          logo_dna_rev = logo_dna_file+".rev.png"
          try:
            if pwm_dna_file  in os.listdir(os.path.join(models_dir,folder)): shutil.copy(os.path.join(models_dir,folder,pwm_dna_file),os.path.join(motif_dir, "PWM",pwm_file))
            if meme_dna_file in os.listdir(os.path.join(models_dir,folder)): shutil.copy(os.path.join(models_dir,folder,meme_dna_file),os.path.join(motif_dir, "PWM",meme_file))
            if msa_dna_file  in os.listdir(os.path.join(models_dir,folder)): shutil.copy(os.path.join(models_dir,folder,msa_dna_file),os.path.join(motif_dir, "PWM",msa_file))
            if logo_dna_fwd  in os.listdir(os.path.join(models_dir,folder)): shutil.copy(os.path.join(models_dir,folder,logo_dna_fwd),os.path.join(motif_dir, "PWM",logo_fwd))
            if logo_dna_rev  in os.listdir(os.path.join(models_dir,folder)): shutil.copy(os.path.join(models_dir,folder,logo_dna_rev),os.path.join(motif_dir, "PWM",logo_rev))
          except Exception as e:
            print("Failed to link PWM of %s (Error %s)"%(pdb_name,e))
            continue
         for html_file in  os.listdir(os.path.join(models_dir,folder)):
           if not html_file.endswith(".html"): continue
           html_file2=html_file
           if pdb_name in html_file:
              shutil.copy(os.path.join(models_dir,folder,html_file),os.path.join(output_dir, pdb_name, "DNA", "binding_profiles",html_file2))
           if html_file.startswith("profile"):
              fragment=folder.lstrip("profile.")
              fragment_set.add(fragment)
              shutil.copy(os.path.join(models_dir,folder,html_file),os.path.join(output_dir, "DNA", "binding_profiles",html_file2.rstrip("html")+fragment+".html"))
         for csv_file in  os.listdir(os.path.join(models_dir,folder)):
           if not csv_file.endswith(".csv"): continue
           csv_file2=csv_file
           if pdb_name in csv_file:
              shutil.copy(os.path.join(models_dir,folder,csv_file),os.path.join(output_dir, pdb_name, "DNA", "binding_profiles",csv_file2))
           if csv_file.startswith("profile"):
              fragment=folder.lstrip("profile.")
              shutil.copy(os.path.join(models_dir,folder,csv_file),os.path.join(output_dir, "DNA", "binding_profiles",csv_file2.rstrip("csv").replace("mean","average")+fragment+".csv"))
    else:
      for pdb_file in pdb_files:
        base_name=pdb_file.rstrip(".pdb")
        pwm_file = base_name+".pwm"
        meme_file= base_name+".meme"
        msa_file = base_name+".msa"
        logo_file= base_name+".logo"
        logo_fwd = logo_file+".fwd.png"
        logo_rev = logo_file+".rev.png"
        pdb_name=os.path.basename(base_name)
        print("\t-- Copy profiles of %s"%os.path.basename(base_name))
        # Create output directories #
        motif_dir = os.path.join(output_dir, pdb_name)
        if os.path.exists(motif_dir) == False:
            os.mkdir(motif_dir)
        if os.path.exists(os.path.join(motif_dir,"PWM")) == False:
            os.mkdir(os.path.join(motif_dir, "PWM"))
        if os.path.exists(os.path.join(motif_dir, "DNA")) == False:
            os.mkdir(os.path.join(motif_dir, "DNA"))    
        if os.path.exists(os.path.join(motif_dir, "DNA", "binding_profiles")) == False:
            os.mkdir(os.path.join(motif_dir, "DNA", "binding_profiles"))
        if os.path.exists(os.path.join(motif_dir, "DNA", "DNA_modeling")) == False:
            os.mkdir(os.path.join(motif_dir, "DNA", "DNA_modeling"))
        # Make accessible to evey user the folders
        os.system("chmod -R 777 %s"%motif_dir)
        # Hard-links to the original PWMs 
        
        try:
           shutil.copy(pwm_file,os.path.join(motif_dir, "PWM",os.path.basename(pwm_file)))
           shutil.copy(meme_file,os.path.join(motif_dir, "PWM",os.path.basename(meme_file)))
           shutil.copy(msa_file,os.path.join(motif_dir, "PWM",os.path.basename(msa_file)))
           shutil.copy(logo_fwd,os.path.join(motif_dir, "PWM",os.path.basename(logo_fwd)))
           shutil.copy(logo_rev,os.path.join(motif_dir, "PWM",os.path.basename(logo_rev)))
        except Exception as e:
           print("Failed to link PWM of %s (Error %s)"%(pdb_name,e))
           continue

    os.system("chmod -R 777 %s"%output_dir)

    fragment_list=[frg for frg in fragment_set]

    return fragment_list


def tf_modeling(json_dict, input_file, output_dir,  remote,dummy_dir):

    pbm_dir = config.get("Paths", "pbm_dir")
    pdb_dir = config.get("Paths", "pdb_dir")

    # Check the number of input fasta sequences in order to set if the output should only contain dimers #
    seqs = 0
    for line in open(input_file, "r").readlines():
        if line.startswith(">"):
            seqs += 1
    if seqs == 2:
        json_dict["monomers"] = "false"
        json_dict["dimers"] = "true"
    # Define the options for the modeling #
    tags = " -f  --renumerate -v --twilight --restrictive "
    # Select the number of templates #
    if (json_dict["templates"] == "1") and (len(json_dict["input"]["_id"].split(",")) == 1):
        tags += "--best "
    elif json_dict["templates"] == "all":
        tags += "-a --unrestrictive "
    else:
        tags += "--n-temp=" + str(json_dict["templates"]) + " "
    # Select the number of models per template #
    tags += "--n-model=" + str(json_dict["models"]) + " "
    # Select whether to model monomers or dimers or both #
    if (json_dict["monomers"] == "true") and (json_dict["dimers"] == "false"):
        tags += "--monomers "
    elif (json_dict["monomers"] == "false") and (json_dict["dimers"] == "true"):
        tags += "--dimers "

    parameters= " -i " + input_file + " -o " + os.path.join(output_dir, "TF_modeling") + " -l model -p " + pdb_dir + " --dummy=" + dummy_dir + " " + tags
    logfile = os.path.join(output_dir,"modeling.log")
    command="model_protein.py"
    # Execute modeling in remote#
    if remote:
          print("RUN %s ON REMOTE"%command)
          job = functions.execute_in_remote(command,parameters,config,output_dir,logfile,waiting=True)
          print("FINISHED REMOTE JOB: "+job)
          if str(job).startswith("Error"):
               json_dict["error_msg"] = job
               return []
    else:
    # Execute modeling in local#
          print("RUN %s ON LOCAL \n"%(command))
          os.system("%s %s/%s %s \n"%(python,scripts_path,command,parameters))

    # Get the output files # 
    pdb_files = get_pdb_models(output_dir)
    os.system("chmod -R 777 %s"%output_dir)

    return pdb_files

def get_pdb_models(output_dir):

    # Get the output files # 
    pdb_files = []
    for f in os.listdir(os.path.join(output_dir, "TF_modeling")):
        if (f.startswith("model")) and (f.endswith(".pdb")):
           pdb_files.append(os.path.join(output_dir, "TF_modeling", f))

    return sorted(pdb_files)


def tf_threading(input_file, output_dir, pdb_dir, dummy_dir):

    print(python + " " + scripts_path + "/model_protein.py -i " + input_file + " -o " + os.path.join(output_dir, "TF_modeling") + " -p " + pdb_dir + " --dummy=" + dummy_dir + " -t -v -l model")
    process = subprocess.Popen([python, scripts_path + "/model_protein.py", "-i", input_file, "-o", os.path.join(output_dir, "TF_modeling"), "-p " + pdb_dir, "--dummy=" + dummy_dir, "-t", "-v", "-l", "model"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    pdb_list = []
    for file in os.listdir(os.path.join(output_dir, "TF_modeling")):
        if (file.startswith("model")) and (file.endswith(".pdb")):
            pdb_list.append(os.path.join(os.path.join(output_dir, "TF_modeling"), file))

    return pdb_list



def get_seq_identity(json_dict, pdb_path):

    seq_identity = None
    if json_dict['input']['type'] == "PDB":
        return 100
    else:

        model_name = os.path.basename(pdb_path)
        template = model_name.split("_")[1]
        chain = model_name.split("_")[2]
        num = model_name.split("_")[3].replace(".pdb", "")
        model_dir = "/".join(pdb_path.split("/")[:-1])
        model_summary = open(os.path.join(model_dir, "model_model.summary.txt"), "r").readlines()[1:]
        for line in model_summary:
            line = line.split(";")
            if (str(line[1]) == str(num)) and (str(line[2]) == str(template)):
                if not "," in str(line[5]):
                    # For monomers #
                    if line[5] == chain:
                        seq_identity = float(line[-2])
                else:
                    # For dimers #
                    if chain in line[5]:
                        seq_identity = float(line[-2])

        return seq_identity


def get_mono_dim_hetero(json_dict, pdb_obj):

    protein_chains = []
    for p in pdb_obj.chains:
        if isinstance(p, ChainOfProtein):
            protein_chains.append(p)

    if len(protein_chains) == 1:
        return "Monomer"
    elif len(protein_chains) == 2:
        return "Dimer"
            

def tf_dna_modeling(protein, dna, remote=False):

    # Obtain the interface file #
    #Required files and directories
    dna = dna.upper()
    if "TF_modeling" in protein:
        main_dir = protein.split("/")[-3]
        output_dir = os.path.join(root_path,main_dir)
        input_protein = os.path.join(root_path,main_dir,"TF_modeling",os.path.basename(protein))
    else:
        main_dir = protein.split("/")[-2]
        output_dir = os.path.join(root_path,main_dir)
        input_protein = os.path.join(root_path,main_dir,os.path.basename(protein))

    dummy_dir = os.path.join(output_dir,"dummy")
    interface_file = os.path.join(dummy_dir, os.path.basename(protein) + "_" + dna + "_interface.txt")
    pdb_dir = config.get("Paths", "pdb_dir")
    print("MODEL protein %s with DNA %s"%(input_protein,dna.upper()))

    #MAKE INTERFACE
    parameters = " -i " + input_protein + " --dummy " + dummy_dir + " -d dinucleotides " + " -o " + interface_file
    command="interface.py"
    if not os.path.exists(interface_file):
    # Execute INTERFACE in remote#
      if remote:
          print("RUN %s ON REMOTE"%command)
          job = functions.execute_in_remote(command,parameters,config,output_dir,logfile,waiting=True)
          print("JOB: "+job)
      else:
    # Execute INTERFACE in local#
          print("RUN %s ON LOCAL \n"%(command))
          print("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))
          os.system("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))


    #MAKE MODEL
    output_file = os.path.join(output_dir, "model." + dna + ".pdb")
    parameters = " -i " + interface_file + " --dummy " + dummy_dir + " -p " + input_protein+" -o " + output_dir + " -s " + dna + " --pdb=" + pdb_dir
    command="model_dna.py"
    if not os.path.exists(output_file):
    # Execute INTERFACE in remote#
      if remote:
          print("RUN %s ON REMOTE"%command)
          job = functions.execute_in_remote(command,parameters,config,output_dir,logfile,waiting=True)
          print("JOB: "+job)
      else:
    # Execute INTERFACE in local#
          print("RUN %s ON LOCAL \n"%(command))
          print("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))
          os.system("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))
    
    #Modify path of output file
    os.system("chmod -R 777 %s"%output_file)
    output_file = os.path.join(output_dir, "model." + dna + ".pdb").replace(root_path, "")

    return output_file

def handle_input_pwms(input_pwms):

    if input_pwms == "" or input_pwms == "null":
        return []
    input_pwms_data = []
    pwm_list = input_pwms.split(",")
    for pwm in pwm_list:
        pwm_obj = {}
        pwm_obj["pwm"] = pwm.split(";")[0]
        pwm_obj["cluster"] = pwm.split(";")[1]
        pwm_obj["motif"] = pwm.split(";")[2]
        pwm_obj["name"] = pwm.split(";")[3]


        input_pwms_data.append(pwm_obj)

    return input_pwms_data

def load_input_pwm(input_pwms_data,json_dict,output_dir):
   
    print("LOAD NEW PWMs")
    for k in range(0, len(json_dict["motifs"])):
      for i in range(0, len(json_dict["motifs"][k])):
        number_of_cluster = k+1
        number_of_motif   = i+1
        print("  -- LOAD PWM:  Cluster %d and Motif %s"%(number_of_cluster,number_of_motif))
        for pwm_data in input_pwms_data:
          print("  -- LOAD PWM:  Check INPUTED cluster %d and INPUTED motif %s"%(int(pwm_data["cluster"]),int(pwm_data["motif"])))
          if int(pwm_data["cluster"]) == number_of_cluster and int(pwm_data["motif"]) == number_of_motif:
            motif = json_dict["motifs"][k][i]
            pwm_file=os.path.join(root_path,motif["meme_pwm"].lstrip("/"))
            print("  -- LOAD PWM: MOVE PWM %s to %s "%(pwm_file,pwm_file+".old"))
            shutil.move(pwm_file,pwm_file+".old")
            print("  -- LOAD PWM:  PWM MODIFIED %s"%pwm_file)
            pm = open(pwm_file,"w")
            for line in str(pwm_data["pwm"]).split("\r"):
               pm.write("%s\n"%(line.rstrip()))
               print("  -- LOAD PWM: new line "+line.rstrip())
            pm.close()
            protein_name=os.path.basename(pwm_file).rstrip(".meme")
            print("  -- LOAD PWM: model name is "+protein_name)
            if os.path.exists(os.path.join(output_dir,"TF_modeling")):
              for folder in os.listdir(os.path.join(output_dir,"TF_modeling")):
                print("  -- LOAD PWM: check %s"%(folder))
                if protein_name in folder and folder.endswith("meme"):
                         print("  -- LOAD PWM: copy %s to %s"%(pwm_file,os.path.join(output_dir,"TF_modeling",folder)))
                         shutil.copy(pwm_file,os.path.join(output_dir,"TF_modeling",folder))
                if os.path.isdir(os.path.join(output_dir,"TF_modeling",folder)) and folder.startswith("profile"):
                   for file in os.listdir(os.path.join(output_dir,"TF_modeling",folder)):
                     print("  -- LOAD PWM: check file %s"%(file))
                     if protein_name in file and file.endswith("meme"):
                         print("  -- LOAD PWM: copy %s to %s"%(pwm_file,os.path.join(output_dir,"TF_modeling",folder,file)))
                         shutil.copy(pwm_file,os.path.join(output_dir,"TF_modeling",folder,file))
                   



def clean_up_profiles(json_dict,working_dna):
  # Remove old binding profiles #
  # For Global
  if json_dict.has_key("binding_profiles"):
      for number_of_cluster in json_dict["binding_profiles"].keys():
        for key in json_dict["binding_profiles"][number_of_cluster].keys():
            if key == "comparison": continue
            if not key in working_dna: 
                json_dict["binding_profiles"].pop(key, None)
        if json_dict["binding_profiles"][number_of_cluster].has_key("comparison"):
          for twokey in json_dict["binding_profiles"][number_of_cluster]["comparison"].keys():
            key1,key2 = twokey.split("_")
            if not key1 in working_dna or not key2 in working_dna: 
                json_dict["binding_profiles"]["comparison"].pop(twokey, None)
  # For single Motifs
  if  json_dict.has_key("motifs"):
    for k in range(0, len(json_dict["motifs"])):
     for i in range(0, len(json_dict["motifs"][k])):
       if json_dict["motifs"][k][i].has_key("binding_profiles"):
        for key in json_dict["motifs"][k][i]["binding_profiles"].keys():
            if key == "comparison": continue
            if not key in working_dna:
               json_dict["motifs"][k][i]["binding_profiles"].pop(key, None)
        if json_dict["motifs"][k][i]["binding_profiles"].has_key("comparison"):
          for twokey in json_dict["motifs"][k][i]["binding_profiles"]["comparison"].keys():
            key1,key2 = twokey.split("_")
            if not key1 in working_dna or not key2 in working_dna: 
                json_dict["motifs"][k][i]["binding_profiles"]["comparison"].pop(twokey, None)


def fix_pdb(input_file,output_file,dummy_dir):
  chain_ids=list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
  n=0
  write_file=False
  dummy_list=[]
  with open(input_file,"r") as f:
    for line in f:
        if line.startswith("ATOM") and not write_file:
            dummy_file = open(os.path.join(dummy_dir,"dummy_"+str(n)+".pdb"),"w")
            dummy_list.append(os.path.join(dummy_dir,"dummy_"+str(n)+".pdb"))
            write_file=True
        if line.startswith("END") or line.startswith("TER"):
            dummy_file.close()
            n = n + 1
            write_file=False
            continue
        if write_file:
            dummy_file.write("%s\n"%line.strip())
  dummy_file.close()
  f.close()
  chains=[]
  for d in dummy_list:
    p=PDB(d)
    for c in p.chains:
        chains.append(c)
  q=PDB()

  for i in range(len(chains)):
    c =chains[i]
    c.chain = chain_ids[i]
    q.add_chain(c)


  q.write(output_file,force=True)
  for d in dummy_list:
   os.remove(d)


#---------------#
# Options       #
#---------------#

def parse_options():
    '''
    This function parses the command line arguments and returns an optparse object.
    '''

    parser = optparse.OptionParser("Usage: modcre.py [--dummy=DUMMY_DIR] -j JSON_FILE [-o OUTPUT_DIR -v]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="DUMMY_DIR")
    parser.add_option("-j", "--json", action="store", type="string", dest="json_file", help="Json file (from jsonBuilder.py; default = None)", metavar="JSON_FILE")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="OUTPUT_DIR")
    parser.add_option("-u", "--update",default=False,action="store_true", dest="update",help="Run update mode of protein 2 dna (default=False)")
    parser.add_option("-r", "--remote",default=False,action="store_true", dest="remote",help="Run in remote host server (default=False)")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")
    
    (options, args) = parser.parse_args()

    if options.json_file == None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

#---------------#
# Protein 2 DNA #
#---------------#

def protein_2_dna(json_dict, output_dir, remote=False,verbose=True):


    # Initialize #
    input_file = os.path.join(output_dir, json_dict["input"]["file"])
    dummy_dir = os.path.join(output_dir, "dummy")
    pdb_dir = config.get("Paths", "pdb_dir")
    # Get the name of the TF we are working with #
    if json_dict["uniprotSearch"] == "true":
        # If it comes straight from uniprot #
        json_dict["TF_name"] = json_dict["input"]["_id"].split("|")[-1]
    else:
      if json_dict["input"]["_id"].count("|") >= 1:
        # If it has an uniprot-like name #
         json_dict["input"]["_id"]= json_dict["input"]["_id"].replace("|","_")
    if "," in json_dict["input"]["_id"]:
        # If it is an heterodimer and has two protein names #
        json_dict["TF_name"] = json_dict["input"]["_id"].split(",")[0].lstrip(">") + "," + json_dict["input"]["_id"].split(",")[1].lstrip(">")
    else:
        json_dict["TF_name"] = json_dict["input"]["_id"].lstrip(">")

    # Clean upf DNA names without "_"
    for i in range(len(json_dict["dna"]["seq"])):
        json_dict["dna"]["seq"][i]["_id"] = json_dict["dna"]["seq"][i]["_id"].replace("_",":")

    # Create output directories #
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)
        os.system("chmod 777 %s"%output_dir)
    if os.path.exists(dummy_dir) == False:
        os.mkdir(dummy_dir)
        os.system("chmod 777 %s"%dummy_dir)
    if os.path.exists(os.path.join(output_dir, "TF_modeling")) == False:
        os.mkdir(os.path.join(output_dir, "TF_modeling"))  
        os.system("chmod 777 %s"%os.path.join(output_dir, "TF_modeling"))
    if os.path.exists(os.path.join(output_dir, "DNA")) == False:
        os.mkdir(os.path.join(output_dir, "DNA"))   
        os.system("chmod 777 %s"%os.path.join(output_dir, "DNA"))
    if os.path.exists(os.path.join(output_dir, "DNA", "binding_profiles")) == False:
        os.mkdir(os.path.join(output_dir, "DNA", "binding_profiles")) 
        os.system("chmod 777 %s"%os.path.join(output_dir, "DNA", "binding_profiles"))

    # Make accessible to evey user the folders
    os.system("chmod -R 777 %s"%output_dir)

    # Initialize the dictionary for motifs #
    if not "motifs" in json_dict:
        json_dict["motifs"] = {}

    # If input file is "input.fa"... #
    if input_file.endswith("input.fa") and json_dict['input']['type'] == 'FASTA':
        if verbose == True: sys.stdout.write("\nUser provided a FASTA file...\n")
        # Create PDB model #
        if verbose == True: sys.stdout.write("\t\t-- creating PDB model...\n")
        # Get the length of the input protein #
        protein_length = []
        for line in open(input_file, "r").readlines()[1:]:
            if line.startswith(">"):
                protein_length.append(0)
            else:
                if protein_length == []:
                    protein_length.append(0)
                protein_length[-1] += len(line.rstrip())
        json_dict['protein_length'] = protein_length
        # Check if the input fasta is empty #
        if len(open(input_file, "r").readlines()) <= 1:
            if not json_dict.has_key("error_msg"):
              json_dict["error_msg"] = "The input fasta sequence is empty."
            return json_dict

        # Make the models
        pdb_files = tf_modeling(json_dict, input_file, output_dir, remote,dummy_dir)
        # If no output was produced... #
        if pdb_files == []:
            if not json_dict.has_key("error_msg"):
              json_dict["error_msg"] = "No models could be created!"
            return json_dict

    # If input file is "input.pdb"... #
    if input_file.endswith("input.pdb") and json_dict['input']['type'] == 'PDB':
        if verbose == True: sys.stdout.write("\nUser provided a PDB file...\n")
        # Initialize #
        print("PDB in JOB %s : %s"%(json_dict["_id"],input_file))
        pdb_obj = PDB(input_file)
        # Unless it is a PDB... #
        if len(pdb_obj.chains) == 0:
            if not json_dict.has_key("error_msg"):
              json_dict["error_msg"] = "No chains could be read from input PDB file!"
            return json_dict
        # Unless PDB has protein... #
        if pdb_obj.has_protein == False:
            if not json_dict.has_key("error_msg"):
              json_dict["error_msg"] = "No protein chains could be read from input PDB file!"
            return json_dict
        # Unless PDB has DNA... #   
        if pdb_obj.has_nucleotide == False:
            if not json_dict.has_key("error_msg"):
             json_dict["error_msg"] = "No DNA chains could be read from input PDB file!"
            return json_dict
        # Add PDB file #
        seq=""
        for pdb_chain_obj in pdb_obj.chains:
              if pdb_chain_obj.chaintype =='P':
                 seq= seq + str(pdb_chain_obj.gapped_protein_sequence.replace("-",""))
        print("PDB in JOB %s is %s with length %d"%(json_dict["_id"],input_file,len(seq)))
        pdb_file = "model:1:"+str(len(seq))+"_input_1.pdb"
        print("JOB %s write PDB as %s"%(json_dict["_id"],os.path.join(output_dir,"TF_modeling",pdb_file)))
        if not os.path.exists(os.path.join(output_dir,"TF_modeling",pdb_file)):
           pdb_obj.write(os.path.join(output_dir,"TF_modeling",pdb_file))
        pdb_files = [os.path.join(output_dir,"TF_modeling",pdb_file)]

    if input_file.endswith(".txt") and json_dict['input']['type'] == 'threading':
        pdb_file = tf_threading(input_file, output_dir, pdb_dir, dummy_dir)
        pdb_obj.write(os.path.join(output_dir,"TF_modeling",input_file))
        pdb_files = [os.path.join(output_dir,"TF_modeling",input_file)]
        # If no output was produced... #
        if pdb_files == []:
            if not json_dict.has_key("error_msg"):
              json_dict["error_msg"] = "No models could be created!"
            return json_dict

    # Parse the families file #
    families = {}
    for line in functions.parse_file(os.path.join(pdb_dir, "families.txt")):
        if line.startswith("#"): continue
        pdb_chain, family = line.split(";")
        families[pdb_chain] = family

    #reference of DNA in FastA
    n_dna=0
    dna_fa={}
    for dna in json_dict["dna"]["seq"]:
          n_dna   = n_dna+1
          dna_id  = dna["_id"]
          dna_fa.setdefault(dna_id,dna_id.lstrip(">")+"_"+str(n_dna))

    # Generate profiles
    print("MP2D: GENERATE PROFILES")
    try:
      generate_profiles(json_dict,output_dir,remote)
    except Exception as e:
      json_dict["error_msg"] = "No profiles generated: "+str(e)
      return json_dict

    # Rearrange files hyerarchies for WEB
    print("MP2D: ARRANGE PROFILES")
    try:
      fragment_list=arrange_profiles_for_web(json_dict,output_dir,pdb_files)
      for frg in fragment_list:
        print("MP2D: Domain "+frg)
    except:
      json_dict["error_msg"] = "No reallocation of files "
      return json_dict

    # Sort PDB files #
    pdb_files = sorted(pdb_files)

    # Cluster pdb files
    print("MP2D:  CLUSTER PDBS")
    try:
      #clusters=PROFILER.classify_pdbs(set(pdb_files))
      if len(fragment_list)>0: clusters=update_classified_pdbs(pdb_files)
      else: clusters=PROFILER.classify_pdbs(set(pdb_files))
    except:
      json_dict["error_msg"] = "Wrong domains "
      return json_dict

    # Generate motif_list information 
    print("MP2D:   GENERATE MOTIF LIST")
    #try:
    motif_list = generate_motif_list(pdb_files,families,output_dir,json_dict,dna_fa)
    #except:
    #  json_dict["error_msg"] = "Missing motifs for the analysis"
    #  return json_dict

    #Group motif list  by cluster number         
    print("MP2D:   GROUP MOTIF LIST")
    groups    = group_motifs_by_clusters(clusters,motif_list)
    fragment  = get_cluster_fragments(clusters,fragment_list)
    for clu in fragment.keys():
        print("MP2D: Cluster "+str(clu)+" Domain "+fragment[clu])


    #Motif JSON
    #json_dict["motifs"]     = [groups.get(number_of_cluster)[:min(10,len(groups.get(number_of_cluster)))] for number_of_cluster in clusters.iterkeys()]
    print("MP2D:   MAKE JSON MOTIFS AND INFORMATION")
    json_dict["motifs"]=[]
    try:
      for number_of_cluster in clusters.iterkeys():
         try:
            if number_of_cluster is None: continue
            if not groups.has_key(number_of_cluster): continue
            if groups.get(number_of_cluster)[:] is None: continue
            json_dict["motifs"].append(groups.get(number_of_cluster)[:])
         except:
            continue
      json_dict["motif_info"] = get_motif_info(json_dict["motifs"])
    except:
      json_dict["error_msg"] = "Missing description of domains"
      return json_dict

    #Binding profiles JSON
    print("MP2D:  GENERATE BINDING PROFILES")
    try:
      json_dict["binding_profiles"] = generate_binding_profiles_json(clusters,fragment,output_dir,json_dict,dna_fa)
    except:
      json_dict["error_msg"] = "Missing global profiles"
      return json_dict


    return json_dict

def protein_2_dna_update(json_dict, output_dir, remote=False, verbose=True):

    # Initialize #
    input_file = os.path.join(output_dir, json_dict["input"]["file"])
    dummy_dir = os.path.join(output_dir, "dummy")
    pdb_dir = config.get("Paths", "pdb_dir")
    if not os.path.exists(os.path.join(output_dir, "DNA", "binding_profiles")): os.mkdir(os.path.join(output_dir, "DNA", "binding_profiles")) 
    # Make accessible to evey user the folders
    os.system("chmod -R 777 %s"%output_dir)

    # Clean upf DNA names without "_"
    for i in range(len(json_dict["dna"]["seq"])):
        json_dict["dna"]["seq"][i]["_id"] = json_dict["dna"]["seq"][i]["_id"].replace("_",":")

    #reference of DNA in FastA and new working DNA sequences
    n_dna=0
    dna_fa={}
    for dna in json_dict["dna"]["seq"]:
          n_dna   = n_dna+1
          dna_id  = dna["_id"]
          dna_fa.setdefault(dna_id,dna_id.lstrip(">")+"_"+str(n_dna))
    working_dna = []
    for dna in json_dict["dna"]["seq"]:
        working_dna.append(dna["_id"])

    # Parse the families file #
    families = {}
    for line in functions.parse_file(os.path.join(pdb_dir, "families.txt")):
        if line.startswith("#"): continue
        pdb_chain, family = line.split(";")
        families[pdb_chain] = family

    # Remove binding profiles of previous DNA sequences
    try:
      clean_up_profiles(json_dict,working_dna)
    except:
      json_dict["error_msg"] = "Previous data was unclean"
      return json_dict

    # Handle the input pwms if provided #
    input_pwms_data = []
    try:
      if "input_pwms" in json_dict.keys(): input_pwms_data = handle_input_pwms(json_dict["input_pwms"])
    except:
      json_dict["error_msg"] = "error handling PWMs"
      return json_dict

    # Load PWM and remove previous PWMs
    try:
      print("MP2D: LOAD PWM")
      load_input_pwm(input_pwms_data,json_dict,output_dir)
    except:
      json_dict["error_msg"] = "Wrong format of external PWMs"
      return json_dict

    # Update new profiles
    try:
      print("MP2D: UPDATE PROFILES")
      update_profiles(json_dict,output_dir,remote)
    except Exception as e:
      json_dict["error_msg"] = "No updated profiles: "+str(e)
      return json_dict
      

    # Get the output files # 
    try:
      print("MP2D: READ PDBS")
      pdb_files = get_pdb_models(output_dir)
    except:
      json_dict["error_msg"] = "No models"
      return json_dict

    # Rearrange files hyerarchies for WEB
    try:
      print("MP2D: ARRANGE PROFILES")
      fragment_list=arrange_profiles_for_web(json_dict,output_dir,pdb_files)
      for frg in fragment_list:
        print("MP2D: Domain "+frg)
    except:
      json_dict["error_msg"] = "No reallocation of updated files "
      return json_dict

    # Cluster pdb files
    try:
      print("MP2D: CLASSIFY PDB")
      if len(fragment_list)>0: clusters=update_classified_pdbs(pdb_files)
      else: clusters=PROFILER.classify_pdbs(set(pdb_files))
    except:
      print "ERROR"
      json_dict["error_msg"] = "Wrong updated domains "
      return json_dict

    # Generate motif_list information 
    try:
      print("MP2D: GENERATE MOTIF LIST")
      motif_list = generate_motif_list(pdb_files,families,output_dir,json_dict,dna_fa)
    except:
      json_dict["error_msg"] = "Missing updated motifs for the analysis"
      return json_dict

    #Group motif list  by cluster number   
    print "MOTIF LIST"
    print motif_list
    print("MP2D: LIST AND REORDER MOTIFS")
    groups    = group_motifs_by_clusters(clusters,motif_list)
    fragment  = get_cluster_fragments(clusters,fragment_list)
    for clu in fragment.keys():
        print("MP2D: Cluster "+str(clu)+" Domain "+fragment[clu])

    #Motif JSON
    print("MP2D:  COPY MOTIFS")
    json_dict["motifs"] = [groups.get(number_of_cluster)[:] for number_of_cluster in clusters.iterkeys()]

    #Binding profiles JSON
    try:
      print("MP2D: GENERATE BINDING PROFILES")
      json_dict["binding_profiles"] = generate_binding_profiles_json(clusters,fragment,output_dir,json_dict,dna_fa)
    except:
      json_dict["error_msg"] = "Missing global updated profiles"
      return json_dict


    return json_dict


#---------------#
# Main          #
#---------------#

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Define which python to be used #
python = os.path.join(config.get("Paths", "python_path"), "python")

# Directories #
pbm_dir = config.get("Paths", "pbm_dir")
pdb_dir = config.get("Paths", "pdb_dir")
root_path = config.get("Paths", "root_path") 

# Load the bashrc file from apache #

if __name__ == "__main__":


    # Arguments & Options #
    options = parse_options()
    print("Root path is")
    print(root_path)

    # Execute protein 2 DNA #
    if options.update:
        json_dict = protein_2_dna_update(json_dict=json.load(open(os.path.abspath(options.json_file))), output_dir=os.path.abspath(options.output_dir), remote=options.remote, verbose=options.verbose)
    else:
        json_dict = protein_2_dna(json_dict=json.load(open(os.path.abspath(options.json_file))), output_dir=os.path.abspath(options.output_dir), remote=options.remote, verbose=options.verbose)
    # Save json_dict #
    if json_dict != None:
        print "DONE %s"%str(json_dict)
        # Create json file #
        out = open(options.json_file, "wt")
        out.write(json.dumps(json_dict, separators=(',', ':'), indent=2))
        out.close()



 
