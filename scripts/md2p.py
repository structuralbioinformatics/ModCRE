import os, sys, re
import ConfigParser
import json
import optparse
import shutil
import subprocess
import imp
import math
import pickle
import numpy as np
from HTMLParser import HTMLParser
from bs4 import BeautifulSoup
import urllib
import operator

# Import sklearn modules #
from sklearn.cluster import SpectralClustering
from sklearn import metrics

# Add "scripts" to sys.path #
scripts_path = os.path.join(os.path.dirname(os.path.abspath(__file__)))

# Imports jbonet's module #
from SBI.structure          import PDB
from SBI.structure.chain    import Chain
from SBI.structure.chain    import ChainOfProtein
from SBI.structure.chain    import ChainOfNucleotide
from SBI.data               import aminoacids1to3
from SBI.external.blast     import blast_parser

sys.path.append(scripts_path)

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Define which python to be used #
python = os.path.join(config.get("Paths", "python_path"), "python")

# Imports my functions #
import functions

# Import modules #
import blast, contacts, dssp, interface, model_protein, model_dna, x3dna, triads, pwm_pbm, scorer, functions, fimo, tomtom 
import pwm_pbm as PWM

from JSONfile import dumpJSON, readJSON
#import tf_complex


# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Define which python to be used #
python = os.path.join(config.get("Paths", "python_path"), "python")

# Directories #
pbm_dir = config.get("Paths", "pbm_dir")
pdb_dir = config.get("Paths", "pdb_dir")
jaspar_dir = config.get("Paths", "jaspar_dir")
cisbp_dir  = config.get("Paths", "cisbp_dir")
root_path = config.get("Paths", "root_path") 


#-------------#
#  Functions  #
#-------------#

def get_specie_code(specie):

    if (specie == "null") or (specie == None) :
        return None 
    else:
        # Search for the specie code in a file #
        files_dir = config.get("Paths", "files_path")
        specielist = config.get("Paths", "species")
        for line in open(os.path.join(files_dir, specielist), "r").readlines():
            if specie in line:
                main_fields = line.strip().split("=")
                if specie != main_fields[1]:continue
                fields = line.split()
                if len(fields)<3: continue
                specie_code = fields[2].replace(":", "")
                return specie_code


def get_userdb(folder,fasta_file,mdl2fasta,uid_set,uid_db,label):

    #Initialize

    src_path = config.get("Paths", "src_path")
    blast_path = os.path.join(src_path, config.get("Paths", "blast_path"))
    tz_parameter = int(config.get("Parameters", "twilight_zone_parameter"))
    tz_type = config.get("Parameters", "twilight_zone_type")
    mkblastdb = os.path.join(blast_path,"makeblastdb")
    if not os.path.exists(os.path.join(folder,"pwms_DB")):
        os.mkdir(os.path.join(folder,"pwms_DB"))
        os.system("chmod -R 777 %s"%os.path.join(folder,"pwms_DB"))
    if not os.path.exists(os.path.join(folder,"results")):
        os.mkdir(os.path.join(folder,"results"))
        os.system("chmod -R 777 %s"%os.path.join(folder,"results"))
    if not os.path.exists(os.path.join(folder,"homologs")):
        os.mkdir(os.path.join(folder,"homologs"))
        os.system("chmod -R 777 %s"%os.path.join(folder,"homologs"))
    homologs_dir = os.path.join(folder,"homologs")
    pwm_dir      = os.path.join(folder,"pwms_DB")
    database_pwm = os.path.join(pwm_dir,"database.txt")
    if os.path.exists(database_pwm): os.rename(database_pwm,database_pwm+".previous")
    # Built datasets
    # PWMs
    association={}
    pdb_db=set()
    if os.path.exists(database_pwm): shutil.rmtree(database_pwm)
    if os.path.exists(uid_db):
        os.system("cat %s >> %s"%(uid_db,database_pwm))
    for pwm_file in os.listdir(folder):
        if not pwm_file.endswith(".meme"):continue
        done=False
        for uid in uid_set:
            if uid in pwm_file: done=True
        pdb_chain = pwm_file.rstrip(".meme").split("_")[-3]+"_"+pwm_file.rstrip(".meme").split("_")[-2]
        if done:
           association.setdefault(uid,pdb_chain)
           pdb_db.add(pdb_chain)
        else:
           msa_obj = PWM.nMSA(os.path.join(folder,pwm_file),pwm_file.rstrip(".meme"),"meme")
           if msa_obj.get_binding_site_length() >0 :
              os.system("cat %s >> %s"%(os.path.join(folder,pwm_file),database_pwm))
              association.setdefault(pwm_file.rstrip(".meme"),pdb_chain)
              pdb_db.add(pdb_chain)
           else:
              print("-- Skip PWM %s"%pwm_file)
    association_file = open(os.path.join(folder,"results","userdb2pdb.txt"),"w")
    for pwm in association.keys():
        association_file.write("%s\t%s\n"%(pwm,association[pwm]))
    association_file.close()
    # Sequences
    sequence_db = os.path.join(homologs_dir,"sequences_db.fa")
    shutil.copy(fasta_file,sequence_db)
    #dbseq = open(sequence_db,"w")
    #for name,sequence in functions.parse_fasta_file(fasta_file):
    #    done=False
    #    for uid in uid_set:
    #        if uid in name: done=True
    #    if done: protein = name.replace("|","_").split()[0]
    #    else: protein = label+"_"+name.replace("|","_").split()[0]
    #    dbseq.write(">%s\n%s\n"%(protein.lstrip(">"),sequence))
    #dbseq.close()
    # homologs#
    print("MAKE BLASTDB: %s -in %s -dbtype prot -out %s"%(mkblastdb,sequence_db,sequence_db))
    os.system("%s -in %s -dbtype prot -out %s"%(mkblastdb,sequence_db,sequence_db))
    for pdb_chain in pdb_db:
        print("\t-- SEARCH HOMOLOGS OF: %s"%pdb_chain)
        homolog_file = os.path.join(homologs_dir,pdb_chain+".txt")
        blast_obj = blast.get_blast_obj( sequence_db ,os.path.join(pdb_dir,"split",pdb_chain+".fasta"))
        functions.write(homolog_file, blast_obj.str_compacted_blast(tz_parameter=tz_parameter, tz_type=tz_type))


    # return external DB data
    os.system("chmod -R 777 %s"%folder)
    return database_pwm, os.path.join(folder,"results","userdb2pdb.txt"), homologs_dir


def get_pwm_database(json_dict, output_dir, dummy_dir, remote=True):

    print("Construction of the PWM database")
    if not os.path.exists(os.path.join(output_dir, "pwm_database")): 
        os.mkdir(os.path.join(output_dir, "pwm_database"))
        os.system("chmod -R 777 %s"%os.path.join(output_dir, "pwm_database"))
    if not os.path.exists(os.path.join(output_dir, "pwm_database","pwms")): 
        os.mkdir(os.path.join(output_dir, "pwm_database","pwms"))
        os.system("chmod -R 777 %s"%os.path.join(output_dir, "pwm_database","pwms"))


    database_file = os.path.join(output_dir, "pwm_database", "pwm_database.txt")
    database_dir  = os.path.join(output_dir, "pwm_database","pwms")
    uid_set=set()
    database = False

    # jaspar DB #
    print("MD2P: jaspar_pwm_db = " + str(json_dict["jaspar_pwm_db"]))
    print("MD2P: Jaspar DB data = "+os.path.join(jaspar_dir, "results", "species_families.out"))
    if not database and str(json_dict["jaspar_pwm_db"]) == "true":
        print("MD2P: CHECK SPECIE ... ")
        print(json_dict["specie"])
        if str(json_dict["specie"]) == "null" or str(json_dict["specie"])=="None":
            json_dict["specie"] = "Homo sapiens"
        print("MD2P: JASPAR DB SPECIE IS .... "+ json_dict["specie"])
        # Create your own database selecting only those pwms from the specie you are interested in #
        pwm_files = []
        for line in open(os.path.join(jaspar_dir, "results", "species_families.out")).readlines():
            fields = line.split("\t")
            #print("MD2P: Check family "+str(fields[2].rstrip().upper())+" to select family "+str(json_dict["specie"]).upper())
            if len(json_dict["specie"].split())>1:
              if str(fields[2].rstrip().upper()) == str(json_dict["specie"]).upper():
                shutil.copy(os.path.join(jaspar_dir, "pwms", fields[0] + ".meme.s"),os.path.join(database_dir, fields[0] + ".meme"))
                #print("MD2P: COPY "+fields[0] + ".meme.s"+" TO "+os.path.join(database_dir, fields[0] + ".meme"))
                pwm_files.append(os.path.join(database_dir, fields[0] + ".meme"))
            else:
              genus = fields[2].rstrip().split()[0]
              epithet = " ".join(fields[2].rstrip().split()[1:])
              name_specie = genus.title()+" "+epithet.lower()
              code = get_specie_code(name_specie)
              if code == str(json_dict["specie"]).upper():
                shutil.copy(os.path.join(jaspar_dir, "pwms", fields[0] + ".meme.s"),os.path.join(database_dir, fields[0] + ".meme"))
                #print("MD2P: COPY "+fields[0] + ".meme.s"+" TO "+os.path.join(database_dir, fields[0] + ".meme"))
                pwm_files.append(os.path.join(database_dir, fields[0] + ".meme"))


        if pwm_files == []:
            print("no pwms files in JASPAR database ")
            database = False
            database_file = None
            database_dir  = None
            association_file = None
            homologs_dir = None
            error_msg = "no pwms files in JASPAR database"
            return json_dict, database, database_dir, database_file,association_file,homologs_dir,error_msg
        else:
            database = True
            for file in pwm_files:
                if os.path.isfile(database_file):
                    os.system("cat " + file + " >> " + database_file)
                else:
                    os.system("cat " + file + " > " + database_file)
            association_file = os.path.join(jaspar_dir, "results","jaspar2pdb.out")
            homologs_dir = None
            error_msg = None

    # cisbp DB #
    if not database and str(json_dict["cisbp_pwm_db"]) == "true":
        print("MD2P: CHECK SPECIE ... ")
        print(json_dict["specie"])
        if str(json_dict["specie"]) == "null" or str(json_dict["specie"]) == "None":
            json_dict["specie"] = "Homo sapiens"
        print("MD2P: CisBP DB SPECIE IS .... "+ json_dict["specie"])

        # Create your own database selecting only those pwms from the specie you are interested in #
        pwm_files = []
        for line in open(os.path.join(cisbp_dir, "results", "species_families.out")).readlines():
            fields = line.split("\t")
            #print("MD2P: cisBP Fields species "+str(fields))
            if len(json_dict["specie"].split())>1:
              if str(fields[2].rstrip().upper()) == str(json_dict["specie"]).upper():
                shutil.copy(os.path.join(cisbp_dir, "pwms", fields[0] + ".meme"),os.path.join(database_dir, fields[0] + ".meme"))
                print("MD2P: COPY "+fields[0] + ".meme"+" TO "+os.path.join(database_dir, fields[0] + ".meme"))
                pwm_files.append(os.path.join(database_dir, fields[0] + ".meme"))
            else:
              genus = fields[2].rstrip().split()[0]
              epithet = " ".join(fields[2].rstrip().split()[1:])
              name_specie = genus.title()+" "+epithet.lower()
              code = get_specie_code(name_specie)
              if code == str(json_dict["specie"]).upper():
                shutil.copy(os.path.join(cisbp_dir, "pwms", fields[0] + ".meme"),os.path.join(database_dir, fields[0] + ".meme"))
                print("MD2P: COPY "+fields[0] + ".meme"+" TO "+os.path.join(database_dir, fields[0] + ".meme"))
                pwm_files.append(os.path.join(database_dir, fields[0] + ".meme"))

        if pwm_files == []:
            print("no pwms files in CISBP database")
            database = False
            database_file = None
            database_dir  = None
            association_file = None
            homologs_dir = None
            error_msg = "no pwms files in CISBP database"
            return json_dict, database, database_dir, database_file,association_file,homologs_dir,error_msg

        else:
            database = True
            for file in pwm_files:
                if os.path.isfile(database_file):
                    os.system("cat " + file + " >> " + database_file)
                else:
                    os.system("cat " + file + " > " + database_file)
            association_file = os.path.join(cisbp_dir, "results","cisbp2pdb.out")
            homologs_dir = None
            error_msg = None

    #If user selects PWMS and proteins
    if json_dict.has_key("custom_pwms") and json_dict.has_key("custom_proteins"):
       if  str(json_dict["custom_pwms"]) != "null" and str(json_dict["custom_proteins"]) != "null":
            if len(str(json_dict["custom_pwms"]).split(",")) == len(str(json_dict["custom_proteins"]).split(",")):
                database = True
                database_file_up = database_file + ".custom"
                pwms = str(json_dict["custom_pwms"]).split(",")
                uni_ids = str(json_dict["custom_proteins"]).split(",")
                uid_set = uid_set.update(set(uni_ids))
                for i in range(0, len(str(json_dict["custom_pwms"]).split(","))):
                    pwm_content = pwms[i].splitlines()
                    uni_motif = uni_ids[i].rstrip().replace("\r", "").replace("\n", "").replace(" ", "")
                    uni_pwm = open( os.path.join(database_dir,uni_motif + ".meme"),"w")
                    with open(database_file_up, 'ab') as db_file:
                        for line in pwm_content:
                            if line.startswith(","):
                                continue
                            if line.startswith("MOTIF"):
                                line = "MOTIF " + uni_ids[i].rstrip().replace("\r", "").replace("\n", "").replace(" ", "") + " " + uni_ids[i].rstrip().replace("\r", "").replace("\n", "").replace(" ", "")
                            db_file.write(line.rstrip() + "\n")
                            uni_pwm.write(line.rstrip() + "\n")
                        db_file.close()
                        uni_pwm.close()
            else:
                database = True
                database_file = None
                database_dir  = None
                association_file = None
                homologs_dir = None
                error_msg = "The number of input PWMs and input Uniprot IDs is different"
                return json_dict, database, database_dir, database_file,association_file,homologs_dir,error_msg

    if not json_dict.has_key("custom_proteins"): json_dict["custom_proteins"] = "null"


    selected_proteins = False
    if json_dict.has_key("prot"):
        if len(json_dict["prot"]["seq"])>0: selected_proteins=True

    if str(json_dict["custom_proteins"]) != "null" or selected_proteins:  
            # Here we have to execute a whole protein2DNA workflow to get pwms from uniprot IDs #
            # Create subdirectory structure #
            if not os.path.exists(os.path.join(output_dir, "pwm_database", "sequences")):
                    os.mkdir(os.path.join(output_dir, "pwm_database", "sequences"))
                    os.system("chmod -R 777 %s"%os.path.join(output_dir, "pwm_database", "sequences"))
            if not os.path.exists(os.path.join(output_dir, "pwm_database", "models")):
                    os.mkdir(os.path.join(output_dir, "pwm_database", "models"))
                    os.system("chmod -R 777 %s"%os.path.join(output_dir, "pwm_database", "models"))
            # Parse the input file #
            sequence_db = open(os.path.join(output_dir, "pwm_database", "sequences", "sequences.fasta"), "w")
            model_db    = open(os.path.join(output_dir, "pwm_database", "sequences", "models.fasta"), "w")
            fasta_file  = os.path.join(output_dir, "pwm_database", "sequences", "sequences.fasta")
            mdl2fasta   = {}
            uniprot_ids = []
            if str(json_dict["custom_proteins"]) != "null":
                uniprot_ids = str(json_dict["custom_proteins"]).split(",")
                print("uniprot_ids: " + str(uniprot_ids))
                # Iterate over uniprot IDs #
                for uid in uniprot_ids:
                    # Get the protein sequence #
                    uid = uid.rstrip().replace("\r", "").replace("\n", "").replace(" ", "")
                    url = 'http://www.uniprot.org/uniprot/' + uid + '.fasta'
                    f = urllib.urlopen(url)
                    seq = f.read()
                    m = re.search("^(tr|sp)\|(\S+)\|\S+\_(\S+)(.+)",seq)
                    if m:
                       mdl2fasta.setdefault(m.group(2),uid)
                       mdl2fasta.setdefault(uid,m.group(2))
                    m = re.search("^(tr|sp|TR|SP)\|(\S+)\|\S+\_(\S+)",seq)
                    if m:
                       mdl2fasta.setdefault(m.group(2),uid)
                       mdl2fasta.setdefault(uid,m.group(2))
                    # Save the protein in a dummy fasta file #
                    fasta = open(os.path.join(output_dir, "pwm_database", "sequences", uid + ".fasta"), "w")
                    fasta.write(seq)
                    sequence_db.write(seq)
                    fasta.close()
            if len(json_dict["prot"]["seq"])>0:
                for data in json_dict["prot"]["seq"]:
                    name=data["_id"].lstrip(">")
                    seq =data["full_sequence"]
                    m = re.search("^(tr|sp|TR|SP)\|(\S+)\|\S+\_(\S+)",name)
                    if m:
                        sequence_db.write(">%s\n%s\n"%(name,seq))
                        model_db.write(">%s\n%s\n"%(m.group(2),seq))
                        mdl2fasta.setdefault(m.group(2),name)
                        mdl2fasta.setdefault(name,m.group(2))
                    else:
                        print("SKIP SEQUENCE %s :  Wrong format sp|ACCESSION|GENE_SPECIE"%name)
            for uid in uniprot_ids:
                  uid = uid.rstrip().replace("\r", "").replace("\n", "").replace(" ", "")
                  for name,seq in functions.parse_fasta_file(os.path.join(output_dir, "pwm_database", "sequences", uid + ".fasta")):
                    model_db.write(">%s\n%s\n"%(uid.lstrip(">"),seq))
            sequence_db.close()
            model_db.close()

            # Model the proteins #
            dummy_dir  = os.path.join(output_dir,"dummy")
            label      = "model"
            parameters = " -i " +  os.path.join(output_dir, "pwm_database", "sequences", "models.fasta")
            parameters = parameters + " -o " + os.path.join(output_dir,"pwm_database", "models" ) 
            parameters = parameters + " -l " + label
            parameters = parameters + " -p " + pdb_dir 
            parameters = parameters + " --dummy=" + dummy_dir 
            parameters = parameters + " -f -a -v --renumerate --chains_fixed --parallel "
            logfile = os.path.join(output_dir,"modelling_db.log")
            command="model_multiple_proteins.py"
            # Execute modeling in remote#
            if remote:
               print("RUN %s ON REMOTE"%command)
               job = functions.execute_in_remote(command,parameters,config,output_dir,logfile,waiting=True)
               print("FINISHED REMOTE JOB: "+job)
               if str(job).startswith("Error"):
                  json_dict["error_msg"] =job
                  return json_dict
            else:
            # Execute modeling in local#
               print("RUN %s ON LOCAL \n"%(command))
               print("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))
               os.system("%s %s/%s %s \n"%(python,scripts_path,command,parameters))

            if not database or len(json_dict["prot"]["seq"])>0:
              # Model PWMs
              dummy_dir  = os.path.join(output_dir,"dummy")
              parameters = " -i " + os.path.join(output_dir,"pwm_database", "models" )
              parameters = parameters + " -o " + os.path.join(output_dir,"pwm_database", "pwms" )
              parameters = parameters + " --info " + os.path.join(output_dir,"pwm_database", "pwms","pwms_list.log")
              parameters = parameters + " --pdb " + pdb_dir 
              parameters = parameters + " --pbm " + pbm_dir 
              parameters = parameters + " --dummy=" + dummy_dir 
              parameters = parameters + " -a -v --reuse --parallel "
              logfile    = os.path.join(output_dir,"pwm_database", "pwms","pwm_by_models.log")
              command="pwm_pbm.py"
              # Execute PWM in remote#
              if remote:
                  print("RUN %s ON REMOTE"%command)
                  job = functions.execute_in_remote(command,parameters,config,output_dir,logfile,waiting=True)
                  print("FINISHED REMOTE JOB: "+job)
                  if str(job).startswith("Error"):
                    json_dict["error_msg"] =job
                    return json_dict
              else:
              # Execute PWM in local#
                  print("RUN %s ON LOCAL \n"%(command))
                  print("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))
                  os.system("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))
            uniprot_db = database_file + ".custom"
            folder=os.path.join(output_dir,"pwm_database" ,"pwms")
            if not os.path.exists(folder):
                os.mkdir(folder)
                os.system("chmod -R 777 %s"%(folder))
            try:
              database_pwm, association_pwm, homologs_dir = get_userdb(folder,fasta_file,mdl2fasta,uid_set,uniprot_db,label) 
              shutil.copy(database_pwm,database_file)
              database_dir = folder
              association_file = association_pwm
              error_msg=None
              database = True
            except:
              error_msg="Fail when constructing database of PWMs: check the format of the protein sequences"
              database_dir = folder
              association_file = association_pwm
              database = True 

    if not database:
        if str(json_dict["specie"]) == "null" or str(json_dict["specie"]) == "None":
            json_dict["specie"] = "9606"
        shutil.copy(os.path.join(pbm_dir, "pwms", "database.txt"), database_file)
        database_dir = os.path.join(pbm_dir, "pwms")
        association_file = None
        homologs_dir = None
        error_msg =  None

    return json_dict, database, database_dir, database_file,association_file,homologs_dir,error_msg


def run_scan(json_dict, output_dir, database_external, database_dir, database_file, association_file, homologs_dir , scan_dir,  dummy_dir, remote=True):

    input_file = os.path.join(output_dir, json_dict["input"]["file"])
    scan_dir   = os.path.join(output_dir,"scan")
    dummy_dir  = os.path.join(output_dir,"dummy_scan")
    parameters  = " -i " + input_file
    parameters  = parameters + " -o " + scan_dir
    parameters = parameters + " --pdb " + pdb_dir 
    parameters = parameters + " --pbm " + pbm_dir 
    parameters = parameters + " --dummy " + dummy_dir 
    parameters = parameters + " --parallel --reuse "
    parameters = parameters + " --verbose "
    # Check if all orthologs are requested, then ranking is unnecessary
    rank=True
    if json_dict.has_key("all_orthologs"):
        if json_dict["all_orthologs"]=="true": rank=False
    if rank:  parameters = parameters + " --rank "
 
    # Check orthologs of a specie
    if json_dict.has_key("specie"):
        if str(json_dict["specie"]) not in ["null", "None"]:
           if len(json_dict["specie"].split())>1:
              genus = json_dict["specie"].split()[0]
              epithet = " ".join(json_dict["specie"].split()[1:])
              name_specie = genus.title()+" "+epithet.lower()
              #json_dict["specie"] = name_specie
              json_dict["specie"] = get_specie_code(name_specie)
           #print("MD2P: RUN SCAN " + " --specie " + "'" + str(json_dict["specie"])) + "'"
           #parameters  = parameters  + " --specie " + "'" + str(json_dict["specie"]) + "'"
           print("MD2P: RUN SCAN " + " --specie " +  str(json_dict["specie"])) 
           parameters  = parameters  + " --specie " + str(json_dict["specie"]) 

    if json_dict.has_key("fimo_threshold"):
        if str(json_dict["fimo_threshold"]) not in ["null", "None"]:
           print("fimo_threshold: " + str(json_dict["fimo_threshold"]))
           parameters  = parameters  + " --ft=" + str(json_dict["fimo_threshold"])

    if database_external:
        if homologs_dir is None:
           parameters  = parameters  + " --external " + ",".join([database_dir, database_file, association_file])
        else:
           parameters  = parameters  + " --external " + ",".join([database_dir, database_file, association_file,homologs_dir])

    logfile = os.path.join(output_dir,"run_scanner.log")
    command = "scanner.py"

    # Execute PWM in remote#
    if remote:
          print("RUN %s ON REMOTE for %s PARAMETERS %s "%(command,os.path.basename(output_dir),parameters))
          job = functions.execute_in_remote(command,parameters,config,output_dir,logfile,waiting=True)
          print("FINISHED JOB %s for %s"%(str(job),os.path.basename(output_dir)))
          if str(job).startswith("Error"):
              json_dict["error_msg"] =job
              return json_dict

    else:
    # Execute PWM in local#
          print("RUN %s ON LOCAL \n"%(command))
          print("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))
          os.system("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))

    # Get the output files #
    if os.path.isfile(os.path.join(scan_dir, "orthologs_with_best_templates.json")):
        json_dict["orthologs_with_best_templates"] = os.path.join(scan_dir, "orthologs_with_best_templates.json")
    else:
        json_dict["error_msg"] = "Scan without orthologs"

    # Get the output files #
    if os.path.isfile(os.path.join(scan_dir, "orthologs.json")):
        json_dict["orthologs"] = os.path.join(scan_dir, "orthologs.json")
    else:
        json_dict["error_msg"] = "Scan witout orthologs"

    return json_dict


def get_protein_name(uniprot_id):

    print(str(uniprot_id))
    prot_name = "Not found"
    try:
        html = urllib.urlopen("http://www.uniprot.org/uniprot/" + uniprot_id + ".txt")
        soup = str(BeautifulSoup(html.read())).split("\n")
        for line in soup:
            if "Name: Full=" in line:
                prot_name = line.split("Name: Full=")[1].replace(";", "")
                if "{" in prot_name:
                    prot_name = prot_name.split("{")[0]
                break
    except:
        prot_name = "Not found"

    if prot_name == "Not found":
        print(uniprot_id + " not found\n")
    return prot_name


def process_scan_outputs(json_dict,  output_dir, dummy_dir):

    # Check if all orthologs are requested, then ranking is unnecessary
    rank=True
    if json_dict.has_key("all_orthologs"):
        if json_dict["all_orthologs"]=="true": rank=False

    # Parse the output of scan#
    scan_dir       = os.path.join(output_dir,"scan")
    if os.path.exists(os.path.join(scan_dir, "orthologs_with_best_templates.json")):
       orthologs_list = json.loads(''.join([line for line in functions.parse_file(os.path.join(scan_dir, "orthologs_with_best_templates.json"))]))
    else:
       orthologs_list = functions.parse_best_orthologs(os.path.join(scan_dir,"orthologs_with_best_templates.txt"),pdb_dir,config,rank)
       output_json    = os.path.join(scan_dir, "orthologs_with_best_templates.json")
       out = open(output_json, "wt")
       out.write(json.dumps(orthologs_list, separators=(',', ':'), indent=2))
       out.close()


    cluster_number = 0
    start_list=[]
    cluster_id=[]
    # the ortholog_dict is redefined as "cluster"
    new_list =[]
    for cluster in orthologs_list:
        cluster_number = cluster_number + 1
        start_list.append(cluster["start"])
        cluster["ID"] = str(cluster["start"])
        cluster_id.append(cluster["ID"])
        cluster["logpval"] = round(-np.log10(float(cluster["pval"])), 3)
        cluster["protein"] = ";".join([uid for uid,gene in cluster["proteins"]])
        cluster["protein_name"] = ";".join([gene for uid,gene in cluster["proteins"]])
        print("MD2P: PROTEIN = "+cluster["protein"])
        print("MD2P: PROTEIN_NAME = "+cluster["protein_name"])
        if len(cluster["families"])>0: cluster["family"] = ";".join([fam for fam in cluster["families"]])
        else:                          cluster["family"] ="unknown"
        cluster_number = cluster_number + 1
        if int(cluster["start"]) == 1:
           cluster["x_values"] = [x for x in range(cluster["start"], cluster["end"]+1)]
           cluster["y_values"] = [cluster["logpval"]]*(cluster["end"]-cluster["start"]) + [0]
        else:
           cluster["x_values"] = [x for x in range(cluster["start"]-1, cluster["end"]+1)]
           cluster["y_values"] =[0] + [cluster["logpval"]]*(cluster["end"]-cluster["start"]) + [0]
        cluster["zero_values"] = [0] * len(cluster["y_values"])
        cluster["orthologs"]=[[],[]]
        ortholog_threads={}
        for thread in cluster["monomer"]:
            thread_file, score, d_score = thread
            #Vector ortholog
            # thread_file, score, d_score, uid, gene/protein name, start, end, index cluster, logpval, accept to model (boolean),monomer(0) or dimer(1),dummy
            protein= ".".join(thread_file.lstrip("aux_files/").split(".")[:-3])
            #protein=thread_file.lstrip("aux_files/").split(".")[0]
            ortholog_threads.setdefault(protein,set()).add(thread_file)
            for i in range(len(cluster["protein"].split(";"))):
                if protein == cluster["protein"].split(";")[i]:
                    protein_name = cluster["protein_name"].split(";")[i]
            vector=[thread_file,round(float(score),2),round(float(d_score),2),protein,protein_name,cluster["start"],cluster["end"],cluster_number,round(float(cluster["logpval"]), 2),"",0,"dummy_cluster"]
            cluster["orthologs"][0].append(vector)
        for thread_dimer in cluster["dimer"]:
          for m in range(2):
            thread = thread_dimer[m]
            thread_file, score, d_score = thread
            #Vector ortholog
            # thread_file, score, d_score, uid, gene/protein name, start, end, index cluster, logpval, accept to model (boolean),monomer(0) or dimer(1),dummy
            #protein=thread_file.lstrip("aux_files/").split(".")[0]
            protein= ".".join(thread_file.lstrip("aux_files/").split(".")[:-3])
            ortholog_threads.setdefault(protein,set()).add(thread_file)
            for i in range(len(cluster["protein"].split(";"))):
                if protein == cluster["protein"].split(";")[i]:
                    protein_name = cluster["protein_name"].split(";")[i]
            vector=[thread_file,round(float(score),2),round(float(d_score),2),protein,protein_name,cluster["start"],cluster["end"],cluster_number,round(float(cluster["logpval"]), 2),"",m,"dummy_cluster"]
            cluster["orthologs"][m].append(vector)
        new_list.append(cluster)


    json_dict["output_dir"] = output_dir.replace(root_path,"")

    if new_list == []:
        json_dict["error_msg"] = "No transcription factor binding sites were found!"
    else:
        new_list = sorted(new_list, key=lambda k: k['start'])

    full_json_dict = json_dict
    full_json_dict["clusters"] = new_list

    return full_json_dict
    
def build_dna_structure(dna_seq, dna_str, output_dir,  begining=None, ending=None, verbose=False):

    x3dna_path = config.get("Paths", "x3dna_path")
    os.environ['X3DNA'] = x3dna_path[:-4]
    print("MD2P: build_dna_structure ")
    print("MD2P: X3DNA =  " + x3dna_path[:-4])

    if (begining != None) and (ending != None):
        dna_file = os.path.join(output_dir, "dna_" + str(begining) + "-" + str(ending) + ".pdb")
        modeling_seq = dna_seq[int(begining)-1:int(ending)]
    else:
        dna_file = os.path.join(output_dir, "dna_" + str(os.getpid()) + ".pdb")
        modeling_seq = dna_seq

    print("MD2P:\tDNA =>   begining: " + str(begining) + "  ending: " + str(ending) + "\n\tSEQUENCE: " + modeling_seq )


    if dna_str == "B_bent":
        histone_alphabet="CDEFGHIJKLMNOPQRSTUVWXYZ"
        # Check the length of the sequence #
        if len(modeling_seq) > 146:
            modeling_seq = modeling_seq[:146]
        # Open the reference structure of nucleosomal DNA #
        files_dir = config.get("Paths", "files_path")
        nuc_dna = PDB(os.path.join(files_dir, "nucleosomal_dna.pdb"))
        # Generate a new pdb with the nucleotide length that we are interested in #
        # This structure has two chains: I and J #
        new_nuc_dna = PDB()
        chain_1 = ChainOfNucleotide("dna", "A")
        chain_2 = ChainOfNucleotide("dna", "B")
        for i in range(0, len(modeling_seq)):
            res_for = nuc_dna.get_chain_by_id("I").nucleotides[i]
            #res_for.number = i + 1
            res_rev = nuc_dna.get_chain_by_id("J").nucleotides[-(i+1)]
            #res_rev.number = len(modeling_seq) - i
            chain_1.add_residue(res_for)
            chain_2.add_residue(res_rev)
        new_nuc_dna.add_chain(chain_1)
        new_nuc_dna.add_chain(chain_2)
        print("MD2P: new dna structure created at: " + os.path.join(output_dir, "nucleosomal_dna.pdb") + "\n\n")
        new_nuc_dna.write(output_file=os.path.join(output_dir, "nucleosomal_dna.pdb"), force=True)
        # Now modify the sequence of this pdb structure with mutate_bases #
        start_1 = chain_1.first_nucleotide.number
        end_1 = chain_1.last_nucleotide.number
        start_2 = chain_2.first_nucleotide.number
        end_2 = chain_2.last_nucleotide.number
        # Create the mutations file #
        complementary = {"A": "T", "C": "G", "G": "C", "T": "A"}
        mut_file = open(os.path.join(output_dir, str(os.getpid()) + ".mutations.dat"), "w")
        for i in range(0, len(modeling_seq)):
            mut_file.write("chain=%s seqnum=%s mutation=%s\n" % ("A", str(start_1+i), modeling_seq[i]))
            mut_file.write("chain=%s seqnum=%s mutation=%s\n" % ("B", str(start_2-i), complementary[modeling_seq[i]]))
        mut_file.close()
        print("MD2P: Mutation file created in: " + os.path.join(output_dir, str(os.getpid()) + ".mutations.dat") + "\n\n")
        # Execute mutate_bases #    
        print("MD2P: run X3DNA command "+os.path.join(x3dna_path, "mutate_bases") + " -l " + os.path.join(output_dir, str(os.getpid()) + ".mutations.dat") + " " + os.path.join(output_dir, "nucleosomal_dna.pdb") + " " + os.path.join(output_dir, "nucleosomal_dna_modified.pdb"))
        process = subprocess.check_output([os.path.join(x3dna_path, "mutate_bases"), "-l", os.path.join(output_dir, str(os.getpid()) + ".mutations.dat"), os.path.join(output_dir, "nucleosomal_dna.pdb"), os.path.join(output_dir, "nucleosomal_dna_modified.pdb")])
        print("MD2P: new dna PDB created in: " + os.path.join(output_dir, "nucleosomal_dna_modified.pdb") + "\n\n")
        dna_pdb = PDB(os.path.join(output_dir, "nucleosomal_dna_modified.pdb"))
        histone_number=0
        for chain in nuc_dna.chains:
            if chain.chaintype=="P":
                print("MD2P: add HISTONE chain: " + histone_alphabet[histone_number] + "\n\n")
                chain.chain=histone_alphabet[histone_number]
                histone_number = histone_number +1
                dna_pdb.add_chain(chain)
        dna_pdb.clean()
        dna_pdb.write(dna_file,force=True)

    elif dna_str == "B":
        print("MD2P run X3DNA command: " + os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -b " + dna_file)
        os.system(os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -b " + dna_file)

    elif dna_str == "A":
        if verbose: print(os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -a " + dna_file)
        process = subprocess.check_output([os.path.join(x3dna_path, "fiber"), "-seq=" + modeling_seq, "-a", dna_file], stderr=subprocess.STDOUT, env=os.environ)
    
    elif dna_str == "C":
        if verbose: print(os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -c " + dna_file)
        process = subprocess.check_output([os.path.join(x3dna_path, "fiber"), "-seq=" + modeling_seq, "-c", dna_file], stderr=subprocess.STDOUT, env=os.environ)

    elif dna_str == "D":
        if verbose: print(os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -d " + dna_file)
        process = subprocess.check_output([os.path.join(x3dna_path, "fiber"), "-seq=" + modeling_seq, "-d", dna_file], stderr=subprocess.STDOUT, env=os.environ)

    elif dna_str == "Z":
        if verbose: print(os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -z " + dna_file)
        process = subprocess.check_output([os.path.join(x3dna_path, "fiber"), "-seq=" + modeling_seq, "-z", dna_file], stderr=subprocess.STDOUT, env=os.environ)

    else:
        #by default use B dna straigth line
        print("MD2P run X3DNA command: " + os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -b " + dna_file)
        os.system(os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -b " + dna_file)

    return dna_file


def get_models(json_dict,job_dir,verbose=True,remote=True):


    job_id     = os.path.basename(job_dir)
    scan_dir   = os.path.join(job_dir,"scan")
    output_dir = os.path.join(scan_dir, "models")
    dummy_dir  = os.path.join(output_dir,"dummy_modeling")
    if not os.path.exists(output_dir): 
        os.mkdir(output_dir)
        os.system("chmod -R 777 %s"%output_dir)

    #get edges to model
    edges2model=json_dict["edges2model"].split(",")
    ppi=set()
    pdi=set()
    if verbose: print("Reading edges in MD2P")
    for edge in edges2model:
        edge_type=None
        for data in edge.split(";"):
            if len(data.split("="))>0:
                if data.split("=")[0] == "source": alfa = data.split("=")[1].replace("binding site ","BS-")
                if data.split("=")[0] == "target": beta = data.split("=")[1].replace("binding site ","BS-")
                if data.split("=")[0] == "id": 
                   interaction = data.split("=")[1]
                   n=0
                   for node in interaction.split("-"):
                       if node.startswith("binding"): n=n+1
                   if n==1: edge_type = "PDI"
                   if n==0: edge_type = "PPI"
        if edge_type is None: continue
        if edge_type == "PDI": 
            if   alfa.startswith("BS-"): pdi.add((alfa,beta))
            elif beta.startswith("BS-"): pdi.add((beta,alfa))
            else: continue
        if edge_type == "PPI": 
            if (alfa,beta) not in ppi and (beta,alfa) not in ppi: ppi.add((alfa,beta))

        #Add forcing homo-dimers and let the modelling decide wether they may occur 
        # comment next 4 lines to apply only to co-factors
        if edge_type == "PDI": 
            if   alfa.startswith("BS-"): ppi.add((beta,beta))
            elif beta.startswith("BS-"): ppi.add((alfa,alfa))
            else: continue
        if edge_type == "PPI": 
            if (alfa,alfa) not in ppi: ppi.add((alfa,alfa))
            if (beta,beta) not in ppi: ppi.add((beta,beta))
        


    #get threading files of TF and TF dimers
    threads_to_model=set()
    for edge in pdi:
      print("Model TF binding: %s %s"%(edge[0], edge[1]))
      start = edge[0].lstrip("BS-")
      tf    = edge[1]
      for bs in json_dict["clusters"]:
          if bs["ID"] != start: continue
          if tf not in bs["protein"].split(";"): continue
          if not bs.has_key("monomer"): continue
          for threads in bs["monomer"]:
            if tf in threads[0]:
                  if verbose: print("\t-- model monomer "+threads[0])
                  threads_to_model.add(threads[0])
          for pairs in bs["dimer"]:
            if tf  in pairs[0][0] or  tf in pairs[1][0]:
                alfa = pairs[0][0].lstrip("aux_files/").split(".")[0]
                beta = pairs[1][0].lstrip("aux_files/").split(".")[0]
                threads_to_model.add(pairs[0][0])
                threads_to_model.add(pairs[1][0])
                if verbose: print("\t-- model dimer "+pairs[0][0]+" "+pairs[1][0])
                if (alfa,beta) not in ppi and (beta,alfa) not in ppi: ppi.add((alfa,beta))
                
    #model TFs
    modeling_file = open(os.path.join(scan_dir,"threading_models.txt"),"w")
    modeling_file.write("#threading files\n")
    for thread_file in threads_to_model:
        modeling_file.write("%s\n"%os.path.join(scan_dir,thread_file))
    modeling_file.close()
    logfile    = os.path.join(scan_dir,"modeling_by_threading.log")
    parameters = " -i " + os.path.join(scan_dir,"threading_models.txt")
    parameters = parameters + " -p " + pdb_dir 
    parameters = parameters + " --dummy=" + dummy_dir 
    parameters = parameters + " -v --parallel --dna "
    parameters = parameters + " --threading "
    parameters = parameters + " -o " + output_dir
    command    = "model_multiple_proteins.py"

    # Execute MODELING in remote#
    if remote:
          print("\nRUN %s ON REMOTE with PARAMETERS %s"%(command,parameters))
          job = functions.execute_in_remote(command,parameters,config,output_dir,logfile,waiting=True)
          print("FINISHED REMOTE JOB: "+job)
          if str(job).startswith("Error"):
              json_dict["error_msg"] =job
              return json_dict
    else:
    # Execute MODELING in local#
          print("\nRUN %s ON LOCAL \n"%(command))
          print("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))
          os.system("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))

    # model dna
    dna_seq = json_dict["dna"]["seq"][0]["full_sequence"]
    if json_dict.has_key("dna_str"):
       dna_str = json_dict["dna_str"]
    else:
       dna_str = "B"
       json_dict["dna_str"] = dna_str

    begining = 1
    ending    = len(dna_seq)
    dna_files=[]
    window  = 250
    overlap = 25
    step    = max( (window - overlap),2)
    number_of_fragments = 1+int(ending)/step
    print("MD2P: NUMBER OF FRAGMENTS <= %d"%number_of_fragments)
    end_interval=0
    try:
      for i in range(number_of_fragments):
         begin_interval        = max((end_interval - overlap),0)
         if begin_interval    >= ending: continue
#         print("MD2P: BUILD PARTIAL DNA OUT OF SIZE "+str(ending)+" IN WINDOWS OF "+str(window))
#         print("MD2P: BUILD PARTIAL DNA FRAGMWENT "+str(i)+" OUT OF "+str(number_of_fragments)+" WITH OVERLAP "+str(overlap))
#         print("MD2P: BUILD PARTIAL DNA BEGIN INTERVAL "+str(begin_interval))
         end_interval          = min((begin_interval + window),ending)
         #if next interval is smaller than twice the window finish with the last fragment
         if (end_interval + window) >=  ending : end_interval=ending
         print("MD2P: BUILD PARTIAL DNA "+str(begin_interval+1)+ " - " + str(end_interval))
         end=end_interval
         begin=begin_interval+1
         dna_file     = build_dna_structure(dna_seq, dna_str, output_dir,  begin, end, verbose=False)
         print("MD2P: STORAGE OF "+dna_file)
         dna_files.append(dna_file)
         if end_interval==ending:break
    except:
      json_dict["error_msg"]="Failed to build DNA structure"
      return json_dict

    #Model large complexes of TF-DNA interactions
    if json_dict["combined_models"]=="true":
       # Get the output files #
       if os.path.isfile(os.path.join(scan_dir, "orthologs_with_best_templates.json")):
        orthologs_file = os.path.join(scan_dir, "orthologs_with_best_templates.json")
       else:
        orthologs_file = None
       print("MD2P: Modelling the complexes TF-DNA by fragments of 250bp or less ")
       logfile    = os.path.join(scan_dir,"complexbuilder.log")
       parameters = " -d " + output_dir
       #Use information of orthologs file to select dimers
       if orthologs_file is not None: 
          parameters = parameters + " -ortho " + orthologs_file
       #Maximum number of iterations is 7!, if more than 7 TF are given the model will be heuristic
       parameters = parameters + " -maxit 5040 " 
       parameters = parameters + " -e 3 " 
       parameters = parameters + " -o " + os.path.join(output_dir,"TF_DNA_FRAGMENTS")
       command    = "complexbuilder.py"

       # Execute MODELING in remote#
       if remote:
          print("\nRUN %s ON REMOTE with PARAMETERS %s"%(command,parameters))
          job = functions.execute_in_remote3(command,parameters,config,output_dir,logfile,waiting=True)
          print("FINISHED REMOTE JOB: "+job)
          if str(job).startswith("Error"):
              json_dict["error_msg"] =job
              return json_dict
       else:
       # Execute MODELING in local#
          print("\nRUN %s ON LOCAL \n"%(command))
          print("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))
          print("\nWARNING PYTHON3 IS MANDATORY \n")
          os.system("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))

       #Optimize the TF-DNA complexes
       logfile    = os.path.join(scan_dir,"optimize_complex.log")
       parameters = " -d " + os.path.join(output_dir,"TF_DNA_FRAGMENTS")
       parameters = parameters + " --parallel "
       command    = "optimize_complex.py"

       # Execute OPTIMIZE in remote#
       if remote:
          print("\nRUN %s ON REMOTE with PARAMETERS %s"%(command,parameters))
          job = functions.execute_in_remote(command,parameters,config,output_dir,logfile,waiting=True)
          print("FINISHED REMOTE JOB: "+job)
          if str(job).startswith("Error"):
              json_dict["error_msg"] =job
              return json_dict
       else:
       # Execute OPTIMIZE in local#
          print("\nRUN %s ON LOCAL \n"%(command))
          print("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))
          os.system("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))


    # model PPIs
    sequences = os.path.join(job_dir,"pwm_database", "sequences","sequences.fasta")
    files_dir = config.get("Paths", "files_path")
    TcoFseq   = config.get("Paths", "TcoFseq")
    TcoF_seq  = os.path.join(files_dir, TcoFseq)
    ppi_fasta = os.path.join(output_dir,"ppi.fa")
    try:
        if os.path.exists(sequences) and os.path.isfile(sequences): os.system("cat %s > %s"%(sequences,ppi_fasta))
        os.system("cat %s >> %s"%(TcoF_seq,ppi_fasta))
    except:
      json_dict["error_msg"]="Failed to include co-factors"
      return json_dict
    ppi_input = open(os.path.join(output_dir,"ppi.ppi"),"w")
    ppi_input.write("# PPIs to  model\n")
    for pp in ppi:
        ppi_input.write("%s\t%s\n"%(pp[0],pp[1]))
    ppi_input.close()
    logfile    = os.path.join(scan_dir,"modpin.log")
    parameters = " -i " + os.path.join(output_dir,"ppi.ppi")
    parameters = parameters + " -seq " + ppi_fasta
    parameters = parameters + " --dummy=" + dummy_dir 
    parameters = parameters + " -o " + output_dir
    parameters = parameters + " --complete 0.95 "
    parameters = parameters + " -v --parallel --hydrogens --renumerate "
    command    = "execute_modpin.py"
    # Execute MODPIN in remote#
    if remote: 
          print("\nRUN %s ON REMOTE"%command)
          job = functions.execute_in_remote(command,parameters,config,output_dir,logfile,waiting=True)
          print("FINISHED REMOTE JOB: "+job)
          if str(job).startswith("Error"):
              json_dict["error_msg"] =job
              return json_dict
    else:
    # Execute MODPIN in local#
          print("\nRUN %s ON LOCAL \n"%(command))
          print("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))
          os.system("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))

    # Compact models into a TAR file
    pairs_dir  = os.path.join(output_dir,"pair_interactions")
    if not os.path.exists(pairs_dir): 
        os.mkdir(pairs_dir)
        os.system("chmod -R 777 %s"%pairs_dir)
    n_models=0
    for pdb_files in os.listdir(output_dir):
        if not pdb_files.endswith(".pdb"): continue
        n_models = n_models + 1
        shutil.copy(os.path.join(output_dir,pdb_files),os.path.join(pairs_dir,pdb_files))
    if n_models==0: json_dict["error_msg"]="Failed to build TF-DNA interaction models"

    ppi_input = open(os.path.join(output_dir,"ppi.ppi"),"r")
    for line in ppi_input:
        if line.startswith("#"):continue
        pa,pb=line.strip().split()
        ppi_dir=pa+"::"+pb
        if not os.path.isdir(os.path.join(output_dir,ppi_dir)):continue
        for pdb_files in os.listdir(os.path.join(output_dir,ppi_dir)):
            if not pdb_files.endswith(".pdb"): continue
            pdb_name = ppi_dir+"_"+pdb_files
            shutil.copy(os.path.join(output_dir,ppi_dir,pdb_files),os.path.join(pairs_dir,pdb_name))
    
    model_file = os.path.join(output_dir, job_id +"_models.tar.gz")

    print("MD2P: COMPRESSED FILE WITH BINARY INTERACTION MODELS ... "+model_file)
    os.system("tar -cv --gzip --directory=%s -f %s %s "%(os.path.dirname(pairs_dir), model_file,os.path.basename(pairs_dir)))

    print("MD2P: Modify dictionary to include the file and remove root path")
    json_dict["tar_file_models"] = model_file.replace(root_path,"") 

    if json_dict["combined_models"]=="true" and json_dict["dna_str"]!="B_bent":
       ############################
       # Combine models with IMP
       ############################
       #Create the input files for IMP
       print("MD2P: Models with IMP ... "+job_id)
       print("MD2P: Input files in Integrative_Modelling ")
       imp_folder = os.path.join(output_dir,"Integrative_Modelling")
       logfile    = os.path.join(scan_dir,"regulatory_complex.log")
       #clean folder
       if os.path.exists(imp_folder):
          shutil.rmtree(imp_folder)
       #get parameters
       parameters = " -i " + output_dir
       parameters = parameters + " -d " + imp_folder
       parameters = parameters + " --dummy "+ os.path.join(scan_dir,"dummy_regulatory_complex")
       parameters = parameters + " --max_iterations 1 " 
       parameters = parameters + " --max_branches 1 " 
       parameters = parameters + " --trim_restraints 0.2 " 
       parameters = parameters + " -o model " 
       command    = "regulatory_complex.py"
       # Execute MODELING in remote#
       if remote:
          print("\nRUN %s ON REMOTE with PARAMETERS %s"%(command,parameters))
          job = functions.execute_in_remote(command,parameters,config,output_dir,logfile,waiting=True)
          print("FINISHED REMOTE JOB: "+job)
          if str(job).startswith("Error"):
              json_dict["error_msg"] =job
              return json_dict
       else:
       # Execute MODELING in local#
          print("\nRUN %s ON LOCAL \n"%(command))
          print("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))
          os.system("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))
       # Execute multiple IMP modelings (limitted to 3 models)
       print("MD2P: Multiple models ")
       files_dir  = config.get("Paths", "files_path")
       imp_script = os.path.join(files_dir,"model_IMP.py")
       logfile    = os.path.join(scan_dir,"IMP_modelling.log")
       parameters = " -i " + imp_folder
       parameters = parameters + " -o " + os.path.join(imp_folder,"IMP")
       parameters = parameters + " --dummy "+ os.path.join(scan_dir,"dummy_IMP")
       parameters = parameters + " --parallel " 
       parameters = parameters + " --verbose " 
       parameters = parameters + " --max_models 3 " 
       parameters = parameters + " --info " +  os.path.join(imp_folder,"IMP_execution.log")
       command    = "model_multiple_IMP.py"
       # Execute MODELING in remote#
       if remote:
          print("\nRUN %s ON REMOTE with PARAMETERS %s"%(command,parameters))
          job = functions.execute_in_remote(command,parameters,config,output_dir,logfile,waiting=True)
          print("FINISHED REMOTE JOB: "+job)
          if str(job).startswith("Error"):
              json_dict["error_msg"] =job
              return json_dict
       else:
       # Execute MODELING in local#
          print("\nRUN %s ON LOCAL \n"%(command))
          print("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))
          os.system("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))

    if json_dict["combined_models"]=="true":
       imp_script = os.path.join(files_dir,"model_IMP.py")
       print("MD2P: copy %s %s "%(imp_script,os.path.join(imp_folder,"model_IMP.py")))
       shutil.copy(imp_script,os.path.join(imp_folder,"model_IMP.py"))
       logfiles= [x for x in os.listdir(imp_folder) if x.endswith("log")]
       for logfile in logfiles:
           os.rename(os.path.join(imp_folder,logfile),os.path.join(scan_dir,logfile))
       model_file = os.path.join(output_dir, job_id +"_complexes.tar.gz")
       print("MD2P: COMPRESSED FILE WITH FRAGMENTS OF COMPLEXES ... "+model_file)
       os.system("tar -cv --gzip --directory=%s -f %s %s %s "%(output_dir, model_file,"TF_DNA_FRAGMENTS","Integrative_Modelling"))
       print("MD2P: Modify dictionary to include the file and remove root path")
       json_dict["tar_file_complex"] = model_file.replace(root_path,"")

    print("MD2P: RESULTING DICTIONARY WITH MODELS")
    print(json_dict)

    return json_dict


def handle_networks(full_json_dict, output_dir, dummy_dir,remote=False):

    # Create input file for create_network.py #
    
    
    # Parse data from json dictionary #
    input_prot_dict = {}
    for bs in full_json_dict["clusters"]:
        key = (str(bs["ID"]) , str(bs["hit_name"]),  str(bs["start"]) , str(bs["end"]) )
        uni_ids = bs["protein"].split(";")
        input_prot_dict[key] = uni_ids

    # Execute create_network.py #

    # 1) Make input file  #
    job_id = full_json_dict["_id"]
    network_input = os.path.join(dummy_dir, job_id + "_network_input.txt")
    ni = open(network_input, "w")
    for bs in input_prot_dict.keys():
        ni.write("binding site " + bs[0] + "\t" + bs[1] + "\t" + bs[2] + "\t" + bs[3] + "\n")
        ni.write("\t".join(input_prot_dict[bs]) + "\n\n")
    ni.close()

    # 2) Define output
    network_output = os.path.join(dummy_dir, job_id + "_network_output.json")

    # 3) Get database of co-factors
    files_dir = config.get("Paths", "files_path")
    TcoF      = config.get("Paths", "TcoF")
    database_file = os.path.join(files_dir, TcoF)

    # parameters
    parameters = " -i " +  network_input
    parameters = parameters + " -d " + database_file
    parameters = parameters +  " -o " + network_output
    if full_json_dict["cofactors"] == "true": parameters = parameters +  " -c "
    logfile    = os.path.join(output_dir,"network.log")
    command    = "create_network.py"

        # Execute PWM in remote#
    if remote:
          print("RUN %s ON REMOTE"%command)
          job = functions.execute_in_remote(command,parameters,config,output_dir,logfile,waiting=True)
          print("FINISHED REMOTE JOB: "+job)
          if str(job).startswith("Error"):
              json_dict["error_msg"] =job
              return json_dict
    else:
    # Execute PWM in local#
          print("RUN %s ON LOCAL \n"%(command))
          print("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))
          os.system("%s %s/%s %s  \n"%(python,scripts_path,command,parameters))


    # Read the output json #
    print("network output is: " + str(os.path.exists(network_output)))
    network_data = readJSON(network_output)
    print("network_data: " + str(network_data))

    return network_data


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
    parser.add_option("-m", "--modeling", default=False, action="store_true", dest="modeling",help="Get the models selected (default=False)")
    parser.add_option("-r", "--remote", default=False, action="store_true", dest="remote",help="Submit to remote cluster (default=False)")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")
    
    (options, args) = parser.parse_args()

    if options.json_file == None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options


#---------------#
# DNA 2 protein #
#---------------#
def dna_2_protein(json_dict, output_dir, verbose=True, remote=True):

    print("\n\n" + str(json_dict) + "\n\n")
    #get data
    dna = json_dict["dna"]["seq"][0]
    fasta_file = os.path.join(output_dir, "DNA", dna["file"])
    # Create directories #
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)
    dummy_dir = os.path.join(output_dir, "dummy")
    if os.path.exists(dummy_dir) == False:
        os.mkdir(dummy_dir)
    scan_dir = os.path.join(output_dir, "scan")
    if os.path.exists(scan_dir) == False:
        os.mkdir(scan_dir)
    #Make available the output to everybody
    print("MD2P: OUTPUT_DIR "+output_dir)
    os.system("chmod -R 777 %s"%output_dir)
    print("MD2P: DUMMY_DIR "+dummy_dir)
    os.system("chmod -R 777 %s"%dummy_dir)
    print("MD2P: SCAN_DIR "+scan_dir)
    os.system("chmod -R 777 %s"%scan_dir)
    print("MD2P: SPECIE "+json_dict["specie"])
    # Get the PWM database to perform the scanning #
    json_dict, database_external,database_dir, database_file, association_file, homologs_dir,error_msg = get_pwm_database(json_dict, output_dir, dummy_dir, remote)
    print("MD2P pwm_database: " )
    print(database_file)
    print("MD2P database_dir: ")
    print(database_dir)
    print("MD2P homologs_dir: ")
    print(homologs_dir)
    print("MD2P association_file: ")
    print(association_file)
    if error_msg is not None: 
        json_dict["error_msg"] = error_msg
        return json_dict

    # Execute scan.py #
    json_dict = run_scan( json_dict, output_dir, database_external, database_dir, database_file, association_file, homologs_dir , scan_dir,  dummy_dir, remote)

    #Make available the output to everybody
    os.system("chmod -R 777 %s"%output_dir)

    # Parse output files #
    print("json_dict: " + str(json_dict))    
    full_json_dict = process_scan_outputs(json_dict,  output_dir, dummy_dir)

    # Create networks #
    remote_network=False
    print("HANDLE NETWORK %s"%os.path.basename(output_dir))
    network_data = handle_networks(full_json_dict, output_dir, dummy_dir, remote_network)

    #Make available the output to everybody
    os.system("chmod -R 777 %s"%output_dir)

    #store the network in the dictionary for the output
    full_json_dict["network"] = network_data


    return full_json_dict


#---------------#
# Main          #
#---------------#


if __name__ == "__main__":


    # Arguments & Options #
    options = parse_options()

    print("Runing MD2P")
    
    # Execute protein 2 DNA #
    if options.modeling:
       print("Runing MODELING MD2P")
       json_dict = get_models(json_dict=json.load(open(os.path.abspath(options.json_file))),job_dir=os.path.abspath(options.output_dir), verbose=options.verbose,remote=options.remote)
    else:
       json_dict = dna_2_protein(json_dict=json.load(open(os.path.abspath(options.json_file))), output_dir=os.path.abspath(options.output_dir),  verbose=options.verbose, remote=options.remote)
    # Save json_dict #
    json_file = os.path.join(options.output_dir,os.path.basename(options.output_dir)+"_full.json")
    if json_dict != None:
        # Create json file #
        out = open(json_file, "wt")
        out.write(json.dumps(json_dict, separators=(',', ':'), indent=2))
        out.close()


