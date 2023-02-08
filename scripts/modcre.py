import json
import os
import sys
import time
import ConfigParser
import bottle
import shutil
import hashlib

from lib import ExperimentFabric, Experiment
from lib import text2file
from lib import readJSON, dumpJSON
#from modcrelib.scripts.mp2d import protein_2_dna
#from modcrelib.scripts.mp2d import protein_2_dna_update
#from modcrelib import dna_2_protein
#from modcrelib.scripts.mp2d import tf_dna_modeling
#from modcrelib import get_tf_dna_models


app = bottle.Bottle()

# Read configuration file #
scripts_path = "/interactomix/sbi/modcrebackend/modcrelib/scripts/"
sys.path.append(scripts_path)
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Imports jbonet's module #
from SBI.structure          import PDB
from SBI.structure.chain    import Chain
from SBI.structure.chain    import ChainOfProtein
from SBI.data               import aminoacids1to3
from SBI.external.blast     import blast_parser

import mp2d
import md2p

#root_path = "/interactomix/sbi/modcrebackend/result/"
root_path = config.get("Paths", "root_path")


routes = {
    'submit'  : '/submit',
    'update'  : '/update',
    'download': '/download',
    'job'     : '/job',
    'session' : '/session',
    'files'   : '/results'
}


def copy_experiment(exp0):

            if os.path.exists(exp0.main_dir):
                # Handle running jobs #
                job_id_o = exp0.__dict__["_id"]
                job_time = get_time()
                if not  exp0.__dict__.has_key("job_name"):
                    name = "NEW"
                else:
                    name = exp0.__dict__["job_name"]
                n_id=0
                print("UPDATE: Check results folder "+os.path.dirname(exp0.main_dir))
                for folder in os.listdir(os.path.dirname(exp0.main_dir)):
                    print("UPDATE: file/folder "+folder)
                    if job_id_o in folder:
                        n_id = n_id + 1
                        print("UPDATE: Found updated folder "+str(n_id))
                print("\nJOB ID "+job_id_o)
                print("Update number "+str(n_id))
                name_data=name.split(".")
                name_data.append(str(n_id))
                name=".".join(name_data)
                job_type = bottle.request.forms.get('type')
                # write the job in the running_jobs_file #
                print("\nread_job: " + str([job_id_o, job_time, name, job_type]))
                #remove_job(job_id, name)
                job_id_data=job_id_o.split(".")
                print job_id_data
                print ("UPDATE: length of JOB_ID file: "+str(len(job_id_data)))
                job_id_data.append(str(n_id))
                job_id=".".join(job_id_data)
                print("\nCONFIRM JOB_ID "+job_id)
                new_dir=os.path.join(os.path.dirname(exp0.main_dir),job_id)
                if not os.path.exists(new_dir):  os.mkdir(new_dir)
                for file in os.listdir(exp0.main_dir):
                       if exp0.__dict__["_id"] in file:
                          os.system("cp -r %s %s"%(os.path.join(exp0.main_dir,file),os.path.join(new_dir,file.replace(exp0.__dict__["_id"],job_id))))
                       else:
                          os.system("cp -r %s %s"%(os.path.join(exp0.main_dir,file),os.path.join(new_dir,file)))
                os.system("chmod -R 777 %s"%new_dir)
            print("\n\nUpdate job "+job_id) 
            exp=Experiment(job_id)
            for key in exp0.__dict__.keys():
                try:
                  if key=="motifs":
                     exp.__dict__[key]=exp0.__dict__[key][:]
                     for k in range(0, len(exp0.__dict__[key])):
                       exp.__dict__[key][k]=exp0.__dict__[key][k][:]
                       for i in range(0, len(exp0.__dict__[key][k])):
                           motif = exp0.__dict__[key][k][i]
                           for key_motif in motif.keys():
                             try:
                               if key_motif=="pwm_frequencies":
                                 exp.__dict__[key][k][i][key_motif]= exp0.__dict__[key][k][i][key_motif][:]
                               elif key_motif=="logos":
                                 exp.__dict__[key][k][i][key_motif][0]= exp0.__dict__[key][k][i][key_motif][0].replace(job_id_o,job_id)
                                 exp.__dict__[key][k][i][key_motif][1]= exp0.__dict__[key][k][i][key_motif][1].replace(job_id_o,job_id)
                               else:
                                 exp.__dict__[key][k][i][key_motif]= exp0.__dict__[key][k][i][key_motif].replace(job_id_o,job_id)
                             except:
                               print("Unmodified key "+key)
                  elif key=="network":
                       print("Copy network")
                       exp.__dict__[key]=exp0.__dict__[key][:]
                  elif key=="clusters":
                       print("Copy clusters")
                       exp.__dict__[key]=exp0.__dict__[key][:]
                  else:
                       print("Copy KEY: "+key+ " => ")
                       print(str(exp0.__dict__[key]))
                       exp.__dict__[key]=exp0.__dict__[key].replace(job_id_o,job_id)
                except:
                   print("Unmodified key "+key)
            exp.__dict__["job_name"]=name

            return exp


def filter_input_json(input_dict, mode, output_dir, job_id):

    input_file = os.path.join(output_dir, job_id + "_input_data.txt")
    input_json_file = os.path.join(output_dir, input_dict["input"]["file"])

    in_file = open(input_file, "w")

    if (mode == "protein2dna") or (mode == "protein2dna_update"):
        # Get the protein fasta sequence #
        if input_json_file.endswith("input.fa") and input_dict['input']['type'] == 'FASTA':
         try:
          header = open(os.path.join(output_dir, "input.fa"), "r").readlines()[0].rstrip()
          seq = ""
          for line in open(os.path.join(output_dir, "input.fa"), "r").readlines()[1:]:
            seq += line.rstrip()
          seq = seq.replace("\n", "")
         except:
          header="EMPTY"
          seq=""
        if input_json_file.endswith("input.pdb") and input_dict['input']['type'] == 'PDB':
         try:
          header="PDB_INPUT"
          mp2d.fix_pdb(os.path.join(output_dir,"input.pdb"),os.path.join(output_dir,"input.fix.pdb"),output_dir)
          pdb_obj = PDB(os.path.join(output_dir,"input.fix.pdb"))
          pdb_obj.clean()
          shutil.move(os.path.join(output_dir,"input.pdb"),os.path.join(output_dir,"input.original.pdb"))
          new_pdb = PDB()
          seq=""
          chains=[]
          for pdb_chain_obj in pdb_obj.chains:
              chains.append(pdb_chain_obj.chain)
          for chain_id in chains:
              pdb_chain_obj = pdb_obj.get_chain_by_id(chain_id)
              if pdb_chain_obj.chaintype =='P':
                 seq= seq + str(pdb_chain_obj.gapped_protein_sequence.replace("-",""))
              if pdb_chain_obj.chaintype =='N':
                 pdb_chain_obj.renumerate_residues(init=1)
              new_pdb.add_chain(pdb_chain_obj)
          new_pdb.write(os.path.join(output_dir,"input.pdb"))
         except:
          header="EMPTY"
          seq=""

        in_file.write("Input protein: " + header + "\n")
        in_file.write("Input sequence: " + seq + "\n")
        #in_file.write("Specie: " + input_dict["specie"] + "\n")
        in_file.write("Templates: " + input_dict["templates"] + "\n")
        in_file.write("Models: " + input_dict["models"] + "\n")
        in_file.write("Monomers: " + input_dict["monomers"] + "\n")
        in_file.write("Dimers: " + input_dict["dimers"] + "\n")
        if len(input_dict["dna"]["seq"]) > 0:
            dna_ids = []
            dna_seq = []
            for seq in input_dict["dna"]["seq"]:
                dna_ids.append(seq["_id"])
                dna_seq.append(seq["full_sequence"])
            in_file.write("DNA IDs: " + ";".join(dna_ids) + "\n")
            in_file.write("DNA sequences: " + ";".join(dna_seq) + "\n")


    elif (mode == "dna2protein") or (mode == "dna2protein_modeling"):
        # Get the DNA sequence #
        in_file.write("Input dna: " + input_dict["dna"]["seq"][0]["_id"] + "\n")
        in_file.write("Input sequence: " + input_dict["dna"]["seq"][0]["full_sequence"] + "\n")
        if input_dict.has_key("specie"):
           if input_dict["specie"] is not None: 
               in_file.write("Specie: " + input_dict["specie"] + "\n")
           else:
               in_file.write("Specie: defined as None\n")
        else:
           in_file.write("Specie: Not defined\n") 
        if input_dict.has_key("family"):
           if input_dict["family"] is not None: in_file.write("Family: " + input_dict["family"] + "\n")
        if str(input_dict["fimo_threshold"]) == "null":
            in_file.write("FIMO threshold: 0.001" + "\n")
        else:
            in_file.write("FIMO threshold: " + str(input_dict["fimo_threshold"]) + "\n")
        in_file.write("Using default PWMs: " + str(input_dict["default_pwm_db"]) + "\n")
        if input_dict.has_key("cisbp_pwm_db"):
           if input_dict["cisbp_pwm_db"] is not None: in_file.write("Using Cis-BP PWMs: " + str(input_dict["cisbp_pwm_db"]) + "\n")
        if input_dict.has_key("jaspar_pwm_db"):
           if input_dict["jaspar_pwm_db"] is not None: in_file.write("Using JASPAR PWMs: " + str(input_dict["jaspar_pwm_db"]) + "\n")
        if input_dict.has_key("custom_proteins"):
           if input_dict["custom_proteins"] is not None: in_file.write("Using UniProt " + str(input_dict["custom_proteins"]) + "\n")
        if input_dict.has_key("prot"):
           if len(input_dict["prot"]["seq"])>0:
               for data in input_dict["prot"]["seq"]:
                    name=data["_id"].lstrip(">")
                    seq =data["full_sequence"]
                    in_file.write("User provided TF:  " + name + " SEQUENCE: " + seq + "\n")
        if input_dict.has_key("tar_file_models"):
           if  input_dict["tar_file_models"] is not None: in_file.write("Binary models: " + input_dict["tar_file_models"] + "\n")
        if input_dict.has_key("combined_models"):
            in_file.write("Combined models: " + input_dict["combined_models"] + "\n")
        if input_dict.has_key("dna_str"):
            in_file.write("DNA structure: " + input_dict["dna_str"] + "\n")
        if input_dict.has_key("output_type"):
            in_file.write("Output type: " + input_dict["output_type"] + "\n")
        if input_dict.has_key("edges2model"):
            in_file.write("Interactions to model: \n")
            ppi=set()
            pdi=set()
            for edge in input_dict.get("edges2model").split(","):
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
            in_file.write("\tPDI: ")
            for pair in pdi:
                alfa,beta=pair
                in_file.write("\t%s::%s"%(alfa,beta))
            if len(ppi)>0:
               in_file.write("\n\tPPI: ")
               for pair in ppi:
                  alfa,beta=pair
                  in_file.write("\t%s::%s"%(alfa,beta))
            in_file.write("\n")

    in_file.close()
    
    return input_file


def get_time():

    from datetime import datetime
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y  %H:%M:%S")
    return dt_string


def erase_results(result_dir):

    # First remove the directories #
    for dir in os.listdir(result_dir):
        if os.path.isdir(os.path.join(result_dir, dir)):
            os.system("rm -rf " + os.path.join(result_dir, dir))
        if dir.endswith(".txt"):
            os.system("rm -f " + os.path.join(result_dir, dir))
    # Then, create empty files #
    for file in ["finished_jobs.txt", "running_jobs.txt", "failed_jobs.txt"]:
        complete_file = os.path.join(result_dir, file)
        op = open(complete_file, "w")
        op.close()

    return 


#
# CONFIGURATION
#
@app.hook('after_request')
def enable_cors():
    """
    You need to add some headers to each request.
    Don't use the wildcard '*' for Access-Control-Allow-Origin in production.
    """
    bottle.response.headers['Access-Control-Allow-Origin']  = '*'
    bottle.response.headers['Access-Control-Allow-Methods'] = 'PUT, GET, POST, DELETE, OPTIONS'
    bottle.response.headers['Access-Control-Allow-Headers'] = 'Origin, Accept, Content-Type, X-Requested-With, X-CSRF-Token'

#
# PATHS
#
@app.get('/')
def api_routes():
    return routes

@app.route(routes['submit'], method=['POST', 'OPTIONS'])
def api_submit():
    print("Running the api_submit function!!!")
    print("ROOT IS "+root_path)
    if bottle.request.method == "OPTIONS":
        pass
    if bottle.request.method == "POST":
        from lib import ExperimentFabric, Experiment
        #-----------------------#
        #      Protein2DNA      #
        #-----------------------#
        if bottle.request.forms.get('type') == "protein2dna":
            #update_jobs()
            # Create an experiment object #
            exp = ExperimentFabric.build(
                bottle.request.environ.get('REMOTE_ADDR'),
                bottle.request.forms.get('type'),
                {
                    'type': bottle.request.forms.get('input'),
                    'file': None,
                    '_id' : bottle.request.forms.get('inid')
                }, job_id=bottle.request.forms.get('job_id'))
            if exp.input['type'] == 'PDB':
                pdb = bottle.request.files.get('file')
                pdb.save(exp.file_path('input.pdb'))
                exp.input['file'] = 'input.pdb'
            elif exp.input['type'] == 'FASTA':
                fasta = bottle.request.forms.get('fasta')
                text2file(fasta, exp.file_path('input.fa'))
                exp.input['file'] = 'input.fa'
            elif exp.input['type'] == 'threading':
                exp.input['file'] = bottle.request.forms.get('threading_file')
            # Iterate over DNA sequences #
            for x in bottle.request.forms.get('dna').split(','):
                if len(x) > 3:
                    exp.add_dna(x)
            # Create a DNA file containing the DNA sequences #
            exp.make_dna_file()
            # Include the options of execution into the experiment #
            exp.include_options(
                templates = bottle.request.forms.get('templates'),
                models = bottle.request.forms.get('models'),
                monomers = bottle.request.forms.get('monomers'),
                dimers = bottle.request.forms.get('dimers'),
                uniprotSearch = bottle.request.forms.get('uniprotSearch')
            )
            job_id = bottle.request.forms.get('job_id')
            job_time = get_time()
            name = bottle.request.forms.get('job_name')
            job_type = bottle.request.forms.get('type')
            # write the job in the running_jobs_file #
            write_send_job(job_id, job_time, name, job_type)
            # Create input json #
            output_dir = (exp.main_dir)
            input_json_file = os.path.join(output_dir, job_id + "_input.json")
            input_file = filter_input_json(exp.__dict__, "protein2dna", output_dir, job_id)
            dumpJSON(input_json_file, exp.__dict__)
            exp.__dict__["job_time"]=job_time
            exp.__dict__["input_json"] = input_file.replace(root_path, "")
            # Execute programs from the server-side 
            print("\n\nbefore protein2dna the experiment dict has " + str(len(exp.__dict__.keys())) + " keys\n")
            print("EXECUTE MP2D IN "+output_dir)
            exp.__dict__ = mp2d.protein_2_dna(exp.__dict__, output_dir, remote=True, verbose=True)
            print("\n\nafter protein2dna the experiment dict has " + str(len(exp.__dict__.keys())) + " keys\n")
            exp.__dict__["job_name"] = name
            # Include the ID into the experiment dict #
            exp.include_id()
            # Create a json file from the experiment dict #
            exp.toJSON()
            # Remove the job from running jobs #
            write_finish_job(job_id, job_time, name, job_type)
            return exp.confirm()
        #-----------------------#
        #   Protein2DNA update  #
        #-----------------------#
        if bottle.request.forms.get('type') == "protein2dna_update":
            # Create an experiment object #
            exp0 = Experiment(bottle.request.forms.get('jobid'))
            #Copy folder into new one
            exp = copy_experiment(exp0)
            exp.restart()
            print("\n\ninid: " + str(bottle.request.forms.get('inid')))
            print("\n\ninid_list: " + str([x for x in bottle.request.forms.get('inid').split(',')]))
            exp.start("protein2dna", 
            {
                'type': bottle.request.forms.get('input_type'),
                'file': bottle.request.forms.get('file'),
                '_id': [x for x in bottle.request.forms.get('inid').split(',')]
            })
            print("\n\nExperiment input: " + str(bottle.request.forms.get('input')))
            job_id = exp.__dict__["_id"]
            print("\n\nrestart the updated job "+job_id)
            for  key in exp0.__dict__.keys():
               print(" -- Key "+key+" : ")
               print exp.__dict__[key]
               print("\n")
            # Remove all files in the output directory related with the DNA sequence #
            if os.path.exists(exp.main_dir):
                # Iterate over DNA sequences #
                for x in bottle.request.forms.get('dna').split(','):
                    if len(x) > 3:
                        print("\n\nDNA in update: " + str(x))
                        exp.add_dna(x)
                # Create a DNA file containing the DNA sequences #
                exp.make_dna_file()
                # Include the input pwms #
                exp.__dict__["input_pwms"] = bottle.request.forms.get('pwms')
                # Handle running jobs #
                job_id = exp.__dict__["_id"]
                print("\n\nUpdate job "+job_id) 
                job_time = get_time()
                if not "job_name" in exp.__dict__:
                    name = "None"
                else:
                    name = exp.__dict__["job_name"]
                job_type = bottle.request.forms.get('type')
                # write the job in the running_jobs_file #
                print("\n\nwrite_send_job UPDATE: " + str([job_id, job_time, name, job_type]))
                write_send_job(job_id, job_time, name, job_type)
                # Create input json #
                output_dir = (exp.main_dir)
                input_json_file = os.path.join(output_dir, job_id + "_input.json")
                input_file = filter_input_json(exp.__dict__, "protein2dna", output_dir, job_id)
                if os.path.exists(input_json_file): shutil.copy(input_json_file,input_json_file+".start")
                input_file = filter_input_json(exp.__dict__, "protein2dna_update", output_dir, job_id)
                dumpJSON(input_json_file, exp.__dict__)
                exp.__dict__["job_time"]=job_time
                exp.__dict__["input_json"] = input_file.replace(root_path, "")
                # Execute programs from the server-side #
                exp.__dict__ = mp2d.protein_2_dna_update(exp.__dict__, output_dir, remote=True, verbose=True)
                # Include the ID into the experiment dict #
                exp.include_id()
                # Create a json file from the experiment dict #
                exp.toJSON()
                write_finish_job(job_id, job_time, name, job_type)
                return exp.confirm()
            else:
                return None
        #-----------------------#
        #      DNA2protein      #
        #-----------------------#
        if bottle.request.forms.get('type') == "dna2protein":
            # Create an experiment object #
            #update_jobs()
            #if str(bottle.request.forms.get('modeling')) == "false":
            exp = ExperimentFabric.build(
                bottle.request.environ.get('REMOTE_ADDR'),
                bottle.request.forms.get('type'),
                {
                    'type': bottle.request.forms.get('input'),
                    'file': None,
                    '_id' : bottle.request.forms.get('inid')
                }, job_id=bottle.request.forms.get('job_id'))
            #### Get DNA sequence
            exp.input['type'] = 'FASTA'
            fasta = bottle.request.forms.get('fasta')
            text2file(fasta, exp.file_path('input.fa'))
            exp.input['file'] = 'input.fa'
            ### Iterate over DNA sequences #
            for x in bottle.request.forms.get('dna').split(','):
                if len(x) > 3:
                    exp.add_dna(x)
            ### Create a DNA file containing the DNA sequences #
            exp.make_dna_file()
            ### Iterate over TF sequences #
            for x in bottle.request.forms.get('prot').split(','):
                if len(x) > 3:
                    exp.add_prot(x)
            ### Create a PROTEINS file containing the TF sequences #
            exp.make_prot_file()
            ####
            # Include the options of execution into the experiment #
            exp.include_options(
                specie = bottle.request.forms.get('specie'),
                family = bottle.request.forms.get('family'),
                fimo_threshold = bottle.request.forms.get('fimo_threshold'), 
                default_pwm_db = bottle.request.forms.get('default_pwm_db'), 
                jaspar_pwm_db = bottle.request.forms.get('jaspar_pwm_db'),
                cisbp_pwm_db = bottle.request.forms.get('cisbp_pwm_db'),
                output_type = bottle.request.forms.get('output_type'),
                cofactors  =  bottle.request.forms.get('cofactors'),
                all_orthologs  =  bottle.request.forms.get('all_orthologs')
                )
            # Handle running jobs #
            job_id = bottle.request.forms.get('job_id')
            job_time = get_time()
            name = bottle.request.forms.get('job_name')
            job_type = bottle.request.forms.get('type')
            # write the job in the running_jobs_file #
            write_send_job(job_id, job_time, name, job_type)
            # Create input json #
            output_dir = (exp.main_dir)
            print("output_dir: " + output_dir)
            print("job_id: " + job_id)
            print("job_time: " + str(job_time))
            #exit(0)
            if os.path.exists(output_dir) == False:
               os.mkdir(output_dir)
            #Make available the output directory to everybody
            os.system("chmod -R 777 %s"%output_dir)
            print("MODCRE: CHMOD 777 of "+output_dir)

            input_json_file = os.path.join(output_dir, job_id + "_input.json")
            dumpJSON(input_json_file, exp.__dict__)
            exp.__dict__["job_time"]=job_time
            
            # Execute programs from the server-side #
            print("\n\n\n\nbefore dna2protein the experiment dict has " + str(len(exp.__dict__.keys())) + " keys\n\n\n\n")
            exp.__dict__["show_models"] = False
            full_json_dict = md2p.dna_2_protein(exp.__dict__, output_dir,remote=True, verbose=True)
            print("\n\n\n\nafter dna2protein the experiment dict has " + str(len(exp.__dict__.keys())) + " keys\n\n\n\n")
            full_json_dict["job_name"] = name
            full_json_dict['jobid'] = exp._id

            input_file = filter_input_json(exp.__dict__, "dna2protein", output_dir, job_id)
            full_json_dict["input_json"] = input_file.replace(root_path, "")
            job_id = os.path.basename(output_dir)
            full_json_file = os.path.join(output_dir, job_id + "_full.json")
            dumpJSON(full_json_file, full_json_dict)

            write_finish_job(job_id, job_time, name, job_type)

            return exp.confirm()
        #-----------------------#
        #        Examples       #
        #-----------------------#
        if bottle.request.forms.get('type') == "example":
            # Select the json file for the example #
            examples_dir = config.get("Paths", "example_dir")
            example = bottle.request.forms.get('example')

            if example == "TF; 1 model":
                ex_dir = os.path.join(examples_dir, "p2d_1")
                ex_type = "protein2dna"

            if example == "TF; >1 models":
                ex_dir = os.path.join(examples_dir, "p2d_2")
                ex_type = "protein2dna"

            if example == "TF; >1 DBD":
                ex_dir = os.path.join(examples_dir, "p2d_3")
                ex_type = "protein2dna"

            if example == "TF; >1 models; >1 DBD":
                ex_dir = os.path.join(examples_dir, "p2d_4")
                ex_type = "protein2dna"

            if example == "TF; Heterodimeric":
                ex_dir = os.path.join(examples_dir, "p2d_5")
                ex_type = "protein2dna"

            if example == "TF + DNA":
                ex_dir = os.path.join(examples_dir, "p2d_6")
                ex_type = "protein2dna"

            if example == "TF + 2DNA":
                ex_dir = os.path.join(examples_dir, "p2d_7")
                ex_type = "protein2dna"

            if example == "Interferon beta enhancer":
                ex_dir = os.path.join(examples_dir, "d2p_1")
                ex_type = "dna2protein"

            if example == "ZRS sonic hedgehog enhancer":
                ex_dir = os.path.join(examples_dir, "d2p_2")
                ex_type = "dna2protein"

            exp = ExperimentFabric.build(
                bottle.request.environ.get('REMOTE_ADDR'),
                ex_type,
                {
                    'type': None,
                    'file': None,
                    '_id' : None
                },
                job_id=bottle.request.forms.get('job_id'))
            # Get the main dir of the experiment #
            for file in os.listdir(ex_dir):
                if file.endswith("_full.json"):
                    dir_id = file.replace("_full.json", "")
                    break
            exp._id = dir_id
            # Remove the previous content from the experiment directory #
            print("rm -rf " + os.path.join(exp.main_dir, "*"))
            os.system("rm -rf " + os.path.join(exp.main_dir, "*"))
            # If the destination directory doesn't exists, create a new one #
            if os.path.exists(exp.main_dir) == False:
                print("mkdir " + exp.main_dir)
                os.system("mkdir " + exp.main_dir)
            # Copy the content of the example directory to the experiment directory #
            print("cp -r " + os.path.join(ex_dir, "*") + " " + exp.main_dir)
            os.system("cp -r " + os.path.join(ex_dir, "*") + " " + exp.main_dir)
            # Get the json file from the experiment directory #
            for file in os.listdir(exp.main_dir):
                if file.endswith("_full.json"):
                    json_file = os.path.join(exp.main_dir, file)
                    break
            # Read the json file and return it #
            #print(json_file)
            return readJSON(json_file)

        #-----------------------#
        #    Modeling  (P2D)    #
        #-----------------------#
        if bottle.request.forms.get('type') == "modelDNA":
            model_dict = {}
            protein = bottle.request.forms.get('protein')
            dna = bottle.request.forms.get('dna')
            print("MODEL FORM "+protein+" DNA "+dna)
            model_dict["file"] = mp2d.tf_dna_modeling(protein, dna)
            
            return model_dict

        #-----------------------#
        #    Modeling  (D2P)    #
        #-----------------------#
        
        if bottle.request.forms.get('type') == "dna2protein_modeling":

            # Create an experiment object #
            print("\n JOB ID "+bottle.request.forms.get('jobid'))
            print("\n JOB NaMe "+bottle.request.forms.get('job_name'))

            exp0 = Experiment(bottle.request.forms.get('jobid'))
            #Copy folder into new one
            print("\n JOB_NAME "+exp0.__dict__["job_name"])
            exp = copy_experiment(exp0)
            print("\n\nNEW DATA")
            print exp.__dict__
            #exp.restart(remove_dna=False)
            print("\n\ninid: " + str(bottle.request.forms.get('inid')))
            print("\n\ninid_list: " + str([x for x in bottle.request.forms.get('inid').split(',')]))
            exp.start("dna2protein_modeling", 
            {
                'type': bottle.request.forms.get('input_type'),
                'file': bottle.request.forms.get('file'),
                '_id': [x for x in bottle.request.forms.get('inid').split(',')]
            })
            job_id = exp.__dict__["_id"]
            print("\n\nrestart the updated job "+job_id)
            for  key in exp0.__dict__.keys():
               print(" -- Key "+key+" : ")
               print exp.__dict__[key]
               print("\n")

            print("EDGES")
            #print bottle.request.forms.get('edges2model')
            print("NODES")
            #print bottle.request.forms.get('nodes2model')
            print("COMBINED MODELS")
            #print bottle.request.forms.get('combined_models')

            exp.include_options(
                dna_seq=bottle.request.forms.get('dna_seq'),
                dna_str=bottle.request.forms.get('dna_str'),
                individual_models=bottle.request.forms.get('individual_models'),
                combined_models=bottle.request.forms.get('combined_models'),
                edges2model=bottle.request.forms.get('edges2model')
            )

            # Handle running jobs #
            job_type = "dna2protein_modeling"
            job_time = get_time()
            name =  exp.__dict__["job_name"]
            # write the job in the running_jobs_file #
            write_send_job(job_id, job_time, name, job_type)
            # Handle information regarding the models #
            dna_seq = bottle.request.forms.get('dna_seq')
            dna_str = bottle.request.forms.get('dna_str')
            if not  exp.__dict__.has_key("dna_str"):  exp.__dict__["dna_str"]=dna_str
            # Create input json #
            output_dir = os.path.abspath(exp.main_dir)
            input_json_file = os.path.join(output_dir, job_id + "_input.json")
            dumpJSON(input_json_file, exp.__dict__)
            input_file = filter_input_json(exp.__dict__, "dna2protein_modelling", output_dir, job_id)
            exp.__dict__["input_json"] = input_file.replace(root_path, "")
            # rename output full.json as original to restart the job while keeping previous inputs
            json_file = os.path.join(output_dir, job_id + "_full.json")
            if os.path.exists(json_file): shutil.move(json_file,json_file+".original")
            ### get_tf_dna_models should return a dict which will be written as a JSON ###
            print("MD2P: GET MODELS IN "+output_dir)
            exp.__dict__ = md2p.get_models(exp.__dict__,output_dir,verbose=True, remote=True)
            print("MD2P: CURSOR RECOVERED FROM GET MODELS")
                        
            # Rearrange the type for additional modeling runs
            exp.__dict__["type"]="dna2protein_modelling"
            # store the information in input_data
            input_file = filter_input_json(exp.__dict__, "dna2protein", output_dir, job_id)
            exp.__dict__["input_json"] = input_file.replace(root_path, "")

            # Create a json file from the experiment dict #
            print("MD2P: MAKE DICTIONARY OUTPUT FULL.JSON")
            exp.toJSON()
            print(exp.__dict__)
            # Write the job with the finished jobs #
            write_finish_job(job_id, job_time, name, job_type)

            return exp.confirm()

        #-----------------------#
        #   Get cluster data    #
        #-----------------------#
        if bottle.request.forms.get('type') == "get_cluster_data":
            print(str(bottle.request.forms.get('output_dir')))
            print(str(bottle.request.forms.get('cluster_id')))
            cluster_file = os.path.join(bottle.request.forms.get('output_dir'), "jsons", "full_cluster_" + str(bottle.request.forms.get('cluster_id')) + ".json")
            cluster_dict = readJSON(cluster_file)

            return cluster_dict


        #-----------------------#
        #      Load species     #
        #-----------------------#
        if bottle.request.forms.get('type') == "species":
            files_dir = config.get("Paths", "files_path")
            specie_dict = {}
            specie_dict["specie_list"] = set()
            for line in open(os.path.join(files_dir, "species.txt"), "r").readlines():
                specie_dict["specie_list"].add(line.rstrip())
            specie_dict["specie_list"] = sorted(list(specie_dict["specie_list"]))

            return specie_dict

        #-----------------------#
        #      Load families    #
        #-----------------------#
        if bottle.request.forms.get('type') == "families":
            family_dict = {}
            family_dict["family_list"] = []
            pdb_dir = config.get("Paths", "pdb_dir")
            for line in open(os.path.join(pdb_dir, "families.txt"), "r").readlines()[1:]:
                family = line.split(";")[1].rstrip()
                if not family in family_dict["family_list"]:
                    family_dict["family_list"].append(family)
            family_dict["family_list"] = ["All"] + sorted(family_dict["family_list"])

            return family_dict


        #-----------------------#
        #        Get Jobs       #
        #-----------------------#
        if bottle.request.forms.get('type') == "jobs":

            result_dir = config.get("Paths", "result_dir")
            
            ##
            #### In case that you want to remove all the jobs and clean the results directory, execute the next command #
            #### erase_results(result_dir)
            ##

            update_jobs()
            jobs_dict = {}
            jobs_dict["jobs"] = []
            # Open finished_jobs_file #
            if os.path.isfile(os.path.join(result_dir, "finished_jobs.txt")):
               try:
                with open(os.path.join(result_dir, "finished_jobs.txt"), "r") as finished:
                 for line in finished:
                    fields = line.split("\t")
                    #print("finished: " + str(fields))
                    job_obj = {"id":fields[0], "time":fields[1], "name":fields[2], "status": "Finished", "type":fields[3]}
                    jobs_dict["jobs"].append(job_obj)
                finished.close()
               except:
                print("Failed to open finished_jobs.txt")
 
            # Open running_jobs_file #
            if os.path.isfile(os.path.join(result_dir, "running_jobs.txt")):
               try: 
                with open(os.path.join(result_dir, "running_jobs.txt"), "r") as running:
                 for line in running:
                    fields = line.split("\t")
                    #print("running: " + str(fields))
                    job_obj = {"id":fields[0], "time":fields[1], "name":fields[2], "status": "Running", "type":fields[3]}
                    jobs_dict["jobs"].append(job_obj)
                running.close()
               except:
                print("Failed to open running_jobs.txt")

            if os.path.isfile(os.path.join(result_dir, "failed_jobs.txt")):
               try:
                with open(os.path.join(result_dir, "failed_jobs.txt"), "r") as failed:
                 for line in failed:
                    fields = line.split("\t")
                    if len(fields) == 5:
                        job_obj = {"id":fields[0], "time":fields[1], "name":fields[2], "status": "Failed", "type":fields[3]}
                    elif len(fields) == 6:
                        job_obj = {"id":fields[0], "time":fields[1], "name":fields[2], "status": "Failed: " + str(fields[5]), "type":fields[3]}
                    jobs_dict["jobs"].append(job_obj)
                failed.close()
               except:
                print("Failed to open failed_jobs.txt")


            return jobs_dict

        
    return exp.confirm()

@app.get(routes['job'])
def api_job():
    print("Running api_job function!!!")
    search = dict((k, bottle.request.query.getall(k)) for k in bottle.request.query.keys())
    print("search: " + str(search))
    if '_id' not in search: 
        return json.dumps({'success': False, 'error': 'Job _id needed'})
    exp = Experiment(search['_id'][0])
    if not os.path.exists(exp.main_dir): 
        return json.dumps({'success': False, 'error': 'Job {0} does not exist'.format(search['_id'][0])})

    return exp.confirm()


# Supplier for result files
@app.route('/results/:path#.+#', name='results')
def result(path):
    return bottle.static_file(path, root='results')



def write_send_job(job_id, job_time, name, job_type):

    print("you are in write_send_job")
    result_dir = config.get("Paths", "result_dir")
    running_jobs_file = os.path.join(result_dir, "running_jobs.txt")
    if os.path.exists(running_jobs_file):
        running = open(running_jobs_file, "a")
    else:
        running = open(running_jobs_file, "w")
    print("RUNNING JOBS "+str(job_id) + "\t" + str(job_time) + "\t" + str(name) + "\t" + str(job_type) + "\t" + str(time.time()) + "\n")
    running.write(str(job_id) + "\t" + str(job_time) + "\t" + str(name) + "\t" + str(job_type) + "\t" + str(time.time()) + "\n")
    running.close()


def write_finish_job(job_id, job_time, name, job_type):

    print("write_finish_job")
    # Check for errors in the job, and what type of error in case there has been one #
    error = False
    error_type = "Unknown"
    result_dir = config.get("Paths", "result_dir")
    if os.path.isfile(os.path.join(result_dir, job_id, job_id + "_full.json")):
        json_file = os.path.join(result_dir, job_id, job_id + "_full.json")
        json_content = readJSON(json_file)
        if "error_msg" in json_content.keys():
            # This means that there has been a error in this execution #
            print("finished job goes to failed jobs, there is an error key")
            error = True
            error_type = json_content["error_msg"]
        if len(json_content.keys()) == 0:
            # This means that there has been a error in this execution #
            print("finished job goes to failed jobs, the json dict is empty")
            error = True
    else:
        print("finished job goes to failed jobs, there is no json dict")
        error = True
    # If there is no error, write the new file in the finished_jobs file #
    removed_jobs=[]
    if error == False:
        print("job didn't fail")
        finished_jobs_file = os.path.join(result_dir, "finished_jobs.txt")
        finished_jobs_file2 = os.path.join(result_dir, "finished_jobs.txt.previous")
        if not os.path.exists(finished_jobs_file):
            with open(finished_jobs_file, "w") as finished:
               finished.write(str(job_id) + "\t" + str(job_time) + "\t" + str(name) + "\t" + str(job_type) + "\t" + str(time.time()) + "\n")
            finished.close()
        else:
            done_lines = []
            for line in open(finished_jobs_file, "r").readlines():
                    if len(line.split()) == 0: continue
                    line_id = line.split("\t")[0]
                    job_start_time = float(line.split("\t")[-1])
                    if os.path.isfile(os.path.join(result_dir, line_id, line_id + "_full.json")):
                        done_lines.append(line)
                    else:
                        removed_jobs.append(line_id)
            shutil.copy(finished_jobs_file,finished_jobs_file2)
            with open(finished_jobs_file, "w") as finished:
              if str(job_id) not in removed_jobs:
                finished.write(str(job_id) + "\t" + str(job_time) + "\t" + str(name) + "\t" + str(job_type) + "\t" + str(time.time()) + "\n")
              for line in done_lines:
                line_id = line.split("\t")[0]
                if line_id in removed_jobs:continue
                finished.write(line.rstrip() + "\n")
            finished.close()
    # If there has been an error, write the information of the job into failed jobs #
    else:
        print("there was an error")
        failed_jobs_file = os.path.join(result_dir, "failed_jobs.txt")
        if os.path.isfile(failed_jobs_file):
            with open(failed_jobs_file, "a") as failed:
                failed.write(str(job_id) + "\t" + str(job_time) + "\t" + str(name) + "\t" + str(job_type) + "\t" + str(time.time()) + "\t" + error_type + "\n")
            failed.close()
            removed_jobs.append(job_id)
        else:
            with open(failed_jobs_file, "w") as failed:
                failed.write(str(job_id) + "\t" + str(job_time) + "\t" + str(name) + "\t" + str(job_type) + "\t" + str(time.time()) + "\t" + error_type + "\n")
            failed.close()
            removed_jobs.append(job_id)

    # Remove the job from the running_jobs_file #
    running_jobs_file = os.path.join(result_dir, "running_jobs.txt")
    keep_lines = []
    current_time = float(time.time())
    for line in open(running_jobs_file, "r").readlines():
        # First, make sure than the executed job has finished #
        if str(job_id) not in removed_jobs:
           if line.startswith(str(job_id) + "\t" + str(job_time) + "\t" + str(name) + "\t" + str(job_type)):
              continue
        keep_lines.append(line)
    # Re-write the running jobs file #
    if len(keep_lines)>0:
       os.system("rm -f " + running_jobs_file)
       with open(running_jobs_file, "w") as running:
         for line in keep_lines:
             running.write(line.rstrip() + "\n")
       running.close()

def update_jobs():

    from datetime import datetime
    
    print("\nUpdating jobs")
    result_dir = config.get("Paths", "result_dir")
    running_jobs_file = os.path.join(result_dir, "running_jobs.txt")
    running_jobs_file2 = os.path.join(result_dir, "running_jobs.txt.previous")
    finished_jobs_file = os.path.join(result_dir, "finished_jobs.txt")
    finished_jobs_file2 = os.path.join(result_dir, "finished_jobs.txt.previous")
    failed_jobs_file = os.path.join(result_dir, "failed_jobs.txt")
    failed_jobs_file2 = os.path.join(result_dir, "failed_jobs.txt.previous")
    failed_jobs = []
    # Update running jobs #
    keep_lines = []
    finished_lines = []
    current_time = float(time.time())
    if os.path.isfile(running_jobs_file):
        shutil.copy(running_jobs_file,running_jobs_file2)
        with open(running_jobs_file2, "r") as running:
         for line in running:
            #print("Running job line: " + line)
            if len(line.split()) == 0: continue
            # This first filter cleans the running file from files that got stuck there #
            line_id = line.split("\t")[0]
            job_start_time = float(line.split("\t")[-1])
            job_type = line.split("\t")[-2]
            time_difference = float(current_time) - float(job_start_time)
            # If .json file has been created, parse it looking for errors #
            if os.path.isfile(os.path.join(result_dir, line_id, line_id + "_full.json")):
                json_file = os.path.join(result_dir, line_id, line_id + "_full.json")
                try:
                    json_content = readJSON(json_file)
                except:
                    # This means that there has been a error in this file (it may be old file) #
                    print("Error of old file "+line_id)
                    error_msg="FAILED FOLDER INCOMPLETED"
                    failed_jobs.append(line.rstrip()+"\t"+error_msg+"\n") 
                    continue
                if "error_msg" in json_content.keys():
                    # This means that there has been a error in this execution #
                    print("Written error "+line_id)
                    failed_jobs.append(line.rstrip()+"\t"+json_content["error_msg"]+"\n")
                    continue
                elif len(json_content.keys()) == 0:
                    # This means that there has been a error in this execution #
                    print("empty json "+line_id)
                    failed_jobs.append(line.rstrip()+"\t"+"Wrong input\n")
                    continue
                # The initial number of keys of the json dict for protein2dna is not known yet #
                elif job_type == "protein2dna" and len(json_content.keys()) <= 15:
                    print("unfinished json "+line_id)
                    failed_jobs.append(line.rstrip()+"\t"+"Wrong input\n")
                    continue
                elif job_type == "dna2protein" and len(json_content.keys()) <= 27:
                    print("unfinished json "+line_id)
                    failed_jobs.append(line)
                    continue
                else:
                    # This means that the job finished but that couldn't be included into the finished jobs section #
                    print("finished job "+line_id)
                    finished_lines.append(line)
                    continue
            # Check execution walltime, (3 or) 30 x 24 hours = 3600x30x24=2592000=(1month), 3600x5x24=432000=(5days) 3600x2x24=172800(2days)
            #if (time_difference > float(2592000) and job_type == "protein2dna") or (time_difference > float(25920000) and job_type == "dna2protein"):
            if time_difference > float(432000):
                print("job walltime")
                failed_jobs.append(line.rstrip()+"\t"+"Exceeding time\n")
                continue
            keep_lines.append(line)
        # Rewrite the running jobs #
        running.close()
        os.system("rm -f " + running_jobs_file)
        running = open(running_jobs_file, "w")
        for line in keep_lines:
            running.write(line.rstrip() + "\n")
        running.close()
    # Update finished jobs #
    #keep_lines=[]
    #keep_lines = finished_lines[:]
    done_lines = []
    run_lines = []
    removed_jobs=[]
    current_time = float(time.time())
    n_jobs=0
    if os.path.isfile(finished_jobs_file):
        for line in open(finished_jobs_file, "r").readlines():
            n_jobs+=1
    if n_jobs == 0:
        finished_time={}
        print("Empty file of finished_jobs. RECALCULATING ....")
        for folder in os.listdir(result_dir):
            if os.path.isfile(os.path.join(result_dir, folder, folder + "_full.json")):
                json_file = os.path.join(result_dir, folder, folder + "_full.json")
                json_content = readJSON(json_file)
                job_id = folder
                job_type=json_content["type"]
                if "job_name" in  json_content.keys():
                    job_name=json_content["job_name"]
                else:
                    job_name=job_id
                if "job_time" in json_content.keys():
                    job_time = json_content["job_time"]
                    date,hours=job_time.split()
                    day,month,year=date.split("/")
                    hour,minute,second=hours.split(":")
                    date_time_obj = datetime(int(year),int(month),int(day),int(hour),int(minute),int(second))
                    #epoch_time=date_time_obj.strftime('%s')
                    epoch_time=os.path.getmtime(os.path.join(result_dir, folder))
                    print("JOB %s TIME %s EPOCH %s"%(job_name,job_time,epoch_time))
                else:
                    epoch_time=os.path.getmtime(os.path.join(result_dir, folder))
                    time_folder=datetime.fromtimestamp(epoch_time)
                    job_time = time_folder.strftime("%d/%m/%Y  %H:%M:%S")
                    print("JOB %s EPOCH %s TIME %s"%(job_name,epoch_time,job_time))
                line = str(job_id) + "\t" + str(job_time) + "\t" + str(job_name) + "\t" + str(job_type) + "\t" + str(epoch_time)
                if "error_msg" in json_content.keys():
                    # This means that there has been a error in this execution #
                    print("Written error")
                    failed_jobs.append(line.rstrip()+"\t"+json_content["error_msg"]+"\n")
                    continue
                elif len(json_content.keys()) == 0:
                    # This means that there has been a error in this execution #
                    print("empty json")
                    failed_jobs.append(line.rstrip()+"\t"+"Wrong input\n")
                    continue
                else:
                    # This means that the job finished but that couldn't be included into the finished jobs section #
                    print("finished job")
                    finished_time.setdefault(float(epoch_time),[]).append(line)
                    #finished_lines.append(line)
                    continue
        print("REORDERING ....")
        for epoch_time in sorted(finished_time.keys(),reverse=True):
            print("Time "+str(epoch_time))
            for line in finished_time[epoch_time]:
                print("  -- ADD "+line)
                finished_lines.append(line)
    if os.path.isfile(finished_jobs_file):
        shutil.copy(finished_jobs_file,finished_jobs_file2)
        with open(finished_jobs_file2, "r") as finished:
         for line in finished:
            if len(line.split()) == 0: continue
            line_id = line.split("\t")[0]
            job_start_time = float(line.split("\t")[-1])
            if line_id in [ids.split("\t")[0] for ids in keep_lines]: 
                run_lines.append(line_id)
            if os.path.isfile(os.path.join(result_dir, line_id, line_id + "_full.json")):
                done_lines.append(line)
            else:
                removed_jobs.append(line_id)
        finished.close()
        if len(finished_lines)>0: os.system("rm -f " + finished_jobs_file)
    finished = open(finished_jobs_file, "w")
    #for line in keep_lines:
    for line in finished_lines:
        line_id = line.split("\t")[0]
        if line_id in removed_jobs:continue
        finished.write(line.rstrip() + "\n")
        removed_jobs.append(line_id)
    for line in done_lines:
        line_id = line.split("\t")[0]
        if line_id in removed_jobs:continue
        finished.write(line.rstrip() + "\n")
    finished.close()

    # Remove old jobs #
    current_time = time.time()
    for dire in os.listdir(result_dir):
      try:
        time_difference = float(current_time) - float(os.path.getmtime(os.path.join(result_dir, dire)))
        if time_difference > float(7776000):
            # The directory is 3 months old, erase it #
            print("\n\n\nRemove the results directory: "+dire)
            print("rm -rf " + os.path.join(result_dir, dire))
            os.system("rm -rf " + os.path.join(result_dir, dire))
            print("\n\n")
      except Exception as e:
        print("Failed to check %s"%dire)
        print("Error: %s"%(str(e)))
        continue

    # Remove old failed jobs #
    failed_folders={}
    remove_failed=[]
    for line in failed_jobs:
           failed_folders.setdefault(line.split()[0],line)
    current_time = time.time()
    for dire in os.listdir(result_dir):
      try:
        time_difference = float(current_time) - float(os.path.getmtime(os.path.join(result_dir, dire)))
        if time_difference > float(2592000) and dire in failed_folders.keys():
            # The directory is 15 days old, erase it #
            print("\n\n\nRemove the results of a failed directory: "+dire)
            print("rm -rf " + os.path.join(result_dir, dire))
            os.system("rm -rf " + os.path.join(result_dir, dire))
            print("\n\n")
            remove_failed.append(failed_folders[dire])
      except Exception as e:
        print("Failed to check %s"%dire)
        print("Error: %s"%(str(e)))
        continue
    # Update failed jobs #
    failed_jobs_done = []
    if os.path.isfile(failed_jobs_file):
        with open(failed_jobs_file, 'r') as f_file:
            for line in f_file:
                line_id = line.split("\t")[0]
                failed_jobs_done.append(line_id)
            f_file.close()
    if os.path.isfile(failed_jobs_file) and len(remove_failed)<=0:
        with open(failed_jobs_file, 'ab') as f_file:
            for line in failed_jobs:
                line_id = line.split("\t")[0]
                if line_id in failed_jobs_done: continue
                f_file.write(line.rstrip() + "\n")
            f_file.close()
    else:
        with open(failed_jobs_file, 'w') as f_file:
            for line in failed_jobs:
                print("Failed job "+line)
                if line in remove_failed: 
                    print("Hide obsolete job "+line)
                    continue
                f_file.write(line.rstrip() + "\n")
            f_file.close()




def remove_job(job_id, name):

    print("you are in remove job: " + job_id)
    result_dir = config.get("Paths", "result_dir")
    running_jobs_file = os.path.join(result_dir, "running_jobs.txt")
    finished_jobs_file = os.path.join(result_dir, "finished_jobs.txt")
    print("finished_jobs_file: " + finished_jobs_file)
    failed_jobs_file = os.path.join(result_dir, "failed_jobs.txt")
    # Iterate the three files looking for the job that we want to erase #
    for file in [running_jobs_file, finished_jobs_file, failed_jobs_file]:
        f_open = open(file, "r")
        f = f_open.readlines()
        final_lines = []
        for line in f:
            if line.startswith(str(job_id)) and (name in line):
                print("removing line: " + str(line))
                continue
            else:
                final_lines.append(line)
        f_open.close()
        os.system("rm -f " + file)
        new_f = open(file, "w")
        for line in final_lines:
            new_f.write(line.rstrip() + "\n")
        new_f.close()


if __name__ == '__main__':
    app.run(host='localhost', port=8082, debug=True)
