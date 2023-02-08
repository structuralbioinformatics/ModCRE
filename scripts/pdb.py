import os, sys, re
import ConfigParser
from collections import Counter
import gzip
import numpy
import optparse
import shutil
import subprocess

# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.dirname(__file__))

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Imports my functions #
import functions

# Imports jbonet's module #
from SBI.structure import PDB
from SBI.structure.chain import Chain

# Import my modules #
import clean, contacts, dimers, dssp, folds, helices, interface, triads, split, nr, spotentials, x3dna, hmmer

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python pdb.py -p pdb_dir -t tfs_file [--dummy=dummy_dir -o output_dir --start=start_step --stop=stop_step -v]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="{directory}")
    parser.add_option("-p", "--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (where \"all\" and \"biounits\" directories are placed)", metavar="{directory}")
    parser.add_option("--start", default=1, action="store", type="int", dest="start_step", help="Start at a given step (default = 1; first)", metavar="{int}")
    parser.add_option("--stop", default=8, action="store", type="int", dest="stop_step", help="Stop at a given step (default = 8; last)", metavar="{int}")
    parser.add_option("-j","--parallel", default=False, action="store_true", dest="parallel", help="Submit JOBS to Queues in parallel. The program stops, you need to re-start with the requested value (default = False)", metavar="{boolean}")
    parser.add_option("-t", action="store", default=None, type="string", dest="tfs_file", help="TFs file or any set of DNA binding proteins (from tfinder2.py)", metavar="{filename}")
    parser.add_option("--cisbp_family", action="store",  type="string", default=None,dest="family_cisbp_file", help="TFs Family file to check with PFAM the missing families (from CIS-BP; i.e. cisbp_1.02.tf_families.sql, if used the families in file of option '-t' will be updated)", metavar="{filename}")
    parser.add_option("-f", action="store",  default=None, type="string", dest="family_file", help="TFs file with the names of known families of transcription factors (from tfinder2.py, i.e. tfs.txt)", metavar="{filename}")
    parser.add_option("--radius",default=0, action="store", type="float", dest="radius", help="Maximum contact distance to calculate interactions to create PWMs (default=0 implies the use of 'max_contact_distance' from configuration", metavar="{string}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)", metavar="{boolean}")


    (options, args) = parser.parse_args()

    if options.pdb_dir is None or (options.tfs_file is None and options.start_step<8) :
        parser.error("missing arguments: type option \"-h\" for help")

    return options

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Create output "main" subdirs #
    for subdir in ["clean", "dssp", "x3dna", "split", "contacts", "helices", "interfaces", "triads", "database", "folds", "nr", "potentials", "pwms"]:
        if not os.path.exists(os.path.join(options.output_dir, subdir)):
            os.makedirs(os.path.join(options.output_dir, subdir))

    if not os.path.exists(options.dummy_dir): os.makedirs(options.dummy_dir)
            
    ##################################
    # 1. Create associated PDB files #
    ##################################
    if options.verbose: sys.stdout.write("\nCreate associated PDB files...\n\n")
    # Skip if starts later #
    if options.start_step <= 1:

      if options.parallel:
            if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
            else: cluster_queue=config.get("Cluster", "cluster_queue")
            functions.submit_command_to_queue("%s %s --pdb=%s -o %s --dummy=%s -t %s --start=1 --stop=1 -v" % (os.path.join(config.get("Paths", "python_path"), "python"),  os.path.join(scripts_path,os.path.basename(__file__)) , os.path.abspath(options.pdb_dir), os.path.abspath(options.output_dir) , os.path.abspath(options.dummy_dir), os.path.abspath(options.tfs_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
      else:

        ##############################
        # 1.1 Parse TFs file         #
        ##############################
        tfs = set()
        tfs_fam ={}
        # For each line... #
        tfs_file=os.path.abspath(options.tfs_file)
        for line in functions.parse_file(os.path.abspath(options.tfs_file)):
            if line.startswith("#"): continue
            # Add PDB to TFinDit #
            tfs.add(line[:4])
            pdb_tfs=line.split(";")[0]
            chain_tfs=line.split(";")[1]
            fam_tfs=line.split(";")[2]
            tfs_fam.setdefault(pdb_tfs+"_"+chain_tfs,fam_tfs)
        # For each PDB... #
        for pdb_name in sorted(tfs):
            # Verbose mode... #
            if options.verbose: sys.stdout.write("\t%s...\n" % pdb_name)

            ##############################
            # 1.2 Clean PDB              #
            ##############################
            if options.verbose: sys.stdout.write("\t\t-- cleaning PDB...\n")
            # Skip if formatted PDB file already exists #
            pdb_file = os.path.abspath(os.path.join(options.output_dir, "clean", pdb_name + ".pdb"))
            if not os.path.exists(pdb_file):
                try:
                    original_pdb_file = os.path.abspath(os.path.join(options.pdb_dir, "biounits", pdb_name + ".pdb1"))
                    if not os.path.exists(original_pdb_file): original_pdb_file = os.path.abspath(os.path.join(options.pdb_dir, "biounits", pdb_name.upper() + ".pdb1"))
                    if not os.path.exists(original_pdb_file): original_pdb_file = os.path.abspath(os.path.join(options.pdb_dir, "biounits", pdb_name.lower() + ".pdb1"))
                    if not os.path.exists(original_pdb_file): original_pdb_file = os.path.abspath(os.path.join(options.pdb_dir, "all", pdb_name + ".pdb"))
                    if not os.path.exists(original_pdb_file): original_pdb_file = os.path.abspath(os.path.join(options.pdb_dir, "all", pdb_name.upper() + ".pdb"))
                    if not os.path.exists(original_pdb_file): original_pdb_file = os.path.abspath(os.path.join(options.pdb_dir, "all", pdb_name.lower() + ".pdb"))
                    if not os.path.exists(original_pdb_file): original_pdb_file = os.path.abspath(os.path.join(options.pdb_dir, "all" , "pdb" + pdb_name + ".ent"))
                    pdb_obj = clean.get_clean_pdb_obj(original_pdb_file)
                    pdb_obj.write(pdb_file, force=True, clean=True)
                except:
                    sys.stderr.write("\t\t\t%s PDB does not exist!\n" % pdb_name)
                    continue
            # Get PDB object #
            else:
                try:
                    pdb_obj = PDB(pdb_file)
                except:
                    sys.stderr.write("\t\t\t%s PDB cannot be read!\n" % pdb_name)
                    continue

            ##############################
            # 1.3 Get DSSP               #
            ##############################
            if options.verbose: sys.stdout.write("\t\t-- getting DSSP info...\n")
            # Skip if DSSP file already exists #
            dssp_file = os.path.abspath(os.path.join(options.output_dir, "dssp", pdb_name + ".txt"))
            if not os.path.exists(dssp_file):
                try:
                    dssp_obj = dssp.get_dssp_obj(pdb_file,dummy_dir=options.dummy_dir)
                    dssp_obj.write(dssp_file)
                except:
                    sys.stderr.write("\t\t\tcould not retrieve DSSP info!\n")
                    continue
            # Get DSSP object #
            else:
                dssp_obj = dssp.DSSP(dssp_file)

            ##############################
            # 1.4 Get 3DNA               #
            ##############################
            if options.verbose: sys.stdout.write("\t\t-- getting 3DNA info...\n")
            # Skip if X3DNA file already exists #
            x3dna_file = os.path.abspath(os.path.join(options.output_dir, "x3dna", pdb_name + ".txt"))
            if not os.path.exists(x3dna_file):
                try:
                    x3dna_obj = x3dna.get_x3dna_obj(pdb_file,dummy_dir=options.dummy_dir)
                    x3dna_obj.write(x3dna_file)
                except:
                    sys.stderr.write("\t\t\tcould not retrieve 3DNA info!\n")
                    continue
            # Get X3DNA object #
            else:
                x3dna_obj = x3dna.X3DNA(x3dna_file)

            ##############################
            # 1.5 Split PDB to chains    #
            ##############################
            if options.verbose: sys.stdout.write("\t\t-- splitting PDB into chains...\n")
            # Split chains #
            try:
              split.split_pdb_to_chains(pdb_obj, dssp_obj, x3dna_obj, output_dir=os.path.abspath(os.path.join(options.output_dir, "split")))
            except:
              sys.stderr.write("\t\t\tcould not split the whole structure\n")

            ##############################
            # 1.6 Get prot-DNA contacts  #
            ##############################
            if options.verbose: sys.stdout.write("\t\t-- getting protein-DNA contacts...\n")
            # Skip if contacts file already exists #
            contacts_file = os.path.abspath(os.path.join(options.output_dir, "contacts", pdb_name + ".txt"))
            if not os.path.exists(contacts_file):
                contacts_obj = contacts.get_contacts_obj(pdb_obj, x3dna_obj)
                contacts_obj.write(contacts_file)
            # Get contacts object #
            else:
                contacts_obj = contacts.Contacts(contacts_file)

            ##############################
            # 1.7 Get prot-DNA helix     #
            ##############################
            if options.verbose: sys.stdout.write("\t\t-- getting protein-DNA helices...\n")
            # For protein chain... #
            for protein_chain_obj in pdb_obj.proteins:
                # Initialize #
                protein_chain_pdb_file = os.path.join(options.output_dir, "split", pdb_name + "_" + protein_chain_obj.chain + ".pdb")
                # Skip if PDB chain file does not exist #
                if not os.path.exists(protein_chain_pdb_file): continue
                # Get PDB chain object #
                try:
                 protein_chain_pdb_obj = PDB(protein_chain_pdb_file)
                except:
                 sys.stderr.write("\t\t\tcould not get protein-DNA helices\n")
                 continue
                # Modify PDB id... #
                protein_chain_pdb_obj._chains[0]._pdb = pdb_name
                # Get protein chains DNA helices #
                helices.get_protein_chains_dna_helices(protein_chain_pdb_obj, x3dna_obj, contacts_obj, output_dir=os.path.abspath(os.path.join(options.output_dir, "helices")))

            ##############################
            # 1.8 Get prot-DNA interface #
            ##############################
            if options.verbose: sys.stdout.write("\t\t-- getting protein-DNA interface...\n")
            # Initialize #
            min_basepairs = int(config.get("Parameters", "min_basepairs"))
            tf_accepted = config.get("Parameters", "tf_accepted").split(",")
            tf_rejected = config.get("Parameters", "tf_rejected").split(",")
            accepted=set()
            rejected=set()
            for tf_pdb_chain in tf_accepted:
               if tf_pdb_chain != "None": accepted.add(( tf_pdb_chain.split("_")[0] , tf_pdb_chain.split("_")[1] ))
            for tf_pdb_chain in tf_rejected:
               if tf_pdb_chain != "None": rejected.add(( tf_pdb_chain.split("_")[0] , tf_pdb_chain.split("_")[1] ))
            # For protein chain... #
            for protein_chain_obj in pdb_obj.proteins:
                # Initialize #
                helix_file = os.path.abspath(os.path.join(options.output_dir, "helices", pdb_name + "_" + protein_chain_obj.chain + ".txt"))
                # Skip if helix file does not exist #
                if not os.path.exists(helix_file): continue
                # Skip if interface file already exists #
                interface_file = os.path.join(options.output_dir, "interfaces", pdb_name + "_" + protein_chain_obj.chain + ".txt")
                if not os.path.exists(interface_file):
                    # Initialize #
                    protein_chain_pdb_obj = PDB(os.path.join(options.output_dir, "split", pdb_name + "_" + protein_chain_obj.chain + ".pdb"))
                    # For each helix... #
                    for helix in functions.parse_file(helix_file):
                        # Initialize #
                        dna_pdb_file = os.path.join(options.output_dir, "split", pdb_name + ".dna." + helix + ".pdb")
                        # Skip if helix file does not exist #
                        if not os.path.exists(dna_pdb_file): continue
                        try:
                         dna_pdb_obj = PDB(os.path.abspath(dna_pdb_file))
                        except:
                         sys.stderr.write("\t\t\tcould not get protein-DNA interface\n")
                         break
                        for dna_chain_obj in dna_pdb_obj.nucleotides:
                            protein_chain_pdb_obj.add_chain(dna_chain_obj)
                        interface_obj = interface.get_interface_obj(protein_chain_pdb_obj, x3dna_obj, contacts_obj)
                        # Skip if not enough basepairs #
                        if interface_obj.get_interface_basepairs() is None: continue
                        if len(interface_obj.get_interface_basepairs()) < min_basepairs:
                            if not (pdb_name, protein_chain_obj.chain) in accepted: continue
                        if (pdb_name, protein_chain_obj.chain) in rejected: continue
                        interface_obj.write(interface_file)
                        break

            ##############################
            # 1.9 Get prot-DNA triads    #
            ##############################
            if options.verbose: sys.stdout.write("\t\t-- getting protein-DNA triads...\n")
            # For protein chain... #
            for protein_chain_obj in pdb_obj.proteins:
                # Initialize #
                helix_file = os.path.abspath(os.path.join(options.output_dir, "helices", pdb_name + "_" + protein_chain_obj.chain + ".txt"))
                interface_file = os.path.join(options.output_dir, "interfaces", pdb_name + "_" + protein_chain_obj.chain + ".txt")
                # Skip if interface file does not exist #
                if not os.path.exists(interface_file): continue
                # Skip if triads file already exists #
                triads_file = os.path.join(options.output_dir, "triads", pdb_name + "_" + protein_chain_obj.chain + ".txt")
                if not os.path.exists(triads_file):
                    # Initialize #
                    protein_chain_pdb_obj = PDB(os.path.join(options.output_dir, "split", pdb_name + "_" + protein_chain_obj.chain + ".pdb"))
                    # Modify PDB id... #
                    protein_chain_pdb_obj._chains[0]._pdb = pdb_name
                    # For each helix... #
                    for helix in functions.parse_file(helix_file):
                        # Initialize #
                        dna_pdb_file = os.path.join(options.output_dir, "split", pdb_name + ".dna." + helix + ".pdb")
                        # Skip if helix file does not exist #
                        if not os.path.exists(dna_pdb_file): continue
                        dna_pdb_obj = PDB(os.path.abspath(dna_pdb_file))
                        for dna_chain_obj in dna_pdb_obj.nucleotides:
                            protein_chain_pdb_obj.add_chain(dna_chain_obj)
                        triads_obj = triads.get_triads_obj(protein_chain_pdb_obj, dssp_obj, x3dna_obj, contacts_obj, complementary=True)
                        triads_obj.write(triads_file)
                        break

        #################################################
        # 1.10 Modify tfs_file with missing families    #
        #################################################
        if options.verbose: sys.stdout.write("\n\tmodify TFs file %s...\n"%(options.tfs_file))
        #Parse all potential families
        family_id={}
        family_name={}
        family_dbd={}
        tfs_fam_new={}
        src_path = config.get("Paths", "src_path")
        pfam = os.path.join(src_path, config.get("Paths", "pfam_path"))
        if options.verbose: sys.stdout.write("\t-- parsing %s ...\n"%(os.path.abspath(options.family_cisbp_file)))
        if options.family_cisbp_file is not None:
          for line in functions.parse_file(os.path.abspath(options.family_cisbp_file)):
            m = re.search("\('(.+)', '(.+)', '(.+)', .+, .+\)", line)
            if m:           
               family_id.setdefault(m.group(1).upper(),set()).add(m.group(2))
               family_dbd.setdefault(m.group(3).upper(),set()).add(m.group(2))
               family_name.setdefault(m.group(2).upper(),set()).add((m.group(1),m.group(2),m.group(3)))
          for pdb_fasta in os.listdir( os.path.abspath(os.path.join(options.output_dir, "split"))):
            # For protein chain... #
            if not pdb_fasta.endswith(".fasta"):continue
            pdb_codes=pdb_fasta.split(".")
            if pdb_codes[1]=="dna":continue
            pdb_split=pdb_fasta.strip(".fasta")
            pdb_code=pdb_split.split("_")[0]
            pdb_chain=pdb_split.split("_")[1]
            if not tfs_fam.has_key(pdb_split):continue
            if options.verbose: sys.stdout.write("\t\t-- check family of %s (%s)...\n"%(pdb_split,tfs_fam[pdb_split]))
            input_fasta=os.path.abspath(os.path.join(options.output_dir, "split",pdb_fasta))
            if tfs_fam.has_key(pdb_split):
             if tfs_fam[pdb_split]=="Unknown":
                tfs_fam_set=set()
                try:
                 families=hmmer.get_pfam(pfam,input_fasta,dummy_dir=options.dummy_dir)
                except ValueError as e:
                 sys.stderr.write("%s\n"%e)
                 exit(0)
                for fam in families:
                 if family_dbd.has_key(fam.upper()): tfs_fam_set.update(family_dbd[fam.upper()])
                 if family_name.has_key(fam.upper()): tfs_fam_set.update(set([x[1] for x in  family_name[fam.upper()]]))
                if len(tfs_fam_set)>0:
                 tfs_fam_new.setdefault((pdb_code,pdb_chain),",".join([x for x in tfs_fam_set]))
                else:
                 for fam in families:
                     if fam.startswith("DUF"):continue
                     tfs_fam_set.add(fam)
                 if len(tfs_fam_set)>0:
                  tfs_fam_new.setdefault((pdb_code,pdb_chain),",".join([x for x in tfs_fam_set]))
                 else:
                  tfs_fam_new.setdefault((pdb_code,pdb_chain),"Unknown")
             else:
                tfs_fam_new.setdefault((pdb_code,pdb_chain), tfs_fam[pdb_split])
            if options.verbose: sys.stdout.write("\t\t\t-- family is %s ...\n"%(tfs_fam_new[(pdb_code,pdb_chain)]))

          tfs_file_previous=tfs_file+".previous"
          shutil.move(tfs_file,tfs_file_previous)
          functions.write(options.tfs_file, "#%s;%s;%s" %("pdb","chain","family"))
          for pair,families in tfs_fam_new.iteritems():
               pdb,chain=pair
               functions.write(options.tfs_file, "%s;%s;%s"%(pdb,chain,families))           
                    


    if options.parallel and options.start_step<2:
        sys.stdout.write("Exiting: wait until all submitted runs have finnished and restart again with start=2\n")
        exit(0)
    # Exit if stops here #
    if options.stop_step == 1:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ##############################fileName#
    # 2. Build blast formatted db #
    ###############################
    if options.start_step < 2:
        if options.verbose: sys.stdout.write("\n")
    if options.verbose: sys.stdout.write("Build blast formatted database...\n\n")
    # Skip if starts later #
    if options.start_step <= 2:
      #Get the TFs to use
      tfs = set()
      # For each line... #
      for line in functions.parse_file(os.path.abspath(options.tfs_file)):
            if line.startswith("#"): continue
            # Add PDB to TFinDit #
            tfs.add(line[:4])
      if options.parallel:
            if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
            else: cluster_queue=config.get("Cluster", "cluster_queue")
            functions.submit_command_to_queue("%s %s --pdb=%s -o %s --dummy=%s -t %s --start=2 --stop=2 -v" % (os.path.join(config.get("Paths", "python_path"), "python"),  os.path.join(scripts_path,os.path.basename(__file__)) , os.path.abspath(options.pdb_dir), os.path.abspath(options.output_dir) , os.path.abspath(options.dummy_dir), os.path.abspath(options.tfs_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
      else:
        # Skip if blast formatted database already exists #
        database_file = os.path.join(options.output_dir, "database", "database.fasta")
        if not os.path.exists(database_file):
            fasta_file_contents = []
            # For each PDB chain fasta file... #
            for fasta_file in sorted(os.listdir(os.path.join(options.output_dir, "split"))):
                # Skip if not PDB chain fasta file #
                if not fasta_file.endswith(".fasta"): continue
                if "dna" in fasta_file: continue
                #Skip if not in TFs
                if fasta_file[:4] not in tfs: continue
                # Skip if does not interact with DNA #
                triads_file = os.path.join(options.output_dir, "triads", fasta_file[:-6] + ".txt")
                if os.path.exists(triads_file):
                    # Add PDB chain fasta file to list #
                    fasta_file_contents.append(functions.parse_file(os.path.join(options.output_dir, "split", fasta_file)))
            # For each FASTA file... #
            for fasta_file_content in fasta_file_contents:
                for line in fasta_file_content:
                    # Create output #
                    functions.write(database_file, line)
        try:
            # Initialize #
            src_path = config.get("Paths", "src_path")
            blast_path = os.path.join(src_path, config.get("Paths", "blast_path"))
            # Get path and exec process #
            process = subprocess.check_output([os.path.join(blast_path, "makeblastdb"), "-in", database_file, "-dbtype", "prot"], stderr=subprocess.STDOUT)
        except:
            raise ValueError("Could not exec makeblastdb")

    if options.parallel and options.start_step<3:
        sys.stdout.write("Exiting: wait until all submitted runs have finnished and restart again with start=3\n")
        exit(0)
    # Exit if stops here #
    if options.stop_step == 2:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ################################
    # 3. Join PDB chains by fold   #
    ################################
    if options.verbose: sys.stdout.write("Cluster PDB chains by TM-score...\n\n")
    # Skip if starts later #
    if options.start_step <= 3:
        # Initialize #
        min_tm_score = float(config.get("Parameters", "min_tm_score"))
        #Get the TFs to use
        tfs = set()
        # For each line... #
        for line in functions.parse_file(os.path.abspath(options.tfs_file)):
            if line.startswith("#"): continue
            # Add PDB to TFinDit #
            tfs.add(line[:4])
        # prepare parallel jobs
        if options.parallel:
            jobs=[]
        # For each PDB chain file... #
        for pdb_file in sorted(os.listdir(os.path.join(options.output_dir, "split"))):
            # Skip if not PDB chain file #
            if not pdb_file.endswith(".pdb"): continue
            if "dna" in pdb_file: continue
            #Skip if not in TFs
            if pdb_file[:4] not in tfs: continue
            # Skip if does not interact with DNA #
            triads_file = os.path.join(options.output_dir, "triads", pdb_file[:-4] + ".txt")
            if not os.path.exists(triads_file): continue
            # Skip if formatted PDB file already exists #
            folds_file = os.path.join(options.output_dir, "folds", pdb_file[:-4] + ".txt")
            if options.verbose: sys.stdout.write("\t-- check file %s\n"%pdb_file[:-4])
            if not os.path.exists(folds_file):
              if options.parallel:
                python=os.path.join(config.get("Paths","python_path"),"python")
                jobs.append('%s %s -i %s -o %s --pdb="%s" --dummy="%s"'%(python,os.path.join(scripts_path,"folds.py"),pdb_file[:-4],os.path.join(os.path.abspath(options.output_dir), "folds"),os.path.abspath(options.output_dir),os.path.abspath(options.dummy_dir)))
              else:
                # At TM-score > 0.7 the probability of 2 proteins to be in the same SCOP fold is > 90% #
                folds.get_fold(pdb_file[:-4], os.path.abspath(options.output_dir), os.path.join(os.path.abspath(options.output_dir), "folds"),options.dummy_dir)
        if options.parallel and len(jobs)>0:
            # Create tmp directory #
            tmp_dir = os.path.join(options.dummy_dir, "sh")
            if os.path.exists(tmp_dir) == False:
                os.mkdir(tmp_dir)
            # Initialize #
            min_jobs_in_queue=int(config.get("Cluster","min_jobs_in_queue"))
            if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
            else: cluster_queue=config.get("Cluster", "cluster_queue")
            command_queue=config.get("Cluster","command_queue")
            submit=config.get("Cluster","cluster_submit")
            n_run=0
            n_jobs=0
            
            for job in jobs:
              if n_jobs == min_jobs_in_queue:
                 n_jobs=0
                 n_run=n_run+1
                 out.close()
              if n_jobs==0: 
                 process="p_%d.sh"%(n_run)
                 dummy_file = os.path.join(tmp_dir, process)
                 shutil.copy(os.path.join(scripts_path,config.get("Cluster","command_queue")),dummy_file)
                 out = open(dummy_file , "a")
              out.write("%s\n"%(job))
              n_jobs=n_jobs+1
            out.close()
            # Get current working directory #
            cwd = os.getcwd()
            # Change directory #
            os.chdir(tmp_dir)
            # Submit all jobs
            for n in xrange(0,n_run):
             sh_file="p_%d.sh"%(n)
             if cluster_queue is None:
               os.system("%s  %s" % (submit,sh_file)) 
             else:
               os.system("%s -q %s %s" % (submit,cluster_queue,sh_file)) 
            # Return to original directory #
            os.chdir(cwd)
            # Exit and wait for end of submissions
            sys.stdout.write("Exiting: wait until submitted jobs are completed. Continue at start=4\n\n")
            exit(0)

    if options.parallel and options.start_step<4:
        sys.stdout.write("Exiting: wait until all submitted runs have finnished and restart again with start=4\n")
        exit(0)

    # Exit if stops here #
    if options.stop_step == 3:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ################################
    # 4. Create non-redundant sets #
    ################################
    if options.verbose: sys.stdout.write("Create non-redundant sets...\n")
    # Skip if starts later #
    if options.start_step <= 4:
        #Get the TFs to use
        tfs = set()
        # For each line... #
        for line in functions.parse_file(os.path.abspath(options.tfs_file)):
            if line.startswith("#"): continue
            # Add PDB to TFinDit #
            tfs.add(line[:4])
        # Initialize #
        triads_objects_files = []
        # For each contacts file... #
        for triads_file in sorted(os.listdir(os.path.join(os.path.abspath(options.output_dir), "triads"))):
            #Skip if not in TFs
            if triads_file[:4] not in tfs: continue
            triads_objects_files.append(os.path.join(os.path.abspath(options.output_dir), "triads", triads_file))
            #triads_objects.append(triads.Triads(os.path.join(os.path.abspath(options.output_dir), "triads", triads_file)))

        ##############################
        # 4.1 General nr set         #
        ##############################
        if options.verbose: sys.stdout.write("\t\t-- nr general set...\n")
        # Skip if nr file already exists #
        nr_file = os.path.join(options.output_dir, "nr", "general.txt")
        if not os.path.exists(nr_file):
          # If parallelized #
          if options.parallel:
                # Submit to queue #
                if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
                else: cluster_queue=config.get("Cluster", "cluster_queue")
                functions.submit_command_to_queue("%s %s --pdb=%s -r %s -o %s " % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "nr.py"), options.output_dir, os.path.abspath(os.getcwd()), nr_file), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
          # Else... #
          else:
            # Write output #
            functions.write(nr_file, "#filename")
            # For each nr triads object... #
            for nr_triads_obj in nr.get_nr_triads(options.output_dir, "general", threshold=float(config.get("Parameters", "max_redundancy_general"))):
                functions.write(nr_file, nr_triads_obj._file)

        ##############################
        # 4.2 Family nr sets         #
        ##############################
        if options.verbose: sys.stdout.write("\t\t-- nr family sets...\n")
        # For each triads object... #
        folds_check_done={}
        for triads_obj_file in triads_objects_files:
            triads_obj=triads.Triads(triads_obj_file)
            m = re.search("(\S{4}_\S).txt$", triads_obj._file)
            nr_file = os.path.join(options.output_dir, "nr", m.group(1) + ".txt")
            # Skip if folds file is empty or missing #
            folds_file = os.path.join(options.output_dir, "folds", m.group(1) + ".txt")
            if folds_check_done.has_key(folds_file):
              skip_nr=folds_check_done.get(folds_file)
            else:
              if os.path.exists(folds_file):
                folds_check=[]
                skip_nr=True
                # For each line... #
                for line in functions.parse_file(folds_file):
                   if line.startswith("#"): continue
                   if len(folds_check)>0: break
                   pdb_chain, tm_score = line.strip().split(";")
                   # Add PDB chain to fold
                   folds_check.append(pdb_chain)
                   if len(folds_check)>0: 
                      skip_nr=False
                folds_check_done.setdefault(folds_file,skip_nr)
              else:
                  sys.stdout.write("\t\t\t-- Missing folds %s\n"% m.group(1))
                  continue
            if skip_nr:
                  sys.stdout.write("\t\t\t-- Empty folds %s\n"% m.group(1))
                  continue
            # Skip if nr file already exists #
            if not os.path.exists(nr_file):
              # If parallelized #
              if options.parallel:
                  # Submit to queue #
                  if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
                  else: cluster_queue=config.get("Cluster", "cluster_queue")
                  functions.submit_command_to_queue("%s %s --pdb=%s -r %s -i %s -o %s  -t %s" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "nr.py"), options.output_dir, os.path.abspath(os.getcwd()), m.group(1), nr_file, config.get("Parameters", "max_redundancy_family")), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
              # Else... #
              else:
                # Write output #
                functions.write(nr_file, "#filename")
                # For each nr triads object... #
                for nr_triads_obj in nr.get_nr_triads(options.output_dir, pdb_chain=m.group(1), threshold=float(config.get("Parameters", "max_redundancy_family"))):
                    functions.write(nr_file, nr_triads_obj._file)
    # Verbose mode #
    if options.verbose: sys.stdout.write("\n")
    # Exit if stops here #
    if options.parallel and options.start_step<5:
        sys.stdout.write("Exiting: wait until all submitted runs have finnished and restart again with start=5\n")
        exit(0)
    if options.stop_step == 4:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ################################
    # 5. Derive stat. potentials   #
    ################################
    if options.verbose: sys.stdout.write("Derive statistical potentials...\n\n")
    # Skip if starts later #
    if options.start_step <= 5:
        # For each nr file... #
        for nr_file in os.listdir(os.path.join(options.output_dir, "nr")):
            # Initialize #
            m = re.search("(\S+).txt$", nr_file)
            # Set labels' dictionary
            label={}
            label.setdefault((True,True,True),'.taylor.bins')
            label.setdefault((False,True,True),'.bins')
            label.setdefault((False,False,True),'.pmf.bins')
            label.setdefault((True,False,True),'.pmf.taylor.bins')
            label.setdefault((True,True,False),'.taylor.acc')
            label.setdefault((False,True,False),'.acc')
            label.setdefault((False,False,False),'.pmf.acc')
            label.setdefault((True,False,False),'.pmf.taylor.acc')
            if options.verbose: sys.stdout.write("\t-- potentials of %s\n"%nr_file)
            for approach,zscores,bins in label:
               # Skip if potentials file already exists #
               potentials_file = os.path.abspath(os.path.join(options.output_dir, "potentials", m.group(1) + label[(approach,zscores,bins)]+".txt"))
               if not os.path.exists(potentials_file):
                # If parallelized #
                if options.parallel:
                    if options.verbose: sys.stdout.write("\t\t-- ")
                    if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
                    else: cluster_queue=config.get("Cluster", "cluster_queue")
                    # Submit to queue #
                    if bins:
                     if approach and zscores:
                      functions.submit_command_to_queue("%s %s -i %s -o %s -s -a -z -b" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "spotentials.py"), os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), os.path.abspath(potentials_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                     if not approach and zscores:
                      functions.submit_command_to_queue("%s %s -i %s -o %s -s -z -b" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "spotentials.py"), os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), os.path.abspath(potentials_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                     if approach and not zscores:
                      functions.submit_command_to_queue("%s %s -i %s -o %s -s -a -b" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "spotentials.py"), os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), os.path.abspath(potentials_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                     if not approach and not zscores:
                      functions.submit_command_to_queue("%s %s -i %s -o %s -s -b" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "spotentials.py"), os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), os.path.abspath(potentials_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                    else: 
                     if approach and zscores:
                      functions.submit_command_to_queue("%s %s -i %s -o %s -s -a -z" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "spotentials.py"), os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), os.path.abspath(potentials_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                     if not approach and zscores:
                      functions.submit_command_to_queue("%s %s -i %s -o %s -s -z" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "spotentials.py"), os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), os.path.abspath(potentials_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                     if approach and not zscores:
                      functions.submit_command_to_queue("%s %s -i %s -o %s -s -a" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "spotentials.py"), os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), os.path.abspath(potentials_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                     if not approach and not zscores:
                      functions.submit_command_to_queue("%s %s -i %s -o %s -s " % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "spotentials.py"), os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), os.path.abspath(potentials_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                # Else... #
                else:
                   # Derive statistical potentials #
                   if options.verbose: sys.stdout.write("\t\t--Derive statistical potential %s\n"%(potentials_file))
                   pmf_3d, pmf_3dc, pmf_s3dc, pmf_s3dc_dd, pmf_s3dc_di, pmf_local, pmf_pair, distances = spotentials.get_statistical_potentials(os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), approach=approach, smooth=True, zscores=zscores, computation=bins, dummy_dir=options.dummy_dir)
                   # Initialize #
                   potentials_obj = spotentials.Potentials()
                   # Add statistical potentials to object #
                   potentials_obj._pmf_3d = pmf_3d
                   potentials_obj._pmf_3dc = pmf_3dc
                   potentials_obj._pmf_s3dc = pmf_s3dc
                   potentials_obj._pmf_s3dc_dd = pmf_s3dc_dd
                   potentials_obj._pmf_s3dc_di = pmf_s3dc_di
                   potentials_obj._pmf_local = pmf_local
                   potentials_obj._pmf_pair = pmf_pair
                   potentials_obj._distances = distances
                   # Write statistical potentials #
                   potentials_obj.write(potentials_file)
    # Exit if stops here #
    if options.verbose: sys.stdout.write("\n")
    if options.parallel and options.start_step<6:
        sys.stdout.write("Exiting: wait until all submitted runs have finnished and restart again with start=6\n")
        exit(0)
    if options.stop_step == 5:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ################################
    # 6. Identify dimers in PDBs   #
    ################################
    if options.verbose: sys.stdout.write("Identify PDB chain dimers...\n\n")
    # Skip if starts later #
    if options.start_step <= 6:
        # Skip if dimers file already exists #
        dimers_file = os.path.join(options.output_dir, "dimers.txt")
        if os.path.exists(dimers_file):
           previous=dimers_file+".previous"
           if os.path.exists(previous): os.remove(previous)
           shutil.move(dimers_file,previous)
        if not os.path.exists(dimers_file):
         if options.parallel:
            if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
            else: cluster_queue=config.get("Cluster", "cluster_queue")
            functions.submit_command_to_queue("%s %s --pdb=%s -o %s --dummy=%s -t %s --start=6 --stop=6 -v" % (os.path.join(config.get("Paths", "python_path"), "python"),  os.path.join(scripts_path,os.path.basename(__file__)) , os.path.abspath(options.pdb_dir), os.path.abspath(options.output_dir) , os.path.abspath(options.dummy_dir), os.path.abspath(options.tfs_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
         else:
            # Initialize #
            pdb_dimers = set()
            pdb_dimerize = set()
            pdb_chains_dimerize = set()
            families_dimerize = config.get("Parameters", "families_dimerize").split(",")
            # For each line... #
            for line in functions.parse_file(os.path.abspath(options.tfs_file)):
                if line.startswith("#"): continue
                line = line.split(";")
                # Add TF chain if dimerizes #
                if line[-1] in families_dimerize:
                    pdb_chains_dimerize.add("%s_%s" % (line[0], line[1]))
                    if os.path.exists(os.path.join(options.output_dir, "triads", "%s_%s.txt" % (line[0], line[1]))):
                        pdb_dimerize.add(line[0])
            # For each PDB file... #
            if len(pdb_dimerize)>0:
              for pdb_name in sorted(pdb_dimerize):
                # Skip rejected PDB entries #
                for dimer in dimers.get_dimers(os.path.join(os.path.abspath(options.output_dir), "clean", "%s.pdb" % pdb_name), os.path.abspath(options.output_dir)):
                    if dimer[0] in pdb_chains_dimerize and dimer[1] in pdb_chains_dimerize:
                        pdb_dimers.add(dimer)
            else: 
              for pdb_file in sorted(os.listdir(os.path.join(options.output_dir, "clean"))):
                # Skip rejected PDB entries #
                for dimer in dimers.get_dimers(os.path.join(os.path.abspath(options.output_dir), "clean", pdb_file), os.path.abspath(options.output_dir)):
                    pdb_dimers.add(dimer)
            # Write output #
            functions.write(dimers_file, "#monomerA;monomerB;contacts;overlap")
            for dimer in sorted(pdb_dimers):
                functions.write(dimers_file, "%s" % ";".join(map(str, dimer)))

    # Exit if stops here #
    if options.parallel and options.start_step<7:
        sys.stdout.write("Exiting: wait until all submitted runs have finnished and restart again with start=7\n")
        exit(0)
    # Exit if stops here #
    if options.stop_step == 6:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ################################
    # 7. Get PDB chain families    #
    ################################
    if options.verbose: sys.stdout.write("Get PDB chain structural family...\n\n")
    # Skip if starts later #
    if options.start_step <= 7:
        # Skip if families file already exists #
        families_file = os.path.join(options.output_dir, "families.txt")
        if os.path.exists(families_file):
           previous=families_file+".previous"
           if os.path.exists(previous): os.remove(previous)
           shutil.move(families_file,previous)
        if not os.path.exists(families_file):
         if options.family_file is None:
            sys.stderr.write("Missing a file with the names of families of known TFs. Please add option '-f' with the missing file (from tfinder2, i.e. tfs.txt)\n")
            exit(0)
         if not os.path.exists(os.path.abspath(options.family_file)):
            sys.stderr.write("Missing a file with the names of families of known TFs. Please check a correct file (from tfinder2, i.e. tfs.txt)\n")
            exit(0)
         if options.parallel:
            if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
            else: cluster_queue=config.get("Cluster", "cluster_queue")
            functions.submit_command_to_queue("%s %s --pdb=%s -o %s --dummy=%s -t %s -f %s --start=7 --stop=7 -v" % (os.path.join(config.get("Paths", "python_path"), "python"),  os.path.join(scripts_path,os.path.basename(__file__)) , os.path.abspath(options.pdb_dir), os.path.abspath(options.output_dir) , os.path.abspath(options.dummy_dir), os.path.abspath(options.tfs_file), os.path.abspath(options.family_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
         else:
            # Initialize #
            tfs = {}
            # Write output #
            functions.write(families_file, "#pdb_chain;family")
            # For each line... #
            for linen in functions.parse_file(os.path.abspath(options.family_file)):
                if linen.startswith("#"): continue
                line = linen.strip().split(";")
                tfs.setdefault(line[0] + "_" + line[1],[]).extend( line[2].split(","))
            # Check if exist triads
            failed={}
            for pdb_chain,families in tfs.iteritems():
                most_common_family = Counter(families).most_common(2)
                if most_common_family[0][0] != "Unknown" :
                  if not  os.path.exists(os.path.join(options.output_dir, "triads",pdb_chain + ".txt")):
                   pdb,chain=pdb_chain.split("_")
                   failed.setdefault(pdb,set()).add((chain,most_common_family[0][0]))
            failed_list=[pdb for pdb in failed.iterkeys()]
            # For each triads file... #
            for triads_file in os.listdir(os.path.join(options.output_dir, "triads")):
                # Initialize #
                pdb_chain = triads_file[:6]
                folds_file = os.path.join(options.output_dir, "folds", pdb_chain + ".txt")
                # Exist fold, then check family
                families = []
                # For each line... #
                for linen in functions.parse_file(folds_file):
                    if linen.startswith("#"): continue
                    line = linen.split(";")
                    if line[0] not in tfs: continue
                    for family in tfs[line[0]]:
                        families.append(family)
                most_common_family = Counter(families).most_common(2)
                if len(families)>0:
                    pdb,chain=pdb_chain.split("_")
                    if pdb in failed_list:
                       failed_list.remove(pdb)
                if len(most_common_family)>1:
                   second_common_ratio = float (most_common_family[1][1])/(float(most_common_family[1][1])+float(most_common_family[0][1]))
                try:
                    if len(most_common_family)>1:
                      if most_common_family[0][0] != "Unknown" :
                         functions.write(families_file, "%s;%s\n# ----\t\t\t# %s" % (pdb_chain, most_common_family[0][0],str(Counter(families))))
                      elif ( second_common_ratio > 0.25):
                         functions.write(families_file, "%s;%s\n# ----\t\t\t# %s" % (pdb_chain, most_common_family[1][0],str(Counter(families))))
                      else:
                         functions.write(families_file, "%s;%s\n# ----\t\t\t# %s" % (pdb_chain, most_common_family[0][0],str(Counter(families))))
                    else:
                      functions.write(families_file, "%s;%s\n# ----\t\t\t# %s" % (pdb_chain, most_common_family[0][0],str(Counter(families))))
                except:
                    # No family #
                    functions.write(families_file, "%s;Undefined\n# ----\t\t\t# %s" % (pdb_chain,str(families)))
            # Write missing pdb with known family
            for pdb in failed_list:
                 for chain_family in failed.get(pdb):
                     chain=chain_family[0]
                     family=chain_family[1]
                     pdb_chain=pdb+"_"+chain
                     functions.write(families_file, "#%s;%s;Failed (not in a biounit)" % (pdb_chain,family))
    if options.parallel and options.start_step<8:
        sys.stdout.write("Exiting: wait until all submitted runs have finnished, check Failed folds in pdb/families.txt, fix or continue with start=8\n")
        exit(0)

    # Exit if stops here #
    if options.stop_step == 7:
        sys.stdout.write("Exiting...\n")
        exit(0)

    ######################
    # 8. Calculate PWMs #
    ######################
    if options.verbose: sys.stdout.write("Calculate position weight matrices...\n\n")
    # Skip if starts later #
    if options.start_step <= 8:

       #Initialize
       if not os.path.exists(os.path.join(os.path.abspath(options.dummy_dir), "pwm")):
         os.makedirs(os.path.join(os.path.abspath(options.dummy_dir), "pwm"))
       radius=float(options.radius)
       if radius <= 0: radius=float(config.get("Parameters", "max_contact_distance"))
       #Read dimers
       if options.verbose:sys.stdout.write("\t\tget dimers...\n")
       dimers={}
       dimers_file=open(os.path.join(os.path.abspath(options.output_dir),"dimers.txt"),"r")
       for line in dimers_file:
         if line.startswith("#"):continue
         chain_A,chain_B,contacts,overlap=line.strip().split(";")
         dimers.setdefault(chain_A,set().add(chain_B))
         dimers.setdefault(chain_B,set().add(chain_A))
       if options.verbose:sys.stdout.write("\t\tget families...\n")
       families={}
       for line in functions.parse_file(os.path.join(options.output_dir, "families.txt")):
         if line.startswith("#"): continue
         pdb_name, family = line.split(";")
         families.setdefault(pdb_name,family)
       #Read chains in pdb/split folder
       for pdb_file in sorted(os.listdir(os.path.join(options.output_dir, "split"))):
         skip=False
         code = pdb_file.strip().split(".")
         if code[-1] != "pdb": continue
         #get the PDB file
         pdb_code  = pdb_file[:4]
         try:
           pdb_obj   = PDB(os.path.join(options.output_dir,"clean",pdb_code+".pdb"))
         except:
           skip=True
           continue
         complex_chains=set()
         #Select the chains that bind the same dna helix
         if code[1] == "dna": 
            dna_helix = code[2]
            for chain_id in pdb_obj.chain_identifiers:
               if not os.path.exists(os.path.join(os.path.abspath(options.output_dir),"helices",pdb_code+"_"+chain_id+".txt")): continue
               helix_file= open(os.path.join(os.path.abspath(options.output_dir),"helices",pdb_code+"_"+chain_id+".txt"),"r")
               for helix in helix_file:
                 if helix.strip()==dna_helix:  complex_chains.add(chain_id)
               helix_file.close()
         else:
            chain_id=code[0].split("_")[-1]
            if os.path.exists(os.path.join(os.path.abspath(options.output_dir),"helices",pdb_code+"_"+chain_id+".txt")): 
              complex_chains.add(chain_id)
              helix_file= open(os.path.join(os.path.abspath(options.output_dir),"helices",pdb_code+"_"+chain_id+".txt"),"r")
              for helix in helix_file:
                dna_helix=helix.strip()
              helix_file.close()
         #Check conditions to continue: exist protein and exist DNA
         if len( complex_chains ) <= 0: continue 
         if not os.path.exists(os.path.join(os.path.abspath(options.output_dir),"split",pdb_code+".dna."+dna_helix+".pdb")): continue
         #Create a dummy PDB file with the protein-DNA complex
         try:
           dummy_pdb_obj = PDB()
           dna_pdb_obj   = PDB(os.path.join(os.path.abspath(options.output_dir),"split",pdb_code+".dna."+dna_helix+".pdb"))
           for dna_chain in dna_pdb_obj.chains:
             dummy_pdb_obj.add_chain(dna_chain)
           pdb_chain=pdb_code+"_"+dna_helix
           for chain_id in complex_chains:
             pdb_chain +=   chain_id 
             chain      =   pdb_obj.get_chain_by_id(chain_id)
             dummy_pdb_obj.add_chain(chain)
           dummy_file = os.path.join(os.path.abspath(options.dummy_dir), "pwm", pdb_chain + ".pdb" )
           if not os.path.exists(dummy_file): dummy_pdb_obj.write(dummy_file)
         except:
           skip=True
           continue
         #Start the PWM calculation
         if options.verbose: sys.stdout.write("%s...\n" % pdb_chain)
         output_file = os.path.join(os.path.abspath(options.output_dir), "pwms", pdb_chain)
         output_pwm  = output_file+".pwm"
         output_meme = output_file+".meme"
         output_msa  = output_file+".msa"
         output_logo = output_file+".logo"
         output_fwd_logo = output_file+".logo.fwd.png"
         output_rev_logo = output_file+".logo.rev.png"
         if os.path.exists(output_pwm) and os.path.exists(output_meme) and os.path.exists(output_msa) and os.path.exists(output_fwd_logo) and os.path.exists(output_rev_logo):
            skip=True
         if not skip:
            if options.parallel:
              if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
              else: cluster_queue=config.get("Cluster", "cluster_queue")
              functions.submit_command_to_queue("%s %s --radius=%f -i %s -o %s --pdb=%s -v --dummy=%s --known --taylor --family " % ( os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "pwm_pbm.py"), radius, dummy_file, output_file, os.path.abspath(options.output_dir) , os.path.abspath(options.dummy_dir)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
            else:
              if os.path.exists(output_msa):
                  msa_obj=PWM.MSA(output_msa,output_file)
              else:
                  # protein is dummy_pdb_obj
                  # Get DSSP object #
                  if options.verbose:sys.stdout.write("\t\t\t-- calculate secondary structure ...\n")
                  try:
                    dummy_dssp_obj = dssp.get_dssp_obj(dummy_file)
                  except:
                    continue
                  # Get X3DNA object #
                  if options.verbose:sys.stdout.write("\t\tget DNA %s ...\n"%pdb_chain)
                  try:
                    dummy_x3dna_obj = x3dna.get_x3dna_obj(dummy_file)
                  except:
                    continue
                  # Get contacts object #
                  if options.verbose:sys.stdout.write("\t\t\t-- calculate contacts ...\n")
                  try:
                    dummy_contacts_obj = contacts.get_contacts_obj(dummy_pdb_obj, dummy_x3dna_obj)
                  except:
                    continue
                  # Get triads object #
                  if options.verbose:sys.stdout.write("\t\t\t-- calculate protein-dna pairs ...\n")
                  try:
                    dummy_triads_obj = triads.get_triads_obj(dummy_pdb_obj, dummy_dssp_obj, dummy_x3dna_obj, dummy_contacts_obj)
                  except:
                    continue
                  # Load statistical potential #
                  if options.verbose:sys.stdout.write("\t\tload potentials ...\n")
                  try:
                    potentials, thresholds,radii,  structural_homologs_by_chain = PWM.load_statistical_potentials(dummy_pdb_obj, os.path.abspath(options.output_dir), None, families, options.radius, None, "s3dc_dd", False, True, False, None, True, False,False,True ,None, os.path.abspath(options.dummy_dir),options.verbose)
                  except:
                    continue
                  # Get MSA object #
                  if options.verbose:sys.stdout.write("\t\tget MSA object ...\n")
                  try:
                    msa_obj = PWM.get_msa_obj(dummy_triads_obj, dummy_x3dna_obj, potentials, radii, None, None, "s3dc_dd" , thresholds)
                  except:
                    continue
                  motif_name=pdb_chain
                  msa_obj.set_motif(motif_name)
                  # Write PWM #
                  if options.verbose:sys.stdout.write("\t\twrite PWM ...\n")
                  if not os.path.exists(output_pwm):
                     if options.verbose:sys.stdout.write("\t\t\t--PWM...\n")
                     msa_obj.write(output_pwm, option="pwm")
                  if not os.path.exists(output_meme):
                     if options.verbose:sys.stdout.write("\t\t\t--PWM in MEME format...\n")
                     msa_obj.write(output_meme, option="meme")
                     #Skip meme format with uniprobe2meme. 
                     #Remove next commented line if you want both ways of obtaining the PWM
                     #PWM.write_pwm_by_meme(msa_obj,output_meme+".s",options.dummy_dir)
                  if not os.path.exists(output_msa):
                     if options.verbose:sys.stdout.write("\t\t\t--MSA...\n")
                     msa_obj.write(output_msa, option="msa")
                  if not os.path.exists(output_fwd_logo) or  not os.path.exists(output_rev_logo):
                     if options.verbose:sys.stdout.write("\t\t\t--Logos...\n")
                     PWM.write_logo(msa_obj,output_logo,options.dummy_dir)
                  os.remove( dummy_file )
         else:
            if options.verbose: sys.stdout.write("...Done...\n")




    if options.parallel:
        sys.stdout.write("Exiting: wait until your submitted run have finnished and it's done\n")
        exit(0)
    else:
        sys.stdout.write("done\n")
        exit(0)

   ################################


