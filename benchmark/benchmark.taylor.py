import os, sys, re
import ConfigParser
import numpy
import optparse
import pandas
import socket
import subprocess

# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.dirname(__file__))

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Imports my functions #
import functions, nr, spotentials

# Imports jbonet's module #
from SBI.data import aminoacids1to3, aminoacids_polarity_boolean
from SBI.external.blast import blast_parser

# Import my modules #
import contacts, dream5, model_protein, nr, pbm, pwm, spotentials, tomtom, triads, x3dna

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python benchmark.taylor.py --pbm=pbm_dir --pdb=pdb_dir --pwm=pwm_dir [--dummy=dummy_dir -o output_dir --start=start_step --stop=stop_step -v]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="{directory}")
    parser.add_option("--pbm", action="store", type="string", dest="pbm_dir", help="PBM directory (output directory from pbm.py)", metavar="{directory}")
    parser.add_option("--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (output directory from pdb.py)", metavar="{directory}")
    parser.add_option("--pwm", action="store", type="string", dest="pwm_dir", help="PWM directory (from CIS-BP)", metavar="{directory}")
    parser.add_option("--start", default=1, action="store", type="int", dest="start_step", help="Start at a given step (default = 1; first)", metavar="{int}")
    parser.add_option("--stop", default=10, action="store", type="int", dest="stop_step", help="Stop at a given step (default = 10; last)", metavar="{int}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")

    (options, args) = parser.parse_args()

    if options.pbm_dir is None or options.pdb_dir is None or options.pwm_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Get hostname #
    hostname = socket.gethostname()

    # Create output "main" subdirs #
    for subdir in ["files", "nr", "potentials", "threading", "pwms", "meme", "tomtom", "results", "parsed", "figures"]:
        if not os.path.exists(os.path.join(options.output_dir, subdir)):
            os.makedirs(os.path.join(options.output_dir, subdir))

    # Verbose mode #
    if options.verbose: sys.stdout.write("\n")

    # Parse dimers file #
    dimers = {}
    # For each line... #
    for line in functions.parse_file(os.path.join(options.pdb_dir, "dimers.txt")):
        if line.startswith("#"): continue
        line = line.split(";")
        dimers[line[0]] = line[1]
        dimers[line[1]] = line[0]

    ################################
    # 1. Create initial files      #
    ################################
    if options.verbose: sys.stdout.write("Create initial files...\n\n")
    # Skip if starts later #
    if options.start_step <= 1:
        # Skip if sequences file already exists #
        sequences_file = os.path.join(options.output_dir, "files", "sequences.fasta")
        if not os.path.exists(sequences_file):
            # Initialize #
            sequences = set()
            threading = set()
            # For each threading file... #
            for threading_file in os.listdir(os.path.join(options.pbm_dir, "threading")):
                m = re.search("(T\d+_\d+\.\d+)\.(\d+)\.\S{4}_\S+\.txt", threading_file)
                if m:
                    # Get TF object #
                    tf_obj = pbm.TF(file_name=os.path.join(options.pbm_dir, "tfs", m.group(1) + ".txt"))
                    # Add sequence to sequences #
                    sequences.add(("%s.%s" % (m.group(1), m.group(2)), tf_obj._sequences[int(m.group(2))]))
            # For each sequence... #
            for tf_id, sequence in sorted(sequences):
                # Initialize #
                hits = []
                hits_file = os.path.join(options.pbm_dir, "hits", "%s.txt" % tf_id)
                # For each line... #
                for line in functions.parse_file(hits_file):
                    line = line.split("\t")
                    if line[0] == "None": continue
                    hits.append(line[0])
                # Skip if no hits #
                if len(hits) == 0: continue
                # For each hit... #
                for hit in hits:
                    # If TF dimerizes... #
                    if hit in dimers:
                        # Initialize #
                        dimer = None
                        # For each hit... #
                        for next_hit in hits:
                            if next_hit == dimers[hit]:
                                dimer = next_hit
                                break
                        # Continue if dimer is None #
                        if dimer is None: continue
                    # Add sequence to threading #
                    threading.add((tf_id, hit))
                    break
            # For each sequence... #
            for tf_id, sequence in sorted(sequences):
                # For each threaded TF... #
                for i in threading:
                    if tf_id == i[0]:
                        # Write sequences file #
                        functions.write(sequences_file, ">%s\n%s" % (tf_id, sequence))
            # For each sequence... #
            for tf_id, hit in sorted(threading):
                # Write sequences file #
                functions.write(os.path.join(options.output_dir, "files", "threading.txt"), "%s\t%s" % (tf_id, hit))
    # Exit if stops here #
    if options.stop_step == 1:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ################################
    # 2. Create homology clusters  #
    ################################
    if options.verbose: sys.stdout.write("Create homology clusters...\n\n")
    # Skip if starts later #
    if options.start_step <= 2:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        kclust_path = os.path.join(src_path, config.get("Paths", "kclust_path"))
        thresholds = [(70, 3.53)]
        # For each threshold... #
        for identity, score in thresholds:
            # Skip if representatives file already exists #
            representatives_file = os.path.join(cluster_dir, "representatives.txt")
            if not os.path.exists(representatives_file):
                # Initialize #
                representatives = set()
                # Cluster sequences #
                process =  subprocess.check_output([os.path.join(kclust_path, "kClust"), "-i", os.path.join(options.output_dir, "files", "sequences.fasta"), "-d", os.path.join(options.output_dir, "files"), "-s", str(score)], stderr=subprocess.STDOUT)
                # For each line... #
                for line in functions.parse_file(os.path.join(options.output_dir, "files", "representatives.fas")):
                    if line.startswith(">"):
                        representatives.add(line[1:13])
                # For each TF... #
                for tf_id in sorted(representatives):
                    functions.write(representatives_file, tf_id)
    # Exit if stops here #
    if options.stop_step == 2:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ################################
    # 3. Create non-redundant sets #
    ################################
    if options.verbose: sys.stdout.write("Create non-redundant sets...\n\n")
    # Skip if starts later #
    if options.start_step <= 3:
        # Initialize #
        threading = {}
        thresholds = [70]
        # For each line... #
        for line in functions.parse_file(os.path.join(options.output_dir, "files", "threading.txt")):
            tf_id, hit = line.split("\t")
            threading[tf_id] = hit
        # For each threshold... #
        for identity in thresholds:
            # Initialize #
            representatives = []
            nr_dir = os.path.join(options.output_dir, "nr")
            # Make subdirs #
            if not os.path.exists(nr_dir):
                os.makedirs(nr_dir)
            # For each line... #
            for line in functions.parse_file(os.path.join(options.output_dir, "files", "representatives.fas")):
                if line.startswith(">"):
                    representatives.append(line[1:])
            # For each TF... #
            for tf_id in sorted(representatives):
                # Initialize #
                pdb_chains = [threading[tf_id]]
                # If TF dimerizes... #
                if threading[tf_id] in dimers:
                    pdb_chains.append(threading[tf_id])
                # For each PDB chain... #
                for pdb_chain in pdb_chains:
                    # Skip if nr file already exists #
                    nr_file = os.path.join(nr_dir, "%s.%s.txt" % (tf_id, pdb_chain[-1]))
                    if not os.path.exists(nr_file):
                        # If hostname is cluster... #
                        if hostname == config.get("Cluster", "cluster_name"):
                            # Submit to queue #
                            functions.submit_command_to_queue("%s %s --pdb=%s -r %s -f %s -i %s -l %s -o %s --pbm=%s -t %s" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "nr.py"), options.pdb_dir, os.path.abspath(os.getcwd()), tf_id[:12], pdb_chain, os.path.join(options.output_dir, "files", "representatives.txt"), nr_file, options.pbm_dir, config.get("Parameters", "max_redundancy_family"), ), config.get("Cluster", "cluster_queue"), int(config.get("Cluster", "max_jobs_in_queue")))
                        # Else... #
                        else:
                            # Get non-redundant triads #
                            nr_triads = nr.get_nr_triads(options.pdb_dir, pdb_chain=pdb_chain, threshold=float(config.get("Parameters", "max_redundancy_family")), pbm_dir=options.pbm_dir, filter_tf_id=tf_id[:12], list_file=os.path.join(options.output_dir, "files", "representatives.txt"))
                            # For each nr triads object... #
                            for nr_triads_obj in nr_triads:
                                functions.write(nr_file, nr_triads_obj._file)
    # Exit if stops here #
    if options.stop_step == 3:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ################################
    # 4. Derive stat. potentials   #
    ################################
    if options.verbose: sys.stdout.write("Derive statistical potentials...\n\n")
    # Skip if starts later #
    if options.start_step <= 4:
        # Initialize #
        threading = {}
        thresholds = [70]
        # For each threshold... #
        for identity in thresholds:
            # Initialize #
            potentials_dir = os.path.join(options.output_dir, "potentials")
            # Make subdirs #
            if not os.path.exists(potentials_dir):
                os.makedirs(potentials_dir)
            # For each nr file... #
            for nr_file in os.listdir(os.path.join(options.output_dir, "nr")):
                # Initialize #
                m = re.search("(\S+).txt$", nr_file)
                # Skip if potentials file already exists #
                potentials_file = os.path.join(potentials_dir, m.group(1) + ".txt")
                if not os.path.exists(potentials_file):
                    # If hostname is cluster... #
                    if hostname == config.get("Cluster", "cluster_name"):
                        # Submit to queue #
                        functions.submit_command_to_queue("%s %s -i %s -o %s -s -z" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "spotentials.py"), os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), os.path.abspath(potentials_file)), config.get("Cluster", "cluster_queue"), int(config.get("Cluster", "max_jobs_in_queue")))
                    # Else... #
                    else:
                        # Derive statistical potentials #
                        pmf_3d, pmf_3dc, pmf_s3dc, pmf_s3dc_dd, pmf_s3dc_di, pmf_local, pmf_pair, distances = spotentials.get_statistical_potentials(os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), approach=False, smooth=True, zscores=True)
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
                # Skip if potentials file already exists #
                potentials_file = os.path.join(potentials_dir, m.group(1) + ".taylor.txt")
                if not os.path.exists(potentials_file):
                    # If hostname is cluster... #
                    if hostname == config.get("Cluster", "cluster_name"):
                        # Submit to queue #
                        functions.submit_command_to_queue("%s %s -i %s -o %s -a -s -z" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "spotentials.py"), os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), os.path.abspath(potentials_file)), config.get("Cluster", "cluster_queue"), int(config.get("Cluster", "max_jobs_in_queue")))
                    # Else... #
                    else:
                        # Derive statistical potentials #
                        pmf_3d, pmf_3dc, pmf_s3dc, pmf_s3dc_dd, pmf_s3dc_di, pmf_local, pmf_pair, distances = spotentials.get_statistical_potentials(os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), approach=True, smooth=True, zscores=True)
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
    if options.stop_step == 4:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ###############################
    # 5. Get PWM predictions      #
    ###############################
    if options.verbose: sys.stdout.write("Get PWM predictions...\n\n")
    # Skip if starts later #
    if options.start_step <= 5:
        # Initialize #
        threaded_pdb_chains = {}
        thresholds_identities = [70]
        thresholds_scores = numpy.arange(0.9, 1.01, 0.01)
        # For each line... #
        for line in functions.parse_file(os.path.join(options.output_dir, "files", "threading.txt")):
            tf_id, hit = line.split("\t")
            threaded_pdb_chains[tf_id] = hit
        # For each threshold... #
        for identity in thresholds_identities:
            # Initialize #
            representatives = []
            # For each line... #
            for line in functions.parse_file(os.path.join(options.output_dir, "files", "representatives.fas")):
                if line.startswith(">"):
                    representatives.append(line[1:])
            # For each TF... #
            for tf_id in sorted(representatives):
                # Initialize #
                pdb_chains = [threaded_pdb_chains[tf_id]]
                threading_triads_obj = triads.Triads()
                # Get x3dna object #
                x3dna_obj = x3dna.X3DNA(os.path.join(options.pdb_dir, "x3dna", pdb_chains[0][:4] + ".txt"))
                # If TF dimerizes... #
                if threaded_pdb_chains[tf_id] in dimers:
                    pdb_chains.append(threaded_pdb_chains[tf_id])
                # For each PDB chain... #
                for pdb_chain in pdb_chains:
                    # Get threaded protein #
                    protein = {}
                    threading = None
                    # For each line... #
                    for line in functions.parse_file(os.path.join(os.path.join(options.pbm_dir, "threading", "%s.%s.txt" % (tf_id, pdb_chain)))):
                        if line == "": continue
                        if line.startswith("#") or line.startswith("//"): continue
                        elif line.startswith(">"):
                            threading = line[1:]
                        elif threading == "protein":
                            line = line.split(";")
                            protein.setdefault((line[0], int(line[1])), line[2])                  
                    # Get triads object #
                    triads_obj = triads.Triads(os.path.join(options.pdb_dir, "triads", pdb_chain + ".txt"))
                    # For each triad object... #
                    for triad_obj in triads_obj.get_triads():
                        # Get triad #
                        (residue_A, aminoacid_environment), (residue_B, dinucleotide_environment) = triad_obj.get_triad()
                        distance = triad_obj.get_contact_distance()
                        # Get amino acid chain, number #
                        chain, number = residue_A.split("-")
                        # If threaded residue... #
                        if (chain, int(number)) in protein:
                            # Skip gaps or weird amino acids #
                            if protein[(chain, int(number))] not in aminoacids_polarity_boolean: continue
                            # Get amino acid environment #
                            aminoacid, hydrophobicity, degree_of_exposure, secondary_structure = aminoacid_environment.split("-")
                            # Set new amino acid environment #
                            aminoacid = aminoacids1to3[protein[(chain, int(number))]]
                            hydrophobicity = "N"
                            if aminoacids_polarity_boolean[protein[(chain, int(number))]]:
                                hydrophobicity = "P"
                            aminoacid_environment = "%s-%s-%s-%s" % (aminoacid, hydrophobicity, degree_of_exposure, secondary_structure)
                            # Add triad #
                            threading_triads_obj.add_triad(triads.Triad(aminoacid_environment, dinucleotide_environment, distance, residue_A, residue_B))
                # Skip if no triads #
                if len(threading_triads_obj.get_triads()) == 0: continue
                # For each type of potential... #
                for potential in ["pdb.family", "pdb.family.taylor", "pbm.family", "pbm.family.taylor", "pdb.general", "pdb.general.taylor", "pbm.general", "pbm.general.taylor"]:
                    # Initialize #
                    binding_sites = {}
                    all_kmers_scaled_scores = {}
                    # For each PDB chain... #
                    for pdb_chain in pdb_chains:
                        # Get potentials #
                        potentials = {}
                        if "family" in potentials:
                            if "pdb" in potential:
                                if "taylor" in potentials: potentials.setdefault(pdb_chain[-1],spotentials.Potentials(os.path.join(options.pdb_dir, "potentials", "%s.taylor.txt" % pdb_chain), "s3dc_dd"))
                                else: potentials.setdefault(pdb_chain[-1],spotentials.Potentials(os.path.join(options.pdb_dir, "potentials", "%s.txt" % pdb_chain), "s3dc_dd"))
                            else:
                                if "taylor" in potentials: potentials.setdefault(pdb_chain[-1],spotentials.Potentials(options.output_dir, "potentials", "%s.%s.taylor.txt" % (tf_id, pdb_chain[-1]), "s3dc_dd"))
                                else: potentials.setdefault(pdb_chain[-1],spotentials.Potentials(os.path.join(options.output_dir, "potentials", "%s.%s.txt" % (tf_id, pdb_chain[-1])), "s3dc_dd"))
                        else:
                            if "pdb" in potential:
                                if "taylor" in potentials: potentials.setdefault(pdb_chain[-1],spotentials.Potentials(os.path.join(options.pdb_dir, "potentials", "general.taylor.txt"), "s3dc_dd"))
                                else: potentials.setdefault(pdb_chain[-1],spotentials.Potentials(os.path.join(options.pdb_dir, "potentials", "general.txt"), "s3dc_dd"))
                            else:
                                if "taylor" in potentials: potentials.setdefault(pdb_chain[-1],spotentials.Potentials(os.path.join(options.pbm_dir, "potentials", "general.taylor.txt"), "s3dc_dd"))
                                else: potentials.setdefault(pdb_chain[-1],spotentials.Potentials(os.path.join(options.pbm_dir, "potentials", "general.txt"), "s3dc_dd"))
                        # Get dinucleotide raw scores and binding site region #
                        scores, binding_site = pwm.get_scores_and_binding_site(threading_triads_obj, x3dna_obj, potentials, "s3dc_dd", pdb_chain[-1])
                        # Get k-mers scaled scores #
                        all_kmers_scaled_scores.setdefault(pdb_chain, pwm.get_kmers_scaled_scores(scores, binding_site, "s3dc_dd"))
                        # Add binding sites #
                        binding_sites.setdefault(pdb_chain, binding_site)
                    # Get whole binding site #
                    binding_site = set()
                    # For each PDB chain... #
                    for pdb_chain in sorted(binding_sites):
                        # For each site... #
                        for i in binding_sites[pdb_chain].keys():
                            # Add site to binding site #
                            binding_site.add(i)
                    # Sort binding site #
                    binding_site = sorted(binding_site)
                    # Define binding site #
                    binding_site = list(range(binding_site[0], binding_site[-1] + 2))
                    # Inialize #
                    default_sequence = "N" * (binding_site[-1] - binding_site[0] + 1)
                    # For each score threshold... #
                    for score in thresholds_scores:
                        # Initialize #
                        msa_obj = pwm.MSA()
                        pwm_dir = os.path.join(options.output_dir, "pwms", potential, str(score))
                        # Skip if PWM file already exists #
                        pwm_file = os.path.join(pwm_dir, "%s.txt" % tf_id)
                        if os.path.exists(pwm_file): continue
                        # Add binding site length #
                        msa_obj.set_binding_site_length(len(binding_site))
                        # Make subdirs #
                        if not os.path.exists(pwm_dir):
                            os.makedirs(pwm_dir)
                        # For each PDB chain... #
                        for pdb_chain in sorted(all_kmers_scaled_scores):
                            # For each k-mer... #
                            for kmer, start, end in all_kmers_scaled_scores[pdb_chain]:
                                if all_kmers_scaled_scores[pdb_chain][(kmer, start, end)] >= score:
                                    sequence = default_sequence[:start - binding_site[0]] + kmer + default_sequence[end - binding_site[0] + 1:]
                                    msa_obj.add_sequence(sequence, all_kmers_scaled_scores[pdb_chain][(kmer, start, end)])
                        # Write PWM #
                        msa_obj.write(os.path.join(pwm_dir, "%s.txt" % tf_id), "pwm")
    # Exit if stops here #
    if options.stop_step == 5:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ###############################
    # 6. Create Tomtom database   #
    ###############################
    if options.verbose: sys.stdout.write("Create Tomtom database...\n\n")
    # Skip if starts later #
    if options.start_step <= 6:
        # For each line... #
        for line in functions.parse_file(os.path.join(options.output_dir, "files", "representatives.fas")):
            if line.startswith(">"):
                # Skip if PWM file already exists #
                pwm_file = os.path.join(options.output_dir, "meme", "%s.txt" % line[1:])
                if not os.path.exists(pwm_file):
                    # Get TF object #
                    tf_obj = pbm.TF(file_name=os.path.join(options.pbm_dir, "tfs", "%s.txt" % line[1:13]))
                    i = line[14]
                    # Initialize #
                    pwm_list = []
                    dna_string = "ACGT"
                    dummy_file = os.path.join(options.dummy_dir, "pwm.txt")
                    # For each line... #
                    for line in functions.parse_file(os.path.join(options.pwm_dir, tf_obj._motifs[int(i)] + ".txt")):
                        if line.startswith("Pos"): continue
                        line = line.split("\t")
                        line.pop(0)
                        pwm_list.append(line)
                    pwm_list = zip(*pwm_list)
                    # Initialize #
                    functions.write(dummy_file, "%s" % tf_obj.get_id() + "." + str(i))
                    # For each nucleotide... 
                    for j in range(len(dna_string)):
                        functions.write(dummy_file, "%s:\t%s" % (dna_string[j], "\t".join(pwm_list[j])))
                    # Convert PWM to MEME format #
                    dream5.uniprobe_to_meme(dummy_file, pwm_file)
                    # Remove dummy file #
                    os.remove(dummy_file)
        # Skip if hits file already exists #
        database_file = os.path.join(options.output_dir, "database.txt")
        if not os.path.exists(database_file):
            # For each pwm file... #
            for pwm_file in os.listdir(os.path.join(options.output_dir, "meme")):
                # If a MEME file... #
                if pwm_file.endswith(".txt"):
                    # For each line... #
                    for line in functions.parse_file(os.path.join(options.output_dir, "meme", pwm_file)):
                        functions.write(database_file, line)
    # Exit if stops here #
    if options.stop_step == 6:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ###############################
    # 7. Reformat PWM to MEME     #
    ###############################
    if options.verbose: sys.stdout.write("Reformat PWM to MEME...\n\n")
    # Skip if starts later #
    if options.start_step <= 7:
        # For each type of potential... #
        for potential in ["pdb.family", "pdb.family.taylor", "pbm.family", "pbm.family.taylor", "pdb.general", "pdb.general.taylor", "pbm.general", "pbm.general.taylor"]:
            # For each score threshold... #
            for score in os.listdir(os.path.join(options.output_dir, "pwms", potential)):
                # Initialize #
                pwm_dir = os.path.join(options.output_dir, "pwms", potential, str(score))
                meme_dir = os.path.join(options.output_dir, "meme", potential, str(score))
                # Make subdirs #
                if not os.path.exists(meme_dir):
                    os.makedirs(meme_dir)
                # For each PWM file... #
                for pwm_file in os.listdir(pwm_dir):
                    # Skip if MEME file already exists #
                    meme_file = os.path.join(meme_dir, pwm_file)
                    if not os.path.exists(meme_file):
                        # Initialize #
                        pwm_list = []
                        nucleotides = list("ACGT")
                        dummy_file = os.path.join(options.dummy_dir, "pwm.txt")
                        # For each line... #
                        for line in functions.parse_file(os.path.join(pwm_dir, pwm_file)):
                            line = line.split()
                            line.pop(0)
                            pwm_list.append(line)
                        # Write dummy file #
                        functions.write(dummy_file, pwm_file[:-4])
                        # For each nucleotide... 
                        for i in range(len(nucleotides)):
                            functions.write(dummy_file, "%s:\t%s" % (nucleotides[i], "\t".join(pwm_list[i])))
                        # Convert PWM to MEME format #
                        dream5.uniprobe_to_meme(dummy_file, meme_file)
                        # Remove dummy file #
                        os.remove(dummy_file)
    # Exit if stops here #
    if options.stop_step == 7:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ###############################
    # 8. Execute Tomtom           #
    ###############################
    if options.verbose: sys.stdout.write("Execute Tomtom...\n\n")
    # Skip if starts later #
    if options.start_step <= 8:
        # Initialize #
        database_file = os.path.join(options.output_dir, "database.txt")
        # For each type of potential... #
        for potential in ["pdb.family", "pdb.family.taylor", "pbm.family", "pbm.family.taylor", "pdb.general", "pdb.general.taylor", "pbm.general", "pbm.general.taylor"]:
            # For each score threshold... #
            for score in os.listdir(os.path.join(options.output_dir, "pwms", potential)):
                # Initialize #
                meme_dir = os.path.join(options.output_dir, "meme", potential, str(score))
                tomtom_dir = os.path.join(options.output_dir, "tomtom", potential, str(score))
                # Make subdirs #
                if not os.path.exists(tomtom_dir):
                    os.makedirs(tomtom_dir)
                # For each PWM file... #
                for meme_file in os.listdir(meme_dir):
                    # Skip if tomtom file already exists #
                    tomtom_file = os.path.join(tomtom_dir, meme_file)
                    if os.path.exists(tomtom_file): continue
                    # Get Tomtom object #
                    tomtom_obj = tomtom.get_tomtom_obj(database_file, os.path.join(meme_dir, meme_file))
                    # Write object #
                    tomtom_obj.write(tomtom_file)
    # Exit if stops here #
    if options.stop_step == 8:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ###############################
    # 9. Parse Tomtom files       #
    ###############################
    if options.verbose: sys.stdout.write("Parse Tomtom files...\n\n")
    # Skip if starts later #
    if options.start_step <= 9:
        # Initialize #
        data = []
        tomtom = {}
        data_frame_file = os.path.join(options.output_dir, "tomtom", "data.csv")
        # For each type of potential... #
        for potential in ["pdb.family", "pdb.family.taylor", "pbm.family", "pbm.family.taylor", "pdb.general", "pdb.general.taylor", "pbm.general", "pbm.general.taylor"]:
            # For each score threshold... #
            for score in os.listdir(os.path.join(options.output_dir, "tomtom", potential)):
                # For each tomtom file... #
                for tomtom_file in os.listdir(os.path.join(options.output_dir, "tomtom", potential, score)):
                    # Initialize #
                    tf_id = tomtom_file[:12]
                    tf_id_sequence = tomtom_file[:14]
                    tomtom.setdefault(tf_id, {})
                    tomtom[tf_id].setdefault((potential, score), 1)
                    # For each line... #
                    for line in functions.parse_file(os.path.join(options.output_dir, "tomtom", potential, score, tomtom_file)):
                        if line.startswith("#"):    continue
                        line = line.split()
                        if len(line) == 10:
                            if line[1] == tf_id_sequence:
                                if float(line[3]) < tomtom[tf_id][(potential, score)]:
                                    tomtom[tf_id][(potential, score)] = float(line[3])
        # For each tf id... #
        for tf_id in sorted(tomtom):
            # Get TF object #
            tf_obj = pbm.TF(file_name=os.path.join(options.pbm_dir, "tfs", "%s.txt" % tf_id))
            for label in sorted(tomtom[tf_id]):
                potential, score = label
                p_value = tomtom[tf_id][label]
                data.append([tf_id, tf_obj.get_family(), potential, score, p_value])
        # Sort data #
        data.sort(key=lambda x: x[0])
        data.sort(key=lambda x: x[1].upper())
        # Get data frame #
        data_frame = pandas.DataFrame(data)
        data_frame.columns = ["name", "family", "potential", "score", "p-value"]
        data_frame.to_csv(data_frame_file)
    # Exit if stops here #
    if options.stop_step == 9:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ###############################
    # 10. Parse Tomtom files       #
    ###############################
    if options.verbose: sys.stdout.write("Parse Tomtom files...\n\n")
    # Skip if starts later #
    if options.start_step <= 10:
        # Initialize #
        data_frame = pandas.read_csv(os.path.join(options.output_dir, "tomtom", "data.csv"), index_col=0, parse_dates=True)
        # For each family... #
        for family in ["AFT", "AP2", "APSES", "ARID/BRIGHT", "bHLH", "bZIP", "B3", "C2H2 ZF", "C2HC ZF", "CUT", "CxxC", "DM", "E2F", "Ets", "Forkhead", "GATA", "GCM", "Homeodomain", "IRF", "MADS box", "MBD", "Myb/SANT", "NAC/NAM", "Ndt80/PhoG", "Nuclear receptor", "POU", "Paired box", "Prospero", "Rap1", "Rel", "Runt", "SMAD", "Sox", "T-box", "TBP", "THAP finger", "WRKY", "Zinc cluster"]:
            # For each type of potential... #
            for potential in ["pdb.family", "pdb.family.taylor", "pbm.family", "pbm.family.taylor", "pdb.general", "pdb.general.taylor", "pbm.general", "pbm.general.taylor"]:
                # For each score threshold... #
                for score in os.listdir(os.path.join(options.output_dir, "tomtom", potential)):
                    # Get data frame #
                    print(family, potential, score)
                    df = data_frame[(data_frame["family"] == family) & (data_frame["potential"] == potential) & (data_frame["score"] == float(score))]
                    print(df)
                    exit(0)
    # Exit if stops here #
    if options.stop_step == 10:
        sys.stdout.write("Exiting...\n\n")
        exit(0)