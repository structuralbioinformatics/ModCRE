import os, sys, re
import ConfigParser
import optparse
import random
import shutil

# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.dirname(__file__))

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Imports my functions #
import functions, pbm

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python benchmark.parsed.py -n n --pbm=pbm_dir [-o output_dir -r randoms]")

    parser.add_option("-n", action="store", type="int", dest="n", help="N (this number is added to label files/folders)", metavar="{int}")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="{directory}")
    parser.add_option("--pbm", action="store", type="string", dest="pbm_dir", help="PBM directory (output directory from pbm.py)", metavar="{directory}")
    parser.add_option("-r", default=100, action="store", type="int", dest="randoms", help="Number of randomly selected negative k-mers (default = 100)", metavar="{int}")    

    (options, args) = parser.parse_args()

    if options.n is None or options.pbm_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Initialize #
    tfs = {}
    max_escore_negatives = float(config.get("Parameters", "max_escore_negatives"))
    # Create negatives subdir #
    if not os.path.exists(os.path.join(options.output_dir, "parsed", "negatives.%s" % str(options.n))):
        os.makedirs(os.path.join(options.output_dir, "parsed", "negatives.%s" % str(options.n)))
    # For each TF file... #
    for tf_file in os.listdir(os.path.join(options.pbm_dir, "tfs")):
        # Initialize #
        tf_id = tf_file[:-4]
        tfs.setdefault(tf_id, pbm.TF(file_name=os.path.join(options.pbm_dir, "tfs", tf_file)))
    # Initialize #
    done = set()
    # For each negatives file... #
    for negatives_file in os.listdir(os.path.join(options.output_dir, "parsed", "negatives.%s" % str(options.n))):
        m = re.search("(T\S+).(\d+).(\S{6}).([ACGT]{8}).(\d+).txt$", negatives_file)
        # For each line... #
        for line in functions.parse_file(os.path.join(options.output_dir, "parsed", "negatives.%s" % str(options.n), negatives_file)):
            done.add((m.group(1), int(m.group(2)), m.group(3), m.group(4), int(m.group(5)), line))

    ##############################
    # 5.1 PDB general potentials #
    ##############################
    pdb_general = set()
    # If output file already exists... #
    pdb_general_file = os.path.join(options.output_dir, "parsed", "pdb.general.%s.txt" % str(options.n))
    if os.path.exists(pdb_general_file):
        # For each line... #
        for line in functions.parse_file(pdb_general_file):
            if line.startswith("#"): continue
            line = line.split(";")
            if int(line[9]):
                pdb_general.add((line[0], int(line[1]), line[2], line[3], int(line[4])))
    else:
        functions.write(pdb_general_file, "#tf_id;experiment;pdb_chain;kmer;start;end;fold;family;score;class")
    # For each line... #
    for line in functions.parse_file(os.path.join(options.output_dir, "folds", "general.txt")):
        if line.startswith("#"): continue
        # Initialize #
        kmers = {}
        negatives = set()
        positive_kmers = []
        negative_kmers = []
        file_name, fold = line.split(";")
        m = re.search("/(T\S+).(\d+).(\S{6}).([ACGT]{8}).(\d+)-(\d+).txt$", file_name)
        tf_id = m.group(1)
        i = m.group(2)
        pdb_chain = m.group(3)
        kmer = m.group(4)
        start = int(m.group(5))
        end = int(m.group(6))
        # Skip if already done #
        if (tf_id, int(i), pdb_chain, kmer, start) in pdb_general: continue
        # If negatives file exists... #
        negatives_file = os.path.join(options.output_dir, "parsed", "negatives.%s" % str(options.n), tf_id + "." + i + "." + pdb_chain + "." + kmer + "." + str(start) + ".txt")
        if os.path.exists(negatives_file):
            # For each line... #
            for line in functions.parse_file(negatives_file):
                negatives.add(line)
        # For each line... #
        for line in functions.parse_file(os.path.join(options.pbm_dir, "motifs", tfs[tf_id]._motifs[int(i)] + ".txt")):
            if line.startswith("#"): continue
            line = line.split(";")
            try:
                kmers.setdefault(line[0], float(line[2]))
                kmers.setdefault(line[1], float(line[2]))
            except:
                pass
        # For each line... #
        for line in functions.parse_file(os.path.join(options.output_dir, "results", "pdb.general", tf_id + "." + i + "." + pdb_chain + "." + str(start) + "-" + str(end) + "." + fold + ".txt")):
            if line.startswith("#"): continue
            line = line.split(";")
            if line[0] == kmer:
                positive_kmers.append((tf_id, int(i), pdb_chain, line[0], start, end, int(fold), tfs[tf_id].get_family(), float(line[1]), 1))
            else:
                if len(negatives) > 0:
                    if line[0] in negatives:
                        negative_kmers.append((tf_id, int(i), pdb_chain, line[0], start, end, int(fold), tfs[tf_id].get_family(), float(line[1]), 0))
                elif (tf_id, int(i), pdb_chain, kmer, start, line[0]) not in done:
                    if line[0] not in kmers: continue
                    elif kmers[line[0]] <= max_escore_negatives:
                        negative_kmers.append((tf_id, int(i), pdb_chain, line[0], start, end, int(fold), tfs[tf_id].get_family(), float(line[1]), 0))
        # If k-mers for TF... #
        if len(positive_kmers) == 1 and len(negative_kmers) >= options.randoms:
            # Write positive k-mer #
            functions.write(pdb_general_file, "%s" % ";".join(map(str, positive_kmers[0])))
            # Randomly select negative k-mers #
            for random_negative_kmer in random.sample(negative_kmers, options.randoms):
                # Write negative k-mer #
                functions.write(pdb_general_file, "%s" % ";".join(map(str, random_negative_kmer)))
                # Add to done #
                done.add((positive_kmers[0][:5] + tuple([random_negative_kmer[3]])))
                if len(negatives) == 0:
                    # Write negatives file #
                    functions.write(negatives_file, random_negative_kmer[3])

    ###################################################
    # 5.2 PDB general potentials approached by Taylor #
    ###################################################
    pdb_general_taylor = set()
    # If output file already exists... #
    pdb_general_taylor_file = os.path.join(options.output_dir, "parsed", "pdb.general.taylor.%s.txt" % str(options.n))
    if os.path.exists(pdb_general_taylor_file):
        # For each line... #
        for line in functions.parse_file(pdb_general_taylor_file):
            if line.startswith("#"): continue
            line = line.split(";")
            if int(line[9]):
                pdb_general_taylor.add((line[0], int(line[1]), line[2], line[3], int(line[4])))
    else:
        functions.write(pdb_general_taylor_file, "#tf_id;experiment;pdb_chain;kmer;start;end;fold;family;score")
    # For each line... #
    for line in functions.parse_file(os.path.join(options.output_dir, "folds", "general.txt")):
        if line.startswith("#"): continue
        # Initialize #
        negatives = set()
        positive_kmers = []
        negative_kmers = []
        file_name, fold = line.split(";")
        m = re.search("/(T\S+).(\d+).(\S{6}).([ACGT]{8}).(\d+)-(\d+).txt$", file_name)
        tf_id = m.group(1)
        i = m.group(2)
        pdb_chain = m.group(3)
        kmer = m.group(4)
        start = int(m.group(5))
        end = int(m.group(6))
        # Skip if already done #
        if (tf_id, int(i), pdb_chain, kmer, start) in pdb_general_taylor: continue
        # If negatives file exists... #
        negatives_file = os.path.join(options.output_dir, "parsed", "negatives.%s" % str(options.n), tf_id + "." + i + "." + pdb_chain + "." + kmer + "." + str(start) + ".txt")
        # If negatives file exists... #
        if os.path.exists(negatives_file):
            # For each line... #
            for line in functions.parse_file(negatives_file):
                negatives.add(line)
        # For each line... #
        for line in functions.parse_file(os.path.join(options.output_dir, "results", "pdb.general.taylor", tf_id + "." + i + "." + pdb_chain + "." + str(start) + "-" + str(end) + "." + fold + ".txt")):
            if line.startswith("#"): continue
            line = line.split(";")
            if line[0] == kmer:
                positive_kmers.append((tf_id, int(i), pdb_chain, line[0], start, end, int(fold), tfs[tf_id].get_family(), float(line[1]), 1))
            elif line[0] in negatives:
                negative_kmers.append((tf_id, int(i), pdb_chain, line[0], start, end, int(fold), tfs[tf_id].get_family(), float(line[1]), 0))
            # If k-mers for TF... #
            if len(positive_kmers) == 1 and len(negative_kmers) >= options.randoms: break
        # If k-mers for TF... #
        if len(positive_kmers) == 1 and len(negative_kmers) >= options.randoms:
            # Write positive k-mer #
            functions.write(pdb_general_taylor_file, "%s" % ";".join(map(str, positive_kmers[0])))
            # For each negative k-mer... #
            for negative_kmer in negative_kmers:
                # Write negative k-mer #
                functions.write(pdb_general_taylor_file, "%s" % ";".join(map(str, negative_kmer)))

    ####################################
    # 5.3 PDB + PBM general potentials #
    ####################################
    pbm_general = set()
    # If output file already exists... #
    pbm_general_file = os.path.join(options.output_dir, "parsed", "pbm.general.%s.txt" % str(options.n))
    if os.path.exists(pbm_general_file):
        # For each line... #
        for line in functions.parse_file(pbm_general_file):
            if line.startswith("#"): continue
            line = line.split(";")
            if int(line[9]):
                pbm_general.add((line[0], int(line[1]), line[2], line[3], int(line[4])))
    else:
        functions.write(pbm_general_file, "#tf_id;experiment;pdb_chain;kmer;start;end;fold;family;score")
    # For each line... #
    for line in functions.parse_file(os.path.join(options.output_dir, "folds", "general.txt")):
        if line.startswith("#"): continue
        # Initialize #
        negatives = set()
        positive_kmers = []
        negative_kmers = []
        file_name, fold = line.split(";")
        m = re.search("/(T\S+).(\d+).(\S{6}).([ACGT]{8}).(\d+)-(\d+).txt$", file_name)
        tf_id = m.group(1)
        i = m.group(2)
        pdb_chain = m.group(3)
        kmer = m.group(4)
        start = int(m.group(5))
        end = int(m.group(6))
        # Skip if already done #
        if (tf_id, int(i), pdb_chain, kmer, start) in pbm_general: continue
        # If negatives file exists... #
        negatives_file = os.path.join(options.output_dir, "parsed", "negatives.%s" % str(options.n), tf_id + "." + i + "." + pdb_chain + "." + kmer + "." + str(start) + ".txt")
        # If negatives file exists... #
        if os.path.exists(negatives_file):
            # For each line... #
            for line in functions.parse_file(negatives_file):
                negatives.add(line)
        # For each line... #
        for line in functions.parse_file(os.path.join(options.output_dir, "results", "pbm.general", tf_id + "." + i + "." + pdb_chain + "." + str(start) + "-" + str(end) + "." + fold + ".txt")):
            if line.startswith("#"): continue
            line = line.split(";")
            if line[0] == kmer:
                positive_kmers.append((tf_id, int(i), pdb_chain, line[0], start, end, int(fold), tfs[tf_id].get_family(), float(line[1]), 1))
            elif line[0] in negatives:
                negative_kmers.append((tf_id, int(i), pdb_chain, line[0], start, end, int(fold), tfs[tf_id].get_family(), float(line[1]), 0))
            # If k-mers for TF... #
            if len(positive_kmers) == 1 and len(negative_kmers) >= options.randoms: break
        # If k-mers for TF... #
        if len(positive_kmers) == 1 and len(negative_kmers) >= options.randoms:
            # Write positive k-mer #
            functions.write(pbm_general_file, "%s" % ";".join(map(str, positive_kmers[0])))
            # For each negative k-mer... #
            for negative_kmer in negative_kmers:
                # Write negative k-mer #
                functions.write(pbm_general_file, "%s" % ";".join(map(str, negative_kmer)))

    #########################################################
    # 5.4 PDB + PBM general potentials approached by Taylor #
    #########################################################
    pbm_general_taylor = set()
    # If output file already exists... #
    pbm_general_taylor_file = os.path.join(options.output_dir, "parsed", "pbm.general.taylor.%s.txt" % str(options.n))
    if os.path.exists(pbm_general_taylor_file):
        # For each line... #
        for line in functions.parse_file(pbm_general_taylor_file):
            if line.startswith("#"): continue
            line = line.split(";")
            if int(line[9]):
                pbm_general_taylor.add((line[0], int(line[1]), line[2], line[3], int(line[4])))
    else:
        functions.write(pbm_general_taylor_file, "#tf_id;experiment;pdb_chain;kmer;start;end;fold;family;score")
    # For each line... #
    for line in functions.parse_file(os.path.join(options.output_dir, "folds", "general.txt")):
        if line.startswith("#"): continue
        # Initialize #
        negatives = set()
        positive_kmers = []
        negative_kmers = []
        file_name, fold = line.split(";")
        m = re.search("/(T\S+).(\d+).(\S{6}).([ACGT]{8}).(\d+)-(\d+).txt$", file_name)
        tf_id = m.group(1)
        i = m.group(2)
        pdb_chain = m.group(3)
        kmer = m.group(4)
        start = int(m.group(5))
        end = int(m.group(6))
        # Skip if already done #
        if (tf_id, int(i), pdb_chain, kmer, start) in pbm_general_taylor: continue
        # If negatives file exists... #
        negatives_file = os.path.join(options.output_dir, "parsed", "negatives.%s" % str(options.n), tf_id + "." + i + "." + pdb_chain + "." + kmer + "." + str(start) + ".txt")
        # If negatives file exists... #
        if os.path.exists(negatives_file):
            # For each line... #
            for line in functions.parse_file(negatives_file):
                negatives.add(line)
        # For each line... #
        for line in functions.parse_file(os.path.join(options.output_dir, "results", "pbm.general.taylor", tf_id + "." + i + "." + pdb_chain + "." + str(start) + "-" + str(end) + "." + fold + ".txt")):
            if line.startswith("#"): continue
            line = line.split(";")
            if line[0] == kmer:
                positive_kmers.append((tf_id, int(i), pdb_chain, line[0], start, end, int(fold), tfs[tf_id].get_family(), float(line[1]), 1))
            elif line[0] in negatives:
                negative_kmers.append((tf_id, int(i), pdb_chain, line[0], start, end, int(fold), tfs[tf_id].get_family(), float(line[1]), 0))
            # If k-mers for TF... #
            if len(positive_kmers) == 1 and len(negative_kmers) >= options.randoms: break
        # If k-mers for TF... #
        if len(positive_kmers) == 1 and len(negative_kmers) >= options.randoms:
            # Write positive k-mer #
            functions.write(pbm_general_taylor_file, "%s" % ";".join(map(str, positive_kmers[0])))
            # For each negative k-mer... #
            for negative_kmer in negative_kmers:
                # Write negative k-mer #
                functions.write(pbm_general_taylor_file, "%s" % ";".join(map(str, negative_kmer)))

    #############################
    # 5.5 PDB family potentials #
    #############################
    pdb_family = set()
    # If output file already exists... #
    pdb_family_file = os.path.join(options.output_dir, "parsed", "pdb.family.%s.txt" % str(options.n))
    if os.path.exists(pdb_family_file):
        # For each line... #
        for line in functions.parse_file(pdb_family_file):
            if line.startswith("#"): continue
            line = line.split(";")
            if int(line[9]):
                pdb_family.add((line[0], int(line[1]), line[2], line[3], int(line[4])))
    else:
        functions.write(pdb_family_file, "#tf_id;experiment;pdb_chain;kmer;start;end;fold;family;score")
    # For each line... #
    for line in functions.parse_file(os.path.join(options.output_dir, "folds", "family.txt")):
        if line.startswith("#"): continue
        # Initialize #
        kmers = {}
        negatives = set()
        positive_kmers = []
        negative_kmers = []
        file_name, fold = line.split(";")
        m = re.search("/(T\S+).(\d+).(\S{6}).([ACGT]{8}).(\d+)-(\d+).txt$", file_name)
        tf_id = m.group(1)
        i = m.group(2)
        pdb_chain = m.group(3)
        kmer = m.group(4)
        start = int(m.group(5))
        end = int(m.group(6))
        # Skip if already done #
        if (tf_id, int(i), pdb_chain, kmer, start) in pdb_family: continue
        # If negatives file exists... #
        negatives_file = os.path.join(options.output_dir, "parsed", "negatives.%s" % str(options.n), tf_id + "." + i + "." + pdb_chain + "." + kmer + "." + str(start) + ".txt")
        # If negatives file exists... #
        if os.path.exists(negatives_file):
            # For each line... #
            for line in functions.parse_file(negatives_file):
                negatives.add(line)
        # For each line... #
        for line in functions.parse_file(os.path.join(options.pbm_dir, "motifs", tfs[tf_id]._motifs[int(i)] + ".txt")):
            if line.startswith("#"): continue
            line = line.split(";")
            try:
                kmers.setdefault(line[0], float(line[2]))
                kmers.setdefault(line[1], float(line[2]))
            except:
                pass
        # For each line... #
        for line in functions.parse_file(os.path.join(options.output_dir, "results", "pdb.family", tf_id + "." + i + "." + pdb_chain + "." + str(start) + "-" + str(end) + "." + fold + ".txt")):
            if line.startswith("#"): continue
            line = line.split(";")
            if line[0] == kmer:
                positive_kmers.append((tf_id, int(i), pdb_chain, line[0], start, end, int(fold), tfs[tf_id].get_family(), float(line[1]), 1))
            else:
                if len(negatives) > 0:
                    if line[0] in negatives:
                        negative_kmers.append((tf_id, int(i), pdb_chain, line[0], start, end, int(fold), tfs[tf_id].get_family(), float(line[1]), 0))
                elif (tf_id, int(i), pdb_chain, kmer, start, line[0]) not in done:
                    if line[0] not in kmers: continue
                    elif kmers[line[0]] <= max_escore_negatives:
                        negative_kmers.append((tf_id, int(i), pdb_chain, line[0], start, end, int(fold), tfs[tf_id].get_family(), float(line[1]), 0))
        # If k-mers for TF... #
        if len(positive_kmers) == 1 and len(negative_kmers) >= options.randoms:
            # Write positive k-mer #
            functions.write(pdb_family_file, "%s" % ";".join(map(str, positive_kmers[0])))
            # Randomly select negative k-mers #
            for random_negative_kmer in random.sample(negative_kmers, options.randoms):
                # Write negative k-mer #
                functions.write(pdb_family_file, "%s" % ";".join(map(str, random_negative_kmer)))
                # Add to done #
                done.add((positive_kmers[0][:5] + tuple([random_negative_kmer[3]])))
                if len(negatives) == 0:
                    # Write negatives file #
                    functions.write(negatives_file, random_negative_kmer[3])

    ##################################################
    # 5.6 PDB family potentials approached by Taylor #
    ##################################################
    pdb_family_taylor = set()
    # If output file already exists... #
    pdb_family_taylor_file = os.path.join(options.output_dir, "parsed", "pdb.family.taylor.%s.txt" % str(options.n))
    if os.path.exists(pdb_family_taylor_file):
        # For each line... #
        for line in functions.parse_file(pdb_family_taylor_file):
            if line.startswith("#"): continue
            line = line.split(";")
            if int(line[9]):
                pdb_family_taylor.add((line[0], int(line[1]), line[2], line[3], int(line[4])))
    else:
        functions.write(pdb_family_taylor_file, "#tf_id;experiment;pdb_chain;kmer;start;end;fold;family;score")
    # For each line... #
    for line in functions.parse_file(os.path.join(options.output_dir, "folds", "family.txt")):
        if line.startswith("#"): continue
        # Initialize #
        negatives = set()
        positive_kmers = []
        negative_kmers = []
        file_name, fold = line.split(";")
        m = re.search("/(T\S+).(\d+).(\S{6}).([ACGT]{8}).(\d+)-(\d+).txt$", file_name)
        tf_id = m.group(1)
        i = m.group(2)
        pdb_chain = m.group(3)
        kmer = m.group(4)
        start = int(m.group(5))
        end = int(m.group(6))
        # Skip if already done #
        if (tf_id, int(i), pdb_chain, kmer, start) in pdb_family_taylor: continue
        # If negatives file exists... #
        negatives_file = os.path.join(options.output_dir, "parsed", "negatives.%s" % str(options.n), tf_id + "." + i + "." + pdb_chain + "." + kmer + "." + str(start) + ".txt")
        # If negatives file exists... #
        if os.path.exists(negatives_file):
            # For each line... #
            for line in functions.parse_file(negatives_file):
                negatives.add(line)
        # For each line... #
        for line in functions.parse_file(os.path.join(options.output_dir, "results", "pdb.family.taylor", tf_id + "." + i + "." + pdb_chain + "." + str(start) + "-" + str(end) + "." + fold + ".txt")):
            if line.startswith("#"): continue
            line = line.split(";")
            if line[0] == kmer:
                positive_kmers.append((tf_id, int(i), pdb_chain, line[0], start, end, int(fold), tfs[tf_id].get_family(), float(line[1]), 1))
            elif line[0] in negatives:
                negative_kmers.append((tf_id, int(i), pdb_chain, line[0], start, end, int(fold), tfs[tf_id].get_family(), float(line[1]), 0))
            # If k-mers for TF... #
            if len(positive_kmers) == 1 and len(negative_kmers) >= options.randoms: break
        # If k-mers for TF... #
        if len(positive_kmers) == 1 and len(negative_kmers) >= options.randoms:
            # Write positive k-mer #
            functions.write(pdb_family_taylor_file, "%s" % ";".join(map(str, positive_kmers[0])))
            # For each negative k-mer... #
            for negative_kmer in negative_kmers:
                # Write negative k-mer #
                functions.write(pdb_family_taylor_file, "%s" % ";".join(map(str, negative_kmer)))

    ###################################
    # 5.7 PDB + PBM family potentials #
    ###################################
    pbm_family = set()
    # If output file already exists... #
    pbm_family_file = os.path.join(options.output_dir, "parsed", "pbm.family.%s.txt" % str(options.n))
    if os.path.exists(pbm_family_file):
        # For each line... #
        for line in functions.parse_file(pbm_family_file):
            if line.startswith("#"): continue
            line = line.split(";")
            if int(line[9]):
                pbm_family.add((line[0], int(line[1]), line[2], line[3], int(line[4])))
    else:
        functions.write(pbm_family_file, "#tf_id;experiment;pdb_chain;kmer;start;end;fold;family;score")
    # For each line... #
    for line in functions.parse_file(os.path.join(options.output_dir, "folds", "family.txt")):
        if line.startswith("#"): continue
        # Initialize #
        negatives = set()
        positive_kmers = []
        negative_kmers = []
        file_name, fold = line.split(";")
        m = re.search("/(T\S+).(\d+).(\S{6}).([ACGT]{8}).(\d+)-(\d+).txt$", file_name)
        tf_id = m.group(1)
        i = m.group(2)
        pdb_chain = m.group(3)
        kmer = m.group(4)
        start = int(m.group(5))
        end = int(m.group(6))
        # Skip if already done #
        if (tf_id, int(i), pdb_chain, kmer, start) in pbm_family: continue
        # If negatives file exists... #
        negatives_file = os.path.join(options.output_dir, "parsed", "negatives.%s" % str(options.n), tf_id + "." + i + "." + pdb_chain + "." + kmer + "." + str(start) + ".txt")
        # If negatives file exists... #
        if os.path.exists(negatives_file):
            # For each line... #
            for line in functions.parse_file(negatives_file):
                negatives.add(line)
        # For each line... #
        for line in functions.parse_file(os.path.join(options.output_dir, "results", "pbm.family", tf_id + "." + i + "." + pdb_chain + "." + str(start) + "-" + str(end) + "." + fold + ".txt")):
            if line.startswith("#"): continue
            line = line.split(";")
            if line[0] == kmer:
                positive_kmers.append((tf_id, int(i), pdb_chain, line[0], start, end, int(fold), tfs[tf_id].get_family(), float(line[1]), 1))
            elif line[0] in negatives:
                negative_kmers.append((tf_id, int(i), pdb_chain, line[0], start, end, int(fold), tfs[tf_id].get_family(), float(line[1]), 0))
            # If k-mers for TF... #
            if len(positive_kmers) == 1 and len(negative_kmers) >= options.randoms: break
        # If k-mers for TF... #
        if len(positive_kmers) == 1 and len(negative_kmers) >= options.randoms:
            # Write positive k-mer #
            functions.write(pbm_family_file, "%s" % ";".join(map(str, positive_kmers[0])))
            # For each negative k-mer... #
            for negative_kmer in negative_kmers:
                # Write negative k-mer #
                functions.write(pbm_family_file, "%s" % ";".join(map(str, negative_kmer)))

    ########################################################
    # 5.8 PDB + PBM family potentials approached by Taylor #
    ########################################################
    pbm_family_taylor = set()
    # If output file already exists... #
    pbm_family_taylor_file = os.path.join(options.output_dir, "parsed", "pbm.family.taylor.%s.txt" % str(options.n))
    if os.path.exists(pbm_family_taylor_file):
        # For each line... #
        for line in functions.parse_file(pbm_family_taylor_file):
            if line.startswith("#"): continue
            line = line.split(";")
            if int(line[9]):
                pbm_family_taylor.add((line[0], int(line[1]), line[2], line[3], int(line[4])))
    else:
        functions.write(pbm_family_taylor_file, "#tf_id;experiment;pdb_chain;kmer;start;end;fold;family;score")
    # For each line... #
    for line in functions.parse_file(os.path.join(options.output_dir, "folds", "family.txt")):
        if line.startswith("#"): continue
        # Initialize #
        negatives = set()
        positive_kmers = []
        negative_kmers = []
        file_name, fold = line.split(";")
        m = re.search("/(T\S+).(\d+).(\S{6}).([ACGT]{8}).(\d+)-(\d+).txt$", file_name)
        tf_id = m.group(1)
        i = m.group(2)
        pdb_chain = m.group(3)
        kmer = m.group(4)
        start = int(m.group(5))
        end = int(m.group(6))
        # Skip if already done #
        if (tf_id, int(i), pdb_chain, kmer, start) in pbm_family_taylor: continue
        # If negatives file exists... #
        negatives_file = os.path.join(options.output_dir, "parsed", "negatives.%s" % str(options.n), tf_id + "." + i + "." + pdb_chain + "." + kmer + "." + str(start) + ".txt")
        # If negatives file exists... #
        if os.path.exists(negatives_file):
            # For each line... #
            for line in functions.parse_file(negatives_file):
                negatives.add(line)
        # For each line... #
        for line in functions.parse_file(os.path.join(options.output_dir, "results", "pbm.family.taylor", tf_id + "." + i + "." + pdb_chain + "." + str(start) + "-" + str(end) + "." + fold + ".txt")):
            if line.startswith("#"): continue
            line = line.split(";")
            if line[0] == kmer:
                positive_kmers.append((tf_id, int(i), pdb_chain, line[0], start, end, int(fold), tfs[tf_id].get_family(), float(line[1]), 1))
            elif line[0] in negatives:
                negative_kmers.append((tf_id, int(i), pdb_chain, line[0], start, end, int(fold), tfs[tf_id].get_family(), float(line[1]), 0))
            # If k-mers for TF... #
            if len(positive_kmers) == 1 and len(negative_kmers) >= options.randoms: break
        # If k-mers for TF... #
        if len(positive_kmers) == 1 and len(negative_kmers) >= options.randoms:
            # Write positive k-mer #
            functions.write(pbm_family_taylor_file, "%s" % ";".join(map(str, positive_kmers[0])))
            # For each negative k-mer... #
            for negative_kmer in negative_kmers:
                # Write negative k-mer #
                functions.write(pbm_family_taylor_file, "%s" % ";".join(map(str, negative_kmer)))

