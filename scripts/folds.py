import os, sys, re
import ConfigParser
import copy
import optparse

# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.dirname(__file__))

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Imports my functions #
import functions

# Import my modules #
import tmalign

#-------------#
# Options     #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python folds.py -i input_pdb_name_chain -p pdb_dir [--dummy=dummy_dir -o output_dir]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_pdb_name_chain", help="Input PDB name and chain (e.g. chain \"F\" of PDB \"1A02\" would be \"1a02_F\")", metavar="{string}")
    parser.add_option("-o", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="{directory}")
    parser.add_option("-p", "--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (i.e. from pdb.py)", metavar="{directory}")

    (options, args) = parser.parse_args()

    if options.input_pdb_name_chain is None or options.pdb_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_fold(pdb_chain, pdb_dir, output_dir="./", dummy_dir="/tmp"):
    """
    """

    # Initialize #
    fold = []
    folds_file = os.path.join(output_dir, pdb_chain + ".txt")
    min_tm_score = float(config.get("Parameters", "min_tm_score"))

    # For each PDB chain file... #
    for pdb_file in sorted(os.listdir(os.path.join(pdb_dir, "split"))):
        # Skip if not PDB chain file #
        if not pdb_file.endswith(".pdb"): continue
        if "dna" in pdb_file: continue
        # Skip if does not interact with DNA #
        triads_file = os.path.join(pdb_dir, "triads", pdb_file[:-4] + ".txt")
        if not os.path.exists(triads_file): continue
        # Superimpose PDB file over PDB chain #
        tmalign_obj = tmalign.get_tmalign_obj(os.path.join(pdb_dir, "split", pdb_chain + ".pdb"), os.path.join(pdb_dir, "split", pdb_file),dummy_dir)
        if min_tm_score <= min(tmalign_obj.get_tm_scores()):
            fold.append((pdb_file[:-4], min(tmalign_obj.get_tm_scores())))
    # Write output #
    functions.write(folds_file, "#pdb_chain;tm-score")
    for pdb_chain, tm_score in fold:
        functions.write(folds_file, "%s;%s" % (pdb_chain, tm_score))

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Get TMalign object #
    get_fold(options.input_pdb_name_chain, os.path.abspath(options.pdb_dir), os.path.abspath(options.output_dir), os.path.abspath(options.dummy_dir))
