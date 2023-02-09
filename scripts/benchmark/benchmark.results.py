import os, sys, re
import optparse
import shutil

# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.dirname(__file__))

# Imports my functions #
import functions

# Import my modules #
import pwm, spotentials, triads, x3dna

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python benchmark.results.py -c pdb_chain -o output_file -p potentials_file -t triads_file -x x3dna_file")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-c", action="store", type="string", dest="pdb_chain", help="PDB chain", metavar="{string}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file", metavar="{filename}")
    parser.add_option("-p", action="store", type="string", dest="potentials_file", help="Potentials file (from spotentials.py)", metavar="{filename}")
    parser.add_option("-t", action="store", type="string", dest="triads_file", help="Triads file (from triads.py)", metavar="{filename}")
    parser.add_option("-x", action="store", type="string", dest="x3dna_file", help="X3DNA file (from x3dna.py)", metavar="{filename}")

    (options, args) = parser.parse_args()

    if options.pdb_chain is None or options.output_file is None or options.potentials_file is None or options.triads_file is None or options.x3dna_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Get triads object #
    triads_obj = triads.Triads(options.triads_file)

    # For each triad... #
    for i in sorted(range(len(frozenset(triads_obj._triads))), reverse=True):
        # If complementary triad... #
        if i % 2 != 0:
            # Delete triad #
            del triads_obj._triads[i]

    # Get X3DNA object #
    x3dna_obj = x3dna.X3DNA(options.x3dna_file)

    # Load statistical potential #
    potentials = {}
    potentials.setdefault(options.pdb_chain, spotentials.Potentials(options.potentials_file, "s3dc_dd"))

    # Get dinucleotide raw scores and binding site region #
    scores, binding_site = pwm.get_scores_and_binding_site(triads_obj, x3dna_obj, potentials, "s3dc_dd", options.pdb_chain)

    # Get k-mers scaled scores #
    all_kmers_scaled_scores = pwm.get_kmers_scaled_scores(scores, binding_site, "s3dc_dd") 

    if not os.path.exists(options.output_file):
        # Initialize #
        dummy_file = os.path.join(options.dummy_dir, str(os.getpid()) + ".txt")
        # Write output #
        functions.write(dummy_file, "#kmer;score")
        # For each k-mer ... #
        for kmer, start, end in all_kmers_scaled_scores:
            functions.write(dummy_file, "%s;%.3f" % (kmer, all_kmers_scaled_scores[(kmer, start, end)]))
        # Copy dummy file #
        shutil.copy(dummy_file, options.output_file)
        # Remove dummy file #
        os.remove(dummy_file)
