import os, sys, re
import ConfigParser
import optparse

# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.dirname(__file__))

# Append scripts path to python path #
sys.path.append(scripts_path)

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Imports jbonet's module #
from SBI.structure import PDB
from SBI.structure.chain import ChainOfNucleotide

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python clean.py -i input_file [--dummy=dummy_dir -o output_file]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in PDB format)", metavar="{filename}")
    parser.add_option("-o", default="clean.pdb", action="store", type="string", dest="output_file", help="Output file (default = clean.pdb)", metavar="{filename}")

    (options, args) = parser.parse_args()

    if options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_clean_pdb_obj(pdb_file, dummy_dir="/tmp"):
    """
    This function executes cleans a PDB file from HETATMs, overlapped DNA
    chains, etc., and returns a {PDB} object.

    @input:
    pdb_file {string}
    dummy_dir {string}

    @return:
    pdb_obj {PDB}

    """

    try:
        # Get PDB object #
        pdb_obj = PDB(pdb_file)
        # Remove discontinued nucleotides and overlapped DNA chains; it interferes with 3DNA #
        pdb_obj = remove_discontinued_nucleotides_and_overlapped_chains(pdb_obj)
    except:
        raise ValueError("Could not clean PDB %s" % pdb_file)

    return pdb_obj

def remove_discontinued_nucleotides_and_overlapped_chains(pdb_obj):
    """
    This function removes discontinued {Nucleotide}s from {ChainOfNucleotide}
    and any overlapped {ChainOfNucleotide}.

    @input:
    pdb_obj {PDB}

    @return:
    clean_pdb_obj {PDB}

    """

    # Initialize #
    blocks = {}
    chains = {}
    nucleotides = []
    max_overlap_distance = float(config.get("Parameters", "max_overlap_distance"))
    clean_pdb_obj = PDB()

    # Add protein chains #
    for protein_chain in pdb_obj.proteins:
        clean_pdb_obj.add_chain(protein_chain)
    
    # For each DNA chain... #
    for dna_chain in pdb_obj.nucleotides:
        # Initialize #
        blocks.setdefault((dna_chain.pdb, dna_chain.chain), 0)
        chains.setdefault((dna_chain.pdb, dna_chain.chain), [])
        chains[(dna_chain.pdb, dna_chain.chain)].append([])
        chains[(dna_chain.pdb, dna_chain.chain)][-1].append(dna_chain.nucleotides[0])
        # For each residue... #
        for i in range(1, len(dna_chain.nucleotides)):
            if not chains[(dna_chain.pdb, dna_chain.chain)][-1][-1].is_followed(dna_chain.nucleotides[i]):
                # Chain is discontinued! #
                chains[(dna_chain.pdb, dna_chain.chain)].append([])
            # Add nucleotide to chain #
            chains[(dna_chain.pdb, dna_chain.chain)][-1].append(dna_chain.nucleotides[i])
            # Chain largest block length #
            if len(chains[(dna_chain.pdb, dna_chain.chain)][-1]) > blocks[(dna_chain.pdb, dna_chain.chain)]:
                blocks[(dna_chain.pdb, dna_chain.chain)] = len(chains[(dna_chain.pdb, dna_chain.chain)][-1])
    # For each DNA chain... #
    for pdb, pdb_chain in sorted(blocks, key=lambda x: (blocks[x] * -1, x)):
        # For continuous set of nucleotides... #
        for continuous_nucleotides in sorted(chains[(pdb, pdb_chain)], key=lambda x: len(x), reverse=True):
            # Initialize #
            is_overlapped = False
            # For continuous nucleotide... #
            for continuous_nucleotide in continuous_nucleotides:
                # Skip if is overlapped #
                if is_overlapped: continue
                # For nucleotide in DNA chain... #
                for nucleotide in nucleotides:
                    # If residue overlaps... #
                    a, b, distance = continuous_nucleotide.distance(nucleotide, "geometric")
                    if distance < max_overlap_distance:
                        is_overlapped = True
                        break
            # If non-overlapped continuous set of nucleotides... #
            if not is_overlapped:
                chain_obj = ChainOfNucleotide(pdb, pdb_chain)
                for continuous_nucleotide in continuous_nucleotides:
                    chain_obj.add_residue(continuous_nucleotide)
                    nucleotides.append(continuous_nucleotide)
                clean_pdb_obj.add_chain(chain_obj)
                break

    return clean_pdb_obj

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Get clean PDB object #
    pdb_obj = get_clean_pdb_obj(os.path.abspath(options.input_file), os.path.abspath(options.dummy_dir))

    # Create output #
    pdb_obj.write(os.path.abspath(options.output_file), force=True, clean=True)
