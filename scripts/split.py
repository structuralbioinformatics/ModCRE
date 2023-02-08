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

# Import my functions #
import functions

# Import jbonet's module #
from SBI.structure import PDB
from SBI.structure.chain import ChainOfProtein, ChainOfNucleotide

# Import my modules #
import dssp, x3dna

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python split.py -i input_file [--dummy=dummy_dir -o output_dir]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in PDB format)", metavar="{filename}")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="{directory}")

    (options, args) = parser.parse_args()

    if options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def split_pdb_to_chains(pdb_obj, dssp_obj, x3dna_obj, output_dir="./", dummy_dir="/tmp"):
    """
    This function splits a {PDB} into chains. Note that DNA chains will be
    grouped according to continuous helix and not by individual chains.

    @input:
    pdb_obj {PDB}
    dssp_obj {DSSP}
    x3dna_obj {X3DNA}
    dummy_dir {string}

    """

    # Initialize #
    min_aminoacids = int(config.get("Parameters", "min_aminoacids"))
    min_basepairs = int(config.get("Parameters", "min_basepairs"))

    # For protein chain... #
    for protein_chain_obj in pdb_obj.proteins:
        # Initialize #
        pdb_file = os.path.join(output_dir, pdb_obj.id.lower() + "_" + protein_chain_obj.chain + ".pdb")
        fasta_file = os.path.join(output_dir, pdb_obj.id.lower() + "_" + protein_chain_obj.chain + ".fasta")
        # Create PDB object #
        split_chain_pdb_obj = PDB()
        # Create Chain object #
        split_chain_obj = ChainOfProtein(pdb_obj.id, protein_chain_obj.chain)
        # For each aminoacid... #
        for aminoacid_obj in protein_chain_obj.aminoacids:
            if dssp_obj.has_residue(protein_chain_obj.chain, aminoacid_obj.number):
                split_chain_obj.add_residue(aminoacid_obj)
        # Add split chain to split PDB object #
        split_chain_pdb_obj.add_chain(split_chain_obj)
        # Filter chain if not enough amino acids #
        m = re.search("[^x]{%s}" % min_aminoacids, split_chain_pdb_obj.chains[0].gapped_protein_sequence)
        if not m: continue
        # Create output #
        split_chain_pdb_obj.write(pdb_file, force=True)
        # Erase FASTA file if already exists #
        if os.path.exists(fasta_file): os.remove(fasta_file)
        # Get PDB chain object #
        split_chain_pdb_obj = PDB(pdb_file)
        # Create output #
        functions.write(fasta_file, ">%s_%s\n%s" % (pdb_obj.id.lower(), protein_chain_obj.chain, split_chain_pdb_obj.chains[0].gapped_protein_sequence))
    # For each DNA helix... #
    for helix in x3dna_obj.get_dna_helices():
        # Skip if not enough base pairs... #
        if len(x3dna_obj.get_helix_basepairs(helix)) < min_basepairs: continue
        # Initialize #
        pdb_file = os.path.join(output_dir, pdb_obj.id.lower() + ".dna." + helix + ".pdb")
        fasta_file = os.path.join(output_dir, pdb_obj.id.lower() + ".dna." + helix + ".fasta")
        # Create PDB object #
        split_helix_pdb_obj = PDB()
        # For dna chain... #
        for dna_chain_obj in pdb_obj.nucleotides:
            # Create Chain object #
            split_helix_obj = ChainOfNucleotide(pdb_obj.id, dna_chain_obj.chain)
            # For each nucleotide... #
            for nucleotide_obj in dna_chain_obj.nucleotides:
                if x3dna_obj.helix_has_residue(helix, dna_chain_obj.chain, nucleotide_obj.number):
                    split_helix_obj.add_residue(nucleotide_obj)
            # Add split chain to split PDB object #
            split_helix_pdb_obj.add_chain(split_helix_obj)
        # Create output #
        split_helix_pdb_obj.write(pdb_file, force=True)
        # Erase FASTA file if already exists #
        if os.path.exists(fasta_file): os.remove(fasta_file)
        # Initialize #
        dna_sequence = ""
        for basepair in x3dna_obj.get_helix_basepairs(helix):
            for pdb_chain, residue_num in x3dna_obj.get_basepair(basepair):
                if pdb_obj.chain_exists(pdb_chain):
                    if pdb_obj.get_chain_by_id(pdb_chain).residue_exists(str(residue_num)):
                        dna_sequence += pdb_obj.get_chain_by_id(pdb_chain).get_residue_by_identifier(str(residue_num)).single_letter
                        # Take the sequence of the fwd strand #
                        break
                # Basepair cannot be defined #
                dna_sequence += "n"
                # Take the sequence of the fwd strand #
                break
        # Create output #
        functions.write(fasta_file, ">%s.dna.%s\n%s" % (pdb_obj.id.lower(), helix, dna_sequence))

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Get PDB object #
    pdb_obj = PDB(os.path.abspath(options.input_file))

    # Get DSSP object #
    dssp_obj = dssp.get_dssp_obj(os.path.abspath(options.input_file), os.path.abspath(options.dummy_dir))

    # Get X3DNA object #
    x3dna_obj = x3dna.get_x3dna_obj(os.path.abspath(options.input_file), os.path.abspath(options.dummy_dir))

    # PDB split to chains #
    pdb_split_to_chains(pdb_obj, dssp_obj, x3dna_obj, os.path.abspath(options.output_dir), os.path.abspath(options.dummy_dir))
