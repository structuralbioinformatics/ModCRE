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

# Import my modules #
import x3dna, contacts

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python helix.py -i input_file [--dummy=dummy_dir -o output_file]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in PDB format)", metavar="{filename}")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="{directory}")

    (options, args) = parser.parse_args()

    if options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_protein_chains_dna_helices(pdb_obj, x3dna_obj, contacts_obj, output_dir="./", dummy_dir="/tmp"):
    """
    This function identifies the DNA helix that is more likely to interact
    with each {ChainOfProtein} within a {PDB} file. For each {ChainOfProtein},
    it writes the resulting helix in an individual file.

    @input:
    pdb_obj {PDB}
    x3dna_obj {X3DNA}
    contacts_obj {Contacts}
    dummy_dir {string}

    """

    # Initialize #
    max_contact_distance = float(config.get("Parameters", "max_contact_distance"))
    min_contacts = int(config.get("Parameters", "min_contacts"))

    # For protein chain... #
    for protein_chain_obj in pdb_obj.proteins:
        # Initialize #
        contacts = {}
        helix_file = os.path.join(output_dir, pdb_obj.id.lower() + "_" + protein_chain_obj.chain + ".txt")
        # For each helix... #
        for helix in x3dna_obj.get_dna_helices():
            # Initialize #
            contacts.setdefault(helix, set())
            # Get helix dinucleotides #
            dinucleotides = x3dna_obj.get_helix_dinucleotides(helix)
            # For each contact... #
            for contact_obj in contacts_obj.get_contacts():
                # If protein chain... #
                if contact_obj._A_chain == protein_chain_obj.chain:
                    # For each dinucleotide... #
                    for dinucleotide in x3dna_obj.get_basepair_dinucleotides(x3dna_obj.get_residue_basepair(contact_obj._B_chain[0], contact_obj._B_residue_obj[0].number)):
                        # If dinucleotide in helix dinucleotides... #
                        if dinucleotide in dinucleotides:
                            contacts[helix].add(contact_obj)
                            break
        # For each helix... #
        for helix in sorted(contacts, key=lambda x: contacts[x], reverse=True):
            # Skip if not enough contacts #
            if len(contacts[helix]) < min_contacts: continue
            # Erase helix file if already exists #
            if os.path.exists(helix_file): os.remove(helix_file)
            # Create output #
            functions.write(helix_file, str(helix))
            break

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Get PDB object #
    pdb_obj = PDB(os.path.abspath(options.input_file))

    # Get X3DNA object #
    x3dna_obj = x3dna.get_x3dna_obj(os.path.abspath(options.input_file), os.path.abspath(options.dummy_dir))

    # Get contacts object #
    contacts_obj = contacts.get_contacts_obj(pdb_obj, x3dna_obj, "pdi", "dinucleotides", os.path.abspath(options.dummy_dir))

    # Get protein chains DNA helices #
    get_protein_chains_dna_helices(pdb_obj, x3dna_obj, contacts_obj, os.path.abspath(options.output_dir), os.path.abspath(options.dummy_dir))
