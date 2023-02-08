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

# Imports jbonet's module #
from SBI.structure import PDB

# Import my modules #
import interface, contacts

#-------------#
# Options     #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python dimers.py -i pdb_file -p pdb_dir [--dummy=dummy_dir -o output_file]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in PDB format)", metavar="{filename}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")
    parser.add_option("-p", "--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (i.e. from pdb.py)", metavar="{directory}")

    (options, args) = parser.parse_args()

    if options.input_file is None or options.pdb_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_dimers(pdb_file, pdb_dir, dummy_dir="/tmp"):
    """
    """

    # Initialize #
    dimers = set()
    monomers = {}
    min_dimer_contacts = int(config.get("Parameters", "min_dimer_contacts"))
    #sys.stderr.write("\t\t--Dimers: reading file %s\n"%pdb_file)
    try:
     pdb_obj = PDB(pdb_file)
    except:
     sys.stderr.write("\t\t--Skip %s cleaning was unsuccessfull, please remove or modify\n"%pdb_file)
     return dimers

    # For protein chain... #
    for i in range(len(pdb_obj.chains) - 1):
        # Skip if not {ChainOfProtein} #
        if not pdb_obj.chains[i].chaintype == "P": continue
        # Skip if does not interact with DNA #
        triads_file = os.path.join(pdb_dir, "triads", pdb_obj.id.lower() + "_" + pdb_obj.chains[i].chain + ".txt")
        if not os.path.exists(triads_file): continue
        # Get other PDB chains in same fold #
        fold = set()
        for line in functions.parse_file(os.path.join(pdb_dir, "folds", pdb_obj.id.lower() + "_" + pdb_obj.chains[i].chain + ".txt")):
            if line.startswith("#"): continue
            pdb_chain, tm_score = line.split(";")
            if pdb_obj.id.lower() == pdb_chain[:4]:
                if pdb_chain[5] == pdb_obj.chains[i].chain: continue
                fold.add(pdb_chain[5])
        # Skip if fold set is empty #
        if len(fold) == 0: continue
        # Get DNA helix #
        for helix in functions.parse_file(os.path.join(pdb_dir, "helices", pdb_obj.id.lower() + "_" + pdb_obj.chains[i].chain + ".txt")):
            break 
        # Get interface object #
        interface_obj = interface.Interface(os.path.join(pdb_dir, "interfaces", pdb_obj.id.lower() + "_" + pdb_obj.chains[i].chain + ".txt"))
        # For protein chain... #
        for j in range(i + 1, len(pdb_obj.chains)):
            # Skip if not {ChainOfProtein} #
            if not pdb_obj.chains[j].chaintype == "P": continue
            # Skip if does not interact with DNA #
            triads_file = os.path.join(pdb_dir, "triads", pdb_obj.id.lower() + "_" + pdb_obj.chains[j].chain + ".txt")
            if not os.path.exists(triads_file): continue
            # Skip if proteins belong to different folds #
            if pdb_obj.chains[j].chain not in fold: continue
            # Skip if binding sites in different helices #
            for other_helix in functions.parse_file(os.path.join(pdb_dir, "helices", pdb_obj.id.lower() + "_" + pdb_obj.chains[j].chain + ".txt")):
                break
            if other_helix != helix: continue
            # Skip if binding sites do not overlap #
            overlap = interface_obj.get_interface_overlap(interface.Interface(os.path.join(pdb_dir, "interfaces", pdb_obj.id.lower() + "_" + pdb_obj.chains[j].chain + ".txt")))
            if len(overlap) == 0: continue
            # Get residue-residue contacts #
            residue_residue_contacts = set()
            dimer_pdb_obj = PDB()
            dimer_pdb_obj.add_chain(pdb_obj.chains[i])
            dimer_pdb_obj.add_chain(pdb_obj.chains[j])
            contacts_obj = contacts.get_contacts_obj(dimer_pdb_obj, contacts_type="ppi")
            for contact_obj in contacts_obj.get_contacts():
                if contact_obj.is_disulfide_bridge() or contact_obj.is_hydrogen_bond() or contact_obj.is_salt_bridge() or contact_obj.is_van_der_waals():
                    residue_residue_contacts.add((contact_obj._A_residue_obj, contact_obj._B_residue_obj))
            # If enough contacts... #
            if len(residue_residue_contacts) >= min_dimer_contacts:
                dimers.add((pdb_obj.id.lower() + "_" + pdb_obj.chains[i].chain, pdb_obj.id.lower() + "_" + pdb_obj.chains[j].chain, len(residue_residue_contacts), len(overlap)))

    # For each dimer sorted by contacts, overlap... #
    for pdb_chain_A, pdb_chain_B, residue_residue_contacts, overlap in sorted(frozenset(dimers), key=lambda x: x[-2] * x[-1], reverse=True):
        if pdb_chain_A in monomers or pdb_chain_B in monomers:
            sys.stdout.write("\t\t--Dimers: consider as monomers %s %s although they interact by %d AA-pairs\n"%(pdb_chain_A, pdb_chain_B, residue_residue_contacts))
            dimers.remove((pdb_chain_A, pdb_chain_B, residue_residue_contacts, overlap))
        #   continue
        monomers.setdefault(pdb_chain_A, 1)
        monomers.setdefault(pdb_chain_B, 1)

    return dimers

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Get TMalign object #
    dimers = get_dimers(options.input_file, options.pdb_dir, dummy_dir="/tmp")

     # Output #
    if options.output_file is not None:
        functions.write(options.output_file, "#monomerA;monomerB;contacts;overlap")
        for dimer in dimers:
            functions.write(options.output_file, "%s" % ";".join(map(str, dimer)))
    else:
        sys.stdout.write("#monomerA;monomerB;contacts;overlap\n")
        for dimer in dimers:
            sys.stdout.write("%s\n" % ";".join(map(str, dimer)))
