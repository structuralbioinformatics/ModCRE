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
# Classes     #
#-------------#

class Interface(object):
    """
    This class defines an {Interface} object.

    """

    def __init__(self, file_name=None):
        self._file = file_name
        self._interface = {}
        # Initialize #
        if self._file is not None:
            self._parse_file()

    def _parse_file(self):
        for line in functions.parse_file(self._file):
            if line.startswith("#"): continue
            basepair, contacts, hbonds, chbonds, sbridges, vwaals, backbone, nucleobase = line.split(";")
            self._interface.setdefault(int(basepair), {'contacts': int(contacts), 'hbonds': int(hbonds), 'chbonds': int(chbonds), 'sbridges': int(sbridges), 'vwaals': int(vwaals), 'backbone': int(backbone), 'nucleobase': int(nucleobase)})

    def add_basepair_contact(self, basepair, contact_obj=None):
        """
        This function assigns a contact to a basepair in the protein-DNA
        interface.

        @input:
        basepair {int}
        contact_obj {Contact}

        """
        # Initialize #
        self._interface.setdefault(basepair, {'contacts': 0, 'hbonds': 0, 'chbonds': 0, 'sbridges': 0, 'vwaals': 0, 'backbone': 0, 'nucleobase': 0})

        if contact_obj is not None:
            if type(contact_obj._B) is list:
                self._interface[basepair]['contacts'] += 1
            else:
                contact = False
                if contact_obj.is_hydrogen_bond():
                    self._interface[basepair]['hbonds'] += 1
                    contact = True
                elif contact_obj.is_salt_bridge():
                    self._interface[basepair]['sbridges'] += 1
                    contact = True
                elif contact_obj.is_C_hydrogen_bond():
                    self._interface[basepair]['chbonds'] += 1
                    contact = True
                elif contact_obj.is_van_der_waals():
                    self._interface[basepair]['vwaals'] += 1
                    contact = True
                if contact:
                    self._interface[basepair]['contacts'] += 1
                    if contact_obj.contacts_backbone():
                        self._interface[basepair]['backbone'] += 1
                    else:
                        self._interface[basepair]['nucleobase'] += 1

    def has_basepair(self, basepair):
        return basepair in self._interface

    def get_start(self):
        return self.get_interface_start()

    def get_interface_start(self):
        for basepair in sorted(self._interface):
            if self._interface[basepair]['contacts'] > 0:
                return basepair

        return None

    def get_end(self):
        return self.get_interface_end()

    def get_interface_end(self):
        end = None

        for basepair in sorted(self._interface):
            if self._interface[basepair]['contacts'] > 0:
                end = basepair

        return end

    def get_interface_length(self):
        if self.get_interface_basepairs() is not None: return len(self.get_interface_basepairs())
        else: return 0

    def get_basepair_contacts(self, basepair):
        if self.has_basepair(basepair):
            return self._interface[basepair]

        return None

    def get_interface_basepairs(self):
        if self.get_interface_start() is not None and self.get_interface_end() is not None:
            return range(self.get_interface_start(), self.get_interface_end() + 1)

        return None

    def get_interface_overlap(self, interface_obj):
        """
        This function calculates returns the overlap between two
        {Interface objects.

        @input:
        interface_obj {int}
        
        @return: {set}

        """
        
        return set(self.get_interface_basepairs()).intersection(interface_obj.get_interface_basepairs()) 

    def write(self, file_name):
        # Initialize #
        done = set()
        functions.write(file_name, "#basepair;contacts;hbonds;chbonds;sbridges;vwaals;backbone;nucleobase")

        for basepair in sorted(self._interface):
            functions.write(file_name, self.return_as_string(basepair))

    def return_as_string(self, basepair):
        # Initialize #
        contacts = self._interface[basepair]['contacts']
        hbonds = self._interface[basepair]['hbonds']
        chbonds = self._interface[basepair]['chbonds']
        sbridges = self._interface[basepair]['sbridges']
        vwaals = self._interface[basepair]['vwaals']
        backbone = self._interface[basepair]['backbone']
        nucleobase = self._interface[basepair]['nucleobase']
        
        return "%s;%s;%s;%s;%s;%s;%s;%s" % (str(basepair), str(contacts), str(hbonds), str(chbonds), str(sbridges), str(vwaals), str(backbone), str(nucleobase))

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python interface.py -i input_file [-d distance_type --dummy=dummy_dir -o output_file]")

    parser.add_option("-d", default="basepairs", action="store", type="string", dest="distance_type", help="Distance type (i.e. \"basepairs\", \"dinucleotides\" or \"mindist\"; default = dinucleotides)", metavar="{string}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in PDB format)", metavar="{filename}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")

    (options, args) = parser.parse_args()

    if options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")
    if options.distance_type != "basepairs" and options.distance_type != "dinucleotides" and options.distance_type != "mindist":
        parser.error("incorrect distance type: accepted values are \"basepairs\" (for contact distance between amino acid cb and geometric center of a basepair), \"dinucleotides\" (for contact distance between amino acid cb and geometric center of a dinucleotide) and \"mindist\" (for minimum distance)")

    return options

def get_interface_obj(pdb_obj, x3dna_obj, contacts_obj, dummy_dir="/tmp"):
    """
    This function calculates the protein-DNA interface of a interface of a complex
    from a series of provided contacts and returns an {Interface} object.

    @input:
    pdb_obj {PDB}
    x3dna_obj {X3DNA}
    contacts_obj {Contacts}
    distance_type {string} either "basepairs", "dinucleotides" or "mindist"
    dummy_dir {string}

    @return:
    interface_obj {Contacts}

    """

    # Initialize #
    done = set() # Removes redundant contacts from interface
    interface_obj = Interface()

    # For each basepair... #
    for basepair in x3dna_obj.get_basepairs():
        # Initialize #
        interface_obj.add_basepair_contact(basepair)
        # For each PDB chain, residue number... #
        for pdb_chain, residue_num in x3dna_obj.get_basepair(basepair):
            if not pdb_obj.chain_exists(pdb_chain): continue
            if not pdb_obj.get_chain_by_id(pdb_chain).residue_exists(str(residue_num)): continue
            nucleotide_obj = pdb_obj.get_chain_by_id(pdb_chain).get_residue_by_identifier(str(residue_num))
            # For each contact #
            for contact_obj in contacts_obj.get_contacts():
                if type(contact_obj._B) is list:
                    if (basepair, contact_obj._A_residue_obj) in done: continue
                else:
                    if (basepair, contact_obj._A, contact_obj._B) in done: continue
                if not pdb_obj.chain_exists(contact_obj._A_chain): continue
                if not pdb_obj.get_chain_by_id(contact_obj._A_chain).residue_exists(str(contact_obj._A_residue_obj.number)): continue
                # Initialize #
                proceed = True
                if type(contact_obj._B) is list:
                    # For each dinucleotide... #
                    for i in range(len(contact_obj._B)):
                        if not pdb_obj.chain_exists(contact_obj._B_chain[i]):
                            proceed = False
                            break
                        if not pdb_obj.get_chain_by_id(contact_obj._B_chain[i]).residue_exists(str(contact_obj._B_residue_obj[i].number)):
                            proceed = False
                            break
                else:
                    if not pdb_obj.chain_exists(contact_obj._B_chain): continue
                    if not pdb_obj.get_chain_by_id(contact_obj._B_chain).residue_exists(str(contact_obj._B_residue_obj.number)): continue
                # This is a valid contact: proceed... #
                if proceed:
                    # Initialize #
                    same_basepair = False
                    # If residue B is a list of residues... #
                    if type(contact_obj._B) is list:
                        # For each residue object... #
                        for i in range(len(contact_obj._B_residue_obj)):
                            if contact_obj._B_chain[i] == pdb_chain and contact_obj._B_residue_obj[i].number == nucleotide_obj.number:
                                same_basepair = True
                                break
                    # Else... #
                    else:
                        if contact_obj._B_chain == pdb_chain and contact_obj._B_residue_obj.number == nucleotide_obj.number: same_basepair = True
                    # If contacts the same basepair... #
                    if same_basepair:
                        interface_obj.add_basepair_contact(basepair, contact_obj)
                        if type(contact_obj._B) is list:
                            done.add((basepair, contact_obj._A_residue_obj))
                        else:
                            done.add((basepair, contact_obj._A, contact_obj._B))

    return interface_obj

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
    contacts_obj = contacts.get_contacts_obj(pdb_obj, x3dna_obj, "pdi", options.distance_type, os.path.abspath(options.dummy_dir))

    # Get interface object #
    interface_obj = get_interface_obj(pdb_obj, x3dna_obj, contacts_obj, os.path.abspath(options.dummy_dir))

    # Output #
    if options.output_file is not None:
        outx3dna  =options.output_file+".x3dna.out"
        outcontact=options.output_file+".contact.out"
        interface_obj.write(os.path.abspath(options.output_file))
        x3dna_obj.write(os.path.abspath(outx3dna))
        contacts_obj.write(os.path.abspath(outcontact))
    else:
        sys.stdout.write("#basepair;contacts;hbonds;chbonds;sbridges;vwaals;backbone;nucleobase\n")
        for basepair in sorted(interface_obj._interface):
            sys.stdout.write("%s\n" % interface_obj.return_as_string(basepair))
