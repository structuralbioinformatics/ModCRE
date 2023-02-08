import os, sys, re
import ConfigParser
import copy
import numpy
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
from SBI.data import element_dic
from SBI.structure import PDB
from SBI.structure.atom import AtomOfAminoAcid, AtomOfNucleotide
from SBI.structure.residue import ResidueOfAminoAcid, ResidueOfNucleotide

# Import my modules #
import x3dna

#-------------#
# Classes     #
#-------------#

class Contacts(object):
    """
    This class defines a {Contacts} object.

    """

    def __init__(self, file_name=None):
        self._file = file_name
        self._contacts = []
        # Initialize #
        if self._file is not None:
            self._parse_file()

    def _parse_file(self):
        for line in functions.parse_file(self._file):
            if line.startswith("#"): continue
            line = line.split(";")
            atomA = re.findall("<(\S+):\s*\[(\S+),\s*(\S+),\s*(.)\]:\((\S+),\s*(\S+),\s*(\S+)\)>", line[0])
            A = AtomOfAminoAcid(atomA[0][2], atomA[0][1], atomA[0][4], atomA[0][5], atomA[0][6], element=atomA[0][3])
            atomB = re.findall("<(\S+):\s*\[(\S+),\s*(\S+),\s*(.)\]:\((\S+),\s*(\S+),\s*(\S+)\)>", line[1])
            if len(atomB) == 1:
                if atomB[0][0] == "AtomOfAminoAcid":
                    B = AtomOfAminoAcid(atomB[0][2], atomB[0][1], atomB[0][4], atomB[0][5], atomB[0][6], element=atomB[0][3])
                if atomB[0][0] == "AtomOfNucleotide":
                    B = AtomOfNucleotide(atomB[0][2], atomB[0][1], atomB[0][4], atomB[0][5], atomB[0][6], element=atomB[0][3])
            elif len(atomB) > 1:
                B = []
                for i in range(len(atomB)):
                    if atomB[i][0] == "AtomOfAminoAcid":
                        B.append(AtomOfAminoAcid(atomB[i][2], atomB[i][1], atomB[i][4], atomB[i][5], atomB[i][6], element=atomB[i][3]))
                    if atomB[i][0] == "AtomOfNucleotide":
                        B.append(AtomOfNucleotide(atomB[i][2], atomB[i][1], atomB[i][4], atomB[i][5], atomB[i][6], element=atomB[i][3]))
            distance = float(line[2])
            residueA = re.findall("<(\S+):\s*\[(\S+),\s*(\S+),\s*\]>", line[3])
            A_residue_obj = ResidueOfAminoAcid(number=residueA[0][2], Rtype=residueA[0][1], version="", mode="")
            residueB = re.findall("<(\S+):\s*\[(\S+),\s*(\S+),\s*\]>", line[4])
            if len(residueB) == 1:
                if residueB[0][0] == "ResidueOfAminoAcid":
                    B_residue_obj = ResidueOfAminoAcid(number=residueB[0][2], Rtype=residueB[0][1], version="", mode="")
                if residueB[0][0] == "ResidueOfNucleotide":
                    B_residue_obj = ResidueOfNucleotide(number=residueB[0][2], Rtype=residueB[0][1], version="", mode="")
            elif len(residueB) > 1:
                B_residue_obj = []
                for i in range(len(residueB)):
                    if residueB[i][0] == "ResidueOfAminoAcid":
                        B_residue_obj.append(ResidueOfAminoAcid(number=residueB[i][2], Rtype=residueB[i][1], version="", mode=""))
                    if residueB[i][0] == "ResidueOfNucleotide":
                        B_residue_obj.append(ResidueOfNucleotide(number=residueB[i][2], Rtype=residueB[i][1], version="", mode=""))
            A_chain = line[5]
            B_chain = line[6].split(",")
            if len(B_chain) == 1:
                B_chain = B_chain[0]
            self.add_contact(Contact(A, B, distance, A_residue_obj, B_residue_obj, A_chain, B_chain))

    def add_contact(self, contact_obj):
        add = True

        for other_contact_obj in self.get_contacts():
            if other_contact_obj._A == contact_obj._A and other_contact_obj._B == contact_obj._B and other_contact_obj._distance == contact_obj._distance:
                add = False
                break

        if add:
            self._contacts.append(contact_obj)

    def get_contacts(self):
        return copy.copy(self._contacts)

    def write(self, file_name, distance_type=None, filter_contacts=False):
        # Initialize #
        if os.path.exists(file_name): os.remove(file_name)
        functions.write(file_name, "#atomA;atomB(s);distance;residueA;residueB(s);chainA;chainB(s)")

        for contact_obj in self.get_contacts():
            if distance_type == "mindist" and filter_contacts:
                if not contact_obj.is_disulfide_bridge() and not contact_obj.is_hydrogen_bond() and not contact_obj.is_salt_bridge() and not contact_obj.is_van_der_waals():
                    continue
            functions.write(file_name, contact_obj.return_as_string())

class Contact(object):
    """
    This class defines a {Contact} object.

    """

    def __init__(self, A, B, distance, A_residue_obj=None, B_residue_obj=None, A_chain=None, B_chain=None):
        self._A = A
        self._B = B
        self._distance = distance
        self._A_residue_obj = A_residue_obj
        self._B_residue_obj = B_residue_obj
        self._A_chain = A_chain
        self._B_chain = B_chain

    def get_contact(self):
        return copy.copy(self)

    def get_contact_distance(self):
        return copy.copy(self._distance)

    def is_disulfide_bridge(self):
        """
        This function returns whether two atoms form a disulfide bridge or not.
        
        According to definition by Mosca R., Ceol A. & Aloy P., 2013, a disulfide
        bridge is formed by any atom pair S-S from 2 Cys at 2.56A.

        @return: {boolean}

        """
        if type(self._B) is list: return None
        
        try:
            if type(self._A) == type(self._B):
                if self._A_residue_obj.single_letter == "C" and self._B_residue_obj.single_letter == "C":
                    if self._A.element == "S" and self._B.element == "S":
                        return self._distance <= 2.56
            return False
        except:
            return None

    def is_hydrogen_bond(self):
        """
        This function returns whether two atoms form a hydrogen bond or not.
        
        According to definition by Mosca R., Ceol A. & Aloy P., 2013, an hydrogen
        bond is formed by a N-O pair at 3.5A.

        @return: {boolean}

        """
        if type(self._B) is list: return None

        try:
            if not self.atoms_clash():
                if (self._A.element == "O" and self._B.element == "N") or (self._A.element == "N" and self._B.element == "O"): 
                    return self._distance <= 3.5
            return False
        except:
            return None

    def is_C_hydrogen_bond(self):
        """
        This function returns whether two atoms form a C-hydrogen bond or not.
        
        According to definition by Mandel-Gutfreund Y., Margalit H., Jernigan J.L.
        & Zhurkin V.B., 1998, a C-hydrogen bond is formed by a CH-O pair at 3.5A.
        Specifically between the C5 of Cytosine and the C5M of Thymine and an O.

        @return: {boolean}

        """
        # Initialize #
        c_hbonds = {'C': "C5", 'T': "C5M"}

        if type(self._B) is list: return None

        try:
            if not self.atoms_clash():
                if type(self._A) != type(self._B):
                    if self._A.element == "O" and self._B.element == "C":
                        if self._B_residue_obj.single_letter in c_hbonds:
                            if self._B.name == c_hbonds[self._B_residue_obj.single_letter]:
                                return self._distance <= 3.5
            return False
        except:
            return None

    def is_salt_bridge(self):
        """
        This function returns whether two atoms form a salt bridge or not.
        
        According to definition by Mosca R., Ceol A. & Aloy P., 2013, a salt
        bridge is formed by any atom pair N-O and O-N at 5.5A.

        @return: {boolean}

        """
        if type(self._B) is list: return None

        try:
            if not self.atoms_clash():
                if (self._A.element == "O" and self._B.element == "N") or (self._A.element == "N" and self._B.element == "O"):
                    return self._distance <= 5.5
            return False
        except:
            return None

    def is_van_der_waals(self):
        """
        This function returns whether two atoms form a van der waals interaction
        or not.
        
        According to definition by Mosca R., Ceol A. & Aloy P., 2013, a Van der
        Waals interaction is formed by any atom pair C-C at 5.0A.
        
        @return: {boolean}

        """
        if type(self._B) is list: return None

        try:
            if not self.atoms_clash():
                if self._A.element == "C" and self._B.element == "C":
                    return self._distance <= 5.0
            return False
        except:
            return None

    def contacts_backbone(self):
        """
        This function returns wheter protein-DNA contact is through the backbone
        or not.

        @return: {boolean}

        """
        if type(self._B) is list: return None

        try:
            if not self.atoms_clash():
                return self._B.is_backbone
            return False
        except:
            return None

    def contacts_nucleobase(self):
        """
        This function returns wheter protein-DNA contact is through the nucleobase
        or not.

        @return: {boolean}

        """
        if type(self._B) is list: return None

        try:
            return bool(abs(self.contacts_backbone() - 1))
        except:
            return None

    def atoms_clash(self):
        """
        This function returns wheter two atoms clash or not.

        According to definition by Mosca R., Ceol A. & Aloy P., 2013, any atom
        pairs at distance less than the sum of the two covalent radii plus 0.5A
        that are not forming a disulfide bridge are considered clashes.

        """
        if type(self._B) is list: return None

        try:
            if not self.is_disulfide_bridge():            
                if self._distance < ((float(element_dic[self._A.element].radius) / 100) + (float(element_dic[self._B.element].radius) / 100) + 0.5):
                    return True
            return False
        except:
            return None

    def return_as_string(self):
        try:
            if type(self._B) is list:
                return "%s;%s;%s;%s;%s;%s;%s" % ('<{0.__class__.__name__}: [{0.name}, {0.number}, {0.element}]:({0.x:.3f}, {0.y:.3f}, {0.z:.3f})>'.format(self._A), ",".join(['<{0.__class__.__name__}: [{0.name}, {0.number}, {0.element}]:({0.x:.3f}, {0.y:.3f}, {0.z:.3f})>'.format(i) for i in self._B]), self._distance, '<{0.__class__.__name__}: [{0.type}, {0.number}, {0.version}]>'.format(self._A_residue_obj), ",".join(['<{0.__class__.__name__}: [{0.type}, {0.number}, {0.version}]>'.format(i) for i in self._B_residue_obj]), self._A_chain, ",".join(self._B_chain))
            else:
                return "%s;%s;%s;%s;%s;%s;%s" % ('<{0.__class__.__name__}: [{0.name}, {0.number}, {0.element}]:({0.x:.3f}, {0.y:.3f}, {0.z:.3f})>'.format(self._A), '<{0.__class__.__name__}: [{0.name}, {0.number}, {0.element}]:({0.x:.3f}, {0.y:.3f}, {0.z:.3f})>'.format(self._B), self._distance, '<{0.__class__.__name__}: [{0.type}, {0.number}, {0.version}]>'.format(self._A_residue_obj), '<{0.__class__.__name__}: [{0.type}, {0.number}, {0.version}]>'.format(self._B_residue_obj), self._A_chain, self._B_chain)
        except:
            return "%s;%s;%s;%s;%s;%s;%s" % (self._A, self._B, self._distance, self._A_residue_obj, self._B_residue_obj, self._A_chain, self._B_chain)

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python contacts.py -i input_file [-c contact_type -d distance_type --dummy=dummy_dir -o output_file]")

    parser.add_option("-c", default="pdi", action="store", type="string", dest="contact_type", help="Contact type (i.e. \"pdi\" or \"ppi\"; default = pdi)", metavar="{string}")
    parser.add_option("-d", default="dinucleotides", action="store", type="string", dest="distance_type", help="Distance type (i.e. \"basepairs\", \"dinucleotides\" or \"mindist\"; default = dinucleotides)", metavar="{string}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-f", "--filter",  default=False, action="store_true", dest="filter", help="Filter non-standard contacts (only if \"-d\" option is \"mindist\")")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in PDB format)", metavar="{filename}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")

    (options, args) = parser.parse_args()

    if options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")
    if options.contact_type != "pdi" and options.contact_type != "ppi":
        parser.error("incorrect contact type: accepted values are \"pdi\" (for protein-DNA interactions), \"ppi\" (for protein-protein interactions)")
    if options.distance_type != "basepairs" and options.distance_type != "dinucleotides" and options.distance_type != "mindist":
        parser.error("incorrect distance type: accepted values are \"basepairs\" (for contact distance between amino acid cb and geometric center of a basepair), \"dinucleotides\" (for contact distance between amino acid cb and geometric center of a dinucleotide) and \"mindist\" (for minimum distance)")

    return options

def get_contacts_obj(pdb_obj, x3dna_obj=None, contacts_type="pdi", distance_type="dinucleotides", dummy_dir="/tmp"):
    """
    This function extracts all protein-DNA/protein contacts from a PDB file and
    returns a {Contacts} object.

    @input:
    pdb_obj {PDB}
    x3dna_obj {X3DNA}
    contacts_type {string} either "pdi" or "ppi"
    distance_type {string} either "basepairs", "dinucleotides" or "mindist"
    dummy_dir {string}

    @return:
    contacts_obj {Contacts}

    """

    # Initialize #
    done = set()
    contacts_obj = Contacts()
    max_contact_distance = float(config.get("Parameters", "max_contact_distance"))
    if contacts_type == "ppi" or distance_type == "mindist": distance_threshold = max_contact_distance
    else: distance_threshold = max_contact_distance + max_contact_distance / 2.0

    # For each protein chain... #
    for protein_chain_obj in pdb_obj.proteins:
        # For each amino acid... #
        for aminoacid_obj in protein_chain_obj.aminoacids:
            # If protein-DNA interactions... #
            if contacts_type == "pdi":
                # For each nucleotide chain... #
                for dna_chain_obj in pdb_obj.nucleotides:
                    # For each nucleotide... #
                    for nucleotide_obj in dna_chain_obj.nucleotides:
                        # Get geometric center distance between amino acid and nucleotide #
                        a, b, distance = aminoacid_obj.distance(nucleotide_obj, "geometric")
                        # Skip if residues are too far away #
                        if distance >= distance_threshold: continue
                        # If min. distance... #
                        if distance_type == "mindist":
                            # For each atom... #
                            for aminoacid_atom_obj in aminoacid_obj.atoms:
                                # Skip if amino acid atom from peptidic bond #
                                if aminoacid_atom_obj.is_Calpha or aminoacid_atom_obj.is_N or aminoacid_atom_obj.is_O: continue
                                # For each atom... #
                                for dna_atom_obj in nucleotide_obj.atoms:
                                    # Allow contacts #
                                    proceed = False
                                    if aminoacid_atom_obj.element == "O" and dna_atom_obj.element == "N": proceed = True
                                    if aminoacid_atom_obj.element == "O" and dna_atom_obj.element == "C": proceed = True
                                    if aminoacid_atom_obj.element == "N" and dna_atom_obj.element == "O": proceed = True
                                    if aminoacid_atom_obj.element == "C" and dna_atom_obj.element == "C": proceed = True
                                    if not proceed: continue
                                    # Skip if already done #
                                    if ((protein_chain_obj.chain, aminoacid_obj.number, aminoacid_atom_obj), (dna_chain_obj.chain, nucleotide_obj.number, dna_atom_obj)) in done: continue
                                    # Get distance #
                                    distance = aminoacid_atom_obj.distance(dna_atom_obj)
                                    # Add contact #
                                    if distance <= max_contact_distance:
                                        contacts_obj.add_contact(Contact(aminoacid_atom_obj, dna_atom_obj, distance, aminoacid_obj, nucleotide_obj, protein_chain_obj.chain, dna_chain_obj.chain))
                                    done.add(((protein_chain_obj.chain, aminoacid_obj.number, aminoacid_atom_obj), (dna_chain_obj.chain, nucleotide_obj.number, dna_atom_obj)))
                        # Else... #
                        else:
                            # Initialize #
                            aminoacid_atom_obj = get_aminoacid_cb_or_ca(aminoacid_obj)
                            phosphate_atoms = []
                            dna_residues = []
                            dna_chains = []
                            basepair = x3dna_obj.get_residue_basepair(dna_chain_obj.chain, nucleotide_obj.number)
                            # Skip if aminoacid atom is None #
                            if aminoacid_atom_obj is None: continue
                            # Skip if basepair is None #
                            if basepair is None: continue
                            # If using basepairs... #
                            if distance_type == "basepairs":
                                for pdb_chain, residue_num in x3dna_obj.get_basepair(basepair):
                                    phosphate_atoms.append(get_nucleotide_p_or_bb(pdb_obj.get_chain_by_id(pdb_chain).get_residue_by_identifier(str(residue_num))))
                                    dna_residues.append(pdb_obj.get_chain_by_id(pdb_chain).get_residue_by_identifier(str(residue_num)))
                                    dna_chains.append(pdb_chain)
                                # Skip if not enough contacts #
                                if len(phosphate_atoms) != 2: continue
                            # If using dinucleotides... #
                            if distance_type == "dinucleotides":
                                for dinucleotide in x3dna_obj.get_basepair_dinucleotides(basepair):
                                    for basepair in x3dna_obj.get_dinucleotide(dinucleotide):
                                        for chain, residue_num in x3dna_obj.get_basepair(basepair):
                                            phosphate_atoms.append(get_nucleotide_p_or_bb(pdb_obj.get_chain_by_id(chain).get_residue_by_identifier(str(residue_num))))
                                            dna_residues.append(pdb_obj.get_chain_by_id(chain).get_residue_by_identifier(str(residue_num)))
                                            dna_chains.append(chain)
                                    # Do the 1st dinucleotide only (eventually will do all) #
                                    break
                                # Skip if not enough contacts #
                                if len(phosphate_atoms) != 4: continue
                            # Skip if already done #
                            if ((protein_chain_obj.chain, aminoacid_obj.number, aminoacid_atom_obj), (tuple(dna_chains), tuple(dna_residues), tuple(phosphate_atoms))) in done: continue
                            # Get distance to geometric center #
                            distance = aminoacid_atom_obj.distance_to_point(numpy.array([sum([i.x for i in phosphate_atoms]) / len(phosphate_atoms), sum([i.y for i in phosphate_atoms]) / len(phosphate_atoms), sum([i.z for i in phosphate_atoms]) / len(phosphate_atoms)]))
                            # Add contact #
                            if distance <= max_contact_distance:
                                contacts_obj.add_contact(Contact(aminoacid_atom_obj, phosphate_atoms, distance, aminoacid_obj, dna_residues, protein_chain_obj.chain, dna_chains))
                            done.add(((protein_chain_obj.chain, aminoacid_obj.number, aminoacid_atom_obj), (tuple(dna_chains), tuple(dna_residues), tuple(phosphate_atoms))))
            # If protein-protein interactions... #
            if contacts_type == "ppi":
                # For each protein chain... #
                for other_protein_chain_obj in pdb_obj.proteins:
                    # Skip if it is the same protein PDB chain #
                    if other_protein_chain_obj.chain == protein_chain_obj.chain: continue
                    # For each amino acid... #
                    for other_aminoacid_obj in other_protein_chain_obj.aminoacids:
                        # Get geometric center distance between amino acid and nucleotide #
                        a, b, distance = aminoacid_obj.distance(other_aminoacid_obj, "geometric")
                        # Skip if residues are too far away #
                        if distance >= distance_threshold: continue
                        # For each atom... #
                        for aminoacid_atom_obj in aminoacid_obj.atoms:
                            # Skip if amino acid atom from peptidic bond #
                            if aminoacid_atom_obj.is_Calpha or aminoacid_atom_obj.is_N or aminoacid_atom_obj.is_O: continue
                            # For each atom... #
                            for other_aminoacid_atom_obj in other_aminoacid_obj.atoms:
                                # Skip if amino acid atom from peptidic bond #
                                if other_aminoacid_atom_obj.is_Calpha or other_aminoacid_atom_obj.is_N or other_aminoacid_atom_obj.is_O: continue
                                # Allow contacts #
                                proceed = False
                                if aminoacid_atom_obj.element == "S" and other_aminoacid_atom_obj.element == "S": proceed = True
                                if aminoacid_atom_obj.element == "O" and other_aminoacid_atom_obj.element == "N": proceed = True
                                if aminoacid_atom_obj.element == "N" and other_aminoacid_atom_obj.element == "O": proceed = True
                                if aminoacid_atom_obj.element == "C" and other_aminoacid_atom_obj.element == "C": proceed = True
                                if not proceed: continue
                                # Skip if already done #
                                if ((protein_chain_obj.chain, aminoacid_obj.number, aminoacid_atom_obj), (other_protein_chain_obj.chain, other_aminoacid_obj.number, other_aminoacid_atom_obj)) in done: continue
                                # Get distance #
                                distance = aminoacid_atom_obj.distance(other_aminoacid_atom_obj)
                                # Add contact #
                                if distance <= max_contact_distance:
                                    contacts_obj.add_contact(Contact(aminoacid_atom_obj, other_aminoacid_atom_obj, distance, aminoacid_obj, other_aminoacid_obj, protein_chain_obj.chain, other_protein_chain_obj.chain))
                                done.add(((protein_chain_obj.chain, aminoacid_obj.number, aminoacid_atom_obj), (other_protein_chain_obj.chain, other_aminoacid_obj.number, other_aminoacid_atom_obj)))
                                done.add(((other_protein_chain_obj.chain, other_aminoacid_obj.number, other_aminoacid_atom_obj), (protein_chain_obj.chain, aminoacid_obj.number, aminoacid_atom_obj)))
        # Add protein chain to done #
        done.add(protein_chain_obj.chain)
    return contacts_obj

def get_aminoacid_cb_or_ca(aminoacid_obj):
    """
    This function returns the CB atom or the CA atom of an amino acid;
    otherwise returns "None".

    @input:
    aminoacid_obj {ResidueOfAminoAcid}
    @return:
    atom_obj {AtomOfAminoAcid}

    """

    # Initialize #
    aminoacid_atom_obj = None

    if aminoacid_obj.has_cb:
        aminoacid_atom_obj = aminoacid_obj.cb
    elif aminoacid_obj.has_ca: # for glycines
        aminoacid_atom_obj = aminoacid_obj.ca

    return aminoacid_atom_obj

def get_nucleotide_p_or_bb(nucleotide_obj):
    """
    This function returns the phosphate atom or the the 1st backbone atom
    of a nucleotide; otherwise returns "None".

    @input:
    nucleotide_obj {ResidueOfNucleotide}
    @return:
    atom_obj {AtomOfNucleotide}

    """

    # Initialize #
    dna_atom_obj = None

    if nucleotide_obj.has_p:
        dna_atom_obj = nucleotide_obj.p
    else:
        dna_atom_obj = nucleotide_obj._backbone_atoms[0]

    return dna_atom_obj

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
    contacts_obj = get_contacts_obj(pdb_obj, x3dna_obj, options.contact_type, options.distance_type, os.path.abspath(options.dummy_dir))

    # Output #
    if options.output_file is not None:
        contacts_obj.write(os.path.abspath(options.output_file), options.distance_type, options.filter)
    else:
        sys.stdout.write("#atomA;atomB(s);distance;residueA;residueB(s);chainA;chainB(s)\n")
        for contact_obj in contacts_obj.get_contacts():
            if options.distance_type == "mindist" and options.filter:
                if not contact_obj.is_disulfide_bridge() and not contact_obj.is_hydrogen_bond() and not contact_obj.is_salt_bridge() and not contact_obj.is_van_der_waals():
                    continue
            sys.stdout.write("%s\n" % contact_obj.return_as_string())
