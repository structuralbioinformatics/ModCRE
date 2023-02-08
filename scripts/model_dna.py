import os, sys, re
import ConfigParser
import optparse
import subprocess

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

# Imports jbonet's module #
from SBI.data import dna_complementary
from SBI.structure import PDB
from SBI.structure.chain import ChainOfNucleotide
from SBI.structure.residue import ResidueOfNucleotide

# Import my modules #
import x3dna, interface

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python model_dna.py  -p pdb_file ( -i interface_file -p pdb_file and (-s dna_sequence or (-t threading_file --pdb pdb_dir ) ) ][--dummy=dummy_dir -o output_dir]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", default=None, type="string", dest="interface_file", help="Interface file (from interface.py)", metavar="{filename}")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="{directory}")
    parser.add_option("-p", action="store", type="string", dest="pdb_file", help="PDB file (e.g. from model_protein.py)", metavar="{filename}")
    parser.add_option("-s", action="store", type="string", dest="dna_sequence", help="DNA sequence ", metavar="{string}")
    parser.add_option("-t", action="store", default=None, type="string", dest="threading_file", help="Threading file (e.g. from threader.py)", metavar="{string}")
    parser.add_option("--pdb", action="store", default=None, type="string", dest="pdb_dir", help="PDB directory (i.e. output dir from pdb.py)", metavar="PDB_DIR")
 

    (options, args) = parser.parse_args()

    if (options.interface_file is None or options.pdb_file is None or (options.dna_sequence is None and (options.threading_file is None or options.pdb_dir is None)) ):
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_dna_model_pdb_obj(pdb_file, dna_sequence, x3dna_obj, interface_obj, interface_range=None, dummy_dir="/tmp"):

    # Initialize #
    basepairs = set()
    nucleobases = set()
    mutations_file = os.path.join(dummy_dir, str(os.getpid()) + ".mutations.dat")
    model_pdb_file = os.path.join(dummy_dir, str(os.getpid()) + "mutated.pdb")
    # Correct DNA sequence #
    dna_sequence = re.sub("U", "T", dna_sequence)
    dna_sequence = re.sub("[^acgtACGT]", "N", dna_sequence)
    # Get nucleotides #
    nucleotides = list(dna_sequence)
    # For each interface basepair... #
    if  interface_obj.get_interface_basepairs() is not None:
        interface_range=interface_obj.get_interface_basepairs()
    elif interface_range is None:
        sys.stdout.write("[WARNING] No DNA interface to model or unrecognized DNA molecule. Original DNA is preserved\n")
        return PDB(pdb_file)
    for basepair in interface_range:
      if len(nucleotides)>0:
        # Initialize #
        nucleotide = nucleotides.pop(0).upper()
        # Get basepair nucleotide #
        ((pdb_chain_A, residue_num_A), (pdb_chain_B, residue_num_B)) = x3dna_obj.get_basepair(basepair)
        # Mutate the forward strand #
        if nucleotide != "N":
            functions.write(mutations_file, "chain=%s seqnum=%s mutation=%s" % (pdb_chain_A, str(residue_num_A), nucleotide))
            # Keep nucleobase #
            nucleobases.add((pdb_chain_A, residue_num_A))
        # Mutate the reverse strand #
        if nucleotide != "N":
            functions.write(mutations_file, "chain=%s seqnum=%s mutation=%s" % (pdb_chain_B, str(residue_num_B), dna_complementary[nucleotide]))
            # Keep nucleobase #
            nucleobases.add((pdb_chain_B, residue_num_B))
        # Keep basepair #
        basepairs.add(basepair)
    # Execute 3DNA #
    try:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        x3dna_path = os.path.join(src_path, config.get("Paths", "x3dna_path"))
        os.environ['X3DNA'] = x3dna_path[:-4]
        # Exec process #
        process = subprocess.check_output([os.path.join(x3dna_path, "mutate_bases"), "-l", mutations_file, pdb_file, model_pdb_file], stderr=subprocess.STDOUT, env=os.environ)
    except:
        raise ValueError("Could not exec X3DNA for %s" % pdb_file)
    # Get raw model PDB object #
    model_pdb_obj = PDB(model_pdb_file)
    # For each PDB chain... #
    for i in range(len(model_pdb_obj.chains)):
        # Skip if not DNA chain #
        if model_pdb_obj.chains[i].chaintype != "N": continue
        # Initialize #
        pdb_chain_obj = ChainOfNucleotide(model_pdb_obj.chains[i]._pdb, model_pdb_obj.chains[i]._chain)
        # For each nucleotide... #
        for nucleotide_obj in model_pdb_obj.chains[i].nucleotides:
            # If allowed basepair... #
            if x3dna_obj.get_residue_basepair(model_pdb_obj.chains[i].chain, int(nucleotide_obj.number)) in basepairs:
                # If allowed nucleobase... #
                if (model_pdb_obj.chains[i].chain, int(nucleotide_obj.number)) in nucleobases:
                    pdb_chain_obj.add_residue(nucleotide_obj)
                # Else... #
                else:
                    pdb_chain_obj.add_residue(ResidueOfNucleotide(nucleotide_obj._number, nucleotide_obj._version, nucleotide_obj._type, nucleotide_obj._mode))
                    for backbone_atom_obj in nucleotide_obj.backbone:
                        pdb_chain_obj.nucleotides[-1].add_atom(backbone_atom_obj)
        model_pdb_obj.chains[i] = pdb_chain_obj
    # Erase tmp files #
    os.remove(mutations_file)
    os.remove(model_pdb_file)
    return model_pdb_obj

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Create output directory #
    if not os.path.exists(options.output_dir):
        os.makedirs(options.output_dir)
    if not os.path.exists(options.dummy_dir):
        os.makedirs(options.dummy_dir)
    pdb_file = options.pdb_file
    if not pdb_file.startswith("/"): pdb_file = os.path.abspath(pdb_file)
    pdb_dir = options.pdb_dir
    if pdb_dir is not None:
       if not pdb_dir.startswith("/"): pdb_dir = os.path.abspath(options.pdb_dir)
    if options.threading_file is not None:
       threading_file=options.threading
       if not threading_file.startswith("/"): threading_file=os.path.abspath(threading_file)
    if options.interface_file is not None:
       interface_file=options.interface_file
       if not interface_file.startswith("/"): interface_file=os.path.abspath(interface_file) 

    # Define DNA sequences
    dna=set()
    kmers={}

    # Get X3DNA object #
    x3dna_obj = x3dna.get_x3dna_obj(pdb_file)

    # Get interface object for threading object#
    if options.threading_file is not None and options.pdb_dir is not None:
     if os.path.exists(threading_file) and os.path.exists(pdb_dir):
       thread_obj = threader.Threaded(threading_file)
       pdb_name   = thread_obj.get_pdb_name()
       pdb_chain  = thread_obj.get_pdb_chain()
       interface_file = os.path.join(pdb_dir,"interfaces",pdb_name+"_"+ pdb_chain+".txt")
       kmers      = thread_obj.get_kmers()
       kmers_fixed= thread_obj.get_kmers_fixed()
       for dna_sequence in kmers.iterkeys():
           if "N" not in dna_sequence:
               dna.add(dna_sequence)
       for dna_sequence in kmers_fixed.iterkeys():
               dna.add(dna_sequence)
     else:
       sys.stderr.write("Threading file %s is not found\n"%(threading_file))
       exit(0)
    else:
       kmers.setdefault(options.dna_sequence,0)
       dna.add(options.dna_sequence)
    if os.path.exists(interface_file):
       interface_obj = interface.Interface(interface_file)
    else:
       sys.stderr.write("Interface file %s is not found\n"%(interface_file))
       exit(0)


    # For all dna sequences

    for dna_sequence in dna:

      # If DNA sequence does not cover the whole interface... #
      if len(dna_sequence) != interface_obj.get_interface_length():
        raise ValueError("DNA sequence does not match the protein-DNA interface. Try \"-s %s\" (please replace \"N\" by any nucleotide {A, C, G, T} or your choice)." % ("N" * interface_obj.get_interface_length()))

      # Get model object #
      pdb_obj = get_dna_model_pdb_obj(pdb_file, dna_sequence, x3dna_obj, interface_obj, interface_range=None, dummy_dir=options.dummy_dir)

      # Write model object #
      pdb_obj.write(os.path.join(options.output_dir, "model." + options.dna_sequence + ".pdb"))

