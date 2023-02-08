import os, sys, re
import optparse
import functions
import ConfigParser
import threader
import blast, contacts, dssp, hmmer, interface, model_protein, nr, triads, x3dna, homologs
scripts_path = os.path.abspath(os.path.dirname(__file__))
import itertools as it
# Append scripts path to python path #
sys.path.append(scripts_path)

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Imports jbonet's module #
from SBI.data import aminoacids1to3, aminoacids_polarity_boolean, nitrogenous_bases, dna_complementary
from SBI.structure import PDB

def  get_triads_from_folder(thr_obj,pdb_structure,pdb_dir):
    # get pdb_name and chain of the template#
    dna_sequence={}
    pdb_name, pdb_chain = thr_obj._pdb_name, thr_obj._pdb_chain
    pdb_structure = os.path.join(pdb_name+"_"+pdb_chain)
    if not os.path.exists(os.path.join(pdb_dir, "helices", pdb_structure + ".txt")):
       sys.stdout.write("ERROR: file %s is not found. The PDB folder MUST contain the PDB of the thread"%(os.path.join(pdb_dir, "helices", pdb_structure + ".txt")))
    # get helix structure #
    for helix in functions.parse_file(os.path.join(pdb_dir, "helices", pdb_structure + ".txt")):
        # Get PDB object #
        pdb_obj = PDB(os.path.join(pdb_dir, "split", pdb_name + ".dna." + helix + ".pdb"))
        break
    pdb_obj_prot = PDB(os.path.join(pdb_dir, "split", pdb_structure + ".pdb"))
    pdb_obj.add_chains(pdb_obj_prot.chains)
    for chain_id in pdb_obj.chain_identifiers:
        chain=pdb_obj.get_chain_by_id(chain_id)
        if chain.chaintype=="N":
            dna_sequence.setdefault(chain_id,chain.nucleotide_sequence())
    # Get X3DNA object #
    x3dna_obj = x3dna.X3DNA(os.path.join(pdb_dir, "x3dna", pdb_name + ".txt"))

    # Get contacts object #
    contacts_obj = contacts.Contacts(os.path.join(pdb_dir, "contacts",pdb_name + ".txt"))

    # Skip if no contacts #
    if len(contacts_obj._contacts) == 0: raise ValueError("No protein-DNA contacts found!")


    # Get helix #
    dna_helices={}
    for helix in x3dna_obj.get_dna_helices():
        # Get helix dinucleotides #
        dinucleotides = x3dna_obj.get_helix_dinucleotides(helix)
        # For each contact... #
        for contact_obj in contacts_obj.get_contacts():
            if contact_obj._A_chain == pdb_chain:
                dna_helices.setdefault(helix,set()).add(tuple(contact_obj._B_chain))

    # Get interface object #
    #interface_obj = interface.Interface(os.path.join(pdb_dir, "interfaces", pdb_structure + ".txt"))
    #interface_basepairs = interface_obj.get_interface_basepairs()
    #for dna_helix in dna_helices.iterkeys():
    #    dna_sequence_interface = dna_sequence[list(dna_helices[dna_helix])[0][0]][interface_obj.get_start()-1:interface_obj.get_interface_length()+interface_obj.get_start()+1]
    # Get nucleotide positions #

    # Get triads object #
    triads_obj = triads.Triads(os.path.join(pdb_dir, "triads", pdb_name + "_" + pdb_chain + ".txt"))

    return triads_obj, pdb_obj, x3dna_obj

def get_triads_from_template(template,verbose,dummy_dir):

    # Get PDB object #
    if verbose: sys.stdout.write("\t-- Reading PDB file %s ...\n"%(os.path.basename(template)))
    pdb_obj = PDB(template)

    # Get DSSP object #
    if verbose: sys.stdout.write("\t\t-- Get DSSP ...\n")
    dssp_obj = dssp.get_dssp_obj(template)

    # Get X3DNA object #
    if verbose: sys.stdout.write("\t\t-- Get DNA object ...\n")
    x3dna_obj = x3dna.get_x3dna_obj(template, dummy_dir)

    # Get contacts object #
    if verbose: sys.stdout.write("\t\t-- Get contacts ...\n")
    contacts_obj = contacts.get_contacts_obj(pdb_obj, x3dna_obj)

    # Skip if no contacts #
    if len(contacts_obj._contacts) == 0: raise ValueError("No protein-DNA contacts found!")

    # Get triads object #
    triads_obj = triads.get_triads_obj(pdb_obj, dssp_obj, x3dna_obj, contacts_obj)


    return triads_obj, pdb_obj, x3dna_obj


def thread_kmer_triads(kmer,triads_obj):
    # Initialize #
    full_dna_seq = list(it.islice(kmer,0,None))
    full_dna_pos = list(range(1,len(full_dna_seq)+1))
    basepairs    = dict(zip(full_dna_pos,full_dna_seq))
    thread_triads_obj = triads.Triads()
    kmer = list(kmer)
    # For each triad object... #
    for triad_obj in triads_obj.get_triads():
   	    # Initialize #
        aminoacid_environment = []
        dinucleotide_environment = []
        a_oa, b_ob, distance, residue_A, residue_B = triad_obj.return_as_string().split(';')
        # Get chain, number #
        chain, number = residue_A.split("-")
        nucleotide_environment = residue_B.split(",")
        dinucleotide_environment = b_ob.split("-")
        # Information stored in b_ob  
        # b, twobases, dna_strand, dna_groove, dna_chemical_group = b_ob.split('-')
        # Get aa_environment and change its residues according to sequence on strand in F2_aa #
        # Forward strand #
        if dinucleotide_environment[2] == 'F':
           if int(nucleotide_environment[0].split("-")[1]) in basepairs.keys() and int(nucleotide_environment[2].split("-")[1]) <= max(basepairs.keys()):
              aminoacid_environment = a_oa.split("-")
              # Information stored in a_oa  
              # a, hydrophobicity, degree_of_exposure, secondary_structure = a_oa.split("-")
              # is preserved
              dinucleotide_environment[0] = basepairs[int(nucleotide_environment[0].split("-")[1])] + basepairs[int(nucleotide_environment[2].split("-")[1])]
              dinucleotide_environment[1] = "".join(nitrogenous_bases[nucleotide] for nucleotide in dinucleotide_environment[0])
              triad_obj._A_environment,triad_obj._B_environment = "-".join(aminoacid_environment), "-".join(dinucleotide_environment)
              thread_triads_obj.add_triad(triad_obj)
        # Reverse strand #
        if dinucleotide_environment[2] == 'R':
           if int(nucleotide_environment[0].split("-")[1]) in basepairs.keys() and int(nucleotide_environment[2].split("-")[1]) <= max(basepairs.keys()):
              aminoacid_environment = a_oa.split("-")
              # Information stored in a_oa  
              # a, hydrophobicity, degree_of_exposure, secondary_structure = a_oa.split("-")
              # is preserved
              dinucleotide_environment[0] = dna_complementary[basepairs[int(nucleotide_environment[0].split("-")[1])]] + dna_complementary[basepairs[int(nucleotide_environment[2].split("-")[1])]]
              dinucleotide_environment[1] = "".join(nitrogenous_bases[nucleotide] for nucleotide in dinucleotide_environment[0])
              triad_obj._A_environment,triad_obj._B_environment = "-".join(aminoacid_environment), "-".join(dinucleotide_environment)
              thread_triads_obj.add_triad(triad_obj)
    # Return  new threaded triads
    return thread_triads_obj


def threading_triads(threading_file, pdb_dir=None, template=None, verbose=True, dummy_dir="/tmp"):
    thr_obj = threader.Threaded(threading_file=threading_file)
    # add threading data to object #
    thr_obj._check_parsing()
    pdb_name, pdb_chain = thr_obj._pdb_name, thr_obj._pdb_chain
    pdb_structure = os.path.join(pdb_name+"_"+pdb_chain)
    # Get nucleotide positions #
    full_dna_seq = list(it.islice(thr_obj._dna.keys()[0],0,None))
    position=1 + int(thr_obj.get_kmer(thr_obj._dna.keys()[0]))
    full_dna_pos = list(range(position,len(full_dna_seq)+position))
    dna_int_obj = dict(zip(full_dna_pos,full_dna_seq))

    #if template is None, use pdb_dir
    if template is None:
       try:
         triads_obj, pdb_obj, x3dna_obj  = get_triads_from_folder(thr_obj,pdb_structure,pdb_dir)
       except ValueError as e:
         if verbose: sys.stderr.write("Failed get triads from folder %s\n"%e)
         exit(0)
    else:
       try:
         triads_obj, pdb_obj, x3dna_obj  = get_triads_from_template(os.path.abspath(template),verbose,dummy_dir)
       except ValueError as e:
         if verbose: sys.stderr.write("Failed get triads from template %s\n"%e)
         exit(0)


    # Initialize #
    protein = {}
    kmers = []
    for line in functions.parse_file(threading_file):
        if line == "": continue
        if line.startswith("#") or line.startswith("//"): continue
        elif line.startswith(">"):
            threading = line[1:]
        elif threading == "protein":
            line = line.split(";")
            protein.setdefault((line[0], int(line[1])), line[2])
        elif threading == "dna":
            kmers.append(line.split(";"))
    # For each k-mer... #
    for kmer, position in kmers:
        # Initialize #
        thread_triads_obj = triads.Triads()
        basepairs = {}
        kmer = list(kmer)
        # For each contact... #
        for pos in dna_int_obj.keys():
            basepairs.setdefault(pos, dna_int_obj[pos])
        # Skip if triads file already exists #
        # Initialize #
        done = set()
        # For each triad object... #
        for triad_obj in triads_obj.get_triads():
        	# Initialize #
            aminoacid_environment = []
            dinucleotide_environment = []
            a_oa, b_ob, distance, residue_A, residue_B = triad_obj.return_as_string().split(';')
            if residue_A.split("-")[0] == pdb_chain:
                # Get chain, number #
                chain, number = residue_A.split("-")
                nucleotide_environment = residue_B.split(",")
                dinucleotide_environment = b_ob.split("-")
                # Information stored in b_ob  
                # b, twobases, dna_strand, dna_groove, dna_chemical_group = b_ob.split('-')
                # Get aa_environment and change its residues according to sequence on strand in F2_aa #
                # Forward strand #
                if (pdb_chain,int(number)) not in protein.keys(): continue
                else: pass
                if dinucleotide_environment[2] == 'F':
                    if int(nucleotide_environment[0].split("-")[1]) in basepairs.keys() and int(nucleotide_environment[2].split("-")[1]) <= max(basepairs.keys()):
                            aminoacid_environment = a_oa.split("-")
                            # Information stored in a_oa  
                            #a, hydrophobicity, degree_of_exposure, secondary_structure = a_oa.split("-")
                            aminoacid_environment[0] = aminoacids1to3[protein[(pdb_chain,int(number))]]
                            if aminoacids_polarity_boolean[protein[(pdb_chain,int(number))]]:
                                    hydrophobicity=aminoacid_environment[1] = "P"
                            else:
                                    hydrophobicity=aminoacid_environment[1] = "N"
                            dinucleotide_environment[0] = basepairs[int(nucleotide_environment[0].split("-")[1])] + basepairs[int(nucleotide_environment[2].split("-")[1])]
                            dinucleotide_environment[1] = "".join(nitrogenous_bases[nucleotide] for nucleotide in dinucleotide_environment[0])
                            triad_obj._A_environment,triad_obj._B_environment = "-".join(aminoacid_environment), "-".join(dinucleotide_environment)
                            thread_triads_obj.add_triad(triad_obj)

                # Reverse strand #
                if dinucleotide_environment[2] == 'R':
                    if int(nucleotide_environment[0].split("-")[1]) in basepairs.keys() and int(nucleotide_environment[2].split("-")[1]) <= max(basepairs.keys()):
                            aminoacid_environment = a_oa.split("-")
                            # Information stored in a_oa  
                            #a, hydrophobicity, degree_of_exposure, secondary_structure = a_oa.split("-")
                            aminoacid_environment[0] = aminoacids1to3[protein[(pdb_chain,int(number))]]
                            if aminoacids_polarity_boolean[protein[(pdb_chain,int(number))]]:
                                    hydrophobicity=aminoacid_environment[1] = "P"
                            else:
                                    hydrophobicity=aminoacid_environment[1] = "N"
                            dinucleotide_environment[0] = dna_complementary[basepairs[int(nucleotide_environment[0].split("-")[1])]] + dna_complementary[basepairs[int(nucleotide_environment[2].split("-")[1])]]
                            dinucleotide_environment[1] = "".join(nitrogenous_bases[nucleotide] for nucleotide in dinucleotide_environment[0])
                            triad_obj._A_environment,triad_obj._B_environment = "-".join(aminoacid_environment), "-".join(dinucleotide_environment)
                            thread_triads_obj.add_triad(triad_obj)
    return thread_triads_obj, pdb_obj, x3dna_obj


def parse_options():
    parser = optparse.OptionParser("threading_to_triads.py --threading <input_file>")
    parser.add_option("--threading", action="store", type="string", dest="threading_file",  default=None, help="path to threading file", metavar="THREADING_FILE")
    parser.add_option("-o", default="threaded_triads.dat", action="store", type="string", dest="output", help="Output file (default = threaded_triad.dat)", metavar="{file}")
    parser.add_option("--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (i.e. output dir from pdb.py)", metavar="{directory}")

    (options, args) = parser.parse_args()
    return options

if __name__ == "__main__":

    # Arguments & Options #
    options  = parse_options()

    threading_file = options.threading_file
    triad_obj, pdb_obj, x3dna_obj = threading_triads(threading_file, pdb_dir= options.pdb_dir)
    # create the threaded object and get the threading file #
    triad_obj.write(options.output)
