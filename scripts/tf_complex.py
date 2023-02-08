import os, sys, re
import ConfigParser
import json
import optparse
import shutil
import subprocess
import imp
import math
import pickle
import numpy as np

from SBI.structure import PDB
from SBI.structure.chain import ChainOfProtein, ChainOfNucleotide
from SBI.data import dna_complementary

import x3dna, functions

# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.dirname(__file__))

# Append scripts path to python path #
sys.path.append(scripts_path)

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Get python path #
python = os.path.join(config.get("Paths", "python_path"), "python")

def get_cooperative_model(coop, coords, coop_dir, dna_seq, dna_file, dna_str, dummy_dir, beggining=1, model_label="", verbose=False):

    list_of_models = coop.split(";")
    output_dir = os.path.dirname(coop_dir)
    
    nuc_dna = PDB(dna_file)
    only_protein_structures = []
    #print("list of models: " + str(list_of_models))
    for i in range(0, len(list_of_models)):
        # Get the DNA structure that will be superimposed #
        model = list_of_models[i]
        coord_1 = int(coords[i].split(":")[0])
        coord_2 = int(coords[i].split(":")[1])
        new_nuc_dna = PDB()
        chain_1 = ChainOfNucleotide("dna", "A")
        chain_2 = ChainOfNucleotide("dna", "B")
        # Get the forward and the reverse chain from the DNA #
        forward_chain = None
        reverse_chain = None
        for chain in nuc_dna.chains:
            if isinstance(chain, ChainOfNucleotide):
                if forward_chain == None:
                    forward_chain = chain
                elif reverse_chain == None:
                    reverse_chain = chain
                else:
                    break
        # Get nucleotides from both chains #
        bs_length = coord_2 - coord_1
        #first_nuc_rev = reverse_chain.first_nucleotide
        for i in range(coord_1, coord_2+1):
            try:
                chain_1.add_residue(forward_chain.nucleotides[i])
            except:
                continue
        for i in sorted(range(coord_1, coord_2+1), reverse=True):
            try:
                chain_2.add_residue(reverse_chain.nucleotides[-i-1])
            except:
                continue

        new_nuc_dna.add_chain(chain_1)
        new_nuc_dna.add_chain(chain_2)
        new_nuc_dna.write(output_file=os.path.join(dummy_dir, "dna_" + os.path.basename(model)), force=True)
        if verbose: print("\n\nnew dna structure created at: " + os.path.join(dummy_dir, "dna_" + os.path.basename(model)) + "\n\n")
        # Execute DNA superimposer #
        label = model[:-4] + "_superimposed"
        model = os.path.join(output_dir, "binary_interactions", model)
        if verbose: print("python " + os.path.join(scripts_path, "dna_superimposer.py") + " --dummy=" + dummy_dir + " -a " + dna_file + " -b " + model + " -o " + os.path.join(output_dir, "working", "cooperative") + " -l " + label + " --dna=" + dna_seq)
        os.system("python " + os.path.join(scripts_path, "dna_superimposer.py") + " --dummy=" + dummy_dir + " -a " + dna_file + " -b " + model + " -o " + os.path.join(output_dir, "working", "cooperative") + " -l " + label + " --dna=" + dna_seq)
        if os.path.isfile(os.path.join(output_dir, "working", "cooperative", label + ".pdb")):
            superimposition = os.path.join(output_dir, "working", "cooperative", label + ".pdb")
            if verbose: print("\n\nSuperimposed file in: " + superimposition + "\n\n")
            # From the superimposed structure isolate the protein chains #
            sup_pdb = PDB(superimposition)
            new_sup = PDB()
            for chain in sup_pdb.chains:
                if isinstance(chain, ChainOfProtein):
                    new_sup.add_chain(chain)
            new_sup.write(output_file=os.path.join(output_dir, "working", "cooperative", label + "_only_protein.pdb"), force=True)
            only_protein_structures.append(os.path.join(output_dir, "working", "cooperative", label + "_only_protein.pdb"))
    # Check clashes #
    #clashes_dict = get_clashes(only_protein_structures, dummy_dir)
    clashes_dict = {}
    # Write clashes report #
    # Join all structures into the same structure #
    chain_names = list("QWERTYUIOPASDFGHJKLZXCVBNM")
    final_model = PDB(dna_file)
    used_chains = []
    for chain in final_model.chains:
        used_chains.append(chain.chain)
    for ps in only_protein_structures:
        prot = PDB(ps)
        for chain in prot.chains:
            if isinstance(chain, ChainOfProtein):
                for cn in chain_names:
                    if not cn in used_chains:
                        chain.chain = cn
                        final_model.add_chain(chain)
                        used_chains.append(cn)
                        #if verbose: print("\n\nAdding to the model structure " + os.path.basename(ps) + " with chain id " + cn + "\n\n")
                        break
    if model_label == "":
        final_model.write(output_file=os.path.join(output_dir, "combined_models", "model.pdb"), force=True)
        if verbose: 
            print("Model finished in: " + os.path.join(output_dir, "combined_models", "model.pdb"))
    else:
        print("cooperative modeling label: " + str(model_label))
        final_model.write(output_file=os.path.join(output_dir, "combined_models", "model_" + str(model_label) + ".pdb"), force=True)
        if verbose: 
            print("Model finished in: " + os.path.join(output_dir, "combined_models", "model_" + str(model_label) + ".pdb"))


       
def build_dna_structure(dna_seq, dna_str, output_dir, dummy_dir, beggining=None, ending=None, verbose=False):

    x3dna_path = config.get("Paths", "x3dna_path")
    os.environ['X3DNA'] = x3dna_path[:-4]
    print("build_dna_structure in tf_complex")
    print("dna_seq: " + dna_seq + "  beggining: " + str(beggining) + "  ending: " + str(ending))
    if (beggining != None) and (ending != None):
        dna_dummy = os.path.join(dummy_dir, "dna_" + str(beggining) + "-" + str(ending) + ".pdb")
        dna_file = os.path.join(output_dir, "dna_" + str(beggining) + "-" + str(ending) + ".pdb")
        modeling_seq = dna_seq[int(beggining)-1:int(ending)]
    else:
        dna_dummy = os.path.join(dummy_dir, "dna_" + str(os.getpid()) + ".pdb")
        dna_file = os.path.join(output_dir, "dna_" + str(os.getpid()) + ".pdb")
        modeling_seq = dna_seq

    print("modeling_seq: " + str(modeling_seq))

    if dna_str == "B_bent":
        """
        # Check the length of the sequence #
        if len(modeling_seq) > 146:
            modeling_seq = modeling_seq[:146]
        # Build the parameters file #
        original_parameters_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "nucleosomal_x3dna", "reference_parameters_147nuc.par"))
        opf = open(original_parameters_file, "r").readlines()
        new_parameters_file = os.path.join(dummy_dir, "nucleosomal_dna_parameters.par")
        npf = open(new_parameters_file, "w")
        # Write first line including sequence length #
        npf.write(" " + str(len(modeling_seq)) + " # base-pairs\n")
        # Write two lines like in the original file #
        for line in opf[1:3]:
            npf.write(line.rstrip() + "\n")
        # Write a line for each basepair #
        complementary = {"A": "T", "C": "G", "G": "C", "T": "A"}
        for i in range(0, len(modeling_seq)):
            line = opf[4+i].rstrip()
            nuc = modeling_seq[i]
            comp_nuc = complementary[nuc]
            npf.write(nuc + "-" + comp_nuc + line[3:] + "\n")
        npf.close()
        print("\n\nParameters file created in: " + new_parameters_file + "\n\n")
        # Use the new parameters file to create a PDB structure #
        print(os.path.join(x3dna_path, "rebuild") + " -atomic " + new_parameters_file + " " + dna_file)
        process = subprocess.check_output([os.path.join(x3dna_path, "rebuild"), "-atomic", new_parameters_file, dna_file])
        """
        # Check the length of the sequence #
        if len(modeling_seq) > 146:
            modeling_seq = modeling_seq[:146]
        # Open the reference structure of nucleosomal DNA #
        files_dir = config.get("Paths", "files_path")
        nuc_dna = PDB(os.path.join(files_dir, "nucleosomal_dna.pdb"))
        # Generate a new pdb with the nucleotide length that we are interested in #
        # This structure has two chains: I and J #
        new_nuc_dna = PDB()
        chain_1 = ChainOfNucleotide("dna", "A")
        chain_2 = ChainOfNucleotide("dna", "B")
        for i in range(0, len(modeling_seq)):
            res_for = nuc_dna.get_chain_by_id("I").nucleotides[i]
            #res_for.number = i + 1
            res_rev = nuc_dna.get_chain_by_id("J").nucleotides[-(i+1)]
            #res_rev.number = len(modeling_seq) - i
            chain_1.add_residue(res_for)
            chain_2.add_residue(res_rev)
        new_nuc_dna.add_chain(chain_1)
        new_nuc_dna.add_chain(chain_2)
        new_nuc_dna.write(output_file=os.path.join(dummy_dir, "nucleosomal_dna.pdb"), force=True)
        if verbose: print("\n\nnew dna structure created at: " + os.path.join(dummy_dir, "nucleosomal_dna.pdb") + "\n\n")
        # Now modify the sequence of this pdb structure with mutate_bases #
        start_1 = chain_1.first_nucleotide.number
        end_1 = chain_1.last_nucleotide.number
        start_2 = chain_2.first_nucleotide.number
        end_2 = chain_2.last_nucleotide.number
        # Create the mutations file #
        complementary = {"A": "T", "C": "G", "G": "C", "T": "A"}
        mut_file = open(os.path.join(dummy_dir, str(os.getpid()) + ".mutations.dat"), "w")
        for i in range(0, len(modeling_seq)):
            mut_file.write("chain=%s seqnum=%s mutation=%s\n" % ("A", str(start_1+i), modeling_seq[i]))
            mut_file.write("chain=%s seqnum=%s mutation=%s\n" % ("B", str(start_2-i), complementary[modeling_seq[i]]))
        mut_file.close()
        if verbose: print("\n\nMutation file created in: " + os.path.join(dummy_dir, str(os.getpid()) + ".mutations.dat") + "\n\n")
        # Execute mutate_bases #    
        if verbose: print(os.path.join(x3dna_path, "mutate_bases") + " -l " + os.path.join(dummy_dir, str(os.getpid()) + ".mutations.dat") + " " + os.path.join(dummy_dir, "nucleosomal_dna.pdb") + " " + os.path.join(dummy_dir, "nucleosomal_dna_modified.pdb"))
        process = subprocess.check_output([os.path.join(x3dna_path, "mutate_bases"), "-l", os.path.join(dummy_dir, str(os.getpid()) + ".mutations.dat"), os.path.join(dummy_dir, "nucleosomal_dna.pdb"), os.path.join(dummy_dir, "nucleosomal_dna_modified.pdb")])
        if verbose: print("\n\nnew dna PDB created in: " + os.path.join(dummy_dir, "nucleosomal_dna_modified.pdb") + "\n\n")
        return os.path.join(dummy_dir, "nucleosomal_dna_modified.pdb")



    elif dna_str == "B":
        if verbose: print(os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -b " + dna_file)
        #process = subprocess.check_output([os.path.join(x3dna_path, "fiber"), "-seq=" + modeling_seq, "-b", dna_dummy], stderr=subprocess.STDOUT, env=os.environ)
        os.system(os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -b " + dna_file)
        #os.system("cp " + dna_dummy + " " + dna_file)
        return dna_file

    elif dna_str == "A":
        if verbose: print(os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -a " + dna_file)
        process = subprocess.check_output([os.path.join(x3dna_path, "fiber"), "-seq=" + modeling_seq, "-a", dna_dummy], stderr=subprocess.STDOUT, env=os.environ)
        return dna_dummy
    
    elif dna_str == "C":
        if verbose: print(os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -c " + dna_file)
        process = subprocess.check_output([os.path.join(x3dna_path, "fiber"), "-seq=" + modeling_seq, "-c", dna_dummy], stderr=subprocess.STDOUT, env=os.environ)
        return dna_dummy

    elif dna_str == "D":
        if verbose: print(os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -d " + dna_file)
        process = subprocess.check_output([os.path.join(x3dna_path, "fiber"), "-seq=" + modeling_seq, "-d", dna_dummy], stderr=subprocess.STDOUT, env=os.environ)
        return dna_dummy

    elif dna_str == "Z":
        if verbose: print(os.path.join(x3dna_path, "fiber") + " -seq=" + modeling_seq + " -z " + dna_file)
        process = subprocess.check_output([os.path.join(x3dna_path, "fiber"), "-seq=" + modeling_seq, "-z", dna_file], stderr=subprocess.STDOUT, env=os.environ)

    return dna_file


def get_binding_site_coordinates(dna_seq, coop, coords_file=None, get_coords_from_file=False, verbose=False):


    if get_coords_from_file:
        # Parse the files and get the coordinates from the file name #
        coords = []
        list_of_models = coop.split(";")
        for i in range(0, len(list_of_models)):
            # Get the DNA structure that will be superimposed #
            model = list_of_models[i]
            coord_1 = os.path.basename(model).split(".")[2].split(":")[0].split("-")[0]
            coord_2 = os.path.basename(model).split(".")[2].split(":")[0].split("-")[1]
            coords.append(coord_1 + ":" + coord_2)

        return coords

    if os.path.isfile(str(coords_file)):
        # Parse the coordinates file #
        coords = []
        for line in open(coords_file, "r").readlines():
            coords.append(line.rstrip())

        return coords 

    else:
        # We have to get the coordinates from the models #
        model_list = coop.split(";")
        coords = []
        for mod in model_list:
            pdb_obj = PDB(mod)
            for chain in pdb_obj.chains():
                if isinstance(chain, ChainOfNucleotide):
                    nuc_seq = chain.nucleotide_sequence()
                    # Now, fit the DNA sequence of this model in the global DNA sequence we are working with #
                    start = dna_seq.find(nuc_seq)
                    end = start + len(nuc_seq)
                    coords.append(str(start) + ":" + str(end))
                    break

        return coords


def get_individual_models(input_models, verbose):

    coop_list = []
    for line in open(os.path.abspath(input_models), "r").readlines():
        if os.path.isfile(line.rstrip()):
            coop_list.append(line.rstrip())
    coop = ";".join(coop_list)

    return coop


def get_threading_models(input_threading, output_dir, pdb_dir, dummy_dir, verbose):

    # For file in input_threading execute model_protein with a threading input #
    for line in open(os.path.abspath(input_threading), "r").readlines():
        if os.path.isfile(line.rstrip()):
            if verbose: print(python + " " + scripts_path + "/model_protein.py -i " + line.rstrip() + " -o " + os.path.join(output_dir, "binary_interactions") + " --p=" + pdb_dir + " --dummy=" + dummy_dir + " -t -d")
            os.system(python + " " + scripts_path + "/model_protein.py -i " + line.rstrip() + " -o " + os.path.join(output_dir, "binary_interactions") + " --p=" + pdb_dir + " --dummy=" + dummy_dir + " -t -d")

    model_list = []
    for file in os.listdir(os.path.join(output_dir, "binary_interactions")):
        if file.endswith(".pdb"):
            if not file.startswith("dna"):
                model_list.append(os.path.join(output_dir, "binary_interactions", file))
    coop = ";".join(model_list)

    return coop


def model_dna_structure(dna_structure, dna_seq, verbose):

    # Load PDB #
    dna_pdb = PDB(dna_structure)
    forward_chain = None
    for chain in dna_pdb.chains:
        if isinstance(chain, ChainOfNucleotide):
            forward_chain = chain
            break

    # Check if the length of the DNA structure and sequence are the same #
    if len(forward_chain.nucleotide_sequence()) != len(dna_seq):
        if verbose: print("The input DNA structure doesn't fits with the input DNA sequence.\n\n")
        exit(0)

    # Load an x3dna object and initialize #
    mutations_file = os.path.join(dummy_dir, str(os.getpid()) + ".mutations.dat")
    model_pdb_file = os.path.join(dummy_dir, str(os.getpid()) + "mutated.pdb")
    x3dna_obj = x3dna.get_x3dna_obj(dna_structure)

    # Execute mutate bases #
    nucleotides = list(dna_seq)
    for bp in x3dna_obj.get_basepairs():
        nucleotide = nucleotides.pop(0).upper()
        ((pdb_chain_A, residue_num_A), (pdb_chain_B, residue_num_B)) = x3dna_obj.get_basepair(bp)
        # Mutate the forward strand #
        if nucleotide != "N":
            functions.write(mutations_file, "chain=%s seqnum=%s mutation=%s" % (pdb_chain_A, str(residue_num_A), nucleotide))
            
        # Mutate the reverse strand #
        if nucleotide != "N":
            functions.write(mutations_file, "chain=%s seqnum=%s mutation=%s" % (pdb_chain_B, str(residue_num_B), dna_complementary[nucleotide]))
            
    try:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        x3dna_path = os.path.join(src_path, config.get("Paths", "x3dna_path"))
        os.environ['X3DNA'] = x3dna_path[:-4]
        # Exec process #
        process = subprocess.check_output([os.path.join(x3dna_path, "mutate_bases"), "-l", mutations_file, dna_structure, model_pdb_file], stderr=subprocess.STDOUT, env=os.environ)
        if verbose: print("Created file in: " + model_pdb_file)
    except:
        raise ValueError("Could not exec X3DNA for %s" % dna_structure)


    return model_pdb_file


def fill_dna_gaps(dna_seq, coop, coords):

    print("\nmodels")
    for c in coop.split(";"):
        print(c)

    print("\ncoordinates")
    for c in coords:
        print(c)




#-------------#
# Options     #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("Usage: dna_superimposer.py -a pdb_static -b pdb_mobile [--dummy=dummy_dir -i interface_alignment -o output_dir -v]")

    parser.add_option("-m", action="store", type="string", dest="input_models", help="Input file containing the paths to the TF-DNA models that should be merged in a single TF-DNA complex", metavar="{finelame}")
    parser.add_option("-t", action="store", type="string", dest="input_threading", default=None, help="Input file containing the paths to the threading files that will be used to create a TF-DNA complex", metavar="{finelame}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="{directory}")
    parser.add_option("-l", "--label", default="", action="store", type="string", dest="label", help="Label for the output file", metavar="{directory}")
    parser.add_option("-d", default=None, action="store", type="string", dest="dna", help="Reference DNA sequence for building the model", metavar="{directory}")
    parser.add_option("--dna_str", default=None, action="store", type="string", dest="dna_structure", help="PDB structure to be used in the TF-DNA superimposition", metavar="{directory}")
    parser.add_option("--dna_conf", default="B", action="store", type="string", dest="dna_conformation", help="Conformation that the DNA should have in the final model (either B, A or B_bent)", metavar="{directory}")
    parser.add_option("-c", "--coord", default=None, action="store", type="string", dest="coords_file", help="File containing the binding site coordinates for each of the input models", metavar="{directory}")
    parser.add_option("-p", "--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (i.e. output dir from pdb.py)", metavar="{directory}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")


    (options, args) = parser.parse_args()

    if (options.input_models is None) and (options.input_threading is None):
        parser.error("missing arguments: type option \"-h\" for help")

    return options


if __name__ == "__main__":

    # Options and arguments #
    options = parse_options()
    output_dir = os.path.abspath(options.output_dir)
    pdb_dir = os.path.abspath(options.pdb_dir)
    dummy_dir = os.path.abspath(options.dummy_dir)
    verbose= options.verbose
    
    # Handle directories #
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if not os.path.exists(os.path.join(output_dir, "binary_interactions")):
        os.mkdir(os.path.join(output_dir, "binary_interactions"))
    else:
        os.system("rm -rf " + os.path.join(output_dir, "binary_interactions", "*"))

    if not os.path.exists(os.path.join(output_dir, "combined_models")):
        os.mkdir(os.path.join(output_dir, "combined_models"))
    else:
        os.system("rm -rf " + os.path.join(output_dir, "combined_models", "*"))
    
    if not os.path.exists(os.path.join(output_dir, "working")):
        os.mkdir(os.path.join(output_dir, "working"))
    else:
        os.system("rm -rf " + os.path.join(output_dir, "working", "*"))
    
    coop_dir = os.path.join(output_dir, "combined_models")

    # Parse the DNA sequence #
    dna_seq = ""
    if options.dna != None:
        with open(os.path.abspath(options.dna), 'r') as ff:
            for line in ff.readlines()[1:]:
                dna_seq += str(line.rstrip())
        ff.close()


    print("dna_structure: " + str(options.dna_structure))

    # Create the DNA structure #
    if options.dna_structure == None:
        dna_pdb = build_dna_structure(dna_seq, options.dna_conformation, output_dir, dummy_dir, verbose=verbose)

    else:
        if options.dna != None:
            dna_pdb = model_dna_structure(os.path.abspath(options.dna_structure), dna_seq, verbose)
            
        else:
            dna_pdb = os.path.abspath(options.dna_structure)

    # Parse the input file and get the input structures #
    if options.input_threading != None:
        # Use threading files #
        coop = get_threading_models(options.input_threading, output_dir, pdb_dir, dummy_dir, verbose)
    else:
        # Use structural models #
        coop = get_individual_models(options.input_models, verbose)
        
    # Get binding site coordinates #
    if options.input_threading != None:
        coords = get_binding_site_coordinates(dna_seq, coop, options.coords_file, get_coords_from_file=True, verbose=verbose)
    else:
        coords = get_binding_site_coordinates(dna_seq, coop, options.coords_file, verbose=verbose)

    # Correct the models to do not have DNA gaps #
    fill_dna_gaps(dna_seq, coop, coords)

    
    # Build the cooperative model #
    #get_cooperative_model(coop, coords, coop_dir, dna_seq, dna_pdb, options.dna_conformation, dummy_dir, verbose=verbose)

    