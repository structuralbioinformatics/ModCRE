import os, sys, re
from Bio import motifs as mm
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import ConfigParser
import gzip
import numpy
import optparse
import shutil
import socket
import subprocess

# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.dirname(__file__))

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Imports my functions #
import functions

# Imports jbonet's module #
from SBI.data import aminoacids1to3, aminoacids_polarity_boolean, nitrogenous_bases
from SBI.structure import PDB

# Import my modules #
import blast, contacts, dssp, hmmer, interface, model_protein, nr, triads, x3dna, homologs, threader
import pwm_pbm as PWM
def parse_options():
    '''
    Create threading files to perform modeling for a set of homologs using the aligned residues of the binding site with a PDB template.
    '''

    parser = optparse.OptionParser("thread.py -i INPUT_FILE  --pdb=PDB_FILE --chain=CHAIN --dna=DNA_SEQ [--filter --dummy DUMMY_DIR --out OUTPUT_DIRECTORY --specie SPECIE -v]")

    parser.add_option("-i", action="store", type="string", dest="input_file", default=None, help="Input PDB file", metavar="INPUT_FILE")
    parser.add_option("--chain", action="store", default=None, type="string", dest="chain", help="Selected chain of PDB", metavar="CHAIN")
    parser.add_option("--dna", action="store", default=None, type="string", dest="dna_seq", help="DNA sequence to be thread", metavar="DNA_SEQ")
    parser.add_option("--dummy", action="store", type="string", dest="dummy_dir", default="/tmp", help="Dummy directory (default is /tmp)", metavar="DUMMY_DIR")
    parser.add_option("-o","--out", action="store", default="aux_files", type="string", dest="output_dir", help="Output directory (default is 'aux_files')", metavar="OUTPUT_FILE")
    parser.add_option("--specie", action="store", type="string", dest="specie",  default=None, help="Specie to obtain specific orthologs (i.e. taxon/code/common_name as 9606/HUMAN/'Homo sapiens')", metavar="SPECIE")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")
    parser.add_option("-d", default="basepairs", action="store", type="string", dest="distance_type", help="Distance type (i.e. \"basepairs\", \"dinucleotides\" or \"mindist\"; default = dinucleotides)", metavar="{string}")

    (options, args) = parser.parse_args()

    if options.input_file is None  or options.pdb_file is None  or options.chain is None or options.dna_seq is None :
        parser.error("missing arguments: type option \"-h\" for help")

    return options



#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options   = parse_options()

    # Initialize
    output_dir=options.output_dir
    dummy_dir=options.dummy_dir
    if not os.path.exists(output_dir): os.mkdir(output_dir)
    if not os.path.exists(dummy_dir):  os.mkdir(dummy_dir)
    pdb_chain = options.chain
 
           if options.verbose: sys.stdout.write("\t\t-- threading TF sequences...\n")
            # For each sequence... #
            for i in range(len(tf_obj._sequences)):
                # Skip if sequence file does not exists #
                sequence_file = os.path.join(options.output_dir, "sequences", tf_obj.get_id() + "." + str(i) + ".fa")
                if not os.path.exists(sequence_file): continue
                # Initialize #
                hits = []
                hits_file = os.path.join(options.output_dir, "hits", tf_obj.get_id() + "." + str(i) + ".txt")
                if os.path.exists(hits_file):
                    # For each line... #
                    for line in functions.parse_file(hits_file):
                        hits.append(line.split("\t"))
                else:
                    # Blast TF sequence #
                    blast_obj = blast.get_blast_obj(os.path.join(options.pdb_dir, "database", "database.fasta"), sequence_file)
                    # Get unfiltered hits #
                    hits = model_protein.filter_blast_hits_by_interface(os.path.abspath(sequence_file), blast_obj, os.path.abspath(options.pdb_dir), filter_twilight_zone_hits=True)
                # Skip if not enough hits #
                if len(hits) == 0: continue
                # Initialize #
                pwm = []
                dummy_file = os.path.join(options.dummy_dir, "pwm.txt")
                pwm_file = os.path.join(os.path.join(options.pwm_dir, tf_obj._motifs[i] + ".txt"))
                sys.stdout.write("\t\t-- Using PWM %s\n"%pwm_file)
                # Remove hits file if exists #
                if os.path.exists(hits_file): os.remove(hits_file)
                # For each line... #
                for line in functions.parse_file(pwm_file):
                    line = line.split("\t")
                    position = line.pop(0)
                    if position != "Pos":
                        pwm.append(map(str, [float(j) * 100000 for j in line]))
                # For each hit... #
                for hit in hits:
                    functions.write(hits_file, "\t".join(map(str, hit)))
                    # Skip if threading file already exists #
                    threading_file = os.path.join(options.output_dir, "threading", tf_obj.get_id() + "." + str(i) + "." + hit[0] + ".txt")
                    if not os.path.exists(threading_file):
                        # Initialize #
                        kmers_wrk= []
                        # Get PDB object #
                        try:
                          pdb_obj = PDB(os.path.join(options.pdb_dir, "clean", hit[0][:4] + ".pdb"))
                        except:
                          sys.stdout.write("Error: cannot read %s\n"%hit[0][:4])
                          continue
                        # Get X3DNA object #
                        try:
                          x3dna_obj = x3dna.X3DNA(os.path.join(options.pdb_dir, "x3dna", hit[0][:4] + ".txt"))
                        except:
                          sys.stdout.write("Error: cannot read DNA %s\n"%hit[0][:4])
                          continue
                        # Get interface object #
                        interface_obj = interface.Interface(os.path.join(options.pdb_dir, "interfaces", hit[0] + ".txt"))
                        # Get sequence #
                        try:
                          sequence = ""
                          for basepair in interface_obj.get_interface_basepairs():
                            for pdb_chain, residue_num in x3dna_obj.get_basepair(basepair):
                                sequence += pdb_obj.get_chain_by_id(pdb_chain).get_residue_by_identifier(str(residue_num)).single_letter
                                break
                        except:
                          sys.stdout.write("Error: cannot get interface sequence %s\n"%(hit[0] + ".txt"))
                          continue
                        # Format sequence #
                        sequence = Seq(sequence, IUPAC.unambiguous_dna)
                        # For each row... #
                        for j in zip(*pwm):
                            functions.write(dummy_file, "\t".join(j))
                        # Get motif #
                        f = open(dummy_file)
                        motif = mm.read(f, "pfm")
                        f.close()
                        # Remove dummy file #
                        #os.remove(dummy_file)
                        # Get pseudocounts #
                        motif.pseudocounts = mm.jaspar.calculate_pseudocounts(motif)
                        # Get motif max. score, min. score and score threshold #
                        max_score = motif.pssm.max
                        min_score = motif.pssm.min
                        # Get the 80% best matches
                        motif_score_threshold = (max_score - min_score) * 0.8 + min_score
                        # For each matching position, score... #
                        for position, score in motif.pssm.search(sequence, threshold=motif_score_threshold):
                            strand = "+"
                            # If number is negative is reverse complement:
                            if position < 0:
                                # Real position = position + len(sequence) #
                                position += len(sequence)
                                strand = "-"
                            # For each line... #
                            for line in functions.parse_file(os.path.join(options.output_dir, "kmers", tf_obj.get_id() + "." + str(i) + ".txt")):
                                if line.startswith("#"): continue
                                line = line.split(";")
                                if strand == "+":
                                    kmers_wrk.append([line[0], position + int(line[1])])
                                else:
                                    kmers_wrk.append([triads.get_complementary_dna_sequence(line[0]), position + len(motif) - 8 - int(line[1])])
                            break
                        # Skip if k-mers could not be mapped #
                        if len(kmers_wrk) == 0: continue
                        # Get identities #
                        identities = sum(map(int, numpy.frombuffer(hit[1], dtype=numpy.byte) == numpy.frombuffer(hit[2], dtype=numpy.byte)))
                        threading_protein = {}
                        # Get PDB file #
                        protein_chain_pdb_obj = PDB(os.path.join(options.pdb_dir, "split", hit[0] + ".pdb"))
                        # Get positions #
                        first = protein_chain_pdb_obj.chains[0].gapped_protein_sequence.find(hit[2].replace("-", ""))
                        positions = range(first, len(hit[2].replace("-", "")) + first)
                        # Get sequence/PDB correlation #
                        sequence_to_crystal, crystal_to_sequence = model_protein.get_sequence_to_crystal_correlations(protein_chain_pdb_obj, protein_chain_pdb_obj.chains[0].chain)
                        # For each aligned position... #
                        for position in range(len(hit[2])):
                            # Skip if hit alignment position is gapped #
                            if hit[2][position] != "-":
                                sequence_position = positions.pop(0)
                                if sequence_position in sequence_to_crystal:
                                    threading_protein.setdefault( (hit[0][-1], sequence_to_crystal[sequence_position]), hit[1][position] )
                        threading_pdb_name   =hit[0][:4]
                        threading_pdb_chain  =hit[0][-1]
                        threading_dna_helix  =x3dna_obj.get_basepair_helix(interface_obj.get_interface_start())
                        threading_identity   =(float(identities) / (len(hit[2]) - hit[2].count("-")))
                        threading_query_ali  =hit[1]
                        threading_hit_ali    =hit[2]
                        threading_coverage   =100* float(min( (len(hit[2]) - hit[2].count("-")), (len(hit[1]) - hit[1].count("-"))) ) / float(max( (len(hit[2]) - hit[2].count("-")), (len(hit[1]) - hit[1].count("-"))) )
                        threading_dna        ={}
                        threading_dnae       ={}
                        for kmer in kmers_wrk:
                            threading_dna.setdefault(str(kmer[0]),str(kmer[1]))
                            threading_dnae.setdefault(str(kmer[0]),str(kmer[1]))
                        threaded_obj = threader.Threaded(None)
                        threaded_obj.set_pdb_name(threading_pdb_name)
                        threaded_obj.set_pdb_chain(threading_pdb_chain)
                        threaded_obj.set_dna_helix(threading_dna_helix)
                        threaded_obj.set_identity(threading_identity)
                        threaded_obj.set_coverage(threading_coverage)
                        threaded_obj.set_protein(threading_protein)
                        threaded_obj.set_dna(threading_dna)
                        threaded_obj.set_dna_fixed(threading_dnae)
                        threaded_obj.set_query_ali(threading_query_ali)
                        threaded_obj.set_pdb_ali(threading_hit_ali)
                        # Create threaded file #
                        threaded_obj.write(threading_file)


