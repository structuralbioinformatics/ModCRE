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

#-------------#
# Classes     #
#-------------#

class TF(object):
    """
    This class defines an {TF} object.

    """

    def __init__(self, id=None, name=None, specie=None, family=None, motifs=None, sources=None, sequences=None, file_name=None):
        self._file = file_name
        self._id = id
        self._name = name
        self._specie = specie
        self._family = family
        self._motifs = motifs
        self._sources = sources
        self._sequences = sequences

        if self._file is not None:
            self._parse_file()

    def _parse_file(self):
        for line in functions.parse_file(self._file):
            if line.startswith("#"): continue
            line = line.split(";")
            self._id = line[0]
            self._name = line[1]
            self._specie = line[2]
            self._family = line[3].split(",")
            self._motifs = line[4].split(",")
            self._sources = line[5].split(",")
            self._sequences = line[6].split(",")

    def get_id(self):
        return self._id

    def get_name(self):
        return self._name

    def get_family(self):
        return self._family[0]

    def get_sequence(self):
        return self._sequences[0]

    def get_motif(self):
        return self._motifs[0]

    def get_families(self):
        return self._family

    def get_sequences(self):
        return self._sequences

    def get_motifs(self):
        return self._motifs

    def write(self, file_name):
        # Initialize #
        if os.path.exists(file_name): os.remove(file_name)
        functions.write(file_name, "#id;name;specie;family;motifs;sources;sequences")
        functions.write(file_name, "%s;%s;%s;%s;%s;%s;%s" % (self._id, self._name, self._specie, ",".join(self._family), ",".join(self._motifs), ",".join(self._sources), ",".join(self._sequences)))

#-------------#
# Options     #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python pbm.py -e escores_file -f families_file -m motifs_file --pdb=pdb_dir --pwm=pwm_dir -s sources_file -t tfs_file  -u uniprot_file [-c --dummy dummy_dir -o output_dir --start=start_step --stop=stop_step -v]")

    # parser.add_option("-c", "--cluster", default=False, action="store_true", dest="cluster_mode", help="Cluster mode (this submits jobs to queues; default = False)")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-e", action="store", type="string", dest="escores_file", help="E-scores file (from CIS-BP; i.e. Escores.txt)", metavar="{filename}")
    parser.add_option("-f", action="store", type="string", dest="families_file", help="Families file (from CIS-BP; i.e. cisbp_1.02.tf_families.sql)", metavar="{filename}")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="{directory}")
    parser.add_option("-m", action="store", type="string", dest="motifs_file", help="Motifs file (from CIS-BP; i.e. cisbp_1.02.motifs.sql)", metavar="{filename}")
    parser.add_option("--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (output directory from pdb.py)", metavar="{directory}")
    parser.add_option("--pwm", action="store", type="string", dest="pwm_dir", help="PWM directory (from CIS-BP)", metavar="{directory}")
    parser.add_option("-s", action="store", type="string", dest="sources_file", help="Sources file (from CIS-BP; i.e. cisbp_1.02.motif_sources.sql)", metavar="{filename}")    
    parser.add_option("--start", default=1, action="store", type="int", dest="start_step", help="Start at a given step (default = 1; first)", metavar="{int}")
    parser.add_option("--stop", default=11, action="store", type="int", dest="stop_step", help="Stop at a given step (default = 11; last)", metavar="{int}")
    parser.add_option("-t", action="store", type="string", dest="tfs_file", help="TFs file (from CIS-BP; i.e. cisbp_1.02.tfs.sql)", metavar="{filename}")
    parser.add_option("-u", action="store", type="string", dest="uniprot_file", help="UniProt file (i.e. uniprot_sprot.fasta + uniprot_trembl.fasta", metavar="{filename}")
    parser.add_option("-j","--parallel", default=False, action="store_true", dest="parallel", help="Submit JOBS to Queues in parallel. The program stops, you need to re-start with the requested value (default = False)", metavar="{boolean}")
    parser.add_option("--radius",default=0, action="store", type="float", dest="radius", help="Maximum contact distance to calculate interactions to create PWMs (default=0 implies the use of 'max_contact_distance' from configuration", metavar="{string}")
  

    # parser.add_option("-t", action="store", type="string", dest="taxonomy_file", help="Taxonomy file (i.e. speclist.txt)", metavar="{filename}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")

    (options, args) = parser.parse_args()

    # if options.pdb_dir is None or options.tfs_file is None:
    #     parser.error("missing arguments: type option \"-h\" for help")

    return options

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Get hostname #
    hostname = socket.gethostname()

    # Create output "main" subdirs #
    for subdir in ["motifs", "tfs", "kmers", "sequences", "hits", "threading", "triads", "msa", "hmm", "nr", "potentials", "uniprot", "homologs","pwms"]:
        if not os.path.exists(os.path.join(options.output_dir, subdir)):
            os.makedirs(os.path.join(options.output_dir, subdir))
    if not os.path.exists(os.path.abspath(options.dummy_dir)):
            os.makedirs(os.path.abspath(options.dummy_dir))

    # Verbose mode #
    if options.verbose: sys.stdout.write("\n")

    ###############################
    # 1. Parse CIS-BP database    #
    ###############################
    if options.verbose: sys.stdout.write("Parse files from CIS-BP database...\n\n")
    # Skip if starts later #
    if options.start_step <= 1:
       if options.parallel:
        if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
        else: cluster_queue=config.get("Cluster", "cluster_queue")
        functions.submit_command_to_queue("%s %s --pdb=%s -o %s --dummy=%s -t %s -e %s -f %s -m %s --pwm %s -s %s -u %s --start=1 --stop=1 -v" % ( os.path.join(config.get("Paths", "python_path"), "python"),  os.path.join(scripts_path,os.path.basename(__file__)) , os.path.abspath(options.pdb_dir), os.path.abspath(options.output_dir) , os.path.abspath(options.dummy_dir), os.path.abspath(options.tfs_file), os.path.abspath(options.escores_file), os.path.abspath(options.families_file), os.path.abspath(options.motifs_file), os.path.abspath(options.pwm_dir), os.path.abspath(options.sources_file), os.path.abspath(options.uniprot_file)),cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
       else:
        ##############################
        # 1.1 Parse E-scores file    #
        ##############################
        if options.verbose: sys.stdout.write("\t\t-- parsing e-scores...\n")
        # Initialize #
        motif_names = {}
        # For each line... #
        escores=open(os.path.abspath(options.escores_file),"r")
        for linee in escores:
            line = linee.strip().split("\t")
            for i in xrange(len(line)):
                motif_names.setdefault(line[i][:11],[]).append(i)
            break
        escores.close()
        for motif_name,columns in motif_names.iteritems():
            escores=open(os.path.abspath(options.escores_file),"r")
            motif_file = os.path.join(options.output_dir, "motifs", motif_name + ".txt")
            if os.path.exists(motif_file): 
               if options.verbose: sys.stdout.write("\t\t\t-- use file %s\n"%(motif_name+ ".txt"))
               continue
            if options.verbose: sys.stdout.write("\t\t\t-- create file %s\n"%(motif_name+ ".txt"))
            functions.write(motif_file, "#k-mer,complementary;e-score")
            skip=True
            for line in escores:
                line = line.strip().split("\t")
                if skip:
                   skip=False
                   continue
                kmer=line[0]
                values=[]
                for j in columns:
                  try: values.append('{0:.4g}'.format(float(line[j])))
                  except: continue
                try:    complementary=triads.get_complementary_dna_sequence(kmer)
                except: continue
                if len(values) > 0:
                   functions.write(motif_file, "%s;%s;%s" % (kmer, triads.get_complementary_dna_sequence(kmer), '{0:.3g}'.format(numpy.mean(map(float, values)))))
                else:
                   functions.write(motif_file, "%s;%s;None" % (kmer, triads.get_complementary_dna_sequence(kmer)))
            escores.close()
        ##############################
        # 1.2 Parse motif sources    #
        ##############################
        if options.verbose: sys.stdout.write("\t\t-- parsing motif sources...\n")
        # Initialize #
        sources = {}
        # For each line... #
        for line in functions.parse_file(os.path.abspath(options.sources_file)):
            #m = re.search("\('(.+)', '(.+)', '.+', '.+', \d+, '.+', '.+'\),*", line)
            m = re.search("\('(.+)', '(.+)', '.+', '.+', '.+', '.+', '.+'\),*", line)
            if m:
                sources.setdefault(m.group(1), m.group(2))
        ##############################
        # 1.3 Parse motifs info      #
        ##############################
        if options.verbose: sys.stdout.write("\t\t-- parsing motif information...\n")
        # Initialize #
        motifs = {}
        # For each line... #
        for line in functions.parse_file(os.path.abspath(options.motifs_file)):
            m = re.search("\('(.+)', '(.+)', '(.+)', '.+', 'PBM', '(.+)', '.+', '.+'\),*", line)
            if m:
                # Skip if motif file does not exist #
                motif_name=m.group(1)
                if not os.path.exists(os.path.join(options.output_dir, "motifs", motif_name + ".txt")): continue
                # Skip if not characterized protein sequence #
                if m.group(4) == "NULL": continue
                motifs.setdefault(m.group(1), [m.group(2), m.group(3), re.sub("[^A-Z]", "X", m.group(4).upper())])
        ##############################
        # 1.4 Parse TFs families     #
        ##############################
        if options.verbose: sys.stdout.write("\t\t-- parsing TF families...\n")
        # Initialize #
        families = {}
        # For each line... #
        for line in functions.parse_file(os.path.abspath(options.families_file)):
            #m = re.search("\('(.+)', '(.+)', '.+', \d+, .+\),*", line)
            m = re.search("\('(.+)', '(.+)', '.+', '.+', '.+'\),*", line)
            if m:
                families.setdefault(m.group(1), set(m.group(2).split(",")))
        ##############################
        # 1.5 Parse TFs info         #
        ##############################
        if options.verbose: sys.stdout.write("\t\t-- parsing TF information...\n")
        # For each line... #
        for line in functions.parse_file(os.path.abspath(options.tfs_file)):
            m = re.search("\('(.+)', '(.+)', '.+', '.+', '(.+)', '(.+)', '[DIN]'\),*", line)
            if m:
                # Initialize #
                tf_motifs = []
                tf_sources = []
                tf_sequences = []
                # For each motif... #
                for motif in motifs:
                    if motifs[motif][0] == m.group(1):
                        tf_motifs.append(motif)
                        tf_sources.append(sources[motifs[motif][1]])
                        tf_sequences.append(motifs[motif][2])
                if len(tf_motifs) > 0:
                    tf_obj = TF(m.group(1), m.group(3), re.sub("_", " ", m.group(4)), families[m.group(2)], tf_motifs, tf_sources, tf_sequences)
                    # Skip if TF file already exists #
                    tf_file = os.path.join(options.output_dir, "tfs", m.group(1) + ".txt")
                    if not os.path.exists(tf_file):
                        if options.verbose: sys.stdout.write("\t\t\t-- create file %s\n"%( m.group(1)+ ".txt"))
                        tf_obj.write(tf_file)
                    else:
                        if options.verbose: sys.stdout.write("\t\t\t-- use file %s\n"%( m.group(1)+ ".txt"))

    # Exit if stops here #
    if options.parallel and options.start_step<2:
        sys.stdout.write("Exiting: wait until all submitted runs have finnished and restart again with start=2\n")
        exit(0)

    # Exit if stops here #
    if options.stop_step == 1:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ##################################
    # 2. Create associated PBM files #
    ##################################
    if options.verbose: sys.stdout.write("Create associated PBM files...\n\n")
    # Skip if starts later #
    if options.start_step <= 2:
       if options.parallel:
        if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
        else: cluster_queue=config.get("Cluster", "cluster_queue")
        functions.submit_command_to_queue("%s %s --pdb=%s -o %s --dummy=%s -t %s -e %s -f %s -m %s --pwm %s -s %s -u %s --start=2 --stop=2 -v" % ( os.path.join(config.get("Paths", "python_path"), "python"),  os.path.join(scripts_path,os.path.basename(__file__)) , os.path.abspath(options.pdb_dir), os.path.abspath(options.output_dir) , os.path.abspath(options.dummy_dir), os.path.abspath(options.tfs_file), os.path.abspath(options.escores_file), os.path.abspath(options.families_file), os.path.abspath(options.motifs_file), os.path.abspath(options.pwm_dir), os.path.abspath(options.sources_file), os.path.abspath(options.uniprot_file)),cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
       else:
        # For each TF... #
        for tf_file in os.listdir(os.path.join(options.output_dir, "tfs")):
            # Get TF object #
            tf_obj = TF(file_name=os.path.join(options.output_dir, "tfs", tf_file))
            # Verbose mode... #
            if options.verbose: sys.stdout.write("\t%s...\n" % tf_obj.get_id())

            ##############################
            # 2.1 Align positive k-mers  #
            ##############################
            if options.verbose: sys.stdout.write("\t\t-- aligning positive k-mers...\n")
            # For each motif... #
            for i in range(len(tf_obj._motifs)):
                # Skip if k-mers file already exists #
                kmers_file = os.path.join(options.output_dir, "kmers", tf_obj.get_id() + "." + str(i) + ".txt")
                if not os.path.exists(kmers_file):
                    # Initialize #
                    positive_kmers = []
                    # For each line... #
                    for line in functions.parse_file(os.path.join(options.output_dir, "motifs", tf_obj._motifs[i] + ".txt")):
                        if line.startswith("#"): continue
                        line = line.split(";")
                        try:
                            # If positive k-mer #
                            if float(line[2]) >= float(config.get("Parameters", "min_escore_positives")):
                                positive_kmers.append(line)
                        except: pass
                    # If positive k-mers... #
                    if len(positive_kmers) > 0:
                        # Initialize #
                        pwm = []
                        pwm_file = os.path.join(os.path.join(options.pwm_dir, tf_obj._motifs[i] + ".txt"))
                        sys.stdout.write("\t\t-- Using PWM %s\n"%pwm_file)
                        functions.write(kmers_file, "#kmer;position")
                        # For each line... #
                        for line in functions.parse_file(pwm_file):
                            line = line.split("\t")
                            position = line.pop(0)
                            if position != "Pos":
                                pwm.append(map(str, [float(j) * 100000 for j in line]))
                        # For each k-mer... #
                        for kmer in sorted(positive_kmers, key=lambda x: float(x[2]), reverse=True):
                            # Initalize #
                            hits = set()
                            # For each sliding window... #
                            for j in range(len(pwm) - 8 + 1):
                                # Initialize #
                                dummy_file = os.path.join(options.dummy_dir, "pwm.txt")
                                # For each row... #
                                sub_pwm = pwm[j:j + 8]
                                for k in zip(*sub_pwm):
                                    functions.write(dummy_file, "\t".join(k))
                                # Get motif #
                                f = open(dummy_file)
                                motif = mm.read(f, "pfm")
                                f.close()
                                # Remove dummy file #
                                os.remove(dummy_file)
                                # Get pseudocounts #
                                motif.pseudocounts = mm.jaspar.calculate_pseudocounts(motif)
                                # Get motif max. score, min. score and score threshold #
                                max_score = motif.pssm.max
                                min_score = motif.pssm.min
                                motif_score_threshold = (max_score - min_score) * 0.8 + min_score
                                # Format sequence #
                                sequence = Seq("".join(kmer[0]), IUPAC.unambiguous_dna)
                                # For each matching position, score... #
                                for position, score in motif.pssm.search(sequence, threshold=motif_score_threshold):
                                    if position == 0:
                                        hits.add((kmer[0], j, score))
                                    break
                                # Scan k-mers against motif #q
                                sequence = Seq("".join(kmer[1]), IUPAC.unambiguous_dna)
                                # For each matching position, score... #
                                for position, score in motif.pssm.search(sequence, threshold=motif_score_threshold):
                                    if position == 0:
                                        hits.add((kmer[1], j, score))
                                    break
                            # For each hit... #
                            for hit in sorted(hits, key=lambda x: x[2], reverse=True):
                                functions.write(kmers_file, "%s;%s" % (hit[0], hit[1]))
                                break

            ##############################
            # 2.2 Extract TF sequences   #
            ##############################
            if options.verbose: sys.stdout.write("\t\t-- extracting TF sequences...\n")
            # For each sequence... #
            for i in range(len(tf_obj._sequences)):
                # Skip if k-mers file does not exists #
                kmers_file = os.path.join(options.output_dir, "kmers", tf_obj.get_id() + "." + str(i) + ".txt")
                if not os.path.exists(kmers_file): continue
                # Skip if sequence file already exists #
                sequence_file = os.path.join(options.output_dir, "sequences", tf_obj.get_id() + "." + str(i) + ".fa")
                if not os.path.exists(sequence_file):
                    functions.write(sequence_file, ">%s\n%s" % (tf_obj.get_id(), tf_obj._sequences[i]))

            ##############################
            # 2.3 Thread sequences       #
            ##############################
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

            ##############################
            # 2.4 Get prot-DNA triads    #
            ##############################
            if options.verbose: sys.stdout.write("\t\t-- getting protein-DNA triads...\n")
            # For each sequence... #
            for i in range(len(tf_obj._sequences)):
                # Skip if sequence file does not exists #
                sequence_file = os.path.join(options.output_dir, "sequences", tf_obj.get_id() + "." + str(i) + ".fa")
                if not os.path.exists(sequence_file): continue
                # For each threading file... #
                for threading_file in os.listdir(os.path.join(options.output_dir, "threading")):
                    if tf_obj.get_id() in threading_file:
                        m = re.search("T\d+_\d+.\d+.\d+.(\S{4}_\S+).txt", threading_file)
                        if m:
                            pdb_chain = m.group(1)
                        # For each helix... #
                        for helix in functions.parse_file(os.path.join(options.pdb_dir, "helices", pdb_chain + ".txt")):
                            # Get PDB object #
                            pdb_obj = PDB(os.path.join(options.pdb_dir, "split", pdb_chain[:4] + ".dna." + helix + ".pdb"))
                            break
                        # Get X3DNA object #
                        x3dna_obj = x3dna.X3DNA(os.path.join(options.pdb_dir, "x3dna", pdb_chain[:4] + ".txt"))
                        # Get interface object #
                        interface_obj = interface.Interface(os.path.join(options.pdb_dir, "interfaces", pdb_chain + ".txt"))
                        interface_basepairs = interface_obj.get_interface_basepairs()
                        # Get triads object #
                        triads_obj = triads.Triads(os.path.join(options.pdb_dir, "triads", pdb_chain + ".txt"))
                        # Initialize #
                        protein = {}
                        kmers = {}
                        # For each line... #
                        thread_obj = threader.Threaded(os.path.join(os.path.join(options.output_dir, "threading", threading_file)))
                        protein    = thread_obj.get_protein()
                        kmers      = thread_obj.get_kmers()
                        kmers_fixed= thread_obj.get_kmers_fixed()                        
                        # For each k-mer... #
                        for kmer, position in kmers.iteritems():
                            # Initialize #
                            basepairs = {}
                            kmerl = list(kmer)
                            # For each contact... #
                            try:
                             for j in range(int(position), int(position) + 8):
                                if j>=len(interface_basepairs):continue
                                basepairs.setdefault(interface_basepairs[j], kmerl.pop(0))
                            except:
                             if options.verbose: sys.stdout.write("Failed on %s (Basepairs in interface  %s Position %s in KMER  %s  )\n"%(threading_file,str(interface_basepairs),str(position),str(kmer)))
                             continue
                            # Skip if triads file already exists #
                            triads_file = os.path.join(options.output_dir, "triads", tf_obj.get_id() + "." + str(i) + "." + pdb_chain + "." + "".join([basepairs[j] for j in sorted(basepairs)]) + "." + str(min([j for j in basepairs])) + "-" + str(max([j for j in basepairs])) + ".txt")
                            if not os.path.exists(triads_file):
                                # Initialize #
                                done = set()
                                triads_thread = triads.Triads() 
                                # For each triad object... #
                                for triad_obj in triads_obj.get_triads():
                                    # Skip if already done #
                                    if (triad_obj._residue_A, triad_obj._residue_B) in done: continue
                                    # Initialize #
                                    aminoacid_environment = []
                                    dinucleotide_environment = []
                                    # Get chain, number #
                                    chain, number = triad_obj._residue_A.split("-")
                                    # If amino acid exists... #
                                    if (chain, int(number)) in protein:
                                        if protein[(chain, int(number))] == "-" or protein[(chain, int(number))] == "X": continue
                                        aminoacid_environment = triad_obj._A_environment.split("-")
                                        aminoacid_environment[0] = aminoacids1to3[protein[(chain, int(number))]]
                                        if aminoacids_polarity_boolean[protein[(chain, int(number))]]:
                                            aminoacid_environment[1] = "P"
                                        else:
                                            aminoacid_environment[1] = "N"
                                    # Get nucleotides #
                                    nucleotides = triad_obj._residue_B.split(",")
                                    # Initialize #
                                    dinucleotide = []
                                    # For each nucleotide... #
                                    for nucleotide in nucleotides:
                                        # Get chain, number #
                                        chain, number = nucleotide.split("-")
                                        # If nucleotide is among k-mer basepairs... #
                                        if x3dna_obj.get_residue_basepair(chain, int(number)) in basepairs:
                                            if x3dna_obj.get_residue_basepair(chain, int(number)) not in dinucleotide:
                                                dinucleotide.append(x3dna_obj.get_residue_basepair(chain, int(number)))
                                    # If dinucleotide exists... #
                                    if len(dinucleotide) == 2:
                                        dinucleotide_environment = triad_obj._B_environment.split("-")
                                        dinucleotide_environment[0] = basepairs[dinucleotide[0]] + basepairs[dinucleotide[1]]
                                        dinucleotide_environment[1] = "".join(nitrogenous_bases[nucleotide] for nucleotide in dinucleotide_environment[0])
                                    # If environments... #
                                    if len(aminoacid_environment) == 4 and len(dinucleotide_environment) == 5:
                                        triad_obj = triads.Triad("-".join(aminoacid_environment), "-".join(dinucleotide_environment), triad_obj._distance, triad_obj._residue_A, triad_obj._residue_B)
                                        triads_thread.add_triad(triad_obj)
                                        # Get complementary triad #
                                        dinucleotide_environment[0] = triads.get_complementary_dna_sequence(dinucleotide_environment[0])
                                        dinucleotide_environment[1] = "".join(nitrogenous_bases[i] for i in dinucleotide_environment[0])
                                        if dinucleotide_environment[2] == "F": dinucleotide_environment[2] = "R"
                                        else: dinucleotide_environment[2] = "F"
                                        triads_thread.add_triad(triads.Triad("-".join(aminoacid_environment), "-".join(dinucleotide_environment), triad_obj._distance, triad_obj._residue_A, triad_obj._residue_B))
                                    # Add to done #
                                    done.add((triad_obj._residue_A, triad_obj._residue_B))
                                # Write triads in file
                                triads_thread.write(triads_file)
                                
    # Exit if stops here #
    if options.parallel and options.start_step<3:
        sys.stdout.write("Exiting: wait until all submitted runs have finnished and restart again with start=3\n")
        exit(0)
    # Exit if stops here #
    if options.stop_step == 2:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ##########################################
    # 3. Create multiple sequence alignments #
    ##########################################
    if options.verbose: sys.stdout.write("Create multiple sequence alignments...\n\n")
    # Skip if starts later #
    if options.start_step <= 3:
       if options.parallel:
        if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
        else: cluster_queue=config.get("Cluster", "cluster_queue")
        functions.submit_command_to_queue("%s %s --pdb=%s -o %s --dummy=%s -t %s -e %s -f %s -m %s --pwm %s -s %s -u %s --start=3 --stop=3 -v" % ( os.path.join(config.get("Paths", "python_path"), "python"),  os.path.join(scripts_path,os.path.basename(__file__)) , os.path.abspath(options.pdb_dir), os.path.abspath(options.output_dir) , os.path.abspath(options.dummy_dir), os.path.abspath(options.tfs_file), os.path.abspath(options.escores_file), os.path.abspath(options.families_file), os.path.abspath(options.motifs_file), os.path.abspath(options.pwm_dir), os.path.abspath(options.sources_file), os.path.abspath(options.uniprot_file)),cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
       else:
        # For each cluster file... #
        for nr_file in sorted(os.listdir(os.path.join(options.pdb_dir, "nr"))):
            if nr_file == "general.txt": continue
            # Get PDB chain #
            pdb_chain = nr_file[:6]
            # Skip if MSA file already exists #
            msa_file = os.path.join(options.output_dir, "msa", pdb_chain + ".sto")
            if not os.path.exists(msa_file):
                # Initialize #
                sequences = set()
                # For each threading file... #
                for threading_file in os.listdir(os.path.join(options.output_dir, "threading")):
                    if pdb_chain in threading_file:
                        # Initialize #
                        sequence = ""
                        m = re.search("(\S+).%s.txt" % pdb_chain, threading_file)
                        # For each line... #
                        for line in functions.parse_file(os.path.join(os.path.join(options.output_dir, "threading", threading_file))):
                            if line == "": continue
                            if line.startswith("#") or line.startswith("//"): continue
                            elif line.startswith(">"):
                                threading = line[1:]
                            elif threading == "protein":
                                line = line.split(";")
                                sequence += line[2]
                        sequences.add((m.group(1), sequence))
                # If enough sequences... #
                if len(sequences) > 0:
                    # Initialize #
                    dummy_file = os.path.join(options.dummy_dir, "sequences.fa")
                    shutil.copy(os.path.join(options.pdb_dir, "split", pdb_chain + ".fasta"), dummy_file)
                    # For each sequence... #
                    for header, sequence in sequences:
                        functions.write(dummy_file, ">%s\n%s" % (header, sequence.replace("-", "")))
                    # Build MSA #
                    clustalo_path =  config.get("Paths", "clustalo_path")
                    process = subprocess.check_output([os.path.join(clustalo_path, "clustalo"), "-i", dummy_file, "-o", msa_file, "--outfmt=stockholm"], stderr=subprocess.STDOUT)
                    os.remove(dummy_file)
    # Exit if stops here #
    if options.parallel and options.start_step<4:
        sys.stdout.write("Exiting: wait until all submitted runs have finnished and restart again with start=4\n")
        exit(0)
    # Exit if stops here #
    if options.stop_step == 3:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ############################
    # 4. Create HMMER profiles #
    ############################
    if options.verbose: sys.stdout.write("Create hidden-Markov models...\n\n")
    # Skip if starts later #
    if options.start_step <= 4:
       if options.parallel:
        if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
        else: cluster_queue=config.get("Cluster", "cluster_queue")
        for msa_file in sorted(os.listdir(os.path.join(options.output_dir, "msa"))):
            # Get PDB chain #
            pdb_chain = msa_file[:6]
            # Skip if HMM file exists #
            src_path   = config.get("Paths", "src_path")
            hmmer_path = os.path.join(src_path, config.get("Paths", "hmmer_path"))
            hmmbuild   = os.path.join(hmmer_path, "hmmbuild")
            hmmer_file = os.path.join(os.path.abspath(options.output_dir), "hmm", pdb_chain + ".hmm")
            msa_file   = os.path.join(os.path.abspath(options.output_dir), "msa", msa_file)
            if not os.path.exists(hmmer_file):
              functions.submit_command_to_queue(("%s %s  %s "%(hmmbuild,hmmer_file,msa_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
       else:
        # For each msa file... #
        for msa_file in sorted(os.listdir(os.path.join(options.output_dir, "msa"))):
            # Get PDB chain #
            pdb_chain = msa_file[:6]
            # Skip if HMM file exists #
            hmmer_file = os.path.join(options.output_dir, "hmm", pdb_chain + ".hmm")
            if not os.path.exists(hmmer_file):
                # Build HMM #
                src_path = config.get("Paths", "src_path")
                hmmer_path = os.path.join(src_path, config.get("Paths", "hmmer_path"))
                process = subprocess.check_output([os.path.join(hmmer_path, "hmmbuild"), hmmer_file, os.path.join(options.output_dir, "msa", msa_file)], stderr=subprocess.STDOUT)
    # Exit if stops here #
    if options.parallel and options.start_step<5:
        sys.stdout.write("Exiting: wait until all submitted runs have finnished and restart again with start=5\n")
        exit(0)
    # Exit if stops here #
    if options.stop_step == 4:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ################################
    # 5. Create non-redundant sets #
    ################################
    if options.verbose: sys.stdout.write("Create non-redundant sets...\n")
    # Skip if starts later #
    if options.start_step <= 5:
        # Create indices for nr.py #
        threading_index = os.path.join(options.output_dir, "threading", "index.txt")
        triads_index = os.path.join(options.output_dir, "triads", "index.txt")
        if os.path.exists(threading_index):
           previous=threading_index+".previous"
           if os.path.exists(previous): os.remove(previous)
           shutil.move(threading_index,previous)
        if os.path.exists(triads_index):
           previous=triads_index+".previous"
           if os.path.exists(previous): os.remove(previous)
           shutil.move(triads_index,previous)
        if not os.path.exists(threading_index):
            if options.verbose: sys.stdout.write("\t-- create index file for threading...\n")
            n=0
            for file_name in os.listdir(os.path.join(options.output_dir, "threading")):
                if options.verbose:  
                   sys.stdout.write(".")
                   n += 1
                   if n>100:
                      n=0 
                      sys.stdout.write("\n")
                functions.write(threading_index, file_name)
            if options.verbose:sys.stdout.write("\n")
        if not os.path.exists(triads_index):
            if options.verbose: sys.stdout.write("\t-- create index file for triads...\n")
            n=0
            for file_name in os.listdir(os.path.join(options.output_dir, "triads")):
                if options.verbose:  
                   sys.stdout.write(".")
                   n += 1
                   if n>100:
                      n=0 
                      sys.stdout.write("\n")
                functions.write(triads_index, file_name)
            if options.verbose:sys.stdout.write("\n")

        ##############################
        # 5.1 General nr set         #
        ##############################
        if options.verbose: sys.stdout.write("\t\t-- nr general set...\n")
        # Skip if nr file already exists #
        nr_file = os.path.join(options.output_dir, "nr", "general.txt")
        if not os.path.exists(nr_file):
            # If parallelized #
            if options.parallel:
                # Submit to queue #
                if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
                else: cluster_queue=config.get("Cluster", "cluster_queue")
                functions.submit_command_to_queue("%s %s --pdb=%s -r %s -o %s --pbm=%s  -t %s" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "nr.py"), options.pdb_dir, os.path.abspath(os.getcwd()), nr_file, options.output_dir,float(config.get("Parameters", "max_redundancy_general"))), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
            # Else... #
            else:
                # Get non-redundant triads #
                nr_triads = nr.get_nr_triads(options.pdb_dir, threshold=float(config.get("Parameters", "max_redundancy_general")), pbm_dir=options.output_dir)
                # For each nr triads object... #
                for nr_triads_obj in nr_triads:
                    functions.write(nr_file, nr_triads_obj._file)

        ##############################
        # 5.2 Family nr sets         #
        ##############################
        if options.verbose: sys.stdout.write("\t\t-- nr family sets...\n")
        # For each triads file... #
        for triads_file in os.listdir(os.path.join(options.pdb_dir, "triads")):
            # Skip if nr file already exists #
            m = re.search("(\S{4}_\S).txt$", triads_file)
            nr_file = os.path.join(options.output_dir, "nr", m.group(1) + ".txt")
            if not os.path.exists(nr_file):
                # If parallelized #
                if options.parallel:
                    # Submit to queue #
                    if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
                    else: cluster_queue=config.get("Cluster", "cluster_queue")
                    functions.submit_command_to_queue("%s %s --pdb=%s -r %s -i %s -o %s --pbm=%s -t %s" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "nr.py"), options.pdb_dir, os.path.abspath(os.getcwd()), m.group(1), nr_file, options.output_dir, config.get("Parameters", "max_redundancy_family")), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                # Else... #
                else:
                    # Get non-redundant triads #
                    nr_triads = nr.get_nr_triads(options.pdb_dir, pdb_chain=m.group(1), threshold=float(config.get("Parameters", "max_redundancy_family")), pbm_dir=options.output_dir)
                    # For each nr triads object... #
                    for nr_triads_obj in nr_triads:
                        functions.write(nr_file, nr_triads_obj._file)
    # Verbose mode #
    if options.verbose: sys.stdout.write("\n")
    # Exit if stops here #
    if options.parallel and options.start_step<6:
        sys.stdout.write("Exiting: wait until all submitted runs have finnished and restart again with start=6\n")
        exit(0)
    if options.stop_step == 5:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ################################
    # 6. Derive stat. potentials   #
    ################################
    if options.verbose: sys.stdout.write("Derive statistical potentials...\n\n")
    # Skip if starts later #
    if options.start_step <= 6:
        # For each nr file... #
        for nr_file in os.listdir(os.path.join(options.output_dir, "nr")):
            # Initialize #
            m = re.search("(\S+).txt$", nr_file)
            # Set labels' dictionary
            label={}
            label.setdefault((True,True,True),'.taylor.bins')
            label.setdefault((False,True,True),'.bins')
            label.setdefault((False,False,True),'.pmf.bins')
            label.setdefault((True,False,True),'.pmf.taylor.bins')
            label.setdefault((True,True,False),'.taylor.acc')
            label.setdefault((False,True,False),'.acc')
            label.setdefault((False,False,False),'.pmf.acc')
            label.setdefault((True,False,False),'.pmf.taylor.acc')
            for approach,zscores,bins in label:
             if options.verbose: sys.stdout.write("\tcheck %s\n"%( m.group(1) + label[(approach,zscores,bins)] + ".txt" ))
             # Skip if potentials file already exists #
             potentials_file = os.path.join(options.output_dir, "potentials", m.group(1) + label[(approach,zscores,bins)] + ".txt")
             if not os.path.exists(potentials_file):
                # If parallelized #
                if options.parallel:
                    if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
                    else: cluster_queue=config.get("Cluster", "cluster_queue")
                    # Submit to queue #
                    if bins:
                      if approach and zscores:
                       functions.submit_command_to_queue("%s %s -i %s -o %s -s -a -z -b " % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "spotentials.py"), os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), os.path.abspath(potentials_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                      if not approach and zscores:
                       functions.submit_command_to_queue("%s %s -i %s -o %s -s -z -b " % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "spotentials.py"), os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), os.path.abspath(potentials_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                      if approach and not zscores:
                       functions.submit_command_to_queue("%s %s -i %s -o %s -s -a -b " % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "spotentials.py"), os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), os.path.abspath(potentials_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                      if not approach and not zscores:
                       functions.submit_command_to_queue("%s %s -i %s -o %s -s -b " % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "spotentials.py"), os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), os.path.abspath(potentials_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                    else:
                      if approach and zscores:
                       functions.submit_command_to_queue("%s %s -i %s -o %s -s -a -z  " % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "spotentials.py"), os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), os.path.abspath(potentials_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                      if not approach and zscores:
                       functions.submit_command_to_queue("%s %s -i %s -o %s -s -z  " % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "spotentials.py"), os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), os.path.abspath(potentials_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                      if approach and not zscores:
                       functions.submit_command_to_queue("%s %s -i %s -o %s -s -a  " % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "spotentials.py"), os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), os.path.abspath(potentials_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                      if not approach and not zscores:
                       functions.submit_command_to_queue("%s %s -i %s -o %s -s  " % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "spotentials.py"), os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), os.path.abspath(potentials_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                # Else... #
                else:
                    # Derive statistical potentials #
                    pmf_3d, pmf_3dc, pmf_s3dc, pmf_s3dc_dd, pmf_s3dc_di, pmf_local, pmf_pair, distances = spotentials.get_statistical_potentials(os.path.abspath(os.path.join(options.output_dir, "nr", nr_file)), approach=approach, smooth=True, zscores=zscores,computation=bins, dummy_dir=options.dummy_dir)
                    # Initialize #
                    potentials_obj = spotentials.Potentials()
                    # Add statistical potentials to object #
                    potentials_obj._pmf_3d = pmf_3d
                    potentials_obj._pmf_3dc = pmf_3dc
                    potentials_obj._pmf_s3dc = pmf_s3dc
                    potentials_obj._pmf_s3dc_dd = pmf_s3dc_dd
                    potentials_obj._pmf_s3dc_di = pmf_s3dc_di
                    potentials_obj._pmf_local = pmf_local
                    potentials_obj._pmf_pair = pmf_pair
                    potentials_obj._distances = distances
                    # Write statistical potentials #
                    potentials_obj.write(potentials_file)
    # Exit if stops here #
    if options.parallel and options.start_step<7:
        sys.stdout.write("Exiting: wait until all submitted runs have finnished and restart again with start=7\n")
        exit(0)

    if options.stop_step == 6:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ####################################
    # 7. Assign families to PDB chains #
    ####################################
    if options.verbose: sys.stdout.write("Assign families to PDB chains...\n\n")
    # Skip if starts later #
    if options.start_step <= 7:
        # Initialize #
        pdb_chains = {}
        families_file = os.path.join(options.output_dir, "families.txt")
        # Remove families file #
        if not os.path.exists(families_file):
           if options.parallel:
            if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
            else: cluster_queue=config.get("Cluster", "cluster_queue")
            functions.submit_command_to_queue("%s %s --pdb=%s -o %s --dummy=%s -t %s -e %s -f %s -m %s --pwm %s -s %s -u %s --start=7 -v" % ( os.path.join(config.get("Paths", "python_path"), "python"),  os.path.join(scripts_path,os.path.basename(__file__)) , os.path.abspath(options.pdb_dir), os.path.abspath(options.output_dir) , os.path.abspath(options.dummy_dir), os.path.abspath(options.tfs_file), os.path.abspath(options.escores_file), os.path.abspath(options.families_file), os.path.abspath(options.motifs_file), os.path.abspath(options.pwm_dir), os.path.abspath(options.sources_file), os.path.abspath(options.uniprot_file)),cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
           else:
            # For each line... #
            for line in functions.parse_file(os.path.join(options.pdb_dir, "families.txt")):
                if line.startswith("#"): continue
                line = line.split(";")
                pdb_chains.setdefault(line[0], set())
                pdb_chains[line[0]].add(line[-1])
            # For each PDB chain... #
            for pdb_chain in pdb_chains:
                if "Unknown" not in pdb_chains[pdb_chain]: continue
                # Initialize #
                pdb_chains[pdb_chain] = set()
                # For each threading file... #
                for threading_file in os.listdir(os.path.join(options.output_dir, "threading")):
                    if pdb_chain not in threading_file: continue
                    # Get TF id and PDB chain #
                    m = re.search("(.+).\d+..{6}.txt", threading_file)
                    # Initialize #
                    tf_id = m.group(1)
                    # Get TF object #
                    tf_obj = TF(file_name=os.path.join(options.output_dir, "tfs", tf_id + ".txt"))
                    # For each family... #
                    for family in tf_obj._family:
                        pdb_chains[pdb_chain].add(family)
            # For each PDB chain... #
            for pdb_chain in pdb_chains:
                if len(pdb_chains[pdb_chain]) == 0:
                    pdb_chains[pdb_chain].add("Unknown")
            # Create output #
            functions.write(families_file, "#pdb_chain;family")
            # For each PDB chain... #
            for pdb_chain in sorted(pdb_chains):
                # Write #
                functions.write(families_file, "%s;%s" % (pdb_chain, ",".join(sorted(pdb_chains[pdb_chain]))))

    # Exit if stops here #
    if options.parallel and options.start_step<8:
        sys.stdout.write("Exiting: wait until all submitted runs have finnished and restart again with start=8\n")
        exit(0)

    # Exit if stops here #
    if options.stop_step == 7:
        sys.stdout.write("Exiting...  \n\n")
        exit(0)

    ##############################
    # 8. Scan UniProt database   #
    ##############################
    if options.verbose: sys.stdout.write("Scan UniProtKB using blast/HMMer...\n\n")
    # Skip if starts later #
    if options.start_step <= 8:
        # Initialize #
        homologs_hmm =  config.get("Parameters","homologs_hmm") == "yes" or  config.get("Parameters","homologs_hmm") == "YES"
        src_path = config.get("Paths", "src_path")
        blast_path = os.path.join(src_path, config.get("Paths", "blast_path"))
        hmmer_path = os.path.join(src_path, config.get("Paths", "hmmer_path"))
        tz_parameter = int(config.get("Parameters", "twilight_zone_parameter"))
        tz_type = config.get("Parameters", "twilight_zone_type")
        # For each triads file... #
        for triads_file in os.listdir(os.path.join(options.pdb_dir, "triads")):
            # Get PDB chain #
            m = re.search("(.{6}).txt", triads_file)
            pdb_chain = m.group(1)
            # Skip if uniprot file already exists #
            uniprot_file_blast = os.path.join(options.output_dir, "uniprot", pdb_chain + ".blast")
            uniprot_file_hmm = os.path.join(options.output_dir, "uniprot", pdb_chain + ".hmm")
            if not os.path.exists(uniprot_file_hmm) or not  os.path.exists(uniprot_file_blast):
                if options.verbose: sys.stdout.write("\t-- scan homologs of %s in %s\n"%(pdb_chain,uniprot_file_blast))
                # If parallelized #
                if options.parallel:
                    if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
                    else: cluster_queue=config.get("Cluster", "cluster_queue")
                    # If HMMER file exists... #
                    hmm_file = os.path.join(options.output_dir, "hmm", pdb_chain + ".hmm")
                    if os.path.exists(hmm_file) and not os.path.exists(uniprot_file_hmm) and homologs_hmm:
                        # Submit to queue #
                        homolog_file = os.path.join(options.output_dir, "uniprot", pdb_chain + ".hmm")
                        functions.submit_command_to_queue("%s %s -d %s -i %s -f -o %s" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "hmmer.py"), os.path.abspath(options.uniprot_file), os.path.abspath(hmm_file), os.path.abspath(homolog_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                    # Get BLAST object #
                    if not os.path.exists(uniprot_file_blast):
                        fasta_file = os.path.join(options.pdb_dir, "split", pdb_chain + ".fasta")
                        # Submit to queue #
                        homolog_file = os.path.join(options.output_dir, "uniprot", pdb_chain + ".blast")
                        functions.submit_command_to_queue("%s %s -d %s -i %s -f -o %s" % (os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "blast.py"), os.path.abspath(options.uniprot_file), os.path.abspath(fasta_file), os.path.abspath(homolog_file)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
                # Else... #
                else:
                    # If Hmmer file exists... #
                    hmm_file = os.path.join(options.output_dir, "hmm", pdb_chain + ".hmm")
                    if os.path.exists(hmm_file) and not os.path.exists(uniprot_file_hmm) and homologs_hmm:
                        # Get HMMER object #
                        hmmer_obj = hmmer.get_hmmer_obj(options.uniprot_file, hmm_file)
                        # Write uniprot file #
                        homolog_file = os.path.join(options.output_dir, "uniprot", pdb_chain + ".hmm")
                        hmmer_obj.write(homolog_file)
                    # Get BLAST object #
                    if not os.path.exists(uniprot_file_blast):
                        blast_obj = blast.get_blast_obj(options.uniprot_file, os.path.join(options.pdb_dir, "split", pdb_chain + ".fasta"))
                        # Write output #
                        homolog_file = os.path.join(options.output_dir, "uniprot", pdb_chain + ".blast")
                        functions.write(homolog_file, blast_obj.str_compacted_blast(tz_parameter=tz_parameter, tz_type=tz_type))
    # Exit if stops here #
    if options.parallel and options.start_step<9:
        sys.stdout.write("Exiting: wait until all submitted runs have finnished and restart again with start=9\n")
        exit(0)

    if options.stop_step == 8:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ##############################
    # 9. Filter UniProt homologs #
    ##############################
    if options.verbose: sys.stdout.write("Filter UniProtKB by interface and generate Homologs...\n\n")
    # Skip if starts later #
    if options.start_step <= 9:
       # Initialize #
       homologs_hmm =  config.get("Parameters","homologs_hmm") == "yes" or  config.get("Parameters","homologs_hmm") == "YES"
       if options.parallel:
        if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
        else: cluster_queue=config.get("Cluster", "cluster_queue")
        for homologs_file in sorted(os.listdir(os.path.join(options.output_dir, "uniprot"))):
          if (homologs_file.endswith(".blast") and not homologs_hmm) or (homologs_file.endswith(".hmm") and homologs_hmm):
            pdb_chain  = homologs_file[:6]
            if options.verbose: sys.stdout.write("%s...\n" % pdb_chain)
            dummy_file = os.path.join(os.path.abspath(options.output_dir), "homologs", pdb_chain + ".txt")
            if not os.path.exists(dummy_file):
              functions.submit_command_to_queue("%s %s -i %s -o %s --pdb=%s --pbm=%s --dummy=%s --filter" % ( os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "homologs.py"), os.path.join(os.path.abspath(options.output_dir), "uniprot",homologs_file), os.path.join(os.path.abspath(options.output_dir), "homologs", pdb_chain + ".txt"), os.path.abspath(options.pdb_dir), os.path.abspath(options.output_dir) , os.path.abspath(options.dummy_dir)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
            else:
              if options.verbose: sys.stdout.write("...Done...\n")
       else:
        # For each homologs file... #
        for homologs_file in sorted(os.listdir(os.path.join(options.output_dir, "uniprot"))):
           if (homologs_file.endswith(".blast") and not homologs_hmm) or (homologs_file.endswith(".hmm") and homologs_hmm):
            # Initialize #
            pdb_chain  = homologs_file[:6]
            dummy_file = os.path.join(options.output_dir, "homologs", pdb_chain + ".txt")
            # Verbose #
            if options.verbose: sys.stdout.write("%s...\n" % pdb_chain)
            # Skip if dummy file exist #
            if not os.path.exists(dummy_file):
                # Initialize #
                homologs_obj  = homologs.Homologs(os.path.join(options.output_dir, "uniprot",homologs_file))
                filtered      = homologs_obj.filter_hits_by_interface(options.pdb_dir,options.dummy_dir)
                #replace the set of homologs
                homologs_obj._homologs=filtered
                # write output
                homologs_obj.write(dummy_file)
            else:
              if options.verbose: sys.stdout.write("...Done...\n")

    # Exit if stops here #
    if options.parallel and options.start_step < 10:
        sys.stdout.write("Exiting: wait until all submitted jobs are completed and restart again with start=10\n\n")
        exit(0)

    if options.stop_step == 9:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ######################
    # 10. Calculate PWMs #
    ######################
    if options.verbose: sys.stdout.write("Calculate position weight matrices...\n\n")
    # Skip if starts later #
    if options.start_step <= 10:
       #Initialize
       if not os.path.exists(os.path.join(os.path.abspath(options.dummy_dir), "pwms")):
         os.makedirs(os.path.join(os.path.abspath(options.dummy_dir), "pwms"))
       radius=float(options.radius)
       if radius <= 0: radius=float(config.get("Parameters", "max_contact_distance"))
       #Read dimers
       if options.verbose:sys.stdout.write("\t\tget dimers...\n")
       dimers={}
       dimers_file=open(os.path.join(os.path.abspath(options.pdb_dir),"dimers.txt"),"r")
       for line in dimers_file:
         if line.startswith("#"):continue
         chain_A,chain_B,contacts,overlap=line.strip().split(";")
         dimers.setdefault(chain_A,set().add(chain_B))
         dimers.setdefault(chain_B,set().add(chain_A))
       if options.verbose:sys.stdout.write("\t\tget families...\n")
       families={}
       for line in functions.parse_file(os.path.join(options.pdb_dir, "families.txt")):
         if line.startswith("#"): continue
         pdb_name, family = line.split(";")
         families.setdefault(pdb_name,family)
       #Read chains in pdb/split folder
       for pdb_file in sorted(os.listdir(os.path.join(options.pdb_dir, "split"))):
         skip=False
         code = pdb_file.strip().split(".")
         if code[-1] != "pdb": continue
         #get the PDB file
         pdb_code  = pdb_file[:4]
         try:
           pdb_obj   = PDB(os.path.join(options.pdb_dir,"clean",pdb_code+".pdb"))
         except:
           skip=True
           continue
         complex_chains=set()
         #Select the chains that bind the same dna helix
         if code[1] == "dna": 
            dna_helix = code[2]
            for chain_id in pdb_obj.chain_identifiers:
               if not os.path.exists(os.path.join(os.path.abspath(options.pdb_dir),"helices",pdb_code+"_"+chain_id+".txt")): continue
               helix_file= open(os.path.join(os.path.abspath(options.pdb_dir),"helices",pdb_code+"_"+chain_id+".txt"),"r")
               for helix in helix_file:
                 if helix.strip()==dna_helix:  complex_chains.add(chain_id)
               helix_file.close()
         else:
            chain_id=code[0].split("_")[-1]
            if os.path.exists(os.path.join(os.path.abspath(options.pdb_dir),"helices",pdb_code+"_"+chain_id+".txt")): 
              complex_chains.add(chain_id)
              helix_file= open(os.path.join(os.path.abspath(options.pdb_dir),"helices",pdb_code+"_"+chain_id+".txt"),"r")
              for helix in helix_file:
                dna_helix=helix.strip()
              helix_file.close()
         #Check conditions to continue: exist protein and exist DNA
         if len( complex_chains ) <= 0: continue 
         if not os.path.exists(os.path.join(os.path.abspath(options.pdb_dir),"split",pdb_code+".dna."+dna_helix+".pdb")): continue
         #Create a dummy PDB file with the protein-DNA complex
         try:
           dummy_pdb_obj = PDB()
           dna_pdb_obj   = PDB(os.path.join(os.path.abspath(options.pdb_dir),"split",pdb_code+".dna."+dna_helix+".pdb"))
           for dna_chain in dna_pdb_obj.chains:
             dummy_pdb_obj.add_chain(dna_chain)
           pdb_chain=pdb_code+"_"+dna_helix
           for chain_id in complex_chains:
             pdb_chain +=   chain_id 
             chain      =   pdb_obj.get_chain_by_id(chain_id)
             dummy_pdb_obj.add_chain(chain)
           dummy_file = os.path.join(os.path.abspath(options.dummy_dir), "pwms", pdb_chain + ".pdb" )
           if not os.path.exists(dummy_file): dummy_pdb_obj.write(dummy_file)
         except:
           skip=True
           continue
         #Start the PWM calculation
         if options.verbose: sys.stdout.write("%s...\n" % pdb_chain)
         output_file = os.path.join(os.path.abspath(options.output_dir), "pwms", pdb_chain)
         output_pwm  = output_file+".pwm"
         output_meme = output_file+".meme"
         output_msa  = output_file+".msa"
         output_logo = output_file+".logo"
         output_fwd_logo = output_file+".logo.fwd.png"
         output_rev_logo = output_file+".logo.rev.png"
         if os.path.exists(output_pwm) and os.path.exists(output_meme) and os.path.exists(output_msa) and os.path.exists(output_fwd_logo) and os.path.exists(output_rev_logo):
            skip=True
         if not skip:
            if options.parallel:
              if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
              else: cluster_queue=config.get("Cluster", "cluster_queue")
              functions.submit_command_to_queue("%s %s --radius=%f -i %s -o %s --pdb=%s --pbm=%s --auto -v --dummy=%s --known" % ( os.path.join(config.get("Paths", "python_path"), "python"), os.path.join(scripts_path, "pwm_pbm.py"), radius, dummy_file, output_file, os.path.abspath(options.pdb_dir), os.path.abspath(options.output_dir) , os.path.abspath(options.dummy_dir)), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
            else:
              if os.path.exists(output_msa):
                  msa_obj=PWM.nMSA(output_msa,pdb_chain)
              else:
                  # protein is dummy_pdb_obj
                  # Get DSSP object #
                  if options.verbose:sys.stdout.write("\t-- calculate secondary structure ...\n")
                  try:
                    dummy_dssp_obj = dssp.get_dssp_obj(dummy_file)
                  except:
                    continue
                  # Get X3DNA object #
                  if options.verbose:sys.stdout.write("\t-- get DNA %s ...\n"%pdb_chain)
                  try:
                    dummy_x3dna_obj = x3dna.get_x3dna_obj(dummy_file)
                  except:
                    continue
                  # Get contacts object #
                  if options.verbose:sys.stdout.write("\t-- calculate contacts ...\n")
                  try:
                    dummy_contacts_obj = contacts.get_contacts_obj(dummy_pdb_obj, dummy_x3dna_obj)
                  except:
                    continue
                  # Get triads object #
                  if options.verbose:sys.stdout.write("\t-- calculate protein-dna pairs ...\n")
                  try:
                    dummy_triads_obj = triads.get_triads_obj(dummy_pdb_obj, dummy_dssp_obj, dummy_x3dna_obj, dummy_contacts_obj)
                  except:
                    continue
                  # Load statistical potential #
                  if options.verbose:sys.stdout.write("\t-- load potentials ...\n")
                  try:
                    #potentials, thresholds, radii, structural_homologs_by_chain = PWM.load_statistical_potentials(dummy_pdb_obj, os.path.abspath(options.pdb_dir), os.path.abspath(options.output_dir), families, options.radius, potential_file=None, split_potential="s3dc_dd", auto_mode=True, family_potentials=False, pbm_potentials=False,  score_threshold=None, taylor_approach=False, pmf_approach=False, bins_approach=False, known_pdb=True , None, dummy_dir=os.path.abspath(options.dummy_dir),options.verbose)
                    potentials, thresholds, radii, structural_homologs_by_chain = PWM.load_statistical_potentials(dummy_pdb_obj, os.path.abspath(options.pdb_dir), os.path.abspath(options.output_dir), families, options.radius, None, "s3dc_dd", True, False, False,  None, False, False, False, True , None, os.path.abspath(options.dummy_dir),options.verbose)
                  except:
                    continue
                  # Get MSA object #
                  if options.verbose:sys.stdout.write("\t-- get MSA object ...\n")
                  try:
                    msa_obj = PWM.get_msa_obj(dummy_triads_obj, dummy_x3dna_obj, potentials, radii, None, None, "s3dc_dd" , thresholds)
                  except:
                    continue
                  motif_name=pdb_chain
                  msa_obj.set_motif(motif_name)
                  # Write PWM #
              if not os.path.exists(output_pwm):
                     if options.verbose:sys.stdout.write("\twrite PWM ...\n")
                     if options.verbose:sys.stdout.write("\t\t--PWM...\n")
                     msa_obj.write(output_pwm, option="pwm")
              if not os.path.exists(output_meme):
                     if options.verbose:sys.stdout.write("\t\t--PWM in MEME format...\n")
                     msa_obj.write(output_meme, option="meme")
                     #Skip meme format with uniprobe2meme. 
                     #Remove next commented line if you want both ways of obtaining the PWM
                     #PWM.write_pwm_by_meme(msa_obj,output_meme+".s",options.dummy_dir)
              if not os.path.exists(output_msa):
                     if options.verbose:sys.stdout.write("\t\t--MSA...\n")
                     msa_obj.write(output_msa, option="msa")
              if not os.path.exists(output_fwd_logo) or  not os.path.exists(output_rev_logo):
                     if options.verbose:sys.stdout.write("\t\t--Logos...\n")
                     PWM.write_logo(msa_obj,output_logo,options.dummy_dir)
              os.remove( dummy_file )
         else:
            if options.verbose: sys.stdout.write("...Done...\n")


    # Exit if stops here #
    if options.parallel and options.start_step < 11:
        sys.stdout.write("Exiting: wait until all submitted jobs are completed and restart again with start=11\n\n")
        exit(0)

    if options.stop_step == 10:
        sys.stdout.write("Exiting...\n\n")
        exit(0)

    ##########################################
    # 11. Check PWMs and generate a DATABASE #
    ##########################################

    if options.verbose: sys.stdout.write("Check PWM and MEME files and the completeness of the database...\n\n")
    # Skip if starts later #
    if options.start_step <= 11:
       #Initialize
       if not os.path.exists(os.path.join(os.path.abspath(options.output_dir), "pwms")):
         sys.stdout.write("Please run first with start<11 as PWM files are incompleted\n")
         exit(0)
       if not os.path.exists(os.path.join(os.path.abspath(options.dummy_dir), "pwms")):
         os.makedirs(os.path.join(os.path.abspath(options.dummy_dir), "pwms"))
       radius=float(options.radius)
       if radius <= 0: radius=float(config.get("Parameters", "max_contact_distance"))
       database_file = os.path.join(os.path.abspath(options.output_dir), "pwms", "database.txt")
       database_extended_file = os.path.join(os.path.abspath(options.output_dir), "pwms", "database_extended.txt")
       if  os.path.exists(database_extended_file): os.remove(database_extended_file)
       if  os.path.exists(database_file): os.remove(database_file)
       # Check TF families
       tf_codes=set()
       if os.path.exists(os.path.join(options.output_dir, "families.txt")):
          for line in functions.parse_file(os.path.join(options.output_dir, "families.txt")):
              pdb_name, family = line.split(";")
              if family != "Unknown" and family != "Undefined": tf_codes.add(pdb_name[0:4].lower())
       # Check PWM
       pwm_done=set()
       for pwm_file in os.listdir(os.path.join(os.path.abspath(options.output_dir), "pwms")):
           if pwm_file.endswith("pwm") or pwm_file.endswith("msa"):
                if pwm_file.endswith("pwm"): motif_name=pwm_file.rstrip(".pwm")
                if pwm_file.endswith("msa"): motif_name=pwm_file.rstrip(".msa")
                if motif_name in pwm_done: continue
                output_file=os.path.join(os.path.abspath(options.output_dir), "pwms",motif_name)
                output_meme = output_file+".meme"
                output_msa  = output_file+".msa"
                output_pwm  = output_file+".pwm"
                if not os.path.exists(output_pwm) and not os.path.exists(output_msa):
                   if options.verbose: sys.stdout.write("\t\t--Skip\n")
                   continue
                if not os.path.exists(output_msa): 
                   if options.verbose: sys.stdout.write("\t--Load %s \n"%(output_pwm))
                   msa_obj=PWM.nMSA(output_pwm,motif_name,option="pwm")
                else:
                   if options.verbose: sys.stdout.write("\t--Load %s \n"%(output_msa))
                   msa_obj=PWM.nMSA(output_msa,motif_name,option="msa")
                if msa_obj.get_binding_site_length() <= 0: 
                   if options.verbose: sys.stdout.write("\t--Check %s \n"%motif_name)
                   #If MSA does not exist we have to rebuilt as in option 10
                   if not os.path.exists(output_msa):
                      if options.verbose:sys.stdout.write("\t\t--Redo PWM and MSA of %s ...\n"%motif_name)
                      families={}
                      for line in functions.parse_file(os.path.join(options.pdb_dir, "families.txt")):
                        if line.startswith("#"): continue
                        pdb_name, family = line.split(";")
                        families.setdefault(pdb_name,family)
                      m = re.search("(\S+)_(\d+)(\S+)",motif_name)
                      if m:
                         pdb_code          = m.group(1)
                         dna_helix         = m.group(2)
                         complex_chains    = m.group(3)
                      else:
                         pdb_code          = motif_name[0:4]
                         dna_helix         = motif_name[6:6]
                         complex_chains    = motif_name[7:]
                      pdb_obj   = PDB(os.path.join(options.pdb_dir,"clean",pdb_code+".pdb"))
                      if not os.path.exists(os.path.join(os.path.abspath(options.pdb_dir),"split",pdb_code+".dna."+dna_helix+".pdb")): 
                         if options.verbose: sys.stdout.write("\t\t--Skip %s from PWMs\n"%motif_name)
                         continue
                      dummy_pdb_obj = PDB()
                      dna_pdb_obj   = PDB(os.path.join(os.path.abspath(options.pdb_dir),"split",pdb_code+".dna."+dna_helix+".pdb"))
                      for dna_chain in dna_pdb_obj.chains:
                        dummy_pdb_obj.add_chain(dna_chain)
                      pdb_chain=pdb_code+"_"+dna_helix
                      for chain_id in complex_chains:
                        pdb_chain +=   chain_id 
                        chain      =   pdb_obj.get_chain_by_id(chain_id)
                        dummy_pdb_obj.add_chain(chain)
                      dummy_file = os.path.join(os.path.abspath(options.dummy_dir), "pwms", pdb_chain + ".pdb" )
                      if not os.path.exists(dummy_file): dummy_pdb_obj.write(dummy_file)
                      # Get DSSP object #
                      if options.verbose:sys.stdout.write("\t\t-- calculate secondary structure %s ...\n"%pdb_chain)
                      try:
                        dummy_dssp_obj = dssp.get_dssp_obj(dummy_file)
                      except:
                        if options.verbose: sys.stdout.write("\t\t--Skip %s from PWMs\n"%motif_name)
                        continue
                      # Get X3DNA object #
                      if options.verbose:sys.stdout.write("\t\t-- get DNA %s ...\n"%pdb_chain)
                      try:
                        dummy_x3dna_obj = x3dna.get_x3dna_obj(dummy_file)
                      except:
                        if options.verbose: sys.stdout.write("\t\t--Skip %s from PWMs\n"%motif_name)
                        continue
                      # Get contacts object #
                      if options.verbose:sys.stdout.write("\t\t-- calculate contacts %s ...\n"%pdb_chain)
                      try:
                        dummy_contacts_obj = contacts.get_contacts_obj(dummy_pdb_obj, dummy_x3dna_obj)
                      except:
                        if options.verbose: sys.stdout.write("\t\t--Skip %s from PWMs\n"%motif_name)
                        continue
                      # Get triads object #
                      if options.verbose:sys.stdout.write("\t\t-- calculate protein-dna pairs ...\n")
                      try:
                        dummy_triads_obj = triads.get_triads_obj(dummy_pdb_obj, dummy_dssp_obj, dummy_x3dna_obj, dummy_contacts_obj)
                      except:
                        if options.verbose: sys.stdout.write("\t\t--Skip %s from PWMs\n"%motif_name)
                        continue
                      # Load statistical potential #
                      if options.verbose:sys.stdout.write("\t\t-- load potentials ...\n")
                      try:
                        potentials, thresholds, radii, structural_homologs_by_chain = PWM.load_statistical_potentials(dummy_pdb_obj, os.path.abspath(options.pdb_dir), os.path.abspath(options.output_dir), families, options.radius, None, "s3dc_dd", True, False, False,  None, False, False, False, True , None, os.path.abspath(options.dummy_dir),options.verbose)
                      except:
                        if options.verbose: sys.stdout.write("\t\t--Skip %s from PWMs\n"%motif_name)
                        continue
                      # Get MSA object #
                      if options.verbose:sys.stdout.write("\t\t-- get MSA object ...\n")
                      try:
                        msa_obj = PWM.get_msa_obj(dummy_triads_obj, dummy_x3dna_obj, potentials, radii, None, None,"s3dc_dd" , thresholds)
                      except:
                        if options.verbose: sys.stdout.write("\t\t--Skip %s from PWMs\n"%motif_name)
                        continue
                      msa_obj.set_motif(motif_name)
                      msa_obj.write(output_msa, option="msa")
                      output_logo = output_file+".logo"
                      output_fwd_logo = output_file+".logo.fwd.png"
                      output_rev_logo = output_file+".logo.rev.png"
                      if options.verbose:sys.stdout.write("\t\t-- Redo logos...\n")
                      try:
                         PWM.write_logo(msa_obj,output_logo,options.dummy_dir)
                      except:
                         if options.verbose:sys.stdout.write("\t\t--Rewrite MSA %s (msa has %d sequences and binding-site length is %d)\n"%(motif_name,len(msa_obj.get_sequences()),msa_obj.get_binding_site_length()))
                         msa_obj.write(output_msa, option="msa")
                   else:
                      msa_obj=PWM.nMSA(output_msa,motif_name)
                   if msa_obj.get_binding_site_length() > 0 and len(msa_obj.get_sequences())>0:
                      if options.verbose: sys.stdout.write("\t\t--Add %s in database (msa has %d sequences and binding-site length is %d)\n"%(motif_name,len(msa_obj.get_sequences()),msa_obj.get_binding_site_length()))
                      os.remove(output_meme)
                      os.remove(output_pwm)
                      msa_obj.set_motif(motif_name)
                      msa_obj.write(output_meme, option="meme")
                      msa_obj.write(output_pwm, option="pwm")
                      pwm_done.add(motif_name)
                      os.system("cat %s >> %s" %(output_meme,database_extended_file))
                      m = re.search("(\S+)_(\d+)(\S+)",motif_name)
                      if m:
                         chains =  m.group(3)
                         pdb_code = m.group(1).lower()
                         if len(chains) == 1 and pdb_code in tf_codes: os.system("cat %s >> %s" % (output_meme, database_file))
                   else:
                      pwm_done.add(motif_name)
                      if options.verbose: sys.stdout.write("\t\t--Skip %s from PWMs\n"%motif_name)
                      os.remove(output_meme)
                      os.remove(output_pwm)
                      os.remove(output_msa)
                else:
                   if msa_obj.get_binding_site_length() > 0 and len(msa_obj.get_sequences())>0:
                      if options.verbose: sys.stdout.write("\t\t--Add %s in database (msa has %d sequences and binding-site length is %d)\n"%(motif_name,len(msa_obj.get_sequences()),msa_obj.get_binding_site_length()))
                      os.remove(output_meme)
                      os.remove(output_pwm)
                      msa_obj.write(output_meme, option="meme")
                      msa_obj.write(output_pwm, option="pwm")
                      pwm_done.add(motif_name)
                      os.system("cat %s >> %s" %(output_meme,database_extended_file))
                      m = re.search("(\S+)_(\d+)(\S+)",motif_name)
                      if m:
                         chains =  m.group(3)
                         pdb_code = m.group(1).lower()
                         if len(chains) == 1 and pdb_code in tf_codes: os.system("cat %s >> %s" % (output_meme, database_file))
                   else:
                      pwm_done.add(motif_name)
                      if options.verbose: sys.stdout.write("\t\t--Skip %s from PWMs\n"%motif_name)
                      os.remove(output_meme)
                      os.remove(output_pwm)
                      os.remove(output_msa)
                      
                   

    # Exit if stops here #
    if options.parallel and options.start_step <= 11:
        sys.stdout.write("Exiting: wait until all submitted jobs are completed and remove dummy directory\n\n")
        exit(0)

    if options.stop_step == 11:
        sys.stdout.write("Done!!...\n\n")
        exit(0)

