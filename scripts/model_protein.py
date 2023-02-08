import os, sys, re
import ConfigParser
import optparse
import shutil
import subprocess
import socket
import time

#print(socket.gethostname())
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

# Import my functions #
import functions

# Imports jbonet's module #
from SBI.structure import PDB
from SBI.structure.chain import ChainOfProtein, ChainOfNucleotide
from SBI.structure.residue import ResidueOfNucleotide
from SBI.structure.atom import AtomOfNucleotide

# Import my modules #
import blast, contacts, dssp, interface, model_dna, x3dna,threader

#-------------#
# Classes     #
#-------------#

class PirAlignment(object):
    """
    This class defines a {PirAlignment} object.

    """

    def __init__(self, file_name=None):
        self._file = file_name
        self._hit = None
        self._sequence = None
        self._structure = None
        self._sequence_residues = []
        self._structure_residues = []
        self._structure_chains = []
        self._sequence_alignments = []
        self._structure_alignments = []
        self._heterodimer = None
        self._domain = ""
        self._identity = None
        self._similarity = None
        self._main_order = []

        # Initialize #
        if self._file is not None:
            self._parse_file()

    def _parse_file(self):
        # Alignment #
        for line in functions.parse_file(self._file):
            m = re.search("^structureX:(.{4}):(\d+):(.):(\d+):(.)::::", line)
            if m:
                self._structure = m.group(1)
                self._structure_residues.append(m.group(2))
                self._structure_residues.append(m.group(4))
                self._structure_residues.append(m.group(3))
                self._structure_residues.append(m.group(5))
                alignment = "structure"
            m = re.search("^sequence:(\S+)::::::::", line)
            if m:
                self._sequence = m.group(1)
                alignment = "sequence"
            m = re.search("^(.+)[/*]$", line)
            if m:
                if alignment == "structure":
                    self._structure_alignments.append(m.group(1))
                if alignment == "sequence":
                    self._sequence_alignments.append(m.group(1))

    def get_sequence_name(self):
        return self._sequence

    def get_structure_name(self):
        return self._structure

    def get_structure_chains(self):
        return sorted(self._structure_chains)

    def get_identity(self):
        return self._identity

    def get_similarity(self):
        return self._similarity

    def get_alignments(self, atype="sequence"):
        if atype != "sequence" and atype != "structure":
            atype = "sequence"

        if atype == "sequence":
            return self._sequence_alignments

        return self._structure_alignments

    def add_hit_name(self, name):
        self._hit = name

    def add_sequence_name(self, name):
        self._sequence = name

    def add_structure_name(self, name):
        self._structure = name

    def add_sequence_residue(self, residue_num):
        self._sequence_residues.append(residue_num)

    def add_structure_residue(self, residue_num):
        self._structure_residues.append(residue_num)

    def add_structure_chain(self, chain):
        self._structure_chains.append(chain)

    def add_sequence_alignment(self, alignment):
        self._sequence_alignments.append(alignment)

    def add_structure_alignment(self, alignment):
        self._structure_alignments.append(alignment)

    def add_domain(self, domain):
        self._domain = domain

    def get_domain(self):
        return self._domain 

    def add_main_order(self, main_order):
        self._main_order.append(main_order)

    def get_main_order(self):
        return self._main_order 

    def set_heterodimer(self, heterodimer):
        self._heterodimer = heterodimer

    def set_identity(self, identity):
        self._identity = identity

    def set_similarity(self, similarity):
        self._similarity = similarity

    def write(self, file_name):
        functions.write(file_name, ">P1;%s" % self._structure)
        functions.write(file_name, "structureX:%s:%s:%s:%s:%s::::" % (self._structure, self._structure_residues[0], self._structure_chains[0], self._structure_residues[-1], self._structure_chains[-1]))
        functions.write(file_name, "%s*" % "/\n".join([re.sub("[xX]", "-", i) for i in self._structure_alignments]))
        functions.write(file_name, ">P1;%s" % self._sequence)
        functions.write(file_name, "sequence:%s::::::::" % self._sequence)
        functions.write(file_name, "%s*" % "/\n".join([re.sub("[xX]", "-", i) for i in self._sequence_alignments]))

#-------------#
# Functions   #
#-------------#

def fileExist(file):
    '''
    Check existing files
    '''
    if file is not None:
        return os.path.exists(file) and os.path.isfile(file)
    else:
        return False

def make_subdirs(main, subdirs):
    '''
    This function makes all subdirs listed in "subdirs".
    '''

    for subdir in subdirs:
        if not os.path.exists(os.path.join(main, subdir)):
            os.makedirs(os.path.join(main, subdir))

def remove_files(files):
    '''
    This function removes all files listed in "files".
    '''

    for each_file in files:
        if os.path.exists(each_file):
            os.remove(each_file)

def bijection_sequence_to_sequence(sequences_a,sequences_b,dummy_dir="tmp",homodimer=False):
    ''' 
     Define the association of most similar sequences
     sequences_a is a dictionary of sequences {name:sequence}
     sequences_b is a dictionary of sequences {name:sequence}
     dummy_dir  Dummy directory to cerate file
     the output is a dictionary of associations {name_a:name_b} and {name_b:name_a}

    '''
    #Initialize
    from SBI.structure.chain import Chain
    from SBI.structure.chain import ChainOfProtein
    from SBI.structure.chain import ChainOfNucleotide
    from SBI.sequence import Sequence
    from SBI.structure import PDB
    from Bio import SeqIO
    from Bio import ExPASy
    from Bio import AlignIO
    from Bio.Align import Applications
        
    clustal_exe  = os.path.join(config.get('Paths','clustalo_path'),'clustalo')

    if len(sequences_a.keys())>len(sequences_b.keys()):
       main_set   = sequences_a
       second_set = sequences_b
    else: 
       main_set   = sequences_b
       second_set = sequences_a

    pairing={}
    used=set()
    for name_main,sequence_main in main_set.iteritems():
       score={}
       if  len(used) == len(second_set.keys()) and not homodimer:
           pairing.setdefault(name_main,None)
           continue
       for  name_second,sequence_second in second_set.iteritems():
           if name_second in used and not homodimer: continue
           if len(sequence_main) < len(sequence_second): shortest_length=len(sequence_main)
           else: shortest_length=len(sequence_second)
           infile       =dummy_dir+"/tmp_"+name_main+"_"+name_second+".fa"
           infile_query =dummy_dir+"/tmp_"+name_main+".fa"
           infile_hit   =dummy_dir+"/tmp_"+name_second+".fa"
           outfile      =dummy_dir+"/tmp_"+name_main+"_"+name_second+".aln"
           dndfile      =dummy_dir+"/tmp_"+name_main+"_"+name_second+".dnd"
           fd=open(infile,"w")
           fdm=open(infile_query,"w")
           fds=open(infile_hit,"w")
           fd.write(">{0:s}\n{1:s}\n".format(name_main,sequence_main))
           fdm.write(">{0:s}\n{1:s}\n".format(name_main,sequence_main))
           fd.write(">{0:s}\n{1:s}\n".format(name_second,sequence_second))
           fds.write(">{0:s}\n{1:s}\n".format(name_second,sequence_second))
           fd.close()
           fdm.close()
           fds.close()
           try:
             try:  
                # run clustalo
                msa_cline=Applications.ClustalwCommandline(clustal_exe,infile=infile,outfile=outfile)
                child = subprocess.Popen(str(msa_cline),stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell="/bin/bash")
                child.communicate()
                #store alignment in compare
                alignment=AlignIO.read(outfile,'clustal')
                main_aligned=alignment[0].seq
                second_aligned=alignment[1].seq
                alignment_length=alignment.get_alignment_length()
             except:
                 try:  
                    dummy_main_aligned, dummy_second_aligned, dummy_query_start, dummy_query_end, dummy_hit_start, dummy_hit_end, dummy_identity, dummy_similarity = exec_matcher(infile_query, infile_hit)
                    main_aligned=""
                    second_aligned=""
                    for i in range(dummy_query_start-1):
                        main_aligned=main_aligned+sequence_main[i]
                        second_aligned=second_aligned+"-"
                    for i in range(dummy_hit_start-1):
                        main_aligned=main_aligned+"-"
                        second_aligned=second_aligned+sequence_second[i]
                    for i in range(len(dummy_main_aligned)):
                        main_aligned=main_aligned+dummy_main_aligned[i]
                        second_aligned=second_aligned+dummy_second_aligned[i]
                    alignment_length=len(main_aligned)
                 except Exception as e:
                    sys.stdout.write("Error of MATCHER. Confirm EMBOSS package is correctly installed\n")
                    sys.stderr.write("ERROR: %s\n"%e)
                    raise e
             try:
               len_main =len(main_aligned)
               len_second=len(second_aligned)
             except Exception as e:
               sys.stderr.write("ERROR: %s\n"%e)
               raise e
           except Exception as e:
             sys.stderr.write("ERROR: %s\n"%e)
             raise e
           #remove temporary fasta and alignment files
           #remove_files([infile,infile_query,infile_hit,outfile,dndfile])
           #count identities
           identities = 0
           for i in range(alignment_length): 
               if main_aligned[i] == second_aligned[i]: 
                  identities = identities + 1
           score.setdefault(name_second,(float(identities)/shortest_length,shortest_length))
       maximum_id = max([sc[0] for name,sc in score.iteritems()]) 
       maximum_length=0
       for name,sc in score.iteritems():
         if sc[0]==maximum_id:
            if sc[1]>maximum_length:
               maximum_length=sc[1]
               selected_second=name 
       pairing.setdefault(name_main,selected_second)
       pairing.setdefault(selected_second,name_main)
       used.add(selected_second)

    return pairing


def chains_fixed(pdb):
    ''' 
     Rename the chains of PDB file located in path folder 
     A,B,C,D ... for proteins
     a,b,c,d ... for nucleid acids

     path	Folder where PDB file is located
     pdb 	PDB file
     dummy_dir  Dummy directory to cerate files

    '''
    #Initialize
    from SBI.structure.chain import Chain
    from SBI.structure.chain import ChainOfProtein
    from SBI.structure.chain import ChainOfNucleotide
    from SBI.structure import PDB

    protein_chains=list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    nucleic_chains=list("abcdefghijklmnopqrstuvwxyz")
    name_pdb=pdb.id
    new_pdb=PDB()
    protein_ids=[]
    nucleic_ids=[]
    print("\t\t\tFIX CHAINS")
    for chain_id in pdb.chain_identifiers:
        chain=pdb.get_chain_by_id(chain_id)
        print("\t\t\t\tCHAIN ID: %s TYPE %s "%(chain_id,chain.chaintype))
        if chain.chaintype == "N": 
            nucleic_ids.append(chain_id)
        else:
            protein_ids.append(chain_id)
    print("\t\t\t\tPROTEIN CHAINS %s"%str(protein_ids))
    print("\t\t\t\tNUCLEIC CHAINS %s"%str(nucleic_ids))
    for i in range(len(nucleic_ids)):
        chain_id = nucleic_ids[i]
        chain=pdb.get_chain_by_id(chain_id)
        n=int(float(i)/len(nucleic_chains))
        if n>0: new_id=nucleic_chains[i]+str(n)
        else:   new_id=nucleic_chains[i]
        chain.chain=new_id
        new_pdb.add_chain(chain)
    for i in range(len(protein_ids)):
        chain_id = protein_ids[i]
        chain=pdb.get_chain_by_id(chain_id)
        n=int(float(i)/len(protein_chains))
        if n>0: new_id=protein_chains[i]+str(n)
        else:   new_id=protein_chains[i]
        chain.chain=new_id
        new_pdb.add_chain(chain)

    return new_pdb

def renumber_pdb(pdb_file,sequences,dummy_dir="/tmp"):
    ''' 
     Renumber PDB file located in path folder with the real sequences

     path	Folder where PDB file is located
     pdb 	PDB file
     sequences  dictionary of sequence tuples (name,seq) by chain id
                chain identifier is the key of the dictionary
     dummy_dir  Dummy directory to cerate files

    '''

    #Initialize
    from SBI.structure.chain import Chain
    from SBI.structure.chain import ChainOfProtein
    from SBI.structure.chain import ChainOfNucleotide
    from SBI.sequence import Sequence
    from SBI.structure import PDB
    from Bio import SeqIO
    from Bio import ExPASy
    from Bio import AlignIO
    from Bio.Align import Applications

    clustal_exe  = os.path.join(config.get('Paths','clustalo_path'),'clustalo')
    new_pdb=PDB()
    name_pdb = os.path.basename(pdb_file)
    path     = os.path.dirname(pdb_file)
    pdb=PDB(pdb_file)
    pdb.clean()
    for chain_id,chain_seq in sequences.iteritems():
       name_chain = name_pdb+"_"+chain_id
       pdb_chain  = pdb.get_chain_by_id(chain_id)
       if chain_seq[0] is None or chain_seq[1] is None or pdb_chain.chaintype != "P":
         if pdb_chain.chaintype != "P": 
            print("\t\t\t-- Add NON-PROTEIN chain "+chain_id)
            new_pdb.add_chain(pdb_chain)
         continue
    for chain_id,chain_seq in sequences.iteritems():
       name_chain = name_pdb+"_"+chain_id
       pdb_chain  = pdb.get_chain_by_id(chain_id)
       if chain_seq[0] is None or chain_seq[1] is None or pdb_chain.chaintype != "P": continue
       name_seq   = chain_seq[0]
       pdb_chain  = pdb.get_chain_by_id(chain_id)
       new_chain  = ChainOfProtein(name_pdb,chain_id)
       #define/create files
       infile       =dummy_dir+"/tmp_"+name_chain+"_"+name_seq+".fa"
       infile_query =dummy_dir+"/tmp_"+name_chain+".fa"
       infile_hit   =dummy_dir+"/tmp_"+name_seq+".fa"
       outfile      =dummy_dir+"/tmp_"+name_chain+"_"+name_seq+".aln"
       dndfile      =dummy_dir+"/tmp_"+name_chain+"_"+name_seq+".dnd"
       fd=open(infile,"w")
       fdm=open(infile_query,"w")
       fds=open(infile_hit,"w")
       fd.write(">{0:s}\n{1:s}\n".format(name_chain,pdb_chain.protein_sequence))
       fdm.write(">{0:s}\n{1:s}\n".format(name_chain,pdb_chain.protein_sequence))
       fd.write(">{0:s}\n{1:s}\n".format(name_seq,chain_seq[1]))
       fds.write(">{0:s}\n{1:s}\n".format(name_seq,chain_seq[1]))
       fd.close()
       fdm.close()
       fds.close()
       try:
         try:
             # run clustalw2
             msa_cline=Applications.ClustalwCommandline(clustal_exe,infile=infile,outfile=outfile)
             child = subprocess.Popen(str(msa_cline),stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell="/bin/bash")
             child.communicate()
             #store alignment in compare
             alignment=AlignIO.read(outfile,'clustal')
             structure=alignment[0].seq
             reference=alignment[1].seq
         except:
             try:  
                dummy_structure, dummy_reference, dummy_query_start, dummy_query_end, dummy_hit_start, dummy_hit_end, dummy_identity, dummy_similarity = exec_matcher(infile_query, infile_hit)
                structure=""
                reference=""
                for i in range(dummy_query_start-1):
                    structure=structure+pdb_chain.protein_sequence[i]
                    reference=reference+"-"
                for i in range(dummy_hit_start-1):
                    structure=structure+"-"
                    reference=reference+chain_seq[1][i]
                for i in range(len(dummy_structure)):
                    structure=structure+dummy_structure[i]
                    reference=reference+dummy_reference[i]
                alignment_length=len(structure)
             except Exception as e:
                sys.stdout.write("Error of MATCHER. Confirm EMBOSS package is correctly installed\n")
                sys.stderr.write("ERROR: %s\n"%e)
                return e
         try:
           len_3d =len(structure)
           len_ref=len(reference)
         except Exception as e:
           sys.stderr.write("ERROR: %s\n"%e)
           return e
       except Exception as e:
         sys.stderr.write("ERROR: %s\n"%e)
         return e
       #print("!======================== ")
       #print("ALIGNMENT")
       #print("STRUCTURE "+structure)
       #print("SEQUENCE  "+reference)
       #print("=========================! ")
       #remove temporary fasta and alignment files
       remove_files([infile,infile_query,infile_hit,outfile,dndfile])
       #mapping of residues to the original sequence
       mapping=create_mapping(pdb_chain.protein_idx.split(";"),structure,reference)
       #print "MAPPING ", mapping
       #fill the new chain with the correct numbering of residues
       for residue in pdb_chain.aminoacids:
          pair = (str(residue.number),residue.version)
          if not mapping.has_key(pair): continue
          number,version = mapping.get(pair)
          residue.number=number
          residue.version=version
          new_chain.add_residue(residue)
       #fill the new pdb
       new_chain.chaintype=pdb_chain.chaintype
       print("\t\t\t-- Add PROTEIN chain "+chain_id)
       new_pdb.add_chain(new_chain)

    return new_pdb  

def create_mapping(idx,structure,reference):
    mapped={}
    n=0
    m=0
    jump=0
    abc=" ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
    for i in range(len(reference)):
      idx_reference=m+1
      if structure[i]!="-":
         number =idx[n][0:-1]
         version=idx[n][-1]
         if reference[i] == "-":
           if jump < 53: new_version = abc[jump]
           else: new_version="@"
         else:
           new_version = version
         mapped.setdefault((number,version),(idx_reference,new_version))
         n = n +1 
      if reference[i]!="-":
         jump = 0
         m    = m + 1
      else:
         jump = jump + 1
        
    return mapped 




def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python model_protein.py -i input_file -p pdb_dir [--dummy=dummy_dir --n-model=n_model --n-total=n_total -o output_dir] [-a -d -f -m -e -r resolution_file -s --dimer --monomer --unbound_fragments --unrestrictive] [-t] ")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file in FASTA format (sequence) or THREADING format (i.e. from threader.py)", metavar="{filename}")
    parser.add_option("--n-model", default=1, action="store", type="int", dest="n_model", help="Number of models per template (default = 1)", metavar="{int}")
    parser.add_option("--n-total", action="store", type="int", dest="n_total", default=None, help="Total number of models per execution of this program (if not \"None\", automatically sets the \"--n-model\" parameter; supersedes option \"--n-model\"; default = None)", metavar="{int}")
    parser.add_option("--best", action="store_true", dest="best", default=False, help="It only uses the first hit of blast (supersedes n-model and n-total to 1; default = False)", metavar="{boolean}")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="{directory}")
    parser.add_option("-p", "--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (i.e. output dir from pdb.py)", metavar="{directory}")
    parser.add_option("-l", action="store", type="string", dest="label", help="Label to include in the output models name", metavar="{str}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode. If not selected the dummy directory will be removed (default = False)", metavar="{boolean}")
    parser.add_option("--force", default=False, action="store_true", dest="force", help="Force to do the modelling even if the file already exists (default = False)", metavar="{boolean}")
    parser.add_option("--opt", default=False, action="store_true", dest="optimization", help="Run modelling with optimization (default = False). Note: use it with care, this is not advisable as structures can take missconfomations", metavar="{boolean}")
    parser.add_option("--info",default=None,action="store", type="string", dest="info",help="Information LOG file of MODELS that have failed and have been completed")
 

    group = optparse.OptionGroup(parser, "FASTA options", "Provide a FASTA file as input.")
    group.add_option("-a", "--all", default=False, action="store_true", dest="all_models", help="Use ALL possible templates (default = False)", metavar="{boolean}")
    group.add_option("-d", "--dna", default=False, action="store_true", dest="dna", help="DNA modelling (if \"True\", remodels DNA using 3DNA; default = False)", metavar="{boolean}")
    group.add_option("-f","--full", default=False, action="store_true", dest="full_mode", help="Model all protein fragments with template that bind DNA (default = False)", metavar="{boolean}")
    group.add_option("--unbound_fragments", default=False, action="store_true", dest="full_fragments_mode", help="Model all protein fragments with template even unbound to DNA (default = False)", metavar="{boolean}")
    group.add_option("-e", "--twilight", default=False, action="store_true", dest="filter_twilight_zone_hits", help="Select template-hits over the twilight zone unless there are no models (default = False)", metavar="{boolean}")
    group.add_option("--restrictive", default=False, action="store_true", dest="filter_restrictive", help="Restrict modelling to only template-hits over the twilight zone when option --twilight is selected (default = False)", metavar="{boolean}")
    group.add_option("--unrestrictive", default=False, action="store_true", dest="no_restrict", help="Allow models with wrong alignment in the DNA-binding interface (default = False)", metavar="{boolean}")
    group.add_option("-m", action="store", dest="mutation", help="Mutate input sequence (e.g. \"R88C\" mutates the \"R\" at position \"88\" for a \"C\")", metavar="{str}")
    group.add_option("-r", action="store", default=None, dest="resolution_file", help="Resolution file (if provided, sorts similar alignments by RMSD of the template PDB; default = None)", metavar="{filename}")
    group.add_option("-s","--select", default=False, action="store_true", dest="select_mode", help="Select mode (if True, selects either monomers or dimers, whichever has more hits; default = False)", metavar="{boolean}")
    group.add_option("--dimers", default=False, action="store_true", dest="dimer_mode", help="If True forces to work only with dimers; default = False)", metavar="{boolean}")
    group.add_option("--monomers", default=False, action="store_true", dest="monomer_mode", help="If True forces to work only with monomers; default = False)", metavar="{boolean}")
    group.add_option("--renumerate", default=False, dest = 'renumerate' , action = 'store_true', help = 'Flag to renumber the sequences as in the original FastA input (default is False)')
    group.add_option("--chains_fixed",default=False, action="store_true", dest="chains_fixed", help="Replace the names of the protein chains to A-B-C-D... and DNA chains to a-b-c-d....")


    parser.add_option_group(group)

    group = optparse.OptionGroup(parser, "Threading options", "Provide a threading file as input (i.e from threader.py).")
    group.add_option("-t", "--threading", default=False, action="store_true", dest="threading", help="Threading mode (default = False)")
    parser.add_option_group(group)

    (options, args) = parser.parse_args()

    if options.input_file is None  or options.pdb_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")
    if options.mutation is not None:
        m = re.search("^([a-zA-Z]\d+[a-zA-Z])$", options.mutation)
        if not m:
            parser.error("incorrect provided mutation \"%s\": e.g.\"R88C\" mutates the \"R\" at position 88 for a \"C\"" % options.mutation)

    return options

def get_pir_alignments(query_file, pdb_dir, full_mode=False,  filter_twilight_zone_hits=False, resolutions=None, select_mode=False, dummy_dir="/tmp/", original_sequence=None, 
                       thread_obj=None, no_restrict_usage=False, dimer_mode=False, monomer_mode=False, verbose=False, filter_restrictive=False, best_mode=False, iteration=0, start_time=0):
    """
    This function blasts a query sequence against the PDB database and 
    returns a list of {PirAlignment}s, according to some criteria.

    @input:
    query_file {filename}
    pdb_dir {directory}
    filter_twilight_zone_hits {boolean} default = False
    resolutions {dict} default = None
    select_mode {boolean} selects either monomers or dimers, whichever has more hits; default = False
    dimer_mode {boolean} only dimer templates are allowed
    monomer_mode {boolean} only monomer templates are allowed
    dummy_dir {directory} default = /tmp
    verbose {boolean} to print output steps of the work
    filter_restrictive {boolean} default = False to be retsrictive on the twilight_zone_hits criterion
    @return:
    selected_pir_alignments / pir_alignments {list} of {PirAlignment}

    """

    # Initialize #
    hits = []
    pir_alignments = []
    sequences = []
    headers = []
    coords = []
    monomer_hits = set()
    dimer_hits = set()
    done_hits = set()
    min_sequence_length = int(config.get("Parameters", "min_aminoacids"))
    no_restrict = False
    no_restrict_usage_original = no_restrict_usage
    if no_restrict_usage: no_restrict = True
    filter_twilight_zone_hits_original=filter_twilight_zone_hits
    domain_dict = {}
    mod_regions = []
    pdb_inuse = pdb_dir.split("/")[-1]
    start_point = 0
   
    # Check that modelling regions are not exceeding the time of calculation or the number of iterations (   
    current_time = float(time.time())
    iteration    = iteration+1
    if verbose:
       sys.stdout.write("\t-- Iteration to GET ALIGNMENTS: %d (time %.3f h)\n"%(iteration,(current_time- start_time)/3600))
    if iteration >1000 or (current_time - start_time) > 43200:
       sys.stdout.write("Exceeding time (12 hours) or iterations (1000): %.3f hours %d iterations\n"%((current_time-start_time)/3600,iteration))
       sys.stdout.write("STOP new searches of homologs\n")
       return pir_alignments

    if thread_obj == None:
        # For each FASTA sequence... #
        for header, sequence in functions.parse_fasta_file(query_file):
            # Add sequence #
            sequences.append(sequence)
            headers.append(header)
            # Get header in case the original sequence is not None, so we are processing a non modeled region, save the starting point #
            if original_sequence != None:
               for original in original_sequence:
                 header_original,sequence_original=original
                 if header == header_original or len(headers)==1:
                   start_point = int(header.split("_")[1].split(":")[0])
            # Create dummy file #
            dummy_file = os.path.join(dummy_dir, "%s.query.fa" % header+"_"+str(os.getpid()))
            if query_file != dummy_file:
                if fileExist(dummy_file): remove_files([dummy_file])
                functions.write(dummy_file, ">%s\n%s" % (header, sequence))
            # Get BLAST object #
            if verbose: sys.stdout.write("\t-- Blast homologs search of " + header + "... \n")
            blast_obj = blast.get_blast_obj(os.path.join(pdb_dir, "database", "database.fasta"), dummy_file, dummy_dir=dummy_dir)
            # Get unfiltered hits #
            size_new_hits=0
            try:
               accepted_hits = filter_blast_hits_by_interface(dummy_file, blast_obj, pdb_dir, filter_twilight_zone_hits, dummy_dir, no_restrict,best_mode)
               size_new_hits = len(accepted_hits)
               hits.append(accepted_hits)
            except Exception as e:
               sys.stdout.write("ERROR when filtering BLAST hits\n")
               raise e
            if hits[0] != []:
                # if you have found nice hits by using restrictions, then you won't use the no_restrict option to find hits in the non modeled regions #
                no_restrict_usage = True
            if verbose: 
                sys.stdout.write("\t\t-- Found %d hits ...\n"%(size_new_hits))
            # Remove dummy file #
            if os.path.exists(dummy_file)  and not verbose: os.remove(dummy_file)
            # Skip if already did 2 sequences #
            if len(hits) >= 2: break
        #If input is a dimer check if there are results for both sequences such that we can create the dimer model
        dimer_full_test=False
        if len(hits)>1: 
           dimer_full_test=(len(headers)>len(hits) or hits[1]==[])
           if not dimer_full_test:
               # Parse dimers file #
               dimers = {}
               dimer_list = []
               # For each line... #
               for line in functions.parse_file(os.path.join(pdb_dir, "dimers.txt")):
                   if line.startswith("#"): continue
                   line = line.split(";")
                   dimers[line[0]] = line[1]
                   dimers[line[1]] = line[0]
                   dimer_list.append(line[0])
                   dimer_list.append(line[1])
               # For each hit... #
               dimer = None
               for hit in hits[0]:
                   if dimer is not None:continue
                   # If TF dimerizes... #
                   if str(hit[0]) in dimer_list:
                       if dimer is not None:continue
                       # For each hit... 
                       for next_hit in hits[1]:
                           if next_hit[0] == dimers[hit[0]]:
                               dimer = next_hit
                               break
               # Skip if found nothing #
               if dimer is None: dimer_full_test=True
        # In case we don't obtain any result, we repeat the execution of blast but changing the filter parameter #
        if (hits[0] == [] or dimer_full_test ) and filter_twilight_zone_hits == True and not filter_restrictive:
            filter_twilight_zone_hits = False
            no_restrict_usage = False
            headers=[]
            sequences=[]
            hits = []
            for header, sequence in functions.parse_fasta_file(query_file):
                if verbose: sys.stdout.write("\t\t-- Search remote homologs of " + header + " with low similarity ... \n")
                # Add sequence #
                sequences.append(sequence)
                headers.append(header)
                # Create dummy file #
                dummy_file = os.path.join(dummy_dir, "query.%s.fa" % header+"_"+str(os.getpid()))
                if query_file != dummy_file:
                    if fileExist(dummy_file): remove_files([dummy_file])
                    functions.write(dummy_file, ">%s\n%s" % (header, sequence))
                # Get BLAST object #
                blast_obj = blast.get_blast_obj(os.path.join(pdb_dir, "database", "database.fasta"), dummy_file, dummy_dir=dummy_dir)
                # Get unfiltered hits #
                size_new_hits=0
                try:
                  accepted_hits = filter_blast_hits_by_interface(dummy_file, blast_obj, pdb_dir, filter_twilight_zone_hits, dummy_dir, no_restrict,best_mode)
                  size_new_hits = len(accepted_hits)
                  hits.append(accepted_hits)
                except Exception as e:
                  sys.stdout.write("ERROR when filtering BLAST hits\n")
                  raise e
                if hits[0] != []:
                   # if you have found nice hits by using restrictions, then you won't use the no_restrict option to find hits in the non modeled regions #
                   no_restrict_usage = True
                if verbose: 
                   sys.stdout.write("\t\t\t-- Found %d hits  ...\n"%(size_new_hits))
                # Remove dummy file #
                if os.path.exists(dummy_file) and not verbose: os.remove(dummy_file)
                # Skip if already did 2 sequences #
                if len(hits) >= 2: break
        sys.stdout.flush()
        #If input is a dimer check if there are results for both sequences such that we can create the dimer model
        if verbose: 
           sys.stdout.write("\t-- Check dimers\n")
        dimer_full_test=False
        if len(hits)>1: 
           dimer_full_test=(len(headers)>len(hits) or hits[1]==[])
           if not dimer_full_test:
               # Parse dimers file #
               dimers = {}
               dimer_list = []
               # For each line... #
               for line in functions.parse_file(os.path.join(pdb_dir, "dimers.txt")):
                   if line.startswith("#"): continue
                   line = line.split(";")
                   dimers[line[0]] = line[1]
                   dimers[line[1]] = line[0]
                   dimer_list.append(line[0])
                   dimer_list.append(line[1])
               # For each hit... #
               dimer = None
               for hit in hits[0]:
                   if dimer is not None:continue
                   # If TF dimerizes... #
                   if str(hit[0]) in dimer_list:
                       if dimer is not None:continue
                       # For each hit... #
                       for next_hit in hits[1]:
                           if next_hit[0] == dimers[hit[0]]:
                               dimer = next_hit
                               break
               # Skip if found nothing #
               if dimer is None: dimer_full_test=True
        # In case we don't obtain any result, we repeat the execution of blast but changing the filter parameter and we now allow low similarity in the interface sequence and secondary structure #
        if (hits[0] == [] or dimer_full_test ) and (filter_twilight_zone_hits == False) and (no_restrict == False) and (no_restrict_usage == False  or dimer_full_test ) and not filter_restrictive:
            no_restrict = True
            headers=[]
            sequences=[]
            hits = []
            for header, sequence in functions.parse_fasta_file(query_file):
                if verbose: sys.stdout.write("\t\t-- Search remote homologs of " + header + " without restrictions on the interface....\n")
                # Add sequence #
                sequences.append(sequence)
                headers.append(header)
                # Create dummy file #
                dummy_file = os.path.join(dummy_dir, "query.%s.fa" % header+"_"+str(os.getpid()))
                if query_file != dummy_file:
                    if fileExist(dummy_file): remove_files([dummy_file])
                    functions.write(dummy_file, ">%s\n%s" % (header, sequence))
                # Get BLAST object #
                blast_obj = blast.get_blast_obj(os.path.join(pdb_dir, "database", "database.fasta"), dummy_file,  dummy_dir=dummy_dir)
                # Get unfiltered hits #
                size_new_hits=0
                try:
                   accepted_hits = filter_blast_hits_by_interface(dummy_file, blast_obj, pdb_dir, False, dummy_dir, no_restrict,best_mode)
                   size_new_hits = len(accepted_hits)
                   hits.append(accepted_hits)
                except Exception as e:
                   sys.stdout.write("ERROR when filtering BLAST hits\n")
                   raise e
                if verbose: 
                   sys.stdout.write("\t\t\t-- Found %d hits  ...\n"%(size_new_hits))
                # Remove dummy file #
                if os.path.exists(dummy_file)  and not verbose: os.remove(dummy_file)
                # Skip if already did 2 sequences #
                if len(hits) >= 2: break

        # If no hits are found we get out of this function #
        dimer_full_test=False
        if len(hits)>1: dimer_full_test=(len(headers)>len(hits) or hits[1]==[])
        if (hits[0] == [] or dimer_full_test ) and filter_twilight_zone_hits == False:
            if verbose: sys.stdout.write("\t\t-- Skip modeling " + str(headers) +". No matches found....\n")
            return None, None

    else:
        # the hits object is going to be only the threaded pdb #
        hits.append([thread_obj])
        # Actually, when we use threading the query file is not a file, but directly the query sequence #
        # sequences.append(query_file)
        sequences.append(thread_obj[2].replace("-",""))
    # If resolutions are provided #
    if resolutions is not None:
        # Initialize #
        hits_copy = []
        # For each list of hits... #
        for hit_list in hits:
            # Initialize #
            hits_copy.append([])
            sorted_hits = []
            # For each hit... #
            for hit in hit_list:
                sorted_hits.append([hit, (get_identities(hit[1], hit[2]) * 100.0) / get_aligned_residues(hit[1], hit[2]), "RMSD"])
                if hit[0][:4] in resolutions:
                    sorted_hits[-1][-1] = resolutions[hit[0][:4]]
            # Sort #
            sorted_hits.sort(key=lambda x: (-x[1], x[2]))
            # For each hit... #
            for hit in sorted_hits:
                # Copy hit to hits #
                hits_copy[-1].append(hit[0])
        # Copy hits #
        hits = hits_copy
    # Parse dimers file #
    dimers = {}
    dimer_list = []
    if thread_obj is None:
      # For each line... #
      for line in functions.parse_file(os.path.join(pdb_dir, "dimers.txt")):
        if line.startswith("#"): continue
        line = line.split(";")
        dimers[line[0]] = line[1]
        dimers[line[1]] = line[0]
        dimer_list.append(line[0])
        dimer_list.append(line[1])
    # For each hit... #
    for hit in hits[0]:
        # Skip if already done #
        if hit[0] in done_hits: continue
        # Initialize #
        pir_alignment_obj = PirAlignment()
        done_hits.add(hit[0])
        # Get the boundaries of the modeled sequences #
        coords.append(hit[3])
        coords.append(hit[4])
        # If TF dimerizes... #
        if str(hit[0]) in dimer_list:
            # Initialize #
            i = 0
            if len(hits) >= 2: i = 1
            dimer = None
            # For each hit... #
            for next_hit in hits[i]:
                if next_hit[0] == dimers[hit[0]]:
                    dimer = next_hit
                    break
            # Skip if found nothing #
            #if dimer is None or monomer_mode:
            if dimer is None:
                # If hit sequence identity is above the twilight zone and restrictions over the interface have been considered then we model the monomer alone #  
                if (bool(hit[-4]) == True and bool(hit[-3]) == False): 
                    pir_alignment_obj.set_heterodimer(dimers[hit[0]])
                # To allow the modeling of the monomer we remove the hit from the dimers dictionary #
                if verbose: sys.stdout.write("\t\t\t-- Template %s cannot be used for dimers\n"%(str(hit[0])))
                dimers.pop(str(hit[0]), None)
            else:
                # Add dimer unless we use all monomers#
                if not monomer_mode:
                   done_hits.add(dimer[0])
        # Add hit/sequence/structure name #
        pir_alignment_obj.add_hit_name(hit[0])
        pir_alignment_obj.add_sequence_name("query")
        pir_alignment_obj.add_structure_name(hit[0][:4])
        pir_alignment_obj.set_identity(hit[-2])
        pir_alignment_obj.set_similarity(hit[-1])
        # Get PDB object #
        pdb_obj = PDB(os.path.join(pdb_dir, "split",  "%s.pdb" % hit[0]))
        # If TF dimerizes... #
        if hit[0] in dimers:
            dimer_hits.add(hit[0])
            dimer_hits.add(dimer[0])
            if not monomer_mode:
               protein_pdb_obj = PDB(os.path.join(pdb_dir, "split", dimer[0] + ".pdb"))
               pdb_obj.add_chain(protein_pdb_obj.chains[0])
        else:
            monomer_hits.add(hit[0])
        # For each helix... #
        for line in functions.parse_file(os.path.join(pdb_dir, "helices", hit[0] + ".txt")):
            dna_pdb_obj = PDB(os.path.join(pdb_dir, "split", hit[0][:4] + ".dna." + line + ".pdb"))
            for dna_pdb_chain_obj in dna_pdb_obj.chains:
                pdb_obj.add_chain(dna_pdb_chain_obj)
            break
        # For each PDB chain... #
        for pdb_chain_obj in sorted(pdb_obj.chains, key=lambda x: x.chain):
            # Add structure chain #
            pir_alignment_obj.add_structure_chain(pdb_chain_obj.chain)
            # Add structure residues #
            pir_alignment_obj.add_structure_residue(pdb_chain_obj._all_residues[0].number)
            pir_alignment_obj.add_structure_residue(pdb_chain_obj._all_residues[-1].number)
            # Add alignments #
            add_alignment = False
            if pdb_chain_obj.chaintype == "P":
                # If hit PDB chain... #
                if pdb_chain_obj.chain == hit[0][-1]:
                    sequence_alignment, structure_alignment = readjust_alignments_to_crystal_sequence(pdb_chain_obj, hit[1], hit[2], hit[5], hit[6])
                    add_alignment = True
                    # Initialize #
                    sub_sequence = re.sub("\W", "", sequence_alignment)
                    # Get sub-sequence coordinates #
                    if original_sequence is None: start = sequences[0].find(sub_sequence)
                    else:  start = original_sequence[0][1].find(sub_sequence)
                    # Add sequence coordinates to PIR alignment object #
                    pir_alignment_obj.add_sequence_residue(start + 1)
                    pir_alignment_obj.add_sequence_residue(start + len(sub_sequence))
                    pir_alignment_obj.add_main_order(1)
                    pir_alignment_obj.add_main_order(1)
                # If dimer PDB chain... #
                else:
                    sequence_alignment, structure_alignment = readjust_alignments_to_crystal_sequence(pdb_chain_obj, dimer[1], dimer[2], dimer[5], dimer[6])
                    add_alignment = True
                    # Initialize #
                    sub_sequence = re.sub("\W", "", sequence_alignment)
                    # Get sub-sequence coordinates #
                    if original_sequence is None:
                        if len(sequences) == 1: start = sequences[0].find(sub_sequence)
                        else: start = sequences[1].find(sub_sequence)
                    else:  
                        if len(original_sequence) == 1: start = original_sequence[0][1].find(sub_sequence)
                        else: start = original_sequence[1][1].find(sub_sequence)
                    # Add sequence coordinates to PIR alignment object #
                    pir_alignment_obj.add_sequence_residue(start + 1)
                    pir_alignment_obj.add_sequence_residue(start + len(sub_sequence))
                    pir_alignment_obj.add_main_order(0)
                    pir_alignment_obj.add_main_order(0)
                # If dimer PDB chain... #
                # Add alignments to PIR alignment object #
                pir_alignment_obj.add_sequence_alignment(sequence_alignment)
                pir_alignment_obj.add_structure_alignment(structure_alignment)
            else:
                #Use flexible DNA l=DG e=DA t=DT j=DC
                #dna_seq=pdb_chain_obj.nucleotide_sequence().replace("G","l").replace("A","e").replace("T","t").replace("C","j")
                #pir_alignment_obj.add_sequence_alignment(dna_seq)
                #pir_alignment_obj.add_structure_alignment(dna_seq)
                pir_alignment_obj.add_sequence_alignment("." * len(pdb_chain_obj.nucleotides))
                pir_alignment_obj.add_structure_alignment("." * len(pdb_chain_obj.nucleotides))
        # Add PIR alignment to PIR alignments #
        if pir_alignment_obj != None:
            pir_alignments.append(pir_alignment_obj)
        #if full_mode and len(hits) == 1 :
        if full_mode:
            # Processing the modeled coordinates # 
            nt = int(hit[3])
            ct = int(hit[4])
            if len(mod_regions) == 0:
                mod_regions.append([nt, ct])
            elif len(mod_regions) >= 1:
                new = True
                for reg in mod_regions:
                    # check if the two regions overlap #
                    p_set = set(range(nt, ct))
                    reg_set = set(range(reg[0], reg[1]))
                    if len(p_set.intersection(reg_set)) >= len(p_set)/2:
                        # if they overlap we extend the limits of the alignment. The overlap must consider at least half of the protein sequence #
                        reg[0] = min([nt, reg[0]])
                        reg[1] = max([ct, reg[1]])
                        
    # identify the different modeled regions # 
    try:
        non_mod_regions = []
        if len(mod_regions) == 1:
            non_mod_regions.append([0, mod_regions[0][0] - 1])
            non_mod_regions.append([mod_regions[0][1] + 1, len(sequences[0])])
            for p_obj in pir_alignments:
                if original_sequence == None:
                    p_obj.add_domain(str(mod_regions[0][0]) + ":" + str(mod_regions[0][1]))
                else:
                    p_obj.add_domain(str(mod_regions[0][0] + start_point) + ":" + str(mod_regions[0][1] + start_point))
        
    except:
        if verbose: sys.stdout.write("\t\t\t-- No models for the region: " + str(mod_regions) + "\n")


    sys.stdout.flush()
    if (full_mode) :
        # Extract the unmodeled sequences and run blast on them #
        for reg in non_mod_regions:
            # Initialize #
            #if original_sequence is None: original_sequence = "".join(sequences[0])
            if original_sequence is None: 
               original_headers=headers
               original_sequences=sequences
               original_sequence=[]
               for index_header in range(len(sequences)):
                  original_sequence.append((headers[index_header],sequences[index_header]))
            else:
               original_sequences=[]
               original_headers=[]
               for original in original_sequence:
                 header_original,sequence_original=original
                 original_headers.append(header_original)
                 original_sequences.append(sequence_original)
            # Get N-terminal end #
            tail = sequences[0][reg[0]:reg[1]]
            # Skip if N-terminal end is too short #
            if len(tail) < min_sequence_length: continue
            # Create dummy file #
            dummy_name=os.path.basename(query_file).split(".")[0].split("_")[0] + "_" + str(reg[0]) + ":" + str(reg[1])+"_"+str(os.getpid())
            dummy_file = os.path.join(dummy_dir, "%s.tails.fa" %dummy_name)
            if os.path.exists(dummy_file):
                os.system("rm -f " + dummy_file)
            # Write FASTA file #
            if len(hits)>1:
               functions.write(dummy_file, ">%s\n%s\n>%s\n%s" % (os.path.basename(query_file).split(".")[0].split("_")[0] + "_" + str(reg[0]) + ":" + str(reg[1]), tail,original_headers[1]+ "_" + str(0) + ":" + str(len(original_sequences[1])), original_sequences[1]))
            else:
               functions.write(dummy_file, ">%s\n%s" % (os.path.basename(query_file).split(".")[0].split("_")[0] + "_" + str(reg[0]) + ":" + str(reg[1]), tail))
            # Get list of PIR alignments #
            if verbose: sys.stdout.write("\t-- Redo the search of hits in the region " + str(reg[0]) + ":" + str(reg[1]) +  " \n")
            try:
               pir_alignments_list = get_pir_alignments(dummy_file, pdb_dir, full_mode, filter_twilight_zone_hits_original, resolutions, select_mode, dummy_dir, original_sequence, \
                                                        no_restrict_usage=no_restrict_usage_original, dimer_mode=dimer_mode, monomer_mode=monomer_mode, verbose=verbose, \
                                                        filter_restrictive=filter_restrictive,best_mode=best_mode, iteration=iteration, start_time=start_time)
            except Exception as e:
               sys.stdout.write("ERROR when getting PIR alignments\n")
               raise e
            # For each  PIR alignment list... #
            sys.stdout.flush()
            if pir_alignments_list != None:
                for i in range(len(pir_alignments_list)):
                    # For each PIR alignment... #
                    if pir_alignments_list[i] != None:
                        pir_alignments.append(pir_alignments_list[i])
    
    # If selection mode... #
    if dimer_mode or monomer_mode or select_mode:
        # Initialize #
        selected_pir_alignments = []
        # For each list of PIR alignments... #
        for pir_alignment_obj  in pir_alignments:
            # For each PIR alignment object... #
                # Initalize #
                dimer = False
                # For each PDB chain... #
                for pdb_chain in pir_alignment_obj.get_structure_chains():
                    if dimer: continue
                    if pir_alignment_obj.get_structure_name() + "_" + pdb_chain in dimers:
                        dimer = True
                if select_mode:
                  # If as many dimers as there are monomers... #
                  if len(dimer_hits) == len(monomer_hits):
                    #selected_pir_alignments[-1].append(pir_alignment_obj)
                    selected_pir_alignments.append(pir_alignment_obj)
                  # If dimer and more dimers than monomers... #
                  elif dimer and len(dimer_hits) > len(monomer_hits):
                    #selected_pir_alignments[-1].append(pir_alignment_obj)
                    selected_pir_alignments.append(pir_alignment_obj)
                  # If monomer and more monomers than dimers... #
                  elif not dimer and len(monomer_hits) > len(dimer_hits):
                    selected_pir_alignments.append(pir_alignment_obj)
                else:
                  if dimer_mode and dimer:
                    selected_pir_alignments.append(pir_alignment_obj)
                  #uncomment next line if you want to restrict to only monomers  
                  #if monomer_mode and not dimer: 
                  #comment next line if you want to restrict to only monomers  
                  if monomer_mode:
                    selected_pir_alignments.append(pir_alignment_obj)
        if select_mode and (len(dimer_hits)>len(monomer_hits) or len(dimer_hits)==len(monomer_hits)) and verbose: sys.stdout.write("\t-- Found dimers %d\n"%(len(dimer_hits)))
        if select_mode and (len(dimer_hits)<len(monomer_hits) or len(dimer_hits)==len(monomer_hits)) and verbose: sys.stdout.write("\t-- Found monomers %d\n"%(len(monomer_hits)))
        if verbose: sys.stdout.write("\t-- Found dimers %d\n"%(len(dimer_hits)))
        if verbose: sys.stdout.write("\t-- Found monomers %d\n"%(len(monomer_hits)))
        return selected_pir_alignments
    else:
        return pir_alignments


def filter_blast_hits_by_interface(query_file, blast_obj, pdb_dir, filter_twilight_zone_hits=False, dummy_dir="/tmp/", no_restrict=False,best_mode=False):
    """
    This function filters hits if the query does not cover all interface
    interface residues or the alignments contain insertions or deletions in
    structured regions containing any interface residues.

    @input:
    query_file {filename}
    blast_obj {BlastOutput}
    pdb_dir {directory}
    filter_twilight_zone_hits {boolean} default = False
    dummy_dir {directory} default = /tmp

    @return:
    unfiltered_hits {list}
    
    """
    
    # Initialize #
    unfiltered_hits = []
    tz_parameter = 0
    tz_type = None
    if filter_twilight_zone_hits:
        tz_parameter = int(config.get("Parameters", "twilight_zone_parameter"))
        tz_type = config.get("Parameters", "twilight_zone_type")
    # For each hit object... #
    for hit_obj in blast_obj.get_hits(tz_parameter=tz_parameter, tz_type=tz_type):
        #print hit_obj.sequenceID,hit_obj.e_value
        # Skip if already done hit #
        done = False
        for unfiltered_hit in unfiltered_hits:
            if hit_obj.sequenceID in unfiltered_hit: done = True
        if done: continue
        hit_file = os.path.join(pdb_dir, "split", hit_obj.sequenceID + ".fasta")
        # Skip if  e-value is larger than a certain threshold #
        e_thres = float(config.get("Parameters", "e-value_threshold"))
        if (float(str(hit_obj).split("\t")[5]) > float(e_thres)) :
            sys.stdout.write("\t\t-- Skip alignment %s %s : limitted by e-value threshold %.4f > %.4f \n"%(os.path.basename(hit_file),os.path.basename(query_file),float(str(hit_obj).split("\t")[5]),float(e_thres)))
            sys.stdout.flush()
            continue
        # Refine alignment #   
        try:
           sys.stdout.write("\t-- Get alignment with matcher %s %s \n"%(hit_file,query_file))
           query_alignment, hit_alignment, query_start, query_end, hit_start, hit_end, identity, similarity = exec_matcher(query_file, hit_file)
           if query_alignment is None:
              sys.stdout.write("\t-- Skip alignment %s %s \n"%(os.path.basename(hit_file),os.path.basename(query_file)))
              continue
           else:
              sys.stdout.write("\t-- Aligned %s %s \n"%(os.path.basename(hit_file),os.path.basename(query_file)))
           sys.stdout.flush()
        except Exception as e:
           #if query_alignment is None or hit_alignment is None:
           sys.stdout.write("Error of MATCHER. Confirm EMBOSS package is correctly installed\n")
           sys.stdout.flush()
           raise  e
        # Get hit PDB complex with DNA #
        protein_pdb_obj = PDB(os.path.join(pdb_dir, "split", hit_obj.sequenceID + ".pdb"))
        # Get sequence/PDB correlation #
        sequence_to_crystal, crystal_to_sequence = get_sequence_to_crystal_correlations(protein_pdb_obj, protein_pdb_obj.chains[0].chain)
        # Get crystal sequence #
        crystal_sequence = protein_pdb_obj.chains[0].gapped_protein_sequence
        # Get interface residues #
        interface_residues = set()
        for contact_obj in contacts.Contacts(os.path.join(pdb_dir, "contacts", hit_obj.sequenceID[:4] + ".txt")).get_contacts():
            if contact_obj._A_chain == protein_pdb_obj.chains[0].chain:
                interface_residues.add(contact_obj._A_residue_obj.number)
        # Get secondary structure #
        secondary_structure = ""
        dssp_obj = dssp.DSSP(os.path.join(pdb_dir, "dssp", hit_obj.sequenceID[:4] + ".txt"))
        for i in range(len(crystal_sequence)):
            if crystal_sequence[i] == "x":
              secondary_structure += "C"
            else:
              try:
                j = dssp_obj.get_secondary_structure(protein_pdb_obj.chains[0].chain, sequence_to_crystal[i])
                if j == "E" or j == "H":
                    if i + 1 in interface_residues:
                        j = "*"
                secondary_structure += j
              except:
                secondary_structure += "C"
        secondary_structure = re.sub("[^EH*]", "C", secondary_structure)
        # Get query secondary structure #
        query_secondary_structure = ""
        hit_positions = range(hit_start - 1, hit_end)
        for i in range(len(hit_alignment)):
            if hit_alignment[i] != "-":
                j = hit_positions.pop(0)
                if query_alignment[i] != "-":
                    query_secondary_structure += secondary_structure[j]
                else:
                    query_secondary_structure += "-"
            else:
                query_secondary_structure += "-"
        # Skip if interface residues are not conserved #
        if (len(re.findall("\*", query_secondary_structure)) != len(re.findall("\*", secondary_structure))) and (no_restrict == False): 
            sys.stdout.write("\t\t-- Skip: interface has different sequence\n")
            sys.stdout.flush()
            continue
        # Skip if insertions/deletions in structured region holding interface residues #
        if (re.search("[\*|E]+\-+[\*|E]+", query_secondary_structure) or re.search("[\*|H]+\-+[\*|H]+", query_secondary_structure)) and (no_restrict == False):
            sys.stdout.write("\t\t-- Skip: insertions/deletions in interface secondary structure\n")
            sys.stdout.flush()
            continue
        # Add unfiltered hit #
        unfiltered_hits.append((hit_obj.sequenceID, query_alignment, hit_alignment, query_start, query_end, hit_start, hit_end, filter_twilight_zone_hits, no_restrict, identity, similarity))
        # Stop hits if best_mode has been set 
        if best_mode: break 

    return unfiltered_hits

def exec_matcher(A, B):
    """
    This function aligns a pair of sequences "A" and "B" using
    matcher from the EMBOSS package.

    @input:
    query_file {filename}
    hit_file {filename}

    @return:
    query_alignment {string} or None
    hit_alignment {string} or None
    query_start {int} or None
    query_end {int} or None
    hit_start {int} or None
    hit_end {int} or None

    """


    try:
        # Initialize #
        alignment = []
        positions = []
        emboss_path = config.get("Paths", "emboss_path")
        matcher=os.path.join(emboss_path,"matcher")
        query=A.replace(":","_")
        if query != A: shutil.copy(A,query)
        template=B.replace(":","_")
        if template != B: shutil.copy(B,template)
        # Exec process #
        try:
           process = subprocess.Popen([matcher, "-asequence", query, "-bsequence", template, "-outfile", "stdout"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        except Exception as e:
           print("ERROR: Failed MATCHER execution")
           print("matcher -asequence "+query+" -bsequence "+template+" -outfile stdout")
           raise e
        # For each line... #
        count=0
        for line in process.stdout:
            count = count +1
            if count>100000:
               return None, None, None, None, None, None, None, None
            line = line.strip()
            # Skip header #
            if "Waterman-Eggert local alignment of two sequences" in line: continue
            # Get % of identity and of similarity #
            if line.startswith("# Identity:"):
                identity = float(line.split("(")[1][:-2])
            # Get % of identity and of similarity #
            elif line.startswith("# Similarity:"):
                similarity = float(line.split("(")[1][:-2])
            # Skip commented lines #
            elif line.startswith("#"): continue
            # Capture alignment #
            m = re.search("^\S+\s+([a-zA-Z-]+)$", line)
            if m:
                alignment.append(m.group(1))
        # Get alignments #
        A_alignment = "".join(alignment[0::2])
        B_alignment = "".join(alignment[1::2])
        # Get start and tails #
        for header, sequence in functions.parse_fasta_file(A):
            A_start = sequence.find(A_alignment.replace("-", ""))
            A_end = A_start + len(A_alignment.replace("-", "")) - 1
        for header, sequence in functions.parse_fasta_file(B):
            B_start = sequence.find(B_alignment.replace("-", ""))
            B_end = B_start + len(B_alignment.replace("-", "")) - 1
        # Correct starts and tails by 1 #
        return A_alignment, B_alignment, A_start + 1, A_end + 1, B_start + 1, B_end + 1, identity, similarity
    except Exception as e:
        raise e
        #return None, None, None, None, None, None, None, None

def get_sequence_to_crystal_correlations(pdb_obj, pdb_chain, gapped=True):
    """
    This function returns the 1 to 1 correlation between sequence and crystal
    positions.

    @input:
    pdb_obj {PDB}
    pdb_chain {string}
    gapped {boolean} default = True

    @return:
    sequence_to_crystal {dictionary} sequence position, crystal position
    crystal_to_sequence {dictionary} crystal position, sequence position
    
    """

    # Initialize #
    sequence_to_crystal = {}
    crystal_to_sequence = {}

    # For each protein chain object... #
    for protein_chain_obj in pdb_obj.proteins:
        # Skip if not the chain of interest #
        if protein_chain_obj.chain != pdb_chain: continue
        # Get crystal sequence #
        crystal_sequence = protein_chain_obj.protein_sequence
        if gapped:
            crystal_sequence = protein_chain_obj.gapped_protein_sequence
        # Get protein sequence positions #
        sequence_positions = range(len(crystal_sequence))
        # Get crystal positions #
        crystal_positions = [aminoacid.number for aminoacid in protein_chain_obj.aminoacids]
        # For each nucleotide... #
        for aminoacid in crystal_sequence:
            if aminoacid == "x":
                sequence_positions.pop(0)
                continue
            sequence_to_crystal[sequence_positions[0]] = crystal_positions[0]
            crystal_to_sequence[crystal_positions[0]] = sequence_positions[0]
            # Remove last #
            sequence_positions.pop(0)
            crystal_positions.pop(0)

    return sequence_to_crystal, crystal_to_sequence

def readjust_alignments_to_crystal_sequence(pdb_chain_obj, query_alignment, hit_alignment, hit_start, hit_end):
    # Initialize #
    sequence_alignment = "-" * len(pdb_chain_obj.gapped_protein_sequence[:hit_start - 1]) + query_alignment + "-" * len(pdb_chain_obj.gapped_protein_sequence[hit_end:])
    structure_alignment = pdb_chain_obj.gapped_protein_sequence[:hit_start - 1] + hit_alignment + pdb_chain_obj.gapped_protein_sequence[hit_end:]

    return sequence_alignment, structure_alignment

def get_protein_models(pir_alignment_obj, pdb_dir, n=1, dummy_dir="/tmp/", verbose=False, optimization=False):

    # Initialize #
    models = []
    modpy_path    =  config.get("Paths", "modpy_path")
    if not os.path.exists(modpy_path):
       modpy_path    =  os.path.join(config.get("Paths","src_path"),config.get("Paths", "modpy_path"))
    modeller_path =  config.get("Paths", "modeller_path")
    # Get current working directory #
    cwd = os.getcwd()
    # Create tmp directory #
    tmp = os.path.join(dummy_dir, str(os.getpid()))
    if not os.path.exists(tmp): os.makedirs(tmp)
    if verbose: sys.stdout.write("\t-- working directory %s\n"%tmp)
    # Change directory #
    os.chdir(tmp)
    # Create PIR alignment file #
    pir_file = pir_alignment_obj.get_structure_name() + "_"+pir_alignment_obj.get_structure_chains()[0]+"_alignment.pir"
    n_pir_file=0
    while fileExist(pir_file): 
      n_pir_file= n_pir_file+1
      pir_file  = pir_alignment_obj.get_structure_name() + "_"+pir_alignment_obj.get_structure_chains()[0]+"_alignment.pir"+"_"+str(n_pir_file)
    pir_alignment_obj.write(pir_file)
    # Create dummy PDB file #
    dummy_pdb_obj = PDB()
    dummy_name = pir_alignment_obj.get_structure_name()
    dummy_file = dummy_name + ".pdb"
    # For each PDB chain... #
    for pdb_chain_obj in sorted(PDB(os.path.join(pdb_dir, "clean", pir_alignment_obj.get_structure_name() + ".pdb")).chains, key=lambda x: x.chain):
        if pdb_chain_obj.chaintype == "P":
            if pdb_chain_obj.chain in pir_alignment_obj._structure_chains:
                if verbose: sys.stdout.write("\t\t-- add chain %s \n"%( pir_alignment_obj.get_structure_name() + "_" + pdb_chain_obj.chain))
                dummy_pdb_obj.add_chain(PDB(os.path.join(pdb_dir, "split", pir_alignment_obj.get_structure_name() + "_" + pdb_chain_obj.chain + ".pdb")).chains[0])
    # For each helix... #
    for line in functions.parse_file(os.path.join(pdb_dir, "helices", pir_alignment_obj.get_structure_name() + "_" + dummy_pdb_obj.chains[0].chain + ".txt")):
        dna_pdb_obj = PDB(os.path.join(pdb_dir, "split", pir_alignment_obj.get_structure_name() + ".dna." + line + ".pdb"))
        for dna_pdb_chain_obj in dna_pdb_obj.chains:
            dummy_pdb_obj.add_chain(dna_pdb_chain_obj)
        break
    # Sort PDB chains #
    dummy_pdb_obj.chains.sort(key=lambda x: x.chain)
    # Create dummy PDB file #
    dummy_pdb_obj.write(dummy_file, force=True)
    if verbose: sys.stdout.write("\t\t-- Using template: " + pir_alignment_obj.get_structure_name() + "; Chains: " + str(pir_alignment_obj._structure_chains) + " \n")
    # Execute modeller #
    flags=" "
    if  verbose:       flags = flags + " -v "
    if  optimization:  flags = flags + " --optimize "
    if  verbose: sys.stdout.write("\t\t\t-- %s\n"%(os.path.join(modeller_path, "modpy.sh") + " " + python + " " + os.path.join(modpy_path, "simpleModel.py") + " --pir=" + pir_file + " --models=" + str(n) + " -v " + " --out " + dummy_name + flags))
    try:
      os.system(os.path.join(modeller_path, "modpy.sh") + " " + python + " " + os.path.join(modpy_path, "simpleModel.py") + " --pir=" + pir_file + " --models=" + str(n) + " -v " + " --out " + dummy_name + flags)
    except Exceptions as e:
      raise("Error on modeling %s"%e)
    # For each model... #
    for pdb_file in os.listdir("."):
        # If model PDB file... #
        m = re.search("query.B(\d+).pdb", pdb_file)
        if m:
            # Get raw model PDB object #
            model_pdb_obj = PDB(pdb_file)
            # Initialize clean model PDB object #
            clean_model_pdb_obj = PDB()
            # For each PDB chain... #
            skip_failure=0
            for i in range(len(model_pdb_obj.chains)):
                # Correct chains #
                if model_pdb_obj.chains[i].chaintype == "P":
                   try:
                    pdb_chain_obj = ChainOfProtein(dummy_pdb_obj.chains[i]._pdb, dummy_pdb_obj.chains[i]._chain)
                    for aminoacid_obj in model_pdb_obj.chains[i].aminoacids:
                        pdb_chain_obj.add_residue(aminoacid_obj)
                   except Exception as e:
                    skip_failure = skip_failure +1
                    if verbose: sys.stdout.write("\t\t\t-- Fail to add protein chain: " + repr(e) + "\n")
                else:
                   try:
                      pdb_chain_obj = ChainOfNucleotide(dummy_pdb_obj.chains[i]._pdb, dummy_pdb_obj.chains[i]._chain)
                      for j in range(len(model_pdb_obj.chains[i].heteroatoms)):
                        pdb_chain_obj.add_residue(ResidueOfNucleotide(model_pdb_obj.chains[i].heteroatoms[j]._number, model_pdb_obj.chains[i].heteroatoms[j]._version, dummy_pdb_obj.chains[i].nucleotides[j]._type, dummy_pdb_obj.chains[i].nucleotides[j]._mode))
                        for atom_obj in model_pdb_obj.chains[i].heteroatoms[j].atoms:
                            pdb_chain_obj.nucleotides[-1].add_atom(AtomOfNucleotide(atom_obj._number, atom_obj._name, atom_obj._coordinates[0], atom_obj._coordinates[1], atom_obj._coordinates[2], atom_obj._occupancy, atom_obj._tempFactor, atom_obj._element, atom_obj._charge))
                   except Exception as e:
                      skip_failure = skip_failure +1
                      if verbose: sys.stdout.write("\t\t\t-- Fail to add non-protein chain: " + repr(e) + "\n")
                if skip_failure==0: clean_model_pdb_obj.add_chain(pdb_chain_obj)
            # Add model to models #
            if skip_failure==0: models.append((int(m.group(1)), clean_model_pdb_obj))
    # Return to original directory #
    os.chdir(cwd)
    # Erase tmp directory #
    if not verbose: shutil.rmtree(tmp)
    # Sort models #
    if len(models)<0: raise("No models accepted")

    models.sort(key=lambda x: x[0])

    return [i[-1] for i in models]

def get_identities(A, B):

    # Initialize #
    identities = 0

    for i in range(len(A)):
        if A[i].isalpha() and B[i].isalpha():
            if A[i] == B[i]:
                identities += 1

    return identities

def get_aligned_residues(A, B):

    # Initialize #
    aligned_residues = 0

    for i in range(len(A)):
        if A[i].isalpha() and B[i].isalpha():
            if A[i] != "x" and  A[i] != "X" and B[i] != "x" and  B[i] != "X":
                aligned_residues += 1

    return aligned_residues

def get_sequence_residues(sequence):

    # Initialize #
    residues = 0

    for i in range(len(sequence)):
        if sequence[i].isalpha():
            if sequence[i] != "x" and  sequence[i] != "X":
                residues += 1

    return residues    

def get_non_overlapping_pir_alignments(alignments_list):
    """
    This function returns a non-overlapping list of PIR alignments. 

    @input:
    alignments_list {list} of {PirAlignment} objects

    @return:
    non_overlapping_alignments_list {list} of non-overlapping {PirAlignment} objects
    
    """

    # Initialize #
    non_overlapping_alignments_list = []
    max_overlap_percentage = int(config.get("Parameters", "max_overlap_percentage"))

    # For each PIR alignment object... #
    for pir_alignment_obj in alignments_list:
        if type(pir_alignment_obj) == "list":
            pir_alignment_obj = pir_alignment_obj[0]
        # Initialize #
        non_overlapping = True
        residues = range(pir_alignment_obj._sequence_residues[0], pir_alignment_obj._sequence_residues[1] + 1)
        # For each PIR alignment object... #
        for non_overlapping_pir_alignment_obj in non_overlapping_alignments_list:
            # Initialize #
            non_overlapping_residues = range(non_overlapping_pir_alignment_obj._sequence_residues[0], non_overlapping_pir_alignment_obj._sequence_residues[1] + 1)
            # Get intersection #
            intersection = [i for i in residues if i in non_overlapping_residues]
            # If overlapping... #
            if (float(len(intersection)) / len(residues)) * 100 > max_overlap_percentage or (float(len(intersection)) / len(non_overlapping_residues)) * 100 > max_overlap_percentage:
                non_overlapping = False
                break
        # If non-overlapping... #
        if non_overlapping: non_overlapping_alignments_list.append(pir_alignment_obj)

    return non_overlapping_alignments_list

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()
    dummy_dir = options.dummy_dir
    if not dummy_dir.startswith("/"): dummy_dir = os.path.abspath(options.dummy_dir)
    if options.input_file  is not None:
     input_file = options.input_file
     if not input_file.startswith("/"): input_file = os.path.abspath(options.input_file)
    output_dir = options.output_dir
    if not output_dir.startswith("/"): output_dir = os.path.abspath(options.output_dir)
    pdb_dir = options.pdb_dir
    if not pdb_dir.startswith("/"): pdb_dir = os.path.abspath(options.pdb_dir)
    resolution_file=options.resolution_file
    if resolution_file is not None:
      if not resolution_file.startswith("/"):resolution_file=os.path.abspath(resolution_file)
    input_name = ".".join(os.path.basename(input_file).split(".")[0:-1])
    info_file = options.info
    if info_file is None: info_file = os.path.join(output_dir,"modelling_results.log")
    if not info_file.startswith("/"): info_file =  os.path.abspath( info_file )

    iteration  = 0
    start_time = float(time.time())

    # Create output and dummy directory #
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(dummy_dir):
        os.makedirs(dummy_dir)

    if not fileExist(info_file):
        info=open(info_file,"w")
        info.write("# Modelling results\n")
        info.close()
    if not fileExist(input_file) or not os.path.exists(pdb_dir):
       sys.stdout.write("Missing input files, please check help with option -h\n")
       info=open(info_file,"a")
       info.write("%s\tFAIL\t%s\n"%(os.path.basename(input_file).rstrip(".fa"),socket.gethostname()))
       info.close()
       exit(0)

    # Initialize #
    resolutions = None
    # Read resolutions file #
    if options.resolution_file is not None:
        # Initialize #
        resolutions = {}
        # For each line... #
        for line in functions.parse_file(resolution_file):
            if line.startswith("#"): continue
            line = line.split(";")
            resolutions.setdefault(line[0], float(line[1]))

    # Get data on families
    families = {}
    if options.verbose:sys.stdout.write("Check families...\n")
    if not os.path.exists(os.path.join(pdb_dir, "families.txt")):
       sys.stdout.write("Families file %s is not used\n"%(os.path.join(pdb_dir, "families.txt")))
    else:
      for line in functions.parse_file(os.path.join(options.pdb_dir, "families.txt")):
        if line.startswith("#"): continue
        pdb_chain, family = line.split(";")
        families[pdb_chain] = family


    # Get PIR alignment object #

    # Initialize #
    sequence = ""
    sequence_file = os.path.join(dummy_dir, "%s.fasta" % input_name+"_"+str(os.getpid()))
    if options.label is None: label=input_name
    else: label=options.label
    # Remove  sequence file if already exists #
    if os.path.exists(sequence_file): os.remove(sequence_file)

    if options.threading:
        # Reading the threading file #
        full_mode=False
        thread_object   = threader.Threaded(input_file)
        hit_ID          = thread_object.get_pdb_name()+"_"+thread_object.get_pdb_chain()
        query_alignment = thread_object.get_alignment("query")
        hit_alignment   = thread_object.get_alignment("pdb")
        hit_start       = [key[1] for key in sorted(thread_object.get_protein(), key=lambda x: x[-1])][0]
        hit_end         = [key[1] for key in sorted(thread_object.get_protein(), key=lambda x: x[-1])][-1]
        query_start     = 1
        query_end       = len(query_alignment.replace("-",""))
        kmers_thread    = thread_object.get_kmers_fixed()
        #Generate sequence input file
        header= os.path.basename(input_name)
        functions.write(sequence_file, ">%s\n%s" % (header, query_alignment.replace("-","")))
        aligned_residues= get_aligned_residues(query_alignment, hit_alignment)
        identity   =int((get_identities(query_alignment, hit_alignment) * 100) / (float(aligned_residues) + 0.5))
        similarity =identity
        hit_object = [hit_ID, query_alignment, hit_alignment, query_start, query_end, hit_start, hit_end, identity, similarity]
        try:
               pir_alignments_list = get_pir_alignments(sequence_file, pdb_dir, full_mode, options.filter_twilight_zone_hits, resolutions, options.select_mode, dummy_dir, thread_obj=hit_object, \
                                                        no_restrict_usage=options.no_restrict,dimer_mode=options.dimer_mode, monomer_mode=options.monomer_mode, verbose=options.verbose, \
                                                        filter_restrictive=options.filter_restrictive,best_mode=options.best, iteration=iteration, start_time=start_time)
        except Exception as e:
               sys.stdout.write("ERROR when getting PIR alignments\n")
               sys.stdout.write("ERROR: %s\n"%e)
               info=open(info_file,"a")
               info.write("%s\tFAIL\t%s\n"%(os.path.basename(input_file).rstrip(".fa"),socket.gethostname()))
               info.close()
               exit(0)
    else:   
        # For each line in input sequence... #
        for header_title, sequence in functions.parse_fasta_file(input_file):
            header_words=header_title.strip().split()
            header=header_words[0].lstrip(">").split("|")[-1]
            if options.mutation is not None:
                m = re.search("^([a-zA-Z])(\d+)([a-zA-Z])$", options.mutation)
                if sequence[int(m.group(2)) - 1] == m.group(1):
                    sequence = list(sequence)
                    sequence[int(m.group(2)) - 1] = m.group(3)
                    sequence = "".join(sequence)
                else:
                    raise ValueError("ERROR: Fail to introduce mutation %s! Amino acid %s is \"%s\" and not \"%s\"" % (options.mutation, m.group(2), sequence[int(m.group(2)) - 1], m.group(1)))
            functions.write(sequence_file, ">%s\n%s" % (header, sequence))
            if options.mutation is not None:
                break
        # Get list of PIR alignments #
        if options.verbose: sys.stdout.write("Homology Modelling of "+header+" ...\n")
        try:
           pir_alignments_list = get_pir_alignments(sequence_file, pdb_dir, options.full_mode, options.filter_twilight_zone_hits, resolutions, options.select_mode, dummy_dir,  \
                                                    no_restrict_usage=options.no_restrict,dimer_mode=options.dimer_mode, monomer_mode=options.monomer_mode, verbose=options.verbose, \
                                                    filter_restrictive=options.filter_restrictive,best_mode=options.best, iteration=iteration, start_time=start_time)
        except Exception as e:
           sys.stdout.write("ERROR when getting PIR alignments\n")
           sys.stdout.write("ERROR: %s\n"%e)
           info=open(info_file,"a")
           info.write("%s\tFAIL\t%s\n"%(os.path.basename(input_file).rstrip(".fa"),socket.gethostname()))
           info.close()
           exit(0)
    if pir_alignments_list is None:
        if options.verbose: sys.stdout.write("No alignments found. Exit\n")
        info=open(info_file,"a")
        info.write("%s\tFAIL\t%s\n"%(os.path.basename(input_file).rstrip(".fa"),socket.gethostname()))
        info.close()
        exit(0)
    # For each PIR alignment list... #
    numb = 0
    summary_file = os.path.join(output_dir, label + "_model.summary.txt")
    if os.path.isfile(summary_file):
        os.system("rm -f " + summary_file)
    functions.write(summary_file, "#model;n;template;N-tails;C-tails;protein-chain;DNA-chain;#identities;#coverage;template_by_RMSD;domain;%_identity;%_similarity")
    domain_d = {}
    for pir_alignment_obj in pir_alignments_list:
        if pir_alignment_obj is None:continue
        domain_d.setdefault(pir_alignment_obj.get_domain(),[]).append(pir_alignment_obj)
    for domain in domain_d.iterkeys():
        for pir_alignment_obj in domain_d[domain]:
            # Initialize #
            protein_chains = []
            dna_chains = []
            identities = []
            coverage = []
            # For each chain... #
            if options.verbose: sys.stdout.write("\t-- Modelling with template " + pir_alignment_obj._hit + " ...\n")
            for j in range(len(pir_alignment_obj.get_structure_chains())):
                pdb_chain = pir_alignment_obj._structure_chains[j]
                sequence_alignment = pir_alignment_obj._sequence_alignments[j]
                structure_alignment = pir_alignment_obj._structure_alignments[j]
                if re.search("^\.+$", sequence_alignment):
                    dna_chains.append(pdb_chain)
                else:
                    protein_chains.append(pdb_chain)
                    aligned_residues = get_aligned_residues(sequence_alignment, structure_alignment)
                    identities.append(str(int(((get_identities(sequence_alignment, structure_alignment) * 100) / float(aligned_residues)) + 0.5)))
                    coverage.append(str(int(((aligned_residues * 100) / float(get_sequence_residues(structure_alignment))) + 0.5)))
            # Initialize "n" parameter #
            n = options.n_model
            if n==1 and options.n_total is not None:
                if not options.all_models: n = options.n_total
                elif options.n_total < len(domain_d[domain]): n = 1 
                elif options.n_total % len(domain_d[domain]) == 0: n = options.n_total / len(domain_d[domain])
                elif options.n_total % len(domain_d[domain]) != 0: n = options.n_total / len(domain_d[domain]) + 1
                else: n = 1

                
            # Get models only if they haven't been created previously #
            exists = True
            model_list = []
            for i in range(1, n+1):
                if pir_alignment_obj.get_main_order()[0]==1: pdb_file = os.path.join(output_dir, label + ":" + str(pir_alignment_obj._sequence_residues[0]) + ":" + str(pir_alignment_obj._sequence_residues[1]) + "_" + pir_alignment_obj._hit + "_" + str(i) + ".pdb")
                else: pdb_file = os.path.join(output_dir, label + ":" + str(pir_alignment_obj._sequence_residues[2]) + ":" + str(pir_alignment_obj._sequence_residues[3]) + "_" + pir_alignment_obj._hit + "_" + str(i) + ".pdb")
                if fileExist(pdb_file):
                    model_list.append(pdb_file)
                else:
                    exists = False


            if exists == False or options.force:
              try:
                models = get_protein_models(pir_alignment_obj, pdb_dir, n, dummy_dir, verbose=options.verbose, optimization = options.optimization)
                # For each model... #
                for j in range(n):
                    if j+1>len(models): break
                    # Initialize #
                    if options.verbose:sys.stdout.write("\t\t\t-- Modeling %d ...\n"%j)
                    pdb_obj = models[j]
                    # Write model #
                    if pir_alignment_obj.get_main_order()[0]==1: pdb_filename=label + ":" + str(pir_alignment_obj._sequence_residues[0]) + ":" + str(pir_alignment_obj._sequence_residues[1]) + "_" + pir_alignment_obj._hit + "_" + str(j+1) + ".pdb"
                    else: pdb_filename=label + ":" + str(pir_alignment_obj._sequence_residues[2]) + ":" + str(pir_alignment_obj._sequence_residues[3]) + "_" + pir_alignment_obj._hit + "_" + str(j+1) + ".pdb"
                    dummy_name = "dummy_" + pdb_filename + "_" + str(os.getpid())
                    dummy_file = os.path.join(dummy_dir, "%s.pdb" % ( dummy_name ) )
                    if options.verbose: sys.stdout.write("\t\t\t-- Write dummy %s\n"%dummy_name)
                    pdb_obj.write(dummy_file, force=True)
                    # Test model protein-DNA interface #
                    x3dna_obj = x3dna.get_x3dna_obj(dummy_file)
                    contacts_obj = contacts.get_contacts_obj(pdb_obj, x3dna_obj)
                    interface_obj = interface.get_interface_obj(pdb_obj, x3dna_obj, contacts_obj)
                    if interface_obj.get_interface_length() == 0 and not options.full_fragments_mode and not options.threading:
                        # Protein model does not interact w/ DNA #
                        if options.verbose: sys.stdout.write("\t\t-- Skip model without Protein-DNA contacts "+ pdb_filename + " ....\n")
                        if os.path.exists(dummy_file)and not options.verbose: os.remove(dummy_file)
                        continue
                    # If model DNA... #
                    if options.dna:
                        # Initialize #
                        dummy_file = os.path.join(dummy_dir, "%s.pdb" % dummy_name )
                        # Write model #
                        pdb_obj.write(dummy_file, force=True)
                        if options.threading:
                           dna_sequence    = [dna_seq for dna_seq in kmers_thread.iterkeys()][0]
                           interface_start =  int(kmers_thread[dna_sequence])
                           interface_range =  range(interface_start, interface_start + len(dna_sequence) )
                           if options.verbose: sys.stdout.write("\t\t-- Model DNA sequence %s\n"%dna_sequence)
                           pdb_obj = model_dna.get_dna_model_pdb_obj(dummy_file, dna_sequence, x3dna_obj, interface_obj, interface_range, dummy_dir=dummy_dir)
                        else:
                           # Get DNA sequence #
                           dna_sequence = x3dna_obj.get_nucleotide_sequence(interface_obj.get_interface_start(), interface_obj.get_interface_end())
                           # Get model object #
                           if options.verbose: sys.stdout.write("\t\t-- Model DNA sequence %s\n"%dna_sequence)
                           if interface_obj.get_interface_start() is not None:
                              interface_range =  range(interface_obj.get_interface_start(),interface_obj.get_interface_start()+len(dna_sequence) )
                              pdb_obj = model_dna.get_dna_model_pdb_obj(dummy_file, dna_sequence, x3dna_obj, interface_obj, interface_range, dummy_dir=dummy_dir)
                    # If resolutions file... #
                    if options.resolution_file is not None:
                        # Write summary #
                        if pir_alignment_obj.get_structure_name() in resolutions: functions.write(summary_file, "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s" % (numb, (j + 1), pir_alignment_obj.get_structure_name(), ",".join(map(str, pir_alignment_obj._sequence_residues[0::2])), ",".join(map(str, pir_alignment_obj._sequence_residues[1::2])), ",".join(protein_chains), ",".join(dna_chains), ",".join(identities), ",".join(coverage), resolutions[pir_alignment_obj.get_structure_name()], pir_alignment_obj.get_domain(), str(pir_alignment_obj.get_identity()), str(pir_alignment_obj.get_similarity())))
                        else: functions.write(summary_file, "%s;%s;%s;%s;%s;%s;%s;%s;%s;NA;%s;%s;%s" % (numb, (j + 1), pir_alignment_obj.get_structure_name(), ",".join(map(str, pir_alignment_obj._sequence_residues[0::2])), ",".join(map(str, pir_alignment_obj._sequence_residues[1::2])), ",".join(protein_chains), ",".join(dna_chains), ",".join(identities), ",".join(coverage), pir_alignment_obj.get_domain(), str(pir_alignment_obj.get_identity()), str(pir_alignment_obj.get_similarity())))
                    else:
                        numb += 1
                        functions.write(summary_file, "%s;%s;%s;%s;%s;%s;%s;%s;%s;NA;%s;%s;%s" % ((numb), (j + 1), pir_alignment_obj.get_structure_name(), ",".join(map(str, pir_alignment_obj._sequence_residues[0::2])), ",".join(map(str, pir_alignment_obj._sequence_residues[1::2])), ",".join(protein_chains), ",".join(dna_chains), ",".join(identities), ",".join(coverage), pir_alignment_obj.get_domain(), str(pir_alignment_obj.get_identity()), str(pir_alignment_obj.get_similarity())))
                        # Remove file #
                        if os.path.exists(dummy_file) and not options.verbose: os.remove(dummy_file)
                        # Initialize #
                        pdb_file = os.path.join(output_dir, pdb_filename)
                        if options.verbose: sys.stdout.write("\t\t-- Model created as: " + pdb_filename + " ...\n")
                        # Write model #
                        if options.chains_fixed: 
                            pdb_obj=chains_fixed(pdb_obj)
                        pdb_obj.write(pdb_file, force=True)
                        if options.renumerate:
                           if options.verbose: sys.stdout.write("\t\t-- Renumerate residues as in original sequence %s\n"%input_name)
                           sequences_modelled={}
                           sequences_input={}
                           for header, sequence in functions.parse_fasta_file(sequence_file):
                               sequences_input.setdefault(header,sequence)
                           for chain_id in pdb_obj.chain_identifiers:
                               chain=pdb_obj.get_chain_by_id(chain_id)
                               if chain.chaintype == "P":
                                  sequences_modelled.setdefault(chain_id,chain.protein_sequence)    
                           homodimer=False
                           if (len(sequences_input.keys())==1 and len(sequences_modelled.keys())>len(sequences_input.keys())) and not options.monomer_mode: homodimer=True
                           if options.verbose: sys.stdout.write("\t\t-- Homodimer %s\n"%homodimer)
                           output_file=os.path.basename(pdb_file)
                           path_file=os.path.dirname(pdb_file)
                           try:
                              chain_name_pairs=bijection_sequence_to_sequence(sequences_modelled,sequences_input,dummy_dir,homodimer=homodimer)
                              #print "PAIRS ",chain_name_pairs
                              sequences_complex = {}
                              for chain_id in pdb_obj.chain_identifiers:
                                header=chain_name_pairs.get(chain_id)
                                sequences_complex.setdefault(chain_id,(header,sequences_input.get(header))) 
                              #print "COMPLEXES ",sequences_complex
                              pdb_renumber=PDB()
                              pdb_renumber=renumber_pdb(pdb_file,sequences_complex,dummy_dir)
                              if options.chains_fixed: 
                                  pdb_renumber=chains_fixed(pdb_renumber)
                              if len(pdb_renumber.chain_identifiers)>1:
                                  pdb_renumber.write(pdb_file,force=True)
                              else:
                                  if verbose: sys.stdout.write("\t\t\t-- Failed to renumerate %s\n"%output_file)
                           except Exception as e:
                              if options.verbose: sys.stdout.write("\t\t\t-- Failed to renumerate %s\n"%output_file)
                              if options.verbose: sys.stderr.write("\t\t\t-- Error while renumbering %s\n"%str(repr(e)))
              except Exception as e:
                if options.verbose: sys.stdout.write("\t\t\t-- Failed Modelling with template " + pir_alignment_obj._hit + " ...\n")
                if options.verbose: sys.stdout.write("\t\t\t-- Error while getting models %s\n"%str(repr(e)))
                continue
            else:
                template_name=" "
                for model in model_list:
                    numb += 1
                    if template_name==pir_alignment_obj.get_structure_name(): j=j+1
                    else: j=1
                    if j==1: template_name=pir_alignment_obj.get_structure_name()
                    functions.write(summary_file, "%s;%s;%s;%s;%s;%s;%s;%s;%s;NA;%s;%s;%s" % ((numb), (j), pir_alignment_obj.get_structure_name(), ",".join(map(str, pir_alignment_obj._sequence_residues[0::2])), ",".join(map(str, pir_alignment_obj._sequence_residues[1::2])), ",".join(protein_chains), ",".join(dna_chains), ",".join(identities), ",".join(coverage), pir_alignment_obj.get_domain(), str(pir_alignment_obj.get_identity()), str(pir_alignment_obj.get_similarity())))
                    # Initialize #
                    if options.verbose: sys.stdout.write("\t\t-- Reuse PDB model: " + os.path.basename(model) + "\n")
                
        # Stop modelling #
        if not options.all_models:
            break
    # Finalize
    info=open(info_file,"a")
    info.write("%s\tDONE\n"%(os.path.basename(input_file).rstrip(".fa")))
    info.close()
    if options.verbose: print("Done")

    # Remove sequence file #
    try:
        if os.path.exists(dummy_dir) and not options.verbose: shutil.rmtree(dummy_dir)
    except:
        if options.verbose: sys.stdout.write("Dummy directory %s cannot be removed\n"%dummy_dir)
        exit(0)
