import os, sys, re
from collections import Counter
import ConfigParser
import itertools
import numpy
import optparse
import subprocess
import time

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
import random
# Imports jbonet's module #
from SBI.data import aminoacids1to3, aminoacids_polarity_boolean, nitrogenous_bases
from SBI.structure import PDB

# Import my modules #
import contacts, dssp, interface, spotentials, tmalign, triads, x3dna, threading_to_triads

#-------------#
# Classes     #
#-------------#

class MSA(object):
    """
    This class defines a {MSA} object.

    """
    residues = []
    def __init__(self, file_name=None, motif_name=None, option="msa"):
        self._file = file_name
        if motif_name is not None:
           self._motif= motif_name
        elif self._file is not None:
           name=file_name.split("/")[-1]
           if name.endswith(".msa"):self._motif=name.rstrip(".msa")
           elif name.endswith(".pwm"):self._motif=name.rstrip(".pwm")
           elif name.endswith(".txt"):self._motif=name.rstrip(".txt")
           elif name.endswith(".meme"):self._motif=name.rstrip(".meme")
           else:self._motif=name
        else:
           self._motif="UNDEFINED"
        self._option = option
        self._binding_site_length = 0
        self._sequences = set()
        self._pwm = []
        self._max_sequences=(0,0,0,1000)
        self._min_binding_site_length=(0,0,0)
        # Initialize #
        if file_name is not None:
            self._parse_file()

    def _parse_file(self):
        option=self._option
        if option != "msa" and option != "pwm" and option != "meme" and option != "txt":
           option = "msa"
        if option=="msa": 
           self.clean_pwm()
           for line in functions.parse_file(self._file):
               if line.startswith("#"): continue
               sequence,score = line.strip().split(";")
               self.add_sequence(sequence, float(score))
           if len(self._sequences) <1000:
              repeat=int(1000.0/len(self._sequences))
              original_set=set()
              for sequence,score in  self._sequences:
                  original_set.add((sequence,score))
              for i in xrange(repeat):
                  for sequence,score in original_set:
                      self.add_sequence(sequence, float(score-0.01/float(i+1)))
           self.set_binding_site_length()
           self.set_pwm()
        if option=="pwm":
           matrix={}
           for line in functions.parse_file(self._file):
               data=line.strip().split("\t")
               if not data[0].endswith(":"): continue
               residue=data[0].strip(":")
               self.set_binding_site_length(len(data)-1)
               matrix.setdefault(residue,data[1:])
           for i in xrange(self.get_binding_site_length()):
               vector = []
               for j in xrange(len(self.get_residues())):
                 residue=self.residues[j]
                 vector.append(matrix[residue][i])
               self._pwm.append(vector)
           self.set_binding_site_length(len(self._pwm))
        if option=="meme":
           header=True
           for line in functions.parse_file(self._file):
             data=line.strip().split()
             if len(data)<=0: continue
             if data[0]=="ALPHABET=":
                header=True
                self.clean_pwm()
             if header:
               if data[0]=="ALPHABET=":
                  residues=list(data[1])
               if data[0]=="MOTIF":
                   if self._motif=="UNDEFINED":
                      self.set_motif(data[1])
               if data[0]=="letter-probability":
                   self._binding_site_length = int(data[5])
                   header=False
             else:
               vector=[]
               if len(data) == len(self.get_residues()):
                 for j in xrange(len(self.get_residues())):
                   for i in xrange(len(residues)):
                      if  self.residues[j]==residues[i]:
                          vector.append(data[i])
                 self._pwm.append(vector)
           self.set_binding_site_length(len(self._pwm))
        if option=="txt":
           header=True
           for line in functions.parse_file(self._file):
             if header:
               data=line.strip().split("\t")
               residues=data[1:]
               header=False
             else:
               data=line.strip().split("\t")
               self.set_binding_site_length(int(data[0]))
               vector=[]
               for j in xrange(len(self.get_residues())):
                   for i in xrange(len(residues)):
                      if  self.residues[j]==residues[i]:
                          vector.append(data[i+1])
               self._pwm.append(vector)
           self.set_binding_site_length(len(self._pwm))

    def __add__(self,other):
        x=self.__class__()
        x.set_motif(self.get_motif()+"_"+other.get_motif())
        x._pwm.extend(self.get_pwm())
        x._pwm.extend(other.get_pwm())
        x.set_binding_site_length(len(x._pwm))
        x.enhance()
        return x

    def extend_residues(self,new_list):
        self.residues.extend(new_list)

    def enhance(self):
        pwm=[]
        for vector in self._pwm:
            maxi=max([float(x) for x in vector])
            mini=min([float(x) for x in vector])
            v=[]
            for j in xrange(len(self.get_residues())): v.append(vector[j])
            if maxi>mini and mini>0:
               if (maxi-mini)/mini > 0.05:
                  v=[]
                  suma=sum([float(x)-mini for x in vector])
                  for j in xrange(len(self.get_residues())):
                      v.append("%8.4f"%((float(vector[j]) - mini)/suma))
            pwm.append(v)
        self._pwm=pwm 

    def combine(self,other,overlap):
        x=self.__class__()
        x.set_motif(self.get_motif()+"_"+other.get_motif())
        pwm_a=self.get_pwm()
        pwm_b=other.get_pwm()
        x._pwm=pwm_a[:-overlap]
        for i in xrange(overlap):
            vector=[]
            for j in xrange(len(self.get_residues())):
                vector.append("%8.4f"%((float(pwm_b[i][j])+float(pwm_a[i-overlap][j]))/2.0))
            x._pwm.append(vector)
        x._pwm.extend(pwm_b[overlap:])
        x.set_binding_site_length(len(x._pwm))
        x.enhance()
        return x

    def section(self,start=0,end=0):
        x=self.__class__()
        if end==0:end=len(self._pwm)
        x.set_motif(self.get_motif())
        x._pwm=self.get_pwm()[start:end]
        x.set_binding_site_length(len(x._pwm))
        return x

    def add_sequence(self, sequence, score):
        self._sequences.add((sequence, score))

    def set_binding_site_length(self,length=0):
        if length == 0:
           binding = zip(*[i for i in self.get_sequences()])
           self._binding_site_length = len(binding)
        else:
           self._binding_site_length = length

    def set_motif(self,name):
        if name.endswith(".msa"):name.rstrip(".msa")
        if name.endswith(".pwm"):name.rstrip(".pwm")
        if name.endswith(".txt"):name.rstrip(".txt")
        if name.endswith(".meme"):name.rstrip(".meme")
        self._motif=name

    def set_pwm(self):
        self._pwm=[]
        if len(self._sequences)>0:
           pfm = []
           binding = zip(*[i[0] for i in self.get_sequences()])
           self._binding_site_length = len(binding)
           sequences=self.get_sequences()
           # Count residues instances #
           for i in range(self.get_binding_site_length()):
             vector=[]
             for res in self.get_residues():
               vector.append(binding[i].count(res)+1)
             pfm.append(vector)
           # For each position... #
           for i in range(self.get_binding_site_length()):
             vector=[]
             for j in range(len(self.get_residues())):
               vector.append("%.3f" % (float(pfm[i][j]) / sum(pfm[i])))
             self._pwm.append(vector)                           

    def set_sequences(self):
        if self._binding_site_length <self._min_binding_site_length[0]:
           number_of_sequences=self._max_sequences[0]
        elif self._binding_site_length >=self._min_binding_site_length[0] and self._binding_site_length < self._min_binding_site_length[1]:
           number_of_sequences=self._max_sequences[1]
        elif self._binding_site_length >=self._min_binding_site_length[1] and self._binding_site_length < self._min_binding_site_length[2]:
           number_of_sequences=self._max_sequences[2]
        else:
           number_of_sequences=self._max_sequences[3]
        if len(self._sequences) <= 0 and len(self._pwm)>0:
           sequence_residue=[]
           for i in xrange(self.get_binding_site_length()):
               vector=[]
               for j in xrange(len(self.get_residues())):
                   vector.append(int(number_of_sequences*float(self._pwm[i][j])))
               sequence_residue.append(vector)
           done=(sum([sum([n for n in sequence_residue[i]]) for i in xrange(self.get_binding_site_length())]) == 0 )
           for s in xrange(number_of_sequences):
               seq=""
               if done: break
               for i in xrange(self.get_binding_site_length()):
                   add_one=True
                   nr=int(1000*self.get_binding_site_length()*random.random())
                   for jj in xrange(len(self.get_residues())):
                       if not add_one: continue
                       jjj=float(jj+nr)/float(len(self.get_residues()))
                       j=jj+nr-int(jjj)*len(self.get_residues())
                       base = self.residues[j]
                       if sequence_residue[i][j]>0:
                          seq = seq + base
                          sequence_residue[i][j] = sequence_residue[i][j] - 1
                          add_one = False
               for i in xrange(self.get_binding_site_length()):
                 if done: continue
                 if not done:
                  done = (sum([n for n in sequence_residue[i]]) == 0)
               if (len(seq)==self.get_binding_site_length()):
                 if (seq,1.0) not in self.get_sequences():
                    self.add_sequence(seq,1.0) 
                 else:
                    sc=min([float(x[1]) for x in self.get_sequences()])
                    self.add_sequence(seq,sc-0.001/(s+1))

    def get_main_sequence(self):
        if len(self._pwm)<=0:
            self.set_pwm()
        main=""
        for vector in self._pwm:
            maxim=max([float(x) for x in vector])
            skip=False
            for i in xrange(len(vector)):
                if float(vector[i]) == maxim:
                  if not skip:
                    main = main + self.residues[i]
                    skip = True
            if not skip:
               for i in xrange(len(vector)):
                if float(vector[i]) > (maxim - 1.0e-8):
                  if not skip:
                    main = main + self.residues[i]
                    skip = True
        return main

    def clean_pwm(self):
        self._pwm=[]

    def get_pwm(self):
        return self._pwm

    def get_motif(self):
        name=self._motif
        if name.endswith(".msa"):name.rstrip(".msa")
        if name.endswith(".pwm"):name.rstrip(".pwm")
        if name.endswith(".txt"):name.rstrip(".txt")
        if name.endswith(".meme"):name.rstrip(".meme")
        return name

    def get_option(self):
        return self._option

    def get_residues(self):
        return self.residues

    def get_sequences(self):
        return self._sequences

    def clean_sequences(self):
        self._sequences=set()

    def get_binding_site_length(self):
        return  self._binding_site_length

    def trim_pwm(self):
        trim = True
        freq=1.0/len(self.residues)
        vector=[]
        for j in range(len(self.residues)):
            vector.append(freq)
        while trim == True:
            if (self._pwm[0] == vector) or (self._pwm[-1] == vector):
                if self._pwm[0] == vector:
                    self._pwm = self._pwm[1:]
                if self._pwm[-1] == vector:
                    self._pwm = self._pwm[:-1]
            else:
                trim = False

    def information_content(self):
        try:
          x=numpy.array(self._pwm)
          pwm=x.astype(numpy.float)
          lpwm=numpy.log2(pwm)
          pwm[numpy.isinf(pwm)]=0
          pwm[numpy.isnan(pwm)]=0
          lpwm[numpy.isinf(lpwm)]=0
          lpwm[numpy.isnan(lpwm)]=0
          ic = sum(sum( -1. * pwm * lpwm) )
        except:
          ic = 0.0
        return ic

    def write(self, file_name=None, option="msa", overwrite=False):
        #If overwrite a file
        if overwrite == True:
            os.system("rm -f " + file_name)
        if self._binding_site_length == 0:
            self.set_binding_site_length()
        # If binding site... #
        if self.get_binding_site_length() > 0:
            if option != "msa" and option != "pwm" and option != "meme":
                option = "msa"
            if option == "msa":
                # Write output #
                if len(self._sequences) > 0 :
                    functions.write(file_name, "#sequence;score")
                    for sequence, score in sorted(self._sequences, key=lambda x: x[1], reverse=True):
                        functions.write(file_name, "%s;%s" % (sequence, score))
            if option == "pwm" or  option == "meme":
               # Initialize #
               residues=self.get_residues()
               if len(self._pwm) > 0:
                    pwm = self.get_pwm()
               else:
                    print "Warning: set up PWM of ",len(residues)," residues"
                    self.set_pwm()
                    pwm = self.get_pwm()
               if option == "pwm":
                    # For residue... #
                    functions.write(file_name,self.get_motif())
                    for i in range(len(residues)):
                        functions.write(file_name, "%s:\t%s" % (residues[i], "\t".join([row[i] for row in pwm])))
               else:
                    alphabet="".join(self.residues)
                    freq=1.0/len(self.residues)
                    frequencies=""
                    for res in self.residues:
                      frequencies += " %s %7.5f"%(res,freq)
                    header="MEME version 4.4\n\nALPHABET= %s\n\nstrands: + -\n\nBackground letter frequencies (from uniform background):\n%s\n\n"%(alphabet,frequencies)
                    header+="MOTIF %s %s\n\n"%(self.get_motif(), self.get_motif())
                    header+="letter-probability matrix: alength= %d w= %d nsites= 20 E= 0"%(len(self.residues),self.get_binding_site_length())
                    functions.write(file_name, "%s"%header)
                    for row in pwm:
                      functions.write(file_name, "\t".join(["%10.6f"%(float(x)) for x in row]))
                    functions.write(file_name, "")
        else:
            raise ValueError("No DNA sequences bound!")

class nMSA(MSA):
    residues = list("ACGT")
    def __init__(self,file_name=None, motif_name=None, option="msa"):
      self._nucleotides = self.residues
      self._max_sequences=(100,250,500,1000)
      self._min_binding_site_length=(4,6,8)
      super(nMSA,self).__init__(file_name, motif_name, option)
    def get_nucleotides(self):
      return self._nucleotides


class pMSA(MSA):
    residues = list("ACDEFGHIKLMNPQRSTVWY")
    def __init__(self,file_name=None, motif_name=None, option="msa"):
      self._aminoacids = self.residues
      self._max_sequences=(100,500,1000,5000)
      self._min_binding_site_length=(4,8,10)
      self._option=option
      super(pMSA,self).__init__(file_name, motif_name,option)
    def get_aminoacids(self):
        return self._aminoacids



#-------------#
# Options     #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python pwm.py -i input_files --pbm=pbm_dir --pdb=pdb_dir [-o output_dir] [-a -f -p -s potential -t --threading] [--parallel --info LOG_FILE  --complete COMPLETE ]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_files", help="PDB file or THREADING file or directory of files (e.g. from model_protein.py)", metavar="{filename}")
    parser.add_option("-o", action="store", default="output_pwm", type="string", dest="output_file", help="Output file (if single file as input) or directory name (default = 'output_pwm')", metavar="{filename}")
    parser.add_option("--complete",default=1.00, action="store", type="float", dest="complete", help="Ratio of completness over the total number of profiles top be done(default= 0.95). This is useful in the server to stop the profiler when the time of execution exceeds more than 48 hours ", metavar="RATIO")
    parser.add_option("--pbm", action="store", type="string", dest="pbm_dir", default=None, help="PBM directory (i.e. output dir from pbm.py)", metavar="{directory}")
    parser.add_option("--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (i.e. output dir from pdb.py)", metavar="{directory}")
    parser.add_option("-k","--known", default=False, action="store_true", dest="known", help="The name is of a known PDB file, with format 'code_chain' (default = False)")
    parser.add_option("--meme", default=False, action="store_true", dest="meme", help="Use 'uniprobe2meme' to calculate the PWM matrix for 'FIMO' (default = False)")
    parser.add_option("-r","--reset", default=False, action="store_true", dest="reset", help="Clean the sequences of the original MSA and reset them by a random selection in accordance with the PWM (default = False)")
    parser.add_option("--fragment", default=None, action="store", type="string", dest="fragment_restrict", help="Fragment of protein to apply the potential. Format is 'a-b;c-d': two regions between residues a-b and c-d. Example: '45_A-48_A;50_A-51_A' (Default is None it applies to all amino-acids)")
    parser.add_option("--binding", default=None, action="store", type="string", dest="binding_restrict", help="Binding site of DNA to apply the potential. Format is 'a-b;c-d': two regions between residues a-b and c-d of the forward chain (first in PDB). (Default is None it applies to all nucleotides)")
    parser.add_option("--threading", default=False, action="store_true", dest="threading", help="Input file is a threading file of a PDB structure that (MUST!) exist in the PDB folder of ModCRE (default = False)")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")
    parser.add_option("--parallel", default=False, action="store_true", dest="parallel", help="Run in parallel if the input is a directory (default = False)")
    parser.add_option("--info",default=None,action="store", type="string", dest="info",help="Information LOG file of PWMs that have failed and have been completed")
    parser.add_option("--reuse",default=False, action="store_true", dest="reuse", help="Reuse the information files. If the flag is used then profiles that had failed will remain as FAILED, otherwise it tries to redo them (default=False)")

    group = optparse.OptionGroup(parser, "Statistical potentials", "Select your statistical potentials of choice. By default it uses S3DCdd general potential derived from PDB (the simplest one). In \"--auto\" mode, the program uses S3DCdd family potentials derived from both PDB and/or PBM data and/or approached by Taylor as selected in Potentials configuration file. In case family potentials cannot be applied, the program uses general potentials derived from both PDB and PBM data and approached by Taylor. \"-a\" option overrides options \"-f\", \"-p\" and \"-t\".")
    group.add_option("-a", "--auto", default=False, action="store_true", dest="auto_mode", help="Automate the selection of statistical potentials (default = False)")
    group.add_option("-f", "--family", default=False, action="store_true", dest="family_potentials", help="Use family potentials (default = False)")
    group.add_option("-p", default=False, action="store_true", dest="pbm_potentials", help="Use potentials derived from both PBM + PDB data (default = False)")
    group.add_option("-s", default="s3dc_dd", action="store", type="string", dest="split_potential", help="Split-potential to be used (3d, 3dc, s3dc, s3dc_dd, s3dc_di, pair; default = s3dc_dd)", metavar="{string}")
    group.add_option("-m", "--pmf", default=False, action="store_true", dest="pmf", help="Use of raw mean-force potentials with no Z-scoring (default = False)")
    group.add_option("-t", action="store", default=None, type="float", dest="score_threshold", help="Threshold on the scaled score to consider positive k-mers (default is taken from config file)", metavar="{float}")
    group.add_option("--taylor", default=False, action="store_true", dest="taylor_approach", help="Approach PMF by Taylor (default = False)")
    group.add_option("-b", "--bins", default=False, action="store_true",  dest="computation", help="Computate the potentials: by bins (if selected) or accumulative (default).", metavar="{boolean}")
    group.add_option("--file",default=None, action="store", type="string", dest="potential_file", help="Use potentials from specific file (default = None)", metavar="{string}")
    group.add_option("--radius",default=0, action="store", type="float", dest="radius", help="Maximum contact distance to calculate interactions (default=0 implies the use of 'max_contact_distance' or family-specific radius from configuration", metavar="{string}")
    group.add_option("--methylation",default=False, action="store_true", dest="methylation", help="Flag to use methylated cytosines with binding/non-binding specificity (default=False)")

  

    parser.add_option_group(group)
    
    (options, args) = parser.parse_args()

    if (options.input_files is None and not os.path.exists(options.output_file+".msa")) or options.pdb_dir is None :
        parser.error("missing arguments: type option \"-h\" for help")
    if not re.search("^3d|3dc|s3dc$|^s3dc_dd$|^s3dc_di|pair$", options.split_potential):
        parser.error("incorrect value for -s argument: type option \"-h\" for help")

    return options


#-------------#
# Functions   #
#-------------#


def load_statistical_potentials(pdb_obj, pdb_dir, pbm_dir, families, radius_default, potential_file_entry=None, split_potential="s3dc_dd", auto_mode=False, family_potentials=False, pbm_potentials=False,  score_threshold=None, taylor_approach=False, pmf_approach=False, bins_approach=False, known_pdb=False, structural_homologs_by_chain=None, dummy_dir="/tmp", verbose=False):
    """
    This function loads a specific statistical potential for each {ChainOfProtein}
    in a {PDB}.

    @input:
    pdb_obj {PDB}
    pdb_dir {string}
    pbm_dir {string}
    families {dict} contains each PDB chain family
    split_potential {string}
    auto_mode {boolean} if yes, it automatically selects the used potential
    family_potentials {boolean} if yes, it uses family potentials
    pbm_potentials {boolean} if yes, it uses PBM + PDB derived potentials
    score_threshold {float} threshold on the scaled score to consider positive k-mers
    taylor_approach {boolean} if yes, it uses potentials approached by Taylor  
    pmf_approach {boolean} if yes, it uses scores instead of zscores
    bins_approach  {boolean} if yes, it uses samples by bins of ditances instead of accumulated samples
    known_pdb {boolean} if yes, it uses the id of pdb_obj to find the potential
    structural_homologs_by_chain {dict} contains the 'list of structural homologs' (value) per 'chain of the pdb_obj' (key)
    dummy_dir {string}
    verbose {boolean}

    @return:
    potentials {dict}

    """

    # Initialize #
    potentials = {}
    thresholds = {}
    radii={}
    min_tm_score = float(config.get("Parameters", "min_tm_score"))
    if pmf_approach:
       pmf=".pmf"
    else:
       pmf=""
    # For protein chain... #
    for protein_chain_obj in pdb_obj.proteins:
      if potential_file_entry is None:
        if verbose:sys.stdout.write("\t\t-- Search potential for chain %s ... \n"%protein_chain_obj.chain)
        family = None
        pdb_chain = None
        # If automatic mode or family potentials are selected... #
        if auto_mode or family_potentials:
            # Initialize #
            structural_homologs = []
            protein_chain_pdb_obj = PDB()
            # Get current working directory #
            cwd = os.getcwd()
            # Create tmp directory #
            tmp = os.path.join(dummy_dir, str(os.getpid()))
            if not os.path.exists(tmp): os.makedirs(tmp)
            # Change directory #
            os.chdir(tmp)
            # Create PDB chain file #
            protein_chain_pdb_obj.add_chain(protein_chain_obj)
            protein_chain_pdb_obj.write(os.path.join(tmp, "pdb_chain.pdb"), force=True)
            # For each PDB chain file... #
            check_close=True
            if known_pdb:
             pdb_name=pdb_obj.id.lower()[0:4]+"_"+protein_chain_obj.chain
             fold_file=os.path.join(pdb_dir, "folds",pdb_name+".txt")
             if verbose:sys.stdout.write("\t\t-- Read known %s ...\n"%pdb_name)
             if os.path.exists(fold_file):
                list_structures=open(fold_file,"r")
                for struct_score in list_structures:
                    if struct_score.startswith("#"):continue
                    structural_homologs.append(tuple(struct_score.strip().split(";")))
                list_structures.close()
                check_close=False
            if structural_homologs_by_chain is not None:
                if structural_homologs_by_chain.has_key(protein_chain_obj.chain):
                    structural_homologs=structural_homologs_by_chain.get(protein_chain_obj.chain)
                    if verbose:sys.stdout.write("\t\t-- Use known homologs for chain %s ...\n"%protein_chain_obj.chain)
                    check_close=False
            if check_close:
             if verbose:sys.stdout.write("\t\t-- Check PDB structural homologs for chain %s ...\n"%protein_chain_obj.chain)
             for pdb_file in sorted(os.listdir(os.path.join(pdb_dir, "split"))):
                # Skip if not PDB chain file #
                if not pdb_file.endswith(".pdb"): continue
                if "dna" in pdb_file: continue
                # Skip if does not interact with DNA #
                triads_file = os.path.join(pdb_dir, "triads", pdb_file[:-4] + ".txt")
                if not os.path.exists(triads_file): continue
                # Superimpose PDB file over PDB chain #
                tmalign_obj = tmalign.get_tmalign_obj(os.path.join(tmp, "pdb_chain.pdb"), os.path.join(pdb_dir, "split", pdb_file))
                if min_tm_score <= max(tmalign_obj.get_tm_scores()):
                    structural_homologs.append((pdb_file[:-4], max(tmalign_obj.get_tm_scores())))
             if structural_homologs_by_chain is not None and not structural_homologs_by_chain.has_key(protein_chain_obj.chain) and structural_homologs != []:
                structural_homologs_by_chain.setdefault(protein_chain_obj.chain,structural_homologs)
             if structural_homologs_by_chain is None and structural_homologs != []:
                structural_homologs_by_chain={}
                structural_homologs_by_chain.setdefault(protein_chain_obj.chain,structural_homologs)
            # Sort homologs by tm-score #
            structural_homologs.sort(key=lambda x: x[-1], reverse=True)
            # Get PDB chain #
            if len(structural_homologs)>0:
               family_potentials_found=False
               for check in structural_homologs:
                 if not family_potentials_found:
                   pdb_chain = check[0]
                   if pdb_chain in families:
                     family = families[pdb_chain]
                   else:
                     family = None
                   if family  is not None and family != "Unknown" and family != "Undefined" :
                     try:
                      family_defined=config.get("Potentials", family)
                      if verbose:sys.stdout.write("\t\t-- Family %s has potential %s \n"%(family,family_defined))
                      family_potentials_found=True
                     except:
                      if verbose:sys.stdout.write("\t\t-- Family %s has not defined potential\n"%(family))
                      family_potentials_found=False
            # Check family_potentials are accessible
            if pdb_chain not in families and family_potentials:
               sys.stdout.write("\t\tfamily potentials cannot be used...\n")
               family_potentials_found=False
            # If structural homolog in families... #
            if pdb_chain in families:
                # Get family #
                family = families[pdb_chain]
                # Set family to None if family is Unknown #
                if (family == "unknown") or (family == "Unknown"): family = None
                if (family == "undefined") or (family == "Undefined"): family = None
            # Return to original directory #
            os.chdir(cwd)
        # If auto mode is selected... #
        if auto_mode and pbm_dir is not None:
            # Initialize #
            potential = config.get("Potentials", "general").split(",")
            # If structural homolog was found and its family has been benchmarked... #
            if pdb_chain is not None and family is not None:
              try:
               potential = config.get("Potentials", family).split(",")
              except:
               if verbose:sys.stdout.write("\t\t-- Family %s has undefined potentials, continue with general\n"%family)
            # If potentials derived from both PDB and PBM... #
            if potential[0] == "pdb": potential_file = os.path.join(pdb_dir, "potentials")
            # Else... #
            else: potential_file = os.path.join(pbm_dir, "potentials")
            # If general potentials... #
            if potential[1] == "general" or pdb_chain is None: potential_file = os.path.join(potential_file, "general")
            # Else... #
            else: potential_file = os.path.join(potential_file, pdb_chain)
            # Add "pmf" type
            potential_file+=pmf
            # If potentials approached by Taylor... #
            if potential[2] == "taylor": potential_file += ".taylor"
            # If potentials computed by bins... #
            if potential[3] == "bins": 
               potential_file += ".bins"
            else:
               potential_file += ".acc"
            # radius for contacts
            if radius_default >0:
                radii.setdefault(protein_chain_obj.chain,radius_default)
            else:
                radii.setdefault(protein_chain_obj.chain,float(potential[-1]))
            # Set threshold #
            thresholds.setdefault(protein_chain_obj.chain, float(potential[-2]))
       # Else... #
        else:
            # If potentials derived from both PDB and PBM are selected... #
            if pbm_potentials and pbm_dir is not None: potential_file = os.path.join(pbm_dir, "potentials")
            # Else... #
            else: potential_file = os.path.join(pdb_dir, "potentials")
            if family_potentials and family_potentials_found:
              # If structural homolog was found... #
                if pdb_chain is not None: potential_file = os.path.join(potential_file, pdb_chain)
                # No structural homolog was found... #
                else: raise ValueError("Could not use family potentials for chain %s. Please unselect option \"-f\" or select option \"-a\" instead" % protein_chain_obj.chain)
            # Else... #
            else: potential_file = os.path.join(potential_file, "general")
            # Add "pmf" type
            potential_file+=pmf
            # If potentials approached by Taylor are selected... #
            if taylor_approach: potential_file += ".taylor"
            # If potentials computed by bins... #
            if bins_approach: 
               potential_file += ".bins"
            else:
               potential_file += ".acc"
            # radius for contacts
            radii.setdefault(protein_chain_obj.chain,radius_default)
            # Set threshold #
            if score_threshold is None: thresholds.setdefault(protein_chain_obj.chain, float(config.get("Potentials", "default_score_threshold")))
            else: thresholds.setdefault(protein_chain_obj.chain, score_threshold)
        # Load potential file #
        potentials[protein_chain_obj.chain] = spotentials.Potentials(potential_file + ".txt", split_potential)
        if verbose:sys.stdout.write("\t\t-- Load %s  ...\n"%potential_file);
        if verbose:sys.stdout.write("\t\t-- Threshold MSA %f ...\n"%thresholds[protein_chain_obj.chain]);
      else:
        # Load potential file #
        potentials[protein_chain_obj.chain] = spotentials.Potentials(potential_file_entry, split_potential)
        if verbose:sys.stdout.write("\t\t-- Load %s  ...\n"%potential_file_entry);
        # Set threshold #
        if score_threshold is None: thresholds.setdefault(protein_chain_obj.chain, float(config.get("Potentials", "default_score_threshold")))
        else: thresholds.setdefault(protein_chain_obj.chain, score_threshold)
        if verbose:sys.stdout.write("\t\t-- Threshold MSA %f ...\n"%thresholds[protein_chain_obj.chain]);
        radii.setdefault(protein_chain_obj.chain,radius_default)

    return potentials, thresholds, radii, structural_homologs_by_chain

def get_msa_obj(triads_obj, x3dna_obj, potentials, radii, fragment_restrict, binding_restrict, split_potential, thresholds, methylation=False):

    # Initialize #
    msa_obj = nMSA()
    binding_sites = {}
    all_kmers_scaled_scores = {}
    # extend methylated nucleotide residues
    if methylation:
       msa_obj.extend_residues(list("XJOQ"))

    # For each PDB chain... #
    for pdb_chain in sorted(potentials):
        # Get dinucleotide raw scores and binding site region #
        scores, binding_site = get_scores_and_binding_site(triads_obj, x3dna_obj, potentials, radii, fragment_restrict, binding_restrict, split_potential, pdb_chain,methylation)
        # Get k-mers scaled scores #
        if len(binding_site)<=0: continue
        all_kmers_scaled_scores.setdefault(pdb_chain, get_kmers_scaled_scores(scores, binding_site, split_potential,methylation))
        # Add binding site to binding sites #
        binding_sites.setdefault(pdb_chain, binding_site)
    # Get whole binding site #
    binding_site = set()
    # For each PDB chain... #
    for pdb_chain in sorted(binding_sites):
        # For each site... #
        for i in binding_sites[pdb_chain].keys():
            # Add site to binding site #
            binding_site.add(i)
    # Sort binding site #
    binding_site = sorted(binding_site)
    if len(binding_site)>0:
      # Add last basepair to binding site#
      binding_site.append(binding_site[-1] + 1)
      # Add binding site length #
      msa_obj.set_binding_site_length(len(binding_site))
      # # Add original PDB sequence to MSA #
      # msa_obj.add_sequence(x3dna_obj.get_nucleotide_sequence(binding_site[0], binding_site[-1]), "wt")
      # Inialize #
      default_sequence = "N" * (binding_site[-1] - binding_site[0] + 1)
      # For each PDB chain... #
      for pdb_chain in sorted(all_kmers_scaled_scores):
        # For each k-mer... #
        for kmer, start, end in all_kmers_scaled_scores[pdb_chain]:
            #print "score",float(all_kmers_scaled_scores[pdb_chain][(kmer, start, end)]), "thr",float(thresholds[pdb_chain]),"kmer",kmer,"chain",pdb_chain
            if float(all_kmers_scaled_scores[pdb_chain][(kmer, start, end)]) >= float(thresholds[pdb_chain]):
                sequence = default_sequence[:start - binding_site[0]] + kmer + default_sequence[end - binding_site[0] + 1:]
                msa_obj.add_sequence(sequence, all_kmers_scaled_scores[pdb_chain][(kmer, start, end)])

    # If the msa is empty, create a non-informative pwm #
    if len(msa_obj.get_sequences()) == 0 and len(msa_obj.get_pwm()) > 0:
        msa_obj.set_sequences()
    if len(msa_obj.get_sequences()) == 0:
        #Setup as a non-informative  matrix
        msa_obj.add_sequence("A"*len(binding_site), 0.1) 
        msa_obj.add_sequence("C"*len(binding_site), 0.1) 
        msa_obj.add_sequence("G"*len(binding_site), 0.1) 
        msa_obj.add_sequence("T"*len(binding_site), 0.1) 
        if methylation:
           msa_obj.add_sequence("X"*len(binding_site), 0.1)
           msa_obj.add_sequence("J"*len(binding_site), 0.1)
           msa_obj.add_sequence("O"*len(binding_site), 0.1)
           msa_obj.add_sequence("Q"*len(binding_site), 0.1)
        msa_obj.add_sequence("N"*len(binding_site), 1.0)

    # Check the size of the multiple alignment
    if len(msa_obj.get_sequences()) <1000:
       repeat=int(1000.0/len(msa_obj.get_sequences()))
       original_set=set()
       for sequence,score in  msa_obj.get_sequences():
          original_set.add((sequence,score))
       for i in xrange(repeat):
          for sequence,score in original_set:
              msa_obj.add_sequence(sequence, float(score-0.01/float(i+1))) 

    # Create a PWM inside the msa_obj #
    msa_obj.set_pwm()

    return msa_obj

def get_scores_and_binding_site(triads_obj, x3dna_obj, potentials, radii, fragment_restrict, binding_restrict, split_potential, pdb_chain,methylation=False):
    """
    """

    # Initialize #
    scores = {}
    binding_site = {}
    nucleotides = list("ACGT")
    if methylation:
       nucleotides.extend(list("XJOQ"))
    radius = 0
    # Get original DNA sequence indexes
    basepairs=x3dna_obj.get_basepairs()
    dna_idx={}
    for basepair in basepairs.iterkeys():
        (fwd_pdb_chain, fwd_residue_num), (rev_pdb_chain, rev_residue_num) = basepairs.get(basepair)
        dna_idx.setdefault((fwd_pdb_chain,fwd_residue_num),int(basepair))
    #print "DNA_IDX",dna_idx
    if radii.has_key(pdb_chain): radius=radii.get(pdb_chain)
    if radius <= 0: radius=float(config.get("Parameters", "max_contact_distance"))
    # For each contact... #
    for triad_obj in triads_obj.get_triads():
        # Initialize #
        a_oa, b_ob, distance, residue_A, residue_B = triad_obj.return_as_string().split(";")
        chain, residue_num = residue_A.split("-")
        dinucleotide = get_triad_dinucleotide(triad_obj, x3dna_obj)
        # Skip if not the right chain... #
        if chain != pdb_chain: continue
        # Skip if far contact... #
        dab = numpy.floor(float(distance))
        if dab > radius: continue
        # Skip if amino-acid is not in the restricted fragment... #
        #print "CHECK",dinucleotide, a_oa, b_ob, distance, residue_A, residue_B
        if fragment_restrict is not None:
           belongs=False
           if fragment_restrict.has_key(chain):  
               for interval in fragment_restrict.get(chain):
                   if int(residue_num) >= int(interval[0]) and int(residue_num) <= int(interval[1]): belongs=True
                   if belongs: break
               if not belongs: continue
        # Skip if nucleotide is not in the restricted binding... #get_scores_and_binding_site
        bp_1F,bp_1R,bp_2F,bp_2R=residue_B.split(",")
        # Add dinucleotide to binding site #
        chain_dna,base_number_1=bp_1F.split("-")
        chain_dna,base_number_2=bp_2F.split("-")
        #print "NUCLEOTIDES", bp_1F,bp_1R,bp_2F,bp_2R,"CHAIN",chain_dna,"BASES",base_number_1,base_number_2
        if not dna_idx.has_key((chain_dna,int(base_number_1))): continue
        if not dna_idx.has_key((chain_dna,int(base_number_2))): continue
        binucleotide = (dna_idx[(chain_dna,int(base_number_1))],dna_idx[(chain_dna,int(base_number_2))])
        if binding_restrict is not None:
               belongs=False
               for interval in binding_restrict:
                   if int(binucleotide[0]) >= interval[0] and int(binucleotide[1]) <= interval[1]: belongs=True
                   #if int(binucleotide[1]) >= interval[0] and int(binucleotide[0]) <= interval[1]: belongs=True
                   if belongs: break
               if not belongs: continue
        # Add dinucleotide to binding site #
        #sys.stdout.write("\t Add %s %s %s %s %s %s \n"%(dinucleotide,a_oa, b_ob, distance, residue_A, residue_B))
        binding_site.setdefault(dinucleotide, []).append(triad_obj)
        # Skip if environment could not be obtained for amino acid #
        if "None" in a_oa: continue
        # Skip if environment could not be obtained for dinucleotide #
        if "None" in b_ob: continue
        # The following arguments arguments are defined as in the papers #
        dab = numpy.floor(float(distance))
        a, hydrophobicity, degree_of_exposure, secondary_structure = a_oa.split("-")
        oa = "%s-%s-%s" % (hydrophobicity, degree_of_exposure, secondary_structure)
        b, nitrogenous_bases, dna_strand, dna_groove, dna_chemical_group = b_ob.split("-")
        # For each dinucleotide step... #
        for dinucleotide_sequence in itertools.product(nucleotides, repeat=2):
            # Initialize #
            dinucleotide_sequence = "".join(dinucleotide_sequence)
            if "O" in dinucleotide_sequence and ("X" in dinucleotide_sequence or "J" in dinucleotide_sequence): continue
            if "Q" in dinucleotide_sequence and ("X" in dinucleotide_sequence or "J" in dinucleotide_sequence): continue
            if "X" in dinucleotide_sequence and ("O" in dinucleotide_sequence or "Q" in dinucleotide_sequence): continue
            if "J" in dinucleotide_sequence and ("O" in dinucleotide_sequence or "Q" in dinucleotide_sequence): continue
            scores.setdefault((dinucleotide_sequence, dinucleotide), 0.0)
            # Get dinucleotide object #
            dinucleotide_obj = triads.DinucleotideEnvironment(dinucleotide_sequence, dna_strand, dna_groove, dna_chemical_group)
            # The following arguments arguments are defined as in the papers #
            b_ob = dinucleotide_obj.return_as_string()
            b, nitrogenous_bases, dna_strand, dna_groove, dna_chemical_group = b_ob.split("-")
            ob = "%s-%s-%s-%s" % (nitrogenous_bases, dna_strand, dna_groove, dna_chemical_group)
            a_b = "%s;%s" % (a, b)
            a_b_oa_ob = "%s;%s" % (a_oa, b_ob)
            oa_ob = "%s;%s" % (oa, ob)
            # Apply potential... #
            if split_potential == "3d":
                if potentials[chain].get_score(split_potential, distance=dab) is not None:
                    scores[(dinucleotide_sequence, dinucleotide)] += potentials[chain].get_score(split_potential, distance=dab)
            if split_potential == "3dc":
                if potentials[chain].get_score(split_potential, key=oa_ob, distance=dab) is not None:
                    scores[(dinucleotide_sequence, dinucleotide)] += potentials[chain].get_score(split_potential, key=oa_ob, distance=dab)
            if split_potential == "s3dc":
                if potentials[chain].get_score(split_potential, key=a_b_oa_ob, distance=dab) is not None:
                    scores[(dinucleotide_sequence, dinucleotide)] += potentials[chain].get_score(split_potential, key=a_b_oa_ob, distance=dab)
            if split_potential == "s3dc_dd":
                if potentials[chain].get_score(split_potential, key=a_b_oa_ob, distance=dab) is not None:
                    scores[(dinucleotide_sequence, dinucleotide)] += potentials[chain].get_score(split_potential, key=a_b_oa_ob, distance=dab)
            if split_potential == "s3dc_di":
                if potentials[chain].get_score(split_potential, key=a_b_oa_ob) is not None:
                    scores[(dinucleotide_sequence, dinucleotide)] += potentials[chain].get_score(split_potential, key=a_b_oa_ob)
            if split_potential == "pair":
                if potentials[chain].get_score(split_potential, key=a_b, distance=dab) is not None:
                    scores[(dinucleotide_sequence, dinucleotide)] += potentials[chain].get_score(split_potential, key=a_b, distance=dab)
    return scores, binding_site

def get_triad_dinucleotide(triad_obj, x3dna_obj):
    """
    This function returns the dinucleotide involved in the input {Triad}.

    @input:
    triad_obj {Triad}
    x3dna_obj {X3DNA}

    @return {int}

    """

    # Initalize #
    dinucleotides = []
    residues = triad_obj._residue_B.split(",")

    # For each residue... #
    for residue in residues:
        # Initialize #
        pdb_chain, residue_num = residue.split("-")
        # Get basepair #
        basepair = x3dna_obj.get_residue_basepair(pdb_chain, int(residue_num))
        # For each dinucleotide... #
        for dinucleotide in x3dna_obj.get_basepair_dinucleotides(basepair):
            dinucleotides.append(dinucleotide)
    # Get most common dinucleotide #
    most_common_dinucleotide = Counter(dinucleotides).most_common(1)

    return most_common_dinucleotide[0][0]

def get_kmers_scaled_scores(scores, binding_site, split_potential,methylation=False):
    """
    """

    # Initialize #
    binding_site_basepairs = set()
    all_kmers = {}
    all_kmers_scaled_scores = {}
    nucleotides = list("ACGT")
    if methylation:
       nucleotides.extend(list("XJOQ"))
    max_scores = {}
    min_scores = {}
    min_motif_size=int(config.get("Parameters", "min_motif_size"))
    min_basepairs=int(config.get("Parameters", "min_basepairs"))
    # For each dinucleotide... #
    for dinucleotide in binding_site:
        binding_site_basepairs.add(dinucleotide)
        # A dinucleotide involve two basepairs #
        binding_site_basepairs.add(dinucleotide + 1)
    # Convert binding site to a list #
    binding_site_basepairs = sorted(binding_site_basepairs) 
    # Get default sequence #
    default_sequence = "N" * (binding_site_basepairs[-1] - binding_site_basepairs[0] + 1)
    if len(binding_site_basepairs) < max(min_basepairs,min_motif_size):
       min_length=len(binding_site_basepairs)
    else:
       min_length=max(min_basepairs,min_motif_size)
    round=0
    # For each basepair... #
    for basepair in binding_site_basepairs:
        start = basepair
        end = start + min_length - 1
        if methylation:
           kmers = dict(zip(nucleotides, [0, 0, 0, 0, 0, 0, 0, 0]))
        else:
           kmers = dict(zip(nucleotides, [0, 0, 0, 0]))
        if end > binding_site_basepairs[-1]:
            if round < 1:
              sys.stdout.write("Binding site too short (%d nucleotides): reduce 'min_motif_size' in config.ini \n"%(len(binding_site_basepairs)))
            break
        round = round + 1
        for dinucleotide in range(start, end):
            tmp_scores = {}
            # For each k-mer... #
            for kmer in kmers:
                # For each nucleotide:
                for nucleotide in nucleotides:
                    dinucleotide_sequence = kmer[-1] + nucleotide
                    if methylation:
                       if "O" in kmer and (nucleotide == "X" or nucleotide == "J"): continue
                       if "Q" in kmer and (nucleotide == "X" or nucleotide == "J"): continue
                       if "X" in kmer and (nucleotide == "O" or nucleotide == "Q"): continue
                       if "J" in kmer and (nucleotide == "O" or nucleotide == "Q"): continue
                       methylC = 0
                       if Counter(kmer).get("X") is not None: methylC = methylC + Counter(kmer).get("X")
                       if Counter(kmer).get("J") is not None: methylC = methylC + Counter(kmer).get("J")
                       if Counter(kmer).get("O") is not None: methylC = methylC + Counter(kmer).get("O")
                       if Counter(kmer).get("Q") is not None: methylC = methylC + Counter(kmer).get("Q")
                       if methylC > 5: continue
                    if (dinucleotide_sequence, dinucleotide) in scores:
                        tmp_scores[kmer + nucleotide] = kmers[kmer] + scores[(dinucleotide_sequence, dinucleotide)]
                    else:
                        tmp_scores[kmer + nucleotide] = kmers[kmer] + 0
            kmers = tmp_scores.copy()
        for kmer in kmers:
            all_kmers[(kmer, start, end)] = kmers[kmer]
    # If split potential is either "S3DC" or "S3DCdd"... #
    if split_potential == "s3dc" or split_potential == "s3dc_dd" or split_potential == "pair":
        # For each k-mer ... #
        for key in all_kmers:
            all_kmers[key] *= -1
    # For each k-mer ... #
    for kmer, start, end in all_kmers:
        # Get max. and min. scores for scaling #
        max_scores.setdefault((start, end), all_kmers[(kmer, start, end)])
        min_scores.setdefault((start, end), all_kmers[(kmer, start, end)])
        if all_kmers[(kmer, start, end)] > max_scores[(start, end)]:
            max_scores[(start, end)] = all_kmers[(kmer, start, end)]
        if all_kmers[(kmer, start, end)] < min_scores[(start, end)]:
            min_scores[(start, end)] = all_kmers[(kmer, start, end)]
    # For each k-mer ... #
    for kmer, start, end in all_kmers:
        all_kmers_scaled_scores.setdefault((kmer, start, end), 0)
        if max_scores[(start, end)] != min_scores[(start, end)]:
            all_kmers_scaled_scores[(kmer, start, end)] = scale(all_kmers[(kmer, start, end)], max_score=max_scores[(start, end)], min_score=min_scores[(start, end)])
    
    return all_kmers_scaled_scores

def scale(score, max_score, min_score, max_scaled_score=1.0, min_scaled_score=0.0):
    """
    This function scales a score in a range x0..x1 to a new range y0..y1:

    @input:
    score {float}
    max_score {float}
    min_score {float}
    max_scaled_score {float}
    min_scaled_score {float}

    @return {float}
    
    """
    if (max_score - min_score) > 1.0e-6:
      return min_scaled_score + (max_scaled_score - min_scaled_score) * (score - min_score) / (max_score - min_score)
    else:
      return 0.0

def write_pwm_by_meme(msa_obj,output_file,dummy_dir="/tmp"):
    src_path=config.get("Paths","src_path")
    meme_path = os.path.join(src_path, config.get("Paths", "meme_path"))
    label=output_file.split("/")[-1]
    dummy_file= os.path.join(dummy_dir,label+".pwm")
    if len(msa_obj.get_pwm())<0:
       msa_obj.set_pwm()
    if os.path.exists(dummy_file): os.remove(dummy_file)
    msa_obj.write(dummy_file,option="pwm")
    try: 
      process = subprocess.check_output([os.path.join(meme_path, "uniprobe2meme"), dummy_file], stderr=subprocess.PIPE)
      # Create output file #
      out = open(output_file, "wt")
      out.write(process)
      out.close()
      os.remove(dummy_file)
    except:
        # Copy the meme format #
        msa_obj.write(file_name=output_file, option="meme")
        



def write_logo(msa_obj,logo_file,dummy_dir="/tmp"):
    #Initialize programs
    src_path=config.get("Paths","src_path")
    meme_path = os.path.join(src_path, config.get("Paths", "meme_path"))
    weblogo_path = os.path.join(src_path, config.get("Paths", "weblogo_path"))
    ghostscript_path = os.path.join(src_path, config.get("Paths", "ghostscript_path"))
    python = os.path.join(config.get("Paths", "python_path"), "python")
    os.environ['COREBIOPATH'] = ghostscript_path
    # Fill dummy sequences if necessary
    if len(msa_obj.get_sequences()) <=0:
       msa_obj.set_sequences()
    # Define files
    fasta_file   = os.path.join(dummy_dir,logo_file.split("/")[-1]+".msa.fa")
    forward_file = logo_file+".fwd.png"
    reverse_file = logo_file+".rev.png"
    # Create FASTA file #
    if not os.path.exists(fasta_file):
     out = open(fasta_file, "wt")
     count = 0
     # For each sequence... #
     for sequence, score in msa_obj.get_sequences():
      count += 1
      out.write(">%d %s (%.3f)\n%s\n" % (count, msa_obj.get_motif(), score, sequence))
     out.close()
    #Create weblogo
    try:
      process = subprocess.check_output([python, os.path.join(weblogo_path, "weblogo"), "-c", "classic", "-F", "png", "-s", "large", "--scale-width=no", "-f", fasta_file, "-o", forward_file], env=os.environ, stderr=subprocess.PIPE)
      process = subprocess.check_output([python, os.path.join(weblogo_path, "weblogo"), "-c", "classic", "-F", "png", "--revcomp", "-s", "large", "--scale-width=no", "-f", fasta_file, "-o", reverse_file], env=os.environ, stderr=subprocess.PIPE)
    except:
     try:
         os.system("%s %s -c classic  -F png -s large --scale-width=no -f %s -o %s"%(python,os.path.join(weblogo_path, "weblogo"),fasta_file,forward_file))
         os.system("%s %s -c classic  -F png --revcomp -s large --scale-width=no -f %s -o %s"%(python,os.path.join(weblogo_path, "weblogo"),fasta_file,reverse_file))
     except:
         sys.stdout.write("WARNING: \"weblogo\" execution failed!")
    # Remove dummy FASTA file
    if os.path.exists(fasta_file): os.remove(fasta_file)

def write_protein_logo(msa_obj,logo_file, dummy_dir="/tmp"):
    #Initialize programs
    src_path=config.get("Paths","src_path")
    meme_path = os.path.join(src_path, config.get("Paths", "meme_path"))
    weblogo_path = os.path.join(src_path, config.get("Paths", "weblogo_path"))
    ghostscript_path = os.path.join(src_path, config.get("Paths", "ghostscript_path"))
    python = os.path.join(config.get("Paths", "python_path"), "python")
    os.environ['COREBIOPATH'] = ghostscript_path
    # Fill dummy sequences if necessary
    if len(msa_obj.get_sequences()) <=0:
       msa_obj.set_sequences()
    # Define files
    fasta_file   = os.path.join(dummy_dir,logo_file.split("/")[-1]+".msa.fa")
    protein_file = logo_file+".png"

    # Create FASTA file #
    if not os.path.exists(fasta_file):
     out = open(fasta_file, "wt")
     count = 0
     # For each sequence... #
     for sequence, score in msa_obj.get_sequences():
      count += 1
      out.write(">%d %s (%.3f)\n%s\n" % (count, msa_obj.get_motif(), score, sequence))
     out.close()
    #Create weblogo
    try:
      process = subprocess.check_output([python, os.path.join(weblogo_path, "weblogo"), "-A", "protein", "-a", "ARNDCEQGHILKMFPSTWYV", "-c", "chemistry", "-F", "png_print", "-s", "large", "--scale-width=no", "--errorbars=no", "-f", fasta_file, "-o", protein_file], env=os.environ, stderr=subprocess.PIPE)
    except:
     try:
         os.system("%s %s -A protein -a ARNDCEQGHILKMFPSTWYV  -c chemistry  -F png_print -s large --scale-width=no --errorbars no --scale-width=no -f %s -o %s"%(python,os.path.join(weblogo_path, "weblogo"),fasta_file,protein_file))
     except:
         sys.stdout.write("WARNING: \"weblogo\" execution failed!")
    # Remove dummy FASTA file
    os.remove(fasta_file)




def get_score_for_subseq(binding_site, x3dna_obj, dna_sequence, potentials,  split_potential):
    # Initialize #
    final_score = 0.0
    #scores_per_nucleotide = {}
    dn=0
    dna_numbers={}
    # Get original DNA sequence indexes
    for dinucleotide in binding_site.iterkeys():
        if x3dna_obj.has_dinucleotide(dinucleotide):
           bp1,bp2=x3dna_obj.get_dinucleotide(dinucleotide)
           if x3dna_obj.has_basepair(bp1) and x3dna_obj.has_basepair(bp2):
              ((fwd_pdb_chain1, fwd_base_num1), (rev_pdb_chain1, rev_base_num1)) = x3dna_obj.get_basepair(bp1)
              ((fwd_pdb_chain2, fwd_base_num2), (rev_pdb_chain2, rev_base_num2)) = x3dna_obj.get_basepair(bp2)
              dn = dn+1
              dna_numbers.setdefault(dn,(fwd_pdb_chain1,fwd_base_num1,fwd_pdb_chain2,fwd_base_num2,dinucleotide))

    
    # Iterate the protein chains in the potential #
    for pc in potentials.keys():
        potential = potentials[pc]
        # Iterate the DNA sequence #
        for dn in dna_numbers.iterkeys():
            #scores_per_nucleotide.setdefault(dn, 0.0)
            fwd_chain1,fwd_base_num1,fwd_chain2, fwd_base_num2,dinucleotide = dna_numbers.get(dn)
            # Iterate the triads #
            for triad_obj in binding_site.get(dinucleotide):
                # Initialize #
                a_oa, b_ob, distance, residue_A, residue_B = triad_obj.return_as_string().split(";")
                # Skip if environment could not be obtained for amino acid #
                if "None" in a_oa: continue
                # Skip if environment could not be obtained for dinucleotide #
                if "None" in b_ob: continue
                bp_1F,bp_1R,bp_2F,bp_2R=residue_B.split(",")
                # Get dinucleotide of binding site #
                chain_dna,base_number_1=bp_1F.split("-")
                chain_dna,base_number_2=bp_2F.split("-")
                protein = a_oa
                dna = b_ob
                chain, residue_pos = residue_A.split("-")
                residue_num=int(residue_pos)
                # The chain in the triad and the chain in the potentials object must be coincident #
                if pc != chain: continue
                # The base number and the position in the DNA sequence must be coincident #
                if str(base_number_1) == str(fwd_base_num1) and chain_dna == str(fwd_chain1) and str(base_number_2) == str(fwd_base_num2) and chain_dna == str(fwd_chain2):
                        # The following arguments arguments are defined as in the papers #
                        # Set the distance #
                        dab = numpy.floor(float(distance))
                        # Set properties of the protein #
                        a, hydrophobicity, degree_of_exposure, secondary_structure = a_oa.split("-")
                        aminoacid_three_letters = a
                        oa = "%s-%s-%s" % (hydrophobicity, degree_of_exposure, secondary_structure)
                        # Set and modify properties of the DNA #
                        b, nitrogenous_bases, dna_strand, dna_groove, dna_chemical_group = b_ob.split("-")
                        dinucleotide_sequence = dna_sequence[dn-1] + dna_sequence[dn]
                        nucleotide_environment = triads.DinucleotideEnvironment(dinucleotide_sequence, dna_strand, dna_groove, dna_chemical_group)
                        nitrogenous_bases = nucleotide_environment.get_nitrogenous_bases()
                        b = dinucleotide_sequence
                        ob = "%s-%s-%s-%s" % (nitrogenous_bases, dna_strand, dna_groove, dna_chemical_group)
                        dna = dinucleotide_sequence + "-" + nitrogenous_bases + "-" + dna_strand + "-" + dna_groove + "-" + dna_chemical_group
                        a_b = "%s;%s" % (a, b)
                        oa_ob = "%s;%s" % (oa, ob)
                        a_b_oa_ob = "%s;%s" % (protein, dna)
                        # Use different data sources to get the scores depending on the split potentials #
                        if split_potential == "3d":
                            score = potential.get_score(split_potential, distance=dab)
                        if split_potential == "3dc":
                            score = potential.get_score(split_potential, key=oa_ob, distance=dab)
                        if split_potential == "s3dc":
                            score = potential.get_score(split_potential, key=a_b_oa_ob, distance=dab)
                        if split_potential == "s3dc_dd":
                            score = potential.get_score(split_potential, key=a_b_oa_ob, distance=dab)
                        if split_potential == "s3dc_di":
                            score = potential.get_score(split_potential, key=a_b_oa_ob)
                        if split_potential == "pair":
                            score = potential.get_score(split_potential, key=a_b, distance=dab)
                        if split_potential == "local":
                            score = potential.get_score(split_potential, key=a_oa)
                        if score is not None:
                            final_score += score
    if split_potential=="pair" or split_potential=="local" or split_potential=="s3dc" or split_potential=="s3dc_dd":
       final_score = -final_score
    return final_score



def get_min_max_scores(msa_obj, binding_site,  x3dna_obj, potentials,  split_potential, dummy_dir,methylation=False):

    nucleotides=["A","C","G","T"]
    #extend for methylation
    if methylation:
       nucleotides.extend(list("XJOQ"))
    for sc in ["best_max", "worse_min"]:
        if sc == "best_max":
            reverse = True
        elif sc == "worse_min":
            reverse = False
        good_sequences = set()
        pwm = msa_obj.get_pwm()
        best_seq = ""
        for position in pwm:
            check=zip(nucleotides,position)
            best_seq += sorted(check, key=lambda x:x[1], reverse=reverse)[0][0]
        good_sequences.add(best_seq)
        # Get variations of the best sequence # 
        counter = 0
        while len(good_sequences) < max(100,4*len(pwm)):
         counter += 1
         if counter == 1000:
            break
         test_sequences=set()
         for seq in good_sequences:
             test_sequences.add(seq)
         for check_seq in test_sequences:
          for delta in [x * 0.01 for x in range(5, 20, 5)]:
            for i in range(0, len(pwm)):
                position = pwm[i]
                check=zip(nucleotides,position)
                selected_freq = float(sorted(check, key=lambda x:x[1], reverse=reverse)[0][1])
                selected_res  = sorted(check, key=lambda x:x[1], reverse=reverse)[0][0]
                for k in range(1, 4):
                    try_freq  = float(sorted(check, key=lambda x:x[1], reverse=reverse)[k][1])
                    try_res   = sorted(check, key=lambda x:x[1], reverse=reverse)[k][0]
                    if abs(selected_freq - try_freq) <= delta:
                        new_seq    = list(check_seq)
                        new_seq[i] = try_res
                        new_seq = "".join(new_seq)
                        good_sequences.add(new_seq)
            if len(good_sequences) > max(100,4*len(pwm)):
                break

        # score each one of the good sequences #
        top_score = None
        for seq in good_sequences:
            #subseq_score, scores_per_nucleotide = get_score_for_subseq(binding_site,  x3dna_obj, seq, potentials,  split_potential)
            subseq_score = get_score_for_subseq(binding_site,  x3dna_obj, seq, potentials,  split_potential)
            if top_score == None:
                top_score = subseq_score
            else:
                if sc == "best_max":
                    if top_score < subseq_score:
                        top_score = subseq_score

                elif sc == "worse_min":
                    if top_score > subseq_score:
                        top_score = subseq_score
                 
        if sc == "best_max":
            max_score = top_score            
        elif sc == "worse_min":
            min_score = top_score
        

    return min_score, max_score





def get_single_pwm(input_pdb_file,input_threading_file,threading,output_file,pbm_dir,pdb_dir,families,potential_file, radius, fragment_restrict, binding_restrict, split_potential,auto_mode,family_potentials,pbm_potentials,score_threshold,taylor_approach,pmf,bins_approach,known,meme,reset,dummy_dir,verbose,methylation=False):

    # Initialize output/input files
    pwm_file = output_file+".pwm"
    pwm_meme = output_file+".meme"
    meme_file= output_file+".meme.s"
    msa_file = output_file+".msa"
    logo_file= output_file+".logo"


    if os.path.exists(msa_file):

       msa_obj=nMSA(msa_file, output_file)

    elif os.path.exists(pwm_meme):

       msa_obj=nMSA(pwm_meme,output_file,"meme")

    elif os.path.exists(pwm_file):

       msa_obj=nMSA(pwm_meme,output_file,"pwm")

    else:

      if threading:

         #Get triads, pdb etc from threading file #
         triads_obj, pdb_obj, x3dna_obj = threading_to_triads.threading_triads(threading_file=input_threading_file, pdb_dir= pdb_dir)

      else:
        try:
         # Get PDB object #
         if  verbose:sys.stdout.write("\t--Get protein %s ...\n"% input_pdb_file)
         pdb_obj = PDB( input_pdb_file)

         # Get DSSP object #
         if  verbose:sys.stdout.write("\t\t-- calculate secondary structure ...\n")
         dssp_obj = dssp.get_dssp_obj(os.path.abspath( input_pdb_file))

         # Get X3DNA object #
         if  verbose:sys.stdout.write("\t--Get DNA %s ...\n"% input_pdb_file)
         x3dna_obj = x3dna.get_x3dna_obj(os.path.abspath( input_pdb_file),dummy_dir=dummy_dir)
         if len(x3dna_obj.get_dinucleotides().keys() )<1:
            sys.stdout.write("Missing DNA ...\n")
        
         # Get contacts object #
         if  verbose:sys.stdout.write("\t\t-- calculate contacts ...\n")
         contacts_obj = contacts.get_contacts_obj(pdb_obj, x3dna_obj)
         if len(contacts_obj.get_contacts())<1:
            sys.stdout.write("Missing Protein-DNA contacts ...\n")

         # Get triads object #
         if  verbose:sys.stdout.write("\t\t-- calculate protein-dna pairs ...\n")
         triads_obj = triads.get_triads_obj(pdb_obj, dssp_obj, x3dna_obj, contacts_obj)

        except Exception as e:
         raise Exception("Failed to get protein features")

      # Load statistical potential #
      if  verbose:sys.stdout.write("\t--Load potentials ...\n")
      try:
        potentials, thresholds , radii, structural_homologs_by_chain = load_statistical_potentials(pdb_obj,  pdb_dir, pbm_dir, families, radius, potential_file, split_potential,  auto_mode,  family_potentials,  pbm_potentials,  score_threshold,  taylor_approach,  pmf, bins_approach, known, None,  os.path.abspath( dummy_dir), verbose)
      except Exception as e:
         raise Exception("Failed to get potentials ")
      
      # Get MSA object #
      if  verbose:sys.stdout.write("\t--Get MSA object ...\n")
      try:
        msa_obj = get_msa_obj(triads_obj, x3dna_obj, potentials, radii, fragment_restrict, binding_restrict, split_potential, thresholds, methylation)
        motif_name= output_file.split("/")[-1]
        msa_obj.set_motif(motif_name)
      except Exception as e:  
        raise Exception("Failed to get MSA ")

    # Check msa_obj contain sequences, otherwise generate a dummy set
    if reset: msa_obj.clean_sequences()
    if len(msa_obj.get_sequences()) == 0 and len(msa_obj.get_pwm()) > 0:
        msa_obj.set_sequences()
 
    # Write PWM #
    if  verbose:sys.stdout.write("\t--Write PWM ...\n")

    if not os.path.exists(pwm_file):
     if  verbose:sys.stdout.write("\t\t--PWM...\n")
     try:
       msa_obj.write(pwm_file, option="pwm")
     except Exception as e:
       if  verbose:sys.stdout.write("Failed %s\n"%e)
    if not os.path.exists(pwm_meme):
     if  verbose:sys.stdout.write("\t\t--PWM in MEME format...\n")
     try:
       msa_obj.write(pwm_meme, option="meme")
     except Exception as e:
       if  verbose:sys.stdout.write("Failed %s\n"%e)
    if not os.path.exists(msa_file):
     if  verbose:sys.stdout.write("\t\t--MSA...\n")
     try:
       msa_obj.write(msa_file, option="msa")
     except Exception as e:
       if  verbose:sys.stdout.write("Failed %s\n"%e)

    if not os.path.exists(meme_file) and  meme:
     if  verbose:sys.stdout.write("\t\t--PWM by MEME...\n")
     try:
       write_pwm_by_meme(msa_obj,meme_file, dummy_dir)
     except Exception as e:
       if  verbose:sys.stdout.write("Failed %s\n"%e)

    if not os.path.exists(logo_file+".fwd.png") or  not os.path.exists(logo_file+".rev.png"):
     if  verbose:sys.stdout.write("\t\t--Logos...\n")
     try:
       write_logo(msa_obj,logo_file, dummy_dir)
     except Exception as e:
       if  verbose:sys.stdout.write("Failed %s\n"%e)
    


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Initialize #
    if not os.path.exists(options.dummy_dir): os.makedirs(options.dummy_dir)
    pbm_dir=None
    if options.pbm_dir is not None: pbm_dir = os.path.abspath(options.pbm_dir)
    pdb_dir=None
    if options.threading:
       input_threading_files=options.input_files
       input_pdb_files=None
       if not input_threading_files.startswith("/"): input_threading_files=os.path.abspath(options.input_files)
    else:
       input_pdb_files=options.input_files
       input_threading_files=None
       if not input_pdb_files.startswith("/"): input_pdb_files=os.path.abspath(options.input_files)
    if options.pdb_dir is not None: pdb_dir = os.path.abspath(options.pdb_dir)
    complete=float(options.complete)

    families = {}
    if options.verbose:sys.stdout.write("Check families...\n")
    if not os.path.exists(os.path.join(pdb_dir, "families.txt")):
       sys.stdout.write("Families file %s is not used\n"%(os.path.join(pdb_dir, "families.txt")))
    else:
      for line in functions.parse_file(os.path.join(options.pdb_dir, "families.txt")):
        if line.startswith("#"): continue
        pdb_chain, family = line.split(";")
        families[pdb_chain] = family

    fragment_restrict=None
    if options.fragment_restrict is not None:
       fragment_restrict={}
       segment_fragment=options.fragment_restrict.split(";")
       for interval in segment_fragment:
           ac,bc = interval.split("-")
           chain1=chain2=" "
           if len(ac.split("_"))==1: a=ac
           if len(bc.split("_"))==1: b=bc
           if len(ac.split("_"))>1: a,chain1=ac.split("_")
           if len(bc.split("_"))>1: b,chain2=bc.split("_")
           if chain1==chain2:
              fragment_restrict.setdefault(chain1,[]).append((int(a),int(b)))

    binding_restrict=None
    if options.binding_restrict is not None:
       binding_restrict=[]
       segment_binding=options.binding_restrict.split(";")
       for interval in segment_binding:
           a,b = interval.split("-")
           aa = int(a)
           bb = int(b)
           binding_restrict.append((aa,bb))
    # Execute for single file PWM
    if os.path.isfile(options.input_files):
       if options.verbose:sys.stdout.write("Generate file-specific PWM %s ...\n"%(options.output_file))
       try:
         get_single_pwm(input_pdb_files,input_threading_files,options.threading,options.output_file,pbm_dir,pdb_dir,families,options.potential_file, options.radius, fragment_restrict, binding_restrict, options.split_potential,options.auto_mode,options.family_potentials,options.pbm_potentials,options.score_threshold,options.taylor_approach,options.pmf,options.computation,options.known,options.meme,options.reset, options.dummy_dir,options.verbose,options.methylation)
         if options.info is not None:
             log_file=open(options.info,"a")
             log_file.write("%s\tDONE\n"%(os.path.basename(options.output_file)+".meme"))
             log_file.close()
       except Exception as e:
         print("Failed PWM: %s\n"%e)
         if options.info is not None:
             log_file=open(options.info,"a")
             log_file.write("%s\tFAIL\n"%(os.path.basename(options.output_file)+".meme"))
             log_file.close()

    # Execute for all PDB files of a directory
    if os.path.isdir(options.input_files):
     if not os.path.exists(options.output_file): os.makedirs(options.output_file)
     if options.threading:
         pdb_files = [x for x in os.listdir(options.input_files) if x.endswith(".txt")]
     else:
         pdb_files = [x for x in os.listdir(options.input_files) if x.endswith(".pdb")]
     database_file=os.path.join( options.output_file,"database.txt")
     if os.path.exists(database_file): os.remove(database_file)
     name_of_pwms=set()
     for pdb_file in pdb_files:
         if options.threading: output_file = os.path.join( options.output_file,pdb_file.rstrip(".txt") )
         else:                 output_file = os.path.join( options.output_file,pdb_file.rstrip(".pdb") )
         name_of_pwms.add(os.path.basename(output_file)+".meme")
     # Iterate untill all protein profiles are done
     submitted=set()
     n_done=0
     if options.verbose: print("Start iteration to check and run profiles")
     info_file=options.info
     if info_file is None: info_file=os.path.join(options.output_file,"pwm_list_execution.log")
     if options.reuse and functions.fileExist(info_file):
         if options.verbose: print ("Reuse previous information of runs from file %s"%info_file)
         for pdb_file in pdb_files:
           if options.threading:
              input_pdb_file = None
              input_threading_file = os.path.join(options.input_files,pdb_file)
              output_file    = os.path.join( options.output_file,pdb_file.rstrip(".txt") )
           else:
              input_pdb_file = os.path.join(options.input_files,pdb_file)
              input_threading_file = None
              output_file    = os.path.join( options.output_file,pdb_file.rstrip(".pdb") )
           if functions.fileExist(output_file+".meme"):
              if options.verbose:print("\t-- Found PWM %s"%os.path.basename(output_file))
              if options.verbose:print("\t\t-- rewrite the logos...")
              msa_obj=nMSA(output_file+".meme",output_file,"meme")
              write_logo(msa_obj,output_file+".logo")
     else:
         if options.verbose: print ("Open to write %s"%info_file)
         log_file = open(info_file,"w")
         log_file.write("#List of PWMs\n")
         log_file.close()
     done=functions.done_jobs(info_file)
     iterate = functions.check_done(done,name_of_pwms)
     maxtime = 3600 * 2
     start_time = d_time = 0
     while( iterate ):
       for pdb_file in pdb_files:
           if options.threading:
              input_pdb_file = None
              input_threading_file = os.path.join(options.input_files,pdb_file)
              output_file    = os.path.join( options.output_file,pdb_file.rstrip(".txt") )
           else:
              input_pdb_file = os.path.join(options.input_files,pdb_file)
              input_threading_file = None
              output_file    = os.path.join( options.output_file,pdb_file.rstrip(".pdb") )
           if pdb_file in submitted: continue
           submitted.add(pdb_file)
           if options.verbose:sys.stdout.write("Generate PWM %s ...\n"%(os.path.basename(output_file)))
           if functions.fileExist(output_file+".meme"):
              if options.verbose:print("\t-- Found PWM %s"%os.path.basename(output_file))
              if options.reuse:
                 if options.verbose:print("\t\t-- rewrite the logos...")
                 msa_obj=nMSA(output_file+".meme",output_file,"meme")
                 write_logo(msa_obj,output_file+".logo")
              log_file=open(info_file,"r")
              skip_adding=False
              for line in log_file:
                  if os.path.basename(output_file)+".meme" in line.split(): skip_adding=True
              if skip_adding and options.verbose: print("\t\t-- Already in the information file as DONE")
              submitted.add(pdb_file)
              if skip_adding: continue
              if options.verbose: print("\t\t-- Add in the list of done")
              log_file.close()
              log_file=open(info_file,"a")
              log_file.write("%s\tDONE\n"%(os.path.basename(output_file)+".meme"))
              log_file.flush()
              log_file.close()
              continue
           if options.reuse:
              if options.verbose:print("\t-- Not found protein profile %s"%os.path.basename(output_file))
              log_file = open(info_file,"r")
              skip_adding=False
              for line in log_file:
                       if os.path.basename(output_file)+".meme" in line.split(): skip_adding=True
              if skip_adding and options.verbose: print("\t-- Already in the information file but FAILED")
              if skip_adding: submitted.add(pdb_file)
              if skip_adding: continue
              log_file.close()
           if options.parallel:
              if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
              else: cluster_queue=config.get("Cluster", "cluster_queue")
              program=os.path.join(scripts_path,"pwm_pbm.py")
              python = os.path.join(config.get("Paths", "python_path"), "python")
              if options.threading:
                  parameters =   " -i %s "%input_threading_file
                  parameters =   parameters +" --threading "
              else:
                  parameters =   " -i %s "%input_pdb_file
              parameters =  parameters + " -o %s "%output_file
              parameters =  parameters + " --pdb=%s "%options.pdb_dir
              parameters =  parameters + " --pbm=%s "%options.pbm_dir
              parameters =  parameters + " --dummy=%s "%options.dummy_dir
              parameters =  parameters + " -s %s "%options.split_potential
              parameters =  parameters + "--info %s "%info_file
              if options.potential_file  is not None: parameters = parameters + " --file %s "%options.potential_file
              if options.score_threshold is not None: parameters =  parameters + " -t %s "%options.score_threshold
              if options.radius is not None    : parameters = parameters + " --radius %s "%options.radius
              if binding_restrict is not None  : parameters = parameters + " --binding %s "%options.binding_restrict
              if fragment_restrict is  not None: parameters = parameters + " --fragment %s "%options.fragment_restrict
              if options.verbose:                parameters = parameters + " --verbose "
              if options.auto_mode        :      parameters = parameters + " --auto "
              if options.family_potentials:      parameters = parameters + " --family "
              if options.pmf              :      parameters = parameters + " --pmf "
              if options.pbm_potentials   :      parameters = parameters + " -p "
              if options.known            :      parameters = parameters + " --known "
              if options.taylor_approach  :      parameters = parameters + " --taylor "
              if options.computation      :      parameters = parameters + " --bins "
              if options.reset            :      parameters += " --reset "
              if options.methylation      :      parameters = parameters + " --methylation "
              if options.verbose:print("\t-- Submmit %s %s"%(program,parameters))
              functions.submit_command_to_queue("%s %s %s" % (python,program,parameters), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
              submitted.add(pdb_file)
           else:
            try:
              get_single_pwm(input_pdb_file,input_threading_file,options.threading,output_file,pbm_dir,pdb_dir,families,options.potential_file, options.radius,fragment_restrict, binding_restrict, options.split_potential,options.auto_mode,options.family_potentials,options.pbm_potentials,options.score_threshold,options.taylor_approach,options.pmf,options.computation,options.known,options.meme,options.reset, options.dummy_dir,options.verbose,options.methylation)
              log_file=open(info_file,"a")
              log_file.write("%s\tDONE\n"%(os.path.basename(output_file)+".meme"))
              log_file.flush()
              log_file.close()
            except:
              log_file=open(info_file,"a")
              log_file.write("%s\tFAIL\n"%(os.path.basename(output_file)+".meme"))
              log_file.flush()
              log_file.close()
            submitted.add(pdb_file)

       #Check next iteration, profiles submitted and profiles done
       done=functions.done_jobs(info_file)
       iterate= functions.check_done(done,name_of_pwms) 
       if len(done) > n_done:
           n_done=len(done)
           if n_done>1 and start_time==0: start_time = float(time.time())
           if options.verbose: 
                sys.stdout.write("Number of profiles already done %d\n"%n_done)
                if d_time>0: sys.stdout.write("Time: %f\n"%d_time)
                sys.stdout.write("\t-- Check files done ...\n")
                log_file = open(info_file,"r")
                for line in log_file:
                   print("\t\t-- %s"%line.strip())
                log_file.flush()
                log_file.close()
                sys.stdout.write("\t-- Still running protein profiles  %s ...\n"%functions.check_done(done,name_of_pwms))
                sys.stdout.write("\t-- Continue iteration %s ...\n"%iterate)
                sys.stdout.flush() 
       #Check next iteration, if exceeding time and enough profiles stop iteration
       if start_time>0:
          current_time = float(time.time())
          d_time = current_time - start_time
          if float(  len(done) ) / len(name_of_pwms) > complete  and d_time > maxtime: iterate=False
          if d_time > maxtime and iterate and complete < 1.0: 
               complete = complete - 0.01 
               maxtime  = maxtime  + 120
               if options.verbose: sys.stdout.write("Time: %f Done: %d (%d) Ratio to end: %f\n"%(d_time,len(done),len(name_of_pwms),complete))


     #Generate DATABASE
     for name_pwm in name_of_pwms:
        os.system("cat %s >> %s" % (os.path.join(options.output_file,name_pwm), database_file))
    #Done
    if options.verbose:sys.stdout.write("Done\n")
