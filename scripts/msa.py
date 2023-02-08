import os, sys, re
import ConfigParser

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
           self._motif=file_name.split("/")[-1]
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
        if option != "msa" and option != "pwm" and option != "meme":
           option = "msa"
        if option=="msa":
           self.clean_pwm()
           for line in functions.parse_file(self._file):
               if line.startswith("#"): continue
               sequence,score = line.strip().split(";")
               self.add_sequence(sequence, float(score))
           self.set_binding_site_length()
           self.set_pwm()
        if option=="pwm":
           matrix={}
           for line in functions.parse_file(self._file):
               data=line.strip().split("\t")
               residue=data[0].split(":")[0]
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
        x.set_binding_site_length(len(self.get_pwm())+len(other.get_pwm()))
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
        self._motif=name

    def set_pwm(self):
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
                    self.add_sequence(seq,sc-0.001)

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
        return self._motif

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
      self._max_sequences=(10,50,250,500)
      self._min_binding_site_length=(4,6,8)
      super(nMSA,self).__init__(file_name, motif_name, option)
    def get_nucleotides(self):
      return self._nucleotides


class pMSA(MSA):
    residues = list("ACDEFGHIKLMNPQRSTVWY")
    def __init__(self,file_name=None, motif_name=None, option="msa"):
      self._aminoacids = self.residues
      self._short_sequences =set()
      self._max_sequences=(100,500,1000,5000)
      self._min_binding_site_length=(4,8,10)
      self._gapped_sequences =set()
      self._short_pwm = []
      self._gapped_pwm = []
      self._option=option
      if file_name is not None:
         self.parse_file()
      super(pMSA,self).__init__(file_name, motif_name,option)
    def get_aminoacids(self):
        return self._aminoacids
    def parse_file(self):
      if self._option=="msa_short":
         self.clean_pwm()
         for line in functions.parse_file(self._file):
             if line.startswith("#"): continue
             short_sequence, score = line.strip().split(";")
             self.add_short_sequence(short_sequence, float(score))
         self.set_binding_site_length()
         self.set_short_pwm()
      if self._option=="msa_gapped":
         self.clean_pwm()
         for line in functions.parse_file(self._file):
             if line.startswith("#"): continue
             gapped_sequence, score = line.strip().split(";")
             self.add_gapped_sequence(gapped_sequence, float(score))
         self.set_binding_site_length()
         self.set_gapped_pwm()

    def set_gapped_pwm(self):
        if len(self._gapped_sequences)>0:
           if len(self._pwm)>0: self.clean_pwm()
           pfm = []
           binding = zip(*[i[0] for i in self.get_gapped_sequences()])
           self._binding_site_length = len(binding)
           sequences=self.get_gapped_sequences()
           # Count residues instances #
           for i in range(self.get_binding_site_length()):
             vector=[]
             for res in self.residues:
               vector.append(binding[i].count(res)+1)
             pfm.append(vector)
           # For each position... #
           for i in range(self.get_binding_site_length()):
             vector=[]
             for j in range(len(self.residues)):
               vector.append("%.3f" % (float(pfm[i][j]) / sum(pfm[i])))
             self._pwm.append(vector)        

    def set_short_pwm(self):
        if len(self._short_sequences)>0:
           if len(self._pwm)>0: self.clean_pwm()
           pfm = []
           binding = zip(*[i[0] for i in self.get_short_sequences()])
           self._binding_site_length = len(binding)
           sequences=self.get_short_sequences()
           # Count residues instances #
           for i in range(self.get_binding_site_length()):
             vector=[]
             for res in self.residues:
               vector.append(binding[i].count(res)+1)
             pfm.append(vector)
           # For each position... #
           for i in range(self.get_binding_site_length()):
             vector=[]
             for j in range(len(self.residues)):
               vector.append("%.3f" % (float(pfm[i][j]) / sum(pfm[i])))
             self._pwm.append(vector)          


    def add_short_sequence(self, short_sequence, score):
      self._short_sequences.add((short_sequence, score))

    def add_gapped_sequence(self, gapped_sequence, score):
      self._gapped_sequences.add((gapped_sequence, score))

    def write(self,file_name=None, option="msa", overwrite=False):
      if option == "meme_short": 
         self.set_short_pwm()
         super(pMSA,self).write(file_name,"meme",overwrite)
      elif option == "meme_gapped": 
         self.set_gapped_pwm()
         super(pMSA,self).write(file_name,"meme",overwrite)
      else:
         super(pMSA,self).write(file_name,option,overwrite)


