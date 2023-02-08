import sys
import os
import shutil
import subprocess
import numpy as np
import pandas as pd
import argparse
import ConfigParser
import shutil
import optparse
import re
#from scipy.spatial.transform import Rotation
from Bio.PDB import PDBIO, Superimposer, PDBParser,MMCIFParser,MMCIFIO


#"/users/sbi/boliva/ModCRE/scripts"

scripts_path = os.path.abspath(os.path.dirname(__file__))
sys.path.append(scripts_path)

# Import my functions #
from functions import parse_fasta_file


# Imports jbonet's module #
from SBI.structure import PDB
from SBI.structure.chain import ChainOfProtein, ChainOfNucleotide
from SBI.structure.residue import ResidueOfNucleotide
from SBI.structure.atom import AtomOfNucleotide
from SBI.structure.contacts import Complex

#-------------#
# Classes     #
#-------------#

class macro_complex(object):
   """
   This is a class with information of all chains, proteins and DNA forming a complex (i.e. transcription complex)
   It is constructed with PDB files of fragmented DNA and a set of PDB complexes

   Initialize:
   name	= name of the complex (mandatory)
   fragments = DNA fragments of maximum length 250bp (empty list otherwise [])
   dna_files = pdb files with DNA coordiantes for each fragment (empty list otherwise [])
   complex_files = dictionary of PDB coordinates and chain assignement of a complex with DNA for each fragment (empty dictionary otherwise {})

   Additional PDB files can be incorporated to define a protein after the core of the complex is defined. 
   Added files will be defined by the _chain_correspondence_plus dictionary for new and previous proteins of the core
   The size of DNA bead is used as thgeshold to calculate the limit of maximum strength with th erestraints

   Variables

   	_id		is the name given to the complex (initialized with "name")
	
	_fragments	list of tuples for fragments of DNA (i.e 1-250, 226-500, 475- etc.)
	
	_dna_files 	list of PDB files with DNA structures, each corresponds with a fragment	(i.e dna_1-250.pdb is the structure of DNA in fragtament 1-250)

	_complexes	dictionary of { fragment:tuple([PDB,remarks]) } PDB,remarks are two unique files on a DNA fragment to built the complex, 
			The PDB file contains proteins and DNA in the fragment, the remarks show the name of the protein associated with each chain in the PDB
			Fore example:

				REMARK	a = DNA
				REMARK	b = DNA
				REMARK	c = P11474.3cbb_B.96-115:0:75_3cbb_B	(Model name of P11474 modelled with 3cbb_B as template)
				REMARK	d = P11474.4cn3_C.30-43:0:77_4cn3_C
				REMARK	e = P11474.1cit_A.24-39:0:86_1cit_A
				REMARK	f = P11474.4hn6_B.186-201:0:67_4hn6
				REMARK	g = P15976.4hc7_A.129-140:0:104_4hc7
 

	_protein_chains		list of protein chains. Proteins with the same code of UniProt are numbered ( :1, :2, :3, etc.)	to distinguish
	
	_protein_sequences 	list of protein sequences in the same order as _protein_chains

	_dna_chains		list of DNA chains. 
	
	_dna_sequences 	 	list of DNA sequences in the same order as _dna_chains

	_restraints	dictionary of restraints of the form { tuple( (pos_1,pos_2) ): (distance,1.0) }, where pos_1 is also a tuple: (molecule, number), with molecule a protein or DNA chain

	_interfaces	dictionary of interfaces  of the form { molecule: set([numbers])}, with molecule a protein or DNA chain and the set of residue numbers with contacts with other molecules




   """
   def __init__(self,name,fragments=None,complex_files=None,PPI_distance=12):
       self._id=name
       if fragments is not None: self._fragments = fragments
       else:                     self._fragments = []
       if complex_files is not None: self._complexes = complex_files
       else:                         self._complexes = {}
       self._dna_chains=[]
       self._protein_chains=[]
       self._protein_sequences=[]
       self._restraints={}
       self._interfaces={}
       self._chain_correspondence_in_core={}
       self._chain_correspondence_plus={}
       self._protein_reference_site={}
       self._topology={}

       #Internal parameters
           
       # Restraints for protein
       self._protein_dsize     = 50
       self._protein_sfactor   = 2
       self._protein_threshold = 15
       self._protein_maxdist   = 20.0
       self._domain_connection = 2
       self._interface_overlap = 3
       self._minimum_interface = 3
       self._minimum_addition  = 10
       self._PPI_distance      = PPI_distance

       # restrictive condition of equal complex (if True, all restraints must be the same and it increases the number of complexes and of the calculation)
       self._restrict_equal = True


   def initialize(self,fasta_file):
       print("\t--Get protein chains")
       self.set_proteins_in_core(fasta_file)
       self.set_interfaces_and_restraints_in_core()
       self.set_topology_core()

   def copy_protein_chains(self,protein_chains):
       self._protein_chains = protein_chains[:]

   def copy_protein_sequences(self,protein_sequences):
       self._protein_sequences = protein_sequences[:]

   def copy_restraints(self,restraints):
       for k,p in restraints.iteritems():
           self._restraints.setdefault(k,p)

   def copy_interfaces(self,interfaces):
       for k,p in interfaces.iteritems():
           q = set([x for x in p])
           self._interfaces.setdefault(k,q)

   def copy_chain_correspondence_in_core(self,chain_correspondence_in_core):
       for k,p in chain_correspondence_in_core.iteritems():
           self._chain_correspondence_in_core.setdefault(k,p)

   def copy_chain_correspondence_plus(self,chain_correspondence_plus):
       for k,p in chain_correspondence_plus.iteritems():
           self._chain_correspondence_plus.setdefault(k,p)

   def copy_topology(self,topology):
       for k,p in topology.iteritems():
           q=p[:]
           self._topology.setdefault(k,q)

   def copy_fragments(self,fragments):
       self._fragments = fragments[:]
   
   def copy_complexes(self,complexes):
       for k,p in complexes.iteritems():
           self._complexes.setdefault(k,p)

   def copy_protein_reference_site(self,protein_reference_site):
       for k,p in protein_reference_site.iteritems():
           q = set([x for x in p])
           self._protein_reference_site.setdefault(k,q)
  
   def set_protein_dsize(self,protein_dsize):
       self._protein_dsize = protein_dsize

   def set_protein_sfactor(self,protein_sfactor):
       self._protein_sfactor = protein_sfactor

   def set_protein_threshold(self,protein_threshold):
       self._protein_threshold = protein_threshold

   def set_protein_maxdist(self,protein_maxdist):
       self._protein_maxdist =  protein_maxdist

   def set_domain_connection(self,domain_connection):
       self._domain_connection = domain_connection

   def set_interface_overlap(self,interface_overlap):
       self._interface_overlap=interface_overlap

   def set_minimum_interface(self,minimum_interface):
       self._minimum_interface=minimum_interface

   def set_minimum_addition(self,minimum_addition):
       self._minimum_addition=minimum_addition

   def set_restrict_equal(self,restrict_equal):
       self._restrict_equal=restrict_equal

   def set_proteins_in_core(self,fasta_file):
       sequences={}
       for protein,sequence in parse_fasta_file(fasta_file,gz=False,clean=True):
          sequences.setdefault(protein,sequence)
       npr=1
       for fragment in self._complexes.keys():
          pdb_file,remark_file= self._complexes[fragment]
          if not os.path.exists(pdb_file) or not os.path.exists(remark_file): continue
          fo = open(remark_file,"r")
          for line in fo:
              if line.startswith("REMARK"):
                 chain   = line.split("=")[0].split()[-1]
                 protein = ".".join(line.split("=")[1].split()[0].split(":")[0].split(".")[:-2])
                 binding = line.split("=")[1].split()[0].split(":")[0].split(".")[-1]
                 if protein == "" or protein == " ": continue
                 if not sequences.has_key(protein): continue
                 fragment_folder="fragment_"+str(fragment[0])+"-"+str(fragment[1])
                 print("\t\t-- Protein "+protein+" number "+str(npr)+" in chain "+chain+" in "+fragment_folder)
                 self._protein_chains.append(protein+":"+str(npr))
                 self._protein_sequences.append(sequences[protein])
                 self._chain_correspondence_in_core.setdefault(protein+":"+str(npr),tuple((fragment,chain)))
                 self._chain_correspondence_in_core.setdefault(tuple((fragment,chain)),protein+":"+str(npr))
                 npr=npr+1
          print("\tDone fragment "+str(fragment))
          fo.close()

    
   def set_interfaces_and_restraints_in_core(self):
       #this method initializes the set of restraints
       for fragment in self._fragments:
           print("\t--Get restraints PPI and TF-DNA in fragment "+str(fragment))
           pdb_file    = self.get_pdb_by_fragment(fragment)
           print("\t--PDB File "+os.path.basename(pdb_file))
           pdb_obj     = PDB(pdb_file)
           pdb_complex = Complex(pdb_obj,PPI_distance=self._PPI_distance)
           start = fragment[0]
           end   = fragment[1]
           for interface in pdb_complex.PPInterfaces:
               protein_chain_id_1 = interface.protein_chain.chain
               protein_chain_id_2 = interface.protein_interactor.chain
               print("\t\t-- Interface protein-protein of chains "+protein_chain_id_1+" "+protein_chain_id_2)
               protein_num_1      = self.get_protein_by_fragment_and_chain(fragment,protein_chain_id_1)
               protein_num_2      = self.get_protein_by_fragment_and_chain(fragment,protein_chain_id_2)
               if protein_num_1 is None: continue
               if protein_num_2 is None: continue
               self._protein_reference_site.setdefault(protein_num_1,set([x for x in range(start,end,1)]))
               self._protein_reference_site.setdefault(protein_num_2,set([x for x in range(start,end,1)]))
               if self.get_reference_site(protein_num_1) is not None and self.get_reference_site(protein_num_2) is None:
                  self._protein_reference_site.setdefault(protein_num_2,self.get_reference_site(protein_num_1))
               if self.get_reference_site(protein_num_2) is not None and self.get_reference_site(protein_num_1) is None:
                  self._protein_reference_site.setdefault(protein_num_1,self.get_reference_site(protein_num_2))
               print("\t\t-- PPI "+protein_num_1+" with "+protein_num_2)
               for contact in interface.contacts:
                   aminoacid1  = contact.aminoacid1.identifier
                   aminoacid2  = contact.aminoacid2.identifier
                   self._interfaces.setdefault(protein_num_1,set()).add(aminoacid1)
                   self._interfaces.setdefault(protein_num_2,set()).add(aminoacid2)
                   vector     = contact.aminoacid2.ca.coordinates - contact.aminoacid1.ca.coordinates
                   distance   = np.sqrt(vector.dot(vector))
                   pos_prot_1 = tuple( (protein_num_1,aminoacid1) )
                   pos_prot_2 = tuple( (protein_num_2,aminoacid2) )
                   self._restraints.setdefault(tuple( (pos_prot_1,pos_prot_2) ) ,(distance,1.0) )

   def set_protein_inner_restrains_in_core(self,protein,dsize=None,sfactor=None,threshold=None,maxdist=None):
       if dsize is None: dsize=self._protein_dsize
       if sfactor is None: sfactor=self._protein_sfactor
       if threshold is None: threshold=self._protein_threshold
       if maxdist is None: maxdist=self._protein_maxdist
       #Generate forces dictionary
       kforce={}
       for i in range(dsize+1):
         if i<threshold:  
            kforce.setdefault(i,1.0)
         else:
            kforce.setdefault((i),1/float(np.power(sfactor,i-threshold)))
       #Get the PDB file of the core to built constraints
       pdb_file,chain_id = self.get_pdb_and_chain_in_core(protein)
       pdb_obj = PDB(pdb_file)
       #Get coordiantes of CA
       ca_pos=[]
       ca_xyz={}
       chain = pdb_obj.get_chain_by_id(chain_id)
       for n in chain.aminoacids:
           ca_pos.append(n.identifier)
           ca_xyz.setdefault(tuple((protein,n.identifier)),n.ca.coordinates)
       #Get distance restraints
       done=[]
       for i in range(len(ca_pos)):
           a= max(0,i - dsize)
           b= min(i + dsize,len(ca_pos))
           pos_i = tuple((protein,ca_pos[i]))
           for j in range(a,b):
             if j==i:continue
             pos_j = tuple((protein,ca_pos[j]))
             if tuple((pos_i,pos_j)) in done:continue
             if ca_xyz.has_key(pos_i) and ca_xyz.has_key(pos_j):
                diff = ca_xyz[pos_i] - ca_xyz[pos_j]
                dist = np.sqrt(diff.dot(diff))
                k    = kforce[np.abs(j-i)]
                if dist>maxdist: continue
                #print("CHECK DISTANCE ",len(phos_pos[c]),a,b,c,i,j,pos_i,pos_j, phos_xyz[pos_i],phos_xyz[pos_j],diff)
                self._restraints.setdefault(tuple((pos_i,pos_j)),tuple((dist,k)))
                done.append(tuple((pos_i,pos_j)))
                done.append(tuple((pos_j,pos_i)))

   def extend_protein_inner_restrains(self,protein,dsize=None,sfactor=None,threshold=None,maxdist=None):
       if dsize is None: dsize=self._protein_dsize
       if sfactor is None: sfactor=self._protein_sfactor
       if threshold is None: threshold=self._protein_threshold
       if maxdist is None: maxdist=self._protein_maxdist
       #Generate forces dictionary
       kforce={}
       for i in range(dsize+1):
         if i<threshold:  
            kforce.setdefault(i,1.0)
         else:
            kforce.setdefault((i),1/float(np.power(sfactor,i-threshold)))
       #Get the PDB files added on the complex to built constraints
       if self._chain_correspondence_plus.has_key(protein):
           for pdb_chain in self._chain_correspondence_plus[protein]:
               pdb_file,chain_id = pdb_chain
               pdb_obj = PDB(pdb_file)
               #Get coordinates of CA
               ca_pos=[]
               ca_xyz={}
               chain = pdb_obj.get_chain_by_id(chain_id)
               for n in chain.aminoacids:
                   ca_pos.append(n.identifier)
                   ca_xyz.setdefault(tuple((protein,n.identifier)),n.ca.coordinates)
               #Get distance restraints
               done=[]
               for i in range(len(ca_pos)):
                   a= max(0,i - dsize)
                   b= min(i + dsize,len(ca_pos))
                   pos_i = tuple((protein,ca_pos[i]))
                   for j in range(a,b):
                     if j==i:continue
                     pos_j = tuple((protein,ca_pos[j]))
                     if tuple((pos_i,pos_j)) in done:continue
                     if ca_xyz.has_key(pos_i) and ca_xyz.has_key(pos_j):
                        diff = ca_xyz[pos_i] - ca_xyz[pos_j]
                        dist = np.sqrt(diff.dot(diff))
                        k    = kforce[np.abs(j-i)]
                        if dist>maxdist: continue
                        #print("CHECK DISTANCE ",len(phos_pos[c]),a,b,c,i,j,pos_i,pos_j, phos_xyz[pos_i],phos_xyz[pos_j],diff)
                        self._restraints.setdefault(tuple((pos_i,pos_j)),tuple((dist,k)))
                        done.append(tuple((pos_i,pos_j)))
                        done.append(tuple((pos_j,pos_i)))


   def set_all_inner_restraints_in_core(self,dsize=None,sfactor=None,threshold=None,maxdist=None):
       if dsize is None: dsize=self._protein_dsize
       if sfactor is None: sfactor=self._protein_sfactor
       if threshold is None: threshold=self._protein_threshold
       if maxdist is None: maxdist=self._protein_maxdist
       for protein in self._protein_chains:
           self.set_protein_inner_restrains_in_core(protein,dsize,sfactor,threshold,maxdist)

   def extend_all_inner_restraints(self,dsize=None,sfactor=None,threshold=None,maxdist=None):
       if dsize is None: dsize=self._protein_dsize
       if sfactor is None: sfactor=self._protein_sfactor
       if threshold is None: threshold=self._protein_threshold
       if maxdist is None: maxdist=self._protein_maxdist
       for protein in self._protein_chains:
           self.extend_protein_inner_restrains(protein,dsize,sfactor,threshold,maxdist)

   def counter_restraints(self):
       counter={}
       for pair,restr in self._restraints.iteritems():
           pos_1,pos_2 = pair
           mol_1,at_1  = pos_1
           mol_2,at_2  = pos_2
           if counter.has_key((mol_1,mol_2)):
              counter[(mol_1,mol_2)] += 1 
           else:
              counter[(mol_1,mol_2)] = 1
       return counter

   def trim_restraints(self,ratio=0.2,minimum=50,threshold=50,interaction_type="ppi"):
       new_restraints={}
       counter = self.counter_restraints()
       new_counter={}
       for pair,restr in self._restraints.iteritems():
           pos_1,pos_2 = pair
           mol_1,at_1  = pos_1
           mol_2,at_2  = pos_2
           if not counter.has_key((mol_1,mol_2)): continue
           if counter[(mol_1,mol_2)] < max(minimum,threshold): 
              new_restraints.setdefault(pair,restr)
              continue
           if new_counter.has_key((mol_1,mol_2)):
              if new_counter[(mol_1,mol_2)][0] < minimum and  sum(new_counter[(mol_1,mol_2)]) + minimum >= counter[(mol_1,mol_2)]: 
                    new_restraints.setdefault(pair,restr)
                    new_counter[(mol_1,mol_2)][0] += 1
                    continue
           if mol_1 in self._protein_chains and mol_2 in self._protein_chains: itype ="ppi"
           if mol_1 in self._dna_chains and mol_2 in self._dna_chains:         itype ="nni"
           if mol_1 in self._protein_chains and mol_2 in self._dna_chains:     itype ="pni"
           if mol_1 in self._dna_chains and mol_2 in self._protein_chains:     itype ="pni"
           if itype == interaction_type or interaction_type=="all":
              r = np.random.random_sample()
              if r < ratio:
                 new_restraints.setdefault(pair,restr)
                 if new_counter.has_key((mol_1,mol_2)):
                    new_counter[(mol_1,mol_2)][0] += 1
                 else:
                    new_counter[(mol_1,mol_2)] = [1,0]
              else:
                 if new_counter.has_key((mol_1,mol_2)):
                    new_counter[(mol_1,mol_2)][1] += 1
                 else:
                    new_counter[(mol_1,mol_2)] = [0,1]
           else:
              new_restraints.setdefault(pair,restr)
       self._restraints ={}
       for pair,restr in new_restraints.iteritems():
           self._restraints.setdefault(pair,restr)
       

   def set_name(self,name):
       self._id = name

   def get_name(self):
       return self._id

   def get_complexes(self):
       return self._complexes

   def get_fragments(self):
       return self._fragments

   def get_complexes(self):
       return self._complexes

   def get_protein_chains(self):
       return self._protein_chains

   def get_protein_sequences(self):
       return self._protein_sequences

   def get_restraints(self):
       return self._restraints

   def get_interfaces(self):
       return self._interfaces

   def get_interface(self,protein):
       if self._interfaces.has_key(protein):
           return self._interfaces.get(protein)
       else:
           return None

   def get_chain_correspondence_in_core(self):
       return self._chain_correspondence_in_core

   def get_chain_correspondence_plus(self):
       return self._chain_correspondence_plus

   def get_topology(self):
       return self._topology

   def get_protein_reference_site(self):
       return self._protein_reference_site

   def get_reference_site(self,protein):
       if not self._protein_reference_site.has_key(protein): return None
       return self._protein_reference_site[protein]

   def get_protein_dsize(self):
       return self._protein_dsize 

   def get_protein_sfactor(self):
       return self._protein_sfactor 

   def get_protein_threshold(self):
       return self._protein_threshold 

   def get_protein_maxdist(self):
       return self._protein_maxdist

   def get_interface_overlap(self):
       return self._interface_overlap

   def get_minimum_interface(self):
       return self._minimum_interface

   def get_minimum_addition(self):
       return self._minimum_addition

   def get_restrict_equal(self):
       return self._restrict_equal

   def get_pdb_by_fragment(self,fragment):
       if self._complexes.has_key(fragment): return self._complexes[fragment][0]
       else: return None

   def get_remarks_by_fragment(self,fragment):
       if self._complexes.has_key(fragment): return self._complexes[fragment][1]
       else: return None

   def get_pdb_and_chain_in_core(self,protein):
       if protein not in self._protein_chains: return None
       if not self._chain_correspondence_in_core.has_key(protein): return None
       fragment,chain = self._chain_correspondence_in_core[protein]
       pdb_file = self._complexes[fragment][0]
       return (pdb_file,chain)

   def get_pdb_and_chain(self,protein):
       if protein not in self._protein_chains: return None
       pdb_chain=set()
       if self._chain_correspondence_in_core.has_key(protein):
          fragment,chain = self._chain_correspondence_in_core[protein]
          pdb_file = self._complexes[fragment][0]
          pdb_chain.add(tuple((pdb_file,chain)))
       if self._chain_correspondence_plus.has_key(protein):
          pdb_chain.update(self._chain_correspondence_plus.get(protein))
       return pdb_chain

   def get_protein_by_fragment_and_chain(self,fragment,chain):
       if not self._chain_correspondence_in_core.has_key(tuple((fragment,chain))):return None
       return self._chain_correspondence_in_core.get(tuple((fragment,chain)))

   def get_proteins_by_fragment(self,fragment):
       proteins=[]
       for p in self._protein_chains:
           if not self._chain_correspondence_in_core.has_key(p): continue
           fragment_p,chain_p = self._chain_correspondence_in_core.get(p)
           if fragment_p == fragment: proteins.append(p)
       return proteins

   def get_pdbchains_by_fragment(self,fragment):
       chains=[]
       for p in self._protein_chains:
           if not self._chain_correspondence_in_core.has_key(p): continue
           fragment_p,chain_p = self._chain_correspondence_in_core.get(p)
           if fragment_p == fragment: chains.append(chain_p)
       return chains

   def get_protein_by_pdb_and_chain_in_core(self,pdb_file,chain):
       fragment_select = self.get_fragment_by_pdb_in_core(pdb_file)
       if fragment_select is None: return None
       if not self._chain_correspondence_in_core.has_key(tuple((fragment_select,chain))): return None
       return self._chain_correspondence_in_core(tuple((fragment_select,chain)))

   def get_protein_by_pdb_and_chain(self,pdb_file,chain):
       fragment_select = self.get_fragment_by_pdb_in_core(pdb_file)
       if fragment_select is not None:
         if self._chain_correspondence_in_core.has_key(tuple((fragment_select,chain))):
            return self._chain_correspondence_in_core(tuple((fragment_select,chain)))
         else:
            failed=True
       if fragment_select is None or failed:
           if not self._chain_correspondence_plus.has_key(tuple((pdb_file,chain))): return None
           return self._chain_correspondence_plus.get(tuple((pdb_file,chain)))

   def get_fragment_by_pdb_in_core(self,pdb_file):
       fragment_select=None
       for fragment in self._complexes.keys():
           pdb,remark = self._complexes[fragment]
           if pdb == pdb_file or pdb_file in pdb or pdb in pdb_file:
               fragment_select = fragment
               break
       if fragment_select is None: return None
       return fragment_select

   def copy(self,name):
       other = self.__class__(name)
       other.copy_fragments(self.get_fragments())
       other.copy_complexes(self.get_complexes())
       other.copy_protein_chains(self.get_protein_chains())
       other.copy_protein_sequences(self.get_protein_sequences())
       other.copy_restraints(self.get_restraints())
       other.copy_interfaces(self.get_interfaces())
       other.copy_chain_correspondence_in_core(self.get_chain_correspondence_in_core())
       other.copy_chain_correspondence_plus(self.get_chain_correspondence_plus())
       other.copy_topology(self.get_topology())
       other.copy_protein_reference_site(self.get_protein_reference_site())
       other.set_protein_dsize(      self.get_protein_dsize())
       other.set_protein_sfactor(    self.get_protein_sfactor())
       other.set_protein_threshold(  self.get_protein_threshold())
       other.set_protein_maxdist(    self.get_protein_maxdist())
       other.set_interface_overlap(  self.get_interface_overlap())
       other.set_minimum_interface(  self.get_minimum_interface())
       other.set_minimum_addition(   self.get_minimum_addition())
       other.set_restrict_equal(     self.get_restrict_equal())
       return other

   def add_chain_correspondence_plus(self,protein,pdb,chain):
       self._chain_correspondence_plus.setdefault(protein,set()).add(tuple([pdb,chain]))
       self._chain_correspondence_plus.setdefault(tuple([pdb,chain]),protein)

   def accept_addition(self,pdb_obj,protein_num,protein_chain,partner,partner_chain,restraints,condition):
       # Add the new protein
       add_ppi_restraints=True
       print("\t\t--check topology of the potential new interaction "+protein_num+" "+partner)
       ppi_protein_1 = protein_num.split(":")[0]
       ppi_protein_2 = partner.split(":")[0]
       #Check the new protein residues (restraints and interface)
       chain1      = pdb_obj.get_chain_by_id(protein_chain)
       start       = int(chain1.first_aminoacid.identifier)
       end         = int(chain1.last_aminoacid.identifier)
       residues=set()
       for i in range(len(self._topology["molecule_name"])):
           if self._topology["molecule_name"][i] == protein_num:
            try:
              residue_range = self._topology["residue_range"][i]
              rr_start,rr_end = residue_range.split(",")
              for aa in range(int(rr_start),int(rr_end)+1):
                  residues.add(aa)
            except Exception as e:
              print("Skip protein "+str(protein_num))
              print("ERROR: %s"%str(e))
              continue
       all_residues = [x for x in residues]
       all_residues.sort()
       original_range="0-0"
       if len(all_residues)>0: 
          original_range=str(all_residues[0])+" - "+str(all_residues[-1])
       original_gaps=[]
       for i in range(len(all_residues)):
           if i+1 <len(all_residues):
              if all_residues[i+1] - all_residues[i] > 1:
                 ostart = all_residues[i]
                 oend   = all_residues[i+1]
                 original_gaps.append(tuple((ostart,oend)))
       new_residues = set([x for x in range(int(start),int(end)+1)])
       new_residues = [x for x in new_residues.difference(residues)]
       #Check that restraints are from new residues
       contact_residues=set()
       for pair in restraints.keys():
           pos_prot_1,pos_prot_2 = pair
           p_1,aminoacid1  = pos_prot_1
           p_2,aminoacid2  = pos_prot_2
           if p_1 == ppi_protein_1 :
              contact_residues.add(int(aminoacid1))
           if p_2 == ppi_protein_1 :
              contact_residues.add(int(aminoacid2))
       if not contact_residues.issubset(new_residues) and condition == "old":
          print("\t\t\t--no new contacts available for  "+protein_num)
          add_ppi_restraints=False
       if len(new_residues) < self._minimum_addition:
          print("\t\t\t--too few residues to consider adding data for "+protein_num)
          add_ppi_restraints=False

       if not add_ppi_restraints:
          return add_ppi_restraints

       #Get new ranges and gaps
       all_residues = [x for x in residues.union(new_residues)]
       all_residues.sort()
       if add_ppi_restraints:
          print("\t\t\t\t--ADDING Range "+protein_num+" changes from "+original_range+" to "+str(all_residues[0])+" - "+str(all_residues[-1]))
       new_residues.sort()
       new_ranges=[]
       start=0
       for i in range(len(new_residues)):
           if start==0: start=new_residues[i]
           if i+1 <len(new_residues):
              if new_residues[i+1] - new_residues[i] > 1:
                 end   = new_residues[i]
                 new_ranges.append(tuple((start,end)))
                 start = new_residues[i+1]
           else:
                 end   = new_residues[i]
                 new_ranges.append(tuple((start,end)))
       gaps=[]
       for i in range(len(all_residues)):
           if i+1 <len(all_residues):
              if all_residues[i+1] - all_residues[i] > 1:
                 start = all_residues[i]
                 end   = all_residues[i+1]
                 gaps.append(tuple((start,end)))
       print("\t\t\t\t--ADDING Ranges "+str(new_ranges))

       # Check the modification of the topology of the partner
       chain_partner      = pdb_obj.get_chain_by_id(partner_chain)
       partner_start      = int(chain_partner.first_aminoacid.identifier)
       partner_end        = int(chain_partner.last_aminoacid.identifier)
       residues=set()
       for i in range(len(self._topology["molecule_name"])):
           if self._topology["molecule_name"][i] == partner:
              try:
                residue_range = self._topology["residue_range"][i]
                rr_start,rr_end = residue_range.split(",")
                for aa in range(int(rr_start),int(rr_end)+1):
                  residues.add(aa)
              except Exception as e:
                print("Skip protein "+str(partner))
                print("ERROR: %s"%str(e))
                continue
       all_residues = [x for x in residues]
       all_residues.sort()
       original_range="0-0"
       if len(all_residues)>0: 
          original_range=str(all_residues[0])+" - "+str(all_residues[-1])
       original_gaps=[]
       for i in range(len(all_residues)):
           if i+1 <len(all_residues):
              if all_residues[i+1] - all_residues[i] > 1:
                 ostart = all_residues[i]
                 oend   = all_residues[i+1]
                 original_gaps.append(tuple((ostart,oend)))
       new_residues = set([x for x in range(partner_start,partner_end+1)])
       new_residues = [x for x in new_residues.difference(residues)]
       new_residues.sort()
       #Check that restraints are from new residues
       contact_residues=set()
       for pair in restraints.keys():
           pos_prot_1,pos_prot_2 = pair
           p_1,aminoacid1  = pos_prot_1
           p_2,aminoacid2  = pos_prot_2
           if p_2 == ppi_protein_2 :
              contact_residues.add(int(aminoacid2))
           if p_1 == ppi_protein_2 :
              contact_residues.add(int(aminoacid1))
       if not contact_residues.issubset(set(new_residues)) and condition=="old": 
          print("\t\t\t--no new contacts available with the partner in the complex "+partner)
          add_ppi_restraints=False
       if len(new_residues) < self._minimum_addition:
          print("\t\t\t--too few residues to consider adding another stretch of "+partner)
          add_ppi_restraints=False

       if not add_ppi_restraints:
          return add_ppi_restraints


       #Get new ranges and gaps
       all_residues = [x for x in residues.union(new_residues)]
       all_residues.sort()
       if add_ppi_restraints:
          print("\t\t\t\t--ADDING Range "+partner+" changes from "+original_range+" to "+str(all_residues[0])+" - "+str(all_residues[-1]))
       new_residues.sort()
       new_ranges=[]
       start=0
       for i in range(len(new_residues)):
           if start==0: start=new_residues[i]
           if i+1 <len(new_residues):
              if new_residues[i+1] - new_residues[i] > 1:
                 end   = new_residues[i]
                 new_ranges.append(tuple((start,end)))
                 start = new_residues[i+1]
           else:
                 end   = new_residues[i]
                 new_ranges.append(tuple((start,end)))
       gaps=[]
       for i in range(len(all_residues)):
           if i+1 <len(all_residues):
              if all_residues[i+1] - all_residues[i] > 1:
                 start = all_residues[i]
                 end   = all_residues[i+1]
                 gaps.append(tuple((start,end)))
       print("\t\t\t\t--ADDING Ranges "+str(new_ranges))


       return add_ppi_restraints


   def add_topology_and_connect(self,number_of_addition,protein_num,pdb_file,protein_chain,partner,partner_chain,restraints,condition=None):
       # Add the new protein
       add_ppi_restraints=True
       print("\t\t\t--check topology for interaction "+protein_num+" "+partner)
       ppi_protein_1 = protein_num.split(":")[0]
       ppi_protein_2 = partner.split(":")[0]
       fasta_name = self._id+".topology.fasta"
       pdb_obj    = PDB(pdb_file)
       chain_id   = protein_chain
       chain      = pdb_obj.get_chain_by_id(chain_id)
       start      = int(chain.first_aminoacid.identifier)
       end        = int(chain.last_aminoacid.identifier)
       offset     = 0
       residues=set()
       for i in range(len(self._topology["molecule_name"])):
           if self._topology["molecule_name"][i] == protein_num:
              try:
                residue_range = self._topology["residue_range"][i]
                rr_start,rr_end = residue_range.split(",")
                for aa in range(int(rr_start),int(rr_end)+1):
                  residues.add(aa)
              except Exception as e:
                print("Skip protein "+str(protein_num))
                print("ERROR: %s"%str(e))
                continue
       all_residues = [x for x in residues]
       all_residues.sort()
       original_range="0-0"
       if len(all_residues)>0: 
          original_range=str(all_residues[0])+" - "+str(all_residues[-1])
       original_gaps=[]
       for i in range(len(all_residues)):
           if i+1 <len(all_residues):
              if all_residues[i+1] - all_residues[i] > 1:
                 ostart = all_residues[i]
                 oend   = all_residues[i+1]
                 original_gaps.append(tuple((ostart,oend)))
       new_residues = set([x for x in range(int(start),int(end)+1)])
       new_residues = [x for x in new_residues.difference(residues)]
       #Check that restraints are from new residues
       contact_residues=set()
       for pair in restraints.keys():
           pos_prot_1,pos_prot_2 = pair
           p_1,aminoacid1  = pos_prot_1
           p_2,aminoacid2  = pos_prot_2
           if p_1 == ppi_protein_1 :
              contact_residues.add(int(aminoacid1))
           if p_2 == ppi_protein_1 :
              contact_residues.add(int(aminoacid2))
       if not contact_residues.issubset(new_residues) and condition == "old":
          print("\t\t\t\t--DOUBLE CHECK FAILED: no new contacts available for "+protein_num+" for addition "+str(number_of_addition))
          print("\t\t\t\t-- new residues "+str(new_residues))
          print("\t\t\t\t-- contact residues "+str(contact_residues))
          add_ppi_restraints=False
       if len(new_residues) < self._minimum_addition:
          print("\t\t\t\t--DOUBLE CHECK FAILED: too few residues to consider adding "+protein_num+" for addition "+str(number_of_addition))
          print("\t\t\t\t-- new residues "+str(new_residues))
          add_ppi_restraints=False
       #Get new ranges and gaps
       all_residues = [x for x in residues.union(new_residues)]
       all_residues.sort()
       if add_ppi_restraints:
          print("\t\t\t\t--Range "+protein_num+" changes from "+original_range+" to "+str(all_residues[0])+" - "+str(all_residues[-1]))
       new_residues.sort()
       new_ranges=[]
       start=0
       for i in range(len(new_residues)):
           if start==0: start=new_residues[i]
           if i+1 <len(new_residues):
              if new_residues[i+1] - new_residues[i] > 1:
                 end   = new_residues[i]
                 new_ranges.append(tuple((start,end)))
                 start = new_residues[i+1]
           else:
                 end   = new_residues[i]
                 new_ranges.append(tuple((start,end)))
       gaps=[]
       for i in range(len(all_residues)):
           if i+1 <len(all_residues):
              if all_residues[i+1] - all_residues[i] > 1:
                 start = all_residues[i]
                 end   = all_residues[i+1]
                 gaps.append(tuple((start,end)))
       #print("\t\t\t\t--NEW Ranges "+str(new_ranges))
       for new_range in new_ranges:
           start,end = new_range
           if int(end)-int(start) < self._minimum_addition: continue
           self._topology.setdefault("molecule_name",[]).append(protein_num)
           self._topology.setdefault("color",[]).append("green")
           self._topology.setdefault("fasta_fn",[]).append(os.path.basename(fasta_name))
           self._topology.setdefault("fasta_id",[]).append(protein_num)
           self._topology.setdefault("pdb_fn",[]).append(os.path.basename(pdb_file))
           self._topology.setdefault("chain",[]).append(protein_chain)
           self._topology.setdefault("residue_range",[]).append(str(start)+","+str(end))
           self._topology.setdefault("pdb_offset",[]).append(offset)
           self._topology.setdefault("bead_size",[]).append(1)
           self._topology.setdefault("em_residues_per_gaussian",[]).append(0)
           if self._topology.has_key("rigid_body"): rb=len(self._topology.get("rigid_body")) + 1
           else: rb=1
           self._topology.setdefault("rigid_body",[]).append(rb)
           self._topology.setdefault("super_rigid_body",[]).append(" ")
           self._topology.setdefault("chain_of_super_rigid_bodies",[]).append(" ")
       # Add distance restraint between domains of the new protein 
       maximum_gap=0
       for gap in gaps:
           aminoacid1,aminoacid2 = gap
           distance  = min(15.0, float(int(aminoacid2)-int(aminoacid1) ))
           if distance>maximum_gap and gap not in original_gaps: maximum_gap=distance
           restraint = tuple((distance,10.0))
           pos_prot_1 = tuple((protein_num,aminoacid1))
           pos_prot_2 = tuple((protein_num,aminoacid2))
           self._restraints.setdefault(tuple((pos_prot_1,pos_prot_2)),restraint)

       #Check that new PPI restraints are possible
       if self._chain_correspondence_in_core.has_key(protein_num) and condition=="old":
           if maximum_gap < self._domain_connection:
              print("\t\t\t\t--The new domain for "+protein_num+" has not enough flexibility to be added") 
              add_ppi_restraints=False


       # Complete the topology of the partner
       chain_id           = partner_chain
       chain_partner      = pdb_obj.get_chain_by_id(chain_id)
       partner_start      = int(chain_partner.first_aminoacid.identifier)
       partner_end        = int(chain_partner.last_aminoacid.identifier)
       residues=set()
       for i in range(len(self._topology["molecule_name"])):
           if self._topology["molecule_name"][i] == partner:
              try:
                residue_range = self._topology["residue_range"][i]
                rr_start,rr_end = residue_range.split(",")
                for aa in range(int(rr_start),int(rr_end)+1):
                  residues.add(aa)
              except Exception as e:
                print("Skip protein "+str(partner))
                print("ERROR: %s"%str(e))
                continue
       all_residues = [x for x in residues]
       all_residues.sort()
       original_range="0-0"
       if len(all_residues)>0: 
          original_range=str(all_residues[0])+" - "+str(all_residues[-1])
       original_gaps=[]
       for i in range(len(all_residues)):
           if i+1 <len(all_residues):
              if all_residues[i+1] - all_residues[i] > 1:
                 ostart = all_residues[i]
                 oend   = all_residues[i+1]
                 original_gaps.append(tuple((ostart,oend)))
       new_residues = set([x for x in range(partner_start,partner_end+1)])
       new_residues = [x for x in new_residues.difference(residues)]
       new_residues.sort()
       #Check that restraints are from new residues
       contact_residues=set()
       for pair in restraints.keys():
           pos_prot_1,pos_prot_2 = pair
           p_2,aminoacid2  = pos_prot_2
           p_1,aminoacid1  = pos_prot_1
           if p_2 == ppi_protein_2 :
              contact_residues.add(int(aminoacid2))
           if p_1 == ppi_protein_2 :
              contact_residues.add(int(aminoacid1))
       if not contact_residues.issubset(set(new_residues)) and condition=="old": 
          print("\t\t\t\t--DOUBLE CHECK FAILED: no new contacts available for "+partner+" for addition "+str(number_of_addition))
          print("\t\t\t\t-- new residues "+str(new_residues))
          print("\t\t\t\t-- contact residues "+str(contact_residues))
          add_ppi_restraints=False
       if len(new_residues) < self._minimum_addition:
          print("\t\t\t\t--DOUBLE CHECK FAILED: too few residues to consider adding "+partner+" for addition "+str(number_of_addition))
          print("\t\t\t\t-- new residues "+str(new_residues))
          add_ppi_restraints=False
       #Get new ranges and gaps
       all_residues = [x for x in residues.union(new_residues)]
       all_residues.sort()
       if add_ppi_restraints:
          print("\t\t\t\t--Range "+partner+" changes from "+original_range+" to "+str(all_residues[0])+" - "+str(all_residues[-1]))
       new_residues.sort()
       new_ranges=[]
       start=0
       for i in range(len(new_residues)):
           if start==0: start=new_residues[i]
           if i+1 <len(new_residues):
              if new_residues[i+1] - new_residues[i] > 1:
                 end   = new_residues[i]
                 new_ranges.append(tuple((start,end)))
                 start = new_residues[i+1]
           else:
                 end   = new_residues[i]
                 new_ranges.append(tuple((start,end)))
       gaps=[]
       for i in range(len(all_residues)):
           if i+1 <len(all_residues):
              if all_residues[i+1] - all_residues[i] > 1:
                 start = all_residues[i]
                 end   = all_residues[i+1]
                 gaps.append(tuple((start,end)))
       for new_range in new_ranges:
           start,end = new_range
           if int(end)-int(start) < self._minimum_addition: continue
           self._topology.setdefault("molecule_name",[]).append(partner)
           if self.get_fragment_by_pdb_in_core(os.path.basename(pdb_file)) is None: color_partner="green"
           else:  color_partner="blue"
           self._topology.setdefault("color",[]).append(color_partner)
           self._topology.setdefault("fasta_fn",[]).append(os.path.basename(fasta_name))
           self._topology.setdefault("fasta_id",[]).append(partner)
           self._topology.setdefault("pdb_fn",[]).append(os.path.basename(pdb_file))
           self._topology.setdefault("chain",[]).append(partner_chain)
           self._topology.setdefault("residue_range",[]).append(str(start)+","+str(end))
           self._topology.setdefault("pdb_offset",[]).append(offset)
           self._topology.setdefault("bead_size",[]).append(1)
           self._topology.setdefault("em_residues_per_gaussian",[]).append(0)
           if self._topology.has_key("rigid_body"): rb=len(self._topology.get("rigid_body")) + 1
           else: rb=1
           self._topology.setdefault("rigid_body",[]).append(rb)
           self._topology.setdefault("super_rigid_body",[]).append(" ")
           self._topology.setdefault("chain_of_super_rigid_bodies",[]).append(" ")

       # Add distance restraint between domains of the partner
       maximum_gap=0
       for gap in gaps:
           aminoacid1,aminoacid2 = gap
           distance  = min(15.0, float(int(aminoacid2)-int(aminoacid1) ))
           if distance>maximum_gap and gap not in original_gaps: maximum_gap=distance
           restraint = tuple((distance,10.0))
           pos_prot_1 = tuple((partner,aminoacid1))
           pos_prot_2 = tuple((partner,aminoacid2))
           self._restraints.setdefault(tuple((pos_prot_1,pos_prot_2)),restraint)
       #Check that new PPI restraints are possible
       if self._chain_correspondence_in_core.has_key(partner) and condition=="old" and add_ppi_restraints:
           if maximum_gap < self._domain_connection : 
              print("\t\t\t\t--The new domain for "+partner+" has not enough flexibility to be added BUT we still keep it for the new protein") 
              #add_ppi_restraints=False
       if add_ppi_restraints: print("\t\t\t\t--ADD RESTRAINT ")
       else:print("\t\t\t\t--SKIP RESTRAINT (no new contacts)")

       return add_ppi_restraints


   def add_restraints(self,new_restraints,new_protein,partner_in_complex):
       protein_1 = new_protein.split(":")[0]
       protein_2 = partner_in_complex.split(":")[0]
       for pair,restr in new_restraints.iteritems():
           pos_prot_1,pos_prot_2 = pair
           p_1,aminoacid1  = pos_prot_1
           p_2,aminoacid2  = pos_prot_2
           if p_1 == protein_1 and p_2 == protein_2:
              new_pos_prot_1 = tuple((new_protein,aminoacid1))
              new_pos_prot_2 = tuple((partner_in_complex,aminoacid2))
              self._restraints.setdefault(tuple((new_pos_prot_1,new_pos_prot_2)),restr)
           else:
              print("ERROR CHECK PROGRAM ON restraints "+p_1+" != "+protein_1+" or "+p_2+" != "+protein_2)
           
            
   def add_interfaces(self,interface_1,interface_2,new_protein,partner_in_complex):
       protein_1 = new_protein.split(":")[0]
       protein_2 = partner_in_complex.split(":")[0]
       for p,aas in interface_1.iteritems():
           if p==protein_1:
              self._interfaces.setdefault(new_protein,set()).update(aas)
           else:
              print("ERROR CHECK PROGRAM ON  interface_1 "+p+" != "+protein_1)
       for p,aas in interface_2.iteritems():
           if p==protein_2:
              self._interfaces.setdefault(partner_in_complex,set()).update(aas)
           else:
              print("ERROR CHECK PROGRAM ON  interface_2 "+p+" != "+protein_2)
    

   def add_protein(self,protein,sequence,partner_with_reference):
       number = len(self._protein_chains) + 1
       protein_num=protein+":"+str(number)
       self._protein_chains.append(protein_num)
       self._protein_sequences.append(sequence)
       if self.get_reference_site(partner_with_reference) is None:
          print("ERROR: Missing reference site for "+partner_with_reference)
       self.set_reference_site(protein_num,self.get_reference_site(partner_with_reference))
       return protein_num

   def add_ppi(self,pdb_file,protein_1,protein_2,sequence_1,sequence_2,condition=None,maximum=None):
       # Additional interactions will be colored green
       print("\t   check addition of PPI "+os.path.basename(pdb_file))
       protein_codes = set([x.split(":")[0] for x in self._protein_chains])
       if protein_1 not in protein_codes and protein_2 not in protein_codes: return (list(),0)
       try:
         pdb_obj = PDB(pdb_file)
       except:
         print("\t\t--skip unreadable file "+pdb_file)
         return (list(),0)
       pdb_name=os.path.basename(pdb_file)
       chain_ids=[]
       protein_chain={}
       for chain in pdb_obj.chains:
           if chain.chaintype == "P":
              chain_ids.append(chain.chain)
       chain_ids.sort()
       protein_chain.setdefault(protein_1+":A",chain_ids[0])
       protein_chain.setdefault(protein_2+":B",chain_ids[1])
       protein_chain.setdefault(chain_ids[0],protein_1)
       protein_chain.setdefault(chain_ids[1],protein_2)
       set_1=set_2=[]
       if protein_1 in protein_codes:
          set_1=[protein for protein in self.get_protein_chains() if protein.split(":")[0] == protein_1]
       if protein_2 in protein_codes:
          set_2=[protein for protein in self.get_protein_chains() if protein.split(":")[0] == protein_2]
       pdb_complex = Complex(pdb_obj,PPI_distance=self._PPI_distance)
       interface_1={}
       interface_2={}
       restraints_1={}
       restraints_2={}
       for interface in pdb_complex.PPInterfaces:
               protein_chain_id_1 = interface.protein_chain.chain
               protein_chain_id_2 = interface.protein_interactor.chain
               #print("\t\t-- Interface "+ os.path.basename(pdb_file)+" protein-protein of chains "+protein_chain_id_1+" "+protein_chain_id_2)
               protein_num_1 = None
               protein_num_2 = None
               if protein_chain.has_key(protein_chain_id_1): protein_num_1      = protein_chain[protein_chain_id_1]
               if protein_chain.has_key(protein_chain_id_2): protein_num_2      = protein_chain[protein_chain_id_2]
               #print("\t\t-- PPI "+protein_num_1+" with "+protein_num_2)
               if protein_num_1 is None: return (list(),0)
               if protein_num_2 is None: return (list(),0)
               if len([x for x in interface.contacts])<=0: return (list(),0)
               for contact in interface.contacts:
                   aminoacid1  = contact.aminoacid1.identifier
                   aminoacid2  = contact.aminoacid2.identifier
                   vector     = contact.aminoacid2.ca.coordinates - contact.aminoacid1.ca.coordinates
                   distance   = np.sqrt(vector.dot(vector))
                   if protein_num_1 == protein_1 and protein_num_2 == protein_2:
                      pos_prot_1 = tuple( (protein_1,aminoacid1) )
                      pos_prot_2 = tuple( (protein_2,aminoacid2) )
                      restraints_1.setdefault(tuple( (pos_prot_1,pos_prot_2) ) ,(distance,1.0) )
                      restraints_2.setdefault(tuple( (pos_prot_2,pos_prot_1) ) ,(distance,1.0) )
                      interface_1.setdefault(protein_1,set()).add(aminoacid1)
                      interface_2.setdefault(protein_2,set()).add(aminoacid2)
                   elif protein_num_1 == protein_2 and protein_num_2 == protein_1:
                      pos_prot_1 = tuple( (protein_1,aminoacid2) )
                      pos_prot_2 = tuple( (protein_2,aminoacid1) )
                      restraints_1.setdefault(tuple( (pos_prot_1,pos_prot_2) ) ,(distance,1.0) )
                      restraints_2.setdefault(tuple( (pos_prot_2,pos_prot_1) ) ,(distance,1.0) )
                      interface_1.setdefault(protein_1,set()).add(aminoacid2)
                      interface_2.setdefault(protein_2,set()).add(aminoacid1)

       adding=set()
       for protein in set_1:
           interface = self.get_interface(protein)
           if interface is None: continue
           #print("\t\t-- Interface overlap %d"%(len(interface.intersection(interface_1))))
           if len(interface.intersection(interface_1[protein_1]))>self._interface_overlap:
              print("\t\t--Interface for PPI "+protein+"-"+protein_2+" is already occupied")
              continue
           if condition == "old" and protein_2 not in protein_codes: 
              #print("\t\t-- Condition "+condition+" protein "+protein_2+" NOT IN "+str(protein_codes))
              continue
           if condition == "new" and protein_2 in protein_codes: 
              #print("\t\t-- Condition "+condition+" protein "+protein_2+" IN "+str(protein_codes))
              continue
           if condition == "old":
             right_oriented = True
             for p in set_2:
                 if p==protein:continue
                 if tuple([protein,p,right_oriented]) in adding: continue
                 if len(interface_1[protein_1]) < self._minimum_interface: continue
                 if len(interface_2[protein_2]) < self._minimum_interface: continue
                 if not protein_chain.has_key(protein_1+":A"):continue
                 if not protein_chain.has_key(protein_2+":B"):continue
                 if not self.accept_addition(pdb_obj,protein,protein_chain[protein_1+":A"],p,protein_chain[protein_2+":B"],restraints_1,condition): continue
                 try:
                    print("\t\t--Add interaction "+protein+" with "+p)
                    original_interface ="%s-%s"%(min([int(x) for x in interface]),max([int(x) for x in interface]))
                    added_interface    ="%s-%s"%(min([int(x) for x in interface_1[protein_1]] ),max([int(x) for x in interface_1[protein_1]]))
                    print("\t\t--changes of interface in "+protein+" : from "+original_interface+" add "+added_interface)
                 except:
                    print("\t\t--changes of interface in "+protein+" cannot be read, some amino-acids are defined with letter positions")
                 adding.add(tuple([protein,p,right_oriented]))
                 #print("OLD P-set1 P-set2 "+str(adding))
           else:
             right_oriented=False
             if tuple([protein_2,protein,right_oriented]) in adding: continue
             if len(interface_1[protein_1]) < self._minimum_interface: continue
             if len(interface_2[protein_2]) < self._minimum_interface: continue
             if not protein_chain.has_key(protein_1+":A"):continue
             if not protein_chain.has_key(protein_2+":B"):continue
             if not self.accept_addition(pdb_obj,protein_2,protein_chain[protein_2+":B"],protein,protein_chain[protein_1+":A"],restraints_2,condition): continue
             print("\t\t--Add interaction "+protein+" with "+protein_2)
             try:
                original_interface ="%s-%s"%(min([int(x) for x in interface]),max([int(x) for x in interface]))
                added_interface    ="%s-%s"%(min([int(x) for x in interface_1[protein_1]] ),max([int(x) for x in interface_1[protein_1]]))
                print("\t\t--changes of interface in "+protein+" : from "+original_interface+" add "+added_interface)
             except:
                print("\t\t--changes of interface in "+protein+" cannot be read, some amino-acids are defined with letter positions")
             adding.add(tuple([protein_2,protein,right_oriented]))
             #print("new or none P-set1 prot_2 "+str(adding))
       for protein in set_2:
           interface = self.get_interface(protein)
           if interface is None: continue
           #print("\t\t-- Interface overlap %d"%(len(interface.intersection(interface_2))))
           if len(interface.intersection(interface_2[protein_2]))>self._interface_overlap: 
              print("\t\t--Interface for PPI "+protein+"-"+protein_1+" is already occupied")
              continue
           if condition == "old" and protein_1 not in protein_codes: 
              #print("\t\t-- Condition "+condition+" protein "+protein_1+" NOT IN "+str(protein_codes))
              continue
           if condition == "new" and protein_1 in protein_codes: 
              #print("\t\t-- Condition "+condition+" protein "+protein_1+" IN "+str(protein_codes))
              continue
           if condition == "old":
             for p in set_1:
                 right_oriented=False
                 if p==protein:continue
                 if tuple([protein,p]) in adding: continue
                 if len(interface_1[protein_1]) < self._minimum_interface: continue
                 if len(interface_2[protein_2]) < self._minimum_interface: continue
                 if not protein_chain.has_key(protein_1+":A"):continue
                 if not protein_chain.has_key(protein_2+":B"):continue
                 if not self.accept_addition(pdb_obj,protein,protein_chain[protein_2+":B"],p,protein_chain[protein_1+":A"],restraints_2,condition): continue
                 try:
                    print("\t\t--Add interaction "+protein+" with "+p)
                    original_interface ="%s-%s"%(min([int(x) for x in interface]),max([int(x) for x in interface]))
                    added_interface    ="%s-%s"%(min([int(x) for x in interface_2[protein_2]] ),max([int(x) for x in interface_2[protein_2]]))
                    print("\t\t--changes of interface in "+protein+" : from "+original_interface+" add "+added_interface)
                 except:
                    print("\t\t--changes of interface in "+protein+" cannot be read, some amino-acids are defined with letter positions")
                 adding.add(tuple([protein,p,right_oriented]))
                 #print("OLD p-set2 P-set1 "+str(adding))
           else:
             right_oriented=True
             if tuple([protein_1,protein,right_oriented]) in adding: continue
             if len(interface_1[protein_1]) < self._minimum_interface: continue
             if len(interface_2[protein_2]) < self._minimum_interface: continue
             if not protein_chain.has_key(protein_1+":A"):continue
             if not protein_chain.has_key(protein_2+":B"):continue
             if not self.accept_addition(pdb_obj,protein_1,protein_chain[protein_1+":A"],protein,protein_chain[protein_2+":B"],restraints_1,condition): continue
             try:
                print("\t\t--Add interaction "+protein+" with "+protein_1)
                original_interface ="%s-%s"%(min([int(x) for x in interface]),max([int(x) for x in interface]))
                added_interface    ="%s-%s"%(min([int(x) for x in interface_2[protein_2]] ),max([int(x) for x in interface_2[protein_2]]))
                print("\t\t--changes of interface in "+protein+" : from "+original_interface+" add "+added_interface)
             except:
                print("\t\t--changes of interface in "+protein+" cannot be read, some amino-acids are defined with letter positions")
             adding.add(tuple([protein_1,protein,right_oriented]))
             #print("new or none P-set2 prot1 "+str(adding))

      
       if len(adding) == 0: 
          print("\t\t--Nothing to add")
          return (list(),0)


       print("\t\t--Potential protein interactions to add: %d (between %s-%s )"%(len(adding),protein_1,protein_2))
       if maximum is not None:  new_adding = [x for x in adding][:int(maximum)]
       else:                    new_adding = [x for x in adding]
       adding = new_adding
       for i in range(len(adding)):
           print("\t\t\t-- number of addition %d add  %s  "%(i,str(adding[i]) ) )
       new_restraints = 0
       seed_name = self.get_name() + ".0"
       #preserve the original to add other combinations of PPIs
       seed_complex = self.copy(seed_name)

       new_protein,partner_in_complex,right_oriented  = adding[0]
       print("\t\t  --check adding %s  - %s to %s"%(new_protein,partner_in_complex,self.get_name()))
       check = self.get_interface(partner_in_complex)
       if condition == "old":
         if right_oriented:
          #new_protein.split(":")[0] == protein_1:
          protein_num   = new_protein
          chain_id      = protein_chain[protein_1+":A"]
          partner_chain = protein_chain[protein_2+":B"]
          add_ppi_restraints=self.add_topology_and_connect(0,protein_num,pdb_file,chain_id,partner_in_complex,partner_chain,restraints_1,condition)
          if add_ppi_restraints:
             print("\t\t\t-- Add restrain contacts (%s %s) %d"%(protein_num,partner_in_complex,len([x for x in restraints_1.keys()])))
             new_restraints  = new_restraints +len([x for x in restraints_1.keys()])
             self.add_restraints(restraints_1,protein_num,partner_in_complex)
             self.add_interfaces(interface_1,interface_2,protein_num,partner_in_complex)
         else:
          #new_protein.split(":")[0] == protein_2: 
          protein_num   = new_protein
          chain_id      = protein_chain[protein_2+":B"]
          partner_chain = protein_chain[protein_1+":A"]
          add_ppi_restraints=self.add_topology_and_connect(0,protein_num,pdb_file,chain_id,partner_in_complex,partner_chain,restraints_2,condition)
          if add_ppi_restraints:
             print("\t\t\t-- Add restrain contacts (%s %s) %d"%(protein_num,partner_in_complex,len([x for x in restraints_2.keys()])))
             new_restraints  = new_restraints +len([x for x in restraints_2.keys()])
             self.add_restraints(restraints_2,protein_num,partner_in_complex)
             self.add_interfaces(interface_2,interface_1,protein_num,partner_in_complex)
       else:
         if right_oriented:
          #new_protein.split(":")[0] == protein_1:
          protein_num=self.add_protein(protein_1,sequence_1,partner_in_complex)
          print("\t\t\t-- Add new protein: "+protein_num)
          chain_id      = protein_chain[protein_1+":A"]
          partner_chain = protein_chain[protein_2+":B"]
          add_ppi_restraints=self.add_topology_and_connect(0,protein_num,pdb_file,chain_id,partner_in_complex,partner_chain,restraints_1)
          print("\t\t\t-- Add restrain contacts (%s %s) %d"%(protein_num,partner_in_complex,len([x for x in restraints_1.keys()])))
          new_restraints  = new_restraints +len([x for x in restraints_1.keys()])
          self.add_restraints(restraints_1,protein_num,partner_in_complex)
          self.add_interfaces(interface_1,interface_2,protein_num,partner_in_complex)
          self.add_chain_correspondence_plus(protein_num,pdb_file,chain_id)
         else:
          #new_protein.split(":")[0] == protein_2: 
          protein_num=self.add_protein(protein_2,sequence_2,partner_in_complex)
          print("\t\t\t-- Add new protein: "+protein_num)
          chain_id      = protein_chain[protein_2+":B"]
          partner_chain = protein_chain[protein_1+":A"]
          add_ppi_restraints=self.add_topology_and_connect(0,protein_num,pdb_file,chain_id,partner_in_complex,partner_chain,restraints_2)
          print("\t\t\t-- Add restrain contacts (%s %s) %d"%(protein_num,partner_in_complex,len([x for x in restraints_2.keys()])))
          new_restraints  = new_restraints +len([x for x in restraints_2.keys()])
          self.add_restraints(restraints_2,protein_num,partner_in_complex)
          self.add_interfaces(interface_2,interface_1,protein_num,partner_in_complex)
          self.add_chain_correspondence_plus(protein_num,pdb_file,chain_id)
       new_complex_list=[]
       for i in range(1,len(adding)):
           new_protein,partner_in_complex,right_oriented  = adding[i]
           new_name = self.get_name() + "."+str(i)
           new_complex = seed_complex.copy(new_name)
           #check = self.get_interface(partner_in_complex)
           print("\t\t  --check adding %s - %s on a new complex %s"%(new_protein,partner_in_complex,new_name))
           if condition == "old":
            #check = self.get_interface(partner_in_complex)
            if right_oriented:
              #new_protein.split(":")[0] == protein_1:
              protein_num   = new_protein
              chain_id      = protein_chain[protein_1+":A"]
              partner_chain = protein_chain[protein_2+":B"]
              add_ppi_restraints=new_complex.add_topology_and_connect(i,protein_num,pdb_file,chain_id,partner_in_complex,partner_chain,restraints_1,condition)
              if add_ppi_restraints:
                 print("\t\t\t-- Add restrain contacts (%s %s) %d"%(protein_num,partner_in_complex,len([x for x in restraints_1.keys()])))
                 new_restraints  = new_restraints +len([x for x in restraints_1.keys()])
                 new_complex.add_restraints(restraints_1,protein_num,partner_in_complex)
                 new_complex.add_interfaces(interface_1,interface_2,protein_num,partner_in_complex)
            else:
              #new_protein.split(":")[0] == protein_2: 
              protein_num   = new_protein
              chain_id      = protein_chain[protein_2+":B"]
              partner_chain = protein_chain[protein_1+":A"]
              add_ppi_restraints=new_complex.add_topology_and_connect(i,protein_num,pdb_file,chain_id,partner_in_complex,partner_chain,restraints_2,condition)
              if add_ppi_restraints:
                 print("\t\t\t-- Add restrain contacts (%s %s) %d"%(protein_num,partner_in_complex,len([x for x in restraints_2.keys()])))
                 new_restraints  = new_restraints +len([x for x in restraints_2.keys()])
                 new_complex.add_restraints(restraints_2,protein_num,partner_in_complex)
                 new_complex.add_interfaces(interface_2,interface_1,protein_num,partner_in_complex)
           else:
            if right_oriented:
              #new_protein.split(":")[0] == protein_1:
              protein_num=new_complex.add_protein(protein_1,sequence_1,partner_in_complex)
              print("\t\t\t-- Add new protein: "+protein_num)
              chain_id      = protein_chain[protein_1+":A"]
              partner_chain = protein_chain[protein_2+":B"]
              add_ppi_restraints=new_complex.add_topology_and_connect(i,protein_num,pdb_file,chain_id,partner_in_complex,partner_chain,restraints_1)
              print("\t\t\t-- Add restrain contacts (%s %s) %d"%(protein_num,partner_in_complex,len([x for x in restraints_1.keys()])))
              new_restraints  = new_restraints +len([x for x in restraints_1.keys()])
              new_complex.add_restraints(restraints_1,protein_num,partner_in_complex)
              new_complex.add_interfaces(interface_1,interface_2,protein_num,partner_in_complex)
              new_complex.add_chain_correspondence_plus(protein_num,pdb_file,chain_id)
            else:
              # new_protein.split(":")[0] == protein_2: 
              protein_num=new_complex.add_protein(protein_2,sequence_2,partner_in_complex)
              print("\t\t\t-- Add new protein: "+protein_num)
              chain_id      = protein_chain[protein_2+":B"]
              partner_chain = protein_chain[protein_1+":A"]
              add_ppi_restraints=new_complex.add_topology_and_connect(i,protein_num,pdb_file,chain_id,partner_in_complex,partner_chain,restraints_2)
              print("\t\t\t-- Add restrain contacts (%s %s) %d"%(protein_num,partner_in_complex,len([x for x in restraints_2.keys()])))
              new_restraints  = new_restraints +len([x for x in restraints_2.keys()])
              new_complex.add_restraints(restraints_2,protein_num,partner_in_complex)
              new_complex.add_interfaces(interface_2,interface_1,protein_num,partner_in_complex)
              new_complex.add_chain_correspondence_plus(protein_num,pdb_file,chain_id)
           new_complex_list.append( new_complex )

       done_complex=[]
       for c in new_complex_list:
           skip = False
           for d in done_complex:
               if skip: continue
               if c==d: skip=True
           if skip: continue
           done_complex.append(c)
               
       return (done_complex,new_restraints)

   def restraints_to_csv(self,interaction_type="all",kd=20.0,weight=1.0,sd=1.0):
       data={}
       for key in self._restraints:
         posi,posj = key
         dist,ratio= self._restraints[key]
         name_1,res_1=posi
         name_2,res_2=posj
         if name_1 in self._protein_chains: mtype_1="prot"
         if name_2 in self._protein_chains: mtype_2="prot"
         if interaction_type == mtype_1+mtype_2 or interaction_type == mtype_2+mtype_1 or interaction_type=="all":
            data.setdefault("prot1",[]).append(name_1)
            data.setdefault("prot2",[]).append(name_2)
            data.setdefault("prot1_res",[]).append(res_1)
            data.setdefault("prot2_res",[]).append(res_2)
            data.setdefault("distance",[]).append(float(dist))
            data.setdefault("kd",[]).append(float(ratio)*kd)
            data.setdefault("weight",[]).append(float(weight))
            data.setdefault("sd",[]).append(float(sd))
       table=pd.DataFrame(data)
       return table

   def contains(self,other):
       proteins_a = set([protein.split(":")[0] for protein in self.get_protein_chains()])
       proteins_b = set([protein.split(":")[0] for protein in other.get_protein_chains()])
       # chek all proteins in 'other' are in 'self'
       if not proteins_b.issubset(proteins_a): return False
       set_a = {}
       set_b = {}
       # chek the number of a proteins of a protein type in 'other' is the same or less than in 'self'
       skip=False
       for p in proteins_a: 
           set_a.setdefault(p,set([protein for protein in self.get_protein_chains() if protein.split(":")[0] == p]))
       for p in proteins_b: 
           set_b.setdefault(p,set([protein for protein in other.get_protein_chains() if protein.split(":")[0] == p]))
           if not set_a.has_key(p):
              skip=True
           else:
              if len(set_a[p])<len(set_b[p]): skip=True
       if skip: return False
       # check the correspondence of all proteins in other with the same reference_site in self
       correspondence_b={}
       correspondence_a={}
       for p in proteins_b:
           for b in set_b[p]:
               site_b = other.get_reference_site(b)
               if site_b is None: continue
               all_correspondence={}
               for a in set_a[p]:
                   site_a = self.get_reference_site(a)
                   if site_a is None: continue
                   common = site_a.intersection(site_b)
                   if len(common) == 0: continue
                   all_correspondence.setdefault(a, len(common))
               a_rank = [a for a,y in sorted(all_correspondence.iteritems(), key=lambda x: x[1],reverse=True )]
               if len(a_rank) <=0: return False
               a_max=a_rank[0]
               correspondence_b.setdefault(b,a_max)
               correspondence_a.setdefault(a_max,b)
       # check other has at least the same restraints as self
       contain_restraints = True
       if self._restrict_equal:
         restraints_a  = set([x for x in self.get_restraints().keys()])
         orestraints_b = set([x for x in other.get_restraints().keys()])
         restraints_b  = set()
         for x,y in orestraints_b:
           protx,aax = x
           if correspondence_b.has_key(protx): 
               xx = correspondence_b[protx]
           else:
               continue
           proty,aay = y
           if correspondence_b.has_key(proty): 
               yy = correspondence_b[proty]
           else:
               continue
           restraints_b.add( tuple( (tuple( (xx,aax) ),tuple( (yy,aay) )) ) )
         for x,y in restraints_a:
           if contain_restraints:
              if (x,y) not in restraints_b and (y,x) not in restraints_b: 
                 contain_restraints = False
       # Check conditions
       if  set(correspondence_b.keys()) == set(other.get_protein_chains()) and len(correspondence_a.keys()) == len(self.get_protein_chains()) and set([x for k,x in correspondence_a.iteritems()]) == set(other.get_protein_chains()) and contain_restraints: 
              return True
       return False

   def __hash__(self):
       return self._id.__hash__()
                   
   def __eq__(self,other):
       return self.contains(other) and other.contains(self)

   def __lt__(self,other):
       return len(self._protein_chains) < len(other._protein_chains)
       
   def __gt__(self,other):
       return len(self._protein_chains) > len(other._protein_chains)
       
   def __cmp__(self,other):
      if   self == other: return 0
      elif self < other:  return -1
      else:               return 1
      
   def write_restraints(self,output,interaction_type="all",kd=20.0,weight=1.0,sd=1.0):
      if os.path.exists(output) and os.path.isfile(output):
       fo=open(output,"a")
      else:
       fo=open(output,"w")
       fo.write("prot1,prot1_res,prot2,prot2_res,distance,kd,weight,sd\n")
      for key in self._restraints:
        posi,posj = key
        dist,ratio= self._restraints[key]
        name_1,res_1=posi
        name_2,res_2=posj
        if name_1 in self._protein_chains: mtype_1="prot"
        if name_2 in self._protein_chains: mtype_2="prot"
        if interaction_type == mtype_1+mtype_2 or interaction_type == mtype_2+mtype_1 or interaction_type=="all":
           fo.write("%s,%s,%s,%s,%10.3f,%10.3f,%10.3f,%10.3f\n"%(name_1,res_1,name_2,res_2,float(dist),float(ratio)*kd,weight,sd))
      fo.close()

   def get_topology_csv(self):
      table=pd.DataFrame(self._topology)
      return table

   def set_topology_core(self):
      data={}
      #red_color   = ["darkred", "firebrick","tomato","red", "orangered", "coral", "orange", "goldenrod","gold","siena","yellow"]
      #blue_color  = ["navy","darkblue","blue","midnightblue","dodgerblue","steelblue","deepskyblue","skyblue","lightblue","cyan","lightcyan"]
      #green_color = ["darkgreen","forestgreen","green","lawngreen","springgreen","lightgreen","limegreen","turquoise","aquamarine"]
      #color = red_color + blue_color + green_color
      fasta_name = self._id+".topology.fasta"
      color_molecule = 0
      for protein in self._protein_chains:
          color_molecule = color_molecule + 1
          pdb_chain  = self.get_pdb_and_chain_in_core(protein)
          if pdb_chain is None: continue
          pdb_file,c = pdb_chain
          pdb_cif    = pdb_file.rstrip("pdb")+"cif"
          pdb_obj    = PDB(pdb_file)
          chain      = pdb_obj.get_chain_by_id(c)
          start      = chain.first_aminoacid.identifier
          end        = chain.last_aminoacid.identifier
          offset     = 0
          data.setdefault("molecule_name",[]).append(protein)
          #data.setdefault("color",[]).append(color[np.mod(color_molecule, len(color))])
          data.setdefault("color",[]).append("blue")
          data.setdefault("fasta_fn",[]).append(os.path.basename(fasta_name))
          data.setdefault("fasta_id",[]).append(protein)
          data.setdefault("pdb_fn",[]).append(os.path.basename(pdb_cif))
          data.setdefault("chain",[]).append(c)
          data.setdefault("residue_range",[]).append(str(start)+","+str(end))
          data.setdefault("pdb_offset",[]).append(offset)
          data.setdefault("bead_size",[]).append(1)
          data.setdefault("em_residues_per_gaussian",[]).append(0)
          if data.has_key("rigid_body"): rb=len(data.get("rigid_body")) + 1
          else: rb=1
          data.setdefault("rigid_body",[]).append(rb)
          data.setdefault("super_rigid_body",[]).append(" ")
          data.setdefault("chain_of_super_rigid_bodies",[]).append(" ")
      self._topology=data

   def clean_topology(self,fasta_file):
       """
       clean topology by removing isolated proteins unconnected with the complex and without restraints
       missing protein sequences are recovered from the original FastA file with all the sequences 
       """
       sequences={}
       for protein,sequence in parse_fasta_file(fasta_file,gz=False,clean=True):
          sequences.setdefault(protein,sequence)

       molecule_name                = self._topology["molecule_name"]
       color                        = self._topology["color"]
       fasta_fn                     = self._topology["fasta_fn"]
       fasta_id                     = self._topology["fasta_id"]
       pdb_fn                       = self._topology["pdb_fn"]
       chain                        = self._topology["chain"]
       residue_range                = self._topology["residue_range"]
       pdb_offset                   = self._topology["pdb_offset"]
       bead_size                    = self._topology["bead_size"]
       em_residues_per_gaussian     = self._topology["em_residues_per_gaussian"]
       rigid_body                   = self._topology["rigid_body"]
       super_rigid_body             = self._topology["super_rigid_body"]
       chain_of_super_rigid_bodies  = self._topology["chain_of_super_rigid_bodies"]

       data={}
       for i in range(len(molecule_name)):
           skip = True
           skip_protein=False
           if molecule_name[i] not in self._dna_chains:
             if molecule_name[i] not in self._protein_chains:
               code = molecule_name[i].split(":")[0]
               print("\t\t-- Missing protein "+molecule_name[i])
               if sequences.has_key(code):
                  self._protein_chains.append(molecule_name[i])
                  self._protein_sequences.append(sequences[code])
               else:
                  skip_protein=True
             if skip_protein:
               print("\t\t-- SKIP protein "+molecule_name[i])
               delete_restraints=[]
               for pair,restr in self._restraints.iteritems():
                 prot_pos1,prot_pos2 = pair
                 p1,aa1 = prot_pos1
                 p2,aa2 = prot_pos2
                 if skip_protein:
                  if (p1 == molecule_name[i] or p2 == molecule_name[i]): delete_restraints.append(pair)
               for pair in delete_restraints:
                 if pair in self._restraints.keys():
                   del self._restraints[pair]
           for pair,restr in self._restraints.iteritems():
               prot_pos1,prot_pos2 = pair
               p1,aa1 = prot_pos1
               p2,aa2 = prot_pos2
               if not skip: continue
               if (p1 == molecule_name[i] and p2 != molecule_name[i]) or (p1 != molecule_name[i] and p2 == molecule_name[i]):
                   skip = False
                   break
           if skip: continue
           if skip_protein: continue 
           data.setdefault("molecule_name",[]).append(molecule_name[i])
           data.setdefault("color",[]).append(color[i])
           data.setdefault("fasta_fn",[]).append(fasta_fn[i])
           data.setdefault("fasta_id",[]).append(fasta_id[i])
           data.setdefault("pdb_fn",[]).append(pdb_fn[i])
           data.setdefault("chain",[]).append(chain[i])
           data.setdefault("residue_range",[]).append(residue_range[i])
           data.setdefault("pdb_offset",[]).append(pdb_offset[i])
           data.setdefault("bead_size",[]).append(bead_size[i])
           data.setdefault("em_residues_per_gaussian",[]).append(em_residues_per_gaussian[i])
           if data.has_key("rigid_body"): rb=len(data.get("rigid_body")) + 1
           else: rb=1
           data.setdefault("rigid_body",[]).append(rb)
           data.setdefault("super_rigid_body",[]).append(super_rigid_body[i])
           data.setdefault("chain_of_super_rigid_bodies",[]).append(chain_of_super_rigid_bodies[i])


       self._topology=data

       

   def complete_topology(self,bead_size=5,domain=None):
       """
       Complete the misisng regions (loops and tails) without PDB known structure using BEADS
       'domain' defines the missing regions:
              loops   missing regions of the protein chain wthout considering Nt and Ct tails
              tails   missing regions at Nt and Ct tails of the protein sequence
              full    all missing regions (loops and tails)
       if domain is None (default) then all missing regions are left unknown without representation in the model.
       returns the total of regions modelled as beads without known structure in PDB
       """
       if domain is None: return 0
       fasta_name = self._id+".topology.fasta"

       if domain == "loops":
          sum_gaps=0
          for protein in self._protein_chains:
             residues=set()
             for i in range(len(self._topology["molecule_name"])):
               if self._topology["molecule_name"][i] == protein:
                  try:
                    residue_range = self._topology["residue_range"][i]
                    rr_start,rr_end = residue_range.split(",")
                    for aa in range(int(rr_start),int(rr_end)+1):
                      residues.add(aa)
                  except Exception as e:
                    print("Skip protein "+str(protein))
                    print("ERROR: %s"%str(e))
                    continue
             all_residues = [x for x in residues]
             all_residues.sort()
             gaps=[]
             for i in range(len(all_residues)):
               if i+1 <len(all_residues):
                  if all_residues[i+1] - all_residues[i] > 1:
                     start = all_residues[i]+1
                     end   = all_residues[i+1]-1
                     gaps.append(tuple((start,end)))
             offset=0
             for gap in gaps:
                 sum_gaps = sum_gaps +1
                 gstart,gend = gap
                 gbead_size  = bead_size
                 if gend-gstart <bead_size: gbead_size = int(float(gend-gstart)/2)
                 if gbead_size < 1: gbead_size=1
                 for j in range(gstart-1,gend,gbead_size):
                     start= j+1
                     end  = j+gbead_size
                     if end > gend: end=gend
                     self._topology.setdefault("molecule_name",[]).append(protein)
                     self._topology.setdefault("color",[]).append("green")
                     self._topology.setdefault("fasta_fn",[]).append(os.path.basename(fasta_name))
                     self._topology.setdefault("fasta_id",[]).append(protein)
                     self._topology.setdefault("pdb_fn",[]).append("BEADS")
                     self._topology.setdefault("chain",[]).append(" ")
                     self._topology.setdefault("residue_range",[]).append(str(start)+","+str(end))
                     self._topology.setdefault("pdb_offset",[]).append(offset)
                     self._topology.setdefault("bead_size",[]).append(gbead_size)
                     self._topology.setdefault("em_residues_per_gaussian",[]).append(0)
                     if self._topology.has_key("rigid_body"): rb=len(self._topology.get("rigid_body")) + 1
                     else: rb=1
                     self._topology.setdefault("rigid_body",[]).append(rb)
                     self._topology.setdefault("super_rigid_body",[]).append(" ")
                     self._topology.setdefault("chain_of_super_rigid_bodies",[]).append(" ")
          return (sum_gaps) 

       elif  domain == "tails":
          sum_gaps=0
          for j in range(len(self._protein_chains)):
             sequence = self._protein_sequences[j]
             protein  = self._protein_chains[j]
             residues=set()
             for i in range(len(self._topology["molecule_name"])):
               if self._topology["molecule_name"][i] == protein:
                  try:
                    residue_range = self._topology["residue_range"][i]
                    rr_start,rr_end = residue_range.split(",")
                    for aa in range(rr_start,rr_end+1):
                      residues.add(aa)
                  except Exception as e:
                    print("Skip protein "+str(protein))
                    print("ERROR: %s"%str(e))
                    continue
             all_residues = [x for x in residues]
             all_residues.sort()
             gaps=[]
             gaps.append(tuple((1,min([int(x) for x in residues])-1)))
             gaps.append(tuple((len(sequence)+1,max([int(x) for x in residues]))))
             offset=0
             for gap in gaps:
                 sum_gaps = sum_gaps +1
                 gstart,gend = gap
                 gbead_size  = bead_size
                 if gend-gstart <bead_size: gbead_size = int(float(gend-gstart)/2)
                 if gbead_size < 1: gbead_size=1
                 for j in range(gstart-1,gend,gbead_size):
                     start= j+1
                     end  = j+gbead_size
                     if end > gend: end=gend
                     self._topology.setdefault("molecule_name",[]).append(protein)
                     self._topology.setdefault("color",[]).append("green")
                     self._topology.setdefault("fasta_fn",[]).append(os.path.basename(fasta_name))
                     self._topology.setdefault("fasta_id",[]).append(protein)
                     self._topology.setdefault("pdb_fn",[]).append("BEADS")
                     self._topology.setdefault("chain",[]).append(" ")
                     self._topology.setdefault("residue_range",[]).append(str(start)+","+str(end))
                     self._topology.setdefault("pdb_offset",[]).append(offset)
                     self._topology.setdefault("bead_size",[]).append(gbead_size)
                     self._topology.setdefault("em_residues_per_gaussian",[]).append(0)
                     if self._topology.has_key("rigid_body"): rb=len(self._topology.get("rigid_body")) + 1
                     else: rb=1
                     self._topology.setdefault("rigid_body",[]).append(rb)
                     self._topology.setdefault("super_rigid_body",[]).append(" ")
                     self._topology.setdefault("chain_of_super_rigid_bodies",[]).append(" ")
          return (sum_gaps) 

       elif  domain == "full":
          sum_gaps=0
          for j in range(len(self._protein_chains)):
             #Tails
             sequence = self._protein_sequences[j]
             protein  = self._protein_chains[j]
             residues=set()
             for i in range(len(self._topology["molecule_name"])):
               if self._topology["molecule_name"][i] == protein:
                  try:
                    residue_range = self._topology["residue_range"][i]
                    rr_start,rr_end = residue_range.split(",")
                    for aa in range(int(rr_start),int(rr_end)+1):
                      residues.add(aa)
                  except Exception as e:
                    print("Skip protein "+str(protein))
                    print("ERROR: %s"%str(e))
                    continue
             all_residues = [x for x in residues]
             all_residues.sort()
             gaps=[]
             gaps.append(tuple((1,min([int(x) for x in residues])-1)))
             gaps.append(tuple((len(sequence)+1,max([int(x) for x in residues]))))
             #Loops
             for i in range(len(all_residues)):
               if i+1 <len(all_residues):
                  if all_residues[i+1] - all_residues[i] > 1:
                     start = all_residues[i]+1
                     end   = all_residues[i+1]-1
                     gaps.append(tuple((start,end)))
             offset=0
             for gap in gaps:
                 sum_gaps = sum_gaps +1
                 gstart,gend = gap
                 gbead_size  = bead_size
                 if gend-gstart <bead_size: gbead_size = int(float(gend-gstart)/2) 
                 if gbead_size < 1: gbead_size=1
                 for j in range(gstart-1,gend,gbead_size):
                     start= j+1
                     end  = j+gbead_size
                     if end > gend: end=gend
                     self._topology.setdefault("molecule_name",[]).append(protein)
                     self._topology.setdefault("color",[]).append("green")
                     self._topology.setdefault("fasta_fn",[]).append(os.path.basename(fasta_name))
                     self._topology.setdefault("fasta_id",[]).append(protein)
                     self._topology.setdefault("pdb_fn",[]).append("BEADS")
                     self._topology.setdefault("chain",[]).append(" ")
                     self._topology.setdefault("residue_range",[]).append(str(start)+","+str(end))
                     self._topology.setdefault("pdb_offset",[]).append(offset)
                     self._topology.setdefault("bead_size",[]).append(gbead_size)
                     self._topology.setdefault("em_residues_per_gaussian",[]).append(0)
                     if self._topology.has_key("rigid_body"): rb=len(self._topology.get("rigid_body")) + 1
                     else: rb=1
                     self._topology.setdefault("rigid_body",[]).append(rb)
                     self._topology.setdefault("super_rigid_body",[]).append(" ")
                     self._topology.setdefault("chain_of_super_rigid_bodies",[]).append(" ")
          return (sum_gaps) 
       else:
           return 0
       
   def write_topology(self,output_dir,fasta_file,new=False):
      #clean the topology by trimming unconnected proteins
      self.clean_topology(fasta_file)
      fasta_name=self._id+".topology.fasta"
      if not os.path.exists(os.path.join(output_dir,fasta_name)):
         fo = open(os.path.join(output_dir,fasta_name),"w")
         for i in range(len(self._dna_chains)):
             fo.write(">%s\n%s\n"%(self._dna_chains[i],self._dna_sequences[i]))
         for i in range(len(self._protein_chains)):
             fo.write(">%s\n%s\n"%(self._protein_chains[i],self._protein_sequences[i]))
         fo.close()
      #rename the fasta file accordingly in the topology 
      new_fasta_list=[]
      for t in self._topology["fasta_fn"]:
          new_fasta_list.append(fasta_name)
      self._topology["fasta_fn"]=new_fasta_list
      #get topology in CSV table format
      topology = self.get_topology_csv()
      output   = os.path.join(output_dir,self.get_name()+".topology.txt")
      rows     = topology.shape[0]
      columns  =["molecule_name","color","fasta_fn","fasta_id","pdb_fn","chain","residue_range","pdb_offset","bead_size","em_residues_per_gaussian","rigid_body","super_rigid_body","chain_of_super_rigid_bodies"]
      maxlen={}
      for col in columns:
          maxlen.setdefault(col, len(str(col))+2)
          if maxlen[col] < max( [len(str(topology.iloc[i][col])) for i in range(rows)] ): maxlen[col] = max( [len(str(topology.iloc[i][col])) for i in range(rows)] )
      if os.path.exists(output) and os.path.isfile(output) and not new:
       fo=open(output,"a")
       for i in range(rows):
           line="|"
           for col in columns:
               #print("CHECK "+col+" "+ str(topology.iloc[i][col])+" size "+str(maxlen[col])+" in line "+line)
               if maxlen[col]< 20: line=line+"{:20}|".format(topology.iloc[i][col])
               if maxlen[col]>=20 and maxlen[col] < 30: line=line+"{:30}|".format(topology.iloc[i][col])
               if maxlen[col]>=30 and maxlen[col] < 50: line=line+"{:50}|".format(topology.iloc[i][col])
               if maxlen[col]>=50 and maxlen[col] < 70: line=line+"{:70}|".format(topology.iloc[i][col])
               if maxlen[col]>=70: line=line+"{:100}|".format(topology.iloc[i][col])
           fo.write("%s\n"%line)
      else:
       fo=open(output,"w")
       line="|"
       for col in columns:
           if maxlen[col]< 20: line=line+"{:20}|".format(col)
           if maxlen[col]>=20 and maxlen[col] < 30: line=line+"{:30}|".format(col)
           if maxlen[col]>=30 and maxlen[col] < 50: line=line+"{:50}|".format(col)
           if maxlen[col]>=50 and maxlen[col] < 70: line=line+"{:70}|".format(col)
           if maxlen[col]>=70: line=line+"{:100}|".format(col)
       fo.write("%s\n"%line)
       for i in range(rows):
           line="|"
           for col in columns:
               if maxlen[col] < 20: line=line+"{:20}|".format(topology.iloc[i][col])
               if maxlen[col]>=20 and maxlen[col] < 30: line=line+"{:30}|".format(topology.iloc[i][col])
               if maxlen[col]>=30 and maxlen[col] < 50: line=line+"{:50}|".format(topology.iloc[i][col])
               if maxlen[col]>=50 and maxlen[col] < 70: line=line+"{:70}|".format(topology.iloc[i][col])
               if maxlen[col]>=70: line=line+"{:100}|".format(topology.iloc[i][col])
           fo.write("%s\n"%line)

      fo.close()


class TF_complex(macro_complex):
   """
   Macro-complex class that inlcudes DNA

   This is a class with information of all chains, proteins and DNA forming a complex (i.e. transcription complex)
   It is constructed with PDB files of fragmented DNA and a set of PDB complexes

   Initialize:
   name	= name of the complex (mandatory)
   fragments = DNA fragments of maximum length 250bp (empty list otherwise [])
   dna_files = pdb files with DNA coordiantes for each fragment (empty list otherwise [])
   complex_files = dictionary of PDB coordinates and chain assignement of a complex with DNA for each fragment (empty dictionary otherwise {})

   Additional PDB files can be incorporated to define a protein after the core of the complex is defined. 
   Added files will be defined by the _chain_correspondence_plus dictionary for new and previous proteins of the core
   The size of DNA bead is used as thgeshold to calculate the limit of maximum strength with th erestraints

   Variables

   	_id		is the name given to the complex (initialized with "name")
	
	_fragments	list of tuples for fragments of DNA (i.e 1-250, 226-500, 475- etc.)
	
	_dna_files 	list of PDB files with DNA structures, each corresponds with a fragment	(i.e dna_1-250.pdb is the structure of DNA in fragtament 1-250)

	_complexes	dictionary of { fragment:tuple([PDB,remarks]) } PDB,remarks are two unique files on a DNA fragment to built the complex, 
			The PDB file contains proteins and DNA in the fragment, the remarks show the name of the protein associated with each chain in the PDB
			Fore example:

				REMARK	a = DNA
				REMARK	b = DNA
				REMARK	c = P11474.3cbb_B.96-115:0:75_3cbb_B	(Model name of P11474 modelled with 3cbb_B as template)
				REMARK	d = P11474.4cn3_C.30-43:0:77_4cn3_C
				REMARK	e = P11474.1cit_A.24-39:0:86_1cit_A
				REMARK	f = P11474.4hn6_B.186-201:0:67_4hn6
				REMARK	g = P15976.4hc7_A.129-140:0:104_4hc7
 

	_protein_chains		list of protein chains. Proteins with the same code of UniProt are numbered ( :1, :2, :3, etc.)	to distinguish
	
	_protein_sequences 	list of protein sequences in the same order as _protein_chains

	_dna_chains		list of DNA chains. 
	
	_dna_sequences 	 	list of DNA sequences in the same order as _dna_chains

	_restraints	dictionary of restraints of the form { tuple( (pos_1,pos_2) ): (distance,1.0) }, where pos_1 is also a tuple: (molecule, number), with molecule a protein or DNA chain

	_interfaces	dictionary of interfaces  of the form { molecule: set([numbers])}, with molecule a protein or DNA chain and the set of residue numbers with contacts with other molecules





   """
   def __init__(self,name,fragments=None,dna_files=None,complex_files=None,PPI_distance=12):

       super(TF_complex,self).__init__(name,fragments,complex_files,PPI_distance)
       if dna_files is not None: self._dna_files = dna_files
       else:                     self._dna_files = []
       self._dna_chains=[]    
       self._dna_sequences=[]

       #Internal parameters
           
       # restraints for DNA

       self._dna_sdsize  = 15
       self._dna_xdsize  = 15
       self._dna_sfactor = 3.0
       self._dna_xfactor = 2.0
       self._dna_bead    = 10
       self._dna_maxdist = 100.0

   def initialize(self,fasta_file):
       print("\t--Get protein chains")
       self.set_proteins_in_core(fasta_file)
       print("\t--Get DNA chains")
       self.set_dna()
       self.set_interfaces_and_restraints_in_core()
       self.set_DNA_restraints()
       self.set_topology_core()


   
   def copy_dna_files(self,dna_files):
       self.dna_files=dna_files[:]


   def copy_dna_chains(self,dna_chains):
       self._dna_chains = dna_chains[:]

   def copy_dna_sequences(self,dna_sequences):
       self._dna_sequences = dna_sequences[:]

   def set_dna_sdsize(self,dna_sdsize):
       self._dna_sdsize = dna_sdsize

   def set_dna_xdsize(self,dna_xdsize):
       self._dna_xdsize = dna_xdsize

   def set_dna_xfactor(self,dna_xfactor):
       self._dna_xfactor = dna_xfactor

   def set_dna_sfactor(self,dna_sfactor):
       self._dna_sfactor = dna_sfactor

   def set_dna_bead(self,dna_bead):
       self._dna_bead = dna_bead

   def set_dna_maxdist(self,dna_maxdist):
       self._dna_maxdist = dna_maxdist

   def set_dna(self):
       dna_chains=[]
       dna_sequences={}
       for i in range(len(self._fragments)):
         dna_file = self._dna_files[i]
         fragment_a,fragment_b = self._fragments[i]
         #print("Read DNA in MACRO "+dna_file)
         pdb_obj = PDB(dna_file)
         chain_ids=[]
         for chain in pdb_obj.chains:
            if chain.chaintype=="N":
             chain_id = chain.chain
             chain_ids.append(chain_id)
         chain_ids.sort()
         forward = chain_ids[0]
         reverse = chain_ids[1]
         for c in chain_ids:
             dna_code = "DNA"+c
             print("\t\t--get chain "+c+" code "+dna_code)
             chain=pdb_obj.get_chain_by_id(c)
             if dna_sequences.has_key(dna_code):
              overlap = len( dna_sequences[dna_code] ) - fragment_a + 1
              if c == forward:
                 dna_sequences[dna_code] = dna_sequences[dna_code] + chain.nucleotide_sequence()[overlap:]
              if c == reverse:
                 dna_sequences[dna_code] = chain.nucleotide_sequence()[:-overlap] +  dna_sequences[dna_code]
             else:
              dna_chains.append(dna_code)
              dna_sequences.setdefault(dna_code,chain.nucleotide_sequence())
       self._dna_chains=dna_chains
       self._dna_sequences=[]
       for dna_code in dna_chains:
           self._dna_sequences.append(dna_sequences[dna_code])
    
   def set_protein_reference_site(self,reference_site):
       for k,p in reference_site.iteritems():
           self._protein_reference_site.setdefault(k,p)

   def set_reference_site(self,protein,site):
       self._protein_reference_site.setdefault(protein,site)

   def set_interfaces_and_restraints_in_core(self):
       #this method initializes the set of restraints
       for fragment in self._fragments:
           print("\t--Get restraints PPI and TF-DNA in fragment "+str(fragment))
           pdb_file    = self.get_pdb_by_fragment(fragment)
           print("\t--PDB File "+os.path.basename(pdb_file))
           pdb_obj     = PDB(pdb_file)
           pdb_complex = Complex(pdb_obj,PPI_distance=self._PPI_distance)
           for interface in pdb_complex.PNInterfaces:
               protein_chain_id = interface.protein.chain
               print("\t\t-- Interface protein-DNA of chain "+protein_chain_id)
               dna_chain_id     = "DNA"+interface.nucleotide.chain
               protein_num      = self.get_protein_by_fragment_and_chain(fragment,protein_chain_id)
               if protein_num is None: continue
               print("\t\t-- protein "+protein_num)
               for contact in interface.contacts:
                   aminoacid  = contact.aminoacid.identifier
                   self._interfaces.setdefault(protein_num,set()).add(aminoacid)
                   nucleotide = contact.nucleotide.identifier
                   vector     = contact.aminoacid.ca.coordinates - contact.nucleotide.p.coordinates
                   distance   = np.sqrt(vector.dot(vector))
                   pos_prot   = tuple( (protein_num,aminoacid) )
                   pos_dna    = tuple( (dna_chain_id,nucleotide) )
                   self._protein_reference_site.setdefault(protein_num,set()).add(pos_dna)
                   self._restraints.setdefault(tuple( (pos_prot,pos_dna) ) ,(distance,1.0) )
           for interface in pdb_complex.PPInterfaces:
               protein_chain_id_1 = interface.protein_chain.chain
               protein_chain_id_2 = interface.protein_interactor.chain
               print("\t\t-- Interface protein-protein of chains "+protein_chain_id_1+" "+protein_chain_id_2)
               protein_num_1      = self.get_protein_by_fragment_and_chain(fragment,protein_chain_id_1)
               protein_num_2      = self.get_protein_by_fragment_and_chain(fragment,protein_chain_id_2)
               if protein_num_1 is None: continue
               if protein_num_2 is None: continue
               if self.get_reference_site(protein_num_1) is not None and self.get_reference_site(protein_num_2) is None:
                  self._protein_reference_site.setdefault(protein_num_2,self.get_reference_site(protein_num_1))
               if self.get_reference_site(protein_num_2) is not None and self.get_reference_site(protein_num_1) is None:
                  self._protein_reference_site.setdefault(protein_num_1,self.get_reference_site(protein_num_2))
               print("\t\t-- PPI "+protein_num_1+" with "+protein_num_2)
               for contact in interface.contacts:
                   aminoacid1  = contact.aminoacid1.identifier
                   aminoacid2  = contact.aminoacid2.identifier
                   self._interfaces.setdefault(protein_num_1,set()).add(aminoacid1)
                   self._interfaces.setdefault(protein_num_2,set()).add(aminoacid2)
                   vector     = contact.aminoacid2.ca.coordinates - contact.aminoacid1.ca.coordinates
                   distance   = np.sqrt(vector.dot(vector))
                   pos_prot_1 = tuple( (protein_num_1,aminoacid1) )
                   pos_prot_2 = tuple( (protein_num_2,aminoacid2) )
                   self._restraints.setdefault(tuple( (pos_prot_1,pos_prot_2) ) ,(distance,1.0) )

   def set_DNA_restraints(self,fragment=None,sdsize=None,xdsize=None,sfactor=None,xfactor=None,threshold=None,maxdist=None):
       if sdsize is None: sdsize=self._dna_sdsize
       if xdsize is None: xdsize=self._dna_xdsize
       if sfactor is None: sfactor=self._dna_sfactor
       if xfactor is None: xfactor=self._dna_xfactor
       if threshold is None: threshold=self._dna_bead
       if maxdist is None: maxdist=self._dna_maxdist

       #Get the standard DNA file to built constraints
       select = 0
       if fragment is not None:
         i=0
         for f in self._fragments:
             if f == fragment:
                 select = i
                 break
             i = i + 1
       dna_size=0
       for sequence in self._dna_sequences: dna_size = max(dna_size,len(sequence))
       dna_file = self._dna_files[select]
       pdb_obj  = PDB(dna_file)
       #Get coordiantes of phosphates
       phos_pos={}
       phos_xyz={}
       for dna_code in self._dna_chains:
           c=dna_code.lstrip("DNA")
           chain =pdb_obj.get_chain_by_id(c)
           for n in chain.nucleotides:
               phos_pos.setdefault(dna_code,[]).append(n.identifier)
               phos_xyz.setdefault(tuple((dna_code,n.identifier)),n.p.coordinates)
       #Generate forces dictionary
       self_kforce={}
       x_kforce={}
       for i in range(sdsize+1):
         if i<threshold:  
            self_kforce.setdefault(i,1.0)
         else:
            self_kforce.setdefault((i),1/float(np.power(sfactor,i-threshold)))
       for i in range(xdsize+1):
         if i<threshold:  
            x_kforce.setdefault(i,1.0)
         else:
            x_kforce.setdefault((i),1/float(np.power(xfactor,i-threshold)))
       #Get intra chain restraints
       general={}
       done=[]
       kforce=self_kforce
       dsize=sdsize
       for c in self._dna_chains:
        if not phos_pos.has_key(c): continue
        for i in range(len(phos_pos[c])):
          a= max(0,i - dsize)
          b= min(i + dsize,len(phos_pos[c]))
          pos_i = tuple((c,phos_pos[c][i]))
          gen_i = tuple((c,i-i))
          for j in range(a,b):
             if j==i:continue
             pos_j = tuple((c,phos_pos[c][j]))
             gen_j = tuple((c,j-i))
             if tuple((pos_i,pos_j)) in done:continue
             if phos_xyz.has_key(pos_i) and phos_xyz.has_key(pos_j):
                diff = phos_xyz[pos_i] - phos_xyz[pos_j]
                dist = np.sqrt(diff.dot(diff))
                k    = kforce[np.abs(j-i)]
                if dist>maxdist: continue
                #print("CHECK DISTANCE ",len(phos_pos[c]),a,b,c,i,j,pos_i,pos_j, phos_xyz[pos_i],phos_xyz[pos_j],diff)
                general.setdefault(tuple((gen_i,gen_j)),tuple((dist,k)))
                done.append(tuple((pos_i,pos_j)))
                done.append(tuple((pos_j,pos_i)))
       #Get inter chain restraints
       done=[]
       kforce=x_kforce
       dsize=xdsize
       for c in self._dna_chains:
         if not phos_pos.has_key(c): continue
         for d in self._dna_chains:
           if c==d:continue
           if not phos_pos.has_key(d): continue
           for i in range(len(phos_pos[c])):
             pos_i = tuple((c,phos_pos[c][i]))
             gen_i = tuple((c,i-i))
             ii = len(phos_pos[d])-i-1
             a= max(0,ii - dsize)
             b= min(ii + dsize,len(phos_pos[d]))
             for j in range(a,b):
                pos_j = tuple((d,phos_pos[d][j]))
                gen_j = tuple((d,j-ii))
                if tuple((pos_i,pos_j)) in done:continue
                if phos_xyz.has_key(pos_i) and phos_xyz.has_key(pos_j):
                   diff = phos_xyz[pos_i] - phos_xyz[pos_j]
                   dist = np.sqrt(diff.dot(diff))
                   k    = kforce[np.abs(j-ii)]
                   if dist>maxdist: continue
                   general.setdefault(tuple((gen_i,gen_j)),tuple((dist,k)))
                   done.append(tuple((pos_i,pos_j)))
                   done.append(tuple((pos_j,pos_i)))
       #Add DNA restraints
       for i in range(dna_size):
         for key in general:
            posi,posj = key
            dist,ratio= general[key]
            chain_1,res_1=posi
            chain_2,res_2=posj
            pi= i + int(res_1) + 1
            if chain_1 == chain_2:
                 pj = max(1, i + int(res_2) + 1 ) #res_2 = j-i => 
                 pj = min(pj, dna_size )
                 pos_i = tuple( (chain_1,str(pi)+' ') )
                 pos_j = tuple( (chain_2,str(pj)+' ') )
                 self._restraints.setdefault(tuple((pos_i,pos_j)),tuple((dist,ratio)))
            else:
                 ii = dna_size - i  - 1            # ii = len(phos_pos[d])-i-1 => (dna_size -i -1) position
                 pj = max(1, ii + int(res_2) + 1 )   # res_2 = j-ii
                 pj = min(pj, dna_size )
                 pos_i = tuple( (chain_1,str(pi)+' ') )
                 pos_j = tuple( (chain_2,str(pj)+' ') )
                 self._restraints.setdefault(tuple((pos_i,pos_j)),tuple((dist,ratio)))


   def get_dna_chains(self):
       return self._dna_chains

   def get_dna_bead_size(self):
       return self._dna_bead

   def get_dna_files(self):
       return self._dna_files

   def get_dna_sequences(self):
       return self._dna_sequences

   def get_dna_sdsize(self):
       return self._dna_sdsize 

   def get_dna_xdsize(self):
       return self._dna_xdsize 

   def get_dna_xfactor(self):
       return self._dna_xfactor 

   def get_dna_sfactor(self):
       return self._dna_sfactor 

   def get_dna_bead(self):
       return self._dna_bead

   def get_dna_maxdist(self):
       return self._dna_maxdist 

   def copy(self,name):

       other = self.__class__(name)                              
       other.copy_fragments(self.get_fragments())
       other.copy_dna_files(self.get_dna_files())
       other.copy_complexes(self.get_complexes())
       other.copy_protein_chains(self.get_protein_chains())
       other.copy_protein_sequences(self.get_protein_sequences())
       other.copy_dna_chains(self.get_dna_chains())   
       other.copy_dna_sequences(self.get_dna_sequences())
       other.copy_restraints(self.get_restraints())
       other.copy_interfaces(self.get_interfaces())
       other.copy_chain_correspondence_in_core(self.get_chain_correspondence_in_core())
       other.copy_chain_correspondence_plus(self.get_chain_correspondence_plus())
       other.copy_topology(self.get_topology())
       other.copy_protein_reference_site(self.get_protein_reference_site())
       other.set_dna_sdsize(         self.get_dna_sdsize())
       other.set_dna_xdsize(         self.get_dna_xdsize())
       other.set_dna_sfactor(        self.get_dna_sfactor())
       other.set_dna_xfactor(        self.get_dna_xfactor())
       other.set_dna_bead(           self.get_dna_bead())
       other.set_dna_maxdist(        self.get_dna_maxdist())
       other.set_protein_dsize(      self.get_protein_dsize())
       other.set_protein_sfactor(    self.get_protein_sfactor())
       other.set_protein_threshold(  self.get_protein_threshold())
       other.set_protein_maxdist(    self.get_protein_maxdist())
       other.set_interface_overlap(  self.get_interface_overlap())
       other.set_minimum_interface(  self.get_minimum_interface())
       other.set_minimum_addition(   self.get_minimum_addition())
       other.set_restrict_equal(     self.get_restrict_equal())
       return other

   def restraints_to_csv(self,interaction_type="all",kd=20.0,weight=1.0,sd=1.0):
       data={}
       for key in self._restraints:
         posi,posj = key
         dist,ratio= self._restraints[key]
         name_1,res_1=posi
         name_2,res_2=posj
         if name_1 in self._protein_chains: mtype_1="prot"
         if name_2 in self._protein_chains: mtype_2="prot"
         if name_1 in self._dna_chains: mtype_1="dna"
         if name_2 in self._dna_chains: mtype_2="dna"
         if interaction_type == mtype_1+mtype_2 or interaction_type == mtype_2+mtype_1 or interaction_type=="all":
            data.setdefault("prot1",[]).append(name_1)
            data.setdefault("prot2",[]).append(name_2)
            data.setdefault("prot1_res",[]).append(res_1)
            data.setdefault("prot2_res",[]).append(res_2)
            data.setdefault("distance",[]).append(float(dist))
            data.setdefault("kd",[]).append(float(ratio)*kd)
            data.setdefault("weight",[]).append(float(weight))
            data.setdefault("sd",[]).append(float(sd))
       table=pd.DataFrame(data)
       return table

   def write_restraints(self,output,interaction_type="all",kd=20.0,weight=1.0,sd=1.0):
      if os.path.exists(output) and os.path.isfile(output):
       fo=open(output,"a")
      else:
       fo=open(output,"w")
       fo.write("prot1,prot1_res,prot2,prot2_res,distance,kd,weight,sd\n")
      for key in self._restraints:
        posi,posj = key
        dist,ratio= self._restraints[key]
        name_1,res_1=posi
        name_2,res_2=posj
        if name_1 in self._protein_chains: mtype_1="prot"
        if name_2 in self._protein_chains: mtype_2="prot"
        if name_1 in self._dna_chains: mtype_1="dna"
        if name_2 in self._dna_chains: mtype_2="dna"
        if interaction_type == mtype_1+mtype_2 or interaction_type == mtype_2+mtype_1 or interaction_type=="all":
           fo.write("%s,%s,%s,%s,%10.3f,%10.3f,%10.3f,%10.3f\n"%(name_1,res_1,name_2,res_2,float(dist),float(ratio)*kd,weight,sd))
      fo.close()

   def set_topology_core(self):
      data={}
      red_color   = ["darkred", "firebrick","tomato","red", "orangered", "coral", "orange", "goldenrod","gold","siena","yellow"]
      blue_color  = ["navy","darkblue","blue","midnightblue","dodgerblue","steelblue","deepskyblue","skyblue","lightblue","cyan","lightcyan"]
      green_color = ["darkgreen","forestgreen","green","lawngreen","springgreen","lightgreen","limegreen","turquoise","aquamarine"]
      color = red_color + blue_color + green_color
      fasta_name = self._id+".topology.fasta"
      color_molecule = 0
      dna_bead_size  = max(2,self.get_dna_bead_size()) 
      for i in range(len(self._fragments)):
          if i >0: previous_fragment=self._fragments[i-1]
          else: previous_fragment=None
          fragment  = self._fragments[i]
          dna_file  = self._dna_files[i]
          pdb_file  = self.get_pdb_by_fragment(fragment)
          dna_file  = self._dna_files[i]
          dna_cif   = dna_file.rstrip("pdb")+"cif"
          dna_obj   = PDB(dna_file)
          pdb_obj   = PDB(pdb_file)
          chain_ids=[]
          dna_size_chain={}
          for chain in dna_obj.chains:
            if chain.chaintype=="N":
             chain_id = chain.chain
             chain_ids.append(chain_id)
             dna_size_chain.setdefault(chain_id,len(chain.nucleotide_sequence()))
          chain_ids.sort()
          forward  = chain_ids[0]
          reverse   = chain_ids[1]
          start,end = fragment
          fragment_size   = int(end)-int(start) +1 
          for d in range(len(self._dna_chains)):
              chain_code = self._dna_chains[d]
              c = chain_code.lstrip("DNA")
              if c==forward: n_offset = int(start) - 1
              if c==reverse: n_offset = len(self._dna_sequences[d])-int(end)
              #The chains have been already renumbered so offset = 0
              offset = 0
              n = 1
              if previous_fragment is not None and c==forward:
                  last_previous = int(previous_fragment[1])
                  while n+n_offset<=last_previous:
                      n=n+1
              for j in range(0,dna_size_chain[c],dna_bead_size-1):
                if previous_fragment is not None and c==reverse:
                      if n+n_offset+dna_bead_size-1>= len(self._dna_sequences[d])-last_previous:
                         last_residue = len(self._dna_sequences[d])-last_previous - n_offset 
                         if last_residue <=1: break
                         nn=n + n_offset
                         mm=last_residue + n_offset
                         if nn>=mm:break
                         data.setdefault("residue_range",[]).append(str(nn)+","+str(mm))
                         data.setdefault("molecule_name",[]).append(chain_code)
                         #data.setdefault("color",[]).append(color[color_molecule])
                         data.setdefault("color",[]).append("red")
                         data.setdefault("fasta_fn",[]).append(os.path.basename(fasta_name))
                         data.setdefault("fasta_id",[]).append(chain_code+","+"DNA")
                         data.setdefault("pdb_fn",[]).append(os.path.basename(dna_cif))
                         data.setdefault("chain",[]).append(c)
                         data.setdefault("pdb_offset",[]).append(offset)
                         data.setdefault("bead_size",[]).append(1)
                         data.setdefault("em_residues_per_gaussian",[]).append(0)
                         if data.has_key("rigid_body"): rb=len(data.get("rigid_body")) + 1
                         else: rb=1
                         data.setdefault("rigid_body",[]).append(rb)
                         data.setdefault("super_rigid_body",[]).append(" ")
                         data.setdefault("chain_of_super_rigid_bodies",[]).append(" ")
                         break
                if n >= dna_size_chain[c]: break
                nn=n + n_offset
                mm=min(n+dna_bead_size-1,dna_size_chain[c]) + n_offset
                if nn>=mm:break
                data.setdefault("residue_range",[]).append(str(nn)+","+str(mm))
                data.setdefault("molecule_name",[]).append(chain_code)
                #data.setdefault("color",[]).append(color[color_molecule])
                data.setdefault("color",[]).append("red")
                data.setdefault("fasta_fn",[]).append(os.path.basename(fasta_name))
                data.setdefault("fasta_id",[]).append(chain_code+","+"DNA")
                data.setdefault("pdb_fn",[]).append(os.path.basename(dna_cif))
                data.setdefault("chain",[]).append(c)
                data.setdefault("pdb_offset",[]).append(offset)
                data.setdefault("bead_size",[]).append(1)
                data.setdefault("em_residues_per_gaussian",[]).append(0)
                if data.has_key("rigid_body"): rb=len(data.get("rigid_body")) + 1
                else: rb=1
                data.setdefault("rigid_body",[]).append(rb)
                data.setdefault("super_rigid_body",[]).append(" ")
                data.setdefault("chain_of_super_rigid_bodies",[]).append(" ")
                n=n+dna_bead_size
      
      for protein in self._protein_chains:
          color_molecule = color_molecule + 1
          pdb_chain  = self.get_pdb_and_chain_in_core(protein)
          if pdb_chain is None: continue
          pdb_file,c = pdb_chain
          pdb_cif    = pdb_file.rstrip("pdb")+"cif"
          pdb_obj    = PDB(pdb_file)
          chain      = pdb_obj.get_chain_by_id(c)
          start      = chain.first_aminoacid.identifier
          end        = chain.last_aminoacid.identifier
          offset     = 0
          data.setdefault("molecule_name",[]).append(protein)
          #data.setdefault("color",[]).append(color[np.mod(color_molecule, len(color))])
          data.setdefault("color",[]).append("blue")
          data.setdefault("fasta_fn",[]).append(os.path.basename(fasta_name))
          data.setdefault("fasta_id",[]).append(protein)
          data.setdefault("pdb_fn",[]).append(os.path.basename(pdb_cif))
          data.setdefault("chain",[]).append(c)
          data.setdefault("residue_range",[]).append(str(start)+","+str(end))
          data.setdefault("pdb_offset",[]).append(offset)
          data.setdefault("bead_size",[]).append(1)
          data.setdefault("em_residues_per_gaussian",[]).append(0)
          if data.has_key("rigid_body"): rb=len(data.get("rigid_body")) + 1
          else: rb=1
          data.setdefault("rigid_body",[]).append(rb)
          data.setdefault("super_rigid_body",[]).append(" ")
          data.setdefault("chain_of_super_rigid_bodies",[]).append(" ")
      self._topology=data


#-------------#
# Functions   #
#-------------#

def is_dna(chain):
    residues = chain.get_residues()
    for residue in residues:
        if residue.get_resname().strip() in ["DA", "DT", "DG", "DC", "DU"]:
            return True
        else:
            return False

def merge_mmcif_header(header,cif1,cif2):
    if os.path.exists(cif2): os.remove(cif2)
    fo=open(cif2,"w")
    skip=False
    if os.path.exists(header):
       fi=open(header,"r")
       for line in fi:
           fo.write(line)
       fi.close()
       skip=True
    if os.path.exists(cif1):
       fi=open(cif1,"r")
       feature=[]
       for line in fi:
           if line.startswith("#"): 
              fo.write(line)
              skip = False
           if skip:continue
           data = line.split()
           if data[0]=="loop_":
              feature=[]
              fo.write(line)
              continue
           elif data[0].startswith("_atom_site"):
              fo.write(line)
              dummy,promp = data[0].split(".")
              feature.append(promp)
              continue
           else: 
             if len(data) == len(feature):
              atom_dict={}
              for i in range(len(feature)):
                  atom_dict.setdefault(feature[i],data[i])
              for f in feature:
                  if "asym_id" in f:
                     atom_dict[f] = atom_dict["auth_asym_id"]
                  if "seq_id" in f:
                     atom_dict[f] = atom_dict["auth_seq_id"]
              for f in feature:
                  fo.write("%s\t"%str(atom_dict[f]))
              fo.write("\n")
       fi.close()

def renamed_complexes_set(complexes):


    new_complexes = []
    names = {}
    for x in complexes:
        name = x.get_name()
        #x.set_restrict_equal(False)
        root = name.split(".")[0]
        names.setdefault(root,[]).append(x)
    for r_name,list_x in names.iteritems():
        text_of_cmpx=", ".join([cx.get_name() for cx in list_x])
        print("\t\tRename complexes of %s: %s "%(r_name,text_of_cmpx))
        size = len(list_x)
        done = []
        n    = size
        for i in range(size):
            new_name   = r_name+"."+str(i)
            tfcomplex  = list_x[i]
            skip = False
            for t2_complex in done:
              if skip: continue
              if tfcomplex == t2_complex:
                 print("\t\t--Redundancy: %s is already as the new named %s"%(tfcomplex.get_name(),t2_complex.get_name()))
                 skip = True
                 break
            if skip: continue
            for t2_complex in done:
                done_name = t2_complex.get_name()
                if new_name == done_name:
                   new_name = r_name+"."+str(n)
                   n = n +1
            print("\t\t--Rename complex %s as %s"%(tfcomplex.get_name(),new_name))
            tfcomplex.set_name(new_name)
            done.append(tfcomplex)
            new_complexes.append(tfcomplex)
    return new_complexes
 

 
def dna_reverse(nuc):
    dna_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    complement = "".join([dna_complement[n] for n in nuc])
    reverse = "".join([complement[len(complement) - i - 1] for i in range(len(complement))])
    return reverse

def write_pseudo_mmcif(pdb_obj,mmcif_file,header_only=False,force=False):
    """
    preliminar code to write coordiantes of a PDB object from SBI in mmCif format (version for py2 of SBI doesn't have this possibility)
    """

    if force==True and os.path.exists(mmcif_file): os.remove(mmcif_file)

    fo = open(mmcif_file,"w")
    fo.write("data_%s\n"%pdb_obj.id)

    #Write header relevant for Chimera

    fo.write("#\nloop_\n")
    fo.write("_entity_poly_seq.entity_id\n") 
    fo.write("_entity_poly_seq.num\n") 
    fo.write("_entity_poly_seq.mon_id\n") 
    fo.write("_entity_poly_seq.hetero\n") 

    entity_id = 0
    for chain in pdb_obj.chains:
        entity_id = entity_id +1
        num = 0
        if chain.chaintype=="N":
           for residue in chain.nucleotides:
               num = num +1
               fo.write("%s\t%s\t%s\t%s\n"%(str(entity_id),str(num), residue.type,"n"))           
        if chain.chaintype=="P":
           for residue in chain.aminoacids:
               num = num +1
               fo.write("%s\t%s\t%s\t%s\n"%(str(entity_id),str(num), residue.type,"n"))           
    if   header_only:
         fo.close()
         return

    #Write coordinates
  
    data = [
     "group_PDB"
     ,"id"
     ,"type_symbol"
     ,"label_atom_id" 
     ,"label_alt_id" 
     ,"label_comp_id" 
     ,"label_asym_id" 
     ,"label_entity_id" 
     ,"label_seq_id" 
     ,"pdbx_PDB_ins_code" 
     ,"Cartn_x" 
     ,"Cartn_y" 
     ,"Cartn_z" 
     ,"occupancy" 
     ,"B_iso_or_equiv" 
     ,"pdbx_formal_charge" 
     ,"auth_seq_id" 
     ,"auth_comp_id" 
     ,"auth_asym_id" 
     ,"auth_atom_id" 
     ,"pdbx_auth_seq_id" 
     ,"pdbx_auth_comp_id" 
     ,"pdbx_auth_asym_id" 
     ,"pdbx_auth_atom_name" 
     ,"pdbx_PDB_model_num" 
    ]
    data_short=[
     "group_PDB"
     ,"id"
     ,"type_symbol"
     ,"label_atom_id" 
     ,"label_alt_id" 
     ,"label_comp_id" 
     ,"label_asym_id" 
     ,"label_entity_id" 
     ,"label_seq_id" 
     ,"pdbx_PDB_ins_code" 
     ,"Cartn_x" 
     ,"Cartn_y" 
     ,"Cartn_z" 
     ,"occupancy" 
     ,"B_iso_or_equiv" 
     ,"auth_seq_id" 
     ,"auth_asym_id" 
     ,"auth_atom_id" 
     ,"pdbx_PDB_model_num" 
    ]
    n = 0
    m = 0
    pdbx_PDB_ins_code = "?"
    pdbx_PDB_model_num= 1
    group_PDB="ATOM"
    label_alt_id = "."
    fo.write("#\nloop_\n")
    for d in data_short:
        fo.write("_atom_site.%s\n"%d)
    for chain in pdb_obj.chains:
        m                 = m +1
        #label_entity_id   = m
        label_entity_id   = "?"
        chain_id          = chain.chain
        label_asym_id     = chain_id.upper()
        pdbx_auth_asym_id = auth_asym_id  = chain_id 
        label_seq_id      = 0
        if chain.chaintype=="N":
           for residue in chain.nucleotides:
               label_comp_id = residue.standard_type
               pdbx_auth_comp_id = auth_comp_id  = residue.type
               label_seq_id  = label_seq_id + 1
               pdbx_auth_seq_id = auth_seq_id   = residue.number
               for atom in residue.atoms:
                   n        = n +1
                   label_id = n
                   if atom.charge != "": pdbx_formal_charge = atom.charge
                   else: pdbx_formal_charge  = "?"
                   type_symbol = atom.element
                   label_atom_id = atom.name
                   pdbx_auth_atom_name = auth_atom_id = atom.name
                   occupancy = atom.occupancy
                   B_iso_or_equiv =atom.tempFactor
                   Cartn_x = atom.x
                   Cartn_y = atom.y
                   Cartn_z = atom.z
                   atom_dict = {}
                   atom_dict.setdefault("group_PDB",group_PDB)
                   atom_dict.setdefault("id",label_id)
                   atom_dict.setdefault("type_symbol",type_symbol)
                   atom_dict.setdefault("label_atom_id",label_atom_id)
                   atom_dict.setdefault("label_alt_id",label_alt_id)
                   atom_dict.setdefault("label_comp_id",label_comp_id)
                   atom_dict.setdefault("label_asym_id",label_asym_id)
                   atom_dict.setdefault("label_entity_id",label_entity_id)
                   atom_dict.setdefault("label_seq_id",label_seq_id)
                   atom_dict.setdefault("pdbx_PDB_ins_code",pdbx_PDB_ins_code)
                   atom_dict.setdefault("Cartn_x","%.3f"%Cartn_x)
                   atom_dict.setdefault("Cartn_y","%.3f"%Cartn_y)
                   atom_dict.setdefault("Cartn_z","%.3f"%Cartn_z)
                   atom_dict.setdefault("occupancy",occupancy)
                   atom_dict.setdefault("B_iso_or_equiv",B_iso_or_equiv)
                   atom_dict.setdefault("pdbx_formal_charge",pdbx_formal_charge)
                   atom_dict.setdefault("auth_seq_id",auth_seq_id)
                   atom_dict.setdefault("auth_comp_id",auth_comp_id)
                   atom_dict.setdefault("auth_atom_id",auth_atom_id)
                   atom_dict.setdefault("auth_asym_id",auth_asym_id)
                   atom_dict.setdefault("pdbx_auth_seq_id",pdbx_auth_seq_id)
                   atom_dict.setdefault("pdbx_auth_comp_id",pdbx_auth_comp_id)
                   atom_dict.setdefault("pdbx_auth_asym_id",pdbx_auth_asym_id)
                   atom_dict.setdefault("pdbx_auth_atom_name",pdbx_auth_atom_name)
                   atom_dict.setdefault("pdbx_PDB_model_num",pdbx_PDB_model_num)
                   for d in data_short:
                      x=str(atom_dict[d])
                      if "'" in x: fo.write("'%s'\t"%x)
                      else: fo.write("%s\t"%x)
                   fo.write("\n")
        if chain.chaintype=="P":
           for residue in chain.aminoacids:
               label_comp_id = residue.standard_type
               pdbx_auth_comp_id = auth_comp_id  = residue.type
               label_seq_id  = label_seq_id + 1
               pdbx_auth_seq_id = auth_seq_id   = residue.number
               for atom in residue.atoms:
                   n        = n +1
                   label_id = n
                   if atom.charge != "": pdbx_formal_charge = atom.charge
                   else: pdbx_formal_charge  = "?"
                   label_atom_id = atom.name
                   type_symbol = atom.element
                   pdbx_auth_atom_name = auth_atom_id = atom.name
                   B_iso_or_equiv =atom.tempFactor
                   occupancy = atom.occupancy
                   Cartn_x = atom.x
                   Cartn_y = atom.y
                   Cartn_z = atom.z
                   atom_dict = {}
                   atom_dict.setdefault("group_PDB",group_PDB)
                   atom_dict.setdefault("id",label_id)
                   atom_dict.setdefault("type_symbol",type_symbol)
                   atom_dict.setdefault("label_atom_id",label_atom_id)
                   atom_dict.setdefault("label_alt_id",label_alt_id)
                   atom_dict.setdefault("label_comp_id",label_comp_id)
                   atom_dict.setdefault("label_asym_id",label_asym_id)
                   atom_dict.setdefault("label_entity_id",label_entity_id)
                   atom_dict.setdefault("label_seq_id",label_seq_id)
                   atom_dict.setdefault("pdbx_PDB_ins_code",pdbx_PDB_ins_code)
                   atom_dict.setdefault("Cartn_x","%.3f"%Cartn_x)
                   atom_dict.setdefault("Cartn_y","%.3f"%Cartn_y)
                   atom_dict.setdefault("Cartn_z","%.3f"%Cartn_z)
                   atom_dict.setdefault("occupancy",occupancy)
                   atom_dict.setdefault("B_iso_or_equiv",B_iso_or_equiv)
                   atom_dict.setdefault("pdbx_formal_charge",pdbx_formal_charge)
                   atom_dict.setdefault("auth_seq_id",auth_seq_id)
                   atom_dict.setdefault("auth_comp_id",auth_comp_id)
                   atom_dict.setdefault("auth_atom_id",auth_atom_id)
                   atom_dict.setdefault("auth_asym_id",auth_asym_id)
                   atom_dict.setdefault("pdbx_auth_seq_id",pdbx_auth_seq_id)
                   atom_dict.setdefault("pdbx_auth_comp_id",pdbx_auth_comp_id)
                   atom_dict.setdefault("pdbx_auth_asym_id",pdbx_auth_asym_id)
                   atom_dict.setdefault("pdbx_auth_atom_name",pdbx_auth_atom_name)
                   atom_dict.setdefault("pdbx_PDB_model_num",pdbx_PDB_model_num)
                   for d in data_short:
                      x=str(atom_dict[d])
                      if "'" in x: fo.write("'%s'\t"%x)
                      else: fo.write("%s\t"%x)
                   fo.write("\n")
    fo.write("#\n")
    fo.close()


def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python tf_complex.py -i input_folder  -d data_folder [ -o output_name ]")

    parser.add_option("-i", action="store", type="string", dest="input_folder", default=None, help="Input folder, contains folders 'pair_interactions' and 'TF_DNA_FRAGMENTS' as from ModCRE", metavar="{directory}")
    parser.add_option("-d", action="store", type="string", dest="data_folder", default=None, help="Output folder with restraints and topology to run IMP modelling", metavar="{directory}")
    parser.add_option("--dummy", action="store", type="string", dest="dummy_folder", default="/tmp", help="Dummy folder (default is /tmp )", metavar="{directory}")
    parser.add_option("-o", action="store", type="string", dest="output_file", default=None, help="Output rootname (default is the basename of input folder", metavar="{filename}")
    parser.add_option("--all_conformers", default=False, action="store_true", dest="all_conformers", help="Check all complexes with different PPI models (default is false) WARNING: it requieres large computation", metavar="{boolean}")
    parser.add_option("--max_branches", default=None, action="store", dest="max_branches", help="Limit the number of branches of different complexes with different restrains (default is None) WARNING: unsetting a limit may requiere large memory", metavar="{integer}")
    parser.add_option("--max_iterations", default=None, action="store", dest="max_iterations", help="Limit the number of iterations (default is None) WARNING: unsetting a limit may requiere large memory and computation", metavar="{integer}")
    parser.add_option("--compact", default=False, action="store_true", dest="compact", help="Flag to compact the newest complex with all internal PPIs WARNING:it may force to over packing the complex", metavar="{boolean}")
    parser.add_option("--PPI_distance", default=12.0, action="store", dest="PPI_distance", help="Distance between CB atoms to calculate the interface of PPIs", metavar="{float}")    
    parser.add_option("--trim_restraints", default=1.0, action="store", dest="trim_restraints", help="Ratio of acceptance to trim the total of PPI restraints (defult no trimming)", metavar="{float}")    

    (options, args) = parser.parse_args()

    if options.input_folder is None or options.data_folder is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Parameters
    output_folder = options.input_folder
    data_folder   = options.data_folder 
    dummy         = options.dummy_folder 
    all_conformers= options.all_conformers
    max_branches  = options.max_branches
    max_iterations= options.max_iterations
    compact       = options.compact
    ratio_trimming= float(options.trim_restraints)
    PPI_distance  = float(options.PPI_distance)

    if not dummy.startswith("/"): dummy=os.path.abspath(dummy)
    if not os.path.exists(dummy): os.makedirs(dummy)
    output_file   = os.path.basename(options.output_file)
    if output_file is None: output_file   = os.path.basename(output_folder)

    kd_dna = 2.0
    kd_PNI = 3.0
    kd_PPI = 1.0
    w_dna  = 1.0
    w_PNI  = 3.0
    w_PPI  = 2.0
    sd_dna = 2.0
    sd_PNI = 1.0
    sd_PPI = 1.0
    output_DNA_fragments = os.path.join(output_folder,"TF_DNA_FRAGMENTS")
    output_PPI = os.path.join(output_folder,"pair_interactions")
    job_dir    = os.path.join(output_folder,"../..") 
    #Save time using protein_self_restraints=False
    protein_self_restraints=False
    #protein_self_restraints=True


    #Clean output directory
    if not os.path.exists(data_folder): 
        os.makedirs(data_folder)
    else:
        for tfile in os.listdir(data_folder):
            if "topology" in tfile:
                os.remove(os.path.join(data_folder,tfile))
            if "restraints" in tfile:
                os.remove(os.path.join(data_folder,tfile))

    #Make list of complexes (combinations of different fragments)
    new_complex_list=[]
    complex_codes={}
    fragments_a=[]
    fragments_b=[]
    for output_fragment in os.listdir(output_DNA_fragments):
       if output_fragment.startswith("fragment") and os.path.isdir(os.path.join(output_DNA_fragments,output_fragment)):
          complex_list = new_complex_list
          print("\t--Fragment "+output_fragment)
          fragments_a.append(int(output_fragment.split("_")[-1].split("-")[0]))
          fragments_b.append(int(output_fragment.split("_")[-1].split("-")[1]))
          output_dir = os.path.join(output_DNA_fragments,output_fragment)
          new_complex_list = []
          n=0
          # Compose the combination of structures in each fragment, each structure in a fragment being numbered from 1 to M, being M the total of structures 
          for pdb_file in os.listdir(os.path.join(output_DNA_fragments,output_fragment)):
            #print("\t--check file "+pdb_file)
            if pdb_file.endswith(".pdb") and pdb_file.startswith("dna"):
                if pdb_file.endswith(".non_optimized.pdb"):continue
                print("\t-- Add file "+pdb_file)
                complex_codes.setdefault(output_fragment,[]).append(pdb_file.lstrip("dna__").rstrip("pdb"))
                n=n+1
                if len(complex_list)>0:
                  for cmpx in complex_list:
                      new_complex_list.append(cmpx+"_"+output_fragment+":"+str(n))
                else:
                    new_complex_list.append("_"+output_fragment+":"+str(n))
    complex_list = new_complex_list
    fragments_a.sort()
    fragments_b.sort()
    fragments=[x for x in zip(fragments_a,fragments_b)]

    #Make a set for all PPIs
    set_of_ppis=set()
    if os.path.exists(output_PPI):
       for pdb_file in os.listdir(output_PPI):
        if pdb_file.endswith(".pdb") and "::" in pdb_file:
            if len(pdb_file.split("::"))>1:
               try:
                 start_1,end_1 = pdb_file.split("::")[1].split("_")[-1].split("-")
                 chain_id      = pdb_file.split("::")[1].split("_")[-2]
                 pdb_obj       = PDB(os.path.join(output_PPI,pdb_file))
                 chain         = pdb_obj.get_chain_by_id("A")
                 start_a       = int(chain.first_aminoacid.identifier)
                 if start_a != int(start_1): 
                    print("\t-- Skip "+pdb_file+" missmatch chain "+chain_id+": "+str(start_a)+" with "+str(start_1))
                    continue
                 start_2,end_2 = pdb_file.split("::")[2].split("_")[-1].split("-")
                 chain_id2     = pdb_file.split("::")[2].split("_")[-2]
                 chain         = pdb_obj.get_chain_by_id("B")
                 start_b       = int(chain.first_aminoacid.identifier)
                 if start_b != int(start_2): 
                    print("\t-- Skip "+pdb_file+" missmatch chain "+chain_id+": "+str(start_b)+" with "+str(start_2))
                    continue
                 print("\t-- Accept "+pdb_file)
                 set_of_ppis.add(os.path.join(output_PPI,pdb_file))
               except:
                 print("\t-- Read error "+pdb_file)
                 chain_id      = pdb_file.split("::")[1].split("_")[-2]
                 chain_id2     = pdb_file.split("::")[2].split("_")[-2]
                 print("\t-- Skip "+pdb_file+" chains "+chain_id+" "+chain_id2)
                 continue

    #Make data file of sequences 
    sequences_file=os.path.join(output_DNA_fragments,"sequences.fa")
    ppi_file=os.path.join(output_folder,"ppi.fa")
    fasta_file=os.path.join(output_DNA_fragments,"all_sequences.fa")
    if os.path.exists(fasta_file): os.remove(fasta_file)
    sequences={}
    if os.path.exists(sequences_file):
      for protein,sequence in parse_fasta_file(sequences_file,gz=False,clean=True):
          try:
             m = re.search("^(tr|sp|TR|SP)\|(\S+)\|\S+\_(\S+)",protein)
             if m:
                 code=m.group(2)
             else:
                 if len(protein.split("|"))>0: code=protein.split("|")[1]
                 else: code=protein.split("|")[0]
          except:
             code=protein
          sequences.setdefault(code,sequence)
    sequences_file = os.path.join(job_dir,"pwm_database", "sequences","sequences.fasta")
    if os.path.exists(sequences_file) and os.path.isfile(sequences_file):
      for protein,sequence in parse_fasta_file(sequences_file,gz=False,clean=True):
          try:
             m = re.search("^(tr|sp|TR|SP)\|(\S+)\|\S+\_(\S+)",protein)
             if m:
                 code=m.group(2)
             else:
                 if len(protein.split("|"))>0: code=protein.split("|")[1]
                 else: code=protein.split("|")[0]
          except:
             code=protein
          sequences.setdefault(code,sequence)
    if os.path.exists(ppi_file):
      for protein,sequence in parse_fasta_file(ppi_file,gz=False,clean=True):
          try:
             m = re.search("^(tr|sp|TR|SP)\|(\S+)\|\S+\_(\S+)",protein)
             if m:
                 code=m.group(2)
             else:
                 if len(protein.split("|"))>0: code=protein.split("|")[1]
                 else: code=protein.split("|")[0]
          except:
             code=protein
          sequences.setdefault(code,sequence)
    if os.path.exists(os.path.join(job_dir,"DNA")):
      for ff in os.listdir(os.path.join(job_dir,"DNA")):
        if not ff.endswith(".fa"): continue
        for dna_code,sequence in parse_fasta_file(os.path.join(job_dir,"DNA",ff),gz=False,clean=True):
            sequences.setdefault("DNA_fwd",sequence)
            sequences.setdefault("DNA_rev",dna_reverse(sequence))
    fasta = open(os.path.join(output_DNA_fragments,"all_sequences.fa"),"w")
    code_done=set()
    for code,sequence in sequences.iteritems():
        if code in code_done: continue
        if len(sequence)<1: continue
        if len(code)<1: code="DUMMY"
        fasta.write(">%s\n%s\n"%(code,sequence))
        code_done.add(code)
    fasta.close()

    #Make data file of dna_files and modify DNA sequence numbers
    first_dna= min([f[0] for f in fragments])
    last_dna = max([f[1] for f in fragments])
    dna_files=[]
    for i in range(len(fragments)):
        fragment_a,fragment_b = fragments[i]
        fragment = str(fragment_a)+"-"+str(fragment_b)
        output_fragment = "fragment_"+fragment
        out_dna_file = os.path.join(output_DNA_fragments,"dna_"+fragment+".pdb")
        name_pdb = "dna_"+fragment
        if os.path.exists(out_dna_file): os.remove(out_dna_file)
        if not os.path.exists(out_dna_file) :
          print("Make DNA file "+out_dna_file)
          for pdb_file in os.listdir(os.path.join(output_DNA_fragments,output_fragment)):
            #if os.path.exists(out_dna_file): continue
            if pdb_file.endswith("non_optimized.pdb"): continue
            if pdb_file.endswith(".pdb") and pdb_file.startswith("dna"):
               chain_ids=[]
               prot_ids =[]
               print("\t--Read PDB file "+pdb_file)
               pdb_obj  = PDB(os.path.join(output_DNA_fragments,output_fragment,pdb_file))
               for chain in pdb_obj.chains:
                if chain.chaintype=="N":
                   chain_id = chain.chain
                   chain_ids.append(chain_id)
                else:
                   chain_id = chain.chain
                   prot_ids.append(chain_id)
               chain_ids.sort()
               forward  = chain_ids[0]
               reverse  = chain_ids[1]
               new_pdb_obj=PDB()
               offset = fragment_a - 1
               #Renumber PDB file
               chain = pdb_obj.get_chain_by_id(forward)
               new_chain  = ChainOfNucleotide(name_pdb,forward)
               for residue in chain.nucleotides:
                   number  = int(residue.number)
                   version = residue.version
                   residue.number=number+offset
                   residue.version=version
                   new_chain.add_residue(residue)
               new_chain.chaintype=chain.chaintype
               print("\t\t--Add chain "+forward)
               new_pdb_obj.add_chain(new_chain)
               #Add reverse chain of DNA
               offset = last_dna - fragment_b
               #Renumber PDB file
               chain = pdb_obj.get_chain_by_id(reverse)
               new_chain  = ChainOfNucleotide(name_pdb,reverse)
               for residue in chain.nucleotides:
                   number  = int(residue.number)
                   version = residue.version
                   residue.number=number+offset
                   residue.version=version
                   new_chain.add_residue(residue)
               new_chain.chaintype=chain.chaintype
               print("\t\t--Add chain "+reverse)
               new_pdb_obj.add_chain(new_chain)
               # Write DNA file if it doesn't exists yet and store it
               if not os.path.exists(out_dna_file):
                  print("\t\t--Write file " + os.path.basename(out_dna_file))
                  new_pdb_obj.write(out_dna_file)
                  dna_files.append(out_dna_file)
               # Continue with the rest of chains
               for chain_id in prot_ids:
                  chain = pdb_obj.get_chain_by_id(chain_id)
                  new_pdb_obj.add_chain(chain)
               # Rewrite the PDB
               print("\t\t--Write file " + os.path.basename(pdb_file))
               new_pdb_obj.write(os.path.join(data_folder,pdb_file),force=True)



    # DNA CIF FILES ROTATED AND TRANSLATED
    dna_cif_files = []
    for i in range(len(fragments)):
           fragment_a,fragment_b = fragments[i]
           dna_file = dna_files[i]
           dna_obj  = PDB(dna_file)
           name     = os.path.basename(dna_file).rstrip("pdb")
           dna_cif  = os.path.join(data_folder,name+"cif")
           dummy_cif= os.path.join(dummy,name+"cif")
           header   = os.path.join(dummy,name+"header.cif")
           pdbparser= PDBParser()
           struct   = pdbparser.get_structure(name,dna_file)
           if i==0:
              iocif    = MMCIFIO()
              iocif.set_structure(struct)
              iocif.save(dummy_cif)
              write_pseudo_mmcif(dna_obj,header,header_only=True,force=True)
              #write dna_cif
              merge_mmcif_header(header,dummy_cif,dna_cif)
              dna_cif_files.append(dna_cif)
              continue
           template_fragment_a,template_fragment_b = fragments[i-1]
           template_file = dna_cif_files[-1]
           template_name = os.path.basename(template_file).rstrip("cif").rstrip(".")
           cifparser     = MMCIFParser()
           template      = cifparser.get_structure(template_name,template_file)
           chain_ids=[]
           for chain in template.get_chains():
               if is_dna(chain):
                  chain_ids.append(chain.get_id())
           chain_ids.sort()
           forward  = chain_ids[0]
           reverse  = chain_ids[1]
           overlap = max(0,1+int(template_fragment_b -fragment_a))
           if overlap == 0 : overlap=3
           #Add forward residues
           struct_list_fwd=[]
           for atom in struct.get_atoms():
               if atom.get_parent().get_parent().get_id() == forward and atom.element == "P":
                  struct_list_fwd.append(atom)
           template_list_fwd=[]
           for atom in template.get_atoms():
               if atom.get_parent().get_parent().get_id() == forward and atom.element == "P":
                  template_list_fwd.append(atom)
           template_atoms = template_list_fwd[-overlap:]
           moving_atoms   = struct_list_fwd[:overlap] 
           #Add reverse residues
           struct_list_rev=[]
           for atom in struct.get_atoms():
               if atom.get_parent().get_parent().get_id() == reverse and atom.element == "P":
                  struct_list_rev.append(atom)
           template_list_rev=[]
           for atom in template.get_atoms():
               if atom.get_parent().get_parent().get_id() == reverse and atom.element == "P":
                  template_list_rev.append(atom)
           template_atoms.extend(template_list_rev[:overlap])
           moving_atoms.extend(struct_list_rev[-overlap:])
           sup = Superimposer()
           sup.set_atoms(template_atoms,moving_atoms)
           sup.apply(struct.get_atoms())
           iocif    = MMCIFIO()
           iocif.set_structure(struct)
           iocif.save(dummy_cif)
           write_pseudo_mmcif(dna_obj,header,header_only=True,force=True)
           #write dna_cif
           merge_mmcif_header(header,dummy_cif,dna_cif)
           dna_cif_files.append(dna_cif)

    #Make TF_COMPLEX list
    n=0
    tf_complex_list=[]
    done=set()
    for complex_id in complex_list:
      n=n+1
      name=output_file+str(n)
      data=complex_id.split("_fragment_")
      #Make complex_files dictionary of core pdb_files
      complex_files={}
      print("\nREGULATORY COMPLEX WITH PDBs COMBINATION "+str(data))
      for tf_comb in data:
        if ":" not in tf_comb:continue
        #decompose the combination of fragments
        fragment_code,cmp_num  = tf_comb.split(":")
        fragment_folder = "fragment_"+fragment_code
        fragment = tuple( ([int(x) for x in fragment_code.split("-") ]))
        if complex_codes.has_key(fragment_folder):
          idx=int(cmp_num)-1
          pdb_name    = "dna__"+complex_codes[fragment_folder][idx]+"pdb"
          pdb_file    = os.path.join(data_folder,pdb_name)
          remark_name = "remarks_"+complex_codes[fragment_folder][idx].lstrip(fragment_code)+"txt"
          remark_file = os.path.join(output_DNA_fragments,fragment_folder,remark_name)
          data_remark_name = "remarks_"+complex_codes[fragment_folder][idx]+"txt"
          data_remark_file = os.path.join(data_folder,data_remark_name)
          if not os.path.exists(remark_file):
             remark_name = "remarks"+complex_codes[fragment_folder][idx].lstrip(fragment_code)+"txt"
             remark_file = os.path.join(output_DNA_fragments,fragment_folder,remark_name)
             data_remark_name = "remarks"+complex_codes[fragment_folder][idx]+"txt"
             data_remark_file = os.path.join(data_folder,data_remark_name)
          if not os.path.exists(pdb_file): continue
          if not os.path.exists(remark_file): continue
          print("\t-- PDB file :   "+pdb_file)
          print("\t-- CHAINS from: "+remark_file)
          if not os.path.exists(data_remark_file):
             try:
                shutil.copy(remark_file,data_remark_file)
                complex_files.setdefault(fragment,tuple( (pdb_file,data_remark_file) ))
             except:
                print("\t-- FAILED TO COPY %s TO %s"%(remark_file,data_remark_file))
                complex_files.setdefault(fragment,tuple( (pdb_file,remark_file) ))
          else:
             complex_files.setdefault(fragment,tuple( (pdb_file,data_remark_file) ))
      tf_complex=TF_complex(name,
                               fragments=fragments,
                               dna_files=dna_files,
                               complex_files=complex_files,
                               PPI_distance = PPI_distance
                               )
      tf_complex.initialize(fasta_file)
      if protein_self_restraints: tf_complex.set_all_inner_restraints_in_core()
      skip=False
      for t2_complex in done:
              if skip: continue
              if tf_complex == t2_complex: 
                 print("\t--Redundant complex "+tf_complex.get_name()+" is already "+t2_complex.get_name())
                 skip=True
      if skip: continue
      done.add(tf_complex)
      tf_complex_list.append(tf_complex)

    print("\nSUPERIMPOSE FRAGMENTS")       
    # COMPLEX CIF FRAGMENT FILES ROTATED AND TRANSLATED
    for tf_complex in tf_complex_list:
        complex_files=tf_complex.get_complexes()
        print("\t--Model "+tf_complex.get_name())
        for fragment in complex_files.keys():
            for i in range(len(fragments)):
               if fragments[i] != fragment: continue
               pdb_file, remark_file = complex_files[fragment]
               pdb_obj  = PDB(pdb_file)
               name     = os.path.basename(pdb_file).rstrip("pdb")
               dummy_cif= os.path.join(dummy,name+"cif")
               header   = os.path.join(dummy,name+"header.cif")
               pdb_cif  = os.path.join(data_folder,name+"cif")
               pdbparser= PDBParser()
               struct   = pdbparser.get_structure(name,pdb_file)
               if i==0:
                  iocif    = MMCIFIO()
                  iocif.set_structure(struct)
                  iocif.save(dummy_cif)
                  write_pseudo_mmcif(pdb_obj,header,header_only=True,force=True)
                  #write pdb_cif
                  merge_mmcif_header(header,dummy_cif,pdb_cif)
                  continue
               template_fragment_a,template_fragment_b = fragments[i-1]
               template_file = dna_cif_files[i-1]
               template_name = os.path.basename(template_file).rstrip("cif").rstrip(".")
               print("\t\t--Superimpose "+name + " on "+template_name)
               cifparser     = MMCIFParser()
               template      = cifparser.get_structure(template_name,template_file)
               chain_ids=[]
               for chain in struct.get_chains():
                   if is_dna(chain):
                      chain_ids.append(chain.get_id())
               chain_ids.sort()
               forward  = chain_ids[0]
               reverse  = chain_ids[1]
               overlap = max(0,1+int(template_fragment_b -fragment_a))
               if overlap == 0 : overlap=3
               #Add forward residues
               struct_list_fwd=[]
               for atom in struct.get_atoms():
                   if atom.get_parent().get_parent().get_id() == forward and atom.element == "P":
                      struct_list_fwd.append(atom)
               template_list_fwd=[]
               for atom in template.get_atoms():
                   if atom.get_parent().get_parent().get_id() == forward and atom.element == "P":
                      template_list_fwd.append(atom)
               template_atoms = template_list_fwd[-overlap:]
               moving_atoms   = struct_list_fwd[:overlap] 
               #Add reverse residues
               struct_list_rev=[]
               for atom in struct.get_atoms():
                   if atom.get_parent().get_parent().get_id() == reverse and atom.element == "P":
                      struct_list_rev.append(atom)
               template_list_rev=[]
               for atom in template.get_atoms():
                   if atom.get_parent().get_parent().get_id() == reverse and atom.element == "P":
                      template_list_rev.append(atom)
               template_atoms.extend(template_list_rev[:overlap])
               moving_atoms.extend(struct_list_rev[-overlap:])
               sup = Superimposer()
               sup.set_atoms(template_atoms,moving_atoms)
               sup.apply(struct.get_atoms())
               iocif    = MMCIFIO()
               iocif.set_structure(struct)
               iocif.save(dummy_cif)
               write_pseudo_mmcif(pdb_obj,header,header_only=True,force=True)
               #write pdb_cif
               merge_mmcif_header(header,dummy_cif,pdb_cif)


    # ITERATE UNTIL NO MORE NEW COMPLEXES ARE BUILT
    print("\n=====================================\nADD PROTEIN-PROTEIN INTERACTIONS\n=====================================")
    new_tf_complex_list=[]
    start_iteration = True
    end_iteration   = False
    number_of_iterations = 1
    while not end_iteration:
        if len(new_tf_complex_list) == len(tf_complex_list): 
            end_iteration = True
            break
        if max_iterations is not None:
            if number_of_iterations > int(max_iterations):
               end_iteration = True
               break
        # Check start
        if start_iteration:
           start_iteration=False
           #The new complex list contains the riginal complexes
           new_tf_complex_list = tf_complex_list[:]
        else:
           #Restart the original list of complexes in the next iterations
           tf_complex_list = new_tf_complex_list[:]
        
        # ADD NEW PROTEINS (NOT YET IN THE COMPLEX)  AND INTERACTIONS  WITH THE COMPLEX
        print("--ADD NEW PROTEINS WITH PROTEIN-PROTEIN INTERACTIONS WITH THE COMPLEX")
        new_tf_complex_set = renamed_complexes_set(new_tf_complex_list)
        new_tf_complex_list_start = new_tf_complex_set[:]
        done=set()
        for pdb_ppi_file in set_of_ppis:
            pdb_ppi=os.path.basename(pdb_ppi_file)
            print("\t--PPI file "+pdb_ppi)
            protein_1  = pdb_ppi.split("::")[0]
            protein_2  = pdb_ppi.split("::")[1].split("_")[0]
            if not all_conformers and tuple([protein_1,protein_2]) in done: continue
            sequence_1 = sequences[protein_1]
            sequence_2 = sequences[protein_2]
            #Get only the new complexes
            branch_list=[]
            for tf_complex in new_tf_complex_list_start:
                print("\t-- Complex "+tf_complex.get_name())
                new_complexes, new_restraints = tf_complex.add_ppi(os.path.join(output_PPI,pdb_ppi),protein_1,protein_2,sequence_1,sequence_2,condition="new",maximum=max_branches)
                if len(new_complexes) > 0 or new_restraints>0:
                   branch_list.extend([x for x in new_complexes])
                   done.add(tuple([protein_1,protein_2]))
                   done.add(tuple([protein_2,protein_1]))
            if len(branch_list)>0:
              #Transform the full list into unique new complexes and extend the list of newest complexes
              new_tf_complex_list.extend(branch_list)
              #Shrink the list to unique complexes and reorder them in a set
              new_tf_complex_set=renamed_complexes_set(new_tf_complex_list)
              #Get the new list of complexes and continiue with next PPIs
              new_tf_complex_list = new_tf_complex_set[:]
            for tf_complex in new_tf_complex_set:
                topology = tf_complex.get_topology()
                for pdb_topo in topology["pdb_fn"]:
                    if pdb_ppi == pdb_topo:
                       print("\t\t--Check topology: add PDB file "+pdb_topo)
                       shutil.copy(pdb_ppi_file,os.path.join(data_folder,pdb_topo))
        # ADD INTERACTIONS WITHIN EACH COMPLEX
        print("--ADD INTERACTIONS WITHIN THE COMPLEX")
        done=set()
        new_tf_complex_list_start = new_tf_complex_list[:]
        for pdb_ppi_file in set_of_ppis:
            pdb_ppi=os.path.basename(pdb_ppi_file)
            print("\t--PPI file "+pdb_ppi)
            protein_1  = pdb_ppi.split("::")[0]
            protein_2  = pdb_ppi.split("::")[1].split("_")[0]
            if not all_conformers and tuple([protein_1,protein_2]) in done: continue
            sequence_1 = sequences[protein_1]
            sequence_2 = sequences[protein_2]
            branch_list=[]
            for tf_complex in new_tf_complex_list_start:
                print("\t-- Complex "+tf_complex.get_name())
                new_complexes, new_restraints = tf_complex.add_ppi(os.path.join(output_PPI,pdb_ppi),protein_1,protein_2,sequence_1,sequence_2,condition="old",maximum=max_branches)
                if len(new_complexes) > 0 or new_restraints>0:
                   branch_list.extend([x for x in new_complexes])
                   done.add(tuple([protein_1,protein_2]))
                   done.add(tuple([protein_2,protein_1]))
            if len(branch_list)>0:
              #Transform the full list into unique new complexes and extend the list of newest complexes
              new_tf_complex_list.extend(branch_list)
              #Shrink the list to unique complexes and reorder them in a set
              new_tf_complex_set=renamed_complexes_set(new_tf_complex_list)
              #Get the new list of complexes and continiue with next PPIs
              new_tf_complex_list = new_tf_complex_set[:]
            for tf_complex in new_tf_complex_set:
                topology = tf_complex.get_topology()
                for pdb_topo in topology["pdb_fn"]:
                    if pdb_ppi == pdb_topo:
                       print("\t\t--Check topology: add PDB file "+pdb_topo)
                       shutil.copy(pdb_ppi_file,os.path.join(data_folder,pdb_topo))
        # ADD NEW PROTEINS (ANY)  AND INTERACTIONS WITH THE COMPLEX
        print("--ADD ALL PROTEINS WITH ANY PROTEIN-PROTEIN INTERACTION IN THE COMPLEX")
        done=set()
        new_tf_complex_list_start = new_tf_complex_list[:]
        for pdb_ppi_file in set_of_ppis:
            pdb_ppi=os.path.basename(pdb_ppi_file)
            print("\t--PPI file "+pdb_ppi)
            protein_1  = pdb_ppi.split("::")[0]
            protein_2  = pdb_ppi.split("::")[1].split("_")[0]
            if not all_conformers and tuple([protein_1,protein_2]) in done: continue
            sequence_1 = sequences[protein_1]
            sequence_2 = sequences[protein_2]
            branch_list=[]
            for tf_complex in new_tf_complex_list_start:
                print("\t-- Complex "+tf_complex.get_name())
                new_complexes, new_restraints = tf_complex.add_ppi(os.path.join(output_PPI,pdb_ppi),protein_1,protein_2,sequence_1,sequence_2,condition=None,maximum=max_branches)
                if len(new_complexes) > 0 or new_restraints>0:
                   branch_list.extend([x for x in new_complexes])
                   done.add(tuple([protein_1,protein_2]))
                   done.add(tuple([protein_2,protein_1]))
            if len(branch_list)>0:
              #Transform the full list into unique new complexes and extend the list of newest complexes
              new_tf_complex_list.extend(branch_list)
              #Shrink the list to unique complexes and reorder them in a set
              new_tf_complex_set=renamed_complexes_set(new_tf_complex_list)
              #Get the new list of complexes and continiue with next PPIs
              new_tf_complex_list = new_tf_complex_set[:]
            for tf_complex in new_tf_complex_set:
                topology = tf_complex.get_topology()
                for pdb_topo in topology["pdb_fn"]:
                    if pdb_ppi == pdb_topo:
                       print("\t\t--Check topology: add PDB file "+pdb_topo)
                       shutil.copy(pdb_ppi_file,os.path.join(data_folder,pdb_topo))
        # ADD INTERACTIONS WITHIN EACH COMPLEX
        done=set()
        new_tf_complex_list_start = new_tf_complex_list[:]
        if compact:
          print("--ADD INTERACTIONS WITHIN THE NEW COMPLEX")
          for pdb_ppi_file in set_of_ppis:
            pdb_ppi=os.path.basename(pdb_ppi_file)
            print("\t--PPI file "+pdb_ppi)
            protein_1  = pdb_ppi.split("::")[0]
            protein_2  = pdb_ppi.split("::")[1].split("_")[0]
            sequence_1 = sequences[protein_1]
            sequence_2 = sequences[protein_2]
            branch_list=[]
            if not all_conformers and tuple([protein_1,protein_2]) in done: continue
            for tf_complex in new_tf_complex_list_start:
                print("\t-- Complex "+tf_complex.get_name())
                new_complexes, new_restraints = tf_complex.add_ppi(os.path.join(output_PPI,pdb_ppi),protein_1,protein_2,sequence_1,sequence_2,condition="old",maximum=max_branches)
                if len(new_complexes) > 0 or new_restraints>0:
                   branch_list.extend([x for x in new_complexes])
                   done.add(tuple([protein_1,protein_2]))
                   done.add(tuple([protein_2,protein_1]))
            if len(branch_list)>0:
              #Transform the full list into unique new complexes and extend the list of newest complexes
              new_tf_complex_list.extend(branch_list)
              #Shrink the list to unique complexes and reorder them in a set
              new_tf_complex_set=renamed_complexes_set(new_tf_complex_list)
              #Get the new list of complexes and continiue with next PPIs
              new_tf_complex_list = new_tf_complex_set[:]
            for tf_complex in new_tf_complex_set:
                topology = tf_complex.get_topology()
                for pdb_topo in topology["pdb_fn"]:
                    if pdb_ppi == pdb_topo:
                       print("\t\t--Check topology: add PDB file "+pdb_topo)
                       shutil.copy(pdb_ppi_file,os.path.join(data_folder,pdb_topo))
        new_tf_complex_set=renamed_complexes_set(new_tf_complex_list)
        new_tf_complex_list = new_tf_complex_set[:]
        if len(new_tf_complex_list) == len(tf_complex_list): 
            tf_complex_list = new_tf_complex_list[:]
            end_iteration = True
        if max_iterations is not None:
            number_of_iterations = number_of_iterations + 1
            if number_of_iterations > int(max_iterations):
                tf_complex_list = new_tf_complex_list[:]
                end_iteration   = True
        
    # CLEAN COMPLEXES
    # RESULTS
    print("WRITE INPUT MODELS FOR INTEGRATIVE MODELLING")
    done=[]
    fasta_file=os.path.join(output_DNA_fragments,"all_sequences.fa")
    for tf_complex in tf_complex_list:

          # Be non retsrictive on restraints
          tf_complex.set_restrict_equal(False)
      
          skip=False
          for t2_complex in done:
              if skip: continue
              if tf_complex == t2_complex: 
                 print("\t--REMOVE REDUNDANT COMPLEX "+tf_complex.get_name()+" WITH "+t2_complex.get_name())
                 skip=True

          if skip: continue

          done.append(tf_complex)

          output=tf_complex.get_name()

          print("\t--Write complex "+output)

          if os.path.exists(output) and os.path.isfile(output):
             os.remove(output+".txt")
          if os.path.exists(output+".dna"+".csv") and os.path.isfile(output+".dna"+".csv"):
             os.remove(output+".dna"+".csv")
          if os.path.exists(output+".protein_dna"+".csv") and os.path.isfile(output+".protein_dna"+".csv"):
             os.remove(output+".protein_dna"+".csv")

          #get protein inner restraints
          if protein_self_restraints:
             tf_complex.extend_all_inner_restraints()
          #Trim the number of PPI restraints
          if ratio_trimming < 1.0:
             print("\t\t--Reduce the number of restraints to %s"%str(ratio_trimming))
             tf_complex.trim_restraints(ratio= ratio_trimming)
          #Write old format
          #tf_complex.write_restraints(output+".dna"+".txt","dnadna",kd_dna,w_dna,sd_dna)
          #Write as CSV table for DNA
          table_dna=tf_complex.restraints_to_csv("dnadna",kd_dna,w_dna,sd_dna)
          #table_dna.to_csv(output+".dna"+".csv")
          
          #Write as CSV table for protein-dna
          table_prot=tf_complex.restraints_to_csv("protdna",kd_PNI,w_PNI,sd_PNI)
          #table_prot.to_csv(output+".protdna.csv")
          
          #Write as CSV table for protein-protein
          table_ppi =tf_complex.restraints_to_csv("protprot",kd_PPI,w_PPI,sd_PPI)
          #table_ppi.to_csv(output+".protein"+".csv")
          
          #table=tf_complex.restraints_to_csv("all",kd_PNI,w_PNI,sd_PNI)
          #merge by the union of tables
          if table_dna.shape[1] == table_prot.shape[1]:
             table_core = table_dna.merge(table_prot,how="outer")
          if table_core.shape[1] == table_ppi.shape[1]:
             table      = table_core.merge(table_ppi,how="outer")
          else:
             table      = table_core
          table.to_csv(os.path.join(data_folder,output+".restraints"+".csv"))
          tf_complex.complete_topology(bead_size=5,domain="loops")
          tf_complex.write_topology(data_folder,fasta_file,new=True)

    try:
        shutil.rmtree(dummy)
    except:
        print("Cannot cleand dummy folder")



