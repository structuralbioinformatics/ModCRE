import os, sys, re
from collections import Counter
import ConfigParser
import itertools
import numpy as np
import optparse
import shutil
import random
import threader
#import matplotlib as mpl 
## agg backend is used to create plot as a .png file
#mpl.use('Agg')
#import matplotlib.pyplot as plt 

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
from SBI.data import aminoacids1to3, aminoacids_polarity_boolean, nitrogenous_bases
from SBI.structure import PDB

# Import my modules #
import contacts, dssp, interface,  spotentials, tmalign, triads, x3dna, threading_to_triads
import pwm_pbm as PWM


#-------------#
# Functions   #
#-------------#


def scale_energies(energies, split_potential, min_score, max_score):

    norm_e = 0.0

    if split_potential == "3d":
        norm_e = PWM.scale(energies["3d"], max_score, min_score)
    if split_potential == "3dc":
        norm_e = PWM.scale(energies["3dc"], max_score, min_score)
    if split_potential == "local":
        norm_e = PWM.scale(-energies["local"], max_score, min_score)
    if split_potential == "pair":
        norm_e = PWM.scale(-energies["pair"], max_score, min_score)
    if split_potential == "s3dc":
        norm_e = PWM.scale(-energies["s3dc"], max_score, min_score)
    if split_potential == "s3dc_dd":
        norm_e = PWM.scale(-energies["s3dc_dd"], max_score, min_score)
    if split_potential == "s3dc_di":
        norm_e = PWM.scale(energies["s3dc_di"], max_score, min_score)


    return norm_e



def parse_data_scoring(output_file,pdb_dir,pbm_dir,dummy_dir,fragment_restrict,binding_restrict,input_pdb_file,input_threading_file,threading,template,random,dna_background,radius, potential_file, split_potential, auto_mode, family_potentials, pbm_potentials, score_threshold, taylor_approach, pmf , bins, known, norm,  verbose):

    # Initialize #
    families = {}
    # For family potentials... #
    if pdb_dir is not None:
     if not os.path.exists(os.path.join(pdb_dir, "families.txt")):
       sys.stdout.write("Families file %s is not used\n"%(os.path.join(pdb_dir, "families.txt")))
     else:
      for line in functions.parse_file(os.path.join(pdb_dir, "families.txt")):
        if line.startswith("#"): continue
        pdb_chain, family = line.split(";")
        families[pdb_chain] = family

    if threading:
       # Get triads, pdb etc from threading file #
       triads_obj, pdb_obj, x3dna_obj = threading_to_triads.threading_triads(threading_file=input_threading_file, pdb_dir= pdb_dir, template=template, verbose=verbose, dummy_dir=dummy_dir)
       # Get data for random scores
       if random >0:
          random_triads=[]
          background_sequence=""
          for line in functions.parse_file(dna_background):
              if line.startswith(">"):continue
              background_sequence += line.strip()
          # Generate a random threading file
          thr_obj = threader.Threaded(threading_file=input_threading_file) 
          thr_obj._check_parsing()
          if verbose: sys.stdout.write("\t-- Generate random DNA sequences ...\n")
          for nran in range(random):
            random_threading_file=os.path.abspath(os.path.join(dummy_dir,"random_threading_file_"+str(nran)+"_from_"+os.path.basename(input_threading_file)))
            if not os.path.exists(random_threading_file):
              dna={}
              kmers=thr_obj.get_kmers()
              for dna_seq,binding in kmers.iteritems():
                start=random.randint(0,len(background_sequence)-len(dna_seq))
                dna.setdefault(background_sequence[start:start+len(dna_seq)],binding)
              thr_obj.set_dna(dna)
              thr_obj.set_dna_fixed(dna)
              random_threading_file=os.path.abspath(os.path.join(dummy_dir,"random_threading_file_"+str(nran)+"_from_"+os.path.basename(input_threading_file)))
              thr_obj.write(random_threading_file)
            if verbose: sys.stdout.write("\t-- %dth threading to get random triads ...\n"%(nran))
            random_triads_obj,dummy_pdb_obj, dummy_x3dna_obj = threading_to_triads.threading_triads(threading_file=random_threading_file, pdb_dir= pdb_dir, template=template, verbose=verbose, dummy_dir=dummy_dir)
            random_triads.append(random_triads_obj)
              
    else:
       # Get PDB object #
       if verbose: sys.stdout.write("\t-- Reading PDB file %s ...\n"%(os.path.basename(input_pdb_file)))
       pdb_obj = PDB(input_pdb_file)
       # Get DSSP object #
       if verbose: sys.stdout.write("\t\t-- Get DSSP ...\n")
       dssp_obj = dssp.get_dssp_obj(input_pdb_file)
       # Get X3DNA object #
       if verbose: sys.stdout.write("\t\t-- Get DNA object ...\n")
       x3dna_obj = x3dna.get_x3dna_obj(input_pdb_file, dummy_dir)
       # Get contacts object #
       if verbose: sys.stdout.write("\t\t-- Get contacts ...\n")
       contacts_obj = contacts.get_contacts_obj(pdb_obj, x3dna_obj)
       # Skip if no contacts #
       if len(contacts_obj._contacts) == 0: raise ValueError("No protein-DNA contacts found!")
       # Get triads object #
       triads_obj = triads.get_triads_obj(pdb_obj, dssp_obj, x3dna_obj, contacts_obj)

    # Load statistical potential
    if radius <= 0: radius=float(config.get("Parameters", "max_contact_distance"))
    if verbose: sys.stdout.write("\t-- Load potentials...\n")
    potentials, thresholds, radii, structural_homologs_by_chain = PWM.load_statistical_potentials(pdb_obj, pdb_dir, pbm_dir, families, radius, potential_file, split_potential, auto_mode, family_potentials, pbm_potentials, score_threshold, taylor_approach, pmf , bins, known, None, dummy_dir,verbose)
    # Calculate all energy scores #
    if verbose: sys.stdout.write("\t-- Calculate scores...\n")
    scr=score()
    scr.calculate_energies_and_binding(triads_obj, x3dna_obj, pdb_obj, potentials,  thresholds, radii, radius, fragment_restrict, binding_restrict, dummy_dir)
    # Generate random scores
    random_scores=[]
    if random>0:
       if verbose: sys.stdout.write("\t-- Calculate random energies...\n")
       for nran in range(random):
           random_triads_obj=random_triads[nran]
           if verbose: sys.stdout.write("\t\t-- %dth threaded triads...\n"%(nran))
           random_scr = score()
           random_scr.calculate_energies_and_binding(random_triads_obj, x3dna_obj, pdb_obj, potentials,  thresholds, radii, radius, fragment_restrict, binding_restrict, dummy_dir)
           random_scores.append(random_scr)
    if not verbose:
       shutil.rmtree(dummy_dir)

    return  scr,random_scores


#-------------#
# Classes     #
#-------------#

class  score(object):
    """
    Class of energy scores and profiles of scores per nucleotide
    """

    def __init__(self,binding_site=None, energies=None, energies_per_nucleotide=None, norm_energies=None, potential=None):
        
        self._energies={}
        self._energies_per_nucleotide={}
        self._norm_energies={}
        self._binding_site={}
        #Use all energies by default, but retrieves s3dc_dd with get_score by default
        self._potential="all"

        if energies_per_nucleotide is not None:
           self._energies=energies
        if energies is not None:
           self._energies_per_nucleotide=energies_per_nucleotide
        if norm_energies is not None:
           self._norm_energies=norm_energies
        if binding_site is not None:
           self._binding_site=binding_site
        if potential is not None:
           self._potential=potential

    def get_energies(self):
        return self._energies

    def get_energies_per_nucleotide(self):
        return self._energies_per_nucleotide

    def get_score_per_nucleotide(self, normal=False, potential=None):
        #Select a score_per_nucleotide, default is defined by class or "s3dc_dd" if "all"
        if potential is None: split_potential = self._potential
        else:                 split_potential = potential
        if split_potential == "all": split_potential = "s3dc_dd"
        skip=False
        if not self._energies.has_key(potential): skip=True
        if normal and abs(float(self._energies.has_key(potential))) < 1.0e-10: skip=True
        score_per_nucleotide={}
        for nuc in sorted(self._energies_per_nucleotide.keys()):
            if normal: 
                       if skip: score_per_nucleotide.setdefault(nuc,0.0)
                       else:    score_per_nucleotide.setdefault(nuc,float(self._energies_per_nucleotide[nuc][potential]/self._energies[potential]))
            else:      score_per_nucleotide.setdefault(nuc,float(self._energies_per_nucleotide[nuc][potential]))
        return score_per_nucleotide

    def get_score(self, normal=False, potential=None):
        #Select a single score, default is defined by class or "s3dc_dd" if "all"
        if potential is None: split_potential = self._potential
        else:                 split_potential = potential
        if split_potential == "all": split_potential = "s3dc_dd"
        if normal:
            if self._normal_energies.has_key(split_potential): return  self._normal_energies[split_potential]
            else: return  0.0
        else:
            if self._energies.has_key(split_potential): return  self._energies[split_potential]
            else: return  0.0

    def get_norm_energies(self):
        return self._norm_energies

    def get_binding_site(self):
        return self._binding_site

    def get_potential(self):
        return self._potential

    def calculate_energies_and_binding(self, triads_obj, x3dna_obj, pdb_obj, potentials, thresholds, radii, radius=0, fragment_restrict=None, binding_restrict=None, dummy_dir="/tmp"):
        # Initialize #
        split_potential = self._potential
        E_3d = 0
        E_3dc = 0
        E_local = 0
        E_pair = 0
        E_s3dc = 0
        E_s3dc_dd = 0
        E_s3dc_di = 0
        binding_site = {}
        nucleotides = list("ACGT")
        dna_idx={}
        energies={}
        energies_per_nucleotide={}
        if radius <= 0: radius=float(config.get("Parameters", "max_contact_distance"))
        # Get original DNA sequence indexes
        basepairs=x3dna_obj.get_basepairs()
        for basepair in basepairs.iterkeys():
              (fwd_pdb_chain, fwd_residue_num), (rev_pdb_chain, rev_residue_num) = basepairs.get(basepair)
              dna_idx.setdefault((fwd_pdb_chain,fwd_residue_num),int(basepair))
        # For each PDB chain... #
        twonucleotides={}
        for pdb_chain in sorted(potentials):
              # For each contact... #
              for triad_obj in triads_obj.get_triads():
                  # Initialize #
                  a_oa, b_ob, distance, residue_A, residue_B = triad_obj.return_as_string().split(";")
                  chain, residue_num = residue_A.split("-")
                  #check restriction in TF
                  if fragment_restrict is not None:
                   if fragment_restrict.has_key(pdb_chain):
                     belongs=False
                     for interval in fragment_restrict.get(pdb_chain):
                         if int(residue_num) >= int(interval[0]) and int(residue_num) <= int(interval[1]): belongs=True
                         if belongs: break
                     if not belongs: continue
                  dinucleotide = PWM.get_triad_dinucleotide(triad_obj, x3dna_obj)
                  bp_1F,bp_1R,bp_2F,bp_2R=residue_B.split(",")
                  # Skip if not the right chain... #
                  if chain != pdb_chain: continue
                  # Add dinucleotide to binding site and PDB numbering to list of twonucleotides#
                  chain_dna,base_number_1=bp_1F.split("-")
                  chain_dna,base_number_2=bp_2F.split("-")
                  binucleotide = (dna_idx[(chain_dna,int(base_number_1))],dna_idx[(chain_dna,int(base_number_2))])
                  #Check restriction on DNA site
                  if binding_restrict is not None:
                     belongs=False
                     for interval in binding_restrict:
                         if int(binucleotide[0]) >= int(interval[0]) and int(binucleotide[1]) <= int(interval[1]): belongs=True
                         if belongs: break
                     if not belongs: continue
                  # confirm base_numbers from 3xdna object
                  confirm_bases=False
                  if x3dna_obj.has_dinucleotide(dinucleotide):
                     bp1,bp2=x3dna_obj.get_dinucleotide(dinucleotide)
                     if x3dna_obj.has_basepair(bp1) and x3dna_obj.has_basepair(bp2):
                        ((fwd_pdb_chain1, fwd_residue_num1), (rev_pdb_chain1, rev_residue_num1)) = x3dna_obj.get_basepair(bp1)
                        ((fwd_pdb_chain2, fwd_residue_num2), (rev_pdb_chain2, rev_residue_num2)) = x3dna_obj.get_basepair(bp2)
                        if fwd_pdb_chain1==chain_dna and fwd_pdb_chain2==chain_dna and int(fwd_residue_num1)==int(base_number_1) and int(fwd_residue_num2)==int(base_number_2):
                           confirm_bases=True
                  if not confirm_bases:
                     sys.stdout.write("Error: PDB residue numbers of nucleotides are different from 3XDNA definition of basepairs. Check your files\n")
                     exit(0)
                  #Check restriction on distance
                  dab = np.floor(float(distance))
                  if dab > radius: continue
                  #selected for binding site
                  twonucleotides.setdefault(dinucleotide,(int(base_number_1),int(base_number_2)))
                  binding_site.setdefault(dinucleotide, []).append(triad_obj)
                  energies_per_nucleotide.setdefault(dinucleotide, {"3d":0.0, "3dc":0.0, "local":0.0, "pair":0.0, "s3dc":0.0, "s3dc_dd":0.0, "s3dc_di":0.0})
                  # Skip if environment could not be obtained for amino acid #
                  if "None" in a_oa: continue
                  # Skip if environment could not be obtained for dinucleotide #
                  if "None" in b_ob: continue
                  # The following arguments arguments are defined as in the papers #
                  a, hydrophobicity, degree_of_exposure, secondary_structure = a_oa.split("-")
                  oa = "%s-%s-%s" % (hydrophobicity, degree_of_exposure, secondary_structure)
                  b, nitrogenous_bases, dna_strand, dna_groove, dna_chemical_group = b_ob.split("-")
                  ob = "%s-%s-%s-%s" % (nitrogenous_bases, dna_strand, dna_groove, dna_chemical_group)
                  a_b = "%s;%s" % (a, b)
                  a_b_oa_ob = "%s;%s" % (a_oa, b_ob)
                  oa_ob = "%s;%s" % (oa, ob)
                  ### Compute potentials for our input structure ### 
                  if (split_potential == "3d") or (split_potential == "all"):
                      if potentials[chain].get_score("3d", distance=dab) is not None:
                          E_3d += potentials[chain].get_score("3d", distance=dab)
                          energies_per_nucleotide[dinucleotide]["3d"] += potentials[chain].get_score("3d", distance=dab)
                  if (split_potential == "3dc") or (split_potential == "all"):
                      if potentials[chain].get_score("3dc", key=oa_ob, distance=dab) is not None:
                          E_3dc += potentials[chain].get_score("3dc", key=oa_ob, distance=dab)
                          energies_per_nucleotide[dinucleotide]["3dc"] += potentials[chain].get_score("3dc", key=oa_ob, distance=dab)
                  if (split_potential == "local") or (split_potential == "all"):
                      if potentials[chain].get_score("local", key=a_oa) is not None:
                          E_local += potentials[chain].get_score("local", key=a_oa)
                          energies_per_nucleotide[dinucleotide]["local"] += potentials[chain].get_score("local", key=a_oa)
                  if (split_potential == "pair") or (split_potential == "all"):
                      if potentials[chain].get_score("pair", key=a_b, distance=dab) is not None:
                          E_pair += potentials[chain].get_score("pair", key=a_b, distance=dab)
                          energies_per_nucleotide[dinucleotide]["pair"] += potentials[chain].get_score("pair", key=a_b, distance=dab)
                  if (split_potential == "s3dc") or (split_potential == "all"):
                      if potentials[chain].get_score("s3dc", key=a_b_oa_ob, distance=dab) is not None:
                          E_s3dc += potentials[chain].get_score("s3dc", key=a_b_oa_ob, distance=dab)
                          energies_per_nucleotide[dinucleotide]["s3dc"] += potentials[chain].get_score("s3dc", key=a_b_oa_ob, distance=dab)
                  if (split_potential == "s3dc_dd") or (split_potential == "all"):
                      if potentials[chain].get_score("s3dc_dd", key=a_b_oa_ob, distance=dab) is not None:
                          E_s3dc_dd += potentials[chain].get_score("s3dc_dd", key=a_b_oa_ob, distance=dab)
                          energies_per_nucleotide[dinucleotide]["s3dc_dd"] += potentials[chain].get_score("s3dc_dd", key=a_b_oa_ob, distance=dab)
                  if (split_potential == "s3dc_di") or (split_potential == "all"):
                      if potentials[chain].get_score("s3dc_di", key=a_b_oa_ob) is not None:
                          E_s3dc_di += potentials[chain].get_score("s3dc_di", key=a_b_oa_ob)
                          energies_per_nucleotide[dinucleotide]["s3dc_di"] += potentials[chain].get_score("s3dc_di", key=a_b_oa_ob)
        # Get profile of energies per nucleotide (numbered as in PDB)#
        energies_per_nucleotide_processed = {}
        nucleotides=[]
        for dinucleotide,bps in twonucleotides.iteritems():
              base1,base2=bps
              for i in (base1,base2):
                    nucleotides.append(i)
                    if not  energies_per_nucleotide_processed.has_key(i):   
                      energies_per_nucleotide_processed[i]={"3d":0.0,"3dc":0.0,"local":0.0,"pair":0.0,"s3dc":0.0,"s3dc_dd":0.0,"s3dc_di":0.0}
                    energies_per_nucleotide_processed[i]["3d"] += energies_per_nucleotide[dinucleotide]["3d"]
                    energies_per_nucleotide_processed[i]["3dc"] += energies_per_nucleotide[dinucleotide]["3dc"]
                    energies_per_nucleotide_processed[i]["local"] += energies_per_nucleotide[dinucleotide]["local"]
                    energies_per_nucleotide_processed[i]["pair"] += energies_per_nucleotide[dinucleotide]["pair"]
                    energies_per_nucleotide_processed[i]["s3dc"] += energies_per_nucleotide[dinucleotide]["s3dc"]
                    energies_per_nucleotide_processed[i]["s3dc_dd"] += energies_per_nucleotide[dinucleotide]["s3dc_dd"]
                    energies_per_nucleotide_processed[i]["s3dc_di"] += energies_per_nucleotide[dinucleotide]["s3dc_di"]
        # Energies per Nucleotide
        repeats=Counter(nucleotides)
        for i,s in repeats.iteritems():
              for sp in energies_per_nucleotide_processed[i].iterkeys():
                  energies_per_nucleotide_processed[i][sp]=energies_per_nucleotide_processed[i][sp]/s
        # Energies
        energies      = {"3d":E_3d, "3dc":E_3dc, "local":E_local, "pair":E_pair, "s3dc":E_s3dc, "s3dc_dd":E_s3dc_dd, "s3dc_di":E_s3dc_di}
        # Normalized energies
        norm_energies = {"3d":0.0, "3dc":0.0, "local":0.0, "pair":0.0, "s3dc":0.0, "s3dc_dd":0.0, "s3dc_di":0.0}
        if split_potential == "all":
            for sp in ["3d", "3dc", "local", "pair", "s3dc", "s3dc_dd", "s3dc_di"]:
                msa_obj = PWM.get_msa_obj(triads_obj, x3dna_obj, potentials, radii, fragment_restrict, binding_restrict, sp, thresholds)
                min_score, max_score = PWM.get_min_max_scores(msa_obj, binding_site, x3dna_obj,  potentials,  sp, dummy_dir)
                norm_e = scale_energies(energies, sp, min_score, max_score)
                norm_energies[sp] = norm_e
        else:
            msa_obj = PWM.get_msa_obj(triads_obj, x3dna_obj, potentials, radii, fragment_restrict, binding_restrict, split_potential, thresholds)
            min_score, max_score = PWM.get_min_max_scores(msa_obj, binding_site, x3dna_obj, potentials,  split_potential, dummy_dir)
            norm_e = scale_energies(energies, split_potential, min_score, max_score)
            norm_energies[split_potential] = norm_e

        # Constructor
        self._energies                = energies
        self._binding_site            = binding_site
        self._energies_per_nucleotide = energies_per_nucleotide_processed
        self._norm_energies           = norm_energies



    def write(self, output_file, protein_name=None, normal=False, overwrite=False):

        # Create output file #
        if os.path.exists(output_file) and overwrite:
            os.remove(output_file)
        out = open(output_file, "a")

        # Get variables to write
        split_potential         = self._potential
        binding_site            = self._binding_site
        energies                = self._energies
        energies_per_nucleotide = self._energies_per_nucleotide
        norm_energies           = self._norm_energies
        if protein_name is None:  protein_name=".".join([x for x in os.path.basename(ouput_file).split(".")[0:-1]])

        # Write energies
        out.write("\n##############################################################################################################################################\n")
        out.write("##############################################################################################################################################\n")
        out.write("\n#%30s\t"%"Protein_scores")
        out.write("\t%20s\t%20s\t"%("Binding_site_length","Starting_site"))
        if split_potential == "3d" or split_potential == "all":       out.write("%15s\t"%"E_3d")
        if split_potential == "3dc" or split_potential == "all":      out.write("%15s\t"%"E_3dc")
        if split_potential == "local" or split_potential == "all":    out.write("%15s\t"%"E_local")
        if split_potential == "pair" or split_potential == "all":     out.write("%15s\t"%"E_pair")
        if split_potential == "s3dc" or split_potential == "all":     out.write("%15s\t"%"E_s3dc")
        if split_potential == "s3dc_dd" or split_potential == "all":  out.write("%15s\t"%"E_s3dc_dd")
        if split_potential == "s3dc_di" or split_potential == "all":  out.write("%15s\t"%"E_s3dc_di")
        out.write("\n")
        out.write("%30s\t"%protein_name)
        out.write("\t%20d\t%20d\t"%(len(binding_site)+1,min([int(x) for x in binding_site.iterkeys()])))
        if split_potential == "3d" or split_potential == "all": out.write("%15.3f\t" % (float(energies["3d"])))
        if split_potential == "3dc" or split_potential == "all": out.write("%15.3f\t" % (float(energies["3dc"])))
        if split_potential == "local" or split_potential == "all": out.write("%15.3f\t" % (float(energies["local"])))
        if split_potential == "pair" or split_potential == "all": out.write("%15.3f\t" % (float(energies["pair"])))
        if split_potential == "s3dc" or split_potential == "all": out.write("%15.3f\t" % (float(energies["s3dc"])))
        if split_potential == "s3dc_dd" or split_potential == "all": out.write("%15.3f\t" % (float(energies["s3dc_dd"])))
        if split_potential == "s3dc_di" or split_potential == "all": out.write("%15.3f\t" % (float(energies["s3dc_di"])))
        out.write("\n")

        # Write Normalized energies
        if norm_energies is not None and len(norm_energies.keys())>0 and normal:
          out.write("#%29s\t"%"Normalized (min-max scale)")
          if split_potential == "3d" or split_potential == "all":       out.write("%15s\t"%"E_3d")
          if split_potential == "3dc" or split_potential == "all":      out.write("%15s\t"%"E_3dc")
          if split_potential == "local" or split_potential == "all":    out.write("%15s\t"%"E_local")
          if split_potential == "pair" or split_potential == "all":     out.write("%15s\t"%"E_pair")
          if split_potential == "s3dc" or split_potential == "all":     out.write("%15s\t"%"E_s3dc")
          if split_potential == "s3dc_dd" or split_potential == "all":  out.write("%15s\t"%"E_s3dc_dd")
          if split_potential == "s3dc_di" or split_potential == "all":  out.write("%15s\t"%"E_s3dc_di")
          out.write("\n")
          out.write("%30s\t"%protein_name)
          if split_potential == "3d" or split_potential == "all": out.write("%15.3f\t" % float(norm_energies["3d"]))
          if split_potential == "3dc" or split_potential == "all": out.write("%15.3f\t" % float(norm_energies["3dc"]))
          if split_potential == "local" or split_potential == "all": out.write("%15.3f\t" % float(norm_energies["local"]))
          if split_potential == "pair" or split_potential == "all": out.write("%15.3f\t" % float(norm_energies["pair"]))
          if split_potential == "s3dc" or split_potential == "all": out.write("%15.3f\t" % float(norm_energies["s3dc"]))
          if split_potential == "s3dc_dd" or split_potential == "all": out.write("%15.3f\t" % float(norm_energies["s3dc_dd"]))
          if split_potential == "s3dc_di" or split_potential == "all": out.write("%15.3f\n" % float(norm_energies["s3dc_di"]))

        # Write Energies per nucleotide
        out.write("\n\n")
        out.write("#Nucleotide scores\n")
        out.write("%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\n"%("#Nucleotide","E_3d","E_3dc","E_local","E_pair","E_s3dc","E_s3dc_dd","E_s3dc_di"))
        for nuc in sorted(energies_per_nucleotide.keys()):
           out.write("%15s\t" % str(nuc))
           out.write("%15.3f\t" % float(energies_per_nucleotide[nuc]["3d"]) )
           out.write("%15.3f\t" % float(energies_per_nucleotide[nuc]["3dc"]))
           out.write("%15.3f\t" % float(energies_per_nucleotide[nuc]["local"]))
           out.write("%15.3f\t" % float(energies_per_nucleotide[nuc]["pair"]))
           out.write("%15.3f\t" % float(energies_per_nucleotide[nuc]["s3dc"]))
           out.write("%15.3f\t" % float(energies_per_nucleotide[nuc]["s3dc_dd"]))
           out.write("%15.3f\t" % float(energies_per_nucleotide[nuc]["s3dc_di"]))
           out.write("\n")

        # Write Normalized energies per nucleotide
        if norm_energies is not None and len(norm_energies.keys())>0 and normal:
            out.write("#Normalized scores per nucleotide (ratio scaling)\n")
            out.write("%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\t%15s\n"%("#Nucleotide","E_3d","E_3dc","E_local","E_pair","E_s3dc","E_s3dc_dd","E_s3dc_di"))
            for nuc in sorted(energies_per_nucleotide.keys()):
                out.write("%15s\t" % str(nuc))
                if split_potential == "3d" or split_potential == "all":
                    try: 
                        out.write("%15.3f\t" % float(energies_per_nucleotide[nuc]["3d"]/energies["3d"])) 
                    except:
                        out.write("%15.3f\t" %(0.0))
                else:
                    out.write("%15.3f\t" %(0.0))
                if split_potential == "3dc" or split_potential == "all":
                    try: 
                        out.write("%15.3f\t" % float(energies_per_nucleotide[nuc]["3dc"]/energies["3dc"]))
                    except:
                        out.write("%15.3f\t" %(0.0))
                else:
                    out.write("%15.3f\t" %(0.0))
                if split_potential == "local" or split_potential == "all":
                    try: 
                        out.write("%15.3f\t" % float(energies_per_nucleotide[nuc]["local"]/energies["local"]))
                    except:
                        out.write("%15.3f\t" %(0.0))
                else:
                    out.write("%15.3f\t" %(0.0))
                if split_potential == "pair" or split_potential == "all":
                    try: 
                        out.write("%15.3f\t" % float(energies_per_nucleotide[nuc]["pair"]/energies["pair"]))
                    except:
                        out.write("%15.3f\t" %(0.0))
                else:
                    out.write("%15.3f\t" %(0.0))
                if split_potential == "s3dc" or split_potential == "all":
                    try: 
                        out.write("%15.3f\t" % float(energies_per_nucleotide[nuc]["s3dc"]/energies["s3dc"]))
                    except:
                        out.write("%15.3f\t" %(0.0))
                else:
                    out.write("%15.3f\t" %(0.0))
                if split_potential == "s3dc_dd" or split_potential == "all":
                    try: 
                        out.write("%15.3f\t" % float(energies_per_nucleotide[nuc]["s3dc_dd"]/energies["s3dc_dd"]))
                    except:
                        out.write("%15.3f\t" %(0.0))
                else:
                    out.write("%15.3f\t" %(0.0))
                if split_potential == "s3dc_di" or split_potential == "all":
                    try: 
                        out.write("%15.3f\t" % float(energies_per_nucleotide[nuc]["s3dc_di"]/energies["s3dc_di"]))
                    except:
                        out.write("%15.3f\t" %(0.0))
                else:
                    out.write("%15.3f\t" %(0.0))
                out.write("\n")
        out.close()




#-------------#
# Options     #
#-------------#

def parse_options():
    '''
    This function parses the command line arguments and returns an optparse object.
    '''

    parser = optparse.OptionParser("Usage: scorer.py [--dummy=DUMMY_DIR] -i INPUT_FILE [-l LABEL -o OUTPUT_DIR --pbm=PBM_dir] --pdb=PDB_DIR [-m -v --threading] [-a -f -p -s SPLIT_POTENTIAL -t THRESHOLD -k -b --taylor --file POTENTIAL --radius RADIUS]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="DUMMY_DIR")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input PDB|THREADED file. Mandatory", metavar="INPUT_FILE")
    parser.add_option("-l", action="store", type="string", dest="label", help="Label to organize the output files as 'label.energies.txt' (default = scores.txt)", metavar="LABEL")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_name", help="Output directory (default = ./)", metavar="OUTPUT_DIR")
    parser.add_option("--pbm", action="store", type="string", default=None, dest="pbm_dir", help="PBM directory (i.e. output dir from pbm.py). This is ,mandatory unless using option --file", metavar="PBM_DIR")
    parser.add_option("--pdb", action="store", type="string", default=None, dest="pdb_dir", help="PDB directory (i.e. output dir from pdb.py). This is ,mandatory unless using option --template", metavar="PDB_DIR")
    parser.add_option("-m", "--pmf", default=False, action="store_true", dest="pmf", help="Use of raw mean-force potentials with no Z-scoring (default = False)")
    parser.add_option("--fragment", default=None, action="store", type="string", dest="fragment", help="Fragment of protein to apply the potential. Format is 'a-b;c-d': two regions between residues a-b and c-d. (Default is None it applies to all amino-acids)")
    parser.add_option("-n", "--norm", default=False, action="store_true", dest="norm", help="Normalize scores using a min-max scaling")
    parser.add_option("--binding", default=None, action="store", type="string", dest="binding_restrict", help="Binding site of DNA to apply the potential (uses basepair order from 3XDNA). Format is 'a-b;c-d': two regions between nucleotides a-b and c-d of the forward chain (first in PDB). (Default is None it applies to all nucleotides)")
    parser.add_option("--threading", default=False, action="store_true", dest="threading", help="Input file is a threading file of a PDB structure that either exist in the PDB folder of ModCRE or is given with option '--template' (default = False)")
    parser.add_option("--random", default=0,type=int,action="store", dest="random",help="Number of random DNA sequences to generate a distribution of random scores (ONLY APPLICABLE ON THREADED FILES)", metavar="INTEGER")
    parser.add_option("--background",default=None, action="store", type="string", dest="dna_background", help="File of a DNA sequence in FASTA format. This is used as background to generate a distribution of random scores", metavar="FASTA")
    parser.add_option("--template", action="store", default=None, type="string", dest="template", help="PDB structure of the template on threading. Option threading must be active", metavar="TEMPLATE")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")

    group = optparse.OptionGroup(parser, "Statistical potentials", "Select your statistical potentials of choice. By default it uses S3DCdd general potential derived from PDB (the simplest one). In \"--auto\" mode, the program uses S3DCdd family potentials derived from both PDB and/or PBM data and/or approached by Taylor as selected in Potentials configuration file. In case family potentials cannot be applied, the program uses general potentials derived from both PDB and PBM data and approached by Taylor. \"-a\" option overrides options \"-f\", \"-p\" and \"-t\".")
    group.add_option("-a", "--auto", default=False, action="store_true", dest="auto_mode", help="Automate the selection of statistical potentials (default = False)")
    group.add_option("-f", "--family", default=False, action="store_true", dest="family_potentials", help="Use family potentials (default = False)")
    group.add_option("-p", default=False, action="store_true", dest="pbm_potentials", help="Use potentials derived from both PBM + PDB data (default = False)")
    group.add_option("-s", default="s3dc_dd", action="store", type="string", dest="split_potential", help="Split-potential to be used (3d, 3dc, s3dc, s3dc_dd, s3dc_di, pair; default = s3dc_dd)", metavar="{string}")
    group.add_option("-t", action="store", type="float", dest="score_threshold", help="Threshold on the scaled score to consider positive k-mers (default = 0.95)", metavar="{float}")
    group.add_option("-k","--known", default=False, action="store_true", dest="known", help="The name is of a known PDB file, with format 'code_chain' (default = False)")
    group.add_option("--taylor", default=False, action="store_true", dest="taylor_approach", help="Approach PMF by Taylor (default = False)")
    group.add_option("-b", "--bins", default=False, action="store_true",  dest="bins", help="Computate the potentials: by bins (if selected) or accumulative (default).", metavar="{boolean}")
    group.add_option("--file",default=None, action="store", type="string", dest="potential_file", help="Use potentials from specific file (default = None)", metavar="{string}")
    group.add_option("--radius",default=0, action="store", type="float", dest="radius", help="Maximum contact distance to calculate interactions (default=0 implies the use of 'max_contact_distance' from configuration", metavar="{string}")

    parser.add_option_group(group)

    
    (options, args) = parser.parse_args()

    if options.input_file is None or (options.pdb_dir is None and options.template is None) or (options.pbm_dir is None and options.potential_file is None):
        parser.error("missing arguments: type option \"-h\" for help")
    if options.random > 0: 
       if options.dna_background is None:
        parser.error("missing DNA background file: type option \"-h\" for help")
       if options.template is None:
        parser.error("missing template structure : type option \"-h\" for help")
       if not os.path.exists(options.dna_background) or not os.path.exists(options.template):
        parser.error("empty DNA background file or template structure file: type option \"-h\" for help")
       if not options.threading:
        parser.error("threading flag is not active: type option \"-h\" for help")
    if not re.search("^3d$|^3dc$|^local$|^pair$|^s3dc$|^s3dc_dd$|^s3dc_di$|^all$", options.split_potential):
        parser.error("incorrect value for -s argument: type option \"-h\" for help")
        exit(0)
    return options


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()
    dummy_dir = options.dummy_dir
    if not dummy_dir.startswith("/"): dummy_dir = os.path.abspath(options.dummy_dir)
    pdb_dir = options.pdb_dir
    if pdb_dir is not None:
       if not pdb_dir.startswith("/"): pdb_dir = os.path.abspath(options.pdb_dir)
    pbm_dir = options.pbm_dir
    if pbm_dir is not None:
       if not pbm_dir.startswith("/"): pbm_dir = os.path.abspath(options.pbm_dir)
    if options.threading:
       input_threading_file=options.input_file
       input_pdb_file=None
       if not input_threading_file.startswith("/"): input_threading_file=os.path.abspath(options.input_file)
    else:
       input_pdb_file=options.input_file
       input_threading_file=None
       if not input_pdb_file.startswith("/"): input_pdb_file=os.path.abspath(options.input_file)
    output_name = options.output_name
    if not output_name.startswith("/"): output_name = os.path.abspath(options.output_name)

    output_file = output_name + ".txt"
    if options.label is not None:
        output_file = output_name +"."+ options.label + ".txt"

    fragment_restrict=None
    if options.fragment is not None:
       fragment_restrict={}
       segment_fragment=options.fragment.split(";")
       for interval in segment_fragment:
           ac,bc = interval.split("-")
           if len(ac.split("_"))>0: a,chain1=ac.split("_")
           if len(bc.split("_"))>0: b,chain2=bc.split("_")
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

    scr, random_scr = parse_data_scoring(output_file,pdb_dir,pbm_dir,dummy_dir,fragment_restrict,binding_restrict,input_pdb_file,input_threading_file,options.threading,options.template,options.random,options.dna_background,options.radius, options.potential_file, options.split_potential, options.auto_mode, options.family_potentials, options.pbm_potentials, options.score_threshold, options.taylor_approach, options.pmf , options.bins, options.known,options.norm,options.verbose)

    # Write the output #
    # Original
    protein_name=".".join([x for x in os.path.basename(options.input_file).split(".")[0:-1]])
    scr.write(output_file, protein_name=protein_name, normal=options.norm, overwrite=True)
    # Random energies
    if options.random>0:
      random_output_file=output_file.strip("txt")+".random"
      if os.path.exists(random_output_file):
       os.remove(random_output_file)
      for random_score in random_scr:
          random_score.write(random_output_file, protein_name=protein_name, normal=options.norm, overwrite=True)

    print("Done")


