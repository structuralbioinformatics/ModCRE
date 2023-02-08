import os, sys, re
from collections import Counter
import ConfigParser
import itertools
import numpy as np
import optparse
import shutil
import random
import cPickle
import pandas as pd
import subprocess
import scipy.stats as stats

#import graphics
import matplotlib as mpl 
# agg backend is used to create plot as a .png file
mpl.use('Agg')
import matplotlib.pyplot as plt 

import plotly
import plotly.graph_objs as go    
from plotly.offline import plot as pltly

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
import contacts, dssp, interface,  spotentials, tmalign, triads, x3dna, scorer, fimo, model_protein, model_dna, threading_to_triads,threader
import pwm_pbm as PWM

#-------------#
# Functions   #
#-------------#

def build_BDNA(dna_sequence,pdb_file):

    try:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        x3dna_path = os.path.join(src_path, config.get("Paths", "x3dna_path"))
        os.environ['X3DNA'] = x3dna_path[:-4]
        # Exec process #
        process = subprocess.check_output([os.path.join(x3dna_path, "fiber"),"-b","-seq=%s"%dna_sequence, pdb_file], stderr=subprocess.STDOUT, env=os.environ)
    except Exception as e:
        raise ValueError("Could not exec X3DNA for %s with error %s" % (pdb_file,e))

    

##############################################################################################################################
#                                                                                                                               
# Fast calculation of scores using a threading file. SAVE models is not allowed
# Calculate the profile by threading the DNA fragment sequences on the PDB native structure to calculate the potentials
#
##############################################################################################################################

def calculate_single_profile_of_thread(dna_file,thresholds,energy_profile,thread_file,output_file,pbm_dir,pdb_dir,families,potential_file, radius,fragment_restrict, binding_restrict, split_potential,auto_mode,family_potentials,pbm_potentials,score_threshold,taylor_approach,pmf,bins,known,meme,reset, dummy_dir,verbose,save,methylation=False):

 try:
#Get the output names
    output_dir  = os.path.dirname(output_file)
    output_name = os.path.basename(output_file)
    if verbose: sys.stdout.write("Output files %s (folder %s)\n"%(output_name,output_dir))

#Get the DNA sequence    
    use_dna=True
    for header, sequence in functions.parse_fasta_file(dna_file):
        if verbose: sys.stdout.write("Read DNA named %s (lenght %d) \n"%(header,len(sequence)))
        if use_dna:
           dna_sequence = sequence
           dna_name = header.split()[0].replace("/","-")
           use_dna=False
    if save:
        pdb_dna_file=os.path.join(output_dir,dna_name+".pdb")
        build_BDNA(dna_sequence, pdb_dna_file)

#Select the format of PWM when methylation
    if meme: methylation_pwm = False
    else:    methylation_pwm = methylation

#Get the PWM
    if not functions.fileExist(output_file+".meme"):
           if verbose:sys.stdout.write("Generate PWM %s ...\n"%(output_file))
           PWM.get_single_pwm(None,thread_file,True,output_file,pbm_dir,pdb_dir,families,potential_file, radius,fragment_restrict, binding_restrict, split_potential,auto_mode,family_potentials,pbm_potentials,score_threshold,taylor_approach,pmf,bins,known,meme,reset, dummy_dir,verbose,methylation_pwm)

#Scan the DNA sequence with the PWM to get FIMO scores
    if verbose: sys.stdout.write("Get FIMO scores\n")
    pwm_file=output_file+".meme"
    max_stored_matches = None
    profile=ProfileProtein(dna_sequence,thresholds)
    #Modify the DNA sequence to standard only nucleotides 
    if methylation_pwm:
       dummy_dna_file = dna_file
    else:
       if not os.path.exists( os.path.join(dummy_dir,"DNA_FILES") ): os.makedirs(os.path.join(dummy_dir,"DNA_FILES"))
       dummy_dna_file = os.path.join(dummy_dir,"DNA_FILES",os.path.basename(dna_file))
       with open(dummy_dna_file,"w") as fo:
         for header_dna, sequence_dna in functions.parse_fasta_file(dna_file):
             sequence_dna_modified = sequence_dna.upper().upper().replace("X","C").replace("O","C").replace("J","G").replace("Q","G")
             fo.write(">%s\n%s\n"%(header_dna,sequence_dna_modified))
       fo.close()
    for threshold in thresholds:
        if verbose: sys.stdout.write("\t-- Threshold %s\n"%threshold)
        fimo_obj = fimo.get_fimo_obj(pwm_file, dummy_dna_file, float(threshold), max_stored_matches, dummy_dir)
        profile.set_fimo_scores_by_threshold(fimo_obj,threshold)
    
#Get triads, pdb etc from threading file #
    if verbose: sys.stdout.write("\t-- Get triads ...\n")
    triads_obj, pdb_obj, x3dna_obj = threading_to_triads.threading_triads(threading_file=input_threading_file, pdb_dir= pdb_dir)
#Get interface object #
    if verbose: sys.stdout.write("\t-- Get interface ...\n")
    interface_obj = interface.get_interface_obj(pdb_obj, x3dna_obj, contacts_obj)


# Load statistical potential
    if radius <= 0: radius=float(config.get("Parameters", "max_contact_distance"))
    if verbose: sys.stdout.write("\t-- Load potentials ...\n")
    potentials, thresholds, radii, structural_homologs_by_chain = PWM.load_statistical_potentials(pdb_obj, pdb_dir, pbm_dir, families, radius, potential_file, energy_profile, auto_mode, family_potentials, pbm_potentials, score_threshold, taylor_approach, pmf , bins, known, None, dummy_dir,verbose)

# Get the dictionary of msa_objs
    msa_objs={}
    split_potentials=["3d","3dc","local","pair","s3dc","s3dc_di","s3dc_dd"]
    for sp in split_potentials:
          if verbose: sys.stdout.write("\t-- Get MSA for potential %s ...\n"%sp)
          msa_obj = PWM.get_msa_obj(triads_obj, x3dna_obj, potentials, radii, fragment_restrict, binding_restrict, sp, thresholds,methylation_pwm)
          msa_objs.setdefault(sp,msa_obj)
# Flush
    sys.stdout.flush()


# Split DNA sequence according to interface
    # All TF-DNA models generated with model_dna are restricted to the DNA interface
    # Get the interface length
    interface_length=interface_obj.get_interface_length()
    if verbose: 
        sys.stdout.write("Interface by PDB\n")
        sys.stdout.write("\t-- Length %d\n"%interface_length)
        sys.stdout.write("\t-- Starting basepairs in structure %d\n"%interface_obj.get_interface_start())
        sys.stdout.write("\t-- Last basepairs in structure %d\n"%interface_obj.get_interface_end())
    if len(dna_sequence)<interface_length:
        sys.stdout.write("Stop. Length of DNA is shorter than interface\n")
        raise ValueError("Error in PROFILE thread: Length of DNA is shorter than interface\n")

    # Get the energy-scores by windows along the DNA sequence
    for i in range(len(dna_sequence)-interface_length+1):

        # Set the DNA fragment
        dna_fragment=dna_sequence[i:i+interface_length]
        if verbose: sys.stdout.write("Score %d: fragment %s\n"%(i+1,dna_fragment))

        # Thread DNA fragment to get new triads 
        if verbose: sys.stdout.write("\t-- Thread DNA fragment on triads ...\n")
        fragment_triads_obj = threading_to_triads.thread_kmer_triads(dna_fragment,triads_obj)

        # Calculate the energy scores and fill the profile#
        if verbose: sys.stdout.write("\t-- Calculate scores with potential: %s ...\n"%energy_profile)
        if energy_profile is None or energy_profile=="all": scr=scorer.score()
        else: scr=scorer.score(potential=energy_profile)
        scr.calculate_energies_and_binding(fragment_triads_obj, x3dna_obj, msa_objs, potentials,  thresholds, radii, radius, fragment_restrict, binding_restrict, dummy_dir,methylation)
        profile.set_energy_score_by_position(i,scr)
        if save:
           output = output_file+".scr_"+str(i+1)+".txt"
           protein_name=os.path.basename(pdb_file)
           scr.write(output,protein_name=protein_name,normal=True,overwrite=True)
        # Flush
        sys.stdout.flush()

    return profile

 except Exception as e:
    
    raise Exception("Failed profile %s"%e)
    



##############################################################################################################################
#                                                                                                                               
# Fast calculation of scores:
# Calculate the profile by threading the DNA fragment sequences on the PDB native structure to calculate the potentials
#
##############################################################################################################################

def calculate_single_profile_by_thread(dna_file,thresholds,energy_profile,pdb_file,output_file,pbm_dir,pdb_dir,families,potential_file, radius,fragment_restrict, binding_restrict, split_potential,auto_mode,family_potentials,pbm_potentials,score_threshold,taylor_approach,pmf,bins,known,meme,reset, dummy_dir,verbose,save,methylation=False):

 try:
#Get the output names
    output_dir  = os.path.dirname(output_file)
    output_name = os.path.basename(output_file)
    if verbose: sys.stdout.write("Output files %s (folder %s)\n"%(output_name,output_dir))

#Get the DNA sequence    
    use_dna=True
    for header, sequence in functions.parse_fasta_file(dna_file):
        if verbose: sys.stdout.write("Read DNA named %s (lenght %d) \n"%(header,len(sequence)))
        if use_dna:
           dna_sequence = sequence
           dna_name = header.split()[0].replace("/","-")
           use_dna=False
    if save:
        pdb_dna_file=os.path.join(output_dir,dna_name+".pdb")
        build_BDNA(dna_sequence, pdb_dna_file)

#Select the format of PWM when methylation
    if meme: methylation_pwm = False
    else:    methylation_pwm = methylation
    
#Get the PWM
    if not functions.fileExist(output_file+".meme"):
           if verbose:sys.stdout.write("Generate PWM %s ...\n"%(output_file))
           PWM.get_single_pwm(pdb_file,None,False,output_file,pbm_dir,pdb_dir,families,potential_file, radius,fragment_restrict, binding_restrict, split_potential,auto_mode,family_potentials,pbm_potentials,score_threshold,taylor_approach,pmf,bins,known,meme,reset, dummy_dir,verbose,methylation_pwm)

    #Scan the DNA sequence with the PWM to get FIMO scores
    if verbose: sys.stdout.write("Get FIMO scores\n")
    pwm_file=output_file+".meme"
    max_stored_matches = None
    profile=ProfileProtein(dna_sequence,thresholds)
    #Modify the DNA sequence to standard only nucleotides 
    if methylation_pwm:
       dummy_dna_file = dna_file
    else:
       if not os.path.exists( os.path.join(dummy_dir,"DNA_FILES") ): os.makedirs(os.path.join(dummy_dir,"DNA_FILES"))
       dummy_dna_file = os.path.join(dummy_dir,"DNA_FILES",os.path.basename(dna_file))
       with open(dummy_dna_file,"w") as fo:
         for header_dna, sequence_dna in functions.parse_fasta_file(dna_file):
             sequence_dna_modified = sequence_dna.upper().upper().replace("X","C").replace("O","C").replace("J","G").replace("Q","G")
             fo.write(">%s\n%s\n"%(header_dna,sequence_dna_modified))
       fo.close()
    for threshold in thresholds:
        if verbose: sys.stdout.write("\t-- Threshold %s\n"%threshold)
        fimo_obj = fimo.get_fimo_obj(pwm_file, dummy_dna_file, float(threshold), max_stored_matches, dummy_dir)
        profile.set_fimo_scores_by_threshold(fimo_obj,threshold)
    

# Get the properties of the structure
    # 1) Get PDB object #
    if verbose: sys.stdout.write("Reading PDB file %s ...\n"%(os.path.basename(pdb_file)))
    pdb_obj = PDB(pdb_file)
    # 2) Get DSSP object #
    if verbose: sys.stdout.write("\t-- Get DSSP ...\n")
    dssp_obj = dssp.get_dssp_obj(pdb_file)
    # 3) Get X3DNA object #
    if verbose: sys.stdout.write("\t-- Get DNA object ...\n")
    x3dna_obj = x3dna.get_x3dna_obj(pdb_file, dummy_dir)
    # 4) Get contacts object #
    if verbose: sys.stdout.write("\t-- Get contacts ...\n")
    contacts_obj = contacts.get_contacts_obj(pdb_obj, x3dna_obj)
    #    Skip if no contacts #
    if len(contacts_obj._contacts) == 0: raise ValueError("No protein-DNA contacts found!")
    # 5) Get interface object #
    if verbose: sys.stdout.write("\t-- Get interface ...\n")
    interface_obj = interface.get_interface_obj(pdb_obj, x3dna_obj, contacts_obj)
    # 6) Get triads object #
    if verbose: sys.stdout.write("\t-- Get triads ...\n")
    triads_obj = triads.get_triads_obj(pdb_obj, dssp_obj, x3dna_obj, contacts_obj)

# Flush
    sys.stdout.flush()

# Load statistical potential
    if radius <= 0: radius=float(config.get("Parameters", "max_contact_distance"))
    if verbose: sys.stdout.write("\t-- Load potentials ...\n")
    potentials, thresholds, radii, structural_homologs_by_chain = PWM.load_statistical_potentials(pdb_obj, pdb_dir, pbm_dir, families, radius, potential_file, energy_profile, auto_mode, family_potentials, pbm_potentials, score_threshold, taylor_approach, pmf , bins, known, None, dummy_dir,verbose)

# Get the dictionary of msa_objs
    msa_objs={}
    split_potentials=["3d","3dc","local","pair","s3dc","s3dc_di","s3dc_dd"]
    for sp in split_potentials:
          if verbose: sys.stdout.write("\t-- Get MSA for potential %s ...\n"%sp)
          msa_obj = PWM.get_msa_obj(triads_obj, x3dna_obj, potentials, radii, fragment_restrict, binding_restrict, sp, thresholds,methylation_pwm)
          msa_objs.setdefault(sp,msa_obj)


# Flush
    sys.stdout.flush()

# Split DNA sequence according to interface
    # All TF-DNA models generated with model_dna are restricted to the DNA interface
    # Get the interface length
    interface_length=interface_obj.get_interface_length()
    if verbose: 
        sys.stdout.write("Interface by PDB\n")
        sys.stdout.write("\t-- Length %d\n"%interface_length)
        sys.stdout.write("\t-- Starting basepairs in structure %d\n"%interface_obj.get_interface_start())
        sys.stdout.write("\t-- Last basepairs in structure %d\n"%interface_obj.get_interface_end())
    if len(dna_sequence)<interface_length:
        sys.stdout.write("Stop. Length of DNA is shorter than interface\n")
        raise ValueError("Error in PROFILE thread: Length of DNA is shorter than interface\n")

    # Get the energy-scores by windows along the DNA sequence
    for i in range(len(dna_sequence)-interface_length+1):

        # Set the DNA fragment
        dna_fragment=dna_sequence[i:i+interface_length]
        if verbose: sys.stdout.write("Score %d: fragment %s\n"%(i+1,dna_fragment))
        pdb_dna_obj = model_dna.get_dna_model_pdb_obj(pdb_file, dna_fragment, x3dna_obj, interface_obj, interface_range=None, dummy_dir=dummy_dir)

        # Create PDB with DNA fragment
        if save: pdb_dna_obj.write(os.path.join(output_dir, output_name + "." + str(i+1) + "." + dna_fragment + ".pdb"),force=True)

        # Thread DNA fragment to get new triads 
        if verbose: sys.stdout.write("\t-- Thread DNA fragment on triads ...\n")
        fragment_triads_obj = threading_to_triads.thread_kmer_triads(dna_fragment,triads_obj)

        # Calculate the energy scores and fill the profile#
        if verbose: sys.stdout.write("\t-- Calculate scores with potential: %s ...\n"%energy_profile)
        if energy_profile is None or energy_profile=="all": scr=scorer.score()
        else: scr=scorer.score(potential=energy_profile)
        scr.calculate_energies_and_binding(fragment_triads_obj, x3dna_obj, msa_objs, potentials,  thresholds, radii, radius, fragment_restrict, binding_restrict, dummy_dir,methylation)
        profile.set_energy_score_by_position(i,scr)
        if save:
           output = output_file+".scr_"+str(i+1)+".txt"
           protein_name=os.path.basename(pdb_file)
           scr.write(output,protein_name=protein_name,normal=True,overwrite=True)
        # Flush
        sys.stdout.flush()


    return profile

 except Exception as e:
    
    raise Exception("Failed profile %s"%e)
    

#######################################################################################################################################
#
# Accurate calculation of scores:
# Calculate the profile by modeling the DNA fragment sequences using the native structure as template to calculate the potentials
#
#######################################################################################################################################

def calculate_single_profile_by_models(dna_file,thresholds,energy_profile,pdb_file,output_file,pbm_dir,pdb_dir,families,potential_file, radius,fragment_restrict, binding_restrict, split_potential,auto_mode,family_potentials,pbm_potentials,score_threshold,taylor_approach,pmf,bins,known,meme,reset, dummy_dir,verbose,save,methylation=False):

 try:   

#Get the output names
    output_dir  = os.path.dirname(output_file)
    output_name = os.path.basename(output_file)
    if verbose: sys.stdout.write("Output files %s (folder %s)"%(output_name,output_dir))

#Get the DNA sequence    
    use_dna=True
    for header, sequence in functions.parse_fasta_file(dna_file):
        if verbose: sys.stdout.write("Read DNA named %s (lenght %d) \n"%(header,len(sequence)))
        if use_dna:
           dna_sequence = sequence
           dna_name = header.split()[0].replace("/","-")
           use_dna=False
    if save:
        pdb_dna_file=os.path.join(output_dir,dna_name+".pdb")
        build_BDNA(dna_sequence, pdb_dna_file)

#Select the format of PWM when methylation
    if meme: methylation_pwm = False
    else:    methylation_pwm = methylation

#Get the PWM
    if not functions.fileExist(output_file+".meme"):
         if verbose:sys.stdout.write("Generate PWM %s ...\n"%(output_file))
         PWM.get_single_pwm(pdb_file,None,False,output_file,pbm_dir,pdb_dir,families,potential_file, radius,fragment_restrict, binding_restrict, split_potential,auto_mode,family_potentials,pbm_potentials,score_threshold,taylor_approach,pmf,bins,known,meme,reset, dummy_dir,verbose,methylation_pwm)

    #Scan the DNA sequence with the PWM to get FIMO scores
    if verbose: sys.stdout.write("Get FIMO scores\n")
    pwm_file=output_file+".meme"
    max_stored_matches = None
    profile=ProfileProtein(dna_sequence,thresholds)
    #Modify the DNA sequence to standard only nucleotides 
    if methylation_pwm:
       dummy_dna_file = dna_file
    else:
       if not os.path.exists( os.path.join(dummy_dir,"DNA_FILES") ): os.makedirs(os.path.join(dummy_dir,"DNA_FILES"))
       dummy_dna_file = os.path.join(dummy_dir,"DNA_FILES",os.path.basename(dna_file))
       with open(dummy_dna_file,"w") as fo:
         for header_dna, sequence_dna in functions.parse_fasta_file(dna_file):
             sequence_dna_modified = sequence_dna.upper().upper().replace("X","C").replace("O","C").replace("J","G").replace("Q","G")
             fo.write(">%s\n%s\n"%(header_dna,sequence_dna_modified))
       fo.close()
    for threshold in thresholds:
        if verbose: sys.stdout.write("\t-- Threshold %s\n"%threshold)
        fimo_obj = fimo.get_fimo_obj(pwm_file, dummy_dna_file, float(threshold), max_stored_matches, dummy_dir)
        profile.set_fimo_scores_by_threshold(fimo_obj,threshold)
    
# Get the properties of the structure
    # 1) Get PDB object #
    if verbose: sys.stdout.write("Reading PDB file %s ...\n"%(os.path.basename(pdb_file)))
    pdb_obj = PDB(pdb_file)
    # 2) Get DSSP object #
    if verbose: sys.stdout.write("\t-- Get DSSP ...\n")
    dssp_obj = dssp.get_dssp_obj(pdb_file)
    # 3) Get X3DNA object #
    if verbose: sys.stdout.write("\t-- Get DNA object ...\n")
    x3dna_obj = x3dna.get_x3dna_obj(pdb_file, dummy_dir)
    # 4) Get contacts object #
    if verbose: sys.stdout.write("\t-- Get contacts ...\n")
    contacts_obj = contacts.get_contacts_obj(pdb_obj, x3dna_obj)
    #    Skip if no contacts #
    if len(contacts_obj._contacts) == 0: raise ValueError("No protein-DNA contacts found!")
    # 5) Get interface object #
    if verbose: sys.stdout.write("\t-- Get interface ...\n")
    interface_obj = interface.get_interface_obj(pdb_obj, x3dna_obj, contacts_obj)
    # 6) Get triads object #
    if verbose: sys.stdout.write("\t-- Get triads ...\n")
    triads_obj = triads.get_triads_obj(pdb_obj, dssp_obj, x3dna_obj, contacts_obj)

# Flush
    sys.stdout.flush()



# Load statistical potential
    if radius <= 0: radius=float(config.get("Parameters", "max_contact_distance"))
    if verbose: sys.stdout.write("\t-- Load potentials ...\n")
    potentials, thresholds, radii, structural_homologs_by_chain = PWM.load_statistical_potentials(pdb_obj, pdb_dir, pbm_dir, families, radius, potential_file, energy_profile, auto_mode, family_potentials, pbm_potentials, score_threshold, taylor_approach, pmf , bins, known, None, dummy_dir,verbose)

# Get the dictionary of msa_objs
    msa_objs={}
    split_potentials=["3d","3dc","local","pair","s3dc","s3dc_di","s3dc_dd"]
    for sp in split_potentials:
          if verbose: sys.stdout.write("\t-- Get MSA for potential %s ...\n"%sp)
          msa_obj = PWM.get_msa_obj(triads_obj, x3dna_obj, potentials, radii, fragment_restrict, binding_restrict, sp, thresholds,methylation_pwm)
          msa_objs.setdefault(sp,msa_obj)

# Flush
    sys.stdout.flush()


# Split DNA sequence according to interface
    # All TF-DNA models generated with model_dna are restricted to the DNA interface
    # Get the interface length
    interface_length=interface_obj.get_interface_length()
    if verbose: 
        sys.stdout.write("Interface by PDB\n")
        sys.stdout.write("\t-- Length %d\n"%interface_length)
        sys.stdout.write("\t-- Starting basepairs in structure %d\n"%interface_obj.get_interface_start())
        sys.stdout.write("\t-- Last basepairs in structure %d\n"%interface_obj.get_interface_end())
    if len(dna_sequence)<interface_length:
        sys.stdout.write("Stop. Length of DNA is shorter than interface\n")
        raise ValueError("Error in PROFILE: Length of DNA is shorter than interface\n")

    dummy_pdb_dir = os.path.join(dummy_dir,os.path.basename(output_dir))
    if not os.path.exists(dummy_pdb_dir): os.makedirs(dummy_pdb_dir)

      # Get the energy-scores by windows along the DNA sequence
    for i in range(len(dna_sequence)-interface_length+1):

        # Set the DNA fragment
        dna_fragment=dna_sequence[i:i+interface_length]
        if verbose: sys.stdout.write("Score %d: fragment %s\n"%(i+1,dna_fragment))
        pdb_dna_obj = model_dna.get_dna_model_pdb_obj(pdb_file, dna_fragment, x3dna_obj, interface_obj, interface_range=None, dummy_dir=dummy_dir)

        # Create PDB with DNA fragment
        if save: pdb_dna_obj.write(os.path.join(output_dir, output_name + "." + str(i+1) + "." + dna_fragment + ".pdb"),force=True)
        pdb_frag_file = os.path.join(dummy_pdb_dir, output_name + "." + str(i+1) + "." + dna_fragment + ".pdb")
        pdb_dna_obj.write(pdb_frag_file,force=True)
        # Get X3DNA object fragment#
        if verbose: sys.stdout.write("\t-- Get DNA-fragment object ...\n")
        x3dna_dna_obj = x3dna.get_x3dna_obj(pdb_frag_file, dummy_dir)
        # Get contacts object #
        if verbose: sys.stdout.write("\t-- Get contacts on fragment ...\n")
        contacts_dna_obj = contacts.get_contacts_obj(pdb_dna_obj, x3dna_dna_obj)
        # Skip if no contacts #
        if len(contacts_dna_obj._contacts) == 0: raise ValueError("No protein-DNA contacts found!")
        # Get interface object #
        if verbose: sys.stdout.write("\t-- Get interface on fragment ...\n")
        interface_dna_obj = interface.get_interface_obj(pdb_dna_obj, x3dna_dna_obj, contacts_dna_obj)
        # Get triads object #
        if verbose: sys.stdout.write("\t-- Get triads on fragment ...\n")
        triads_obj = triads.get_triads_obj(pdb_dna_obj, dssp_obj, x3dna_dna_obj, contacts_dna_obj)

        # Calculate the energy scores and fill the profile#
        if verbose: sys.stdout.write("\t-- Calculate scores on fragment with potential: %s...\n"%energy_profile)
        if energy_profile is None or energy_profile=="all": scr=scorer.score()
        else: scr=scorer.score(potential=energy_profile)
        scr.calculate_energies_and_binding(triads_obj, x3dna_dna_obj, msa_objs, potentials,  thresholds, radii, radius, fragment_restrict, binding_restrict, dummy_dir,methylation)
        profile.set_energy_score_by_position(i,scr)
        if save:
           output = output_file+".scr_"+str(i+1)+".txt"
           protein_name=os.path.basename(pdb_file)
           scr.write(output,protein_name=protein_name,normal=True,overwrite=True)

        # Flush
        sys.stdout.flush()


    return profile

 except Exception as e:
    
     raise Exception("Failed profile: %s"%e)
    
#-------------#
# Classes     #
#-------------#


class Profile(object):
    """
    This class define a PROFILE object of ProteinProfiles
    """

    def __init__(self,dna,thresholds=["0.05", "0.001"], potential=None):
        self._dna             = dna
        self._thresholds      = thresholds[:]
        self._proteins        = []
        self._split_potentials=["3d","3dc","local","pair","s3dc","s3dc_di","s3dc_dd"]
        self._score_types     =["fimo_binding","fimo_score","fimo_log_score","energy","energy_best","energy_per_nucleotide"]
        if potential is not None: self._potential=potential 

    def set_dna(self,dna):
        self._dna=dna

    def get_dna(self):
        return self._dna

    def set_thresholds(self,thresholds):
        self._thresholds=thresholds[:]

    def get_thresholds(self):
        return self._thresholds

    def set_split_potentials(self,split_potentials):
        self._split_potentials=split_potentials[:]

    def get_split_potentials(self):
        return self._split_potentials

    def set_score_types(self,score_types):
        self._score_types=score_types[:]

    def get_score_types(self):
        return self._score_types

    def get_potential(self):
        return self._potential

    def add(self,profile_protein):
        skip = False
        if profile_protein.get_dna() != self._dna: skip = True
        st=set(self._thresholds)
        sp=set(profile_protein.get_thresholds())
        if not skip and st != st.intersection(sp): skip=True
        if skip:
           print("Protein profile is incompatible with the general profile")
        else:
           self._proteins.append(profile_protein)

    def get_profiles(self):
        return self._proteins

    def get_profile(self,i):
        return self._proteins[i]
    

    def get_fimo_binding_by_threshold(self,threshold):
        if threshold not in self._thresholds: return None
        if len(self.get_profiles())<=0: return None
        dim=0
        for profile in self.get_profiles():
            if profile.get_fimo_binding_by_threshold(threshold) is None: continue
            if dim>0: x=np.vstack((x,profile.get_fimo_binding_by_threshold(threshold)))
            else:     x=profile.get_fimo_binding_by_threshold(threshold)
            dim=dim+1
        return x

    def get_fimo_score_by_threshold(self,threshold):
        if threshold not in self._thresholds: return None
        if len(self.get_profiles())<=0: return None
        dim=0
        for profile in self.get_profiles():
            if profile.get_fimo_score_by_threshold(threshold) is None: continue
            if dim>0: x=np.vstack((x,profile.get_fimo_score_by_threshold(threshold)))
            else:     x=profile.get_fimo_score_by_threshold(threshold)
            dim=dim+1
        return x

    def get_fimo_log_score_by_threshold(self,threshold):
        if threshold not in self._thresholds: return None
        if len(self.get_profiles())<=0: return None
        dim=0
        for profile in self.get_profiles():
            if profile.get_fimo_log_score_by_threshold(threshold) is None: continue
            if dim>0: x=np.vstack((x,profile.get_fimo_log_score_by_threshold(threshold)))
            else:     x=profile.get_fimo_log_score_by_threshold(threshold)
            dim=dim+1
        return x

    def get_energy_profile(self,normal=False,potential="s3dc_dd"):
        dim=0
        if len(self.get_profiles())<=0: return None
        for profile in self.get_profiles():
            if dim>0: x=np.vstack((x,profile.get_energy_profile(normal,potential)))
            else:     x=profile.get_energy_profile(normal,potential)
            dim=dim+1
        return x

    def get_energy_weighted_profile(self,threshold,normal=False,potential="s3dc_dd"):
        dim=0
        if len(self.get_profiles())<=0: return None
        for profile in self.get_profiles():
            if dim>0: x=np.vstack((x,profile.get_energy_weighted_profile(threshold,normal,potential)))
            else:     x=profile.get_energy_weighted_profile(threshold,normal,potential)
            dim=dim+1
        return x

    def get_energy_best_profile(self,potential="s3dc_dd"): 
        dim=0
        if len(self.get_profiles())<=0: return None
        for profile in self.get_profiles():
            if dim>0: x=np.vstack((x,profile.get_energy_best_profile(potential)))
            else:     x=profile.get_energy_best_profile(potential)
            dim=dim+1
        return x

    def get_energy_best_weighted_profile(self,threshold,potential="s3dc_dd"): 
        dim=0
        if len(self.get_profiles())<=0: return None
        for profile in self.get_profiles():
            if dim>0: x=np.vstack((x,profile.get_energy_best_weighted_profile(threshold,potential)))
            else:     x=profile.get_energy_best_weighted_profile(threshold,potential)
            dim=dim+1
        return x

    def get_energy_per_nucleotide_profile(self,normal=False,potential="s3dc_dd"): 
        dim=0
        if len(self.get_profiles())<=0: return None
        for profile in self.get_profiles():
            if dim>0: x=np.vstack((x,profile.get_energy_per_nucleotide_profile(normal,potential)))
            else:     x=profile.get_energy_per_nucleotide_profile(normal,potential)
            dim=dim+1
        return x

    def get_energy_per_nucleotide_weighted_profile(self,threshold,normal=False,potential="s3dc_dd"): 
        dim=0
        if len(self.get_profiles())<=0: return None
        for profile in self.get_profiles():
            if dim>0: x=np.vstack((x,profile.get_energy_per_nucleotide_weighted_profile(threshold,normal,potential)))
            else:     x=profile.get_energy_per_nucleotide_weighted_profile(threshold,normal,potential)
            dim=dim+1
        return x

    def get_profile_score(self,threshold=None,score_type="energy",normal=False,potential="s3dc_dd"):
        """
          score_type = 
                          energy, 
                          energy_best, 
                          energy_per_nucleotide,
                          fimo_binding
                          fimo_score
                          fimo_log_score
          if threshold is not None energies are weighted
        """
        dim=0
        if len(self.get_profiles())<=0: return None
        for profile in self.get_profiles():
            if dim>0: x=np.vstack((x,profile.get_profile_score(threshold,score_type,normal,potential)))
            else:     x=profile.get_profile_score(threshold,score_type,normal,potential)
            dim=dim+1
        if dim==1: x=np.vstack((x,x))
        return x


    def to_mean_table(self):
        split_potentials=self.get_split_potentials()
        thresholds=self.get_thresholds()
        tbl_dict={}
        tbl_dict.setdefault("Position",[int(x)+1 for x in range(len(self._dna))])
        tbl_dict.setdefault("Nucleotide",list(self._dna))
        for thr in thresholds:
         if thr is not None:
            tbl_dict.setdefault("binding_"+thr,self.get_fimo_binding_by_threshold(thr).mean(0))
            tbl_dict.setdefault("fimo_score_"+thr,self.get_fimo_score_by_threshold(thr).mean(0))
            tbl_dict.setdefault("fimo_log_score_"+thr,self.get_fimo_log_score_by_threshold(thr).mean(0))
        for potential in split_potentials:
            y=self.get_energy_profile(normal=False,potential=potential).mean(0)
            if len(y[abs(y)>1.0e-10])<1: continue
            tbl_dict.setdefault(potential,self.get_energy_profile(normal=False,potential=potential).mean(0))
            tbl_dict.setdefault("normal_"+potential,self.get_energy_profile(normal=True,potential=potential).mean(0))
            tbl_dict.setdefault(potential+"_x_nucleotide",self.get_energy_per_nucleotide_profile(normal=False,potential=potential).mean(0))
            tbl_dict.setdefault("normal_"+potential+"_x_nucleotide",self.get_energy_per_nucleotide_profile(normal=True,potential=potential).mean(0))
            tbl_dict.setdefault("best_"+potential,self.get_energy_best_profile(potential=potential).mean(0))
        for thr in thresholds:
          if thr is None: continue
          for potential in split_potentials:
            y=self.get_energy_weighted_profile(thr,normal=False,potential=potential).mean(0)
            if len(y[abs(y)>1.0e-10])<1: continue
            tbl_dict.setdefault(potential+"_weighted_"+thr,self.get_energy_weighted_profile(thr,normal=False,potential=potential).mean(0))
            tbl_dict.setdefault("normal_"+potential+"_weighted_"+thr,self.get_energy_weighted_profile(thr,normal=True,potential=potential).mean(0))
            tbl_dict.setdefault(potential+"_x_nucleotide"+"_weighted_"+thr,self.get_energy_per_nucleotide_weighted_profile(thr,normal=False,potential=potential).mean(0))
            tbl_dict.setdefault("normal_"+potential+"_x_nucleotide"+"_weighted_"+thr,self.get_energy_per_nucleotide_weighted_profile(thr,normal=True,potential=potential).mean(0))
            tbl_dict.setdefault("best_"+potential+"_weighted_"+thr,self.get_energy_best_weighted_profile(thr,potential=potential).mean(0))
        table=pd.DataFrame(tbl_dict)
        return table

    def to_rmsd_table(self):
        split_potentials=self._split_potentials
        thresholds=self.get_thresholds()
        tbl_dict={}
        tbl_dict.setdefault("Position",[int(x)+1 for x in range(len(self._dna))])
        tbl_dict.setdefault("Nucleotide",list(self._dna))
        for thr in thresholds:
         if thr is not None:
            tbl_dict.setdefault("binding_"+thr,self.get_fimo_binding_by_threshold(thr).std(0))
            tbl_dict.setdefault("fimo_score_"+thr,self.get_fimo_score_by_threshold(thr).std(0))
            tbl_dict.setdefault("fimo_log_score_"+thr,self.get_fimo_log_score_by_threshold(thr).std(0))
        for potential in split_potentials:
            y=self.get_energy_profile(normal=False,potential=potential).mean(0)
            if len(y[abs(y)>1.0e-10])<1: continue
            tbl_dict.setdefault(potential,self.get_energy_profile(normal=False,potential=potential).std(0))
            tbl_dict.setdefault("normal_"+potential,self.get_energy_profile(normal=True,potential=potential).std(0))
            tbl_dict.setdefault(potential+"_x_nucleotide",self.get_energy_per_nucleotide_profile(normal=False,potential=potential).std(0))
            tbl_dict.setdefault("normal_"+potential+"_x_nucleotide",self.get_energy_per_nucleotide_profile(normal=True,potential=potential).std(0))
            tbl_dict.setdefault("best_"+potential,self.get_energy_best_profile(potential=potential).std(0))
        for thr in thresholds:
          if thr is None: continue
          for potential in split_potentials:
            y=self.get_energy_weighted_profile(thr,normal=False,potential=potential).mean(0)
            if len(y[abs(y)>1.0e-10])<1: continue
            tbl_dict.setdefault(potential+"_weighted_"+thr,self.get_energy_weighted_profile(thr,normal=False,potential=potential).std(0))
            tbl_dict.setdefault("normal_"+potential+"_weighted_"+thr,self.get_energy_weighted_profile(thr,normal=True,potential=potential).std(0))
            tbl_dict.setdefault(potential+"_x_nucleotide"+"_weighted_"+thr,self.get_energy_per_nucleotide_weighted_profile(thr,normal=False,potential=potential).std(0))
            tbl_dict.setdefault("normal_"+potential+"_x_nucleotide"+"_weighted_"+thr,self.get_energy_per_nucleotide_weighted_profile(thr,normal=True,potential=potential).std(0))
            tbl_dict.setdefault("best_"+potential+"_weighted_"+thr,self.get_energy_best_weighted_profile(thr,potential=potential).std(0))
        table=pd.DataFrame(tbl_dict)
        return table

    def write_mean_table(self,output):
        table=self.to_mean_table()
        table.to_csv(output)

    def write_rmsd_table(self,output):
        table=self.to_rmsd_table()
        table.to_csv(output)

    def plot(self,output):
        table=self.to_mean_table()
        error=self.to_rmsd_table()
        columns=table.columns.values.tolist()
        #x = table["Nucleotide"].values
        #x = np.array(zip(table["Position"].values.tolist(),table["Nucleotide"].values.tolist()))
        x = table["Position"].values
        for column in columns:
          if column == "Nucleotide": continue
          if column == "Position": continue
          out   = output+"_"+column+".png"
          title = "Graph of profile %s"%column
          y     = table[column].values
          e     = error[column].values
          if len(y[abs(y)>1.0e-10])<1: continue
          plt.plot(x,y,'-', ms=5, lw=2, alpha=0.7, mfc='cyan')
          plt.xlabel('nucleotide ')
          plt.ylabel('score profile')
          plt.fill_between(x,y-e,y+e,alpha=0.2)
          plt.savefig(out,format="png")
          plt.close()
    
    def plot_all(self,output):
        x=np.array([int(i)+1 for i in range(len(self._dna))])
        set_columns=set()
        for profile in self.get_profiles():
            table=profile.to_table()
            set_columns.update(set(table.columns.values.tolist()))
        columns=list(set_columns)
        for column in columns:
          if column == "Nucleotide": continue
          if column == "Position": continue
          out   = output+"_"+column+".png"
          title = "Graph of profile %s"%column
          for profile in self.get_profiles():
                table=profile.to_table()
                if column not in set(table.columns.values.tolist()):continue
                y = table[column].values
                if len(y[abs(y)>1.0e-10])<1: continue
                plt.plot(x,y,'-', ms=5, lw=2, alpha=0.7)
                plt.xlabel('nucleotide ')
                plt.ylabel('score profile')
                plt.savefig(out,format="png")
          plt.close()
         

    def plot_html(self,score_types,normal,energies,output,header_title=False):
        x=np.array([int(i)+1 for i in range(len(self._dna))])
        name_profile = os.path.basename(output)
        html_file    = output+".html"
        thresholds=[]
        thresholds.append(None)
        thresholds.extend(self.get_thresholds())
        '''
        COLOR LIST CSS
                aliceblue, antiquewhite, aqua, aquamarine, azure,
                beige, bisque, black, blanchedalmond, blue,
                blueviolet, brown, burlywood, cadetblue,
                chartreuse, chocolate, coral, cornflowerblue,
                cornsilk, crimson, cyan, darkblue, darkcyan,
                darkgoldenrod, darkgray, darkgrey, darkgreen,
                darkkhaki, darkmagenta, darkolivegreen, darkorange,
                darkorchid, darkred, darksalmon, darkseagreen,
                darkslateblue, darkslategray, darkslategrey,
                darkturquoise, darkviolet, deeppink, deepskyblue,
                dimgray, dimgrey, dodgerblue, firebrick,
                floralwhite, forestgreen, fuchsia, gainsboro,
                ghostwhite, gold, goldenrod, gray, grey, green,
                greenyellow, honeydew, hotpink, indianred, indigo,
                ivory, khaki, lavender, lavenderblush, lawngreen,
                lemonchiffon, lightblue, lightcoral, lightcyan,
                lightgoldenrodyellow, lightgray, lightgrey,
                lightgreen, lightpink, lightsalmon, lightseagreen,
                lightskyblue, lightslategray, lightslategrey,
                lightsteelblue, lightyellow, lime, limegreen,
                linen, magenta, maroon, mediumaquamarine,
                mediumblue, mediumorchid, mediumpurple,
                mediumseagreen, mediumslateblue, mediumspringgreen,
                mediumturquoise, mediumvioletred, midnightblue,
                mintcream, mistyrose, moccasin, navajowhite, navy,
                oldlace, olive, olivedrab, orange, orangered,
                orchid, palegoldenrod, palegreen, paleturquoise,
                palevioletred, papayawhip, peachpuff, peru, pink,
                plum, powderblue, purple, red, rosybrown,
                royalblue, saddlebrown, salmon, sandybrown,
                seagreen, seashell, sienna, silver, skyblue,
                slateblue, slategray, slategrey, snow, springgreen,
                steelblue, tan, teal, thistle, tomato, turquoise,
                violet, wheat, white, whitesmoke, yellow,
                yellowgreen
        '''
        red_color_long   = ["darkred", "firebrick","tomato","red", "orangered", "coral", "orange", "goldenrod","gold","siena","yellow"]
        blue_color_long  = ["navy","darkblue","blue","midnightblue","dodgerblue","steelblue","deepskyblue","skyblue","lightblue","cyan","lightcyan"]
        green_color_long = ["darkgreen","forestgreen","green","lawngreen","springgreen","lightgreen","limegreen","turquoise","aquamarine"]
        red_color   = ["darkred", "firebrick", "red","orangered", "orange", "goldenrod","gold"]
        blue_color  = ["navy","blue","steelblue","deepskyblue","lightblue","cyan","lightcyan"]
        green_color = ["darkgreen","forestgreen","green","springgreen","lightgreen","limegreen","turquoise","aquamarine"]
        selection=[]
        select_main=set()
        trace_list=[]
        for sc in  score_types:
          name_sc=sc
          for ii in range(len(thresholds)):
             thr = thresholds[ii]
             if thr is None: name_thr=name_sc+"_1.0"
             else:           name_thr=name_sc+"_"+thr
             dna=self.get_dna()
             text_dna = np.array(["".join(dna[k:min(k+6,len(dna))])+"..." for k in range(len(dna))])
             if "fimo" in sc:
               if thr is None: continue
               select_main.add(name_sc)
               name = name_thr
               if self.get_profile_score(thr,sc) is None: continue
               y = self.get_profile_score(thr,sc).mean(0).round(2)
               e = self.get_profile_score(thr,sc).std(0).round(2)
               #upper and lower bound are not necessary, but kept for symmetry of data with e=0
               if "binding" in sc: e=np.zeros(len(y))
               if "binding" in sc and thr=="0.01":
                  #trace = go.Scatter(x = x, y = y, text= text_dna, mode = 'lines', name = thr, legendgroup = thr, line = dict(color=red_color[ii],width = 3), fill='tonexty',  hoverinfo="x+y+text", visible = True)
                  #upper_bound = go.Scatter(x = x, y = y+e, mode = 'lines', name = thr, legendgroup = thr, line = dict(color=red_color[ii],width = 0.8), fill='tonexty', showlegend = False, hoverinfo='none', visible = True)
                  #lower_bound = go.Scatter(x = x, y = y-e, mode = 'lines', name = thr, legendgroup = thr, line = dict(color=red_color[ii],width = 0.8), showlegend = False,  hoverinfo='none', visible = True)
                  trace =       dict(type="scatter", x = x, y = y, text= text_dna, mode = 'lines', name = thr, legendgroup = thr, line = dict(color=red_color[ii],width = 3), fill='tonexty',  hoverinfo="x+y+text", visible = True)
                  upper_bound = dict(type="scatter", x = x, y = y+e, mode = 'lines', name = thr, legendgroup = thr, line = dict(color=red_color[ii],width = 0.8), fill='tonexty', showlegend = False, hoverinfo='none', visible = True)
                  lower_bound = dict(type="scatter", x = x, y = y-e, mode = 'lines', name = thr, legendgroup = thr, line = dict(color=red_color[ii],width = 0.8), showlegend = False,  hoverinfo='none', visible = True)
               else:
                  visible=False
                  if "binding" in sc: visible="legendonly"
                  #trace = go.Scatter(x = x, y = y, text= text_dna, mode = 'lines', name = thr, legendgroup = thr, line = dict(color=red_color[ii],width = 3),  fill='tonexty', hoverinfo="x+y+text", visible = visible)
                  #upper_bound = go.Scatter(x = x, y = y+e, mode = 'lines', name = thr, legendgroup = thr, line = dict(color=red_color[ii],width = 0.8), fill='tonexty', showlegend = False, hoverinfo='none', visible = visible)
                  #lower_bound = go.Scatter(x = x, y = y-e, mode = 'lines', name = thr, legendgroup = thr, line = dict(color=red_color[ii],width = 0.8),  showlegend = False, hoverinfo='none', visible = visible)
                  trace =       dict(type="scatter", x = x, y = y, text= text_dna, mode = 'lines', name = thr, legendgroup = thr, line = dict(color=red_color[ii],width = 3),  fill='tonexty', hoverinfo="x+y+text", visible = visible)
                  upper_bound = dict(type="scatter", x = x, y = y+e, mode = 'lines', name = thr, legendgroup = thr, line = dict(color=red_color[ii],width = 0.8), fill='tonexty', showlegend = False, hoverinfo='none', visible = visible)
                  lower_bound = dict(type="scatter", x = x, y = y-e, mode = 'lines', name = thr, legendgroup = thr, line = dict(color=red_color[ii],width = 0.8),  showlegend = False, hoverinfo='none', visible = visible)
               trace_list.append(lower_bound)
               trace_list.append(trace)
               trace_list.append(upper_bound)
               selection.append(name)
               selection.append(name)
               selection.append(name)
             if "energy"  in sc :
               for energy in energies:
                   if normal: name_ene = "normal_"+energy
                   else:      name_ene = energy
                   name = name_ene + "_" + name_thr
                   select_main.add(name_ene+"_"+name_sc)
                   if self.get_profile_score(thr,sc,normal,energy) is None : continue
                   y = self.get_profile_score(thr,sc,normal,energy).mean(0).round(2)
                   e = self.get_profile_score(thr,sc,normal,energy).std(0).round(2)
                   visible=False
                   if thr is None: thre="1.0"
                   else: thre=thr
                   #trace = go.Scatter(x = x, y = y, text= text_dna, mode = 'lines', name = thre, legendgroup = thre, line = dict(color=blue_color[ii],width = 3),  fill='tonexty', hoverinfo="x+y+text", visible = visible)
                   #upper_bound = go.Scatter(x = x, y = y+e, mode = 'lines', name = thre, legendgroup = thre, line = dict(color=blue_color[ii],width = 0.8), fill='tonexty', showlegend = False, hoverinfo='none', visible = visible)
                   #lower_bound = go.Scatter(x = x, y = y-e, mode = 'lines', name = thre, legendgroup = thre, line = dict(color=blue_color[ii],width = 0.8),  showlegend = False, hoverinfo='none', visible = visible)
                   trace =       dict(type="scatter", x = x, y = y, text= text_dna, mode = 'lines', name = thre, legendgroup = thre, line = dict(color=blue_color[ii],width = 3),  fill='tonexty', hoverinfo="x+y+text", visible = visible)
                   upper_bound = dict(type="scatter", x = x, y = y+e, mode = 'lines', name = thre, legendgroup = thre, line = dict(color=blue_color[ii],width = 0.8), fill='tonexty', showlegend = False, hoverinfo='none', visible = visible)
                   lower_bound = dict(type="scatter", x = x, y = y-e, mode = 'lines', name = thre, legendgroup = thre, line = dict(color=blue_color[ii],width = 0.8),  showlegend = False, hoverinfo='none', visible = visible)
                   trace_list.append(lower_bound)
                   trace_list.append(trace)
                   trace_list.append(upper_bound)
                   selection.append(name)
                   selection.append(name)
                   selection.append(name)
        # Define the layout and axis names #
        config = {'linkText': "", 'scrollZoom': True}
        # Determine a dropdown with the variables that will be shown in this plot #
        fimo_selection=[s for s in select_main if "fimo" in s]
        ener_selection=[s for s in select_main if "energy" in s]
        order_selection=sorted(fimo_selection)[:]
        order_selection.extend(sorted(ener_selection))
        buttons=[]
        for plt_name in order_selection:
            features=plt_name.split("_")
            score=features[:]
            if "fimo" in plt_name: score=features[1:]
            arg_visible=[]
            if "binding" in plt_name:
              for sel in selection:
                  features_selection=sel.split("_")
                  thr_sel=features_selection[-1]
                  test_sel="_".join(features_selection[0:-1])
                  if plt_name==test_sel:
                     if thr_sel=="0.01":arg_visible.append(True)
                     else:arg_visible.append("legendonly")
                  else:
                     arg_visible.append(False)
            else:
              for sel in selection:
                  features_selection=sel.split("_")
                  thr_sel=features_selection[-1]
                  test_sel="_".join(features_selection[0:-1])
                  if plt_name==test_sel:
                     if  thr_sel=="0.5":arg_visible.append(True)
                     else:arg_visible.append("legendonly")
                  else:
                     arg_visible.append(False)
            arg_yaxis={}
            arg_yaxis.setdefault('titlefont',dict(size=12))
            if "fimo" in plt_name:
               if "log" in plt_name: score[0]="-log(P-value)"
               if "binding" in plt_name: score.append("ratio")
               label="FIMO: "+" ".join(score)
            if "energy" in plt_name:
               label="POTENTIAL: "+ " ".join(score)
            arg_yaxis.setdefault('title',label)
            arguments=[]
            arguments.append(dict(visible=arg_visible))
            arguments.append(dict(yaxis=arg_yaxis))
            buttons.append(dict(label=label,method='update',args=arguments))

        #updatemenus = [dict(buttons=buttons, direction='down', showactive=True, pad= {'r': 0, 't':  15}, x = 0.25 , y = 1.85, font = {'size': 12}, bgcolor = "white")]
        updatemenus = [dict(buttons=buttons, direction='down', showactive=True, pad= {'r': 0, 't':  15}, x = 0.25 , y = 1.95, font = {'size': 12}, bgcolor = "white")]
        annotations = [dict(text='P-value\nthresholds:', x=1.08, y=1.15, xref='paper', yref='paper', showarrow=False, font=dict(size=12), align="center")]
        xaxis_dict=dict(
                        title         ='Nucleotide positions',
                        rangeslider   =dict(visible=True, bordercolor="grey", borderwidth=3),
                        titlefont     =dict(size=18), 
                        tickformat    =',d'
                        )
        if header_title:
            layout=dict(title="PROFILE: "+name_profile,updatemenus=updatemenus,annotations=annotations,plot_bgcolor='rgba(0,0,0,0)',paper_bgcolor='rgba(0,0,0,0)',yaxis= dict(title = "Ratio", titlefont=dict(size=16)),xaxis=xaxis_dict)
        else:
            layout=dict(title="",updatemenus=updatemenus,annotations=annotations,plot_bgcolor='rgba(0,0,0,0)',paper_bgcolor='rgba(0,0,0,0)',yaxis= dict(title = "Ratio", titlefont=dict(size=16)),xaxis=xaxis_dict)
        # Create the plot #
        fig = dict(data=trace_list, layout=layout)
        pltly(fig, filename=html_file, show_link=False, auto_open=False)

            

class ProfileProtein(object):
    """
    This class define a single protein PROFILE object
    """
    def __init__(self,dna,thresholds=["0.05", "0.001"],binding=None,fimo_log_score=None,fimo_score=None, energy_scores=None, potential=None):
        self._dna             = dna
        self._thresholds      = thresholds[:]
        self._potential       = "all"
        self._split_potentials=["3d","3dc","local","pair","s3dc","s3dc_di","s3dc_dd"]
        self._score_types     =["fimo_binding","fimo_score","fimo_log_score","energy","energy_best","energy_per_nucleotide"]
        self._binding         = {}
        self._fimo_log_score  = {}
        self._fimo_score      = {}
        self._energy_scores   = {}
        if binding is not None: self._binding=binding
        if fimo_log_score is not None: self._fimo_log_score=fimo_log_score
        if fimo_score is not None: self._fimo_score=fimo_score
        if energy_scores is not None: self._energy_scores=energy_scores
        if potential is not None: self._potential=potential  

    def set_dna(self,dna):
        self._dna=dna

    def get_dna(self):
        return self._dna

    def set_thresholds(self,thresholds):
        self._thresholds=thresholds[:]

    def get_thresholds(self):
        return self._thresholds

    def set_split_potentials(self,split_potentials):
        self._split_potentials=split_potentials[:]

    def get_split_potentials(self):
        return self._split_potentials

    def set_score_types(self,score_types):
        self._score_types=score_types[:]

    def get_score_types(self):
        return self._score_types

    def get_potential(self):
        return self._potential

    def set_fimo_scores_by_threshold(self, fimo_obj, threshold):
        binding_score  = np.zeros(len(self._dna))
        fimo_score     = np.zeros(len(self._dna))
        fimo_log_score = np.zeros(len(self._dna))
        ones           = np.zeros(len(self._dna))
        ones           = ones + 1
        count          = np.zeros(len(self._dna))
        for fimo_hit in fimo_obj.get_hits():
            start = fimo_hit.get_start()
            end   = fimo_hit.get_end()
            score = fimo_hit.get_score()
            logscr= -np.log(fimo_hit.get_p_value())
            binding_score[start-1:end]   =1
            fimo_score[start-1:end]     += score
            fimo_log_score[start-1:end] += logscr
            count[start-1:end]          += ones[start-1:end]
        fimo_log_score=fimo_log_score/count
        fimo_log_score[np.isnan(fimo_log_score)]=0
        fimo_log_score[np.isinf(fimo_log_score)]=0
        fimo_score=fimo_log_score/count
        fimo_score[np.isnan(fimo_score)]=0
        fimo_score[np.isinf(fimo_score)]=0
        if self._binding.has_key(threshold):        self._binding[threshold]=binding_score
        else:                                       self._binding.setdefault(threshold,binding_score)
        if self._fimo_log_score.has_key(threshold): self._fimo_log_score[threshold]=fimo_log_score
        else:                                       self._fimo_log_score.setdefault(threshold,fimo_log_score)
        if self._fimo_score.has_key(threshold):     self._fimo_score[threshold]=fimo_score
        else:                                       self._fimo_score.setdefault(threshold,fimo_score)

    def get_fimo_binding_by_threshold(self,threshold):
        if threshold not in self._thresholds: return None
        if not self._binding.has_key(threshold): return None
        return self._binding.get(threshold)

    def get_fimo_score_by_threshold(self,threshold):
        if threshold not in self._thresholds: return None
        if not self._fimo_score.has_key(threshold): return None
        return self._fimo_score.get(threshold)

    def get_fimo_log_score_by_threshold(self,threshold):
        if threshold not in self._thresholds: return None
        if not self._fimo_log_score.has_key(threshold): return None
        return self._fimo_log_score.get(threshold)

    def set_energy_score_by_position(self,position,score):
        self._energy_scores.setdefault(position,score)

    def get_energy_score_by_position(self,position):
        if position not in range(len(self._dna)): return None
        if not self._energy_scores.has_key(position): return None
        return self._energy_scores.get(position)

    def get_energy_profile(self,normal=False,potential="s3dc_dd"):
        energy = np.zeros(len(self._dna))
        ones   = np.zeros(len(self._dna))
        ones   = ones + 1
        count  = np.zeros(len(self._dna))
        scr=self.get_energy_score_by_position(0)
        binding_length=len(scr.get_binding_site())
        for i in range(len(self._dna)-binding_length+1):
            scr=self.get_energy_score_by_position(i)
            if scr is None: continue
            score_per_nucleotide=scr.get_score_per_nucleotide(normal,potential)
            if len(score_per_nucleotide.keys())<=0:continue
            start=i+min([int(nuc) for nuc in sorted(score_per_nucleotide.keys())])-1
            end  =i+max([int(nuc) for nuc in sorted(score_per_nucleotide.keys())])
            if end>len(self._dna): end=len(self._dna)
            energy[start:end] += scr.get_score(normal,potential)
            count[start:end]  += ones[start:end]
        energy = energy/count
        energy[np.isnan(energy)]=0
        energy[np.isinf(energy)]=0
        return energy

    def get_energy_weighted_profile(self,threshold,normal=False,potential="s3dc_dd"):
        if self.get_fimo_binding_by_threshold(threshold) is None: return None
        binding_score=self.get_fimo_binding_by_threshold(threshold)
        energy=np.zeros(len(self._dna))
        ones   = np.zeros(len(self._dna))
        ones   = ones + 1
        count  = np.zeros(len(self._dna))
        scr=self.get_energy_score_by_position(0)
        binding_length=len(scr.get_binding_site())
        for i in range(len(self._dna)-binding_length+1):
            scr=self.get_energy_score_by_position(i)
            if scr is None:  continue
            score_per_nucleotide=scr.get_score_per_nucleotide(normal,potential)
            if len(score_per_nucleotide.keys())<=0:continue
            start=i+min([int(nuc) for nuc in sorted(score_per_nucleotide.keys())])-1
            end  =i+max([int(nuc) for nuc in sorted(score_per_nucleotide.keys())])
            if end>len(self._dna): end=len(self._dna)
            #if one nucleotide of the interval is not part of a binding the score is skipped
            #otherwise, if all nucleotide positions are in binding add the score
            energy[start:end] += scr.get_score(normal,potential) * binding_score[start:end].prod()
            count[start:end]  += ones[start:end] * binding_score[start:end].prod()
        energy = energy/count
        energy[np.isnan(energy)]=0
        energy[np.isinf(energy)]=0
        return energy

    def get_energy_best_profile(self,potential="s3dc_dd"):
        normal=True
        energy=np.zeros(len(self._dna))
        scr=self.get_energy_score_by_position(0)
        binding_length=len(scr.get_binding_site())
        for i in range(len(self._dna)-binding_length+1):
            scr=self.get_energy_score_by_position(i)
            if scr is None: continue
            score_per_nucleotide=scr.get_score_per_nucleotide(normal,potential)
            if len(score_per_nucleotide.keys())<=0:continue
            start=i+min([int(nuc) for nuc in sorted(score_per_nucleotide.keys())])-1
            end  =i+max([int(nuc) for nuc in sorted(score_per_nucleotide.keys())])
            if end>len(self._dna): end=len(self._dna)
            try:
              energy[start:end] = max(scr.get_score(normal,potential),energy[start:end].max())
            except Exception as e:
              print("\t\t-- Failed to get best score: %s in %d of [%d - %d]: %s "%(potential,i,start,end,e))
              print("\t\t-- Failed to get best score: NUCLEOTIDES %s MIN %d"%(str(sorted(score_per_nucleotide.keys())),min([int(nuc) for nuc in sorted(score_per_nucleotide.keys())])))
              continue
        return energy


    def get_energy_best_weighted_profile(self,threshold,potential="s3dc_dd"):
        if self.get_fimo_binding_by_threshold(threshold) is None: return None
        binding_score=self.get_fimo_binding_by_threshold(threshold)
        normal=True
        energy=np.zeros(len(self._dna))
        scr=self.get_energy_score_by_position(0)
        binding_length=len(scr.get_binding_site())
        for i in range(len(self._dna)-binding_length+1):
            scr=self.get_energy_score_by_position(i)
            if scr is None: continue
            score_per_nucleotide=scr.get_score_per_nucleotide(normal,potential)
            if len(score_per_nucleotide.keys())<=0:continue
            start=i+min([int(nuc) for nuc in sorted(score_per_nucleotide.keys())])-1
            end  =i+max([int(nuc) for nuc in sorted(score_per_nucleotide.keys())])
            if end>len(self._dna): end=len(self._dna)
            try:
              energy[start:end] = max(scr.get_score(normal,potential),energy[start:end].max())* binding_score[start:end].prod()
            except Exception as e:
              print("\t\t-- Failed to get best weighted score: %s in %d of [%d - %d]: %s "%(potential,i,start,end,e))
              print("\t\t-- Failed to get best weighted score: NUCLEOTIDES %s MIN %d"%(str(sorted(score_per_nucleotide.keys())),min([int(nuc) for nuc in sorted(score_per_nucleotide.keys())])))
              continue
        return energy


    def get_energy_per_nucleotide_profile(self,normal=False,potential="s3dc_dd"):
        energy=np.zeros(len(self._dna))
        count  = np.zeros(len(self._dna))
        scr=self.get_energy_score_by_position(0)
        binding_length=len(scr.get_binding_site())
        for i in range(len(self._dna)-binding_length+1):
            scr=self.get_energy_score_by_position(i)
            if scr is None: continue
            score_per_nucleotide=scr.get_score_per_nucleotide(normal,potential)
            if len(score_per_nucleotide.keys())<=0:continue
            for nuc in sorted(score_per_nucleotide.keys()):
                position = int(nuc) - 1
                if i+position <len(self._dna):
                   energy[i+position] += score_per_nucleotide.get(nuc)
                   count[i+position]  += 1.0 
        energy = energy/count
        energy[np.isnan(energy)]=0
        energy[np.isinf(energy)]=0
        return energy


    def get_energy_per_nucleotide_weighted_profile(self,threshold,normal=False,potential="s3dc_dd"):
        if self.get_fimo_binding_by_threshold(threshold) is None: return None
        binding_score=self.get_fimo_binding_by_threshold(threshold)
        energy=np.zeros(len(self._dna))
        count  = np.zeros(len(self._dna))
        scr=self.get_energy_score_by_position(0)
        binding_length=len(scr.get_binding_site())
        for i in range(len(self._dna)-binding_length+1):
            scr=self.get_energy_score_by_position(i)
            if scr is None: continue
            score_per_nucleotide=scr.get_score_per_nucleotide(normal,potential)
            if len(score_per_nucleotide.keys())<=0:continue
            for nuc in sorted(score_per_nucleotide.keys()):
                position = int(nuc) - 1
                if i+position <len(self._dna):
                   energy[i+position] += score_per_nucleotide.get(nuc)* binding_score[i+position].prod()
                   count[i+position]  += 1.0 * binding_score[i+position].prod()
        energy = energy/count
        energy[np.isnan(energy)]=0
        energy[np.isinf(energy)]=0
        return energy

    def get_profile_score(self,threshold=None,score_type="energy",normal=False,potential="s3dc_dd"):
        """
          score_type = 
                          energy, 
                          energy_best, 
                          energy_per_nucleotide,
                          fimo_binding
                          fimo_score
                          fimo_log_score
          if threshold is not None energies are weighted
        """
        weighted=False
        if threshold is not None: weighted=True

        if weighted:
            if    score_type=="energy":
                  return self.get_energy_weighted_profile(threshold,normal,potential)
            elif  score_type=="energy_best":
                  return self.get_energy_best_weighted_profile(threshold,potential)
            elif  score_type=="energy_per_nucleotide":
                  return self.get_energy_per_nucleotide_weighted_profile(threshold,normal,potential)
            elif  score_type=="fimo_binding":
                  return self.get_fimo_binding_by_threshold(threshold)
            elif  score_type=="fimo_score":
                  return self.get_fimo_score_by_threshold(threshold)
            elif  score_type=="fimo_log_score":
                  return self.get_fimo_log_score_by_threshold(threshold)
            else:
                  return self.get_energy_weighted_profile(threshold,normal,potential)
        else:
            if    score_type=="energy":
                  return self.get_energy_profile(normal,potential)
            elif  score_type=="energy_best":
                  return self.get_energy_best_profile(potential)
            elif  score_type=="energy_per_nucleotide":
                  return self.get_energy_per_nucleotide_profile(normal,potential)
            else:
                  return self.get_energy_profile(normal,potential)

    def to_table(self):
        split_potentials=self._split_potentials
        thresholds=self.get_thresholds()
        tbl_dict={}
        tbl_dict.setdefault("Position",[int(x)+1 for x in range(len(self._dna))])
        tbl_dict.setdefault("Nucleotide",list(self._dna))
        for thr in thresholds:
            tbl_dict.setdefault("binding_"+thr,self.get_fimo_binding_by_threshold(thr))
            tbl_dict.setdefault("fimo_score_"+thr,self.get_fimo_score_by_threshold(thr))
            tbl_dict.setdefault("fimo_log_score_"+thr,self.get_fimo_log_score_by_threshold(thr))
        for potential in split_potentials:
            y=self.get_energy_profile(normal=False,potential=potential)
            if len(y[abs(y)>1.0e-10])<1: continue
            tbl_dict.setdefault(potential,self.get_energy_profile(normal=False,potential=potential))
            tbl_dict.setdefault("normal_"+potential,self.get_energy_profile(normal=True,potential=potential))
            tbl_dict.setdefault(potential+"_x_nucleotide",self.get_energy_per_nucleotide_profile(normal=False,potential=potential))
            tbl_dict.setdefault("normal_"+potential+"_x_nucleotide",self.get_energy_per_nucleotide_profile(normal=True,potential=potential))
            tbl_dict.setdefault("best_"+potential,self.get_energy_best_profile(potential=potential))
        for thr in thresholds:
          for potential in split_potentials:
            y=self.get_energy_weighted_profile(thr,normal=False,potential=potential)
            if len(y[abs(y)>1.0e-10])<1: continue
            tbl_dict.setdefault(potential+"_weighted_"+thr,self.get_energy_weighted_profile(thr,normal=False,potential=potential))
            tbl_dict.setdefault("normal_"+potential+"_weighted_"+thr,self.get_energy_weighted_profile(thr,normal=True,potential=potential))
            tbl_dict.setdefault(potential+"_x_nucleotide"+"_weighted_"+thr,self.get_energy_per_nucleotide_weighted_profile(thr,normal=False,potential=potential))
            tbl_dict.setdefault("normal_"+potential+"_x_nucleotide"+"_weighted_"+thr,self.get_energy_per_nucleotide_weighted_profile(thr,normal=True,potential=potential))
            tbl_dict.setdefault("best_"+potential+"_weighted_"+thr,self.get_energy_best_weighted_profile(thr,potential=potential))
        table=pd.DataFrame(tbl_dict)
        return table

    def write_table(self,output):
        table=self.to_table()
        table.to_csv(output)

    def plot(self,output):
        table=self.to_table()
        columns=table.columns.values.tolist()
        #x = table["Nucleotide"].values
        x = table["Position"].values
        for column in columns:
          if column == "Nucleotide": continue
          if column == "Position": continue
          out   = output+"_"+column+".png"
          title = "Graph of profile %s"%column
          y     = table[column].values
          if len(y[abs(y)>1.0e-10])<1: continue
          plt.plot(x,y,'-', ms=5, lw=2, alpha=0.7, mfc='cyan')
          plt.xlabel('nucleotide')
          plt.ylabel('score profile')
          plt.savefig(out,format="png")
          plt.close()

    def plot_html(self,score_types,normal,energies,output,header_title=False):
        x=np.array([int(i)+1 for i in range(len(self._dna))])
        name_profile = os.path.basename(output)
        html_file    = output+".html"
        thresholds=[]
        thresholds.append(None)
        thresholds.extend(self.get_thresholds())
        '''
        COLOR LIST CSS
                aliceblue, antiquewhite, aqua, aquamarine, azure,
                beige, bisque, black, blanchedalmond, blue,
                blueviolet, brown, burlywood, cadetblue,
                chartreuse, chocolate, coral, cornflowerblue,
                cornsilk, crimson, cyan, darkblue, darkcyan,
                darkgoldenrod, darkgray, darkgrey, darkgreen,
                darkkhaki, darkmagenta, darkolivegreen, darkorange,
                darkorchid, darkred, darksalmon, darkseagreen,
                darkslateblue, darkslategray, darkslategrey,
                darkturquoise, darkviolet, deeppink, deepskyblue,
                dimgray, dimgrey, dodgerblue, firebrick,
                floralwhite, forestgreen, fuchsia, gainsboro,
                ghostwhite, gold, goldenrod, gray, grey, green,
                greenyellow, honeydew, hotpink, indianred, indigo,
                ivory, khaki, lavender, lavenderblush, lawngreen,
                lemonchiffon, lightblue, lightcoral, lightcyan,
                lightgoldenrodyellow, lightgray, lightgrey,
                lightgreen, lightpink, lightsalmon, lightseagreen,
                lightskyblue, lightslategray, lightslategrey,
                lightsteelblue, lightyellow, lime, limegreen,
                linen, magenta, maroon, mediumaquamarine,
                mediumblue, mediumorchid, mediumpurple,
                mediumseagreen, mediumslateblue, mediumspringgreen,
                mediumturquoise, mediumvioletred, midnightblue,
                mintcream, mistyrose, moccasin, navajowhite, navy,
                oldlace, olive, olivedrab, orange, orangered,
                orchid, palegoldenrod, palegreen, paleturquoise,
                palevioletred, papayawhip, peachpuff, peru, pink,
                plum, powderblue, purple, red, rosybrown,
                royalblue, saddlebrown, salmon, sandybrown,
                seagreen, seashell, sienna, silver, skyblue,
                slateblue, slategray, slategrey, snow, springgreen,
                steelblue, tan, teal, thistle, tomato, turquoise,
                violet, wheat, white, whitesmoke, yellow,
                yellowgreen
        '''
        red_color_long   = ["darkred", "firebrick","tomato","red", "orangered", "coral", "orange", "goldenrod","gold","siena","yellow"]
        blue_color_long  = ["navy","darkblue","blue","midnightblue","dodgerblue","steelblue","deepskyblue","skyblue","lightblue","cyan","lightcyan"]
        green_color_long = ["darkgreen","forestgreen","green","lawngreen","springgreen","lightgreen","limegreen","turquoise","aquamarine"]
        red_color   = ["darkred", "firebrick", "red","orangered", "orange", "goldenrod","gold"]
        blue_color  = ["navy","blue","steelblue","deepskyblue","lightblue","cyan","lightcyan"]
        green_color = ["darkgreen","forestgreen","green","springgreen","lightgreen","limegreen","turquoise","aquamarine"]
        selection=[]
        select_main=set()
        trace_list=[]
        for sc in  score_types:
          name_sc=sc
          for ii in range(len(thresholds)):
             thr = thresholds[ii]
             if thr is None: name_thr=name_sc+"_1.0"
             else:           name_thr=name_sc+"_"+thr
             dna=self.get_dna()
             text_dna = np.array(["".join(dna[k:min(k+6,len(dna))])+"..." for k in range(len(dna))])
             if "fimo" in sc:
               if thr is None: continue
               select_main.add(name_sc)
               name = name_thr
               if self.get_profile_score(thr,sc) is None: continue
               y = self.get_profile_score(thr,sc).round(2)
               #upper and lower bound are not necessary, but kept for symmetry of data with e=0
               if "binding" in sc: e=np.zeros(len(y))
               if "binding" in sc and thr=="0.01":
                  #trace = go.Scatter(x = x, y = y, text= text_dna, mode = 'lines', name = thr, legendgroup = thr, line = dict(color=red_color[ii],width = 3),   hoverinfo="x+y+text", visible = True)
                  trace = dict(type="scatter", x = x, y = y, text= text_dna, mode = 'lines', name = thr, legendgroup = thr, line = dict(color=red_color[ii],width = 3),   hoverinfo="x+y+text", visible = True)
               else:
                  visible=False
                  if "binding" in sc: visible="legendonly"
                  #trace = go.Scatter(x = x, y = y, text= text_dna, mode = 'lines', name = thr, legendgroup = thr, line = dict(color=red_color[ii],width = 3),   hoverinfo="x+y+text", visible = visible)
                  trace = dict(type="scatter", x = x, y = y, text= text_dna, mode = 'lines', name = thr, legendgroup = thr, line = dict(color=red_color[ii],width = 3),   hoverinfo="x+y+text", visible = visible)
               trace_list.append(trace)
               selection.append(name)
             if "energy"  in sc :
               for energy in energies:
                   if normal: name_ene = "normal_"+energy
                   else:      name_ene = energy
                   name = name_ene + "_" + name_thr
                   select_main.add(name_ene+"_"+name_sc)
                   if self.get_profile_score(thr,sc,normal,energy) is None : continue
                   y = self.get_profile_score(thr,sc,normal,energy).round(2)
                   visible=False
                   if thr is None: thre="1.0"
                   else: thre=thr
                   #trace = go.Scatter(x = x, y = y, text= text_dna, mode = 'lines', name = thre, legendgroup = thre, line = dict(color=blue_color[ii],width = 3),   hoverinfo="x+y+text", visible = visible)
                   trace = dict(type="scatter", x = x, y = y, text= text_dna, mode = 'lines', name = thre, legendgroup = thre, line = dict(color=blue_color[ii],width = 3),   hoverinfo="x+y+text", visible = visible)
                   trace_list.append(trace)
                   selection.append(name)
        # Define the layout and axis names #
        config = {'linkText': "", 'scrollZoom': True}
        # Determine a dropdown with the variables that will be shown in this plot #
        fimo_selection=[s for s in select_main if "fimo" in s]
        ener_selection=[s for s in select_main if "energy" in s]
        order_selection=sorted(fimo_selection)[:]
        order_selection.extend(sorted(ener_selection))
        buttons=[]
        for plt_name in order_selection:
            features=plt_name.split("_")
            score=features[:]
            if "fimo" in plt_name: score=features[1:]
            arg_visible=[]
            if "binding" in plt_name:
              for sel in selection:
                  features_selection=sel.split("_")
                  thr_sel=features_selection[-1]
                  test_sel="_".join(features_selection[0:-1])
                  if plt_name==test_sel:
                     if thr_sel=="0.01":arg_visible.append(True)
                     else:arg_visible.append("legendonly")
                  else:
                     arg_visible.append(False)
            else:
              for sel in selection:
                  features_selection=sel.split("_")
                  thr_sel=features_selection[-1]
                  test_sel="_".join(features_selection[0:-1])
                  if plt_name==test_sel:
                     if  thr_sel=="0.5":arg_visible.append(True)
                     else:arg_visible.append("legendonly")
                  else:
                     arg_visible.append(False)
            arg_yaxis={}
            arg_yaxis.setdefault('titlefont',dict(size=12))
            if "fimo" in plt_name:
               if "log" in plt_name: score[0]="-log(P-value)"
               if "binding" in plt_name: score.append("ratio")
               label="FIMO: "+" ".join(score)
            if "energy" in plt_name:
               label="POTENTIAL: "+ " ".join(score)
            arg_yaxis.setdefault('title',label)
            arguments=[]
            arguments.append(dict(visible=arg_visible))
            arguments.append(dict(yaxis=arg_yaxis))
            buttons.append(dict(label=label,method='update',args=arguments))

        #updatemenus = [dict(buttons=buttons, direction='down', showactive=True, pad= {'r': 0, 't':  15}, x = 0.25 , y = 1.85, font = {'size': 12}, bgcolor = "white")]
        updatemenus = [dict(buttons=buttons, direction='down', showactive=True, pad= {'r': 0, 't':  15}, x = 0.25 , y = 1.95, font = {'size': 12}, bgcolor = "white")]
        annotations = [dict(text='P-value\nthresholds:', x=1.08, y=1.15, xref='paper', yref='paper', showarrow=False, font=dict(size=12), align="center")]
        xaxis_dict=dict(
                        title         ='Nucleotide positions',
                        rangeslider   =dict(visible=True, bordercolor="grey", borderwidth=3),
                        titlefont     =dict(size=18), 
                        tickformat    =',d'
                        )
        if header_title:
           layout=dict(title="PROFILE: "+name_profile,updatemenus=updatemenus,annotations=annotations,plot_bgcolor='rgba(0,0,0,0)',paper_bgcolor='rgba(0,0,0,0)',yaxis= dict(title = "Ratio", titlefont=dict(size=16)),xaxis=xaxis_dict)
        else:
           layout=dict(title="",updatemenus=updatemenus,annotations=annotations,plot_bgcolor='rgba(0,0,0,0)',paper_bgcolor='rgba(0,0,0,0)',yaxis= dict(title = "Ratio", titlefont=dict(size=16)),xaxis=xaxis_dict)

        # Create the plot #
        fig = dict(data=trace_list, layout=layout)
        pltly(fig, filename=html_file, show_link=False, auto_open=False)

        

class twoprofile(object):

    def __init__(self,a,b,thresholds=["0.05", "0.001"], potential=None):
          self._profile_one     = a
          self._profile_two     = b
          self._thresholds      = thresholds[:]
          self._potential       = "all"
          self._split_potentials=["3d","3dc","local","pair","s3dc","s3dc_di","s3dc_dd"]
          self._score_types     =["fimo_binding","fimo_score","fimo_log_score","energy","energy_best","energy_per_nucleotide"]
          if potential is not None: self._potential=potential

    def get_profile_one(self):
          return self._profile_one

    def set_profile_one(self,x):
          self._profile_one = x

    def get_profile_two(self):
          return self._profile_two

    def set_profile_two(self,x):
          self._profile_two = x

    def set_thresholds(self,thresholds):
          self._thresholds=thresholds

    def get_thresholds(self):
          return self._thresholds

    def set_potential(self,potential):
          self._potential = potential

    def get_potential(self):
          return self._potential

    def set_split_potentials(self,split_potentials):
        self._split_potentials=split_potentials

    def get_split_potentials(self):
        return self._split_potentials

    def set_score_types(self,score_types):
        self._score_types=score_types

    def get_score_types(self):
        return self._score_types

    def get_paired_difference(self,threshold=None,score_type="energy",normal=False,potential="s3dc_dd"):
          profile_one=self.get_profile_one().get_profile_score(threshold,score_type,normal,potential)
          profile_two=self.get_profile_two().get_profile_score(threshold,score_type,normal,potential)
          dna_length_one = len(self.get_profile_one().get_dna())
          dna_length_two = len(self.get_profile_two().get_dna())
          if dna_length_one != dna_length_two: return None
          if len(profile_two) == len(profile_one): return (profile_two - profile_one)
          #else: return (profile_two.mean(0) - profile_one.mean(0)
          else: return (profile_two[:min(len(profile_two),len(profile_one)),:] - profile_one[:min(len(profile_two),len(profile_one)),:])

    def get_standard_difference(self,threshold=None,score_type="energy",normal=False,potential="s3dc_dd"):
          profile_one=self.get_profile_one().get_profile_score(threshold,score_type,normal,potential)
          profile_two=self.get_profile_two().get_profile_score(threshold,score_type,normal,potential)
          dna_length_one = len(self.get_profile_one().get_dna())
          dna_length_two = len(self.get_profile_two().get_dna())
          if dna_length_one != dna_length_two: return None
          difference_mean = profile_two.mean(0) - profile_one.mean(0)
          difference_rmsd = (profile_two.std(0) + profile_one.std(0))/2
          difference_lower= difference_mean - difference_rmsd
          difference_upper= difference_mean + difference_rmsd
          x = difference_lower
          x = np.vstack((x,difference_upper))
          return x
           
    def get_wilcoxon_score(self,threshold=None,score_type="energy",normal=False,potential="s3dc_dd"):
          profile_one=self.get_profile_one().get_profile_score(threshold,score_type,normal,potential)
          profile_two=self.get_profile_two().get_profile_score(threshold,score_type,normal,potential)
          dna_length =min(len(self.get_profile_one().get_dna()),len(self.get_profile_two().get_dna()))
          wscore=[]
          for i in range(dna_length):
              x   = profile_one[:,i]
              y   = profile_two[:,i]
              if len(x) == len(y):
                try:
                   if abs(x-y).sum()> 0.001:
                     w,p = stats.wilcoxon(x,y)
                     wscore.append(-np.log(p))
                   else:
                     wscore.append(0.0)
                except Exception as e:
                   wscore.append(0.0)
              else:
                 wscore.append(0.0)
          wscore=np.array(wscore)
          return wscore

    def get_mannwhitney_score(self,threshold=None,score_type="energy",normal=False,potential="s3dc_dd"):
          profile_one=self.get_profile_one().get_profile_score(threshold,score_type,normal,potential)
          profile_two=self.get_profile_two().get_profile_score(threshold,score_type,normal,potential)
          dna_length =min(len(self.get_profile_one().get_dna()),len(self.get_profile_two().get_dna()))
          wscore=[]
          for i in range(dna_length):
              x   = profile_one[:,i]
              y   = profile_two[:,i]
              if len(x) == len(y):
                 try:
                   if abs(x-y).sum()> 0.001:
                      mw,p=stats.mannwhitneyu(x,y)
                      wscore.append(-np.log(p))
                   else:
                      wscore.append(0.0)
                 except Exception as e:
                   wscore.append(0.0)
              else:
                try:
                  mw,p=stats.mannwhitneyu(x,y)
                  wscore.append(-np.log(p))
                except Exception as e:
                   wscore.append(0.0)
          wscore=np.array(wscore)
          return wscore

    def to_stats_table(self):
        split_potentials=self._split_potentials
        thresholds=self.get_thresholds()[:]
        thresholds.append(None)
        score_types=self.get_score_types()[:]
        normalization=[True,False]
        tbl_dict={}
        dna_length_one = len(self.get_profile_one().get_dna())
        dna_length_two = len(self.get_profile_two().get_dna())
        if dna_length_one != dna_length_two: return None
        dna_length = dna_length_one
        tbl_dict.setdefault("Position",[int(x)+1 for x in range(dna_length)])
        for sc in  score_types:
          if "fimo" in sc: 
              name_sc=sc
              for thr in thresholds:
                  if thr is None: continue
                  name=name_sc+"_"+thr
                  tbl_dict.setdefault(name+"_wilcoxon_score",self.get_wilcoxon_score(thr,sc))
                  tbl_dict.setdefault(name+"_mannwhitney_score",self.get_mannwhitney_score(thr,sc))
          if "energy"  in sc :
              for potential in split_potentials:
                for thr in thresholds:
                  if thr is None: namep = potential + "_" + sc
                  else:           namep = potential + "_" + sc +"_weighted_"+thr
                  for normal in normalization:
                     if "best" in sc and normal:continue
                     if normal: name = "normal_"+namep 
                     tbl_dict.setdefault(name+"_wilcoxon_score",self.get_wilcoxon_score(thr,sc,normal,potential))
                     tbl_dict.setdefault(name+"_mannwhitney_score",self.get_wilcoxon_score(thr,sc,normal,potential))
        table=pd.DataFrame(tbl_dict)
        return table

    def to_mean_table(self):
        split_potentials=self._split_potentials
        thresholds=self.get_thresholds()[:]
        thresholds.append(None)
        score_types=self.get_score_types()[:]
        normalization=[True,False]
        tbl_dict={}
        dna_length_one = len(self.get_profile_one().get_dna())
        dna_length_two = len(self.get_profile_two().get_dna())
        if dna_length_one != dna_length_two: return None
        dna_length = dna_length_one
        tbl_dict.setdefault("Position",[int(x)+1 for x in range(dna_length)])
        for sc in  score_types:
          if "fimo" in sc: 
              name_sc=sc
              for thr in thresholds:
                  if thr is None: continue
                  name=name_sc+"_"+thr
                  tbl_dict.setdefault(name+"_paired",self.get_paired_difference(thr,sc).mean(0))
                  tbl_dict.setdefault(name+"_standard",self.get_standard_difference(thr,sc).mean(0))
          if "energy"  in sc :
              for potential in split_potentials:
                for thr in thresholds:
                  if thr is None: namep = potential + "_" + sc
                  else:           namep = potential + "_" + sc +"_weighted_"+thr
                  for normal in normalization:
                     if "best" in sc and normal:continue
                     if normal: name = "normal_"+namep 
                     tbl_dict.setdefault(name+"_paired",self.get_paired_difference(thr,sc,normal,potential).mean(0))
                     tbl_dict.setdefault(name+"_standard",self.get_standard_difference(thr,sc,normal,potential).mean(0))
        table=pd.DataFrame(tbl_dict)
        return table

    def to_rmsd_table(self):
        split_potentials=self._split_potentials
        thresholds=self.get_thresholds()[:]
        thresholds.append(None)
        score_types=self.get_score_types()[:]
        normalization=[True,False]
        tbl_dict={}
        dna_length_one = len(self.get_profile_one().get_dna())
        dna_length_two = len(self.get_profile_two().get_dna())
        if dna_length_one != dna_length_two: return None
        dna_length = dna_length_one
        tbl_dict.setdefault("Position",[int(x)+1 for x in range(dna_length)])
        for sc in  score_types:
          if "fimo" in sc: 
              name_sc=sc
              for thr in thresholds:
                  if thr is None: continue
                  name=name_sc+"_"+thr
                  tbl_dict.setdefault(name+"_paired",self.get_paired_difference(thr,sc).std(0))
                  tbl_dict.setdefault(name+"_standard",self.get_standard_difference(thr,sc).std(0))
          if "energy"  in sc :
              for potential in split_potentials:
                for thr in thresholds:
                  if thr is None: namep = potential + "_" + sc
                  else:           namep = potential + "_" + sc +"_weighted_"+thr
                  for normal in normalization:
                     if "best" in sc and normal:continue
                     if normal: name = "normal_"+namep 
                     tbl_dict.setdefault(name+"_paired",self.get_paired_difference(thr,sc,normal,potential).std(0))
                     tbl_dict.setdefault(name+"_standard",self.get_standard_difference(thr,sc,normal,potential).std(0))
        table=pd.DataFrame(tbl_dict)
        return table
 

    def write_stats_table(self,output):
        table=self.to_stats_table()
        if table is not None: table.to_csv(output)
 
    def write_mean_table(self,output):
        table=self.to_mean_table()
        if table is not None: table.to_csv(output)

    def write_rmsd_table(self,output):
        table=self.to_rmsd_table()
        if table is not None: table.to_csv(output)


    def stplot(self,output):
        table=self.to_stats_table()
        if table is not None:
         columns=table.columns.values.tolist()
         x = table["Position"].values
         for column in columns:
           if column == "Position": continue
           if column == "Nucleotide": continue
           out   = output+"_"+column+".png"
           title = "Graph of profile %s"%column
           y     = table[column].values
           if len(y[abs(y)>1.0e-10])<1: continue
           plt.plot(x,y,'-', ms=5, lw=2, alpha=0.7, mfc='cyan')
           plt.xlabel('DNA position')
           plt.ylabel('statistical score')
           plt.savefig(out,format="png")
           plt.close()
   

    def plot(self,output):
        table=self.to_mean_table()
        error=self.to_rmsd_table()
        if table is not None and error is not None:
         columns=table.columns.values.tolist()
         x = table["Position"].values
         for column in columns:
           if column == "Position": continue
           if column == "Nucleotide": continue
           out   = output+"_"+column+".png"
           title = "Graph of profile %s"%column
           y     = table[column].values
           e     = error[column].values
           if len(y[abs(y)>1.0e-10])<1: continue
           plt.plot(x,y,'-', ms=5, lw=2, alpha=0.7, mfc='cyan')
           plt.xlabel('DNA position ')
           plt.ylabel('score differencial profile')
           plt.fill_between(x,y-e,y+e,alpha=0.2)
           plt.savefig(out,format="png")
           plt.close()

 
    def twoplot(self,output):
        profile_one=self.get_profile_one()
        profile_two=self.get_profile_two()
        table_one=profile_one.to_mean_table()
        error_one=profile_one.to_rmsd_table()
        table_two=profile_two.to_mean_table()
        error_two=profile_two.to_rmsd_table()
        columns_one=table_one.columns.values.tolist()
        columns_two=table_two.columns.values.tolist()
        x = table_one["Position"].values
        for column in columns_one:
          if column == "Position": continue
          if column == "Nucleotide": continue
          if column not in columns_two: continue
          out   = output+"_"+column+".png"
          title = "Graph of profile %s"%column
          y     = table_one[column].values
          ey    = error_one[column].values
          z     = table_two[column].values
          ez    = error_two[column].values
          if len(y[abs(y)>1.0e-10])<1: continue
          if len(z[abs(z)>1.0e-10])<1: continue
          plt.xlabel('DNA position ')
          plt.ylabel('score profile')
          plt.plot(x,y,'-', ms=5, lw=2, alpha=0.7, mfc='cyan', color='blue')
          plt.fill_between(x,y-ey,y+ey,alpha=0.2, color='blue')
          plt.plot(x,z,'-', ms=5, lw=2, alpha=0.7, mfc='red', color='firebrick')
          plt.fill_between(x,z-ez,z+ez,alpha=0.2,color='red')
          plt.savefig(out,format="png")
          plt.close()

    def plot_html(self,score_types,normal,energies,output,header_title=False):
        name_profile = os.path.basename(output)
        html_file    = output+".html"
        thresholds=[]
        thresholds.append(None)
        thresholds.extend(self.get_thresholds())
        red_color_long   = ["darkred", "firebrick","tomato","red", "orangered", "coral", "orange", "goldenrod","gold","siena","yellow"]
        blue_color_long  = ["navy","darkblue","blue","midnightblue","dodgerblue","steelblue","deepskyblue","skyblue","lightblue","cyan","lightcyan"]
        green_color_long = ["darkgreen","forestgreen","green","lawngreen","springgreen","lightgreen","limegreen","turquoise","aquamarine"]
        red_color   = ["darkred", "firebrick", "red","orangered", "orange", "goldenrod","gold"]
        blue_color  = ["navy","blue","steelblue","deepskyblue","lightblue","cyan","lightcyan"]
        green_color = ["darkgreen","forestgreen","green","springgreen","lightgreen","limegreen","turquoise","aquamarine"]
        profile_one=self.get_profile_one()
        profile_two=self.get_profile_two()
        dna_length =min(len(self.get_profile_one().get_dna()),len(self.get_profile_two().get_dna()))
        x=np.array([int(i)+1 for i in range(dna_length)])
        selection=[]
        select_main=set()
        trace_list=[]
        for sc in  score_types:
          name_sc=sc
          for ii in range(len(thresholds)):
             thr = thresholds[ii]
             if thr is None: name_thr=name_sc+"_1.0"
             else:           name_thr=name_sc+"_"+thr
             dna_one=profile_one.get_dna()
             dna_two=profile_two.get_dna()
             text_one = np.array(["".join(dna_one[k:min(k+6,len(dna_one))])+"..." for k in range(len(dna_one))])
             text_two = np.array(["".join(dna_two[k:min(k+6,len(dna_two))])+"..." for k in range(len(dna_two))])
             if "fimo" in sc:
               if thr is None: continue
               select_main.add(name_sc)
               name = name_thr
               if profile_one.get_profile_score(thr,sc) is None or profile_two.get_profile_score(thr,sc) is None: continue
               y = profile_one.get_profile_score(thr,sc).mean(0).round(2)
               e = profile_one.get_profile_score(thr,sc).std(0).round(2)
               z = profile_two.get_profile_score(thr,sc).mean(0).round(2)
               d = profile_two.get_profile_score(thr,sc).std(0).round(2)
               #upper and lower bound are not necessary, but kept for symmetry of data with e=0
               if "binding" in sc: 
                   e=np.zeros(len(y))
                   d=np.zeros(len(z))
               legend_one="P1 P-val<"+thr
               legend_two="P2 P-val<"+thr
               if "binding" in sc and thr=="0.01":
                  #trace_one = go.Scatter(x = x, y = y, text=text_one, mode = 'lines', name = legend_one, legendgroup = thr, line = dict(color=red_color[ii],width = 3), fill='tonexty',  hoverinfo="x+y+text", visible = True)
                  #upper_bound_one = go.Scatter(x = x, y = y+e, mode = 'lines', name = legend_one, legendgroup = thr, line = dict(color=red_color[ii],width = 0.8), fill='tonexty',  showlegend = False, hoverinfo='none', visible = True)
                  #lower_bound_one = go.Scatter(x = x, y = y-e, mode = 'lines', name = legend_one, legendgroup = thr, line = dict(color=red_color[ii],width = 0.8), showlegend = False,  hoverinfo='none',visible = True)
                  trace_one =       dict(type="scatter", x = x, y = y, text=text_one, mode = 'lines', name = legend_one, legendgroup = thr, line = dict(color=red_color[ii],width = 3), fill='tonexty',  hoverinfo="x+y+text", visible = True)
                  upper_bound_one = dict(type="scatter", x = x, y = y+e, mode = 'lines', name = legend_one, legendgroup = thr, line = dict(color=red_color[ii],width = 0.8), fill='tonexty',  showlegend = False, hoverinfo='none', visible = True)
                  lower_bound_one = dict(type="scatter", x = x, y = y-e, mode = 'lines', name = legend_one, legendgroup = thr, line = dict(color=red_color[ii],width = 0.8), showlegend = False,  hoverinfo='none',visible = True)
                  #trace_two = go.Scatter(x = x, y = z, text=text_two, mode = 'lines', name = legend_two, legendgroup = thr, line = dict(color=blue_color[ii],width = 3), fill='tonexty', hoverinfo="x+y+text", visible = True)
                  #upper_bound_two = go.Scatter(x = x, y = z+d, mode = 'lines', name = legend_two, legendgroup = thr, line = dict(color=blue_color[ii],width = 0.8), fill='tonexty',  showlegend = False, hoverinfo='none',visible = True)
                  #lower_bound_two = go.Scatter(x = x, y = z-d, mode = 'lines', name = legend_two, legendgroup = thr, line = dict(color=blue_color[ii],width = 0.8), showlegend = False, hoverinfo='none',visible = True)
                  trace_two =       dict(type="scatter", x = x, y = z, text=text_two, mode = 'lines', name = legend_two, legendgroup = thr, line = dict(color=blue_color[ii],width = 3), fill='tonexty', hoverinfo="x+y+text", visible = True)
                  upper_bound_two = dict(type="scatter", x = x, y = z+d, mode = 'lines', name = legend_two, legendgroup = thr, line = dict(color=blue_color[ii],width = 0.8), fill='tonexty',  showlegend = False, hoverinfo='none',visible = True)
                  lower_bound_two = dict(type="scatter", x = x, y = z-d, mode = 'lines', name = legend_two, legendgroup = thr, line = dict(color=blue_color[ii],width = 0.8), showlegend = False, hoverinfo='none',visible = True)
               else:
                  visible=False
                  if "binding" in sc: visible="legendonly"
                  #trace_one = go.Scatter(x = x, y = y, text=text_one,mode = 'lines', name = legend_one, legendgroup = thr, line = dict(color=red_color[ii],width = 3),  fill='tonexty',  hoverinfo="x+y+text", visible = visible)
                  #upper_bound_one = go.Scatter(x = x, y = y+e, mode = 'lines', name = legend_one, legendgroup = thr, line = dict(color=red_color[ii],width = 0.8), fill='tonexty', showlegend = False, hoverinfo='none', visible = visible)
                  #lower_bound_one = go.Scatter(x = x, y = y-e, mode = 'lines', name = legend_one, legendgroup = thr, line = dict(color=red_color[ii],width = 0.8),  showlegend = False, hoverinfo='none', visible = visible)
                  #trace_two = go.Scatter(x = x, y = z, text=text_two, mode = 'lines', name = legend_two, legendgroup = thr, line = dict(color=blue_color[ii],width = 3),  fill='tonexty',   hoverinfo="x+y+text", visible = visible)
                  #upper_bound_two = go.Scatter(x = x, y = z+d, mode = 'lines', name = legend_two, legendgroup = thr, line = dict(color=blue_color[ii],width = 0.8), fill='tonexty', showlegend = False, hoverinfo='none', visible = visible)
                  #lower_bound_two = go.Scatter(x = x, y = z-d, mode = 'lines', name = legend_two, legendgroup = thr, line = dict(color=blue_color[ii],width = 0.8),  showlegend = False, hoverinfo='none', visible = visible)
                  trace_one =       dict(type="scatter",x = x, y = y, text=text_one,mode = 'lines', name = legend_one, legendgroup = thr, line = dict(color=red_color[ii],width = 3),  fill='tonexty',  hoverinfo="x+y+text", visible = visible)
                  upper_bound_one = dict(type="scatter",x = x, y = y+e, mode = 'lines', name = legend_one, legendgroup = thr, line = dict(color=red_color[ii],width = 0.8), fill='tonexty', showlegend = False, hoverinfo='none', visible = visible)
                  lower_bound_one = dict(type="scatter",x = x, y = y-e, mode = 'lines', name = legend_one, legendgroup = thr, line = dict(color=red_color[ii],width = 0.8),  showlegend = False, hoverinfo='none', visible = visible)
                  trace_two =       dict(type="scatter",x = x, y = z, text=text_two, mode = 'lines', name = legend_two, legendgroup = thr, line = dict(color=blue_color[ii],width = 3),  fill='tonexty',   hoverinfo="x+y+text", visible = visible)
                  upper_bound_two = dict(type="scatter",x = x, y = z+d, mode = 'lines', name = legend_two, legendgroup = thr, line = dict(color=blue_color[ii],width = 0.8), fill='tonexty', showlegend = False, hoverinfo='none', visible = visible)
                  lower_bound_two = dict(type="scatter",x = x, y = z-d, mode = 'lines', name = legend_two, legendgroup = thr, line = dict(color=blue_color[ii],width = 0.8),  showlegend = False, hoverinfo='none', visible = visible)
               trace_list.append(lower_bound_one)
               trace_list.append(trace_one)
               trace_list.append(upper_bound_one)
               trace_list.append(lower_bound_two)
               trace_list.append(trace_two)
               trace_list.append(upper_bound_two)
               selection.append(name)
               selection.append(name)
               selection.append(name)
               selection.append(name)
               selection.append(name)
               selection.append(name)
             if "energy"  in sc :
               for energy in energies:
                   if normal: name_ene = "normal_"+energy
                   else:      name_ene = energy
                   name = name_ene + "_" + name_thr
                   select_main.add(name_ene+"_"+name_sc)
                   if profile_one.get_profile_score(thr,sc,normal,energy) is None or profile_two.get_profile_score(thr,sc,normal,energy) is None: continue
                   y = profile_one.get_profile_score(thr,sc,normal,energy).mean(0).round(2)
                   e = profile_one.get_profile_score(thr,sc,normal,energy).std(0).round(2)
                   z = profile_two.get_profile_score(thr,sc,normal,energy).mean(0).round(2)
                   d = profile_two.get_profile_score(thr,sc,normal,energy).std(0).round(2)
                   visible=False
                   if thr is None: thre="1.0"
                   else: thre=thr
                   legend_one="P1 P-val<"+thre
                   legend_two="P2 P-val<"+thre
                   #trace_one =       go.Scatter(x = x, y = y, text=text_one, mode = 'lines', name = legend_one, legendgroup = thre, line = dict(color=red_color[ii],width = 3),  fill='tonexty', hoverinfo="x+y+text", visible = visible)
                   #upper_bound_one = go.Scatter(x = x, y = y+e, mode = 'lines', name = legend_one, legendgroup = thre, line = dict(color=red_color[ii],width = 0.8), fill='tonexty', showlegend = False, hoverinfo='none', visible = visible)
                   #lower_bound_one = go.Scatter(x = x, y = y-e, mode = 'lines', name = legend_one, legendgroup = thre, line = dict(color=red_color[ii],width = 0.8),  showlegend = False, hoverinfo='none', visible = visible)
                   #trace_two =       go.Scatter(x = x, y = z, text=text_two, mode = 'lines', name = legend_two, legendgroup = thre, line = dict(color=blue_color[ii],width = 3),  fill='tonexty', hoverinfo="x+y+text", visible = visible)
                   #upper_bound_two = go.Scatter(x = x, y = z+d, mode = 'lines', name = legend_two, legendgroup = thre, line = dict(color=blue_color[ii],width = 0.8), fill='tonexty', showlegend = False, hoverinfo='none', visible = visible)
                   #lower_bound_two = go.Scatter(x = x, y = z-d, mode = 'lines', name = legend_two, legendgroup = thre, line = dict(color=blue_color[ii],width = 0.8),  showlegend = False, hoverinfo='none', visible = visible)
                   trace_one =       dict(type="scatter",x = x, y = y, text=text_one, mode = 'lines', name = legend_one, legendgroup = thre, line = dict(color=red_color[ii],width = 3),  fill='tonexty', hoverinfo="x+y+text", visible = visible)
                   upper_bound_one = dict(type="scatter",x = x, y = y+e, mode = 'lines', name = legend_one, legendgroup = thre, line = dict(color=red_color[ii],width = 0.8), fill='tonexty', showlegend = False, hoverinfo='none', visible = visible)
                   lower_bound_one = dict(type="scatter",x = x, y = y-e, mode = 'lines', name = legend_one, legendgroup = thre, line = dict(color=red_color[ii],width = 0.8),  showlegend = False, hoverinfo='none', visible = visible)
                   trace_two =       dict(type="scatter",x = x, y = z, text=text_two, mode = 'lines', name = legend_two, legendgroup = thre, line = dict(color=blue_color[ii],width = 3),  fill='tonexty', hoverinfo="x+y+text", visible = visible)
                   upper_bound_two = dict(type="scatter",x = x, y = z+d, mode = 'lines', name = legend_two, legendgroup = thre, line = dict(color=blue_color[ii],width = 0.8), fill='tonexty', showlegend = False, hoverinfo='none', visible = visible)
                   lower_bound_two = dict(type="scatter",x = x, y = z-d, mode = 'lines', name = legend_two, legendgroup = thre, line = dict(color=blue_color[ii],width = 0.8),  showlegend = False, hoverinfo='none', visible = visible)
                   trace_list.append(lower_bound_one)
                   trace_list.append(trace_one)
                   trace_list.append(upper_bound_one)
                   trace_list.append(lower_bound_two)
                   trace_list.append(trace_two)
                   trace_list.append(upper_bound_two)
                   selection.append(name)
                   selection.append(name)
                   selection.append(name)
                   selection.append(name)
                   selection.append(name)
                   selection.append(name)
        # Define the layout and axis names #
        config = {'linkText': "", 'scrollZoom': True}
        # Determine a dropdown with the variables that will be shown in this plot #
        fimo_selection=[s for s in select_main if "fimo" in s]
        ener_selection=[s for s in select_main if "energy" in s]
        order_selection=sorted(fimo_selection)[:]
        order_selection.extend(sorted(ener_selection))
        buttons=[]
        for plt_name in order_selection:
            features=plt_name.split("_")
            score=features[:]
            if "fimo" in plt_name: score=features[1:]
            arg_visible=[]
            if "binding" in plt_name:
              for sel in selection:
                  features_selection=sel.split("_")
                  thr_sel=features_selection[-1]
                  test_sel="_".join(features_selection[0:-1])
                  if plt_name==test_sel:
                     if thr_sel=="0.01":arg_visible.append(True)
                     else:arg_visible.append("legendonly")
                  else:
                     arg_visible.append(False)
            else:
              for sel in selection:
                  features_selection=sel.split("_")
                  thr_sel=features_selection[-1]
                  test_sel="_".join(features_selection[0:-1])
                  if plt_name==test_sel:
                     if  thr_sel=="0.5":arg_visible.append(True)
                     else:arg_visible.append("legendonly")
                  else:
                     arg_visible.append(False)
            arg_yaxis={}
            arg_yaxis.setdefault('titlefont',dict(size=12))
            if "fimo" in plt_name:
               if "log" in plt_name: score[0]="-log(P-value)"
               if "binding" in plt_name: score.append("ratio")
               label="FIMO: "+" ".join(score)
            if "energy" in plt_name:
               label="POTENTIAL: "+ " ".join(score)
            arg_yaxis.setdefault('title',label)
            arguments=[]
            arguments.append(dict(visible=arg_visible))
            arguments.append(dict(yaxis=arg_yaxis))
            buttons.append(dict(label=label,method='update',args=arguments))

        #updatemenus = [dict(buttons=buttons, direction='down', showactive=True,  pad= {'r': 0, 't':  15}, x = 0.25 , y = 1.85, font = {'size': 12}, bgcolor = "white")]
        updatemenus = [dict(buttons=buttons, direction='down', showactive=True, pad= {'r': 0, 't':  15}, x = 0.25 , y = 1.95, font = {'size': 12}, bgcolor = "white")]
        annotations = [dict(text='P-value\nthresholds:', x=1.08, y=1.15, xref='paper', yref='paper', showarrow=False, font=dict(size=12), align="center")]
        xaxis_dict=dict(
                        title         ='Nucleotide positions',
                        rangeslider   =dict(visible=True, bordercolor="grey", borderwidth=3),
                        titlefont     =dict(size=18), 
                        tickformat    =',d'
                        )
        if header_title:
            layout=dict(title="PROFILE: "+name_profile,updatemenus=updatemenus,annotations=annotations,plot_bgcolor='rgba(0,0,0,0)',paper_bgcolor='rgba(0,0,0,0)',yaxis= dict(title = "Ratio", titlefont=dict(size=16)),xaxis=xaxis_dict)
        else:
            layout=dict(title="",updatemenus=updatemenus,annotations=annotations,plot_bgcolor='rgba(0,0,0,0)',paper_bgcolor='rgba(0,0,0,0)',yaxis= dict(title = "Ratio", titlefont=dict(size=16)),xaxis=xaxis_dict)
        # Create the plot #
        fig = dict(data=trace_list, layout=layout)
        pltly(fig, filename=html_file, show_link=False, auto_open=False)

            

