import os, sys, re
from collections import Counter
import ConfigParser
import itertools
import numpy as np
import optparse
import subprocess
import copy
import shutil
import hashlib

# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.dirname(__file__))

# Append scripts path to python path #
sys.path.append(scripts_path)

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Import functions #
import functions

# Import SBI modules #
from SBI.structure import PDB
from SBI.structure import Chain
from SBI.structure.chain import ChainOfProtein, ChainOfNucleotide
from SBI.data import aminoacids3to1, aminoacids_polarity_boolean, nitrogenous_bases

# Import isPBM2.0 modules #
import dssp, x3dna, contacts, triads,  tmalign
import pwm_pbm as PWM

##########################
# Aminoacid propensities #
##########################
"""
Mean amino acid propensities for alpha-helix and beta-strand conformations

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3495713/pdf/1472-6807-12-18.pdf
doi: 10.1186/1472-6807-12-18
"""

exposed_helix_propensities_dict = {'V': 0.83, 'I': 0.96, 'L': 1.16, 'M': 1.03, 'P': 0.48, 'A': 1.43, 'C': 0.63, 'F': 0.88, 'Y': 0.91, 'W': 0.87, 'Q': 1.34, 'S': 0.74, 'T': 0.72, 'N':0.74, 'H': 0.90, 'D': 0.91, 'K': 1.25, 'E': 1.51, 'R': 1.31, 'G': 0.28 }
buried_helix_propensities_dict = {'V': 0.89, 'I': 1.01, 'L': 1.27, 'M': 1.29, 'P': 0.41, 'A': 1.37, 'C': 0.85, 'F': 0.99, 'Y': 0.98, 'W': 1.09, 'Q': 1.21, 'S': 0.80, 'T': 0.84, 'N':0.77, 'H': 0.85, 'D': 0.73, 'K': 1.13, 'E': 1.25, 'R': 1.13, 'G': 0.59 }
exposed_b_strand_propensities_dict = {'V': 2.31, 'I': 2.02, 'L': 1.18, 'M': 1.01, 'P': 0.49, 'A': 0.48, 'C': 1.24, 'F': 1.50, 'Y': 1.71, 'W': 1.90, 'Q': 0.96, 'S': 0.86, 'T': 1.58, 'N':0.71, 'H': 1.15, 'D': 0.61, 'K': 1.14, 'E': 0.89, 'R': 1.27, 'G': 0.41 }
buried_b_strand_propensities_dict = {'V': 1.57, 'I': 1.39, 'L': 0.93, 'M': 0.84, 'P': 0.42, 'A': 0.72, 'C': 1.07, 'F': 1.10, 'Y': 1.12, 'W': 0.91, 'Q': 0.82, 'S': 0.85, 'T': 1.08, 'N':0.76, 'H': 0.98, 'D': 0.76, 'K': 0.98, 'E': 0.86, 'R': 0.82, 'G': 0.81 }

###############
#  METHODS  #
###############

def get_patches(pdb_obj, contacts_obj, radius, fragment_restrict):
    
    # Set the maximum distance between the residues that define the patch max-dictance between Aas
    distance_start = float(config.get("Parameters","aa_contact_in_patch"))
    distance=distance_start
    interface_distance   =  float(config.get("Parameters","interface_distance"))
    patch_size           =  float(config.get("Parameters","patch_size"))
    if radius <= 0: radius=interface_distance
    patches_complete={}
    #Get patches per chain
    for pdb_chain_obj in pdb_obj.chains:
       if pdb_chain_obj.chaintype =='P':
          patches = []
          big_patch = True
          pdb_chain_id = pdb_chain_obj.chain
          sequence = pdb_chain_obj.gapped_protein_sequence
          indexes  = pdb_chain_obj.protein_idx.split(";")
          done = set()
          while (big_patch == True):
             for contact1 in contacts_obj.get_contacts():
                 big_patch = False
                 aa1 = contact1._A_residue_obj.identifier
                 chain1 = contact1._A_chain
                 if aa1 in done: continue
                 if contact1.get_contact_distance() > interface_distance: continue
                 if contact1.get_contact_distance() > radius: continue
                 # Check  the protein chain
                 if pdb_chain_obj.chain == chain1:
                     patch = []
                     singleton_patch = []
                     residue_num1= 1 + [pos for pos in range(len(indexes)) if str(indexes[pos].strip(" ")) == str(aa1.strip(" "))][0]
                     belongs = True
                     if fragment_restrict is not None:
                      if fragment_restrict.has_key(pdb_chain_obj.chain):
                       if fragment_restrict.get(pdb_chain_obj.chain) is not None:
                        belongs=False
                        for interval in fragment_restrict.get(pdb_chain_obj.chain):
                            if int(residue_num1) >= int(interval[0]) and int(residue_num1) <= int(interval[1]): belongs=True
                            if belongs: break
                      else: continue
                     if belongs:
                        pair1=(str(aa1.strip(" ")),residue_num1)
                        singleton_patch.append(pair1)
                        for contact2 in contacts_obj.get_contacts():
                            aa2 = contact2._A_residue_obj.identifier
                            chain2 = contact2._A_chain
                            if contact2.get_contact_distance() > interface_distance: continue
                            if contact2.get_contact_distance() > radius: continue
                            # Check  the protein chain
                            if pdb_chain_obj.chain == chain2:
                              if aa1 != aa2:
                                  residue_num2 = 1 + [pos for pos in range(len(indexes)) if str(indexes[pos].strip(" ")) == str(aa2.strip(" "))][0]
                                  belongs = True
                                  if fragment_restrict is not None:
                                    if fragment_restrict.has_key(pdb_chain_obj.chain):
                                      if fragment_restrict.get(pdb_chain_obj.chain) is not None:
                                       belongs=False
                                       for interval in fragment_restrict.get(pdb_chain_obj.chain):
                                           if int(residue_num2) >= int(interval[0]) and int(residue_num2) <= int(interval[1]): belongs=True
                                           if belongs: break
                                    else: continue
                                  if belongs:
                                    pair2=(str(aa2.strip(" ")),residue_num2)
                                    # Get the residues
                                    pdb_residue_obj1 =  pdb_chain_obj.get_residue_by_identifier(aa1)
                                    pdb_residue_obj2 =  pdb_chain_obj.get_residue_by_identifier(aa2)
                                    # Get the distance and if it is smaller than the before set distance append the residue to the patch
                                    if pdb_residue_obj1.distance(pdb_residue_obj2,'min')[2] < distance:
                                        if pair1 not in patch:
                                            patch.append(pair1)
                                        if pair2 not in patch:
                                            patch.append(pair2)
                        if singleton_patch[0] not in patch:
                            patches.append(singleton_patch)
                            done.add(aa1)
                        # Set the length of the patch
                        if len(patch) > patch_size:
                            distance = distance - 1
                            big_patch = True
                            break
                        elif len(patch) <= patch_size and len(patch) != 0:
                            p = [x for x in sorted(patch,key=lambda x: x[-1])]
                            distance=distance_start
                            done.add(aa1)
                            patches.append(p)  

          # Remove duplicates
          patches_all = set(map(tuple, patches))
          patches_list = map(list, patches_all)
          patches_final = []
          for p1 in range(len(patches_list)):
              skip = False
              patch1=patches_list[p1]
              for p2 in range(len(patches_list)):
                  patch2=patches_list[p2]
                  if p1==p2:continue
                  #if set(patch1).issubset(set(patch2)) or set(patch1).issuperset(set(patch2)):
                  if set(patch1).issubset(set(patch2)):
                     skip = True
                     break
              if not skip:
                 patches_final.append(patch1)

          # Sort the patches final
          patches_final.sort(key =lambda x: int(x[0][1]))

          # store in dictionary
          patches_complete.setdefault(pdb_chain_id,patches_final)

    return patches_complete


def get_msa_obj(pdb_obj,triads_obj, x3dna_obj, potentials, radius, fragment_restrict, binding_restrict, split_potential, patches,  seq_threshold=None, verbose=False, specificity=False,allowed_sequences=None):

    # Initialize #
    binding_sites = {}
    all_sequence_patches_scaled_scores = {}
    msa_full={}
    msa_short={}
    msa_gapped={}

    
    # For each PDB chain... #
    for pdb_chain in sorted(potentials):
        patches_chain=patches.get(pdb_chain)
        # Get dinucleotide raw scores and binding site region #
        if verbose:sys.stdout.write("\t\t-- Get binding and scores for all amino-acids in chain %s...\n"%pdb_chain)
        residue_scores, binding_site, ss_exposure, native_aa_propensity = get_residue_scores_and_binding_site(pdb_obj,triads_obj,x3dna_obj, potentials, radius, fragment_restrict, binding_restrict, split_potential, pdb_chain, verbose)
        # Get sequence scaled scores #
        if verbose:sys.stdout.write("\t\t-- Score all sequence patches for chain %s ...\n"%pdb_chain)
        all_sequence_patches_scaled_scores.setdefault(pdb_chain, get_all_scaled_patches_scores( pdb_obj,pdb_chain,residue_scores, ss_exposure, native_aa_propensity, patches_chain, binding_site, split_potential, seq_threshold, dummy_dir, verbose,specificity,allowed_sequences))
        # Add binding site to binding sites #
        binding_sites.setdefault(pdb_chain, binding_site)

  
    # For each PDB chain... #
    for pdb_chain,binding_site in binding_sites.iteritems():
        if verbose:sys.stdout.write("\t\t-- Get MSA on binding of chain %s\n"%pdb_chain)
        # Initialize #
        pdb_chain_obj=pdb_obj.get_chain_by_id(pdb_chain)
        native_sequence = pdb_chain_obj.gapped_protein_sequence
        msa_short_obj  = PWM.pMSA()
        msa_gapped_obj = PWM.pMSA()
        msa_obj        = PWM.pMSA()
        # Add binding site length #
        if len(binding_site) <=0 : continue
        msa_short_obj.set_binding_site_length(len(binding_site))
        msa_gapped_obj.set_binding_site_length(len(native_sequence))
        msa_obj.set_binding_site_length(len(native_sequence))
        # For each sequence and patch... #
        if all_sequence_patches_scaled_scores.has_key(pdb_chain):
          for sequence_patch,score in all_sequence_patches_scaled_scores[pdb_chain].iteritems():
            sequence,patch=sequence_patch
            tmp_full_sequence   = list(native_sequence)
            tmp_short_sequence  = list('.' * len(sequence))
            tmp_gapped_sequence = list('.'*len(native_sequence))
            residues = list(sequence)
            patch_list=list(patch)
            ordered_patches=sorted(patch_list,key=lambda x: x[1])
            reorder={}
            reorder_i={}
            for i in range(len(tmp_short_sequence)):
                pair=patch_list[i]
                idx = pair[0]
                jdx = [j for j in range(len(ordered_patches)) if ordered_patches[j]==pair][0]
                reorder.setdefault(i,jdx)
                reorder_i.setdefault(jdx,i)
            short_sequence_table={}
            residue_in_binding=[]
            sequence_in_binding=[]
            for resnum,residx in binding_site.iterkeys():
                residue_in_binding.append(resnum)
                short_sequence_table.setdefault(resnum,False)
            for pair in list(patch):
                residue_num =pair[0]
                sequence_in_binding.append(pair[1]-1)
                if residue_num in residue_in_binding:
                   short_sequence_table[residue_num]=True
            for j in range(len(tmp_short_sequence)):
                i=reorder_i[j]
                pair= ordered_patches[j]
                residue_num=pair[0]
                if short_sequence_table[residue_num]:
                   tmp_short_sequence[j]=residues[i]
            residue_by_idx={}
            for aa,idx in zip(residues,sequence_in_binding):
                residue_by_idx.setdefault(idx,aa)
            for i in range(len(tmp_full_sequence)):
                if residue_by_idx.has_key(i):
                   tmp_full_sequence[i]=residue_by_idx[i]
                   tmp_gapped_sequence[i]=residue_by_idx[i]
            short_sequence  = ''.join(tmp_short_sequence)
            gapped_sequence = ''.join(tmp_gapped_sequence)
            full_sequence   = ''.join(tmp_full_sequence)
            msa_obj.add_sequence(full_sequence, score)
            msa_gapped_obj.add_sequence(gapped_sequence, score)
            msa_short_obj.add_sequence(short_sequence, score)
          #Add the object on the dictionary by chain
          msa_obj.set_pwm()
          msa_gapped_obj.set_pwm()
          msa_short_obj.set_pwm()
          msa_full.setdefault(pdb_chain,msa_obj)
          msa_gapped.setdefault(pdb_chain,msa_gapped_obj)
          msa_short.setdefault(pdb_chain,msa_short_obj)
    
    return msa_full,msa_gapped,msa_short

def get_residue_scores_and_binding_site(pdb_obj, triads_obj, x3dna_obj, potentials, radius, fragment_restrict, binding_restrict, split_potential,pdb_chain,verbose):
    '''
    '''
    # Initialize #
    nucleotides = list("ACGT")
    # Get original DNA sequence indexes
    basepairs=x3dna_obj.get_basepairs()
    dna_idx={}
    for basepair in basepairs.iterkeys():
        (fwd_pdb_chain, fwd_residue_num), (rev_pdb_chain, rev_residue_num) = basepairs.get(basepair)
        dna_idx.setdefault((fwd_pdb_chain,fwd_residue_num),int(basepair))
    if radius <= 0: radius=float(config.get("Parameters", "max_contact_distance"))
    # Get original protein sequence indexes
    aminoacids = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN","GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    residue_scores = {}
    binding_site = {}
    ss_exposure = {}
    native_aa_propensity = {}
    pdb_chain_obj=pdb_obj.get_chain_by_id(pdb_chain)
    sequence = pdb_chain_obj.gapped_protein_sequence
    indexes  = pdb_chain_obj.protein_idx.split(";")
    # For each patch... #
    for triad_obj in triads_obj.get_triads():
        a_oa, b_ob, distance, residue_A, residue_B = triad_obj.return_as_string().split(';')
        a, hydrophobicity, degree_of_exposure, secondary_structure = a_oa.split('-')
        chain, residue_num = residue_A.split('-')
        residue_num=str(residue_num.strip(" "))
        ss_exposure[residue_num] = (aminoacids3to1[a], degree_of_exposure, secondary_structure)
        a1 = aminoacids3to1[a]
        if chain != pdb_chain: continue
        if degree_of_exposure == 'E' and secondary_structure == 'H':
            native_aa_propensity[(residue_num)] = exposed_helix_propensities_dict[a1]
        elif degree_of_exposure == 'B' and secondary_structure == 'H':
            native_aa_propensity[(residue_num)] = buried_helix_propensities_dict[a1]
        elif degree_of_exposure == 'E' and secondary_structure == 'E':
            native_aa_propensity[(residue_num)] = exposed_b_strand_propensities_dict[a1]
        elif degree_of_exposure == 'B' and secondary_structure == 'E':
            native_aa_propensity[(residue_num)] = buried_b_strand_propensities_dict[a1]
    
    for triad_obj in triads_obj.get_triads():
        a_oa, b_ob, distance, residue_A, residue_B = triad_obj.return_as_string().split(';')
        chain, residue_num = residue_A.split('-')
        residue_num=str(residue_num.strip(" "))
        # Skip if not the right chain... #
        if chain != pdb_chain: continue
        # Skip if the environment could not be obtained for the amino acid
        if 'None' in a_oa: continue
        # Skip if the environment could not be obtained for the nucleotide
        if 'None' in b_ob: continue
        # Skip if not in radius contact
        dab = np.floor(float(distance))
        if dab > radius: continue
        # Skip if amino-acid is not in the restricted fragment... #
        residue_num1= 1 + [pos for pos in range(len(indexes)) if str(indexes[pos].strip(" ")) == str(residue_num.strip(" "))][0]
        if fragment_restrict is not None:
         if fragment_restrict.has_key(pdb_chain):
           if fragment_restrict.get(pdb_chain) is not None:
               belongs=False
               for interval in fragment_restrict.get(pdb_chain):
                   if int(residue_num1) >= int(interval[0]) and int(residue_num1) <= int(interval[1]): belongs=True
                   if belongs: break
               if not belongs: continue
         else: continue
        # Skip if nucleotide is not in the restricted binding... #
        bp_1F,bp_1R,bp_2F,bp_2R=residue_B.split(",")
        # Add dinucleotide to binding site #
        chain_dna,base_number_1=bp_1F.split("-")
        chain_dna,base_number_2=bp_2F.split("-")
        binucleotide = (dna_idx[(chain_dna,int(base_number_1))],dna_idx[(chain_dna,int(base_number_2))])
        if binding_restrict is not None:
               belongs=False
               for interval in binding_restrict:
                   if int(binucleotide[0]) >= interval[0] and int(binucleotide[1]) <= interval[1]: belongs=True
                   if belongs: break
               if not belongs: continue
        # Add aminoacid to the binding site #
        binding_site.setdefault((residue_num,residue_num1), []).append(triad_obj)
        # The following arguments are defined in a paper (from pwm.py script)
        a_native, hydrophobicity, degree_of_exposure, secondary_structure = a_oa.split('-')
        b, twobases, dna_strand, dna_groove, dna_chemical_group = b_ob.split('-')
        for aa in aminoacids:
           residue_scores.setdefault((aminoacids3to1[aa],(residue_num,residue_num1)), 0.0)
           # Change the hydrophobicity according to the new residue. We only change the hydrophobicity, but not the other parameters, this is threading#
           if aminoacids_polarity_boolean[aminoacids3to1[aa]]:
                hydrophobicity = 'P'
           else:
                hydrophobicity = 'N' 
           a1 = aminoacids3to1[aa]
           # Apply potential... #   
           residue_scores[(a1,(residue_num,residue_num1))] += zscore_by_dinucleotide(aa,a_oa,b_ob,potentials,chain,dab,split_potential) 

    return residue_scores, binding_site, ss_exposure, native_aa_propensity

def zscore_by_dinucleotide(aa,a_oa_triad,b_ob_triad,potentials,chain,dab,split_potential):
    '''
    '''
    # Initialize #
    nucleotides = list("ACGT")
    scores_list=[]
    zscore=0
    # triad data
    a_native, hydrophobicity, degree_of_exposure, secondary_structure = a_oa_triad.split('-')
    b_triad, twobases, dna_strand, dna_groove, dna_chemical_group = b_ob_triad.split('-')
    a_oa_triad      = '%s-%s-%s-%s' %(aa,hydrophobicity, degree_of_exposure, secondary_structure)
    oa_triad        = '%s-%s-%s' %(hydrophobicity, degree_of_exposure, secondary_structure)
    ob_triad        = '%s-%s-%s-%s' %(twobases, dna_strand, dna_groove, dna_chemical_group)
    a_b_triad       = "%s;%s" % (aa, b_triad)
    a_b_oa_ob_triad = "%s;%s" % (a_oa_triad, b_ob_triad)
    oa_ob_triad     = "%s;%s" % (oa_triad, ob_triad)
    # score triad
    score=0
    if split_potential == "3d":
     if potentials[chain].get_score(split_potential, distance=dab) is not None:
        score= potentials[chain].get_score(split_potential, distance=dab)
    if split_potential == "3dc":
     if potentials[chain].get_score(split_potential, key=oa_ob_triad, distance=dab) is not None:
        score= potentials[chain].get_score(split_potential, key=oa_ob_triad, distance=dab)
    if split_potential == "s3dc":
     if potentials[chain].get_score(split_potential, key=a_b_oa_ob_triad, distance=dab) is not None:
        score= potentials[chain].get_score(split_potential, key=a_b_oa_ob_triad, distance=dab)
    if split_potential == "s3dc_dd":
     if potentials[chain].get_score(split_potential, key=a_b_oa_ob_triad, distance=dab) is not None:
        score= potentials[chain].get_score(split_potential, key=a_b_oa_ob_triad, distance=dab)
    if split_potential == "s3dc_di":
     if potentials[chain].get_score(split_potential, key=a_b_oa_ob_triad) is not None:
        score= potentials[chain].get_score(split_potential, key=a_b_oa_ob_triad)
    if split_potential == "pair":
     if potentials[chain].get_score(split_potential, key=a_b_triad, distance=dab) is not None:
        score= potentials[chain].get_score(split_potential, key=a_b_triad, distance=dab)
    # all other dinucleotides data
    dinucleotide_set=set()
    for n1 in nucleotides:
      for n2 in nucleotides:
        dinucleotide_set.add(n1+n2)
    for dinucleotide in dinucleotide_set:
        newbases="".join([nitrogenous_bases[nucleotide] for nucleotide in dinucleotide])
        b=dinucleotide
        a_oa = '%s-%s-%s-%s' %(aa,hydrophobicity, degree_of_exposure, secondary_structure)
        oa = '%s-%s-%s' %(hydrophobicity, degree_of_exposure, secondary_structure)
        ob = '%s-%s-%s-%s' %(newbases, dna_strand, dna_groove, dna_chemical_group)
        b_ob = '%s-%s-%s-%s-%s' %(dinucleotide,newbases, dna_strand, dna_groove, dna_chemical_group)
        a_b = "%s;%s" % (aa, dinucleotide)
        a_b_oa_ob = "%s;%s" % (a_oa, b_ob)
        oa_ob = "%s;%s" % (oa, ob)
        #  potential... #    
        if split_potential == "3d":
         if potentials[chain].get_score(split_potential, distance=dab) is not None:
            scores_list.append(potentials[chain].get_score(split_potential, distance=dab))
        if split_potential == "3dc":
         if potentials[chain].get_score(split_potential, key=oa_ob, distance=dab) is not None:
            scores_list.append(potentials[chain].get_score(split_potential, key=oa_ob, distance=dab))
        if split_potential == "s3dc":
         if potentials[chain].get_score(split_potential, key=a_b_oa_ob, distance=dab) is not None:
            scores_list.append(potentials[chain].get_score(split_potential, key=a_b_oa_ob, distance=dab))
        if split_potential == "s3dc_dd":
         if potentials[chain].get_score(split_potential, key=a_b_oa_ob, distance=dab) is not None:
            scores_list.append(potentials[chain].get_score(split_potential, key=a_b_oa_ob, distance=dab))
        if split_potential == "s3dc_di":
         if potentials[chain].get_score(split_potential, key=a_b_oa_ob) is not None:
            scores_list.append(potentials[chain].get_score(split_potential, key=a_b_oa_ob))
        if split_potential == "pair":
         if potentials[chain].get_score(split_potential, key=a_b, distance=dab) is not None:
            scores_list.append(potentials[chain].get_score(split_potential, key=a_b, distance=dab))  
    #get zscore     
    average,sigma=average_sigma(scores_list) 
    if sigma>0: zscore=(score-average)/sigma
    #return
    return zscore
        
    
def scale_all_patches_scores(all_patches_scores,ss_exposure,native_aa_propensity,split_potential,threshold=None,apply_structure=True):
    # Scale the SEQUENCES scores #
    max_scores = {}
    min_scores = {}
    all_patches_scaled_scores={}
    # Modify the sign 
    if split_potential == 's3dc' or split_potential == 's3dc_dd' or split_potential == 'pair':
        # For each patch ...#
        for key in all_patches_scores.iterkeys():
            all_patches_scores[key] *= -1
    # Get Max and min per patch... #
    for (aminoacid_seq, patch) in all_patches_scores.iterkeys():
        # Get max and min scores for scaling #
        if max_scores.has_key(patch):
           if all_patches_scores[(aminoacid_seq, patch)] > max_scores[patch]: max_scores[patch] = all_patches_scores[(aminoacid_seq, patch)]
        else:
           max_scores.setdefault(patch, all_patches_scores[aminoacid_seq, patch])
        if min_scores.has_key(patch):
           if all_patches_scores[(aminoacid_seq, patch)] < min_scores[patch]:min_scores[patch] = all_patches_scores[(aminoacid_seq, patch)]
        else:
           min_scores.setdefault(patch, all_patches_scores[aminoacid_seq, patch])

    # Modify the scores by structural weight and scale the scores

    for seq_patch in all_patches_scores.iterkeys():
        aminoacid_seq , patch = seq_patch
        residues = list(aminoacid_seq)
        patches  = list(patch)
        structure_weight = 1
        if apply_structure:
          for pair in patches:
            number=pair[0]
            position=int(pair[1])
            for residue in residues:
                if ss_exposure[number] == ('E', 'H'):
                    structure_weight *= (1 + ((exposed_helix_propensities_dict[residue] - native_aa_propensity[number])/20))
                elif ss_exposure[number] == ('B', 'H'):
                    structure_weight *= (1 + ((buried_helix_propensities_dict[residue] - native_aa_propensity[number])/20))
                elif ss_exposure[number] == ('E', 'E'):
                    structure_weight *= ( 1 + ((exposed_b_strand_propensities_dict[residue] - native_aa_propensity[number])/20))
                elif ss_exposure[number] == ('B', 'E'):
                    structure_weight *= ( 1 + ((buried_b_strand_propensities_dict[residue] - native_aa_propensity[number])/20))
                else:
                    pass
        scaled_score = scale(all_patches_scores[seq_patch], max_score=max_scores[patch], min_score=min_scores[patch]) * structure_weight
        if threshold == None:
           if max_scores[patch] != min_scores[patch]: all_patches_scaled_scores.setdefault(seq_patch, scaled_score)
        else:
           if scaled_score >= threshold: all_patches_scaled_scores.setdefault(seq_patch, scaled_score)
    return all_patches_scaled_scores

def combine_patches(patches_scores,residue_scores,ss_exposure,  binding_site,verbose,specificity):
    # Initialize #
    aminoacids = list("ACDEFGHIKLMNPQRSTVWY")

    # Algorithm to convert into non-overlapping patches
    x={}
    for seq_patch, score in patches_scores.iteritems():
        sequence = seq_patch[0]
        patch  = seq_patch[1]
        x.setdefault(patch,set()).add(sequence)
    z = x.copy()
    done = False
    while not done:
        # create working dictionary y
        y = dict()
        y = z.copy()
        # start again the construction of z
        # empty z dictionary to be reconstructed
        z = dict()
        found = []
        condition = True
        done = True
        for patch_p in y.keys():
            if patch_p in found: continue
            if condition == True:
                patch_p_is_unique=True
                for patch_q in y.keys():
                    if patch_q in found: continue
                    if patch_p != patch_q:
                        if len(set(patch_p).intersection(set(patch_q))) > 0:
                            #CONDITION TO CONTINUE
                            # we are not done yet because one intersection was found
                            done = False
                            # is the first intersection found, skip next intersections
                            condition = False
                            # Append patches to found #
                            found.append(patch_p)
                            found.append(patch_q)
                            patch_p_is_unique=False
                            # DEFINE new "patches" #
                            #intersection
                            intersection_pq = set(patch_p).intersection(set(patch_q))
                            intersection_pq = sorted(intersection_pq, key=lambda x: x[1])
                            intersection_pq = list(intersection_pq)
                            if len(intersection_pq)>0:
                               intersection_pq_key = tuple(intersection_pq)
                            else:
                               intersection_pq_key = None
                            #p-intersection
                            patch_p_unique=[]
                            for n in patch_p:
                                if n in intersection_pq:continue
                                patch_p_unique.append(n)
                            if len(patch_p_unique)>0:
                                patch_p_unique_key = tuple(patch_p_unique)
                            else:
                                patch_p_unique_key = None
                            #q-intersection
                            patch_q_unique=[]
                            for n in patch_q:
                                if n in intersection_pq:continue
                                patch_q_unique.append(n)
                            if len(patch_q_unique)>0:
                                patch_q_unique_key = tuple(patch_q_unique)
                            else:
                                patch_q_unique_key = None
                            # Sequences for New patches in intersection #
                            if intersection_pq_key is not None:
                              for sequence_p in y[patch_p]:
                                pairs=zip(patch_p,list(sequence_p))
                                tmp_sequence=""
                                for item in pairs:
                                    if item[0] in intersection_pq:
                                        tmp_sequence += item[1]
                                z.setdefault(intersection_pq_key, set()).add(tmp_sequence)
                              for sequence_q in y[patch_q]:
                                pairs=zip(patch_q,list(sequence_q))
                                tmp_sequence=""
                                for item in pairs:
                                    if item[0] in intersection_pq:
                                        tmp_sequence += item[1]
                                z.setdefault(intersection_pq_key, set()).add(tmp_sequence)
                            # Sequences for New patches in p-intersection#
                            if patch_p_unique_key is not None:
                              for sequence_p in y[patch_p]:
                                pairs=zip(patch_p,list(sequence_p))
                                tmp_sequence=""
                                for item in pairs:
                                    if item[0] in patch_p_unique:
                                        tmp_sequence += item[1]
                                z.setdefault(patch_p_unique_key, set()).add(tmp_sequence)
                            # Sequences for New patches in q-intersection#
                            if patch_q_unique_key is not None:
                              for sequence_q in y[patch_q]:
                                pairs=zip(patch_q,list(sequence_q))
                                tmp_sequence=""
                                for item in pairs:
                                    if item[0] in patch_q_unique:
                                        tmp_sequence += item[1]
                                z.setdefault(patch_q_unique_key, set()).add(tmp_sequence)
                            break
                # after checking with all patches patch p has no interscetions
                if patch_p_is_unique:
                    if len(patch_p)>1:
                        patch_p_key = tuple(patch_p)
                    elif len(patch_p)==1:
                        patch_p_key = patch_p[0]
                    else:
                        patch_p_key = None
                    if patch_p_key is not None:
                       patch_p_key=tuple(patch_p)
                       for sequence_p in y[patch_p]:
                        z.setdefault(patch_p_key, set()).add(sequence_p)
            # we found already an intersection and we keep copying the rest of patches
            else:
                if len(patch_p)>1:
                    patch_p_key = tuple(patch_p)
                elif len(patch_p)==1:
                    patch_p_key = patch_p[0]
                else:
                    patch_p_key = None
                if patch_p_key is not None:
                  patch_p_key=tuple(patch_p)
                  for sequence_p in y[patch_p]:
                    z.setdefault(patch_p_key, set()).add(sequence_p)



    # Merge patches and prepare scores
    patches_combined={}
    product="1 X"
    print " "
    for patch,sequence_set in z.iteritems():
       if len(patches_combined.keys())>0:
          if verbose: print "\t\t\tSIZE",len(patches_combined.keys()), " =  ",product
          product+= " %d X"%len(sequence_set)
          tmp_patches_combined={}
          for seq_patch in patches_combined.keys():
              sequence_combined= seq_patch[0]
              patch_combined   = list(seq_patch[1])
              adding_patch     = list(patch)
              patch_combined.extend(adding_patch)
              new_patch=tuple(patch_combined)
              for aminoacid_seq in sequence_set:
                  new_sequence=sequence_combined+aminoacid_seq
                  tmp_patches_combined.setdefault((new_sequence,new_patch), 0.0)
          patches_combined={}
          for key,value in tmp_patches_combined.iteritems():
              patches_combined.setdefault(key,value)
       else:
          product+= " %d X"%len(sequence_set)
          for aminoacid_seq in sequence_set:
              patches_combined.setdefault((aminoacid_seq, patch), 0.0)
    #Release memory
    z={}
    y={}
    x={}
    tmp_patches_combined={}
    if verbose: print "\t\t\tSIZE",len(patches_combined.keys()), " =  ",product.rstrip("X")

    # Sum the scores of each patch # 
    if specificity:
      if verbose: print "\t\t-- Calculate the distribution of amino-acids ..."
      distribution_aa={}
      for (aminoacid_seq, patch) in patches_combined.iterkeys():
        aminoacid_seq_l = list(aminoacid_seq)
        patch_l = list(patch)
        skip=False
        for i in range(len(aminoacid_seq_l)):
                aa=aminoacid_seq_l[i]
                pair=tuple(patch_l[i])
                #Check the sequence is in the binding site, otherwise allow only the native aminoacid
                if not binding_site.has_key(pair):
                   if ss_exposure.has_key(pair[0]):
                      native_aa = ss_exposure.get(pair[0])
                      if aa != native_aa: skip=True
                #Add residue on the set
                if not skip: distribution_aa.setdefault(pair,[]).append(aa)


    if verbose: print "\t\t-- Score the sequences of combined patches ..."
    for (aminoacid_seq, patch) in patches_combined.iterkeys():
        aminoacid_seq_l = list(aminoacid_seq)
        patch_l = list(patch)
        tmp_score=0
        skip=False
        for i in range(len(aminoacid_seq_l)):
                aa=aminoacid_seq_l[i]
                pair=tuple(patch_l[i])
                #Check the sequence is in the binding site, otherwise allow only the native aminoacid
                if not binding_site.has_key(pair):
                   if ss_exposure.has_key(pair[0]):
                      native_aa = ss_exposure.get(pair[0])
                      if aa != native_aa: skip=True
                #Add interaction of the residue
                if residue_scores.has_key((aa,pair)):
                   if specificity:
                      if distribution_aa.has_key(pair):
                        scores_pair  =[residue_scores.get((bb,pair)) for bb in distribution_aa.get(pair) if residue_scores.has_key((bb,pair)) ]
                        average,sigma=average_sigma(scores_pair)
                        if sigma>0: tmp_score += (residue_scores.get((aa,pair)) - average)/sigma
                   else:
                       tmp_score += residue_scores.get((aa,pair))
        if not skip: patches_combined[(aminoacid_seq, patch)] += tmp_score

    # Return scores
    return patches_combined

def average_sigma(l):
    import math
    average=0.0
    sigma=0.0
    if len(l) >0:
       average=sum(l)/len(l)
       sigma  =math.sqrt(sum([(x-average)*(x-average) for x in l])/len(l))
    return average,sigma

def get_all_scaled_patches_scores(pdb_obj,pdb_chain,residue_scores, ss_exposure, native_aa_propensity, patches, binding_site, split_potential, threshold=None, dummy_dir="/tmp",verbose=False,specificity=False,allowed_sequences=None):
    """
    """
    # Initialize #
    all_patches_scores = {}
    all_patches_scaled_scores = {}
    if threshold is not None:
       first_threshold  = threshold[0]
       second_threshold = threshold[1]
    else:
       first_threshold  = None
       second_threshold = None
    aminoacids = list("ACDEFGHIKLMNPQRSTVWY")
    # Set a default dictionary for the Scores of the pacthes #
    if allowed_sequences is not None:
      # Scale is not weighted by secondary structure
      apply_structure=False
      # Unify the pactches with respect to the allowed sequences
      patch_sum=set()
      for patch_l in patches:
        for pair in patch_l:
            patch_sum.add(pair)
      patch_global=[x for x in sorted(patch_sum, key=lambda x:int(x[1]))]
      patch = tuple(patch_global)
      size=0
      for aminoacid_seq in allowed_sequences:
        if len(aminoacid_seq) == len(patch_global):
            size+=1
            all_patches_scores.setdefault((aminoacid_seq, patch), 0.0)
      if verbose:sys.stdout.write("\t\t\t-- pre-requested memory by all patches: %d floats\n"%size)
    else:
      # Scale is weighted by secondary structure
      apply_structure=True
      # Construct dictionary for all patches
      size=0
      for patch_l in patches:
        patch = tuple(patch_l)
        sizep=0
        for aminoacid in itertools.product(aminoacids, repeat=len(patch)):
            aminoacid_seq = ''.join(aminoacid)
            size+=1
            sizep+=1
            # Set a default dictionary for the scaled scores of the patches #
            all_patches_scores.setdefault((aminoacid_seq, patch), 0.0)
        if verbose:sys.stdout.write("\t\t\t-- pre-requested memory by patch %s : %d floats\n"%(patch,sizep))
      if verbose:sys.stdout.write("\t\t\t-- pre-requested memory by all patches: %d floats\n"%size)
    # Sum the scores of each patch # 
    if verbose:sys.stdout.write("\t\t\t-- calculate scores for all interface patches\n")
    size=0
    for (aminoacid_seq, patch) in all_patches_scores.iterkeys():
        size += 1
        aminoacid_seq_l = list(aminoacid_seq)
        patch_l = list(patch)
        tmp_score=0
        skip=False
        for i in range(len(aminoacid_seq_l)):
                aa=aminoacid_seq_l[i]
                pair=tuple(patch_l[i])
                #Check the sequence is in the binding site, otherwise allow only the native aminoacid
                if not binding_site.has_key(pair):
                   if ss_exposure.has_key(pair[0]):
                      native_aa = ss_exposure.get(pair[0])
                      if aa != native_aa: skip=True
                #Add interaction of the residue
                if residue_scores.has_key((aa,pair)):
                   tmp_score += residue_scores.get((aa,pair))
        if not skip: all_patches_scores[(aminoacid_seq, patch)] += tmp_score
    if verbose:sys.stdout.write("\t\t\t\t-- memory requested: %d floats\n"%size)


    #scale all_patches and reduce the set
    if verbose:sys.stdout.write("\t\t\t-- scale the scores for all interface patches\n")
    ###WARNING SIGNS CHANGE
    all_patches_scaled_scores=scale_all_patches_scores(all_patches_scores,ss_exposure,native_aa_propensity,split_potential,first_threshold,apply_structure)
    #Release memory
    all_patches_scores={}

    #Combine all patches with only the selected sequences
    if verbose:sys.stdout.write("\t\t\t-- combine all interface patches in sequences\n")
    ###WARNING SIGNS IS RECOVERED
    ###WARNING SPECIFICITY IMPLIES Z-SCORE WITH ITS SIGN BACK
    all_patches_scores=combine_patches(all_patches_scaled_scores,residue_scores,ss_exposure,binding_site,verbose,specificity)
    #Release memory
    all_patches_scaled_scores={}

    if verbose:
       print " "
       print "\t\t\tPATCHES"
       print " "
       chain=pdb_chain
       pdb_chain_obj=pdb_obj.get_chain_by_id(pdb_chain)
       seq= pdb_chain_obj.gapped_protein_sequence
       count=0
       print chain,"%3d"%count,seq
       done=set()
       for sequence_patch, score in all_patches_scores.iteritems():
         sequence,patch = sequence_patch
         residues = list(sequence)
         patch_list = list(patch)
         if patch not in done:
            done.add(patch)
            dummy_seq=list("." * len(seq))
            for pair in patch_list:
                res=pair[0]
                idx=int(pair[1])
                dummy_seq[idx-1]=seq[idx-1]
            count += 1
            print chain,"%3d"%count,"".join(dummy_seq)
       print " "


    #scale selected patches and sequences and reduce the set
    if verbose:sys.stdout.write("\t\t\t-- scale the scores for all interface sequences\n")
    if verbose:sys.stdout.write("\t\t\t\t--  memory requested: %d floats\n"%len(all_patches_scores.keys()))
    all_patches_scaled_scores=scale_all_patches_scores(all_patches_scores,ss_exposure,native_aa_propensity,split_potential,second_threshold,apply_structure)
    #Release memory
    all_patches_scores={}


    #Filter out sequences with wrong structure or not allowed
    if allowed_sequences is not None:
      if verbose:sys.stdout.write("\t\t\t-- Filter out sequences that are not in the allowed file\n")
      all_sequence_patches_scaled_scores=filter_out_by_selection(all_patches_scaled_scores,binding_site,allowed_sequences)
      if verbose:sys.stdout.write("\t\t\t\t--  memory requested: %d floats\n"%len(all_sequence_patches_scaled_scores.keys()))
    else:
      if verbose:sys.stdout.write("\t\t\t-- Filter out sequences incompatible with structure\n")
      all_sequence_patches_scaled_scores=filter_out_by_PROSA(all_patches_scaled_scores,pdb_obj,pdb_chain,dummy_dir)
      if verbose:sys.stdout.write("\t\t\t\t--  memory requested: %d floats\n"%len(all_sequence_patches_scaled_scores.keys()))
    #Release memory
    all_patches_scaled_scores={}
   


    return all_sequence_patches_scaled_scores

def filter_out_by_selection(all_patches_scaled_scores,binding_site,allowed_sequences):

    #Initialize
    all_sequence_patches_scaled_scores={}
    #Check allowed sequences
    for sequence_patch, score in all_patches_scaled_scores.iteritems():
        sequence,patches = sequence_patch
        residues = list(sequence)
        patches_list = list(patches)
        tmp_short_sequence  = list('.' * len(sequence))
        ordered_patches=sorted(patches_list,key=lambda x: x[1])
        reorder={}
        reorder_i={}
        for i in range(len(tmp_short_sequence)):
                pair=patches_list[i]
                idx = pair[0]
                jdx = [j for j in range(len(ordered_patches)) if ordered_patches[j]==pair][0]
                reorder.setdefault(i,jdx)
                reorder_i.setdefault(jdx,i)
        short_sequence_table={}
        residue_in_binding=[]
        sequence_in_binding=[]
        for resnum,residx in binding_site.iterkeys():
                residue_in_binding.append(resnum)
                short_sequence_table.setdefault(resnum,False)
        for pair in list(patches):
                residue_num =pair[0]
                sequence_in_binding.append(pair[1]-1)
                if residue_num in residue_in_binding:
                   short_sequence_table[residue_num]=True
        for j in range(len(tmp_short_sequence)):
                i=reorder_i[j]
                pair= ordered_patches[j]
                residue_num=pair[0]
                if short_sequence_table[residue_num]:
                   tmp_short_sequence[j]=residues[i]
        short_sequence  = ''.join(tmp_short_sequence)
        if short_sequence in allowed_sequences:
           all_sequence_patches_scaled_scores[(sequence,patches)]=score
    #return selected set
    return all_sequence_patches_scaled_scores
    
def filter_out_by_PROSA(all_patches_scaled_scores,pdb_obj,pdb_chain,dummy="/tmp"):
   
    #Initialize
    all_sequence_patches_scaled_scores={}
    if not os.path.exists(dummy): 
       os.makedirs(dummy)
    dummy_dir = os.path.join(dummy, "prosa_"+str(os.getpid()))
    if not os.path.exists(dummy_dir):
       os.makedirs(dummy_dir)
    prosa_dir = config.get("Paths","ProSaData")
    prosa_exe = config.get("Paths","prosa2003")
    pcb = os.path.abspath(os.path.join(prosa_dir,"pII3.0.pair-cb"))
    scb = os.path.abspath(os.path.join(prosa_dir,"pII3.0.surf-cb"))
    pca = os.path.abspath(os.path.join(prosa_dir,"pII3.0.pair-ca"))
    sca = os.path.abspath(os.path.join(prosa_dir,"pII3.0.surf-ca"))

    # Make dummy pdb
    pdb_file = os.path.join(dummy_dir,"dummy_"+pdb_chain+".pdb")
    chain_obj= pdb_obj.get_chain_by_id(pdb_chain)
    pdb_dummy=PDB()
    pdb_dummy.add_chain(chain_obj)
    pdb_dummy.write(pdb_file)

    # Check the tertiary structure with Prosa # CHECK IT DOESN'TWORK
    for sequence_patch, score in all_patches_scaled_scores.iteritems():
        sequence,patches = sequence_patch
        residues = list(sequence)
        patches_list = list(patches)
        residues_num=[]
        residues_idx=[]
        for pair in patches_list:
            residues_num.append(pair[0])
            residues_idx.append(pair[1])
        
        prosa_file = os.path.abspath(os.path.join(dummy_dir,"prosa_"+hashlib.sha224(sequence).hexdigest()+".cmd"))
        prosa = open(prosa_file, 'w')
        prosa.write('pair potential %s pcb\n'%pcb)
        prosa.write('surface potential %s scb\n'%scb)
        prosa.write('pair potential %s pca\n'%pca)
        prosa.write('surface potential %s sca\n'%sca)
        prosa.write('use potential * pcb scb pca sca\n')
        prosa.write('read pdb %s wt\n' %pdb_file)
        counter = 1
        
        tmp = 'wt'
        substitutions=" "
        for num, residue in zip(residues_idx, residues):
            prosa.write('mutate sequence %s %s %s m%s\n' %(tmp, str(num), residue, str(counter)))
            substitutions += str(num)+" "
            tmp = 'm'+str(counter)
            counter += 1

        counter = counter - 1
        prosa.write('winsize * 10\n')
        prosa.write('analyse energy *\n')
        prosa.write('diff m%s wt diff\n'%str(counter))
        prosa.write('draw * * 0\n')
        prosa.write('color pair wt cyan\n')
        prosa.write('color pair diff green\n')
        prosa.write('color pair m%s yellow\n'%str(counter))
        prosa.write('draw pair wt 1\n')
        prosa.write('draw pair m%s 1\n'%str(counter))
        prosa.write('draw pair diff 1\n')
        
        prosa.write('print energy diff diff.ene\n')
        prosa.write('print energy wt wt.ene\n')
        prosa.write('print energy m%s m%s.ene\n'%(str(counter), str(counter)))
        prosa.write('exit\n')
        prosa.close()

        wild_type_energy = 'wt.ene.ana'
        mutant_energy = 'm%s.ene.ana' %str(counter)
        diff_energy = 'diff.ene.ana'
        wt_ene_file = os.path.abspath(os.path.join(dummy_dir,wild_type_energy))
        mt_ene_file = os.path.abspath(os.path.join(dummy_dir,wild_type_energy))
        error_counter = 0

        try:
          process = subprocess.check_output([prosa_exe, prosa_file], stderr=subprocess.STDOUT)
          #os.system('%s %s' %(prosa_exe,prosa_file))
          wt_ene = open(wt_ene_file, 'r').readlines()[2:]
          mt_ene = open(mt_ene_file, 'r').readlines()[2:]
          pair_wt={}
          pair_mt={}
          for linewt in wt_ene:
                    wt = linewt.rstrip('\n').split()
                    pair_wt.setdefault(wt[0],float(wt[1]))
          for linemt in mt_ene:
                    mt = linemt.rstrip('\n').lstrip(' ').split('   ')
                    pair_mt.setdefault(mt[0],float(mt[1]))
          for reside in residue_idx:
            if pair_mt.has_key(residue) and pair_wt.has_key(residue):
              if pair_mt.get(residue) > pair_wt.get(residue):
                   if pair_mt.get(residue) > 0 and pair_wt.get(residue) < 0:
                      if pair_mt.get(residue) > 1.00:
                        error_counter += 1
                   else:
                      diff = (pair_mt.get(residue)/pair_wt.get(residue))
                      if diff > 2:
                         error_counter += 1
        except:
          #sys.stdout.write("PROSA analysis failed, continue without checking the structure on positions %s (substitute by %s)\n"%(substitutions,sequence))
          pass

        if error_counter < 4:
                all_sequence_patches_scaled_scores[(sequence,patches)]=score
        
        if os.path.exists(prosa_file):      os.remove(prosa_file)
        if os.path.exists(wt_ene_file):     os.remove(wt_ene_file)
        if os.path.exists(mt_ene_file):     os.remove(mt_ene_file)
        
    shutil.rmtree(dummy_dir)

    return all_sequence_patches_scaled_scores


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
  if (max_score - min_score) != 0.0:
    return min_scaled_score + (max_scaled_score - min_scaled_score) * (score - min_score) / (max_score - min_score)
  else:
    return 0.0


def parse_options():
    '''
    This function parses the command line arguemnts and returns an optparse object
    '''
    parser = optparse.OptionParser('python TF_design.py -i input_pdb_file --pbm=pbm_dir --pdb=pdb_dir [-o output_file ] [-a -f -p --taylor -s potential -b -t -k -r --fragment a_A-b_A;c_B-d_B --binding a-b;c-d  --first_threshold 0.95 --second_threshold 0.95 --radius  0.0 -v ]')

    parser.add_option('--dummy', default='/tmp/', action='store', type='string', dest='dummy_dir', help='Dummy directory (default = /tmp/', metavar='{directory}')
    parser.add_option('-i', action='store', type='string', dest='input_pdb_file', help='PDB file', metavar='{filename}')
    parser.add_option('-o', default='/tmp/', action="store", type="string", dest="output_file", help="Output dir(default = 'output')", metavar="{string}")
    parser.add_option('-v', '--verbose', default='false', action='store_true', dest='verbose', help='Verbose mode (default = False')
    parser.add_option("--pbm", action="store", type="string", dest="pbm_dir", default=None, help="PBM directory (i.e. output dir from pbm.py)", metavar="{directory}")
    parser.add_option("--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (i.e. output dir from pdb.py)", metavar="{directory}")
    parser.add_option("-k","--known", default=False, action="store_true", dest="known", help="The name is of a known PDB file, with format 'code_chain' (default = False)")
    parser.add_option("-r","--reset", default=False, action="store_true", dest="reset", help="Clean the sequences of the original MSA and reset them by a random selection in accordance with the PWM (default = False)")
    parser.add_option("--fragment", default=None, action="store", type="string", dest="fragment_restrict", help="Fragment of protein to apply the potential. Format is 'a_A-b_A;c_B-d_B': two regions between residues a-b and c-d from chain A and B. (Default is None it applies to all amino-acids)")
    parser.add_option("--binding", default=None, action="store", type="string", dest="binding_restrict", help="Binding site of DNA to apply the potential. Format is 'a-b;c-d': two regions between residues a-b and c-d of the forward chain (first in PDB). (Default is None it applies to all nucleotides)")
    parser.add_option("--sequences", default=None, action="store", type="string", dest="allowed_sequences" , help="File of interface sequences allowed. This reduces the number of accepted sequences. It forces to skip the structure control")

    group = optparse.OptionGroup(parser, "Statistical potentials", "Select your statistical potentials of choice. By default it uses S3DCdd general potential derived from PDB (the simplest one). In \"--auto\" mode, the program uses S3DCdd family potentials derived from both PDB and/or PBM data and/or approached by Taylor as selected in Potentials configuration file. In case family potentials cannot be applied, the program uses general potentials derived from both PDB and PBM data and approached by Taylor. \"-a\" option overrides options \"-f\", \"-p\" and \"-t\".")
    group.add_option("-a", "--auto", default=False, action="store_true", dest="auto_mode", help="Automate the selection of statistical potentials (default = False)")
    group.add_option("-f", "--family", default=False, action="store_true", dest="family_potentials", help="Use family potentials (default = False)")
    group.add_option("-p", default=False, action="store_true", dest="pbm_potentials", help="Use potentials derived from both PBM + PDB data (default = False)")
    group.add_option("-s", default="s3dc_dd", action="store", type="string", dest="split_potential", help="Split-potential to be used (3d, 3dc, s3dc, s3dc_dd, s3dc_di, pair; default = s3dc_dd)", metavar="{string}")
    group.add_option("-m", "--pmf", default=False, action="store_true", dest="pmf", help="Use of raw mean-force potentials with no Z-scoring (default = False)")
    group.add_option("--first_threshold", action="store", type="float",default=0.95, dest="first_score_threshold", help="Threshold on the scaled score to consider positive short fragments (default = 0.95)", metavar="{float}")
    group.add_option("--second_threshold", action="store", type="float",default=0.95, dest="second_score_threshold", help="Threshold on the scaled score to consider positive full sequences (default = 0.95)", metavar="{float}")
    group.add_option("--taylor", default=False, action="store_true", dest="taylor_approach", help="Approach PMF by Taylor (default = False)")
    group.add_option("-b", "--bins", default=False, action="store_true",  dest="bins_approach", help="Computate the potentials: by bins (if selected) or accumulative (default).", metavar="{boolean}")
    group.add_option("--specific", default=False, action="store_true",  dest="specificity", help="Normalize the scores of aminoacids to increase the specific interactions (not used by default).", metavar="{boolean}")
    group.add_option("--file",default=None, action="store", type="string", dest="potential_file", help="Use potentials from specific file (default = None)", metavar="{string}")
    group.add_option("--radius",default=0, action="store", type="float", dest="radius", help="Maximum contact distance to calculate interactions (default=0 implies the use of 'max_contact_distance' from configuration", metavar="{string}")
    parser.add_option_group(group)

    (options, args) = parser.parse_args()

    if options.input_pdb_file is None or (options.pdb_dir is None and options.potential_file is None):
        parser.error('Missing argmuents: type \'-h\' for help')

    if not re.search("^3d|3dc|s3dc$|^s3dc_dd$|^s3dc_di|pair$", options.split_potential):
        parser.error("incorrect value for -s argument: type option \"-h\" for help")

    return options

###############
#  MAIN   #
###############

if __name__ == "__main__":
    # Arguments & Options #
    options = parse_options()

    # Initialize #
    seq_threshold=[]
    seq_threshold.append(options.first_score_threshold)
    seq_threshold.append(options.second_score_threshold)
    if not os.path.exists(options.dummy_dir): os.makedirs(options.dummy_dir)
    dummy_dir=os.path.abspath(options.dummy_dir)
    pbm_dir=None
    if options.pbm_dir is not None: pbm_dir = os.path.abspath(options.pbm_dir)
    pdb_dir=None
    if options.pdb_dir is not None: pdb_dir = os.path.abspath(options.pdb_dir)
    families = {}
    if options.verbose:sys.stdout.write("\t-- Check families...\n")
    if os.path.exists(os.path.join(options.pdb_dir, "families.txt")):
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
    radius=options.radius

    # Get PDB object #
    if options.verbose:sys.stdout.write('\t-- Get protein %s... \n' %options.input_pdb_file)
    pdb_obj = PDB(options.input_pdb_file) 

    # Get DSSP Object #
    if options.verbose:sys.stdout.write('\t\t-- calculate secondary structure... \n')
    dssp_obj = dssp.get_dssp_obj(os.path.abspath(options.input_pdb_file))

    # Get X3DNA Object #
    if options.verbose:sys.stdout.write('\t-- Get DNA %s ... \n' %options.input_pdb_file)
    x3dna_obj = x3dna.get_x3dna_obj(os.path.abspath(options.input_pdb_file))
    if len(x3dna_obj.get_dinucleotides().keys()) < 1:
        sys.stdout.write('Missing DNA . . . \n')
        exit(0)

    # Get contacts object #
    if  options.verbose:sys.stdout.write("\t\t-- calculate contacts ...\n")
    contacts_obj = contacts.get_contacts_obj(pdb_obj, x3dna_obj)
    if len(contacts_obj.get_contacts())<1:
       sys.stdout.write("Missing Protein-DNA contacts ...\n")
       exit(0)

    # Get triads object #
    if  options.verbose:sys.stdout.write("\t\t-- calculate protein-dna pairs ...\n")
    triads_obj = triads.get_triads_obj(pdb_obj, dssp_obj, x3dna_obj, contacts_obj)

    # Load statistical potential #
    if  options.verbose:sys.stdout.write("\t-- Load potentials ...\n")
    potentials, thresholds , radii, structural_homologs_by_chain = PWM.load_statistical_potentials(pdb_obj,  pdb_dir, pbm_dir, families, radius, options.potential_file, options.split_potential,  options.auto_mode,  options.family_potentials,  options.pbm_potentials,  None,  options.taylor_approach,  options.pmf, options.bins_approach, options.known, None,   dummy_dir, options.verbose)

    # Get protein interface patches
    if  options.verbose:sys.stdout.write("\t-- Get interface patches ...\n")
    patches=get_patches(pdb_obj, contacts_obj, radius, fragment_restrict)

    if options.verbose:
     print " "
     print "\t\t PATCHES "
     print " "
     for pdb_chain_obj in pdb_obj.chains:
      if pdb_chain_obj.chaintype =='P':
        chain = pdb_chain_obj.chain
        seq= pdb_chain_obj.gapped_protein_sequence
        patch_list=patches.get(chain)
        count=0
        print chain,"%3d"%count,seq
        for patch in patch_list:
            dummy_seq=list("." * len(seq))
            for res,idx in patch:
                dummy_seq[idx-1]=seq[idx-1]
            count += 1
            print chain,"%3d"%count,"".join(dummy_seq)
        print " "

    # Get allowed interfaces
    allowed_sequences=None
    if options.allowed_sequences is not None:
      try:
       allowed_file=open(options.allowed_sequences,"r")
       if  options.verbose: sys.stdout.write("\t-- Open file of allowed interface sequences %s\n"%(options.allowed_sequences))
       set_of_sequences=set()
       for line in allowed_file:
           if line.startswith("#"): continue
           data=line.strip().split()
           set_of_sequences.add(data[0])
       if len(set_of_sequences)>0:
          allowed_sequences=set_of_sequences
          if  options.verbose:sys.stdout.write("\t\t-- Total of sequences allowed %d ...\n"%len(set_of_sequences))
      except:
       sys.stdout.write("\t-- Cannot open file: %s\n"%(options.allowed_sequences))

    # Get MSA objects per chain #
    if  options.verbose:sys.stdout.write("\t-- Get MSA object ...\n")
    msa_full,msa_gapped,msa_short= get_msa_obj(pdb_obj, triads_obj, x3dna_obj, potentials,  radius, fragment_restrict, binding_restrict, options.split_potential, patches,  seq_threshold , options.verbose,options.specificity,allowed_sequences)

    # Write PWM #
    #FULL SEQUENCE
    for chain,msa_obj in msa_full.iteritems():
            # Check msa_obj contain sequences, otherwise generate a dummy set
            if options.reset: msa_obj.clean_sequences()
            if len(msa_obj.get_sequences()) == 0 and len(msa_obj.get_pwm()) > 0:
                if  options.verbose:sys.stdout.write("\t\t-- Set up sequences ...\n")
                msa_obj.set_sequences()
            if len(msa_obj.get_sequences())<=0:
                if  options.verbose:sys.stdout.write("\t--No available sequences were accepted for chain  %s ...\n"%chain)
                continue
            #set up name of files
            motif_name= options.output_file.split("/")[-1]+"_"+chain
            msa_obj.set_motif(motif_name)
            pwm_file =options.output_file+"_%s.pwm" %chain
            pwm_meme =options.output_file+"_%s.meme" %chain
            msa_file =options.output_file+"_%s.msa" %chain
            logo_file=options.output_file+"_%s.logo" %chain
            #write results
            if  options.verbose:sys.stdout.write("\t--Write full sequence  PWM of chain %s ...\n"%chain)
            if not os.path.exists(pwm_file):
             if  options.verbose:sys.stdout.write("\t\t--PWM...\n")
             msa_obj.write(pwm_file, option="pwm",overwrite=True)
            if not os.path.exists(pwm_meme):
             if  options.verbose:sys.stdout.write("\t\t--PWM in MEME format...\n")
             msa_obj.write(pwm_meme, option="meme",overwrite=True)
            if not os.path.exists(msa_file):
             if  options.verbose:sys.stdout.write("\t\t--MSA...\n")
             msa_obj.write(msa_file, option="msa",overwrite=True)
            if not os.path.exists(logo_file+".png"):
             if  options.verbose:sys.stdout.write("\t\t--Logos...\n")
             PWM.write_protein_logo(msa_obj,logo_file, dummy_dir)
    #SHORT SEQUENCE
    for chain,msa_obj in msa_short.iteritems():
            # Check msa_obj contain sequences, otherwise generate a dummy set
            if options.reset: msa_obj.clean_sequences()
            if len(msa_obj.get_sequences()) == 0 and len(msa_obj.get_pwm()) > 0:
                if  options.verbose:sys.stdout.write("\t\t-- Set up sequences ...\n")
                msa_obj.set_sequences()
            if len(msa_obj.get_sequences())<=0:
                if  options.verbose:sys.stdout.write("\t--No available sequences were accepted for chain  %s ...\n"%chain)
                continue
            #set up name of files
            motif_name= options.output_file.split("/")[-1]+"_"+chain
            msa_obj.set_motif(motif_name)
            pwm_file =options.output_file+"_short_%s.pwm" %chain
            pwm_meme =options.output_file+"_short_%s.meme" %chain
            msa_file =options.output_file+"_short_%s.msa" %chain
            logo_file=options.output_file+"_short_%s.logo" %chain
            #write results
            if  options.verbose:sys.stdout.write("\t--Write REDUCED PWM of chain %s ...\n"%chain)
            if not os.path.exists(pwm_file):
             if  options.verbose:sys.stdout.write("\t\t--PWM...\n")
             msa_obj.write(pwm_file, option="pwm",overwrite=True)
            if not os.path.exists(pwm_meme):
             if  options.verbose:sys.stdout.write("\t\t--PWM in MEME format...\n")
             msa_obj.write(pwm_meme, option="meme",overwrite=True)
            if not os.path.exists(msa_file):
             if  options.verbose:sys.stdout.write("\t\t--MSA...\n")
             msa_obj.write(msa_file, option="msa",overwrite=True)
            if not os.path.exists(logo_file+".png"):
             if  options.verbose:sys.stdout.write("\t\t--Logos...\n")
             PWM.write_protein_logo(msa_obj,logo_file, dummy_dir)
    #GAPPED SEQUENCE
    for chain,msa_obj in msa_gapped.iteritems():
            # Check msa_obj contain sequences, otherwise generate a dummy set
            if options.reset: msa_obj.clean_sequences()
            if len(msa_obj.get_sequences()) == 0 and len(msa_obj.get_pwm()) > 0:
                if  options.verbose:sys.stdout.write("\t\t-- Set up sequences ...\n")
                msa_obj.set_sequences()
            if len(msa_obj.get_sequences())<=0:
                if  options.verbose:sys.stdout.write("\t--No available sequences were accepted for chain  %s ...\n"%chain)
                continue
            #set up name of files
            motif_name= options.output_file.split("/")[-1]+"_"+chain
            msa_obj.set_motif(motif_name)
            pwm_file =options.output_file+"_gapped_%s.pwm" %chain
            pwm_meme =options.output_file+"_gapped_%s.meme" %chain
            msa_file =options.output_file+"_gapped_%s.msa" %chain
            logo_file=options.output_file+"_gapped_%s.logo" %chain
            #write results
            if  options.verbose:sys.stdout.write("\t--Write GAPPED PWM of chain %s ...\n"%chain)
            if not os.path.exists(pwm_file):
             if  options.verbose:sys.stdout.write("\t\t--PWM...\n")
             msa_obj.write(pwm_file, option="pwm",overwrite=True)
            if not os.path.exists(pwm_meme):
             if  options.verbose:sys.stdout.write("\t\t--PWM in MEME format...\n")
             msa_obj.write(pwm_meme, option="meme",overwrite=True)
            if not os.path.exists(msa_file):
             if  options.verbose:sys.stdout.write("\t\t--MSA...\n")
             msa_obj.write(msa_file, option="msa",overwrite=True)
            if not os.path.exists(logo_file+".png"):
             if  options.verbose:sys.stdout.write("\t\t--Logos...\n")
             PWM.write_protein_logo(msa_obj,logo_file, dummy_dir)



    if options.verbose:sys.stdout.write("\nDone\n")
    exit(0)

