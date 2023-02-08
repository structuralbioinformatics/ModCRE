import os, sys, re
from collections import Counter
import ConfigParser
import itertools
import numpy as np
import optparse
import shutil
import random
import threader
import cPickle
import pandas as pd

import matplotlib as mpl 
# agg backend is used to create plot as a .png file
mpl.use('Agg')
import matplotlib.pyplot as plt 

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
import contacts, dssp, interface,  spotentials, tmalign, triads, x3dna, scorer, fimo, model_protein, model_dna, threading_to_triads
import pwm_pbm as PWM


#-------------#
# Functions   #
#-------------#

def calculate_single_profile_by_thread(fimo_thresholds,pdb_file,output_file,pbm_dir,pdb_dir,families,potential_file, radius,fragment_restrict, binding_restrict, split_potential,auto_mode,family_potentials,pbm_potentials,score_threshold,taylor_approach,pmf,bins,known,meme,reset, dummy_dir,verbose,save):

    
#Get the PWM
    if not functions.fileExist(output_file+".meme"):
           if verbose:sys.stdout.write("Generate PWM %s ...\n"%(output_file))
           PWM.get_single_pwm(pdb_file,None,False,output_file,pbm_dir,pdb_dir,families,potential_file, radius,fragment_restrict, binding_restrict, split_potential,auto_mode,family_potentials,pbm_potentials,score_threshold,taylor_approach,pmf,bins,known,meme,reset, dummy_dir,verbose)

    #Scan the DNA sequence with the PWM to get FIMO scores
    if verbose: sys.stdout.write("Get FIMO scores\n")
    pwm_file=output_file+".meme"
    max_stored_matches = None
    profile=ProfileProtein(dna_sequence,fimo_thresholds)
    for threshold in fimo_thresholds:
        if verbose: sys.stdout.write("\t-- Threshold %s\n"%threshold)
        fimo_obj = fimo.get_fimo_obj(pwm_file, dna_file, float(threshold), max_stored_matches, dummy_dir)
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
    if verbose: sys.stdout.write("\t\t-- Get triads ...\n")
    triads_obj = triads.get_triads_obj(pdb_obj, dssp_obj, x3dna_obj, contacts_obj)


# Load statistical potential
    if radius <= 0: radius=float(config.get("Parameters", "max_contact_distance"))
    if verbose: sys.stdout.write("\t-- Load potentials ...\n")
    potentials, thresholds, radii, structural_homologs_by_chain = PWM.load_statistical_potentials(pdb_obj, pdb_dir, pbm_dir, families, radius, potential_file, split_potential, auto_mode, family_potentials, pbm_potentials, score_threshold, taylor_approach, pmf , bins, known, None, dummy_dir,verbose)

# Get the dictionary of msa_objs
    msa_objs={}
    split_potentials=["3d","3dc","local","pair","s3dc","s3dc_di","s3dc_dd"]
    for sp in split_potentials:
          if verbose: sys.stdout.write("\t-- Get MSA for potential %s ...\n"%sp)
          msa_obj = PWM.get_msa_obj(triads_obj, x3dna_obj, potentials, radii, fragment_restrict, binding_restrict, sp, thresholds)
          msa_objs.setdefault(sp,msa_obj)


# Split DNA sequence according to interface
    # All TF-DNA models generated with model_dna are restricted to the DNA interface
    # Get the interface length
    interface_length=interface_obj.get_interface_length()
    if verbose: 
        sys.stdout.write("Interface match of PDB\n")
        sys.stdout.write("\t-- Length %d\n"%interface_length)
        sys.stdout.write("\t-- Starting basepairs in structure %d\n"%interface_obj.get_interface_start())
        sys.stdout.write("\t-- Last basepairs in structure %d\n"%interface_obj.get_interface_end())
    if len(dna_sequence)<interface_length:
        sys.stdout.write("Stop. Length of DNA is shorter than interface\n")
        info.write("%s\tFAILED\n"%output_name)
        exit(0)

    dummy_pdb_dir = os.path.join(dummy_dir,os.path.basename(output_dir))
    if not os.path.exists(dummy_pdb_dir): os.makedirs(dummy_pdb_dir)

    # Get the energy-scores by windows along the DNA sequence
    for i in range(len(dna_sequence)-interface_length+1):

        # Set the DNA fragment
        dna_fragment=dna_sequence[i:i+interface_length]
        if verbose: sys.stdout.write("Score %d st fragment %s\n"%(i+1,dna_fragment))
        pdb_dna_obj = model_dna.get_dna_model_pdb_obj(pdb_file, dna_fragment, x3dna_obj, interface_obj, interface_range=None, dummy_dir=dummy_dir)

        # Create PDB with DNA fragment
        if save: pdb_dna_obj.write(os.path.join(output_dir, output_name + "." + str(i+1) + "." + dna_fragment + ".pdb"),force=True)

        # Thread DNA fragment to get new triads 
        if verbose: sys.stdout.write("\t-- Thread DNA fragment on triads ...\n")
        fragment_triads_obj = threading_to_triads.thread_kmer_triads(dna_fragment,triads_obj)

        # Calculate the energy scores and fill the profile#
        if verbose: sys.stdout.write("\t-- Calculate scores potential: %s...\n"%energy_profile)
        scr=scorer.score(potential=energy_profile)
        scr.calculate_energies_and_binding(fragment_triads_obj, x3dna_obj, msa_objs, potentials,  thresholds, radii, radius, fragment_restrict, binding_restrict, dummy_dir)
        profile.set_energy_score_by_position(i,scr)
        if save:
           output = output_file+".scr_"+str(i+1)+".txt"
           protein_name=os.path.basename(pdb_file)
           scr.write(output,protein_name=protein_name,normal=True,overwrite=True)

    return profile

#-------------#
# Classes     #
#-------------#

class ProfileProtein(object):
    """
    This class define a single protein PROFILE object
    """
    def __init__(self,dna,thresholds=["0.05", "0.001"],binding=None,fimo_log_score=None,fimo_score=None, energy_scores=None, potential=None):
        self._dna             = dna
        self._thresholds      = thresholds
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
        self._thresholds=thresholds

    def get_thresholds(self):
        return self._thresholds

    def set_fimo_scores_by_threshold(self, fimo_obj, threshold):
        binding_score  = np.zeros(len(self._dna))
        fimo_score     = np.zeros(len(self._dna))
        fimo_log_score = np.zeros(len(self._dna))
        for fimo_hit in fimo_obj.get_hits():
            start = fimo_hit.get_start()
            end   = fimo_hit.get_end()
            score = fimo_hit.get_score()
            logscr= -np.log(fimo_hit.get_p_value())
            binding_score[start-1:end]   =1
            fimo_score[start-1:end]     += score
            fimo_log_score[start-1:end] += logscr
        self._binding.setdefault(threshold,binding_score)
        self._fimo_log_score.setdefault(threshold,fimo_log_score)
        self._fimo_score.setdefault(threshold,fimo_score)

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
        energy=np.zeros(len(self._dna))
        scr=self.get_energy_score_by_position(0)
        binding_length=len(scr.get_binding_site())
        for i in range(len(self._dna)-binding_length+1):   
            scr=self.get_energy_score_by_position(i)
            if scr is None: continue
            energy[i:i+binding_length] += scr.get_score(normal,potential)
        if normal: energy = energy/binding_length
        return energy

    def get_energy_weighted_profile(self,threshold,normal=False,potential="s3dc_dd"):
        if self.get_fimo_binding_by_threshold(threshold) is None: return None
        binding_score=self.get_fimo_binding_by_threshold(threshold)
        energy=np.zeros(len(self._dna))
        scr=self.get_energy_score_by_position(0)
        binding_length=len(scr.get_binding_site())
        for i in range(len(self._dna)-binding_length+1):   
            scr=self.get_energy_score_by_position(i)
            if scr is None: continue
            #if one nucleotide of the interval is not part of a binding the score is skipped
            #otherwise, if all nucleotide positions are in binding add the score
            energy[i:i+binding_length] += scr.get_score(normal,potential) * binding_score[i:i+binding_length].prod()
        if normal: energy = energy/binding_length
        return energy

    def get_energy_best_profile(self,potential="s3dc_dd"):
        normal=True
        energy=np.zeros(len(self._dna))
        scr=self.get_energy_score_by_position(0)
        binding_length=len(scr.get_binding_site())
        for i in range(len(self._dna)-binding_length+1):   
            scr=self.get_energy_score_by_position(i)
            if scr is None: continue
            energy[i:i+binding_length] = max(scr.get_score(normal,potential),energy[i:i+binding_length].max())
        return energy

    def get_energy_per_nucleotide_profile(self,normal=False,potential="s3dc_dd"):
        energy=np.zeros(len(self._dna))
        scr=self.get_energy_score_by_position(0)
        binding_length=len(scr.get_binding_site())
        for i in range(len(self._dna)-binding_length+1):   
            scr=self.get_energy_score_by_position(i)
            if scr is None: continue
            score_per_nucleotide=scr.get_score_per_nucleotide(normal,potential)
            start = min([int(nuc) for nuc in sorted(score_per_nucleotide.keys())])
            for nuc in sorted(score_per_nucleotide.keys()):
                position = int(nuc) - start
                energy[i+position] += score_per_nucleotide.get(nuc)
        if normal: energy = energy/len(score_per_nucleotide)
        return energy

    def to_table(self):
        split_potentials=["3d","3dc","local","pair","s3dc","s3dc_di","s3dc_dd"]
        thresholds=self.get_thresholds()
        tbl_dict={}
        tbl_dict.setdefault("Nucleotide",[int(x)+1 for x in range(len(self._dna))])
        for thr in thresholds:
            tbl_dict.setdefault("binding_"+thr,self.get_fimo_binding_by_threshold(thr))
            tbl_dict.setdefault("fimo_score_"+thr,self.get_fimo_score_by_threshold(thr))
            tbl_dict.setdefault("fimo_log_score_"+thr,self.get_fimo_log_score_by_threshold(thr))
        for potential in split_potentials:
            tbl_dict.setdefault(potential,self.get_energy_profile(normal=False,potential=potential))
            tbl_dict.setdefault("normal_"+potential,self.get_energy_profile(normal=True,potential=potential))
            tbl_dict.setdefault(potential+"_x_nucleotide",self.get_energy_per_nucleotide_profile(normal=False,potential=potential))
            tbl_dict.setdefault("normal_"+potential+"_x_nucleotide",self.get_energy_per_nucleotide_profile(normal=True,potential=potential))
            tbl_dict.setdefault("best_"+potential,self.get_energy_best_profile(potential=potential))
        for thr in thresholds:
          for potential in split_potentials:
            tbl_dict.setdefault(potential+"_weighted_"+thr,self.get_energy_weighted_profile(thr,normal=False,potential=potential))
            tbl_dict.setdefault("normal_"+potential+"_weighted_"+thr,self.get_energy_weighted_profile(thr,normal=True,potential=potential))
        table=pd.DataFrame(tbl_dict)
        return table

    def write_table(self,output):
        table=self.to_table()
        table.to_csv(output)

    def plot(self,output):
        table=self.to_table()
        columns=table.columns.to_list()
        x = table["Nucleotide"].values
        for column in columns:
          if column == "Nucleotide": continue
          out   = output+"_"+column+".png"
          title = "Graph of profile %s"%column
          y     = table[column].values
          if len(y[abs(y)>1.0e-10])<1: continue
          plt.plot(x,y,'-o', ms=5, lw=2, alpha=0.7, mfc='cyan')
          plt.xlabel('nucleotide position')
          plt.ylabel('score profile')
          plt.savefig(out,format="png")
          plt.close()


#-------------#
# Options     #
#-------------#

def parse_options():
    '''
    This function parses the command line arguments and returns an optparse object.
    '''

    parser = optparse.OptionParser("Usage: single_profile.py [--dummy=DUMMY_DIR] -i INPUT_FILE -d DNA_FASTA [-l LABEL -o OUTPUT_NAME --output_dir OUTPUT_DIR ] --pbm=PBM_dir --pdb=PDB_DIR [-v --save --plot --meme --reset] [-a -f -p -s SPLIT_POTENTIAL -e ENERGY_PROFILE -t THRESHOLD -k -b --taylor --file POTENTIAL --radius RADIUS --fragment FRAGMENT]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="DUMMY_DIR")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input PDB file. Mandatory", metavar="INPUT_FILE")
    parser.add_option("-l", action="store", type="string", dest="label", help="Label to organize the output files as 'label.energies.txt' ", metavar="LABEL")
    parser.add_option("--meme", default=False, action="store_true", dest="meme", help="Use 'uniprobe2meme' to calculate the PWM matrix for 'FIMO' (default = False)")
    parser.add_option("-r","--reset", default=False, action="store_true", dest="reset", help="Clean the sequences of the original MSA and reset them by a random selection in accordance with the PWM (default = False)")
    parser.add_option("-o", "--output", default=None, action="store", type="string", dest="output_name", help="Output name for tables and compressed PICKLE format (default = as PDB input)", metavar="OUTPUT_NAME")
    parser.add_option("--output_dir", default=None, action="store", type="string", dest="output_dir", help="Output directory (default = as PDB input)", metavar="OUTPUT_DIR")
    parser.add_option("--pbm", action="store", type="string", default=None, dest="pbm_dir", help="PBM directory (i.e. output dir from pbm.py). This is ,mandatory unless using option --file on potentials", metavar="PBM_DIR")
    parser.add_option("--pdb", action="store", type="string", default=None, dest="pdb_dir", help="PDB directory (i.e. output dir from pdb.py). This is ,mandatory unless using option --template", metavar="PDB_DIR")
    parser.add_option("-d","--dna",default=None, action="store", type="string", dest="dna_file", help="File of a DNA sequence in FASTA format to profile", metavar="FASTA")
    parser.add_option("--fragment", default=None, action="store", type="string", dest="fragment", help="Fragment of protein to apply the potential. Format is 'a-b;c-d': two regions between residues a-b and c-d. (Default is None it applies to all amino-acids)")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")
    parser.add_option("--save", default=False, action="store_true", dest="save", help="Save PDB models and scores while scanning the DNA (default = False)")
    parser.add_option("--plot", default=False, action="store_true", dest="plot", help="Plot profiles (default = False)")
    parser.add_option("--info", default=None, action="store",  dest="info_file", help="File to store information of SUCCESS/FAILURE of the run (default = 'Information_of_profiles.txt')")


    group = optparse.OptionGroup(parser, "Statistical potentials", "Select your statistical potentials of choice. By default it uses S3DCdd general potential derived from PDB (the simplest one). In \"--auto\" mode, the program uses S3DCdd family potentials derived from both PDB and/or PBM data and/or approached by Taylor as selected in Potentials configuration file. In case family potentials cannot be applied, the program uses general potentials derived from both PDB and PBM data and approached by Taylor. \"-a\" option overrides options \"-f\", \"-p\" and \"-t\".")
    group.add_option("-a", "--auto", default=False, action="store_true", dest="auto_mode", help="Automate the selection of statistical potentials (default = False)")
    group.add_option("-f", "--family", default=False, action="store_true", dest="family_potentials", help="Use family potentials (default = False)")
    group.add_option("-p", default=False, action="store_true", dest="pbm_potentials", help="Use potentials derived from both PBM + PDB data (default = False)")
    group.add_option("-s", default="s3dc_dd", action="store", type="string", dest="split_potential", help="Split-potential to be used (all, 3d, 3dc, s3dc, s3dc_dd, s3dc_di, pair; default = s3dc_dd)", metavar="{string}")
    group.add_option("-e","--energy", default="all", action="store", type="string", dest="energy_profile", help="Select a specific Split-potential to be used for the profile (all, 3d, 3dc, s3dc, s3dc_dd, s3dc_di, pair; default = all)", metavar="{string}")
    group.add_option("-t", action="store", type="float", dest="score_threshold", help="Threshold on the scaled score to consider positive k-mers (default = 0.95)", metavar="{float}")
    group.add_option("-k","--known", default=False, action="store_true", dest="known", help="The name is of a known PDB file, with format 'code_chain' (default = False)")
    group.add_option("--taylor", default=False, action="store_true", dest="taylor_approach", help="Approach PMF by Taylor (default = False)")
    group.add_option("-b", "--bins", default=False, action="store_true",  dest="bins", help="Computate the potentials: by bins (if selected) or accumulative (default).", metavar="{boolean}")
    group.add_option("--file",default=None, action="store", type="string", dest="potential_file", help="Use potentials from specific file (default = None)", metavar="{string}")
    group.add_option("--radius",default=0, action="store", type="float", dest="radius", help="Maximum contact distance to calculate interactions (default=0 implies the use of 'max_contact_distance' from configuration", metavar="{string}")
    parser.add_option("-m", "--pmf", default=False, action="store_true", dest="pmf", help="Use of raw mean-force potentials with no Z-scoring (default = False)")


    parser.add_option_group(group)

    
    (options, args) = parser.parse_args()

    if options.input_file is None or (options.pbm_dir is None and options.potential_file is None) or options.dna_file is None or options.pdb_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")
    if not functions.fileExist(options.dna_file):
        parser.error("Missing DNA fasta file")
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
    distance_tyep="dinucleotides"
    verbose=options.verbose
    save=options.save
    plot=options.plot
    dummy_dir = options.dummy_dir
    if not dummy_dir.startswith("/"): dummy_dir = os.path.abspath(options.dummy_dir)
    if not os.path.exists(dummy_dir): os.makedirs(dummy_dir)
    pdb_dir = options.pdb_dir
    if pdb_dir is not None:
       if not pdb_dir.startswith("/"): pdb_dir = os.path.abspath(options.pdb_dir)
    pbm_dir = options.pbm_dir
    if pbm_dir is not None:
       if not pbm_dir.startswith("/"): pbm_dir = os.path.abspath(options.pbm_dir)
    pdb_file = options.input_file
    if not pdb_file.startswith("/"): pdb_file = os.path.abspath(pdb_file)
    output_dir=os.path.dirname(pdb_file)
    if options.output_dir is not None: output_dir = options.output_dir
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    output_name=pdb_file.rstrip(".pdb")
    if options.output_name is not None: output_name = options.output_name
    if options.label is not None:
        output_file = output_name +"."+ options.label
    output_file = os.path.join(output_dir,output_name)
    if options.info_file is not None:
       info   = open(os.path.join(output_dir,options.info_file),"a")
    else:
       info   = sys.stdout
    dna_file=options.dna_file
    if dna_file is not None:
      if not dna_file.startswith("/"): dna_file=os.path.abspath(dna_file)
    for header, sequence in functions.parse_fasta_file(dna_file):
        if verbose: sys.stdout.write("Read DNA named %s (lenght %d) \n"%(header,len(sequence)))
        dna_sequence = sequence
    families = {}
    if options.verbose:sys.stdout.write("Check families...\n")
    if not os.path.exists(os.path.join(pdb_dir, "families.txt")):
       sys.stdout.write("Families file %s is not used\n"%(os.path.join(pdb_dir, "families.txt")))
    else:
      for line in functions.parse_file(os.path.join(options.pdb_dir, "families.txt")):
        if line.startswith("#"): continue
        pdb_chain, family = line.split(";")
        families[pdb_chain] = family
    # Potential input parameters
    potential_file    =options.potential_file
    radius            =options.radius
    split_potential   =options.split_potential
    energy_profile    =options.energy_profile
    auto_mode         =options.auto_mode
    family_potentials =options.family_potentials
    pbm_potentials    =options.pbm_potentials
    score_threshold   =options.score_threshold
    taylor_approach   =options.taylor_approach
    pmf               =options.pmf
    bins              =options.bins
    known             =options.known
    meme              =options.meme
    reset             =options.reset
    # Restrictions on TF fragment
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
    # No restrictions on binding are allowed
    binding_restrict=None
    # Define fimo_thresholds
    fimo_thresholds=config.get("Parameters", "fimo_profile_thresholds")

# TRY EXECUTION AND STORE THE RESULT
    #try:
    profile=calculate_single_profile_by_thread(fimo_thresholds,pdb_file,output_file,pbm_dir,pdb_dir,families,potential_file, radius,fragment_restrict, binding_restrict, split_potential,auto_mode,family_potentials,pbm_potentials,score_threshold,taylor_approach,pmf,bins,known,meme,reset, dummy_dir,verbose,save)

# Write outputs
    # Name the files
    output = output_file
    # Write table in CSV format
    if functions.fileExist(output+".csv"): os.remove(output+".csv")
    profile.write_table(output+".csv")
    # Write profile in pickle format
    if functions.fileExist(output+".pickle"): os.remove(output+".pickle")
    out=open(output+".pickle","wb")
    cPickle.dump(profile,out)
    out.close()
    # Plot profiles
    if plot: profile.plot(output)
    # Write information of success
    info.write("%s\tDONE\n"%output_name)

# IF PROFILE FAILS
    #except Exception as e:
      #info.write("%s\tFAILED\n"%output_name)
      #sys.stderr.write("ERROR %s"%e)



      



