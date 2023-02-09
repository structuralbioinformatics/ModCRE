import sys
import os
import shutil
import subprocess
import numpy as np
import pandas as pd
import argparse
import shutil
import optparse
import IMP.pmi.dof
import IMP.atom
import IMP.pmi.macros
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.tools
import IMP.pmi.topology
import IMP.pmi.output
import IMP.pmi.mmcif
import ihm
import ihm.location
import ihm.model
import ihm.dumper
import IMP



#-------------#
# Functions   #
#-------------#



def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python model_IMP.py -i input_folder -m model_name  -o output [--info log_file]")

    parser.add_option("-i", action="store", type="string", dest="input_folder", default=None, help="Input folder, contains the inputs orftopology and restraints to model macro-complexes", metavar="{directory}")
    parser.add_option("-m", action="store", type="string", dest="model_name", default=None, help="Name of the model", metavar="{filename}")
    parser.add_option("--num_frames", action="store", type="int", dest="num_frames", default=100, help="Number of frames", metavar="{integer}")
    parser.add_option("--num_mc_steps", action="store", type="int", dest="num_mc_steps", default=50, help="Number of Monte Carlo steps", metavar="{integer}")
    parser.add_option("-o", action="store", type="string", dest="output", default=None, help="Output folder rootname  ", metavar="{directory}")
    parser.add_option("--info",default=None,action="store", type="string", dest="info",help="Information LOG file of MODELS that have failed and have been completed (default modelling_IMP_execution.log)")
 

    (options, args) = parser.parse_args()

    if options.input_folder is None or options.output is None or options.model_name is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options   = parse_options()

    # ---------------------------
    # 1. Define Input Data and Output Directories
    # ---------------------------

    input_folder     = options.input_folder
    output_directory = options.output
    model            = options.model_name
    info_file        = options.info
    num_frames       = int(options.num_frames)
    num_mc_steps     = int(options.num_mc_steps)

    if not input_folder.startswith("/"):  input_folder = os.path.abspath(input_folder)
    if not output_directory.startswith("/"):   output_directory = os.path.abspath(output_directory )
    if info_file is None: info_file = os.path.join(input_folder,"modelling_IMP_execution.log")

    #-------------------------#
    # General parameters      #
    #-------------------------#


    # --------------------------
    # 2. Scoring Parameters
    # --------------------------

    resolution=1.0  #Medium resolution to identify residues and speed up modeling

    # --------------
    # ----- Sterochemistry and Physical Restraints
    ev_weight = 1.0  # Weight of excluded volume restraint
    connectivity_scale = 1.0  # weight of Connectivity restraint
    ev_resolution = 10 #To speed up this expensive restraint, we evaluate it at resolution 10

    # --------------------
    # 3. Sampling Parameters
    # --------------------
    num_steps_beads = 100 # Quickly move all flexible beads into place
    if num_frames < 1:
       num_frames = 100  #  Number of frames in MC run
    num_best_scoring_models = 1
    if num_mc_steps < 1:
       num_mc_steps = 50  # Number of MC steps per frame
    mc_temperature = 1.0  # Temperature for MC
    # --- Simulated Annealing (sa)
    #  - Alternates between two MC temperatures
    sim_annealing = True  # If true, run simulated annealing
    sa_min_temp_steps = 100  # Steps at min temp
    sa_max_temp_steps = 20  # Steps at max temp
    sa_temps = (1.0, 5.0)  # Sim annealing temperatures

    # Replica Exchange (rex)
    rex_temps = (1.0, 5.0)  # Temperature bounds for replica exchange

    # --------------------
    # 4. Model INPUT parameters
    # --------------------
    restraints_file = os.path.join(input_folder,model+".restraints.csv")
    fasta_file      = os.path.join(input_folder,model+".topology.fasta")
    topology_file   = os.path.join(input_folder,model+".topology.txt")

    print("RESTRAINTS %s"%restraints_file)
    print("FASTA  %s"%fasta_file)
    print("TOPOLOGY  %s"%topology_file)

    # --------------------
    # 5. Modelling
    # --------------------

    # Try modelling or write exception
    try:

      # Initialize model
      m = IMP.Model()
      # Read in the topology file --> to coarse model 
      # Specify the directory with the PDB files and fasta files
      topology = IMP.pmi.topology.TopologyReader(topology_file,
                                               pdb_dir=input_folder,
                                               fasta_dir=input_folder)

      # Use the BuildSystem macro to build states from the topology file
      bs = IMP.pmi.macros.BuildSystem(m)
      bs.add_state(topology)
      system_molecules = bs.get_molecules()

      #UNCOMMENT next lines to Record the modeling protocol to an mmCIF file
      #po = IMP.pmi.mmcif.ProtocolOutput()
      #bs.system.add_protocol_output(po)
      #po.system.title = "IMP model of "+os.path.basename(input_folder)+" label "+model

      # Build the system representation and degrees of freedom
      root_hierarchy, dof = bs.execute_macro(max_rb_trans=4.0,
                                           max_rb_rot=0.3,
                                           max_bead_trans=4.0,
                                           max_srb_trans=4.0,
                                           max_srb_rot=0.3)

      # Randomize the initial configuration before sampling
      IMP.pmi.tools.shuffle_configuration(root_hierarchy,
                                          # excluded_rigid_bodies=dof.get_rigid_bodies(),
                                          max_translation=50,
                                          verbose=False,
                                          cutoff=5.0,
                                          niterations=100)

      # Make a list of IMP object restraints 
      imp_objects = list()  # reporter objects
      # Introduce distance restraints (dr)
      restraints_data      = pd.read_csv(restraints_file)
      number_of_restraints = restraints_data.shape[0]
      restraint_feature    = [c for c in restraints_data.columns]
      for i in range(number_of_restraints):
         kd       = float(restraints_data.iloc[i]["kd"])
         mol1     = restraints_data.iloc[i]["prot1"]
         mol2     = restraints_data.iloc[i]["prot2"]
         res1_idx = (restraints_data.iloc[i]["prot1_res"])
         res2_idx = (restraints_data.iloc[i]["prot2_res"])
         distance = float(restraints_data.iloc[i]["distance"])
         weight   = float(restraints_data.iloc[i]["weight"])
         sd       = float(restraints_data.iloc[i]["sd"])
         label    = "restraint_"+str(i+1)+"_"+mol1+str(res1_idx)+"_"+mol2+str(res2_idx)
         tuple_selection1 = tuple([res1_idx, res1_idx, mol1, 0])
         tuple_selection2 = tuple([res2_idx, res2_idx, mol2, 0])
         print("Add restraint %d:  %s %s with %s %s at distance %f Weight %f Deviation %f KD %f"%(i,mol1,res1_idx,mol2, res2_idx,distance,weight,sd,kd))
         try:
           dr     = IMP.pmi.restraints.basic.DistanceRestraint(root_hier=root_hierarchy,
                                                             tuple_selection1=tuple_selection1,
                                                             tuple_selection2=tuple_selection2,
                                                             distancemin=max(0,distance-sd),
                                                             distancemax=(distance+sd),
                                                             resolution=resolution,
                                                             kappa=kd,
                                                             label=label,
                                                             weight=weight ) 

           dr.add_to_model()
           dr.evaluate()
           imp_objects.append(dr)
         except:
           print("Missing atoms in model, SKIP RESTRAINT %d (%s %s with %s %s)"%(i,mol1,res1_idx,mol2,res2_idx)) 

      # Introduce the connectivity of molecules as restraints (cr)
      all_molecules = IMP.pmi.tools.get_molecules(root_hierarchy)
      for mol in all_molecules:
         mol_name = mol.get_name()
         IMP.pmi.tools.display_bonds(mol)
         cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol, scale=connectivity_scale)
         cr.add_to_model()
         cr.set_label(mol_name)
         imp_objects.append(cr)

      # Excluded volume restraints (evr)
      evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=all_molecules,resolution=ev_resolution)
      evr.set_weight(ev_weight)
      evr.add_to_model()
      imp_objects.append(evr)

      # Quickly move all flexible beads into place
      dof.optimize_flexible_beads(nsteps=num_steps_beads)

      # --------------------------
      # Monte-Carlo Sampling
      # --------------------------

      # Define all components to be sampled and the sampling protocol
      mc1 = IMP.pmi.macros.ReplicaExchange0(m,
                                          root_hier=root_hierarchy,
                                          monte_carlo_sample_objects=dof.get_movers(),
                                          output_objects=imp_objects,
                                          rmf_output_objects=imp_objects,
                                          monte_carlo_temperature=mc_temperature,
                                          simulated_annealing=sim_annealing,
                                          simulated_annealing_minimum_temperature=min(sa_temps),
                                          simulated_annealing_maximum_temperature=max(sa_temps),
                                          simulated_annealing_minimum_temperature_nframes=sa_min_temp_steps,
                                          simulated_annealing_maximum_temperature_nframes=sa_max_temp_steps,
                                          replica_exchange_minimum_temperature=min(rex_temps),
                                          replica_exchange_maximum_temperature=max(rex_temps),
                                          number_of_best_scoring_models=num_best_scoring_models,
                                          monte_carlo_steps=num_mc_steps,
                                          number_of_frames=num_frames,
                                          global_output_directory=output_directory + "_" + str(model))

      # Start Sampling
      mc1.execute_macro()

      #UNCOMMENT next lines to write system protocol in mmcif
      #po.finalize()
      #s = po.system
      #if not os.path.exists( os.path.join(output_directory + "_" + str(model),"mmcif" ):
      #   os.makedirs( os.path.join(output_directory + "_" + str(model),"mmcif" )
      #mmcif_file = os.path.join(output_directory + "_" + str(model),"mmcif","models.cif")
      #with open(mmcif_file, 'w') as fh:
      #     ihm.dumper.write(fh, [s])

    except Exception as e:
      # Write exception if it fails
      print("Failed IMP modeling: %s "%e)
      info=open(info_file,"a")
      info.write("%s\tFAIL\n"%(model))
      info.close()
      exit(0)


    # --------------------
    # 6. Finalize 
    # --------------------

    
    info=open(info_file,"a")
    info.write("%s\tDONE\n"%(model))
    info.close()
    print("Done")








