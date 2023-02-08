# Example for: conjugate_gradients(), molecular_dynamics(), model.switch_trace()

# This will optimize stereochemistry of a given model, including
# non-bonded contacts.

import sys
import os
import warnings
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter


from modpyHelper import basic_parser
from modpyHelper import identify_pir

from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients, molecular_dynamics, actions


# Get scripts path (i.e. ".") #
scripts_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),"../../","scripts")

# Append scripts path to python path #
sys.path.append(scripts_path)

# Imports jbonet's module #
from SBI.structure import PDB

def set_options(*args, **kargs):
    '''
    Set the specific options for this script

    @return: parsed ArgumentParser
    '''
    parser = ArgumentParser()

    parser.add_argument("--pdb", dest="pdb",     action="store",
                        metavar="PDB", default=None,
                        help="PDB file with input structure (mandatory)")
    parser.add_argument("--output",  dest="output", action="store",
                        default=None, help="PDB file with optimized structure")
    parser.add_argument("--log",  dest="info_file", action="store",
                        default="optimized.log", help="Information file with the name of the optimized structures")


    return parser.parse_args()


def optimize(code,chains,output,verbose=False):

  print('STDOUT & STDERR will be redirected to log files.')
  sys.stdout = open(output.rstrip(".pdb") + '.log', 'w')
  sys.stderr = open(output.rstrip(".pdb") + '.err', 'w')

  if verbose: log.verbose()  # Commands MODELLER to display all log output


  env = environ()
  env.io.atom_files_directory = ['../atom_files','.']
  env.edat.dynamic_sphere = True

  # Read in HETATM records from template PDBs
  env.io.hetatm = True
  env.libs.topology.read(file='$(LIB)/top_heav.lib')
  env.libs.parameters.read(file='$(LIB)/par.lib')




  #code = 'model'
  mdl = complete_pdb(env, code+".pdb")
  mdl.write(file=code+'.ini')

  # Select all atoms:
  atmsel = selection(mdl)

  # Generate the restraints:
  mdl.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)
  mdl.restraints.write(file=code+'.rsr')

  mpdf = atmsel.energy()

  # Create optimizer objects and set defaults for all further optimizations
  cg = conjugate_gradients(output='REPORT')
  md = molecular_dynamics(output='REPORT')

  # Open a file to get basic stats on each optimization
  trcfil = open(code+'.track.energies', 'w')

  # Run CG on the all-atom selection; write stats every 5 steps
  cg.optimize(atmsel, max_iterations=200000, actions=actions.trace(5, trcfil))
  # Run MD; 
  md.optimize(atmsel, temperature=5, max_iterations=50, actions=[actions.trace(10, trcfil)])
  cg.optimize(atmsel, max_iterations=200000, actions=actions.trace(5, trcfil))
  md.optimize(atmsel, temperature=10, max_iterations=100, actions=[actions.trace(10, trcfil)])
  cg.optimize(atmsel, max_iterations=200000, actions=actions.trace(5, trcfil))
  md.optimize(atmsel, temperature=50, max_iterations=150, actions=[actions.trace(10, trcfil)])
  cg.optimize(atmsel, max_iterations=200000, actions=actions.trace(5, trcfil))
  md.optimize(atmsel, temperature=100, max_iterations=200, actions=[actions.trace(10, trcfil)])
  cg.optimize(atmsel, max_iterations=200000, actions=actions.trace(5, trcfil))
  md.optimize(atmsel, temperature=150, max_iterations=150, actions=[actions.trace(10, trcfil)])
  cg.optimize(atmsel, max_iterations=200000, actions=actions.trace(5, trcfil))
  md.optimize(atmsel, temperature=200, max_iterations=100, actions=[actions.trace(10, trcfil)])
  cg.optimize(atmsel, max_iterations=200000, actions=actions.trace(5, trcfil))
  md.optimize(atmsel, temperature=250, max_iterations=50, actions=[actions.trace(10, trcfil)])
  cg.optimize(atmsel, max_iterations=200000, actions=actions.trace(5, trcfil))
  md.optimize(atmsel, temperature=300, max_iterations=50, actions=[actions.trace(10, trcfil)])
  # Finish off with some more CG, and write stats every 5 steps
  #last optimization
  cg.optimize(atmsel, max_iterations=2000000, actions=[actions.trace(5, trcfil)])

  mpdf = atmsel.energy()

  # renumber
  numbering=[1 for c in chains]
  mdl.rename_segments(segment_ids=chains,renumber_residues=numbering)
  # write last model
  mdl.write(file=output)


if __name__ == '__main__':

    options = set_options()

    if options.pdb is None or not os.path.exists(options.pdb):
        print("Please introduce a correct PDB file")
        info=open(options.info_file,"a")
        info.write("%s\tFAIL\n"%(os.path.basename(options.pdb)))
        info.close()
        exit()

    if options.output is None:
       output = options.pdb.rstrip(".pdb")+".opt.pdb"
    else:
       output = options.output
    pdb = PDB(options.pdb)
    code= options.pdb.rstrip(".pdb")
    #chains=p.chain_identifiers
    chains=[]
    for chain in pdb.chains:
        chains.append(chain.chain)

    info=open(options.info_file,"a")
    try:
      optimize(code,chains,output)
      if os.path.exists(output):
        info.write("%s\tDONE\n"%(os.path.basename(options.pdb)))
      else:
        info.write("%s\tFAIL\n"%(os.path.basename(options.pdb)))
    except:
      info.write("%s\tFAIL\n"%(os.path.basename(options.pdb)))
    info.close()

    print("Done")


