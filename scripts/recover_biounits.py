import os, sys
import optparse
import ConfigParser
import shutil

# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.dirname(__file__))

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Imports my functions #
import functions

# Imports jbonet's module #
from SBI.structure import PDB
from SBI.structure.chain import Chain



#-------------#
# Options     #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python remove_biounits.py --pdb=pdb_dir")
    parser.add_option("--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (output directory from pdb.py)", metavar="{directory}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")
    (options, args) = parser.parse_args()

    if options.pdb_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()
    families_file = os.path.join(options.pdb_dir, "families.txt")
    remove_biounit=set()
    for linen in functions.parse_file(families_file):
        if linen.startswith("#"): continue
        line = linen.lstrip("#").strip().split(";")
        if line[0] == "pdb": continue
        if len(line)>2:
         if line[2].split()[0] == "Failed":
           pdb_name=line[0].split("_")[0]
           if options.verbose:sys.stdout.write("\t-- recover from biounits %s\n"%(pdb_name))
           remove_biounit.add(pdb_name)

    for pdb_name in remove_biounit:
        original_pdb_file = os.path.abspath(os.path.join(options.pdb_dir, "biounits", pdb_name + ".pdb1.skip"))
        if not os.path.exists(original_pdb_file): original_pdb_file = os.path.abspath(os.path.join(options.pdb_dir, "biounits", pdb_name.upper() + ".pdb1.skip"))
        if not os.path.exists(original_pdb_file): original_pdb_file = os.path.abspath(os.path.join(options.pdb_dir, "biounits", pdb_name.lower() + ".pdb1.skip"))
        if  os.path.exists(original_pdb_file): 
            modified_pdb_file = os.path.abspath(os.path.join(options.pdb_dir, "biounits", pdb_name.upper() + ".pdb1"))
            shutil.move(original_pdb_file,modified_pdb_file)

        original_pdb_file = os.path.abspath(os.path.join(options.pdb_dir, "clean", pdb_name + ".pdb"))
        if os.path.exists(original_pdb_file):
           os.remove(original_pdb_file)



        
        



