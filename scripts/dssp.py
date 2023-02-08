import os, sys, re
import ConfigParser
import copy
import optparse
import subprocess

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

#-------------#
# Classes     #
#-------------#

class DSSP(object):
    """
    This class defines a {DSSP} object.

    """

    def __init__(self, file_name):
        self._file = file_name
        self._file_content = []
        self._residues = {}
        # Initialize #
        self._parse_file()

    def _parse_file(self):
        for line in functions.parse_file(self._file):
            if not line.endswith(".") and "#" not in line and line != "":
                m = re.search("(\d+)\s\S\s[ACDEFGHIKLMNPQRSTVWY]", line)
                # Get DSSP info #
                if m:
                    pdb_chain = line[11:12]
                    residue_num = int(line[5:10])
                    accessible_surface_area = float(line[35:38])
                    secondary_structure = line[16:17]
                    self._residues[(pdb_chain, residue_num)] = (accessible_surface_area, secondary_structure)
            # Add line to file content #
            if line != "":
               self._file_content.append(line)

    def has_residue(self, pdb_chain, residue_num):
        if (pdb_chain, residue_num) in self._residues:
            return True

        return False
            
    def get_accessible_surface_area(self, pdb_chain, residue_num):
        if self.has_residue(pdb_chain, residue_num):
            return copy.copy(self._residues[(pdb_chain, residue_num)][0])

        return None

    def get_secondary_structure(self, pdb_chain, residue_num):
        if self.has_residue(pdb_chain, residue_num):
            return copy.copy(self._residues[(pdb_chain, residue_num)][1])

        return None

    def write(self, file_name):
        for line in self._file_content:
            functions.write(file_name, line)

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python dssp.py -i input_file [--dummy=dummy_dir -o output_file]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in PDB format)", metavar="{filename}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")

    (options, args) = parser.parse_args()

    if options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_dssp_obj(pdb_file, dummy_dir="/tmp"):
    """
    This function executes "find_pair" from X3DNA package and returns a {X3DNA}.

    @input:
    pdb_file {string}
    dummy_dir {string}

    @return:
    x3dna_obj {X3DNA}

    """

    try:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        dssp_path = os.path.join(src_path, config.get("Paths", "dssp_path"))
        # Exec process #
        dssp_out="dssp_"+str(os.getpid())+".out"
        #print("Execute ",os.path.join(dssp_path, "dssp"), pdb_file, os.path.join(dummy_dir,dssp_out))
        process = subprocess.check_output([os.path.join(dssp_path, "dssp"), pdb_file, os.path.join(dummy_dir,dssp_out)], stderr=subprocess.STDOUT)
        # Get DSSP object #
        dssp_obj = DSSP(os.path.join(dummy_dir,dssp_out))
        # Remove DSSP file #
        os.remove(os.path.join(dummy_dir,dssp_out))
    except:
        print "Failed ",os.path.join(dssp_path, "dssp"),pdb_file,os.path.join(dummy_dir,dssp_out)
        raise ValueError("Could not exec DSSP for %s" % pdb_file)

    return dssp_obj

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Get DSSP object #
    dssp_obj = get_dssp_obj(os.path.abspath(options.input_file), os.path.abspath(options.dummy_dir))

    # Output #
    if options.output_file is not None:
        dssp_obj.write(os.path.abspath(options.output_file))
    else:
        for line in dssp_obj._file_content:
            sys.stdout.write("%s\n" % line)
