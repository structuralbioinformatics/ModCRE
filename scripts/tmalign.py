import os, sys, re
import ConfigParser
import copy
import numpy
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

class TMalign(object):
    """
    This class defines a {TMalign} object.

    """

    def __init__(self, file_content):
        self._file_content = file_content
        self._tm_scores = []
        self._matrix = numpy.identity(3, float)
        self._vector = numpy.zeros(3, float)
        self._query_alignment = None
        self._hit_alignment = None
        # Initialize #
        self._parse_file()

    def _parse_file(self):
        for line in self._file_content:
            # Capture TM-scores #
            m = re.search("^TM-score=\s*(\S+)", line)
            if m:
                self._tm_scores.append(float(m.group(1)))
            # Capture vector and matrix #
            m = re.search("^(0|1|2)\t\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)$", line)
            if m:
                self._matrix[int(m.group(1))][0] = float(m.group(3))
                self._matrix[int(m.group(1))][1] = float(m.group(4))
                self._matrix[int(m.group(1))][2] = float(m.group(5))
                self._vector[int(m.group(1))] = float(m.group(2))
            # Capture alignment #
            m = re.search("^([\w-]+)$", line)
            if m:
                if self._query_alignment == None:
                    self._query_alignment = m.group(1)
                else:
                    self._hit_alignment = m.group(1)

    def get_tm_scores(self):
        return self._tm_scores

    def get_matrix(self):
        return self._matrix

    def get_vector(self):
        return self._vector

    def get_query_alignment(self):
        return self._query_alignment

    def get_hit_alignment(self):
        return self._hit_alignment

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

    parser = optparse.OptionParser("python tmalign.py -i input_file [--dummy=dummy_dir -o output_file]")

    parser.add_option("-a", action="store", type="string", dest="pdb_file_a", help="PDB file", metavar="{filename}")
    parser.add_option("-b", action="store", type="string", dest="pdb_file_b", help="PDB file", metavar="{filename}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")

    (options, args) = parser.parse_args()

    if options.pdb_file_a is None or options.pdb_file_b is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_tmalign_obj(pdb_file_a, pdb_file_b, dummy_dir="/tmp"):
    """
    This function executes "tmalign" and returns a {TMalign}. Note that
    PDB B will be superimposed over PDB A.

    @input:
    pdb_file_a {filename}
    pdb_file_b {filename}
    dummy_dir {directory}

    @return:
    tmalign_obj {TMalign}

    """

    try:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        tmalign_path = os.path.join(src_path, config.get("Paths", "tmalign_path"))
        # If Mac OS X... #
        if sys.platform == "darwin":
            process = subprocess.check_output([os.path.join(tmalign_path, "TMalign"), pdb_file_a, pdb_file_b], stderr=subprocess.STDOUT)
        # Else... #
        else:
            process = subprocess.check_output([os.path.join(tmalign_path, "TMalign"), "-A", pdb_file_a, "-B", pdb_file_b], stderr=subprocess.STDOUT)
        # Get TMalign object #
        tmalign_obj = TMalign(process.split("\n"))
    except:
        raise ValueError("Could not exec TMalign for %s %s" % (pdb_file_a, pdb_file_b))
    return tmalign_obj

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Get TMalign object #
    tmalign_obj = get_tmalign_obj(os.path.abspath(options.pdb_file_a), os.path.abspath(options.pdb_file_b), os.path.abspath(options.dummy_dir))

    # Output #
    if options.output_file is not None:
        tmalign_obj.write(os.path.abspath(options.output_file))
    else:
        for line in tmalign_obj._file_content:
            sys.stdout.write("%s\n" % line)
