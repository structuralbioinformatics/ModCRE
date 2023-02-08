import os, sys, re
import ConfigParser
import copy
import optparse
import shutil
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

class X3DNA(object):
    """
    This class defines a {X3DNA} object.

    """

    def __init__(self, file_name):
        self._file = file_name
        self._file_content = []
        self._sequence = ""
        self._residues = {}
        self._basepairs = {}
        self._dinucleotides = {}
        self._inverse_dinucleotides = {}
        self._helix = {}
        self._helix_dinucleotides = {}
        # Initialize #
        self._parse_file()
        self._initialize_dinucleotides()
        self._initialize_helix_dinucleotides()

    def _parse_file(self):
        for line in functions.parse_file(self._file):
            # Get base-pair residues #
            m = re.search("(\d+)\s+\S\s\S{4}>(\S):\D*(\d+)_:\[\S{3}\](\w)\S{5}\w\[\S{3}\]:\D*(\d+)_:(\S)<\S{4}", line)
            if m:
                basepair = int(m.group(1))
                fwd_pdb_chain = m.group(2)
                fwd_residue_num = int(m.group(3))
                self._sequence += m.group(4)
                rev_pdb_chain = m.group(6)
                rev_residue_num = int(m.group(5))
                self._residues[(fwd_pdb_chain, fwd_residue_num)] = basepair
                self._residues[(rev_pdb_chain, rev_residue_num)] = basepair
                self._basepairs[basepair] = ((fwd_pdb_chain, fwd_residue_num), (rev_pdb_chain, rev_residue_num))
            m = re.search("^#####\s+Helix\s+#(\d+)\s+\(\d+\):\s+(\d+)\s+\-\s+(\d+)", line)
            if m:
                helix = m.group(1)
                first = int(m.group(2))
                last = int(m.group(3))
                self._helix.setdefault(helix, range(first, last + 1))
            # Add line to file content #
            if line != "":
               self._file_content.append(line)
            

    def _initialize_dinucleotides(self):
        basepairs = [key for key in sorted(self._basepairs.iterkeys())]
        while len(basepairs) > 1:
            i = basepairs.pop(0)
            if self.has_dinucleotide(i): continue
            self._dinucleotides.setdefault(i, (i, basepairs[0]))
            self._inverse_dinucleotides.setdefault((i, basepairs[0]),i)

    def _initialize_helix_dinucleotides(self):
        for helix in self._helix:
            self._helix_dinucleotides.setdefault(helix, [])
            basepairs = copy.copy(self.get_helix_basepairs(helix))

            while len(basepairs) > 1:
                i = basepairs.pop(0)
                self._helix_dinucleotides[helix].append(i)

    def has_residue(self, pdb_chain, residue_num):
        return (pdb_chain, residue_num) in self._residues

    def has_basepair(self, basepair):
        return basepair in self._basepairs

    def has_dinucleotide(self, dinucleotide):
        return dinucleotide in self._dinucleotides

    def has_helix(self, helix):
        return helix in self._helix

    def helix_has_residue(self, helix, pdb_chain, residue_num):
        if self.has_residue(pdb_chain, residue_num) and self.has_helix(helix):
            return self._residues[(pdb_chain, residue_num)] in self.get_helix_basepairs(helix)

        return False

    def get_residue_basepair(self, pdb_chain, residue_num):
        if self.has_residue(pdb_chain, residue_num):
            return copy.copy(self._residues[(pdb_chain, residue_num)])

        return None

    def get_basepair(self, basepair):
        if self.has_basepair(basepair):
            return copy.copy(self._basepairs[basepair])

        return None

    def get_basepairs(self):
        return copy.copy(self._basepairs)

    def get_helix_basepairs(self, helix):
        if self.has_helix(helix):
            return copy.copy(self._helix[helix])

        return None

    def get_dinucleotide(self, dinucleotide):
        if self.has_dinucleotide(dinucleotide):
            return copy.copy(self._dinucleotides[dinucleotide])

        return None
  
    def get_inverse_dinucleotide(self, basepair_1, basepair_2):
        if self._inverse_dinucleotides.has_key((basepair_1,basepair_2)):
           return copy.copy(self._inverse_dinucleotides[(basepair_1,basepair_2)])
        return None
 
    def get_sequence(self):
        return copy.copy(self._sequence)

    def get_dinucleotides(self):
        return copy.copy(self._dinucleotides)

    def get_inverse_dinucleotides(self):
        return copy.copy(self._inverse_dinucleotides)


    def get_helix_dinucleotides(self, helix):
        if self.has_helix(helix):
            return copy.copy(self._helix_dinucleotides[helix])

        return None

    def get_basepair_dinucleotides(self, basepair):
        dinucleotides = []

        for dinucleotide in self.get_dinucleotides():
            if basepair in self.get_dinucleotide(dinucleotide):
                dinucleotides.append(dinucleotide)
        
        return copy.copy(dinucleotides)

    def get_basepair_helix(self, basepair):
        if self.has_basepair(basepair):
            for helix in self._helix:
                if basepair in self.get_helix_basepairs(helix):
                    return helix

        return None

    def get_dna_helices(self):
        return copy.copy(self._helix)

    def get_nucleotide_sequence(self, A, B):
        if self.has_basepair(A) and self.has_basepair(B):
            return self._sequence[A - 1:B]

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

    parser = optparse.OptionParser("python x3dna.py -i input_file [--dummy=dummy_dir -o output_file]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in PDB format)", metavar="{filename}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")

    (options, args) = parser.parse_args()

    if options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_x3dna_obj(pdb_file, dummy_dir="/tmp"):
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
        x3dna_path = os.path.join(src_path, config.get("Paths", "x3dna_path"))
        os.environ['X3DNA'] = x3dna_path[:-4]
        # Get current working directory #
        cwd = os.getcwd()
        # Create tmp directory #
        tmp = os.path.join(dummy_dir, str(os.getpid()))
        if not os.path.exists(tmp): os.makedirs(tmp)
        # Change directory #
        os.chdir(tmp)
        # Exec process #
        process = subprocess.check_output([os.path.join(x3dna_path, "find_pair"), pdb_file, "3dna.out"], stderr=subprocess.STDOUT, env=os.environ)
        # Get X3DNA object #
        x3dna_obj = X3DNA("3dna.out")
        # Return to original directory #
        os.chdir(cwd)
        # Erase tmp directory #
        shutil.rmtree(tmp)
    except Exception as e:
        raise ValueError("Could not exec X3DNA for %s with error %s" % (pdb_file,e))

    return x3dna_obj

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Get X3DNA object #
    x3dna_obj = get_x3dna_obj(os.path.abspath(options.input_file), os.path.abspath(options.dummy_dir))

    # Output #
    if options.output_file is not None:
        x3dna_obj.write(os.path.abspath(options.output_file))
    else:
        for line in x3dna_obj._file_content:
            sys.stdout.write("%s\n" % line)
