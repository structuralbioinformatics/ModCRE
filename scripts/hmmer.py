import os, sys, re
import ConfigParser
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

# Imports jbonet's module #
from SBI.data import aminoacids1to3

#-------------#
# Classes     #
#-------------#

class Hmmer(object):
    """
    This class defines a {Hmmer} object.
    
    """

    def __init__(self, file_content):
        self._file_content = file_content
        self._query = None
        self._query_length = None
        self._hits = []
        # Initialize #
        self._parse_file()

    def _get_file(self):
        return self._file

    def _parse_file(self):
        for line in self._file_content:
            if line.startswith("#"): continue
            if "[No individual domains that satisfy reporting thresholds (although complete target did)]" in line:
                self._hits.pop(-1)
                continue
            m = re.search("Query:\s+(\S+)\s+\[M=(\d+)\]", line)
            if m:
                self._query = m.group(1)
                self._query_length = int(m.group(2))
            m = re.search(">>\s+(.+)", line)
            if m:
                self._hits.append(HmmerHit(m.group(1)))
            m = re.search("(\d+)\s+(\?|\!)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+.*\s+(\d+)\s+(\d+)\s+.*\s+\d+\s+\d+\S+", line)
            if m:
                domain = int(m.group(1))
                inclusion = m.group(2)
                score = float(m.group(3))
                bias = float(m.group(4))
                e_value = float(m.group(5))
                query_from = int(m.group(6))
                query_to = int(m.group(7))
                hit_from = int(m.group(8))
                hit_to = int(m.group(9))
                self._hits[-1].add_hsp(HmmerHSP(domain, inclusion, score, bias, e_value, query_from, query_to, hit_from, hit_to))
            m = re.search("== domain (\d+)", line)
            if m:
                domain = int(m.group(1))
            m = re.search("\S+\s+(\d+)\s([\w\-\.]+)\s(\d+)", line)
            if m:
                if len(re.findall("\w", m.group(2))) == int(m.group(3)) - int(m.group(1)) + 1:
                    if len(self._hits) > 0:
                        if len(self._hits[-1]._hsps) > 0:
                            self._hits[-1]._hsps[domain - 1].add_sequence(m.group(2).upper().replace(".", "-"))
            if "+" in line:
                if len(self._hits) > 0:
                    if len(self._hits[-1]._hsps) > 0:
                        if len(self._hits[-1]._hsps[domain - 1]._query_sequence) > len(self._hits[-1]._hsps[domain - 1]._hit_sequence):
                            self._hits[-1]._hsps[domain - 1].add_similarities(len(re.findall("\+", line)))

    def get_hits(self):
        return sorted(self._hits, key=lambda x: x.get_best_hsp().get_e_value())

    def write(self, file_name, filter_hits=False):
        for hmmer_hit_obj in self.get_hits():
            for hmmer_hsp_obj in hmmer_hit_obj.get_hsps():
                try:
                    string = "%s\t%s\t%s\t-1\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self._query, self._query_length, hmmer_hit_obj._hit, hmmer_hsp_obj.get_identities(), hmmer_hsp_obj.get_similarities(), hmmer_hsp_obj.get_gaps(), hmmer_hsp_obj.get_e_value(), hmmer_hsp_obj._query_sequence, hmmer_hsp_obj._hit_sequence, hmmer_hsp_obj.get_alignment_blocks())
                    if filter_hits:
                        if hmmer_hsp_obj._inclusion == "?": continue
                    functions.write(file_name, string)
                except:
                    pass

class HmmerHit(object):
    """
    This class defines a {HmmerHit} object.
    
    """

    def __init__(self, hit):
        self._hit = hit
        self._hsps = []
        
    def get_hit(self):
        return self._hit

    def add_hsp(self, hmmer_hsp_obj):
        self._hsps.append(hmmer_hsp_obj)

    def get_hsps(self):
        return sorted(self._hsps, key=lambda x: x.get_e_value())

    def get_best_hsp(self):
        for hmmer_hsp_obj in self.get_hsps():
            return hmmer_hsp_obj

class HmmerHSP(object):
    """
    This class defines a {HmmerHSP} object.
    
    """

    def __init__(self, domain, inclusion, score, bias, e_value, query_from, query_to, hit_from, hit_to):
        self._domain = domain
        self._inclusion = inclusion
        self._score = score
        self._bias = bias
        self._e_value = e_value
        self._query_from = query_from
        self._query_to = query_to
        self._hit_from = hit_from
        self._hit_to = hit_to
        self._query_sequence = ""
        self._hit_sequence = ""
        self._similarities = 0

    def add_sequence(self, sequence):
        if len(self._query_sequence) == len(self._hit_sequence):
            self.add_query_sequence(sequence)
        else:
            self.add_hit_sequence(sequence)
        
    def add_query_sequence(self, sequence):
        self._query_sequence += sequence

    def add_hit_sequence(self, sequence):
        self._hit_sequence += sequence

    def add_similarities(self, similarities):
        self._similarities += similarities

    def get_e_value(self):
        return self._e_value

    def get_query_sequence(self):
        return self._query_sequence

    def get_hit_sequence(self):
        return self._hit_sequence

    def get_identities(self):
        identities = 0

        for i in range(len(self._query_sequence)):
            if self._query_sequence[i] == self._hit_sequence[i]:
                identities += 1

        return identities

    def get_similarities(self):
        return self.get_identities() + self._similarities

    def get_gaps(self):
        return len(re.findall("\-", self._query_sequence)) + len(re.findall("\-", self._hit_sequence))

    def get_alignment_blocks(self):
        query_blocks = [[]]
        hit_blocks = [[]]
        query_positions = range(self._query_from, self._query_to + 1)
        hit_positions = range(self._hit_from, self._hit_to + 1)
        for i in range(len(self._query_sequence)):
            if aminoacids1to3.has_key(self._query_sequence[i]) and aminoacids1to3.has_key(self._hit_sequence[i]):
                query_blocks[-1].append(query_positions[0])
                hit_blocks[-1].append(hit_positions[0])
            else:
                if len(query_blocks[-1]) > 0:
                    query_blocks.append([])
                if len(hit_blocks[-1]) > 0:
                    hit_blocks.append([])
            if aminoacids1to3.has_key(self._query_sequence[i]):
                query_positions.pop(0)
            if aminoacids1to3.has_key(self._hit_sequence[i]):
                hit_positions.pop(0)

        return ";".join(["%s:%s,%s:%s" % (query_blocks[i][0], hit_blocks[i][0], query_blocks[i][-1], hit_blocks[i][-1]) for i in range(len(query_blocks))])

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python hmmer.py -d database_file -i input_file [--dummy=dummy_dir -f -o output_file]")

    parser.add_option("-d", action="store", type="string", dest="database_file", help="Database file (in FASTA format)", metavar="{filename}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-f", "--filter", default=False, action="store_true", dest="filter_hits", help="Filter hits under inclusion threshold (default = False)", metavar="{boolean}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in HMM format)", metavar="{filename}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")

    (options, args) = parser.parse_args()

    if options.database_file == None or options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_pfam(database_file, input_file, dummy_dir="/tmp"):
    """
    This function executes "hmmscan" from HMMER package and returns a the set of families.

    @input:
    database_file {string}
    input_file {string}
    dummy_dir {string}

    @return:
    pfam_families {set}
    """
    try:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        hmmer_path = os.path.join(src_path, config.get("Paths", "hmmer_path"))
        # Exec process #
        process = subprocess.check_output([os.path.join(hmmer_path, "hmmscan"), database_file, input_file], stderr=subprocess.STDOUT)
        # Get BLAST object #
        families=set()
        for line in process.split("\n"):
            if line.startswith(">>"):
                families.add(line.split()[1])
    except:
        os.system("%s %s %s"%(os.path.join(hmmer_path, "hmmscan"), database_file, input_file))
        raise ValueError("Could not exec hmmscan for %s" % input_file)

    return families


def get_hmmer_obj(database_file, input_file, dummy_dir="/tmp"):
    """
    This function executes "hmmsearch" from HMMER package and returns a {HMMER}.

    @input:
    database_file {string}
    input_file {string}
    dummy_dir {string}

    @return:
    hmmer_obj {Hmmer}

    """

    try:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        hmmer_path = os.path.join(src_path, config.get("Paths", "hmmer_path"))
        # Exec process #
        process = subprocess.check_output([os.path.join(hmmer_path, "hmmsearch"), input_file, database_file], stderr=subprocess.STDOUT)
        # Get BLAST object #
        hmmer_obj = Hmmer(process.split("\n"))
    except:
        raise ValueError("Could not exec hmmsearch for %s" % input_file)

    return hmmer_obj

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()
    if not os.path.exists(os.path.abspath(options.dummy_dir)):
       os.makedirs(os.path.abspath(options.dummy_dir))

    # Get HMMER object #
    hmmer_obj = get_hmmer_obj(os.path.abspath(options.database_file), os.path.abspath(options.input_file), os.path.abspath(options.dummy_dir))

    # Output #
    if options.output_file is not None:
        # Write output #
        dummy_file = os.path.abspath(os.path.join(options.dummy_dir, str(os.getpid()) + ".txt"))
        hmmer_obj.write(dummy_file, options.filter_hits)
        shutil.copy(dummy_file, options.output_file)
        os.remove(dummy_file)
    else:
        hmmer_obj.write(options.output_file, options.filter_hits)
