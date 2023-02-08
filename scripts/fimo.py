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

# Import my modules #
import triads

#-------------#
# Classes     #
#-------------#

class Fimo(object):
    """
    This class defines a {FIMO} object.
    
    """

    def __init__(self, file_content):
        self._file_content = file_content
        self._hits = []
        # Initialize #
        self._parse_file()

    def _parse_file(self):

        #if config.get("Parameters", "fimo_pvalue_threshold") is not None:
        #   fimo_pvalue_threshold=float(config.get("Parameters", "fimo_pvalue_threshold"))
        #else:
        #   fimo_pvalue_threshold=1.0
        # We collect all p-values selected by fimo without exception
        fimo_pvalue_threshold=1.0
        for linefile in self._file_content:
          if linefile.startswith("#"): continue
          try:
             line = linefile.split("\t")
             if len(line) == 9:
                 read_condition=False
                 if line[-3] != "":
                    if float(line[-3]) < fimo_pvalue_threshold: read_condition=True
                 if read_condition:
                      strand=line[4]
                      if strand == "+": sequence=line[-1]
                      else: sequence = triads.get_complementary_dna_sequence(line[-1])
                      self._query = line[1]
                      self._hits.append(FimoHit(line[0], int(line[2]), int(line[3]), strand, float(line[5]), float(line[6]), sequence))
             if len(line) == 10:
                 read_condition=False
                 if line[-2] != "":
                    if float(line[-2]) < fimo_qvalue_threshold: read_condition=True
                 if line[-3] != "":
                    if float(line[-3]) < fimo_pvalue_threshold: read_condition=True
                 if read_condition:
                      strand=line[5]
                      if strand == "+": sequence=line[-1]
                      else: sequence = triads.get_complementary_dna_sequence(line[-1])
                      self._query = line[1]
                      self._hits.append(FimoHit(line[0], int(line[3]), int(line[4]), strand, float(line[-4]), float(line[-3]), sequence))
          except:
             sys.stderr.write("WARNING FIMO: skip %s\n"%(linefile))

    def get_query(self):
        return self._query

    def has_hit(self, hit_name):
        for hit_obj in self._hits:
            if hit_name == hit_obj.get_hit():
                return True

        return False

    def get_hits(self, sort=False):
        if sort:
            return sorted(self._hits, key=lambda x: x.get_p_value())

        return self._hits

    def get_hit(self, hit_name):
        if self.has_hit(hit_name):
            for hit_obj in self._hits:
                if hit_name == hit_obj.get_hit():
                    return hit_obj

        return None

    def write(self, file_name):
        for line in self._file_content:
            functions.write(file_name, line)

class FimoHit(object):
    """
    This class defines a {FimoHit} object.
    
    """

    def __init__(self, hit, start, end, strand, score, p_value, sequence):
        self._hit = hit
        self._start = start
        self._end = end
        self._strand = strand
        self._score = score
        self._p_value = p_value
        self._sequence = sequence
        
    def get_hit(self):
            return self._hit

    def get_start(self):
            return self._start

    def get_strand(self):
            return self._strand

    def get_end(self):
            return self._end

    def get_score(self):
            return self._score

    def get_p_value(self):
            return self._p_value

    def get_sequence(self):
            return self._sequence

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python fimo.py -d database_file -i input_file [--dummy=dummy_dir -o output_file]")

    parser.add_option("-d", action="store", type="string", dest="database_file", help="Database file (in MEME format)", metavar="{filename}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in FASTA format)", metavar="{filename}")
    parser.add_option("--ft", action="store", type="float", dest="fimo_pvalue_threshold",  default=None, help="P-value threhold for fimo matches", metavar="{float}")
    parser.add_option("--max", action="store", type="int",dest="max_stored_matches",  default=None, help="Maximum number of matches stored", metavar="{integer}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")

    (options, args) = parser.parse_args()

    if options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_fimo_obj(database_file, fasta_file,fimo_pvalue_threshold=None, max_stored_matches=None, dummy_dir="/tmp"):
    """
    This function executes "fimo" from MEME package and returns a {FIMO}.

    @input:
    database_file {filename}
    fasta_file {filename}
    dummy_dir {directory}

    @return:
    fimo_obj {FIMO}

    """

    try:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        meme_path =  config.get("Paths", "meme_path")
        if not os.path.exists(config.get("Paths", "meme_path")):
           meme_path = os.path.join(src_path, config.get("Paths", "meme_path"))
        if fimo_pvalue_threshold is None:
           if config.get("Parameters", "fimo_pvalue_threshold") is not None:
              fimo_pvalue_threshold=float(config.get("Parameters", "fimo_pvalue_threshold"))
           else:
              fimo_pvalue_threshold=1
        if not os.path.exists(dummy_dir): os.makedirs(dummy_dir)
        # Exec process #
        log_file=os.path.join(dummy_dir,"fimo"+str(os.getpid())+".log")
        output_fimo=os.path.join(dummy_dir,"fimo"+str(os.getpid()))
        n=1
        while os.path.exists(output_fimo):
              output_fimo=output_fimo+str(n)
              n = n + 1

        try:
           fimo_file=os.path.join(output_fimo,"fimo.txt")
           if max_stored_matches is None:
              process = subprocess.check_output([os.path.join(meme_path, "fimo"), "-o",output_fimo, "--thresh", str(fimo_pvalue_threshold) , database_file, fasta_file], stderr=subprocess.STDOUT)
           else:
              process = subprocess.check_output([os.path.join(meme_path, "fimo"), "-o",output_fimo, "--thresh", str(fimo_pvalue_threshold) ,"--max-stored-scores", str(max_stored_matches), database_file, fasta_file], stderr=subprocess.STDOUT)
           fimo_output=open(fimo_file,"r")
           # Get Fimo object #
           fimo_obj = Fimo(fimo_output.read().split("\n"))
           shutil.rmtree(output_fimo)
           return fimo_obj
        except:
         try:
           if not os.path.exists(output_fimo): os.makedirs(output_fimo)
           fimo_file=os.path.join(output_fimo,"fimo.txt")
           sys.stdout.write("\t-- execute system '--text' fimo option\n")
           if max_stored_matches is None:
              os.system("%s -o %s --text --thresh %s %s %s > %s\n"%(os.path.join(meme_path, "fimo"), output_fimo, str(fimo_pvalue_threshold),database_file, fasta_file,fimo_file))
           else:
              os.system("%s -o %s --text --thresh %s --max-stored-scores %s %s %s > %s\n"%(os.path.join(meme_path, "fimo"), output_fimo, str(fimo_pvalue_threshold),str(max_stored_matches),database_file, fasta_file,fimo_file))
           fimo_output=open(fimo_file,"r")
           # Get Fimo object #
           fimo_obj = Fimo(fimo_output.read().split("\n"))
           shutil.rmtree(output_fimo)
           return fimo_obj
         except:
           raise ValueError("Could not exec FIMO for %s" % fasta_file)
    except:
        raise ValueError("Could not exec FIMO for %s" % fasta_file)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Get FIMO object #
    fimo_obj = get_fimo_obj(os.path.abspath(options.database_file), os.path.abspath(options.input_file), options.fimo_pvalue_threshold,options.max_stored_matches, os.path.abspath(options.dummy_dir))

    # Output #
    if options.output_file is not None:
        fimo_obj.write(os.path.abspath(options.output_file))
    else:
        for line in fimo_obj._file_content:
            sys.stdout.write("%s\n" % line)

