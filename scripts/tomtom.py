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

class Tomtom(object):
    """
    This class defines a {Tomtom} object.
    
    """

    def __init__(self, file_content):
        self._file_content = file_content
        self._hits = []
        self._query = None
        self._size  = 0
        # Initialize #
        self._parse_file()

    def _parse_file(self):
        rank=0
        for line in self._file_content:
            if line.startswith("Warning:"): continue
            if line.startswith("Query_ID"): continue
            if line.startswith("Processing"): continue
            if line.startswith("Estimat"): continue
            if line.startswith("Warning"): continue
            if line.startswith("#"): continue
            line = line.split("\t")
            if len(line) == 10:
                rank = rank + 1
                self._query = line.pop(0)
                hit, offset, p_value, e_value, q_value, overlap, query_sequence, hit_sequence, strand = line
                self._hits.append(TomtomHit(hit, offset, p_value, e_value, q_value, overlap, query_sequence, hit_sequence, strand, rank))
            else:
              if rank==0:
                sys.stderr.write("Could not execute TOMTOM first line %s\n"%(str(line)))
        self._size = rank

    def get_query(self):
        return self._query

    def get_size(self):
        return self._size

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

class TomtomHit(object):
    """
    This class defines a {TomtomHit} object.
    
    """

    def __init__(self, hit, offset, p_value, e_value, q_value, overlap, query_sequence, hit_sequence, strand, rank):
        self._hit = hit
        self._offset = int(offset)
        try:
          self._p_value = float(p_value)
        except:
          self._p_value = p_value
        try:
          self._e_value = float(e_value)
        except:
          self._e_value = e_value
        try:
          self._q_value = float(q_value)
        except:
          self._q_value = q_value
        self._overlap = int(overlap)
        self._query_sequence = query_sequence
        self._hit_sequence = hit_sequence
        self._strand = strand
        self._rank   = int(rank)
        
    def get_hit(self):
            return self._hit

    def get_offset(self):
            return self._offset

    def get_p_value(self):
           return self._p_value

    def get_e_value(self):
            return self._e_value

    def get_q_value(self):
           return self._q_value

    def get_overlap(self):
            return self._overlap

    def get_query_sequence(self):
            return self._query_sequence

    def get_hit_sequence(self):
            return self._hit_sequence

    def get_strand(self):
            return self._strand

    def get_rank(self):
            return self._rank

    def write(self,file_name):
            functions.write(file_name,"Hit %s Offset %d pValue %e eValue %e qValue %e Overlap %d Query sequence %s Hit sequence %s Strand %s Rank %d\n"%(self._hit,self._offset,self._p_value,self._e_value,self._q_value,self._overlap,self._query_sequence,self._hit_sequence,self._strand, self._rank))

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python tomtom.py -d database_file -i input_file [--dummy=dummy_dir -o output_file]")

    parser.add_option("-d", action="store", type="string", dest="database_file", help="Database file (in MEME format)", metavar="{filename}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (pwm in MEME format)", metavar="{filename}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")

    (options, args) = parser.parse_args()

    if options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_tomtom_obj(database_file, pwm_file, dummy_dir="/tmp"):
    """
    This function executes "fimo" from MEME package and returns a {Tomtom}.

    @input:
    database_file {filename}
    pwm_file {filename}
    dummy_dir {directory}

    @return:
    tomtom_obj {Tomtom}

    """
    try:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        meme_path =  config.get("Paths", "meme_path")
        if not os.path.exists(config.get("Paths", "meme_path")):
           meme_path = os.path.join(src_path, config.get("Paths", "meme_path"))
        try:
         # Exec process #
         process = subprocess.check_output([os.path.join(meme_path, "tomtom"), "-text", "-thresh", "1.0", pwm_file, database_file], stderr=subprocess.STDOUT)
         # Get Fimo object #
         tomtom_obj = Tomtom(process.split("\n"))
         return tomtom_obj
        except:
         sys.stdout.write("\t\t\t-- execute system tomtom\n")
         try:
           output_tomtom=os.path.join(dummy_dir,"tomtom"+str(os.getpid()))
           if not os.path.exists(dummy_dir): os.makedirs(dummy_dir)
           if not os.path.exists(output_tomtom): os.makedirs(output_tomtom)
           tomtom_file=os.path.join(output_tomtom,"tomtom.txt")
           sys.stdout.write("\t\t\t-- %s -text -thresh 1.0 %s %s >& %s"%(os.path.join(meme_path, "tomtom"),pwm_file,database_file,tomtom_file))
           os.system("%s -text -thresh 1.0 %s %s >& %s"%(os.path.join(meme_path, "tomtom"),pwm_file,database_file,tomtom_file))
           tomtom_output=open(tomtom_file,"r")
           try:
              tomtom_obj = Tomtom(tomtom_output.read().split("\n"))
           except Exception as e:
              raise ValueError("Could not exec Tomtom for  %s" % pwm_file)
           shutil.rmtree(output_tomtom)
           return tomtom_obj
         except:
           raise ValueError("Could not exec Tomtom for %s" % pwm_file)
    except:
        raise ValueError("Could not exec Tomtom for %s" % pwm_file)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Get Tomtom object #
    tomtom_obj = get_tomtom_obj(os.path.abspath(options.database_file), os.path.abspath(options.input_file), os.path.abspath(options.dummy_dir))

    # Output #
    if options.output_file is not None:
        tomtom_obj.write(os.path.abspath(options.output_file))
    else:
        for line in tomtom_obj._file_content:
            sys.stdout.write("%s\n" % line)
