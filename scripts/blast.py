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
from SBI.external.blast import blast_parser

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python blast.py -d database_file -i input_file [--dummy=dummy_dir -f -o output_file]")

    parser.add_option("-d", action="store", type="string", dest="database_file", help="Database file (in FASTA format)", metavar="{filename}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-f", "--filter", default=False, action="store_true", dest="filter_hits", help="Filter twilight zone hits (default = False)", metavar="{boolean}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in FASTA format)", metavar="{filename}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")

    (options, args) = parser.parse_args()

    if options.database_file is None or options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_blast_obj(database_file, input_file, dummy_dir="/tmp"):
    """
    This function executes "blastpgp" from BLAST package and returns a {BlastOutput}.

    @input:
    database_file {string}
    input_file {string}
    dummy_dir {string}

    @return:
    blast_obj {BlastOutput}

    """

    try:
        # Initialize #
        src_path = config.get("Paths", "src_path")
        blast_path = os.path.join(src_path, config.get("Paths", "blast_path"))
        blast_out = os.path.join(dummy_dir, os.path.basename(input_file).strip()+str(os.getpid()) + ".txt")
        # Exec process #
        print "%s -query %s -db %s -out %s -outfmt 5 > %s"%(os.path.join(blast_path, "blastp"),input_file,database_file,blast_out,blast_out+".log")
        process = subprocess.check_output([os.path.join(blast_path, "blastp"), "-query", input_file, "-db", database_file, "-out", blast_out, "-outfmt", "5"], stderr=subprocess.STDOUT)
        # Get BLAST object #
        blast_obj = blast_parser.parse_blast(query_sequence=None, blast_output_file=blast_out, selfHit=True, hitIDformat="all")
        # Remove BLAST file #
        #os.remove(blast_out)
    except:
        raise ValueError("Could not exec blastp for %s" % input_file)

    return blast_obj

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Get BLAST object #
    blast_obj = get_blast_obj(os.path.abspath(options.database_file), os.path.abspath(options.input_file), os.path.abspath(options.dummy_dir))

    # Twilight zone #
    tz_parameter = 0
    tz_type = None
    if options.filter_hits:
        tz_parameter = int(config.get("Parameters", "twilight_zone_parameter"))
        tz_type = config.get("Parameters", "twilight_zone_type")

    # Output #
    if options.output_file is not None:
        # Write output #
        dummy_file = os.path.abspath(os.path.join(options.dummy_dir, str(os.getpid()) + ".txt"))
        functions.write(dummy_file, blast_obj.str_compacted_blast(tz_parameter=tz_parameter, tz_type=tz_type))
        shutil.copy(dummy_file, options.output_file)
        os.remove(dummy_file)
    else:
        print(blast_obj.str_compacted_blast(tz_parameter=tz_parameter, tz_type=tz_type))
