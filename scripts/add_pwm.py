import os, sys, re
import ConfigParser
import optparse

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
import pwm_pbm as PWM

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python add_pwm.py -a input_pwm_file -b  input_pwm_file  [-o output_pwm_file --fa format file A --fb format file B --fmt format output --name motif_name]")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-a", action="store", default=None, type="string", dest="input_pwm_A", help="PWM file (A) in any format 'pwm' 'msa' or 'meme')", metavar="{filename}")
    parser.add_option("-b", action="store", default=None, type="string", dest="input_pwm_B", help="PWM file (B) in any format 'pwm' 'msa' or 'meme')", metavar="{filename}")
    parser.add_option("-o", action="store", default="output_pwm", type="string", dest="output_file", help="Output file (default = 'output_pwm')", metavar="{filename}")
    parser.add_option("--name", action="store", default=None, type="string", dest="name_motif", help="Name of output motif (default is None and it merges the input)", metavar="{string}")
    parser.add_option("--fmt", action="store", default="meme", type="string", dest="format_pwm", help="Format of the output PWM (default is 'meme')", metavar="{string}")
    parser.add_option("--fa", action="store", default="meme", type="string", dest="format_A", help="Format of the input file (A) PWM (default is 'meme')", metavar="{string}")
    parser.add_option("--fb", action="store", default="meme", type="string", dest="format_B", help="Format of the input file (B) PWM (default is 'meme')", metavar="{string}")
    parser.add_option("-p", "--protein", default=False, action="store_true", dest="protein", help="PWMs are from protein sequences (default = False, therefore sequences are polymers of nucleotides)")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")
 
    (options, args) = parser.parse_args()

    if (options.input_pwm_A is None or options.input_pwm_A is None):
        parser.error("missing arguments: type option \"-h\" for help")
     
    return options



#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Read input files
    if  options.verbose:sys.stdout.write("\t--Get PWM file A %s ...\n"% options.input_pwm_A)
    if  options.protein:
        msa_obj_a = PWM.pMSA(options.input_pwm_A,None,options.format_A)
    else:
        msa_obj_a = PWM.nMSA(options.input_pwm_A,None,options.format_A)
    if  options.verbose:sys.stdout.write("\t--Get PWM file B %s ...\n"% options.input_pwm_B)
    if  options.protein:
        msa_obj_b = PWM.pMSA(options.input_pwm_B,None,options.format_B)
    else:
        msa_obj_b = PWM.nMSA(options.input_pwm_B,None,options.format_B)

    # Add PWMs
    msa_obj = msa_obj_a + msa_obj_b
    if options.name_motif is not None:
       msa_obj.set_motif(options.name_motif)

    # Write PWM #
    if  options.verbose:sys.stdout.write("\t--Write PWM ...\n")
    format_file = options.format_pwm
    pwm_file    = options.output_file + "." + format_file
    logo_file   = options.output_file + ".logo"
    logo_gapped_file   = options.output_file + ".gapped.logo"
    if not os.path.exists(pwm_file):
       if  options.verbose:sys.stdout.write("\t\t--PWM in %s format...\n"%format_file)
       msa_obj.write(pwm_file, option=format_file)
    if not os.path.exists(logo_file+".fwd.png") or  not os.path.exists(logo_file+".rev.png"):
       if not os.path.exists(options.dummy_dir): os.makedirs(options.dummy_dir)
       if  options.verbose:sys.stdout.write("\t\t--Logos...\n")
       if  options.protein:
         PWM.write_protein_logo(msa_obj,logo_file, logo_gapped_file, options.dummy_dir)
       else:
         PWM.write_logo(msa_obj,logo_file, options.dummy_dir)


    
