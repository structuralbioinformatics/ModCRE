import os, sys, re
import ConfigParser
import hashlib
import optparse

# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.join(os.path.dirname(__file__),"../scripts"))
#scripts_path ="/homes/users/boliva/baldo_data/ModCRE/scripts"
# Append scripts path to python path #
sys.path.append(scripts_path)



# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Imports my functions #
import functions

#-------------#
# Options     #
#-------------#

def parse_options():
    '''
    This function parses the command line arguments and returns an optparse object.
    '''
    parser = optparse.OptionParser("python -m cisbp_sql_motifs_file  --ps proteins_file [ -o rootname -d folder -v]")

    parser.add_option("-m", action="store", default=None, type="string", dest="sql_file", help="SQL file of motifs (from CIS-BP)", metavar="{file}")
    parser.add_option("--models", action="store", default=None, type="string", dest="models_dir", help="Folder with the PWM of theoretical models", metavar="{directory}")
    parser.add_option("-o", action="store", type="string",default="CisBP",dest="root", help="CisBP file in SQL format with table of motifs (default CisBP)", metavar="{rootname}")
    parser.add_option("-d", action="store", type="string",default="sequences",dest="output_dir", help="Folder to store the sequences of TFs (default sequences)", metavar="{directory}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")

    (options, args) = parser.parse_args()

    if  options.sql_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Create output directory #
    if not os.path.exists(options.output_dir):
        os.makedirs(options.output_dir)


    #Read motifs
    nn_motifs = {}
    nn_tf_motifs = {}
    for line in functions.parse_file(os.path.abspath(options.sql_file)):
         m = re.search("\('(.+)', '(.+)', '.+', '.+', '.+', '(.+)', '.+', '.+'\),*", line)
         if m:
            if m.group(3) != "NULL":
               nn_motifs.setdefault(m.group(1),m.group(3))
               nn_tf_motifs.setdefault(m.group(2),set()).add(m.group(1))
         else:
            print "skip ",line
    nn_tf_mseq={}
    for tf,motifs in nn_tf_motifs.iteritems():         
        i=0
        for motif_id in motifs:
            sequence_motif = nn_motifs.get(motif_id)
            sequence_file = os.path.join(options.output_dir, tf + "." + str(i) + ".fa")
            nn_tf_mseq.setdefault(tf + "." + str(i) + ".fa", motif_id+".txt")
            i = i + 1
            if not os.path.exists(sequence_file):
                    functions.write(sequence_file, ">%s\n%s" % (tf, sequence_motif))

    nn_motifs_file = os.path.join(os.path.abspath(options.root+"_nn_motifs.dat"))
    in_motifs=open(nn_motifs_file,"w")
    in_motifs.write("#sequence_file\tmotif_file\n")

    for tf,pwm in nn_tf_mseq.iteritems():
        in_motifs.write("%s\t%s\n"%(tf,pwm))

    in_motifs.close()

    #Read models
    if options.models_dir is not None:
      nn_tf_model={}
      for  pwm_file in os.listdir(os.path.abspath(options.models_dir)):
           if pwm_file.endswith("meme"):
              data_model=pwm_file.split(":")
              tf=data_model[0]
              nn_tf_model.setdefault(tf,set()).add(pwm_file)
              
      nn_models_file = os.path.join(os.path.abspath(options.root+"_nn_models.dat"))

      in_models=open(nn_models_file,"w")
      in_models.write("#TF_ID\tPWM_theoretical\n")
      for tf,set_pwm in nn_tf_model.iteritems():
        for pwm in set_pwm:
            in_models.write("%s\t%s\n"%(tf,pwm))



      
