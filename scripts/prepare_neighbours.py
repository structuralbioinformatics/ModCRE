import os, sys, re
import ConfigParser
import optparse
import subprocess

# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.dirname(__file__)+"../../scripts")

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


def fileExist(file):
    '''
    Check existing files
    '''
    if file is not None:
        return os.path.exists(file) and os.path.isfile(file)
    else:
        return False

#-------------#
# Options     #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python prepare_neighbours.py -i tf_motifs_file --pwm=pwm_dir --pdb=pbm_dir [-o output_dir]")

    parser.add_option("--dummy", default="./tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="File of correspondences between FastA sequences of TFs in PDB and PWMs of motifs in CisBP", metavar="{filename}")
    parser.add_option("-o", action="store", default="output_neighbours", type="string", dest="output_file", help="Output file (if single file as input) or directory name (default = 'output_pwm')", metavar="{filename}")
    parser.add_option("--pbm", action="store", type="string", dest="pbm_dir", help="PBM directory (i.e. output dir from pbm.py)", metavar="{directory}")
    parser.add_option("--pwm", action="store", type="string", dest="pwm_dir", default=None, help="PWM directory with PWM motifs from CisBP", metavar="{directory}")
    parser.add_option("-t", action="store", type="string", dest="tfs_file", help="TFs file (from CIS-BP; i.e. cisbp_1.02.tfs.sql)", metavar="{filename}")
    parser.add_option("-f", action="store", type="string", dest="family_file", help="TFs Family file (from CIS-BP; i.e. cisbp_1.02.tf_families.sql)", metavar="{filename}")
    parser.add_option("--use_CisBP_family", default=False, action="store_true", dest="cisbp_family",help="Use family codes from CisBP instead of common family names of TFs (default = False)", metavar="{boolean}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")

   
    (options, args) = parser.parse_args()

    if (options.input_file is None ) or options.pbm_dir is None or options.pwm_dir is None or options.family_file is None or options.tfs_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()
    cdhit=config.get("Paths", "cd-hit")
    mmseqs=config.get("Paths", "mmseqs")
    verbose=options.verbose

    # Initialize #
    if not os.path.exists(options.dummy_dir): os.makedirs(options.dummy_dir)
    tmp    =options.dummy_dir
    pwm_dir=options.pwm_dir
    pwm_out=options.output_file + "/pwms"
    seq_dir=options.output_file + "/fasta"
    if not os.path.exists(options.output_file): 
       os.makedirs(options.output_file)
       os.makedirs(pwm_out)
       os.makedirs(seq_dir)
    if not os.path.exists(pwm_out):os.makedirs(pwm_out)
    if not os.path.exists(seq_dir):os.makedirs(seq_dir)
 
    #Read CisBP families
    cisbp_families={}
    for line in functions.parse_file(os.path.abspath(options.family_file)):
        m = re.search("\('(.+)', '(.+)', '.+', .+, .+\)", line)
        if m:           
            cisbp_families.setdefault(m.group(1).upper(),set()).add(m.group(2))
    fo=open("cisbp_families.txt","w")
    for cis,fam in cisbp_families.iteritems():
        fo.write("%s\t%s\n"%(cis, ";".join([str(x) for x in fam])))
    fo.close()

    #Read CisBP TFs
    tf_species  = {}
    tf_families = {}
    tf_families_cisbp = {}
    for line in functions.parse_file(os.path.abspath(options.tfs_file)):
         m = re.search("\('(.+)', '(.+)', '.+', '.+', '.+', '(.+)', '[DIN]'\),*", line)
         if m:
             tf_species.setdefault(m.group(1), set()).add(m.group(3).replace("_", " ").upper())
             tf_families.setdefault(m.group(1).upper(), set()).add(m.group(2).upper())
             tf_families_cisbp.setdefault(m.group(1).upper(), set()).update(cisbp_families[m.group(2).upper()])
    fo=open("families.txt","w")
    for tf,fam in tf_families.iteritems():
        fo.write("%s\t%s\n"%(tf, ";".join([str(x) for x in fam]) ))
    fo.close()
    fo=open("families_cisbp.txt","w")
    for tf,fam in tf_families_cisbp.iteritems():
        fo.write("%s\t%s\n"%(tf, ";".join([str(x) for x in fam])))
    fo.close()

    #Create TF-Motif correspondence   
    tf_has_motif={}
    fd=open(options.input_file,"r")
    for line in fd:
        tf,motif=line.split()
        tf_has_motif.setdefault(tf,motif)
    fd.close()

    #Create sequence database
    fasta_file=seq_dir+"/tf.fa"
    fasta_db  =seq_dir+"/tf.fa.db"
    search_db =seq_dir+"/tf.search.db"
    compared  =seq_dir+"/tf.compared"
    fa=open(fasta_file,"w")
    for tf,motif in tf_has_motif.iteritems():
      try:
        pwm_file = options.pwm_dir+"/"+motif
        tf_file  = options.pbm_dir+"/sequences/"+tf
        if fileExist(pwm_file) and fileExist(tf_file):
           pwm = PWM.nMSA(pwm_file,motif,"txt")
           name= pwm_out+"/"+os.path.basename(pwm_file).rstrip(".txt")+".pwm"
           pwm.write(name,"pwm")
           name= pwm_out+"/"+os.path.basename(pwm_file).rstrip(".txt")+".meme"
           pwm.write(name,"meme")
           seq = open(tf_file,"r")
           for line in seq:
             fa.write("%s\n"%line.strip())
           seq.close()
        else:
           if verbose: sys.stderr.write("Files TF %s MOTIF %s are not found\n"%(tf,motif))
      except :
        sys.stderr.write("Error to parse TF %s MOTIF %s\n"%(tf,motif))
    
    if verbose: sys.stdout.write("mmseqs createdb  %s %s >& createdb.log\n"%(fasta_file, fasta_db))
    os.system("mmseqs createdb  %s %s >& createdb.log"%(fasta_file, fasta_db))
    if verbose: sys.stdout.write("mmseqs createindex %s %s >& createindex.log \n"%( fasta_db,tmp))
    os.system("mmseqs createindex %s %s >& createindex.log"%( fasta_db,tmp))
    if verbose: sys.stdout.write("mmseqs search %s %s %s %s --threads 1 -s  7.5 --max-seq-id 1.0 --num-iterations 4 >& search.log\n"%(fasta_db,fasta_db,search_db,tmp))
    os.system("mmseqs search %s %s %s %s --threads 1 -s 7.5 --max-seq-id 1.0 --num-iterations 4 >& search.log"%(fasta_db,fasta_db,search_db,tmp))
    process = subprocess.Popen( [mmseqs, "convertalis", fasta_db, fasta_db, search_db,compared ],  stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    process.communicate()

    for d in xrange(10,100,5):
        fasta_hit     =seq_dir+"/tf"+str(d)+".fa"
        fasta_hit_db  =seq_dir+"/tf"+str(d)+".fa.db"
        search_hit_db =seq_dir+"/tf"+str(d)+".search.db"
        compared_hit  =seq_dir+"/tf"+str(d)+".compared"
        word_length = str(2)
        threshold   = str("%4.2f"%(float(d)/100)).split()[0]
        if verbose: sys.stdout.write("Generate dataset %s\n"%fasta_hit)
        process = subprocess.Popen( [cdhit, "-i", fasta_file , "-o", fasta_hit, "-n", word_length , "-c", threshold ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        process.communicate()
        #if verbose: sys.stdout.write("mmseqs createdb  %s %s >& createdb.log\n"%(fasta_hit, fasta_hit_db))
        #os.system("mmseqs createdb  %s %s >& createdb.log"%(fasta_hit, fasta_hit_db))
        #if verbose: sys.stdout.write("mmseqs createindex %s %s >& createindex.log \n"%( fasta_hit_db,tmp))
        #os.system("mmseqs createindex %s %s >& createindex.log"%( fasta_hit_db,tmp))
        #if verbose: sys.stdout.write("mmseqs search %s %s %s %s --threads 1 -s  7.5 --max-seq-id 1.0 --num-iterations 4 >& search.log\n"%(fasta_hit_db,fasta_hit_db,search_hit_db,tmp))
        #os.system("mmseqs search %s %s %s %s --threads 1 -s 7.5 --max-seq-id 1.0 --num-iterations 4 >& search.log"%(fasta_hit_db,fasta_hit_db,search_hit_db,tmp))
        #process = subprocess.Popen( [mmseqs, "convertalis", fasta_hit_db, fasta_hit_db, search_hit_db,compared_hit ],  stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        #process.communicate()

        
