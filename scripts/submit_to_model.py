import os, sys, re
import ConfigParser
import optparse
import shutil
import subprocess

# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.dirname(__file__))

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Imports my functions #
import functions
import SeqIO as SEQ

# Imports jbonet's module #
from SBI.structure import PDB
from SBI.structure.chain import Chain

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python submit_to_model.py -p pdb_dir -i fasta_file [--dummy=dummy_dir -o output_dir --parallel -v]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="{directory}")
    parser.add_option("-p", "--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (where \"all\" and \"biounits\" directories are placed)", metavar="{directory}")
    parser.add_option("-j","--parallel", default=False, action="store_true", dest="parallel", help="Submit JOBS to Queues in parallel. The program stops, you need to re-start with the requested value (default = False)", metavar="{boolean}")
    parser.add_option("-i", action="store",  type="string", dest="fasta_file", help="TFs of FASTA file with sequences to be modelled", metavar="{filename}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)", metavar="{boolean}")
    parser.add_option("--n-model", default=1, action="store", type="int", dest="n_model", help="Number of models per template (default = 1)", metavar="{int}")
    parser.add_option("--n-total", action="store", type="int", dest="n_total", default=None, help="Total number of models per execution of this program (if not \"None\", automatically sets the \"--n-model\" parameter; supersedes option \"--n-model\"; default = None)", metavar="{int}")
    parser.add_option("--force", default=False, action="store_true", dest="force", help="Force to do the modelling even if the file already exists (default = False)", metavar="{boolean}")
    parser.add_option("--opt", default=False, action="store_true", dest="optimization", help="Run modelling with optimization (default = False). Note: use it with care, this is not advisable as structures can take missconfomations", metavar="{boolean}")

    group = optparse.OptionGroup(parser, "FASTA options", "Provide a FASTA file as input.")
    group.add_option("-a", "--all", default=False, action="store_true", dest="all_models", help="Use ALL possible templates (default = False)", metavar="{boolean}")
    group.add_option("-d", "--dna", default=False, action="store_true", dest="dna", help="DNA modelling (if \"True\", remodels DNA using 3DNA; default = False)", metavar="{boolean}")
    group.add_option("-f","--full", default=False, action="store_true", dest="full_mode", help="Model all protein fragments with template that bind DNA (default = False)", metavar="{boolean}")
    group.add_option("--unbound_fragments", default=False, action="store_true", dest="full_fragments_mode", help="Model all protein fragments with template even unbound to DNA (default = False)", metavar="{boolean}")
    group.add_option("-e", "--twilight", default=False, action="store_true", dest="filter_twilight_zone_hits", help="Select template-hits over the twilight zone unless there are no models (default = False)", metavar="{boolean}")
    group.add_option("--restrictive", default=False, action="store_true", dest="filter_restrictive", help="Restrict modelling to only template-hits over the twilight zone when option --twilight is selected (default = False)", metavar="{boolean}")
    group.add_option("--unrestrictive", default=False, action="store_true", dest="no_restrict", help="Allow models with wrong alignment in the DNA-binding interface (default = False)", metavar="{boolean}")
    group.add_option("-s","--select", default=False, action="store_true", dest="select_mode", help="Select mode (if True, selects either monomers or dimers, whichever has more hits; default = False)", metavar="{boolean}")
    group.add_option("--dimers", default=False, action="store_true", dest="dimer_mode", help="If True forces to work only with dimers; default = False)", metavar="{boolean}")
    group.add_option("--monomers", default=False, action="store_true", dest="monomer_mode", help="If True forces to work only with monomers; default = False)", metavar="{boolean}")
    group.add_option("--renumerate", default=False, dest = 'renumerate' , action = 'store_true', help = 'Flag to renumber the sequences as in the original FastA input (default is False)')

    parser.add_option_group(group)

    (options, args) = parser.parse_args()

    if options.pdb_dir is None or options.fasta_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options


#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()
    verbose = options.verbose
    command = os.path.join(scripts_path,"model_protein.py")
    parameters = " "
    if options.force:                          parameters += " --force "
    if options.optimization :                  parameters += " --opt "
    if options.all_models :                    parameters += " --all "
    if options.dna :                           parameters += " --dna "
    if options.full_mode :                     parameters += " --full "
    if options.full_fragments_mode :           parameters += " --unbound_fragments "
    if options.filter_twilight_zone_hits :     parameters += " --twilight "
    if options.filter_restrictive :            parameters += " --restrictive "
    if options.no_restrict :                   parameters += " --unrestrictive "
    if options.select_mode :                   parameters += " --select  "
    if options.dimer_mode :                    parameters += " --dimers "
    if options.monomer_mode :                  parameters += " --monomers "
    if options.renumerate :                    parameters += " --renumerate "
    if options.n_total is not None:            parameters += " --n-total %d"%(options.n_total)
    if options.n_model > 1:                    parameters += " --n-model %d"%(options.n_model)

    # Create output subdirs #
    if not os.path.exists(options.dummy_dir): os.makedirs(options.dummy_dir)
    if not os.path.exists(options.output_dir): os.makedirs(os.path.abspath(options.output_dir))

    for tf in SEQ.FASTA_iterator(options.fasta_file):
        name=tf.get_identifier().split()[0].split("|")[0]
        label=name
        out=os.path.abspath(os.path.join(options.output_dir,name))
        input_fasta=os.path.join(options.dummy_dir,name+".fasta")
        if verbose: sys.stdout.write("\t-- model %s\n"%(label))
        fa=open(input_fasta,"w")
        fa.write(">%s\n%s\n"%(name,tf.get_sequence()))
        fa.close()       
        if options.parallel:
            if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
            else: cluster_queue=config.get("Cluster", "cluster_queue")
            functions.submit_command_to_queue("%s %s --pdb=%s -o %s -l %s --dummy=%s -i %s %s -v " % (os.path.join(config.get("Paths", "python_path"), "python"),command,options.pdb_dir,out,label,options.dummy_dir,input_fasta,parameters), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),options.dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
        else:
            os.system("%s %s --pdb=%s -o %s -l %s --dummy=%s -i %s %s -v " % (os.path.join(config.get("Paths", "python_path"), "python"),command,options.pdb_dir,out,label,options.dummy_dir,input_fasta,parameters))
        if not verbose: os.remove(input_fasta)

    
    if verbose: sys.stdout.write("Done\n")

