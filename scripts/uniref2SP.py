import os, sys, re
import ConfigParser
import hashlib
import optparse

# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.dirname(__file__))

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Imports my functions #
import functions
import TFinderSelect

#-------------#
# Options     #
#-------------#

def parse_options():
    '''
    This function parses the command line arguments and returns an optparse object.
    '''

    parser = optparse.OptionParser("python uniref2SP.py -u uniref_file -o output_file -v]")

    parser.add_option("-o", "--output-file", default="SPROT", action="store", type="string", help="Output file (default =SPROT)",dest="uniprot_file")
    parser.add_option("-u", action="store", type="string", dest="uniref_file", help="UniProt file (i.e. uniprot_sprot.fasta + uniprot_trembl.fasta", metavar="{filename}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)", metavar="{boolean}")

    (options, args) = parser.parse_args()

    if  options.uniprot_file is None :
        parser.error("missing arguments: type option \"-h\" for help")

    return options

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()
    uniprot_file=options.uniprot_file
    uniref_file=options.uniref_file
 
    if not os.path.exists(uniprot_file):
        # Initialize #
        gz = False
        if uniref_file.endswith(".gz"): gz = True
        # For FASTA sequence... #
        for header, sequence in functions.parse_fasta_file(os.path.abspath(uniref_file), gz=gz):
            # Initialize #
            #print "  --Parse",header
            uniacc = None
            specie = None
            mm = re.search("(.+) Tax=(.+) TaxID=(.+) RepID=(.+)", header)
            if mm:
               text=mm.group(1).lstrip(">UniRef")
               data=text.split()
               acc=data[0].split("_")[1]
               name=" ".join(data[1:])
               os=mm.group(2)
               ox=mm.group(3) 
               entry=mm.group(4)
               new_header="tr|%s|%s %s OS=%s OX=%s PE= SV= "%(acc,entry,name,os,ox)
               print acc,entry,os
            else:
               print "Failed HEADER %s"%(header)
               new_header=header
            functions.write(uniprot_file,">%s\n%s\n"%(new_header,sequence))
                
