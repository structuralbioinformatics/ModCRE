import os, sys, re
import ConfigParser
import copy
import optparse
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
import functions, contacts, dssp, model_protein

# Imports jbonet's modules #
from SBI.external.blast import BlastResult as BR
from SBI.external.blast import BlastHit as BH
from SBI.structure import PDB
from SBI.structure.chain import Chain

#-------------#
# Functions   #
#-------------#

def filter_hit(query_alignment, hit_alignment, chain_id, pdb_obj, contacts_obj, dssp_obj):
    '''
    This function filters a hit if the query does not cover all of the core
    (i.e. region of the hit that contacts the DNA). 

    @return           = boolean
    '''

    # Get sequence/PDB correlation #
    sequence_to_crystal, crystal_to_sequence = model_protein.get_sequence_to_crystal_correlations(pdb_obj=pdb_obj, pdb_chain=chain_id)
    # Get sequence #
    sequence = pdb_obj.get_chain_by_id(chain_id).gapped_protein_sequence
    # Get core amino-acids #
    core_residues=set()
    # For each contact #
    for contact_obj in contacts_obj.get_contacts():
        if contact_obj._A_chain != chain_id: continue
        if not pdb_obj.chain_exists(contact_obj._A_chain): continue
        if not pdb_obj.get_chain_by_id(contact_obj._A_chain).residue_exists(str(contact_obj._A_residue_obj.number)): continue
        core_residues.add(str(contact_obj._A_residue_obj.number))
   # Get secondary structure #
    secondary_structure = ""
    for i in range(len(sequence)):
        if sequence[i] == "x":
            secondary_structure += "C"
        else:
            j = dssp_obj.get_secondary_structure(pdb_obj.get_chain_by_id(chain_id), sequence_to_crystal[i])
            if j == "E" or j == "H":
                 if sequence_to_crystal[i] in core_residues:
                     j = "*"
            secondary_structure += j
    secondary_structure = re.sub("[^EH*]", "C", secondary_structure)
    # Get start, end hit #
    start = sequence.index(hit_alignment.replace("-", ""))
    end = start + len(hit_alignment.replace("-", "")) - 1
    # Get query secondary structure #
    query_secondary_structure = ""
    positions = range(start, end + 1)
    for i in range(len(hit_alignment)):
        if hit_alignment[i] != "-":
            position = positions.pop(0)
            if query_alignment[i] != "-":
                query_secondary_structure += secondary_structure[position]
            else:
                query_secondary_structure += "-"
        else:
            query_secondary_structure += "-"
    # If contact residues are not conserved, skip #
    if len(re.findall("\*", query_secondary_structure)) != len(re.findall("\*", secondary_structure)):
        return True
    # If insertions/deletions in structured region holding contact residues, skip #
    if re.search("[\*|E]+\-+[\*|E]+", query_secondary_structure) or re.search("[\*|H]+\-+[\*|H]+", query_secondary_structure):
        return True

    return False


def get_homologs(options):

    input_file=options.input_file
    sys.stdout.write("\t- Homologs of %s\n"%input_file)

    if options.output_file is None:
       output_file=os.path.join(options.pbm_dir,"homologs",input_file.split("/")[-1]+"homologs.out")
    else:
       output_file=options.output_file
    
    if os.path.exists(output_file):
       os.remove(output_file)

    homologs_obj=Homologs(input_file,options.filter_hits)   

    # Select orthologs
    if options.specie is not None:
       orthologs=homologs_obj.get_orthologs(options.specie)
    else:
       orthologs=homologs_obj.get_homologs()

    #replace the set of homologs
    homologs_obj._homologs=orthologs

    # Filter by contact interface
    if options.filter_hits:
       if not os.path.exists(options.dummy_dir):
          os.makedirs(options.dummy_dir)
       filtered=homologs_obj.filter_hits_by_interface(options.pdb_dir,options.dummy_dir)
    else:
       filtered=homologs_obj.get_homologs()

    #replace the set of homologs
    homologs_obj._homologs=filtered


    # Write ouput file
    homologs_obj.write(output_file)

    

def parse_options():
    '''
    Create a file of homologs ensuring the ungapped sequence and DNA contacts.
    '''

    parser = optparse.OptionParser("homologs.py -i INPUT_FILE  --pbm=PBM_DIR --pdb=PDB_DIR  [--filter --dummy DUMMY_DIR --out OUTPUT_FILE --specie SPECIE -v]")

    parser.add_option("-i", action="store", type="string", dest="input_file", default=None, help="Input blast/hmm file", metavar="INPUT_FILE")
    parser.add_option("--filter", default=False, action="store_true", dest="filter_hits", help="Filter twilight zone hits and  homologs that do not cover 100% of interface (default = False)", metavar="{boolean}")
    parser.add_option("--pbm", action="store", type="string", dest="pbm_dir", help="PBM directory (i.e. output dir from pbm.py)", metavar="PBM_DIR")
    parser.add_option("--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (i.e. output dir from pdb.py)", metavar="PDB_DIR")
    parser.add_option("--dummy", action="store", type="string", dest="dummy_dir", default="/tmp", help="Dummy directory (default is /tmp)", metavar="DUMMY_DIR")
    parser.add_option("-o","--out", action="store", type="string", dest="output_file", default=None, help="Output file (default uses PBM/homologs directory", metavar="OUTPUT_FILE")
    parser.add_option("--specie", action="store", type="string", dest="specie",  default=None, help="Specie to obtain specific orthologs (i.e. taxon/code/common_name as 9606/HUMAN/'Homo sapiens')", metavar="SPECIE")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")

    (options, args) = parser.parse_args()

    if options.input_file is None  or options.pdb_dir is None or options.pbm_dir is None :
        parser.error("missing arguments: type option \"-h\" for help")

    return options




#-------------#
# Classes     #
#-------------#

class Species(object):
    '''
    This class defines a Species file.
    '''

    def __init__(self, species_file):
        self._file = species_file
        self._taxons = None

        self._parse_file()

    def _get_file(self):
        return self._file

    def _parse_file(self):
        if os.path.exists(self._get_file()):
            read = False
            self._taxons = {}

            fd = open(self._get_file(), "rt")
            for line in fd:
                if read:
                    m = re.search("^(\S+)\s+([ABEV])\s+(\d+)", line.strip())
                    if m:
                        taxon = m.group(3)
                        self._taxons.setdefault(taxon, {'code' : m.group(1), 'kingdom' : m.group(2), 'officialName' : None, 'commonName' : None, 'synonym' : None})
                    m = re.search("([NCS])=(.+)$", line.strip())
                    if m:
                        if m.group(1) == "N":
                            self._taxons[taxon]['officialName'] = m.group(2)
                        if m.group(1) == "C":
                            self._taxons[taxon]['commonName'] = m.group(2)
                        if m.group(1) == "S":
                            self._taxons[taxon]['synonym'] = m.group(2)
                if "_____ _ _______ _____________________________________________________________" in line:
                    read = True
            fd.close()
        else:
            raise ValueError("Could not open species file %s" % self._get_file())

    def get_taxons(self, code=None, kingdom=None, name=None):
        taxons = []

        if code != None:
            for taxon in self._taxons:
                if self._taxons[taxon]['code'] == code:
                    taxons.append(taxon)

        if kingdom != None:
            for taxon in self._taxons:
                if self._taxons[taxon]['kingdom'] == kingdom:
                    taxons.append(taxon)

        if name != None:
            for taxon in self._taxons:
                if self._taxons[taxon]['officialName'] == name or self._taxons[taxon]['commonName'] == name or self._taxons[taxon]['synonym'] == name:
                    taxons.append(taxon)

        return taxons



class Homologs(object):
    '''
    This class defines a Homologs object.
    '''

    def __init__(self, homologs_file=None, use_filter=False ):
        self._file = homologs_file
        self._homologs = BR()
        self._filter   = use_filter
        if homologs_file != None:
            self._parse_file()

    def _get_file(self):
        return self._file

    def _parse_file(self):
        hfile=self._get_file()
        if os.path.exists(hfile):
            if  hfile.endswith(".gz"):
              inp= gzip.open(hfile,"rt")
            else:
              inp= open(hfile,"rt")
            for hit in inp:
                qpos=[]
                hpos=[]
                if len(hit.strip().split("\t"))<10:
                   sys.stdout.write("\t- Failure: no homologs found\n")
                   exit(0)
                self._homologs._query        = hit.strip().split("\t")[0]
                self._homologs._query_length = int(hit.strip().split("\t")[1])
                hit_name                     = hit.strip().split("\t")[2]
                dummy_length                 = int(hit.strip().split("\t")[3])
                identities                   = int(hit.strip().split("\t")[4])
                positives                    = int(hit.strip().split("\t")[5])
                gaps                         = int(hit.strip().split("\t")[6])
                e_value                      = hit.strip().split("\t")[7]
                query_seq                    = hit.strip().split("\t")[8]
                hit_seq                      = hit.strip().split("\t")[9]
                align_length                 = len(query_seq)
                Positions                    = hit.strip().split("\t")[10]
                segments                     = Positions.split(";")
                for segment in segments:
                    seqpos_a         = segment.split(",")[0]
                    seqpos_b         = segment.split(",")[1]
                    qpos_a           = int(seqpos_a.split(":")[0])
                    qpos_b           = int(seqpos_b.split(":")[0])
                    hpos_a           = int(seqpos_a.split(":")[1])
                    hpos_b           = int(seqpos_b.split(":")[1])
                    qpos.append(qpos_a)
                    qpos.append(qpos_b)
                    hpos.append(hpos_a)
                    hpos.append(hpos_b)
                hit_length  = hpos[-1] - hpos[0] + 1
                OutputHit   = BH(name = hit_name, length = hit_length,  iteration  = 1,  e_value = e_value,  align_length = align_length,
                                 identities = identities,  positives = positives, gaps = gaps, qseq = query_seq, hseq = hit_seq,
                                 qpos = qpos[0],  hpos = hpos[0], score_seq = query_seq)
                self._homologs._hits.append(OutputHit)
                


    def add_homolog(self, blast_hit=None):
        if blast_hit is not None:
          self._homologs._hits.append(blast_hit)

    def get_homologs(self):
        return self._homologs

    def filter_hits_by_interface(self,pdb_dir,dummy_dir="/tmp"):
        # Initialize #
        unfiltered_hits = BR()
        # Twilight zone #
        tz_parameter = 0
        tz_type = None
        if self._filter:
            tz_parameter = int(config.get("Parameters", "twilight_zone_parameter"))
            tz_type = config.get("Parameters", "twilight_zone_type")
        # get structure of main sequence
        pdb_name=self._homologs.query.split("_")[0]
        chain_id=self._homologs.query.split("_")[-1]
        pdb_file = os.path.join(pdb_dir, "clean", pdb_name + ".pdb")
        if os.path.exists(pdb_file):
           try:
              pdb_obj  = PDB(pdb_file)
           except:
              return self._homologs
        else:
           #return self._homologs.get_hits(tz_parameter=tz_parameter, tz_type=tz_type)
           return self._homologs
        # get contacts
        contact_file=os.path.join(pdb_dir, "contacts",pdb_name + ".txt")
        if os.path.exists(contact_file):
           contacts_obj=contacts.Contacts(contact_file)
        else:
           #return self._homologs.get_hits(tz_parameter=tz_parameter, tz_type=tz_type)
           return self._homologs
        # Get sequence/PDB correlation #
        sequence_to_crystal, crystal_to_sequence = model_protein.get_sequence_to_crystal_correlations(pdb_obj=pdb_obj, pdb_chain=chain_id)
        # Get main sequence, its fasta file and name #
        template_file                = os.path.join(pdb_dir, "split",self._homologs._query + ".fasta")
        try:
          template_sequence            = pdb_obj.get_chain_by_id(chain_id).gapped_protein_sequence
        except:
          sys.stdout.write("Failed to get homologs with %s chain %s\n"%(pdb_name,chain_id))
          exit(0)
        unfiltered_hits._query       = self._homologs._query
        unfiltered_hits._query_length= len(template_sequence)
        # Get core amino-acids of the interface#
        core_residues=set()
        # For each contact #
        for contact_obj in contacts_obj.get_contacts():
            if contact_obj._A_chain != chain_id: continue
            if not pdb_obj.chain_exists(contact_obj._A_chain): continue
            if not pdb_obj.get_chain_by_id(contact_obj._A_chain).residue_exists(str(contact_obj._A_residue_obj.number)): continue
            core_residues.add(str(contact_obj._A_residue_obj.number))
        # Get template secondary structure #
        dssp_file= os.path.join(pdb_dir, "dssp", pdb_name + ".txt")
        if os.path.exists(dssp_file): 
           dssp_obj = dssp.DSSP(dssp_file)
        else:
           dssp_obj = dssp.get_dssp_obj(pdb_file,dummy_dir=dummy_dir)
        # Combine template secondary structure with interface definition#
        secondary_structure = ""
        for i in range(len(template_sequence)):
            if template_sequence[i] == "x":
                secondary_structure += "C"
            else:
                j = dssp_obj.get_secondary_structure(chain_id, sequence_to_crystal[i])
                if j is None: continue
                if j == "E" or j == "H":
                     if sequence_to_crystal[i] in core_residues:
                         j = "*"
                secondary_structure += j
        secondary_structure = re.sub("[^EH*]", "C", secondary_structure)
        for hit_obj in self._homologs.get_hits(tz_parameter=tz_parameter, tz_type=tz_type):
           # Skip if already done hit #
           #done = False
           #for unfiltered_hit in unfiltered_hits._hits:
           #  if hit_obj.sequenceID == unfiltered_hit.sequenceID: done = True
           #if done: continue
           # Refine alignment #
           m = re.search("^(tr|sp)\|(\S+)\|\S+\_(\S+)", hit_obj.sequenceID)
           if m:
               hit_name = m.group(2)
           else:
               hit_name=hit_obj.sequenceID.strip().split()[0].split("|")[0]
           hit_file=os.path.join(dummy_dir,hit_name+"_"+self._homologs.query+".fa")
           if not os.path.exists(hit_file):
                    functions.write(hit_file, ">%s\n%s" % (hit_obj.sequenceID, hit_obj.ungapped_hit_seq))
           else:
                    os.remove(hit_file)
                    functions.write(hit_file, ">%s\n%s" % (hit_obj.sequenceID, hit_obj.ungapped_hit_seq))
           template_alignment, hit_alignment, template_start, template_end, hit_start, hit_end, identity, similarity = model_protein.exec_matcher(template_file, hit_file)
           # Get hit secondary structure #
           hit_secondary_structure = ""
           template_positions = range(template_start - 1, template_end)
           for i in range(len(template_alignment)):
               if template_alignment[i] != "-":
                   j = template_positions.pop(0)
                   if hit_alignment[i] != "-" and j<len(secondary_structure):
                       hit_secondary_structure += secondary_structure[j]
                   else:
                       hit_secondary_structure += "-"
               else:
                   hit_secondary_structure += "-"
           # Skip if interface residues are not conserved #
           if len(re.findall("\*", hit_secondary_structure)) != len(re.findall("\*", secondary_structure)): continue
           # Skip if insertions/deletions in structured region holding interface residues #
           if re.search("[\*|E]+\-+[\*|E]+", hit_secondary_structure) or re.search("[\*|H]+\-+[\*|H]+", hit_secondary_structure): continue
           # Adjust query/hit alignments to template
           qpos=1
           Positions= hit_obj.formatPositions()
           segments = Positions.split(";")
           seqpos   = segments[0].split(",")[0]
           hpos     = int(seqpos.split(":")[1])
           hit_alignment      = "-" * (template_start - 1) + hit_alignment + "-" * (len(template_sequence) - template_end)
           template_alignment = template_sequence[:(template_start - 1)] + template_alignment + template_sequence[template_end:]
           # Add unfiltered hit #
           UnfilteredHit= BH(name = hit_obj.sequenceID, length = hit_obj.length,  iteration  = 1,  e_value = hit_obj.e_value,  align_length = len(hit_alignment),
                          identities = hit_obj.identities,  positives = hit_obj.positives, gaps = hit_obj.gaps, qseq = template_alignment, hseq = hit_alignment,
                          qpos = qpos,  hpos = hpos, score_seq = template_alignment)
           unfiltered_hits._hits.append(UnfilteredHit)
           os.remove(hit_file)

        return unfiltered_hits

    def get_orthologs(self, specie=None):
        # Initialize #
        orthologs               = BR()
        orthologs._query        = self._homologs.query
        orthologs._query_length = self._homologs.query_length
        orthologs._hits         = []
        species_file =os.path.join(config.get("Paths","files_path"),config.get("Paths","species"))
        species_obj=Species(species_file)
        # Twilight zone #
        tz_parameter = 0
        tz_type = None
        if self._filter:
            tz_parameter = int(config.get("Parameters", "twilight_zone_parameter"))
            tz_type = config.get("Parameters", "twilight_zone_type")
        for hit_obj in self._homologs.get_hits(tz_parameter=tz_parameter, tz_type=tz_type):
           if specie is None:
              orthologs._hits.append(hit_obj)
              continue
           try:
               m = re.search("^(tr|sp|TR|SP)\|(\S+)\|\S+\_(\S+)", hit_obj.sequenceID)
               if m:
                   hit_name = m.group(2)
                   hit_specie_code= m.group(3)
               else:
                   hit_name=hit_obj.sequenceID.strip().split()[0].split("|")[0]
                   if len(hit_name)<3 :
                     if len(hit_obj.sequenceID.strip().split()[0].split("|"))>2:
                      hit_name=hit_obj.sequenceID.strip().split()[0].split("|")[2]
                     elif len(hit_obj.sequenceID.strip().split()[0].split("|"))>1:
                      hit_name=hit_obj.sequenceID.strip().split()[0].split("|")[1]
                     else:
                      hit_name=hit_obj.sequenceID.strip().split()[0].split("|")[0]
                   if len(hit_name.split("_")) >0:
                      hit_specie_code=hit_name.split("_")[1]
               m = re.search("^(tr|sp|TR|SP)\|(\S+)\|\S+\_(\S+)(.+)", hit_obj.sequenceID)
               taxons=[]  
               if m:
                  species_code = m.group(3)
                  taxons.extend(species_obj.get_taxons(species_code))
                  mm = re.search("(\S+) (\S+)",m.group(4).split("OS=")[1])
                  if mm:
                    organism = mm.group(1) + " " + mm.group(2)
                    taxons.extend(species_obj.get_taxons(organism))
               if specie is None:
                  orthologs._hits.append(hit_obj)
               else:
                  if specie in taxons or specie == organism or specie.upper() == species_code.upper() or specie.upper() in organism.upper():
                    orthologs._hits.append(hit_obj)
           except:
                  sys.stderr.write("WARNING HOMOLOGS: skip %s\n"%hit_obj.sequenceID)

        return orthologs

    def write(self, output_file):
        # Twilight zone #
        tz_parameter = 0
        tz_type = None
        if self._filter:
            tz_parameter = int(config.get("Parameters", "twilight_zone_parameter"))
            tz_type = config.get("Parameters", "twilight_zone_type")
        functions.write(output_file, self._homologs.str_compacted_blast(tz_parameter=tz_parameter, tz_type=tz_type))

#-------------#
# Main        #
#-------------#


if __name__ == "__main__":

    # Arguments & Options #
    options       = parse_options()
    homologs_obj  = get_homologs(options)

