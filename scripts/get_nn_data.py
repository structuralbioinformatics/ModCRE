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
    parser = optparse.OptionParser("python -m cisbp_sql_motifs_file  --ps proteins_file [ -o rootname -d folder -v]")

    parser.add_option("-m", action="store", default=None, type="string", dest="sql_file", help="SQL file of motifs (from CIS-BP)", metavar="{file}")
    parser.add_option("--ps", action="store", default=None, type="string", dest="prot_file", help="Protein sequences (from CIS-BP)", metavar="{file}")
    parser.add_option("-o", action="store", type="string",default="CisBP",dest="root", help="CisBP file in SQL format with table of motifs (default CisBP)", metavar="{rootname}")
    parser.add_option("-d", action="store", type="string",default="sequences",dest="root", help="Folder to store the sequences of TFs (default sequences)", metavar="{foldername}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")

    (options, args) = parser.parse_args()

    if options.prot_file is None or options.sql_file is None:
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

    # Initialize #
    map_pdb_sp={}
    idmap=open(options.map_file,"r")
    for line in idmap:
        data_id=line.strip().split()
        if data_id[1]=="PDB":
           map_pdb_sp.setdefault(data_id[2].lower(),set()).add(data_id[0])
    pdb = {}
    pdball = set()
    pdb_dict_file = os.path.join(os.path.abspath(options.output_dir), "pdb.txt")
    # Skip if dict file already exists #
    if not os.path.exists(pdb_dict_file):
        # For each PDB file... #
        for pdb_file in os.listdir(os.path.abspath(options.pdb_dir)):
            # For each line... #
            for line in functions.parse_file(os.path.join(os.path.abspath(options.pdb_dir), pdb_file)):
                if line.startswith("ATOM") or line.startswith("HETATOM"): break
                if not line.startswith("DBREF"): continue
                m = re.search("^DBREF\s+(\S{4})\s+(\S).+UNP\s+(\S{6}).+", line)
                if m:
                    pdb.setdefault((m.group(1).lower(), m.group(2)), set())
                    pdb[(m.group(1).lower(), m.group(2))].add(m.group(3))
                    if "bundle" in pdb_file: 
                        if pdb_file.endswith("pdb"): pdb[pdb_file.rstrip(".pdb").lower(), m.group(2)].add(m.group(3))
                        if pdb_file.endswith("pdb1"): pdb[pdb_file.rstrip(".pdb1").lower(), m.group(2)].add(m.group(3))
                    if map_pdb_sp.has_key(m.group(1).lower()):
                        pdb[(m.group(1).lower(), m.group(2))].update(map_pdb_sp.get(m.group(1).lower()))
                        if "bundle" in pdb_file: 
                           if pdb_file.endswith("pdb"): pdb[pdb_file.rstrip(".pdb").lower(), m.group(2)].update(map_pdb_sp.get(m.group(1).lower()))
                           if pdb_file.endswith("pdb1"): pdb[pdb_file.rstrip(".pdb1").lower(), m.group(2)].update(map_pdb_sp.get(m.group(1).lower()))
        for pdb_name, pdb_chain in sorted(pdb):
            functions.write(pdb_dict_file, "%s;%s" % (",".join([pdb_name, pdb_chain]), ",".join(sorted(pdb[(pdb_name, pdb_chain)]))))
            pdball.add((pdb_name,pdb_chain))
    else:
        for line in functions.parse_file(pdb_dict_file):
            key, value = line.split(";")
            pdball.add(tuple(key.split(",")))
            pdb.setdefault(tuple(key.split(",")), set(value.split(",")))
    for pdb_file in os.listdir(os.path.abspath(options.pdb_dir)):
      if pdb_file[-3:] == "pdb" or  pdb_file[-4:] == "pdb1":
        pdb_name=pdb_file[0:4]
        if pdb_name not in set([key for key,chain in pdball]):
               pdball.add((pdb_name,None))
    print "PDB length",len([k for k in pdb.iterkeys()]) 
    print "PDB codes",[k for k in pdb.iterkeys()]
    print pdb
    for pdb_name, pdb_chain in sorted(pdb): 
        print "PDB SORTED",pdb_name, pdb_chain
    # Initialize #
    uniaccs = set()
    for key in pdb:
        for uniacc in pdb[key]:
            uniaccs.add(uniacc)

    print "Uniprot size", len(uniaccs)
    print "Uniprot codes", uniaccs

    cisbp_families={}
    print options.family_file
    for line in functions.parse_file(os.path.abspath(options.family_file)):
        m = re.search("\('(.+)', '(.+)', '.+', .+, .+\)", line)
        if m:           
            cisbp_families.setdefault(m.group(1).upper(),set()).add(m.group(2))
            print "Add family ",m.group(1).upper(),m.group(2)
    # Initialize #
    cisbp = {}
    cisbp_full = {}
    cisbp_dict_file = os.path.join(os.path.abspath(options.output_dir), "cisbp.txt")
    # Skip if dict file already exists #
    if not os.path.exists(cisbp_dict_file):
        # Initialize
        species = {}
        families = {}
        # For each line... #
        for line in functions.parse_file(os.path.abspath(options.tfs_file)):
            #print line
            m = re.search("\('(.+)', '(.+)', '.+', '.+', '.+', '(.+)', '[DIN]'\),*", line)
            #print m
            if m:
                species.setdefault(m.group(1), set()).add(m.group(3).replace("_", " ").upper())
                if options.cisbp_family:
                  families.setdefault(m.group(1).upper(), set()).add(m.group(2).upper())
                  print "Add TF on CISBP ",m.group(1).upper(),m.group(2).upper(),species[m.group(1)]
                else:
                  families.setdefault(m.group(1).upper(), set()).update(cisbp_families[m.group(2).upper()])
                  print "Add TF on CISBP ",m.group(1).upper(),m.group(2).upper(),cisbp_families[m.group(2).upper()],species[m.group(1)]
            else:
                print "SKIP TF on CISBP ",line
        # For each line... #
        n_pid=0
        missing_p={}
        missing_t={}
        for line in functions.parse_file(os.path.abspath(options.proteins_file)):
            m = re.search("\('(.+)', '(.+)', '.+', '.+', '(.+)'\),*", line)
            if m:
                print "PROTEIN FOR TF", m.group(2)
                n_pid = n_pid + 1
                if m.group(2) not in [x for x in species.iterkeys()]: 
                    missing_p.setdefault(m.group(1),set()).add( m.group(2))
                    missing_t.setdefault(m.group(2),set()).add( m.group(1))
                    continue
                print "  -- accept TF protein ID", m.group(1)
                sequence = m.group(3)
                if sequence.endswith("*"):
                    sequence = sequence[:-1]
                h = hashlib.new('md5')
                h.update(sequence)
                md5 = h.hexdigest() + sequence[:4] + sequence[-4:]
                print "  -- SPECIES", species[m.group(2)]
                for specie in species[m.group(2)]:
                    cisbp.setdefault((md5, specie), set()).add(m.group(2))
                    cisbp_full.setdefault(tuple([md5, specie]), (sequence,m.group(1)))
            else:
                print "SKIP PROTEIN ",line
        print "Total proteins",n_pid
        print "Total TFs",len([x for x in species.iterkeys()])
        print "Total TFs with known specie and protein sequence",len([x for x in cisbp.iterkeys()])
        print "Total missing proteins",len([x for x in missing_p.iterkeys()])
        print "Total missing TF",len([x for x in missing_t.iterkeys()])
        for mp,mt in missing_p.iteritems():
            print "MISSING PROTEIN",mp,mt
        for mt,mp in missing_t.iteritems():
            print "MISSING TF     ",mt,mp
        for md5, specie in sorted(cisbp):
            functions.write(cisbp_dict_file, "%s;%s" % (",".join([md5, specie]), ",".join(sorted(cisbp[(md5, specie)]))))
    else:
        for line in functions.parse_file(cisbp_dict_file):
            key, value = line.split(";")
            cisbp.setdefault(tuple(key.split(",")), set(value.split(",")))
    print "Total of proteins of TFs accepted",len([x for x in cisbp.iterkeys()])
    total_tfs=set()
    for k,v in cisbp.iteritems():
        mark_tf=False
        if len(v)>1:
            print "Sequence shared by TFs",k
            mark_tf=True
        for tf in v:
            total_tfs.add(tf)
            if mark_tf:
                print "  Shared sequence for ",tf
    print "Total of TFs",len(total_tfs)
    # Initialize #
    uniprot = {}
    uniprot_full = {}
    uniprot_dict_file = os.path.join(os.path.abspath(options.output_dir), "uniprot.txt")
    # Skip if dict file already exists #
    print "UNIPROT PARSING"
    if not os.path.exists(uniprot_dict_file):
        # Initialize #
        gz = False
        if options.uniprot_file.endswith(".gz"): gz = True
        # For FASTA sequence... #
        for header, sequence in functions.parse_fasta_file(os.path.abspath(options.uniprot_file), gz=gz):
            # Initialize #
            #print "  --Parse",header
            uniacc = None
            specie = None
            #m = re.search("\S+\|(\S+)\|.+OS=(.+) GN=.+ PE=.+ SV=.+", header)
            mm = re.search("\S+\|(\S+)\|.+OS=(.+) OX=.+ GN=.+ PE=.+ SV=.+", header)
            if mm:
                uniacc = mm.group(1)
                specie_list = mm.group(2).upper().split()
                specie=" ".join(specie_list[0:2])
            else:
              mn = re.search("\S+\|(\S+)\|.+OS=(.+) GN=.+ PE=.+ SV=.+", header)
              if mn:
                uniacc = mn.group(1)
                specie_list = mn.group(2).upper().split()
                specie=" ".join(specie_list[0:2])
              else:
                mx= re.search("\S+\|(\S+)\|.+OS=(.+) OX=.+ PE=.+ SV=.+", header)
                if mx:
                  uniacc = mx.group(1)
                  specie_list = mx.group(2).upper().split()
                  specie=" ".join(specie_list[0:2])
                else:
                  print "SKIP UNIPROT",header
            if uniacc is None or specie is None: continue
            if uniacc in uniaccs:
                if sequence.endswith("*"):
                    sequence = sequence[:-1]
                h = hashlib.new('md5')
                h.update(sequence)
                md5 = h.hexdigest() + sequence[:4] + sequence[-4:]
                uniprot.setdefault(uniacc, tuple([md5, specie]))
                uniprot_full.setdefault(tuple([md5, specie]), tuple([sequence,uniacc]))
                print "Add in UNIPROT",uniacc
        for uniacc in sorted(uniprot):
            functions.write(uniprot_dict_file, "%s;%s" % (uniacc, ",".join(uniprot[uniacc])))
    else:
        for line in functions.parse_file(uniprot_dict_file):
            key, value = line.split(";")
            uniprot.setdefault(key, tuple(value.split(",")))

    
    # Initialize #
    tfs = {}
    tfs_cisbp = set()
    tfs_dict_file = os.path.join(os.path.abspath(options.output_dir), "tfs_cisbp.txt")
    # Skip if dict file already exists #
    if not os.path.exists(tfs_dict_file):
        #Read again TF families if families dictionary is empty
        refill=False
        try:
         if len(families.viewitems()) < 1: refill=True
        except:
         refill=True
        if refill:
         families = {}
         # For each line... #
         for line in functions.parse_file(os.path.abspath(options.tfs_file)):
            #print line
            m = re.search("\('(.+)', '(.+)', '.+', '.+', '.+', '(.+)', '[DIN]'\),*", line)
            #print m
            if m:
                if options.cisbp_family:
                  families.setdefault(m.group(1).upper(), set()).add(m.group(2).upper())
                  print "Add TF on TF-CISBP ",m.group(1).upper(),m.group(2).upper()
                else:
                  families.setdefault(m.group(1).upper(), set()).update(cisbp_families[m.group(2).upper()])          
                  print "Add TF on TF-CISBP ",m.group(1).upper(),m.group(2).upper(),cisbp_families[m.group(2).upper()]
            else:
                print "SKIP TF on TF-CISBP second check",line
        functions.write(tfs_dict_file, "#pdb;chain;family")
        # For FASTA sequence... #
        #for pdb_name, pdb_chain in sorted(pdb):
        for pdb_name, pdb_chain in pdb.iterkeys():
            print "CHECK PDB Accessions",pdb_name,pdb_chain,pdb[(pdb_name, pdb_chain)]
            u_seq=None
            c_seq=None
            for uniacc in pdb[(pdb_name, pdb_chain)]:
                print "Check UniAcc",uniacc
                if uniacc not in uniprot: 
                   print "Not found in UNIPROT",uniacc
                   continue
                md5, specie = uniprot[uniacc]
                specie_list=specie.split()
                print "Species",uniprot[uniacc], specie_set
                u_seq=(uniprot_full[(md5,specie)])
                c_seq=None
                for cisbp_md5, cisbp_specie in cisbp.iterkeys():
                    cisbp_specie_list=cisbp_specie.split()
                    if not cisbp_full.has_key((cisbp_md5,cisbp_specie)):
                        print "CISBP data not found",cisbp_md5,cisbp_specie
                        continue
                    c_seq=(cisbp_full[(cisbp_md5,cisbp_specie)])
                    if md5 == cisbp_md5 or u_seq[0] in c_seq[0] or c_seq[0] in u_seq[0]:
                        print "Found md5",md5,cisbp_md5,u_seq[1],c_seq[1]
                        if md5!=cisbp_md5: 
                           print "Different MD5 but matching sequences", u_seq[0],c_seq[0]
                           print "  -- check specie matching sequences", cisbp_specie,specie
                           if  u_seq[0] in c_seq[0]:
                               print "  -- check code matching seq U<C  ",u_seq[1],c_seq[1] 
                           else:
                               print "  -- check code matching seq C<U  ",u_seq[1],c_seq[1] 
                        if (cisbp_specie_list[0]==specie_list[0] and cisbp_specie_list[1]==specie_list[1]):
                            print "Found specie",cisbp_specie,specie
                            pdb_family_set=set()
                            for tf_id in cisbp[(cisbp_md5,cisbp_specie)]:
                              if not tf_id.upper() in families: 
                                  print "Not family found for TF",tf_id 
                                  continue
                              for family_id in families[tf_id.upper()]:
                                  pdb_family_set.update(set(family_id.split(",")))
                            pdb_family=",".join([str(x) for x in pdb_family_set])
                            functions.write(tfs_dict_file, "%s;%s;%s" % (pdb_name, pdb_chain,pdb_family))
                            tfs_cisbp.add((pdb_name, pdb_chain,pdb_family))
    else:
       fo=open(tfs_dict_file,"r")
       for line in fo:
         if not line.startswith("#"):
           (pdb_name, pdb_chain, pdb_family)=line.strip().split(";")
           tfs_cisbp.add((pdb_name, pdb_chain, pdb_family))
    # Check with all TF definitios by TFinderSelect
    # Initialize #
    pdbinput=pdball
    tfs_large=tfs_cisbp
    tfs_large_dict_file = os.path.join(os.path.abspath(options.output_dir), "tfs.txt")
    n=0
    if not os.path.exists(tfs_large_dict_file): 
     functions.write(tfs_large_dict_file, "%s;%s;%s" %("pdb","chain","family"))
     while len(pdbinput) > 0 and n<10:
      sys.stderr.write("Step %d: Searching TF among %d PDB files\n"%(n,len(pdbinput)))
      tf_new=set()
      retry=set()
      (tf_new,retry)=TFinderSelect.TFinder(pdbinput,tfs_cisbp)
      pdbinput=retry
      tfs_large.update(tf_new)
      n=n+1

     if n>=10:
      sys.stderr.write("Reached more than 10 trial connections. Please check your internet provider\n")
      for pdb_name, pdb_chain, pdb_family in tfs_large:
        functions.write(tfs_large_dict_file, "%s;%s;%s" % (pdb_name, pdb_chain, pdb_family))
     else:
      sys.stderr.write("Total number of re-dials %d\n"%n)
      for pdb_name, pdb_chain, pdb_family in tfs_large:
        functions.write(tfs_large_dict_file, "%s;%s;%s" % (pdb_name, pdb_chain, pdb_family))

      
