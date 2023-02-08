import os, sys, re
import ConfigParser
import copy
import json
import numpy
import optparse

# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.dirname(__file__))

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Imports my functions #
import functions

# Import jbonet's module #
from SBI.data import aminoacids3to1, aminoacids_polarity_boolean

#-------------#
# Class       #
#-------------#

class Potentials(object):
    """
    This class defines a {Potentials} object.
    
    """

    def __init__(self, file_name=None, select_potential=None):
        self._file = file_name
        self._pmf_3d = None
        self._pmf_3dc = None
        self._pmf_s3dc = None
        self._pmf_local = None
        self._pmf_pair = None
        self._pmf_s3dc_dd = None
        self._pmf_s3dc_di = None
        self._distances = {}
        self._potential = select_potential

        if self._file is not None:
            self._parse_file()

    def _parse_file(self):
        for line in functions.parse_file(self._file):
            potential, json_obj = line.strip("\n").split("\t")
            if self._potential is None or self._potential=="all" or self._potential == potential[4:]:
                if potential == "pmf_3d" :
                    self._pmf_3d = json.loads(json_obj)
                if potential == "pmf_3dc"  :
                    self._pmf_3dc = json.loads(json_obj)
                if potential == "pmf_s3dc"  :
                    self._pmf_s3dc = json.loads(json_obj)
                if potential == "pmf_s3dc_dd"  :
                    self._pmf_s3dc_dd = json.loads(json_obj)
                if potential == "pmf_s3dc_di"  :
                    self._pmf_s3dc_di = json.loads(json_obj)
                if potential == "pmf_local"  :
                    self._pmf_local = json.loads(json_obj)
                if potential == "pmf_pair"  :
                    self._pmf_pair = json.loads(json_obj)
            if potential == "distances":
                distances = json.loads(json_obj)
                for position in range(len(distances)):
                    self._distances[distances[position]] = position

    def _get_distance(self, distance):
        for d in self._distances:
            if distance < d:
                return d

        return None

    def _get_distance_position(self, distance):
        if distance in self._distances:
            return self._distances[distance]

        return None

    def get_score(self, potential, key=None, distance=None):
        
        distance = numpy.floor(self._get_distance(distance))
        position = self._get_distance_position(distance)
        try:
          if position is not None:
            if potential == "3d":
                return self._pmf_3d[position]
            if potential == "3dc":
                if key in self._pmf_3dc:
                    return self._pmf_3dc[key][position]
            if potential == "s3dc":
                if key in self._pmf_s3dc:
                    return self._pmf_s3dc[key][position]
            if potential == "s3dc_dd":
                if key in self._pmf_s3dc_dd:
                    return self._pmf_s3dc_dd[key][position]
            if potential == "pair":
                if key in self._pmf_pair:
                    return self._pmf_pair[key][position]
          if potential == "local":
            if key in self._pmf_local:
                return self._pmf_local[key][0]
          if potential == "s3dc_di":
            if key in self._pmf_s3dc_di:
                return self._pmf_s3dc_di[key][0]
        except:
          return None

    def write(self, file_name):
        functions.write(file_name, "pmf_3d\t%s" % json.dumps(self._pmf_3d, separators=(",", ":")))
        functions.write(file_name, "pmf_3dc\t%s" % json.dumps(self._pmf_3dc, separators=(",", ":")))
        functions.write(file_name, "pmf_s3dc\t%s" % json.dumps(self._pmf_s3dc, separators=(",", ":")))
        functions.write(file_name, "pmf_s3dc_dd\t%s" % json.dumps(self._pmf_s3dc_dd, separators=(",", ":")))
        functions.write(file_name, "pmf_s3dc_di\t%s" % json.dumps(self._pmf_s3dc_di, separators=(",", ":")))
        functions.write(file_name, "pmf_local\t%s" % json.dumps(self._pmf_local, separators=(",", ":")))
        functions.write(file_name, "pmf_pair\t%s" % json.dumps(self._pmf_pair, separators=(",", ":")))
        functions.write(file_name, "distances\t%s" % json.dumps(self._distances, separators=(",", ":")))

#-------------#
# Options     #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python spotentials.py -i input_file [-a --dummy=dummy_dir -o output_file -s -v -z -b]")

    parser.add_option("-a", default=False, action="store_true", dest="approach", help="Approach missing contacts with evolutionary transitions from BLOSUM62 (i.e. Taylor's approach; default = False)", metavar="{boolean}")
    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (full path to triads files to be used to derive the potentials)", metavar="{filename}")
    parser.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="{filename}")
    parser.add_option("-s", default=False, action="store_true", dest="smooth", help="Smooth potentials (default = False)", metavar="{boolean}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)", metavar="{boolean}")
    parser.add_option("-z", default=False, action="store_true", dest="zscores", help="Calculate scanning-Zscores (default = False)", metavar="{boolean}")
    parser.add_option("-b", "--bins", default=False, action="store_true",  dest="computation", help="Computate the potentials: by bins (if selected) or accumulative (default).", metavar="{boolean}")

    (options, args) = parser.parse_args()

    if options.input_file is None:
        parser.error('missing arguments: type option \'-h\' for help')

    return options

def get_statistical_potentials(file_name, approach=False, smooth=False, zscores=False, computation=False, dummy_dir="/tmp"):
    """
    This functions derives statistical potentials from a list of triads
    files.

    @input:
    file_name {filename} full path to triads files

    @return:

    """

    # Initialize #
    if computation:
       bin_distance = float(config.get("Parameters", "bin_distance_bins"))
       smooth_sample= int(config.get("Parameters", "smooth_potential_bins"))
    else:
       bin_distance = float(config.get("Parameters", "bin_distance_accu"))
       smooth_sample= int(config.get("Parameters", "smooth_potential_accu"))
       
    max_contact_distance = float(config.get("Parameters", "max_contact_distance"))
    distances = list(numpy.arange(0, max_contact_distance + bin_distance, bin_distance)) # +bin_distance because of how range works
    kB = 8.31441e-3 # Boltzmann constant
    T = 298.0       # standard temperature

    # Get frequencies #
    
    f_dab, f_a_dab, f_a_dab_oa, f_a_b_dab, f_dab_oa_ob, f_a_b_dab_oa_ob = get_frequencies(file_name, computation=computation, approach=approach)
    # Calculate potentials #
    pmf_3d = [None] * len(distances)
    pmf_3dc = {}
    pmf_s3dc = {}
    pmf_s3dc_dd = {}
    pmf_s3dc_di = {}
    pmf_local = {}
    pmf_pair = {}
    # PMF 3D #
    for dr in range(len(distances)):
        if f_dab[dr] > 0:
            pmf_3d[dr] = kB * T * numpy.log(f_dab[dr] / float(sum(f_dab)))
    # PMF 3DC #
    for oa_ob in f_dab_oa_ob:
        pmf_3dc[oa_ob] = [None] * len(distances)
        for dr in range(len(distances)):
            if f_dab_oa_ob[oa_ob][dr] > 0:
                pmf_3dc[oa_ob][dr] = (kB * T * numpy.log(f_dab_oa_ob[oa_ob][dr] / float(sum(f_dab_oa_ob[oa_ob])))) - pmf_3d[dr]
    # PMF S3DC #
    for a_b_oa_ob in f_a_b_dab_oa_ob:
        a_oa, b_ob = a_b_oa_ob.split(";")
        a_list = a_oa.split("-")
        b_list = b_ob.split("-")
        a = a_list.pop(0)
        b = b_list.pop(0)
        oa = "-".join(a_list)
        ob = "-".join(b_list)
        a_b = a + ";" + b
        oa_ob = oa + ";" + ob
        pmf_s3dc[a_b_oa_ob] = [None] * len(distances)
        for dr in range(len(distances)):
            if f_a_b_dab_oa_ob[a_b_oa_ob][dr] > 0:
                pmf_s3dc[a_b_oa_ob][dr] = pmf_3d[dr] + pmf_3dc[oa_ob][dr] - (kB * T * numpy.log(f_a_b_dab_oa_ob[a_b_oa_ob][dr] / float(sum(f_a_b_dab_oa_ob[a_b_oa_ob]))))
    # PMF S3DC dd #
    for a_b_oa_ob in f_a_b_dab_oa_ob:
        pmf_s3dc_dd[a_b_oa_ob] = [None] * len(distances)
        for dr in range(len(distances)):
            if f_a_b_dab_oa_ob[a_b_oa_ob][dr] > 0:
                pmf_s3dc_dd[a_b_oa_ob][dr] = - (kB * T * numpy.log(f_a_b_dab_oa_ob[a_b_oa_ob][dr] / float(f_dab[dr])))
    # PMF S3DC di #
    for a_b_oa_ob in f_a_b_dab_oa_ob:
        pmf_s3dc_di[a_b_oa_ob] = [None]
        if sum(f_a_b_dab_oa_ob[a_b_oa_ob]) > 0:
            pmf_s3dc_di[a_b_oa_ob] = [kB * T * numpy.log(sum(f_a_b_dab_oa_ob[a_b_oa_ob]) / float(sum(f_dab)))]
    # PMF local #
    for a_oa in f_a_dab_oa:
        pmf_local[a_oa] = [None]
        if sum(f_a_dab_oa[a_oa]) > 0:
            a_list = a_oa.split("-")
            a = a_list.pop(0)
            pmf_local[a_oa] = [(kB * T * numpy.log(sum(f_a_dab_oa[a_oa]) / float(sum(f_dab)))) - (kB * T * numpy.log(sum(f_a_dab[a]) / float(sum(f_dab))))]
    # PMF pair #
    for a_b in f_a_b_dab:
        pmf_pair[a_b] = [None] * len(distances)
        for dr in range(len(distances)):
            if f_a_b_dab[a_b][dr] > 0:
                pmf_pair[a_b][dr] = pmf_3d[dr] - (kB * T * numpy.log(f_a_b_dab[a_b][dr] / float(sum(f_a_b_dab[a_b]))))

    # Approach potentials #
    if approach:
        # PMF S3DC #
        approached = []
        for a_b_oa_ob in pmf_s3dc:
            approached.append([a_b_oa_ob, approach_pmf(a_b_oa_ob, pmf_s3dc, "-")])
        for a_b_oa_ob, approached_pmf in approached:
            pmf_s3dc[a_b_oa_ob] = approached_pmf
        # PMF S3DC dd #
        approached = []
        for a_b_oa_ob in pmf_s3dc_dd:
            approached.append([a_b_oa_ob, approach_pmf(a_b_oa_ob, pmf_s3dc_dd, "-")])
        for a_b_oa_ob, approached_pmf in approached:
            pmf_s3dc_dd[a_b_oa_ob] = approached_pmf
        # PMF S3DC di #
        approached = []
        for a_b_oa_ob in pmf_s3dc_di:
            approached.append([a_b_oa_ob, approach_pmf(a_b_oa_ob, pmf_s3dc_di, "+")])
        for a_b_oa_ob, approached_pmf in approached:
            pmf_s3dc_di[a_b_oa_ob] = approached_pmf
        # PMF pair #
        approached = []
        for a_b in pmf_pair:
            approached.append([a_b, approach_pmf(a_b, pmf_pair, "-")])
        for a_b, approached_pmf in approached:
            pmf_pair[a_b] = approached_pmf

    # Calculate Z-scores #
    if zscores:
        pmf_3dc = calculate_zscores(pmf_3dc, "pmf_3dc", bin_distance)
        pmf_s3dc = calculate_zscores(pmf_s3dc, "pmf_s3dc", bin_distance)
        pmf_s3dc_dd = calculate_zscores(pmf_s3dc_dd, "pmf_s3dc_dd", bin_distance)
        pmf_s3dc_di = calculate_zscores(pmf_s3dc_di, "pmf_s3dc_di", bin_distance)
        pmf_local = calculate_zscores(pmf_local, "pmf_local", bin_distance)
        pmf_pair = calculate_zscores(pmf_pair, "pmf_pair", bin_distance)

    # Smooth potentials #
    if smooth:
        # PMF 3D #
        for dr in range(len(distances)):
            pmf_3d[dr] = smooth_bin(pmf_3d, dr,  smooth_sample)
        # PMF 3DC #
        for envpair in pmf_3dc:
            for dr in range(len(distances)):
                pmf_3dc[envpair][dr] = smooth_bin(pmf_3dc[envpair], dr, smooth_sample )
        # PMF S3DC #
        for pair in pmf_s3dc:
            for dr in range(len(distances)):
                pmf_s3dc[pair][dr] = smooth_bin(pmf_s3dc[pair], dr, smooth_sample )
        # PMF S3DC dd #
        for pair in pmf_s3dc_dd:
            for dr in range(len(distances)):
                pmf_s3dc_dd[pair][dr] = smooth_bin(pmf_s3dc_dd[pair], dr, smooth_sample ) 
        # PMF pair #
        for pair in pmf_pair:
            for dr in range(len(distances)):
                pmf_pair[pair][dr] = smooth_bin(pmf_pair[pair], dr, smooth_sample )


    return pmf_3d, pmf_3dc, pmf_s3dc, pmf_s3dc_dd, pmf_s3dc_di, pmf_local, pmf_pair, distances

def get_frequencies(file_name, computation=False, approach=False):
    """
    This function extracts the different frequencies to calculate the
    statistical potentials: i.e. "f_dab", "f_a_dab", "f_a_dab_oa", "f_a_b_dab",
    "f_dab_oa_ob" and "f_a_b_dab_oa_ob".

    @input:
    file_name {filename}

    @return:
    f_dab {list}
    f_a_dab {dict}
    f_a_dab_oa {dict}
    f_a_b_dab {dict}
    f_dab_oa_ob {dict}
    f_a_b_dab_oa_ob {dict}

    """

    # Initialize #
    if computation:
       bin_distance = float(config.get("Parameters", "bin_distance_bins"))
    else:
       bin_distance = float(config.get("Parameters", "bin_distance_accu"))
    max_contact_distance = float(config.get("Parameters", "max_contact_distance"))
    distances = list(numpy.arange(0, max_contact_distance + bin_distance, bin_distance)) # +bin_distance because of how range works
    aminoacids = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
    path = file_name.split("/")

    # The following arguments arguments are defined as in the papers:
    # "a" and "b" are the residues "a" and "b"
    # "oa" and "ob" are the environments of residues "a" and "b"
    # "a_oa" and "b_ob" are the residues "a" and "b" in their respective environments "oa" and "ob"
    # "dab" is the distance between residues "a" and "b"
    # "a_b_oa_ob" is a contact between residue "a_oa" and residue "b_ob"
    f_dab = [0] * len(distances) # contacts at a max. distance "dab"
    f_a_dab = {}                 # contacts at a max. distance "dab" involving residue "a"
    f_a_dab_oa = {}              # contacts at a max. distance "dab" involving residue "a" in environment "oa"
    f_a_b_dab = {}               # contacts at a max. distance "dab" involving residues "a" and "b"
    f_dab_oa_ob = {}             # contacts at a max. distance "dab" involving environments "oa" and "ob"
    f_a_b_dab_oa_ob = {}         # contacts at a max. distance "dab" involving residues "a" and "b" in their respective environments "oa" and "ob"

    # For each line... #
    # sys.stderr.write("\t\t--Frequencies: Read File %s\n"%(file_name))
    for line in functions.parse_file(file_name):
        # Skip line if not a file name #
        if line.startswith("#"): continue
        # Try... #
        try:
            # For each line... #
            # sys.stderr.write("\t\t--Frequencies: add triad %s\n"%(os.path.join("/".join(path[:-3]), line)))
            for line in functions.parse_file(os.path.join("/".join(path[:-3]), line)):
                # Skip line if not a triad #
                if line.startswith("#"): continue
                # Get triad #
                line = line.split(";")
                dab = line[2]
                # Adjust distance #
                dab = adjust_distance_to_bin(float(dab),bin_distance)
                # Skip if distance is too large #
                if dab > max_contact_distance: continue
                a_list = line[0].split("-")
                b_list = line[1].split("-")
                a = a_list.pop(0)
                b = b_list.pop(0)
                oa = "-".join(a_list)
                ob = "-".join(b_list)
                a_oa = a + "-" + oa
                b_ob = b + "-" + ob
                a_b_oa_ob = a_oa + ";" + b_ob
                a_b = a + ";" + b
                oa_ob = oa + ";" + ob
                if approach:
                    # Initialize #
                    for aa in aminoacids:
                        hydrophobicity = "N"
                        if aminoacids_polarity_boolean[aminoacids3to1[aa]]:
                            hydrophobicity = "P"
                        oaa = hydrophobicity + "-" + a_list[1] + "-" + a_list[2]
                        aa_oaa = aa + "-" + oaa
                        aa_b = aa + ";" + b
                        oaa_ob = oaa + ";" + ob
                        aa_b_oaa_ob = aa_oaa + ";" + b_ob
                        # Initialize frequencies #
                        if aa not in f_a_dab:
                            f_a_dab[aa] = [0] * len(distances)
                        # For local only do amino acids #
                        if aa_oaa not in f_a_dab_oa:
                            f_a_dab_oa[aa_oaa] = [0] * len(distances)
                        if aa_b not in f_a_b_dab:
                            f_a_b_dab[aa_b] = [0] * len(distances)
                        if oaa_ob not in f_dab_oa_ob:
                            f_dab_oa_ob[oaa_ob] = [0] * len(distances)
                        if aa_b_oaa_ob not in f_a_b_dab_oa_ob:
                            f_a_b_dab_oa_ob[aa_b_oaa_ob] = [0] * len(distances)
                # Initialize frequencies #
                if a not in f_a_dab:
                    f_a_dab[a] = [0] * len(distances)
                if b not in f_a_dab:
                    f_a_dab[b] = [0] * len(distances)
                # For local only do amino acids #
                if a_oa not in f_a_dab_oa:
                    f_a_dab_oa[a_oa] = [0] * len(distances)
                if a_b not in f_a_b_dab:
                    f_a_b_dab[a_b] = [0] * len(distances)
                if oa_ob not in f_dab_oa_ob:
                    f_dab_oa_ob[oa_ob] = [0] * len(distances)
                if a_b_oa_ob not in f_a_b_dab_oa_ob:
                    f_a_b_dab_oa_ob[a_b_oa_ob] = [0] * len(distances)
                # Update #
                for dr in range(1,len(distances)):
                    if computation :
                        if distances[dr-1] < dab <= distances[dr]:
                            f_dab[dr] += 1
                            f_a_dab[a][dr] += 1
                            f_a_dab[b][dr] += 1
                            f_a_b_dab[a_b][dr] += 1
                            # For local only do amino acids #
                            f_a_dab_oa[a_oa][dr] += 1
                            f_a_b_dab_oa_ob[a_b_oa_ob][dr] += 1
                            f_dab_oa_ob[oa_ob][dr] += 1
                    else:
                        if dab <= distances[dr]:
                            f_dab[dr] += 1
                            f_a_dab[a][dr] += 1
                            f_a_dab[b][dr] += 1
                            f_a_b_dab[a_b][dr] += 1
                            # For local only do amino acids #
                            f_a_dab_oa[a_oa][dr] += 1
                            f_a_b_dab_oa_ob[a_b_oa_ob][dr] += 1
                            f_dab_oa_ob[oa_ob][dr] += 1
        # Except... #
        except: pass

    return f_dab, f_a_dab, f_a_dab_oa, f_a_b_dab, f_dab_oa_ob, f_a_b_dab_oa_ob

def adjust_distance_to_bin(distance,bin_distance):
    """
    This function adjusts distance to best bin.
    
    """

    # Initialize #
    distances = []

    for i in numpy.arange(numpy.floor(distance), numpy.floor(distance) + 1, bin_distance):
        if i <= distance:
            distances.append(i)

    return distances[-1]


def approach_pmf(a_b_oa_ob, pmf, symbol):

    # Initialize #
    kB = 8.31441e-3 # Boltzmann constant
    T = 298.0       # standard temperature
    aminoacids = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
    approached_pmf = copy.copy(pmf[a_b_oa_ob])
    a_oa, b_ob = a_b_oa_ob.split(";")
    if len(a_oa.split("-"))>1:
       a, hydrophobicity, exposure, secondary_structure = a_oa.split("-")
       env=1
    else:
       a= a_oa
       b= b_ob
       env=0

    # For each distance bin... #
    for dr in range(len(approached_pmf)):
            # Approach PMF for amino acid at distance bin #
            # Function to determine best transition amino acid #
            f_max = []
            for aa in aminoacids:
                # Skip if amino acid is a... #
                if aa == a: continue
                if (env):
                  hydrophobicity = "N"
                  if aminoacids_polarity_boolean[aminoacids3to1[aa]]:
                    hydrophobicity = "P"
                  oaa = hydrophobicity + "-" + exposure + "-" + secondary_structure
                  aa_oaa = aa + "-" + oaa
                  aa_b_oaa_ob = aa_oaa + ";" + b_ob
                else:
                  aa_b_oaa_ob = aa + ";" + b
                if aa_b_oaa_ob in pmf:
                    if pmf[aa_b_oaa_ob][dr] is not None:
                        if symbol == "+":
                            f_max.append((aa, get_transition_probability(A=a, B=aa), pmf[aa_b_oaa_ob][dr], get_transition_probability(A=a, B=aa) * numpy.exp(pmf[aa_b_oaa_ob][dr] / (kB * T))))
                        else:
                            f_max.append((aa, get_transition_probability(A=a, B=aa), pmf[aa_b_oaa_ob][dr], get_transition_probability(A=a, B=aa) * numpy.exp(-pmf[aa_b_oaa_ob][dr] / (kB * T))))
            # Skip if no contacts were found #
            if len(f_max) > 0:
                # Get transition probability and energy of amino acid maximizing previous function #
                f_max.sort(key=lambda x: x[-1], reverse=True)
                max_aa = f_max[0][0]
                max_probability = f_max[0][1]
                max_energy = f_max[0][2]
                # Fill the summation using Taylor's approach #
                summation = 0.0
                for aa in aminoacids:
           	         # Skip if amino acid is a... #
                    if aa == a: continue
                    # Skip if amino acid is b max... #
                    if aa == max_aa: continue
                    if (env):
                      hydrophobicity = "N"
                      if aminoacids_polarity_boolean[aminoacids3to1[aa]]:
                        hydrophobicity = "P"
                      oaa = hydrophobicity + "-" + exposure + "-" + secondary_structure
                      aa_oaa = aa + "-" + oaa
                      aa_b_oaa_ob = aa_oaa + ";" + b_ob
                    else:
                      aa_b_oaa_ob = aa + ";" + b
                    if aa_b_oaa_ob in pmf:
                        if pmf[aa_b_oaa_ob][dr] is not None:
                            probability = get_transition_probability(A=a, B=aa)
                            energy = pmf[aa_b_oaa_ob][dr]
                            if symbol == "+":
                                summation += (probability / max_probability) * numpy.exp((energy - max_energy) / (kB * T))
                            else:
                                summation += (probability / max_probability) * numpy.exp(-((energy - max_energy) / (kB * T)))
                if approached_pmf[dr] is None:
                   if symbol == "+":
                    approached_pmf[dr] = kB * T * ((numpy.log(max_probability) + (max_energy / (kB * T)) + summation))
                   else:
                    approached_pmf[dr] = -kB * T * ((numpy.log(max_probability) - (max_energy / (kB * T)) + summation))
                else:
                   if symbol == "+":
                    #if (max_probability/get_transition_probability(A=a,B=a)) * numpy.exp( (pmf[a_b_oa_ob][dr] - max_energy) / (kB * T) ) < 1 :
                       #approached_pmf[dr] +=  kB * T * ( (numpy.log(get_transition_probability(A=a,B=a)) + (max_probability/get_transition_probability(A=a,B=a)) * numpy.exp( (pmf[a_b_oa_ob][dr] - max_energy) / (kB * T) ) * ( 1 + summation ) )
                    if (max_probability * numpy.exp(-(pmf[a_b_oa_ob][dr] - max_energy) / (kB * T) )) < 1 :
                       approached_pmf[dr] +=  kB * T * (  max_probability * numpy.exp(-(pmf[a_b_oa_ob][dr] - max_energy) / (kB * T) ) * ( 1 + summation ) )
                   else:
                    #if (max_probability/get_transition_probability(A=a,B=a)) * numpy.exp(-(pmf[a_b_oa_ob][dr] - max_energy) / (kB * T) ) < 1 :
                       #approached_pmf[dr] += -kB * T * ( (numpy.log(get_transition_probability(A=a,B=a)) + (max_probability/get_transition_probability(A=a,B=a)) * numpy.exp(-(pmf[a_b_oa_ob][dr] - max_energy) / (kB * T) ) * ( 1 + summation ) )
                    if (max_probability * numpy.exp(+(pmf[a_b_oa_ob][dr] - max_energy) / (kB * T) )) < 1 :
                       approached_pmf[dr] += -kB * T * (  max_probability * numpy.exp(+(pmf[a_b_oa_ob][dr] - max_energy) / (kB * T) ) * ( 1 + summation ) )


    return approached_pmf

def get_transition_probability(A, B):
    """
    This function returns the transition probability of an amino acid "A"
    to an amino acid "B" (according to BLOSUM62 scoring matrix).
    
    """
    
    l = 0.347
    aminoacids = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
    frequences = [0.074, 0.052, 0.045, 0.054, 0.025, 0.034, 0.054, 0.074, 0.026, 0.068, 0.099, 0.058, 0.025, 0.047, 0.039, 0.057, 0.051, 0.013, 0.034, 0.073]
    blosum62   = [['4', '-1', '-2', '-2', '0', '-1', '-1', '0', '-2', '-1', '-1', '-1', '-1', '-2', '-1', '1', '0', '-3', '-2', '0'],
                 ['-1', '5', '0', '-2', '-3', '1', '0', '-2', '0', '-3', '-2', '2', '-1', '-3', '-2', '-1', '-1', '-3', '-2', '-3'],
                 ['-2', '0', '6', '1', '-3', '0', '0', '0', '1', '-3', '-3', '0', '-2', '-3', '-2', '1', '0', '-4', '-2', '-3'],
                 ['-2', '-2', '1', '6', '-3', '0', '2', '-1', '-1', '-3', '-4', '-1', '-3', '-3', '-1', '0', '-1', '-4', '-3', '-3'],
                 ['0', '-3', '-3', '-3', '9', '-3', '-4', '-3', '-3', '-1', '-1', '-3', '-1', '-2', '-3', '-1', '-1', '-2', '-2', '-1'],
                 ['-1', '1', '0', '0', '-3', '5', '2', '-2', '0', '-3', '-2', '1', '0', '-3', '-1', '0', '-1', '-2', '-1', '-2'],
                 ['-1', '0', '0', '2', '-4', '2', '5', '-2', '0', '-3', '-3', '1', '-2', '-3', '-1', '0', '-1', '-3', '-2', '-2'],
                 ['0', '-2', '0', '-1', '-3', '-2', '-2', '6', '-2', '-4', '-4', '-2', '-3', '-3', '-2', '0', '-2', '-2', '-3', '-3'],
                 ['-2', '0', '1', '-1', '-3', '0', '0', '-2', '8', '-3', '-3', '-1', '-2', '-1', '-2', '-1', '-2', '-2', '2', '-3'],
                 ['-1', '-3', '-3', '-3', '-1', '-3', '-3', '-4', '-3', '4', '2', '-3', '1', '0', '-3', '-2', '-1', '-3', '-1', '3'],
                 ['-1', '-2', '-3', '-4', '-1', '-2', '-3', '-4', '-3', '2', '4', '-2', '2', '0', '-3', '-2', '-1', '-2', '-1', '1'],
                 ['-1', '2', '0', '-1', '-3', '1', '1', '-2', '-1', '-3', '-2', '5', '-1', '-3', '-1', '0', '-1', '-3', '-2', '-2'],
                 ['-1', '-1', '-2', '-3', '-1', '0', '-2', '-3', '-2', '1', '2', '-1', '5', '0', '-2', '-1', '-1', '-1', '-1', '1'],
                 ['-2', '-3', '-3', '-3', '-2', '-3', '-3', '-3', '-1', '0', '0', '-3', '0', '6', '-4', '-2', '-2', '1', '3', '-1'],
                 ['-1', '-2', '-2', '-1', '-3', '-1', '-1', '-2', '-2', '-3', '-3', '-1', '-2', '-4', '7', '-1', '-1', '-4', '-3', '-2'],
                 ['1', '-1', '1', '0', '-1', '0', '0', '0', '-1', '-2', '-2', '0', '-1', '-2', '-1', '4', '1', '-3', '-2', '-2'],
                 ['0', '-1', '0', '-1', '-1', '-1', '-1', '-2', '-2', '-1', '-1', '-1', '-1', '-2', '-1', '1', '5', '-2', '-2', '0'],
                 ['-3', '-3', '-4', '-4', '-2', '-2', '-3', '-2', '-2', '-3', '-2', '-3', '-1', '1', '-4', '-3', '-2', '11', '2', '-3'],
                 ['-2', '-2', '-2', '-3', '-2', '-1', '-2', '-3', '2', '-1', '-1', '-2', '-1', '3', '-3', '-2', '-2', '2', '7', '-1'],
                 ['0', '-3', '-3', '-3', '-1', '-2', '-2', '-3', '-3', '3', '1', '-2', '1', '-1', '-2', '-2', '0', '-3', '-1', '4']]

    return frequences[aminoacids.index(A)] * frequences[aminoacids.index(B)] * numpy.exp(l * int(blosum62[aminoacids.index(A)][aminoacids.index(B)]))


def smooth_bin(array, position, n):
    """
    This function smooths a given potential PMF bin by averaging it with
    the PMFs from the bins around +/- "n" positions.
    
    """

    values = []
    first = max(0,position-n)
    last  = min(position+n+1,len(array))
    for i in range(first,last):
        if array[i] != None: values.append(array[i])
    if len(values)>0: 
       return numpy.mean(values)
    else:
       return None

def calculate_zscores(pmf, pmf_name, bin_distance):

    # Initialize #
    zpmf = {}
    aminoacids = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
    max_contact_distance = float(config.get("Parameters", "max_contact_distance"))
    distances = list(numpy.arange(0, max_contact_distance + bin_distance, bin_distance)) # +bin_distance because of how range works


    # Z-score PMF 3DC #
    if pmf_name == "pmf_3dc":
        # For each environment pair... #
        for oa_ob in pmf:
            oa, ob = oa_ob.split(";")
            hydrophobicity, exposure, secondary_structure = oa.split("-")
            zpmf[oa_ob] = [None] * len(distances)
            for dr in range(len(distances)):
                if pmf[oa_ob][dr] is not None:
                    energy = []
                    for aa in aminoacids:
                        hydrophobicity = "N"
                        if aminoacids_polarity_boolean[aminoacids3to1[aa]]:
                            hydrophobicity = "P"
                        oaa = hydrophobicity + "-" + exposure + "-" + secondary_structure
                        oaa_ob = oaa + ";" + ob
                        if oaa_ob in pmf:
                            if pmf[oaa_ob][dr] is not None:
                                energy.append(pmf[oaa_ob][dr])
                    mean = numpy.mean(energy)
                    std = numpy.std(energy)
                    if std != 0.0:
                        zpmf[oa_ob][dr] = (pmf[oa_ob][dr] - mean) / std
        return zpmf
    # Z-score PMF S3DC or S3DC dd #
    if pmf_name == "pmf_s3dc" or pmf_name == "pmf_s3dc_dd":
        # For each residue-environment pair... #
        for a_b_oa_ob in pmf:
            a_oa, b_ob = a_b_oa_ob.split(";")
            a, hydrophobicity, exposure, secondary_structure = a_oa.split("-")
            zpmf[a_b_oa_ob] = [None] * len(distances)
            for dr in range(len(distances)):
                if pmf[a_b_oa_ob][dr] is not None:
                    energy = []
                    for aa in aminoacids:
                        hydrophobicity = "N"
                        if aminoacids_polarity_boolean[aminoacids3to1[aa]]:
                            hydrophobicity = "P"
                        oaa = hydrophobicity + "-" + exposure + "-" + secondary_structure
                        aa_oaa = aa + "-" + oaa
                        aa_b_oaa_ob = aa_oaa + ";" + b_ob
                        if aa_b_oaa_ob in pmf:
                            if pmf[aa_b_oaa_ob][dr] is not None:
                                energy.append(pmf[aa_b_oaa_ob][dr])
                    mean = numpy.mean(energy)
                    std = numpy.std(energy)
                    if std != 0.0:
                        zpmf[a_b_oa_ob][dr] = (pmf[a_b_oa_ob][dr] - mean) / std
        return zpmf
    # Z-score PMF S3DC di #
    if pmf_name == "pmf_s3dc_di":
        # For each residue-environment pair... #
        for a_b_oa_ob in pmf:
            a_oa, b_ob = a_b_oa_ob.split(";")
            a, hydrophobicity, exposure, secondary_structure = a_oa.split("-")
            zpmf[a_b_oa_ob] = [None]
            if pmf[a_b_oa_ob][0] is not None:
                energy = []
                for aa in aminoacids:
                    hydrophobicity = "N"
                    if aminoacids_polarity_boolean[aminoacids3to1[aa]]:
                        hydrophobicity = "P"
                    oaa = hydrophobicity + "-" + exposure + "-" + secondary_structure
                    aa_oaa = aa + "-" + oaa
                    aa_b_oaa_ob = aa_oaa + ";" + b_ob
                    if aa_b_oaa_ob in pmf:
                        if pmf[aa_b_oaa_ob][0] is not None:
                            energy.append(pmf[aa_b_oaa_ob])
                mean = numpy.mean(energy)
                std = numpy.std(energy)
                if std != 0.0:
                    zpmf[a_b_oa_ob][0] = (pmf[a_b_oa_ob][0] - mean) / std
        return zpmf
    # Z-score PMF local #
    if pmf_name == "pmf_local":
        # For each residue-environment... #
        for a_oa in pmf:
            a, hydrophobicity, exposure, secondary_structure = a_oa.split("-")
            zpmf[a_oa] = [None]
            if pmf[a_oa][0] is not None:
                energy = []
                for aa in aminoacids:
                    hydrophobicity = "N"
                    if aminoacids_polarity_boolean[aminoacids3to1[aa]]:
                        hydrophobicity = "P"
                    oaa = hydrophobicity + "-" + exposure + "-" + secondary_structure
                    aa_oaa = aa + "-" + oaa
                    if aa_oaa in pmf:
                        if pmf[aa_oaa][0] is not None:
                            energy.append(pmf[aa_oaa])
                mean = numpy.mean(energy)
                std = numpy.std(energy)
                if std != 0.0:
                    zpmf[a_oa][0] = (pmf[a_oa][0] - mean) / std
        return zpmf
    # Z-score PMF pair #
    if pmf_name == "pmf_pair":
        # For each residue pair... #
        for a_b in pmf:
            a, b = a_b.split(";")
            zpmf[a_b] = [None] * len(distances)
            for dr in range(len(distances)):
                if pmf[a_b][dr] is not None:
                    energy = []
                    for aa in aminoacids:
                        aa_b = aa + ";" + b
                        if aa_b in pmf:
                            if pmf[aa_b][dr] is not None:
                                energy.append(pmf[aa_b][dr])
                    mean = numpy.mean(energy)
                    std = numpy.std(energy)
                    if std != 0.0:
                        zpmf[a_b][dr] = (pmf[a_b][dr] - mean) / std
        return zpmf

    return None

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Get statistical potentials #
    pmf_3d, pmf_3dc, pmf_s3dc, pmf_s3dc_dd, pmf_s3dc_di, pmf_local, pmf_pair, distances = get_statistical_potentials(os.path.abspath(options.input_file), options.approach, options.smooth, options.zscores, options.computation, os.path.abspath(options.dummy_dir))

    # Initialize #
    potentials_obj = Potentials()
    # Assign the potentials #
    potentials_obj._pmf_3d = pmf_3d
    potentials_obj._pmf_3dc = pmf_3dc
    potentials_obj._pmf_s3dc = pmf_s3dc
    potentials_obj._pmf_s3dc_dd = pmf_s3dc_dd
    potentials_obj._pmf_s3dc_di = pmf_s3dc_di
    potentials_obj._pmf_local = pmf_local
    potentials_obj._pmf_pair = pmf_pair
    potentials_obj._distances = distances

    # Write output #
    if options.output_file is not None:
        potentials_obj.write(os.path.abspath(options.output_file))
    else:
        sys.stdout.write("pmf_3d\t%s\n" % json.dumps(potentials_obj._pmf_3d, separators=(",", ":")))
        sys.stdout.write("pmf_3dc\t%s\n" % json.dumps(potentials_obj._pmf_3dc, separators=(",", ":")))
        sys.stdout.write("pmf_s3dc\t%s\n" % json.dumps(potentials_obj._pmf_s3dc, separators=(",", ":")))
        sys.stdout.write("pmf_s3dc_dd\t%s\n" % json.dumps(potentials_obj._pmf_s3dc_dd, separators=(",", ":")))
        sys.stdout.write("pmf_s3dc_di\t%s\n" % json.dumps(potentials_obj._pmf_s3dc_di, separators=(",", ":")))
        sys.stdout.write("pmf_local\t%s\n" % json.dumps(potentials_obj._pmf_local, separators=(",", ":")))
        sys.stdout.write("pmf_pair\t%s\n" % json.dumps(potentials_obj._pmf_pair, separators=(",", ":")))
        sys.stdout.write("distances\t%s\n" % json.dumps(potentials_obj._distances, separators=(",", ":")))

