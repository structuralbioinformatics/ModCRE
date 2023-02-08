import os, sys, re
import ConfigParser
import optparse
import copy
import itertools
import matplotlib.patches as mpatches
import matplotlib.pyplot as plot
import numpy
import pandas
import seaborn

# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.dirname(__file__))

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Imports my functions #
import functions

# Import jbonet's module #
from SBI.data import aminoacids3to1, aminoacids_polarity_boolean, nitrogenous_bases

# Import my modules #
import spotentials

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python benchmark.py --pbm=pbm_dir --pdb=pdb_dir [--dummy=dummy_dir -f folds -o output_dir -r randoms --start=start_step --stop=stop_step -v]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="{directory}")
    parser.add_option("--pbm", action="store", type="string", dest="pbm_dir", help="PBM directory (output directory from pbm.py)", metavar="{directory}")
    parser.add_option("--pdb", action="store", type="string", dest="pdb_dir", help="PDB directory (output directory from pdb.py)", metavar="{directory}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")

    (options, args) = parser.parse_args()

    if options.pbm_dir is None or options.pdb_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Create output dir #
    if not os.path.exists(options.output_dir):
        os.makedirs(os.path.join(options.output_dir))
    # Verbose #
    if options.verbose: sys.stdout.write("Create figures...\n\n")
    # Initialize #
    aminoacids = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
    aminoacid_environments = []
    dinucleotides = ["".join(i) for i in itertools.product(list("ACGT"), repeat=2)]
    dinucleotide_environments = []
    # Get aminoacid environment #
    for aminoacid in aminoacids:
        # Get hydrophobicity #
        hydrophobicity = "N"
        if aminoacids_polarity_boolean[aminoacids3to1[aminoacid]]:
            hydrophobicity = "P"
        # For degree of exposure... #
        for degree_of_exposure in ["B", "E"]:
            # For secondary structure... #
            for secondary_structure in ["C", "E", "H"]:
                aminoacid_environments.append("%s-%s-%s-%s" % (aminoacid, hydrophobicity, degree_of_exposure, secondary_structure))
    # Get dinucleotide environment #
    for dinucleotide in dinucleotides:
        # For DNA groove... #
        for dna_groove in ["A", "I"]:
            # For chemical group... #
            for chemical_group in ["B", "N"]:
                dinucleotide_environments.append("%s-%s-%s-%s-%s" % (dinucleotide, "".join(nitrogenous_bases[nucleotide] for nucleotide in dinucleotide), "F", dna_groove, chemical_group))
    # Seaborn's context #
    seaborn.set(context="paper")
    # Initialize #
    heatmap_file = os.path.join(options.output_dir, "general.heatmaps.svg")
    # Skip if heatmap file already exists #
    if not os.path.exists(heatmap_file):
        # Initialize plot #
        fig, axes = plot.subplots(nrows=2, ncols=2, sharex=True, sharey=True, figsize=(7.44094, 7.44094))
        fig.suptitle("General statistical potentials")
        cbar_ax = fig.add_axes([.925, .25, .015, .5])
        # For each row... #
        for i, ax in enumerate(axes.flat):
            # Initialize #
            data = []
            # Get statistical potentials #
            if i == 0:
                ax.set_title("PDB")
                potentials = spotentials.Potentials(os.path.join(options.pdb_dir, "potentials", "general.txt"), "s3dc_dd")
            if i == 1:
                ax.set_title("PDB (Taylor's approach)")
                potentials = spotentials.Potentials(os.path.join(options.pdb_dir, "potentials", "general.taylor.txt"), "s3dc_dd")
            if i == 2:
                ax.set_title("PDB+PBM")
                potentials = spotentials.Potentials(os.path.join(options.pbm_dir, "potentials", "general.txt"), "s3dc_dd")
            if i == 3:
                ax.set_title("PDB+PBM (Taylor's approach)")
                potentials = spotentials.Potentials(os.path.join(options.pbm_dir, "potentials", "general.taylor.txt"), "s3dc_dd")
            # For each aminoacid environment... #
            for aminoacid_environment in aminoacid_environments:
                # Initialize #
                data.append([])
                # For each dinucleotide environment... #
                for dinucleotide_environment in dinucleotide_environments:
                    # Initialize #
                    triad = "%s;%s" % (aminoacid_environment, dinucleotide_environment)
                    # If triad exists... #
                    if triad in potentials._pmf_s3dc_dd:
                        try:
                            data[-1].append(min([k for k in potentials._pmf_s3dc_dd[triad] if k is not None]))
                        except:
                            data[-1].append(numpy.nan)
            # Plot heatmap #
            seaborn.heatmap(pandas.DataFrame(data), ax=ax, cbar=i == 0, vmin=-3, vmax=3, cbar_ax=None if i else cbar_ax, cmap=seaborn.cubehelix_palette(60, hue=1, dark=0.1, light=0.9, as_cmap=True, reverse=True))
            ax.set_xticks([2, 6, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 62])
            ax.set_xticklabels(dinucleotides, rotation="vertical")
            ax.set_yticks([3, 9, 15, 21, 27, 33, 39, 45, 51, 57, 63, 69, 75, 81, 87, 93, 99, 105, 111, 117])
            ax.set_yticklabels(sorted(aminoacids, reverse=True), rotation="horizontal")
            if i == 2 or i == 3:
                ax.set_xlabel("Dinucleotides")
            if i == 0 or i == 2:
                ax.set_ylabel("Amino acids")
        # Save figure #
        fig.savefig(heatmap_file, bbox_inches="tight", pad_inches=0.1, transparent=True)
    # Initialize #
    families = ["AFT", "AP2", "ARID/BRIGHT", "bHLH", "bZIP", "C2H2 ZF", "C2HC ZF", "CxxC",
                "DM", "E2F", "Ets", "Forkhead", "GATA", "GCM", "Homeodomain", "IRF", "Myb/SANT",
                "NAC/NAM", "Ndt80/PhoG", "Nuclear receptor", "Paired box", "POU", "Rap1", "Rel",
                "Runt", "SMAD", "Sox", "T-box", "WRKY", "Zinc cluster"]
    pdb_families = {}
    # For each line... #
    for line in functions.parse_file(os.path.join(options.pbm_dir, "families.txt")):
        if line.startswith("#"): continue
        line = line.split(";")
        if line[-1] in families:
            pdb_families.setdefault(line[-1], [])
            pdb_families[line[-1]].append(line[0])
    # For each family... #
    for family in families:
        # Initialize #
        heatmap_file = os.path.join(options.output_dir, "%s.heatmaps.svg" % family.replace("/", "_"))
        # Skip if heatmap file already exists #
        if not os.path.exists(heatmap_file):
            # Initialize plot #
            fig, axes = plot.subplots(nrows=2, ncols=2, sharex=True, sharey=True, figsize=(7.44094, 7.44094))
            fig.suptitle("Family statistical potentials (%s)" % family)
            cbar_ax = fig.add_axes([.925, .25, .015, .5])
            # For each row... #
            for i, ax in enumerate(axes.flat):
                # Initialize #
                data = []
                # For aminoacid environment... #
                for aminoacid_environment in aminoacid_environments:
                    # Initialize #
                    data.append([])
                    # For dinucleotide environment... #
                    for dinucleotide_environment in dinucleotide_environments:
                        # Initialize #
                        data[-1].append(numpy.nan)
                # For each PDB chain... #
                for pdb_chain in pdb_families[family]:
                    # Get statistical potentials #
                    if i == 0:
                        ax.set_title("PDB")
                        potentials = spotentials.Potentials(os.path.join(options.pdb_dir, "potentials", "%s.txt" % pdb_chain), "s3dc_dd")
                    if i == 1:
                        ax.set_title("PDB (Taylor's approach)")
                        potentials = spotentials.Potentials(os.path.join(options.pdb_dir, "potentials", "%s.taylor.txt" % pdb_chain), "s3dc_dd")
                    if i == 2:
                        ax.set_title("PDB+PBM")
                        potentials = spotentials.Potentials(os.path.join(options.pbm_dir, "potentials", "%s.txt" % pdb_chain), "s3dc_dd")
                    if i == 3:
                        ax.set_title("PDB+PBM (Taylor's approach)")
                        potentials = spotentials.Potentials(os.path.join(options.pbm_dir, "potentials", "%s.taylor.txt" % pdb_chain), "s3dc_dd")
                    # For aminoacid environment... #
                    for j in range(len(aminoacid_environments)):
                        # For dinucleotide environment... #
                        for k in range(len(dinucleotide_environments)):
                            # Initialize #
                            triad = "%s;%s" % (aminoacid_environments[j], dinucleotide_environments[k])
                            # If triad exists... #
                            if triad in potentials._pmf_s3dc_dd:
                                try:
                                    value = min([l for l in potentials._pmf_s3dc_dd[triad] if l is not None])
                                except:
                                    value = None
                                if value is None: continue
                                if data[j][k] is numpy.nan:
                                    data[j][k] = copy.copy(value)
                                else:
                                    data[j][k] = min([data[j][k], value])
                # Plot heatmap #
                seaborn.heatmap(pandas.DataFrame(data), ax=ax, cbar=i == 0, vmin=-3, vmax=3, cbar_ax=None if i else cbar_ax, cmap=seaborn.cubehelix_palette(60, hue=1, dark=0.1, light=0.9, as_cmap=True, reverse=True))
                ax.set_xticks([2, 6, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 62])
                ax.set_xticklabels(dinucleotides, rotation="vertical")
                ax.set_yticks([3, 9, 15, 21, 27, 33, 39, 45, 51, 57, 63, 69, 75, 81, 87, 93, 99, 105, 111, 117])
                ax.set_yticklabels(sorted(aminoacids, reverse=True), rotation="horizontal")
                if i == 2 or i == 3:
                    ax.set_xlabel("Dinucleotides")
                if i == 0 or i == 2:
                    ax.set_ylabel("Amino acids")
            # Save figure #
            fig.savefig(heatmap_file, bbox_inches="tight", pad_inches=0.1, transparent=True)
