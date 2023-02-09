import os, sys, re
import ConfigParser
import matplotlib.pyplot as plot
import numpy
import optparse
from scipy.stats import wilcoxon

#-------------#
# Options     #
#-------------#

def parse_options():
    '''
    This function parses the command line arguments and returns an optparse object.
    '''

    parser = optparse.OptionParser("figures.py [--dummy=DUMMY_DIR] -i INPUT_DIR [-o OUTPUT_DIR -v]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="DUMMY_DIR")
    parser.add_option("-i", action="store", type="string", dest="input_dir", help="Input directory (i.e. from statistics.py)", metavar="INPUT_DIR")
    parser.add_option("-o", "--output-dir", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="OUTPUT_DIR")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")

    (options, args) = parser.parse_args()

    if options.input_dir == None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

#-------------#
# Functions   #
#-------------#

def significance(text, X_location, Y_location, X_mean, Y_mean, X_stderr, Y_stderr):
    # Find max y point #
    y = 0.005 + max(X_mean + X_stderr, Y_mean + Y_stderr)
    # Define props #
    props = {'arrowstyle' : "-", 'connectionstyle' : "bar"}
    # Create labels #
    plot.annotate(text, xy=(X_location, y + 0.03), size=10)
    plot.annotate("", xy=(X_location, y), xytext=(Y_location, y), arrowprops=props)

def significance_family(significant, X_location, Y_location, X_mean, Y_mean, X_stderr, Y_stderr):
    if significant:
        # Find max y point #
        y = max(X_mean + X_stderr, Y_mean + Y_stderr)
        # Define props #
        props = {'arrowstyle' : "-", 'connectionstyle' : "bar"}
        # Create labels #
        plot.annotate("*", xy=((X_location + Y_location) / 2 - 0.03, y), zorder=10)
        
def get_ppv_and_cov(values, threshold):
    # Initialize #
    tp = 0.0
    fp = 0.0
    tn = 0.0
    fn = 0.0
    ppv = 0
    cov = 0
    # For each line... #
    for score, boolean in values:
        if score >= threshold:
            if boolean:
                tp += 1
            else:
                fp += 1
        else:
            if boolean:
                fn += 1
            else:
                tn += 1
    # Calculate ppv and cov #
    if tp > 0 or fp > 0:
        ppv = tp / (tp + fp)
        cov = tp / (tp + fn)

    return (ppv, cov)

#-------------#
# Main        #
#-------------#

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = "%s/config.ini" % os.path.dirname(__file__)
config.read(config_file)

# Add "src" to sys.path #
sys.path.append(config.get("Paths", "src_path"))

# Imports my functions #
import functions as functions

# Arguments & Options #
options = parse_options()
dummy_dir = os.path.abspath(options.dummy_dir)
input_dir = os.path.abspath(options.input_dir)
output_dir = os.path.abspath(options.output_dir)
verbose = options.verbose

# Create output directory #
if not os.path.exists(output_dir): os.makedirs(output_dir)

#################
# Plot: Barplot #
#################
if verbose: sys.stdout.write("\nCreate barplot...\n")
# Initialize #
n = 3
bar_width = 0.35
error_config = {'ecolor': "black", 'width' : 1}
plot.figure(figsize=(3.50394, 5.25591))
prop = {'size' : 8}
# Start plot #
plot.subplot(1,1,1, autoscale_on=False, aspect="equal")
# Adjust grid #
plot.subplots_adjust(left=0.2)
# Get data #
pdb_general = []
pdb_general_taylor = []
pbm_general = []
pbm_general_taylor = []
pbm_family = []
pbm_family_taylor = []
# Read files #
for line in functions.parse_csv_file(os.path.join(input_dir, "pdb.general.txt")):
    # Skip lines #
    if line[0] == "Family" or line[0] == "All": continue
    # Add line #
    if line[2] == "NA": line[2] = None
    else: line[2] = float(line[2])
    line[3] = int(line[3])
    pdb_general.append(line)
# Read files #
for line in functions.parse_csv_file(os.path.join(input_dir, "taylor.pdb.general.txt")):
    # Skip lines #
    if line[0] == "Family" or line[0] == "All": continue
    # Add line #
    if line[2] == "NA": line[2] = None
    else: line[2] = float(line[2])
    line[3] = int(line[3])
    pdb_general_taylor.append(line)
# Read files #
for line in functions.parse_csv_file(os.path.join(input_dir, "pbm.general.txt")):
    # Skip lines #
    if line[0] == "Family" or line[0] == "All": continue
    # Add line #
    if line[2] == "NA": line[2] = None
    else: line[2] = float(line[2])
    line[3] = int(line[3])
    pbm_general.append(line)
# Read files #
for line in functions.parse_csv_file(os.path.join(input_dir, "taylor.pbm.general.txt")):
    # Skip lines #
    if line[0] == "Family" or line[0] == "All": continue
    # Add line #
    if line[2] == "NA": line[2] = None
    else: line[2] = float(line[2])
    line[3] = int(line[3])
    pbm_general_taylor.append(line)
# Read files #
for line in functions.parse_csv_file(os.path.join(input_dir, "pbm.family.txt")):
    # Skip lines #
    if line[0] == "Family" or line[0] == "All": continue
    # Add line #
    if line[2] == "NA": line[2] = None
    else: line[2] = float(line[2])
    line[3] = int(line[3])
    pbm_family.append(line)
# Read files #
for line in functions.parse_csv_file(os.path.join(input_dir, "taylor.pbm.family.txt")):
    # Skip lines #
    if line[0] == "Family" or line[0] == "All": continue
    # Add line #
    if line[2] == "NA": line[2] = None
    else: line[2] = float(line[2])
    line[3] = int(line[3])
    pbm_family_taylor.append(line)
# Get families #
families = [i[0] for i in pbm_family if i[2] is not None]
families_N = [i[-1] for i in pbm_family if i[2] is not None]
families_as_in_cell_paper = []
for i in range(len(families)):
    families_as_in_cell_paper.append(families[i])
    # Adjust names to Weirauch's Cell paper #
    if families[i] == "C2H2_ZF": families_as_in_cell_paper[-1] = "C2H2 ZF"
    if families[i] == "Myb_SANT": families_as_in_cell_paper[-1] = "Myb/SANT"
    if families[i] == "NAC_NAM": families_as_in_cell_paper[-1] = "NAC/NAM"
    if families[i] == "Nuclear_receptor": families_as_in_cell_paper[-1] = "Nuclear receptor"
    if families[i] == "Paired_box": families_as_in_cell_paper[-1] = "Paired box"
    if families[i] == "POU": families_as_in_cell_paper[-1] = "Homeo+Pou"
    if families[i] == "T_box": families_as_in_cell_paper[-1] = "T box"
    if families[i] == "Zinc_cluster": families_as_in_cell_paper[-1] = "Zinc cluster"
# Transform lists to numpy arrays #
pdb_general = numpy.array([i[2] for i in pdb_general if i[2] is not None])
pdb_general_taylor = numpy.array([i[2] for i in pdb_general_taylor if i[2] is not None])
pbm_general = numpy.array([i[2] for i in pbm_general if i[2] is not None])
pbm_general_taylor = numpy.array([i[2] for i in pbm_general_taylor if i[2] is not None])
pbm_family = numpy.array([i[2] for i in pbm_family if i[2] is not None])
pbm_family_taylor = numpy.array([i[2] for i in pbm_family_taylor if i[2] is not None])
# Get means and stds #
means = (pdb_general.mean(), pbm_general.mean(), pbm_family.mean())
stderrs = (pdb_general.std() / numpy.sqrt(len(pdb_general)), pbm_general.std() / numpy.sqrt(len(pbm_general)), pbm_family.std() / numpy.sqrt(len(pbm_family)))
means_taylor = (pdb_general_taylor.mean(), pbm_general_taylor.mean(), pbm_family_taylor.mean())
stderrs_taylor = (pdb_general_taylor.std() / numpy.sqrt(len(pdb_general_taylor)), pbm_general_taylor.std() / numpy.sqrt(len(pbm_general_taylor)), pbm_family_taylor.std() / numpy.sqrt(len(pbm_family_taylor)))
# Get wilcoxon tests #
wilcoxon_tests = (wilcoxon(pdb_general, pdb_general_taylor), wilcoxon(pbm_general, pbm_general_taylor), wilcoxon(pbm_family, pbm_family_taylor))
# Get "x" bars locations #
x_locations = []
for i in range(n):
    x_locations.append(i)
    x_locations.append(i + bar_width)
x_locations = numpy.array(x_locations)
# Get means values #
values = []
for i in range(n):
    values.append(means[i] * n * 2)
    values.append(means_taylor[i] * n * 2)
values = numpy.array(values)
# Get means errorbars  #
errorbars = []
for i in range(n):
    errorbars.append(stderrs[i] * n * 2)
    errorbars.append(stderrs_taylor[i] * n * 2)
errorbars = numpy.array(errorbars)
# Get bars colors #
colors = []
for i in range(n):
    colors.append("white")
    colors.append("grey")
colors = numpy.array(colors)
# Plot bars #
bars = plot.bar(x_locations, values, bar_width, color=colors, yerr=errorbars, error_kw=error_config)
# Plot significance #
for i in range(n):
    if wilcoxon_tests[i][-1] < 0.05:
        plot.text(bar_width, max([values[i] + errorbars[i], values[i + 1] + errorbars[i + 1]]), "*", size=10, horizontalalignment="center")
    if wilcoxon_tests[i][-1] < 0.01:
        plot.text(bar_width, max([values[i] + errorbars[i], values[i + 1] + errorbars[i + 1]]) + 0.30, "*", size=10, horizontalalignment="center")
    if wilcoxon_tests[i][-1] < 0.001:
        plot.text(bar_width, max([values[i] + errorbars[i], values[i + 1] + errorbars[i + 1]]) + 0.30 + 0.30, "*", size=10, horizontalalignment="center")
# Get plot labels texts and ticks #
plot.xticks((numpy.arange(n) + bar_width), ("general (PDB)", "general (PDB+PBM)", "family (PDB+PBM)"), size=10, horizontalalignment="right", rotation=45)
plot.yticks((numpy.array([0, 0.6, 1.2, 1.8, 2.4, 3])), (numpy.arange(0, 0.6, 0.1)), size=10)
plot.ylabel("AUPRC", size=10)
plot.legend([bars[1]], ["Taylor's approach"], loc=2, frameon=False, prop=prop)
plot.xlim(0, 3)
plot.ylim(0, 3)
# Save figure #
plot.savefig(os.path.join(output_dir, "barplot.png"))

########################
# Plot: PPV & Coverage #
########################
if verbose: sys.stdout.write("\nCreate PPV/Cov plot...\n")
# Start plot #
plot.close('all')
plot.figure(figsize=(3.50394, 5.25591))
plot.subplot(1,1,1, autoscale_on=False, aspect="equal")
# Adjust grid #
plot.subplots_adjust(left=0.2)
# Get data #
scores = numpy.arange(0, 1.01, step=0.01)
ppvs = numpy.empty((len(scores), len(families)))
covs = numpy.empty((len(scores), len(families)))
# For each family... #
for i in range(len(families)):
    if verbose: sys.stdout.write("\t%s...\n" % families[i])
    # Get values #
    values = [[float(line[5]), bool(int(line[6]))] for line in functions.parse_csv_file(os.path.join(input_dir, "taylor.pbm.family." + families[i] + ".txt.gz"), gz=True) if line[0] != "Fold"]
    # For each score... #
    for j in range(len(scores)):
        ppv, cov = get_ppv_and_cov(values, int(scores[j] * 100) / 100.0)
        ppvs[j][i] = ppv
        covs[j][i] = cov
# Get means and stds #
means_ppvs = tuple([i.mean() for i in ppvs])
stderrs_means = tuple([i.std() / numpy.sqrt(len(families)) for i in ppvs])
means_covs = tuple([i.mean() for i in covs])
stderrs_covs = tuple([i.std() / numpy.sqrt(len(families)) for i in covs])
# Create lines #
lines1 = plot.plot(scores, means_ppvs, color="black", linestyle="-")
lines2 = plot.plot(scores, means_covs, color="black", linestyle=":")
lines3 = plot.plot([0.95, 0.95], [-0.05, 1.05], color="red", linestyle="-")
# Plot error bars every 0.05 scores #
plot.errorbar([scores[i] for i in range(len(scores)) if scores[i] > 0 and scores[i] < 1 and int(scores[i] * 100) % 5 == 0],
              [means_ppvs[i] for i in range(len(scores)) if scores[i] > 0 and scores[i] < 1 and int(scores[i] * 100) % 5 == 0],
              yerr=[stderrs_means[i] for i in range(len(scores)) if scores[i] > 0 and scores[i] < 1 and int(scores[i] * 100) % 5 == 0], color="black", fmt='.')
plot.errorbar([scores[i] for i in range(len(scores)) if scores[i] > 0 and scores[i] < 1 and int(scores[i] * 100) % 5 == 0],
              [means_covs[i] for i in range(len(scores)) if scores[i] > 0 and scores[i] < 1 and int(scores[i] * 100) % 5 == 0],
              yerr=[stderrs_covs[i] for i in range(len(scores)) if scores[i] > 0 and scores[i] < 1 and int(scores[i] * 100) % 5 == 0], color="black", fmt='.')
# Labels texts, title and axis ticks
plot.xticks(numpy.arange(0, 1.2, 0.2), numpy.arange(0, 1.2, 0.2), size=10)
plot.yticks(numpy.arange(0, 1.2, 0.2), numpy.arange(0, 1.2, 0.2), size=10)
plot.xlabel("score", size=10)
plot.ylabel("PPV & coverage", size=10)
plot.legend((lines1[0], lines2[0]), ("PPV", "coverage"), loc=6, frameon=False, prop=prop)
plot.xlim(-0.05, 1.05)
plot.ylim(-0.05, 1.05)
# Save figure #
plot.savefig(os.path.join(output_dir, "ppv_cov.png"))

##################
# Plot: Families #
##################
if verbose: sys.stdout.write("\nCreate families barplots...\n")
# Initialize #
position = 1
plot.close('all')
plot.figure(figsize=(7.20472, 9.0059))
# For each family... #
for i in range(len(families)):
    # Get values #
    positives = [float(line[5]) for line in functions.parse_csv_file(os.path.join(input_dir, "pbm.family." + families[i] + ".txt.gz"), gz=True) if line[6] == "1"]
    positives_taylor = [float(line[5]) for line in functions.parse_csv_file(os.path.join(input_dir, "taylor.pbm.family." + families[i] + ".txt.gz"), gz=True) if line[6] == "1"]
    # Start plot #
    plot.subplot(5, 4, i + 1, autoscale_on=False, aspect="equal")
    plot.bar(0, pbm_family[i], bar_width, color="white")
    plot.bar(bar_width, pbm_family_taylor[i], bar_width, color="grey")
    # Plot significance #
    wilcoxon_test = wilcoxon(positives, positives_taylor)
    if wilcoxon_test[-1] < 0.05:
        plot.text(bar_width, max([pbm_family[i], pbm_family_taylor[i]]), "*", size=10, horizontalalignment="center")
    if wilcoxon_test[-1] < 0.01:
        plot.text(bar_width, max([pbm_family[i], pbm_family_taylor[i]]) + 0.05, "*", size=10, horizontalalignment="center")
    if wilcoxon_test[-1] < 0.001:
        plot.text(bar_width, max([pbm_family[i], pbm_family_taylor[i]]) + 0.05 + 0.05, "*", size=10, horizontalalignment="center")
    # Get current axis and plot labels texts and ticks #
    axis = plot.gca()
    axis.set_ylim((0, 1))
    axis.tick_params(axis='both', labelsize=10)
    axis.set_xticks([])
    if i == 0 or i % 4 == 0:
        axis.yaxis.set_label_position("left")
        axis.set_ylabel("AUPRC", size=10)
        axis.set_yticks((0, 1))
    else:
        axis.set_yticks([])
    plot.text(0.98, 0.90, "ppv=" + str(int(ppvs[95][i] * 100)) + "%", size=8, horizontalalignment="right")
    plot.text(0.98, 0.81, "cov=" + str(int(covs[95][i] * 100)) + "%", size=8, horizontalalignment="right")
    plot.text(0.98, 0.02, families_N[i], size=8, horizontalalignment="right")
    plot.text(0.5, 1.03, families_as_in_cell_paper[i], size=10, horizontalalignment="center")
# Save figure #
plot.savefig(os.path.join(output_dir, "families.png"))
