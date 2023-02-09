import sys, os
import numpy as np
import matplotlib.pyplot as plt
import itertools
import scipy.stats
from pylab import plot, show, savefig, xlim, figure, hold, ylim, legend, boxplot, setp, axes, text
import pandas as pd
import seaborn as sns
import matplotlib.patches as mpatches

#import pbm

#input_dir = os.path.abspath(sys.argv[1])
input_file = os.path.abspath(sys.argv[1])
output_dir = os.path.abspath(sys.argv[2])
pbm_dir = os.path.abspath(sys.argv[3])
curve = sys.argv[4]

#pots = ["pdb.general", "pdb.general.taylor", "pbm.general", "pbm.general.taylor", "pdb.family", "pdb.family.taylor", "pbm.family", "pbm.family.taylor"]

#families = ["AFT", "AP2", "APSES", "ARID/BRIGHT", "C2H2 ZF", "C2HC ZF", "CxC", "CxxC", "DM", "E2F", "Ets", "Forkhead", "GATA", "GCM", "HSF", "Homeodomain", "IRF", "MADS box", "MBD", "Myb/SANT", "NAC/NAM", "Ndt80/PhoG", "Nuclear receptor", "POU", "Paired box", "Rap1", "SMAD", "Sox", "T-box", "WRKY", "Zinc cluster", "bHLH", "bZIP"]
"""
# Initialize #
families = []
# For each TF... #
for file_name in os.listdir(os.path.join(os.path.abspath(pbm_dir), "tfs")):
    # Get TF object #
    tf_obj = pbm.TF(file_name=os.path.join(os.path.abspath(pbm_dir), "tfs", file_name))
    # For each sequence... #
    for i in range(len(tf_obj._sequences)):
        # Skip if k-mers file does not exist #
        kmers_file = os.path.join(os.path.abspath(pbm_dir), "kmers", tf_obj.get_id() + "." + str(i) + ".txt")
        if not os.path.exists(kmers_file): continue
        # Get formatted family name #
        #family = tf_obj.get_family(format=True)
        family = None
        family_list = tf_obj.get_family()
        if len(family_list) > len(tf_obj._sequences):
            continue
        try:
            if len(family_list) == 1:
                family = family_list[0]
            elif len(family_list) > 1:
                family = family_list[i]
        except:
            family = None
            continue
        # Add motif to family #
        if family == None:
            continue
        #families.setdefault(family, set())
        #families[family].add(("%s.%s" % (tf_obj.get_id(), str(i)), tf_obj._motifs[i]))
        if not family in families:
            families.append(family)


# Create heatmaps #
auprc_dict = {}
for fam in families:
    auprc_dict.setdefault(fam, {})
    for pot in pots:
        print("curve: " + curve)
        if curve == "roc":
            input_file = os.path.join(input_dir, pot + ".auroc.csv")
        elif curve == "pr":
            input_file = os.path.join(input_dir, pot + ".aucpr.csv")
        in_file = open(input_file, "r").readlines()[1:]
        for line in in_file:
            if line.startswith(fam):
                auprc = line.split(",")[1]
                auprc_dict[fam][pot] = float(auprc)
"""
data = open(input_file, "r").readlines()[0]

data_dict = {}
families = []
pots = ['PDB.General.null', 'PDB.General.Taylor', 'PDB+PBM.General.null', 'PDB+PBM.General.Taylor', 'PDB.family.null', 'PDB.family.Taylor', 'PDB+PBM.family.null', 'PDB+PBM.family.Taylor'] 

split_data = data.split("[[")[1:]
for sp in split_data:
    print(sp)
    
    data_p = str(sp.split(",")[0]).replace('''"''', "").replace(" ", "")
    if data_p.startswith("["):
        data_p = data_p[1:]
    fam_p = str(sp.split(",")[1]).replace('''"''', "").replace(" ", "")
    taylor = str(sp.split(",")[2]).replace('''"''', "").replace(" ", "")
    family = str(sp.split(",")[3][:-1]).replace('''"''', "").replace(" ", "")
    auroc = float(str(sp.split(",")[5])[:-1])
    auprc = str(sp.split(",")[7])[:-2]
    if auprc.endswith("]"):
        auprc = float(auprc[:-1])
    else:
        auprc = float(auprc)

    if not family in families:
        families.append(family)
    if not data_p + "." + fam_p + "." + taylor in pots:
        pots.append(data_p + "." + fam_p + "." + taylor)

    data_dict.setdefault(family, {})
    data_dict[family].setdefault(data_p + "." + fam_p + "." + taylor)
    print("curve: " + curve)
    if curve == "roc":
        data_dict[family][data_p + "." + fam_p + "." + taylor] = auroc
    elif curve == "pr":
        data_dict[family][data_p + "." + fam_p + "." + taylor] = auprc
    elif curve == "ratio":
        data_dict[family][data_p + "." + fam_p + "." + taylor] = auprc/auroc


    #print(str([data_p, fam_p, taylor, family, str(auroc), str(auprc)]))

    

print(str(pots))
df = pd.DataFrame({"potentials": pots})
for fam in sorted(families):
    print("fam: " + str(fam))
    value_list = []
    for pot in pots:
        value_list.append(data_dict[fam][pot])
    df[fam] = value_list

df = df.set_index("potentials")
df = df.transpose()
print(df.to_string())
values = df.as_matrix()

fig = figure(figsize=(20, 22))
sns.set(font_scale=3)
mask = df.isnull()
ax = sns.heatmap(df, annot=True, cmap="Reds", cbar=True, annot_kws={"size": 18}, mask=mask, )
#ax.set_xticklabels(potentials)
#ax.set_yticklabels(tf_families)
for item in ax.get_xticklabels():
    item.set_rotation(60)
    item.set_fontsize(20)
for item in ax.get_yticklabels():
    item.set_rotation(0)
    item.set_fontsize(22)
    item.set_fontweight('bold')
    #item.set_fontname("calibri")

plt.subplots_adjust(left=0.15)
fig.savefig(os.path.join(output_dir, "benchmark_heatmap_" + curve + ".png"))
print("plot created at: " + str(os.path.join(output_dir, "benchmark_heatmap_" + curve + ".png")))
plt.close(fig)


"""
# Create precision-recall curves #
color_list = ["firebrick", "orange", "gold", "lawngreen", "forestgreen", "c", "blue", "orchid"]
for fam in families:
    precision = {}
    x_vals = {}
    for pot in pots:
        precision.setdefault(pot, [])
        x_vals.setdefault(pot, [])
        pre_file = os.path.join(input_dir, pot + ".precision.csv")
        for line in open(pre_file, "r").readlines()[1:]:
            if line.startswith(fam):
                fields = line.split(",")
                try:
                    x = float(fields[2])
                    val = float(fields[3])
                    x_vals[pot].append(x)
                    precision[pot].append(val)
                    print(str(x) + "   " + str(val))
                except:
                    continue
            
    fig = figure(figsize=(10, 10))
    ax = axes()
    ax.patch.set_facecolor("white")
    ax.grid(color="lightgray")
    for i in range(len(pots)):
        pot = pots[i]
        col = color_list[i]
        ax.plot(x_vals[pot], precision[pot], color=col, linewidth=2)

    ax.set_xlabel("Threshold", fontsize=12)
    ax.set_ylabel("Precision", fontsize=12)
    ax.set_ylim(0, 1)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.suptitle("Precision " + fam.replace("/", "_").replace(" ", "_"), fontsize=14, fontweight="bold")

    pdb_general = mpatches.Patch(color='firebrick', label='PDB_general')
    pdb_general_taylor = mpatches.Patch(color='orange', label='PDB_general_taylor')
    pbm_general = mpatches.Patch(color='gold', label='PBM_general')
    pbm_general_taylor = mpatches.Patch(color='lawngreen', label='PBM_general_taylor')
    pdb_family = mpatches.Patch(color='forestgreen', label='PDB_family')
    pdb_family_taylor = mpatches.Patch(color='c', label='PDB_family_taylor')
    pbm_family = mpatches.Patch(color='blue', label='PBM_family')
    pbm_family_taylor = mpatches.Patch(color='orchid', label='PBM_family_taylor')
    handle_list = [pdb_general, pdb_general_taylor, pbm_general, pbm_general_taylor, pdb_family, pdb_family_taylor, pbm_family, pbm_family_taylor]
    plt.legend(handles=(handle_list), bbox_to_anchor=(1.00, 0.95), loc=2, borderaxespad=0., prop={'size':12})
    plt.subplots_adjust(right=0.75)


    fig.savefig(os.path.join(output_dir, "precision_" + fam.replace("/", "_").replace(" ", "_") + "_" + curve + ".png"))
    print("plot created at: " + str(os.path.join(output_dir, "precision_" + fam.replace("/", "_").replace(" ", "_") + "_" + curve + ".png")))
    plt.close(fig)

    recall = {}
    x_vals = {}
    for pot in pots:
        recall.setdefault(pot, [])
        x_vals.setdefault(pot, [])
        rec_file = os.path.join(input_dir, pot + ".recall.csv")
        for line in open(rec_file, "r").readlines()[1:]:
            if line.startswith(fam):
                fields = line.split(",")
                try:
                    x = float(fields[2])
                    val = float(fields[3])
                    x_vals[pot].append(x)
                    recall[pot].append(val)
                except:
                    continue

    fig = figure(figsize=(10, 10))
    ax = axes()
    ax.patch.set_facecolor("white")
    ax.grid(color="lightgray")
    for i in range(len(pots)):
        pot = pots[i]
        col = color_list[i]
        ax.plot(x_vals[pot], recall[pot], color=col, linewidth=2)

    ax.set_xlabel("Threshold", fontsize=12)
    ax.set_ylabel("Recall", fontsize=12)
    ax.set_ylim(0, 1)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.suptitle("Recall " + fam.replace("/", "_").replace(" ", "_"), fontsize=14, fontweight="bold")

    pdb_general = mpatches.Patch(color='firebrick', label='PDB_general')
    pdb_general_taylor = mpatches.Patch(color='orange', label='PDB_general_taylor')
    pbm_general = mpatches.Patch(color='gold', label='PBM_general')
    pbm_general_taylor = mpatches.Patch(color='lawngreen', label='PBM_general_taylor')
    pdb_family = mpatches.Patch(color='forestgreen', label='PDB_family')
    pdb_family_taylor = mpatches.Patch(color='c', label='PDB_family_taylor')
    pbm_family = mpatches.Patch(color='blue', label='PBM_family')
    pbm_family_taylor = mpatches.Patch(color='orchid', label='PBM_family_taylor')
    handle_list = [pdb_general, pdb_general_taylor, pbm_general, pbm_general_taylor, pdb_family, pdb_family_taylor, pbm_family, pbm_family_taylor]
    plt.legend(handles=(handle_list), bbox_to_anchor=(1.00, 0.95), loc=2, borderaxespad=0., prop={'size':12})
    plt.subplots_adjust(right=0.75)

    fig.savefig(os.path.join(output_dir, "recall_" + fam.replace("/", "_").replace(" ", "_") + "_" + curve + ".png"))
    print("plot created at: " + str(os.path.join(output_dir, "recall_" + fam.replace("/", "_").replace(" ", "_") + "_" + curve + ".png")))
    plt.close(fig)
        



"""
