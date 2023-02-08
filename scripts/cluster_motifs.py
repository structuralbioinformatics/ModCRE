import os, sys, re
import ConfigParser
import json
import optparse
import shutil
import subprocess
import imp
import math
import numpy as np
import pickle

# Import sklearn modules #
from sklearn.cluster import SpectralClustering
from sklearn import metrics

import tomtom


def merge_fimos_and_nucleotide_data(motif_list, pval_cutoffs):

    # Merge binding site data #
    merged_motif = {}
    merged_motif["binding_sites"] = {}
    for allele in ["ref", "alt"]:
        fimo_obj = fimo.Fimo()
        for motif in motif_list:
            for fh in motif["binding_sites"][allele].get_hits():
                fimo_obj.add_hit(fh)
        merged_motif["binding_sites"][allele] = fimo_obj

    return merged_motif


def get_affinity_between_clusters(cluster1, cluster2):

    max_pval = None
    for mot1 in cluster1:
        for mot2 in cluster2:
            pwm1 = mot1["complete_meme_pwm"]
            pwm2 = mot2["complete_meme_pwm"]
            tomtom_obj = tomtom.get_tomtom_obj(pwm1, pwm2)
            tomtom_hit_obj = tomtom_obj.get_hits()[0]
            pval = float(tomtom_hit_obj._p_value)
            #print(str(mot1["model"]) + "  " + str(mot2["model"]) + "  " + str(pval))
            if max_pval == None:
                max_pval = pval
            elif max_pval < pval:
                max_pval = pval

    affinity = 1.0/max_pval
    return affinity


def check_merged_clusters(merged_clusters):

    for clust in merged_clusters:
        if str(type(clust)) == "<type 'list'>":
            affinity = get_affinity_between_clusters(clust, clust)
        elif str(type(clust)) == "<type 'dict'>":
            affinity = get_affinity_between_clusters([clust], [clust])

        if affinity < 2:
            return False
        
    return True


def plot_adjacency_matrix(adjacency_matrix, output_dir):

    import matplotlib.pyplot as plt
    from pylab import plot, show, savefig, xlim, figure, hold, ylim, legend, boxplot, setp, axes, text
    import pandas as pd
    import seaborn as sns
    import time

    df = pd.DataFrame(data=np.array(adjacency_matrix))
    fig = figure(figsize=(50, 45))
    sns.set(font_scale=3)
    mask = df.isnull()
    try:
        ax = sns.heatmap(df, annot=False, cmap="Reds", cbar=True, annot_kws={"size": 18}, mask=mask, vmin=0.0, )
        plt.subplots_adjust(left=0.15)
        print("Adjacency matrix ploted in: " + os.path.join(output_dir, "adjacency_matrix.png"))
        fig.savefig(os.path.join(output_dir, "adjacency_matrix.png"))
        plt.close(fig)
    except:
        return


def cluster_motifs(motif_list):

    motif_list = sorted(motif_list)
    # Before checking pwm similarity, split the motif list in domains #
    clusters_dict = {}
    for mot in motif_list:
        domain = mot["domain"]
        if domain in clusters_dict.keys():
            clusters_dict[domain].append(mot)
        else:
            similar = False
            domain1 = int(domain.split(":")[0])
            domain2 = int(domain.split(":")[1])
            for dom in clusters_dict.keys():    
                dom1 = int(dom.split(":")[0])
                dom2 = int(dom.split(":")[1])
                if len(set(range(domain1, domain2)).intersection(set(range(dom1, dom2)))) > min([domain2-domain1, dom2-dom1])*0.8:
                    # this means that the two domains are quite similar #
                    similar = True
                    clusters_dict[dom].append(mot)
                    break
            if similar == False:
                clusters_dict[str(domain1) + ":" + str(domain2)] = [mot]

    final_clusters = []
    for domain in clusters_dict.keys():
        motif_list = clusters_dict[domain]
        adjacency_matrix = []
        for i in range(0, len(motif_list)):
            row = []
            for j in range(0, len(motif_list)):
                pwm1 = motif_list[i]["complete_meme_pwm"]
                pwm2 = motif_list[j]["complete_meme_pwm"]
                tomtom_obj = tomtom.get_tomtom_obj(pwm1, pwm2)
                tomtom_hit_obj = tomtom_obj.get_hits()[0]
                affinity = -np.log10(float(tomtom_hit_obj._p_value))
                if affinity > 10.0:
                    affinity = 10.0            
                row.append(affinity)
            adjacency_matrix.append(row)
        n_clusters = 1
        # Make a plot of the adjacency matrix #
        merged_clusters_done = False
        while merged_clusters_done == False:
            try:
                n_clusters += 1
                sc = SpectralClustering(n_clusters, affinity='precomputed', n_init=100)
                sc.fit(adjacency_matrix)
                cluster_labels = sc.labels_
                merged_clusters = []
                for i in range(0, max(cluster_labels)+1):
                    complete_cluster = []
                    for j in range(0, len(motif_list)):
                        if cluster_labels[j] == i:
                            complete_cluster.append(motif_list[j])
                    merged_clusters.append(complete_cluster)
                merged_clusters_done = check_merged_clusters(merged_clusters)
            except:
                print("\n\n\n\nError handling the adjacency matrix!!!!!!\n\n\n\n")
                merged_clusters_done = True

        final_clusters = final_clusters + merged_clusters

    return final_clusters


def equal_binding_site(motif_groups):

    final_groups = []
    # Iterate each group #
    for group in motif_groups:
        current_groups = {}
        # Iterate each motif per group #
        for motif in group:
            # cluster motifs by binding site length #
            bs_length = motif["binding_site_length"]
            if not bs_length in current_groups.keys():
                current_groups[bs_length] = [motif]
            else:
                current_groups[bs_length].append(motif)
        for bs_length in current_groups.keys():
            final_groups.append(current_groups[bs_length])

    return final_groups



def sort_motif_clusters_by_IC(motif_groups):

    final_groups = []
    list_of_groups = []
    # Create a list with motifs and their minimum information content #
    for group in motif_groups:
        IC = None
        for motif in group:
            if IC == None:
                IC = motif["IC"]
            elif IC < motif["IC"]:
                IC = motif["IC"]
        list_of_groups.append({"IC": IC, "group": group})
    # Sort the list according to the information content #
    sorted_list = sorted(list_of_groups, key=lambda k: k["IC"], reverse=True)
    for mot in sorted_list:
        final_groups.append(mot["group"])

    return final_groups



def get_motif_clusters(motif_list, pval_cutoffs, equal_bs=True):

    if len(motif_list) == 1:
        return [motif_list]

    else:
        # Group binding motifs by domains and templates #
        domains = set()
        templates = set()
        for mot in motif_list:
            domains.add(mot["domain"])
            templates.add(mot["template"])
        if (len(list(domains)) == 1) and (len(list(templates)) == 1):
            motif_clusters = []
            motif_clusters.append(motif_list)
            return motif_clusters

        else:
            # If we get here it means that we have to deal with some serious clustering #
            # We should get an adjacency matrix for the PWMs of the different motifs #
            motif_clusters = []
            if len(motif_list) > 2:
                motif_groups = cluster_motifs(motif_list)
            else:
                if get_affinity_between_clusters([motif_list[0]], [motif_list[1]]) >= 2:
                    motif_groups = [motif_list]
                else:
                    motif_groups = [[motif_list[0]], [motif_list[1]]]
            # Make sure that in each cluster we have only motifs with the same binding site length #
            if equal_bs == True:
            	motif_groups = equal_binding_site(motif_groups)
            # Sort motif clusters by information content #
            motif_groups = sort_motif_clusters_by_IC(motif_groups)
            # Renumber the ID for each of the binding motifs #
            for i in range(0, len(motif_groups)):
                for m in range(0, len(motif_groups[i])):
                    motif_groups[i][m]["ID"] = m+1

            return motif_groups
