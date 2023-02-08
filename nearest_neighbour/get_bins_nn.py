import sys,os
import pandas as pd
import numpy as np

folder = sys.argv[1]
output = sys.argv[2]

bins = [15,25,35,45,55,65,75,85,95]

table={}
table_nn={}
table_mdl={}
families=set()
for fi in os.listdir(folder):
    for b in bins:
        bin_nn = str(b)+".dat"
        bin_mdl= str(b)+"M.dat"
        if fi.endswith(bin_nn):
            family = fi.strip(bin_nn).split("graph_")[1].replace("_"," ")
            families.add(family)
            with open(os.path.join(folder,fi),"r") as fo:
                tfs=set()
                pred=0
                for line in fo:
                    data = line.split()
                    tf   = data[0].split("+")[0]
                    pred = pred +1
                    tfs.add(tf)
                table_nn.setdefault((family,b),(len(tfs),pred))
            fo.close()
        if fi.endswith(bin_mdl):
            family = fi.strip(bin_mdl).split("graph_")[1].replace("_"," ")
            families.add(family)
            with open(os.path.join(folder,fi),"r") as fo:
                tfs=set()
                pred=0
                for line in fo:
                    data = line.split()
                    tf   = data[0].split("+")[0]
                    pred = pred +1
                    tfs.add(tf)
                table_mdl.setdefault((family,b),(len(tfs),pred))
            fo.close()

for family in families:
    table.setdefault("Family",[]).append(family)
    for  b in bins:
        if table_nn.has_key((family,b)):
            num_tf,num_pred = table_nn[(family,b)]
            table.setdefault("tfs_nn_"+str(b),[]).append(num_tf)
            table.setdefault("predictions_nn_"+str(b),[]).append(num_pred)
        else:
            table.setdefault("tfs_nn_"+str(b),[]).append(0)
            table.setdefault("predictions_nn_"+str(b),[]).append(0)
            
        if table_mdl.has_key((family,b)):
            num_tf,num_pred = table_mdl[(family,b)]
            table.setdefault("tfs_modcre_"+str(b),[]).append(num_tf)
            table.setdefault("predictions_modcre_"+str(b),[]).append(num_pred)
        else:
            table.setdefault("tfs_modcre_"+str(b),[]).append(0)
            table.setdefault("predictions_modcre_"+str(b),[]).append(0)

df=pd.DataFrame(table)
col = []
col.append("Family")
for b in bins:
    col.append("tfs_nn_"+str(b))
    col.append("predictions_nn_"+str(b))
    col.append("tfs_modcre_"+str(b))
    col.append("predictions_modcre_"+str(b))
do = df[col]
do.to_csv(output+".csv")

