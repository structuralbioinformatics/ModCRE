
                        ##### Library import #####

import argparse # Parsing of arguments.
import os, sys # Main system libraries.
import shutil
import json
import warnings
import difflib
import re

from Bio.PDB import PDBIO, Superimposer, PDBParser

        ##### Initialization of Tool objects and information. #####

parser = PDBParser(PERMISSIVE=1)
sup = Superimposer()
io = PDBIO()

incompatible ={}
unreferenced =set()
wrong =set()
alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWZ1234567890"
alphabet2= []
paired_chains=set()

#for x in alphabet:
#    alphabet2.append(x)
for x in alphabet:
    for y in alphabet:
       alphabet2.append(x+y)
for x in alphabet:
    for y in alphabet:
       for z in alphabet:
           alphabet2.append(x+y+z)
aminoacids = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', "DT":"T", "DG":"G", "DC":"C",
     "DA":"A", "DU":"U", "U":"U", "A":"A", "T":"T", "G":"G", "C":"C", "UNK":"X", "X":"X"}

                            ##### Function zone #####
def make_paired_chains(namemap,exact_dimers):
      for p in exact_dimers.keys():
         pairing = False
         for name in namemap.keys():
             if p in namemap[name]: 
                 chain_1 = name
                 pairing=True
         if not pairing: continue
         for pp in exact_dimers[p]:
           for name in namemap.keys():
             if pp in namemap[name]: 
                paired_chains.add((chain_1,name))
      return paired_chains

def make_dimers_list(index,namemap,filelist,added_proteins,selected,conformers,dimers,args):

      for file in os.listdir(args["folder"]):

        if file.startswith("dna"): continue 
        if not file.endswith(".pdb"): continue
        m=re.search("(\S+).(\S+)_(\S+).(\d+)-(\d+):(\d+):(\d+)_(\S+).pdb",file)
        if m:
            start_model = int(m.group(6))
            end_model   = int(m.group(7))
            prot        = m.group(1)
            pdb         = m.group(2)
            chain       = m.group(3)
            binding     = "%d-%d"%(int(m.group(4)),int(m.group(5)))
            pdb_chain   = pdb+"_"+chain
            code_argument = prot+"."+pdb+"_"+chain+"."+binding
            conformers.setdefault((prot,pdb,binding),set()).add(chain)
            modelsize.setdefault((prot,pdb,binding),[]).append(end_model-start_model)
        else:
            start_model = 0
            end_model   = 0
            continue

        #code_argument = file.split(":")[0]
        #codes         = code_argument.split(".")

        if code_argument not in dimers.keys(): continue

        filepath = os.path.join(args["folder"], file)

        #if len(codes) == 3:
        #   prot,pdb_chain,binding = codes
        #   pdb,chain   = pdb_chain.split("_")
        #else:
        #   continue
        print("\tAnalyzing file ... "+filepath)
        if (prot, binding) not in selected.keys():continue
        if pdb == selected[(prot, binding)]:
            print("\t\t TEMPLATE "+pdb+ " CHAINS %s"%conformers[(prot,pdb,binding)])
            use=False
            #if file[6] == ".":
            if m:
            #WARNING: DOBLE CHECK the PDB files should always have a UNiprotAccession code of 6 characters followed by a '.'
                if prot+"."+pdb_chain+"."+binding in added_proteins:continue
                if end is not None:
                   fragment=binding
                   positions=fragment.split("-")
                   if int(positions[0]) >= start and int(positions[1])<= end:
                     added_proteins.append(prot+"."+pdb_chain+"."+binding)
                     use=True
                else:      
                     added_proteins.append(prot+"."+pdb_chain+"."+binding)
                     use=True
            else:
                continue
                #if file.split("_")[0] in added_proteins:continue
                #added_proteins.append(file.split("_")[0])
            if use:
              namemap[alphabet2[index]] = filepath # This is for knowing the original file.
              filelist.append(filepath)
              print("\tSelect to use ... "+filepath)
              index += 1
              if index == len(alphabet2): index = 0
              #check dimers
              if code_argument in dimers.keys():
                 paired_codes = [x for x in dimers[code_argument] if pdb in x]
                 print("\t\t-- Dimmerize with %s "%str(paired_codes))
                 done=False
                 done_hetero=False
                 done_homo=False
                 for code_pair in paired_codes:
                   try:
                     codes_pair    = code_pair.split(".")
                     if len(codes_pair) >= 3:
                        prot_pair=".".join(codes_pair[:-2])
                        pdb_chain_pair = codes_pair[-2] 
                        binding_pair   = codes_pair[-1]
                        pdb_pair,chain_pair   = pdb_chain_pair.split("_")
                     else:
                        continue
                     if prot_pair!=prot and done_hetero:continue
                     if prot_pair==prot and done_homo:continue
                     if (prot_pair, binding_pair) not in selected.keys():continue
                     if pdb_pair ==  selected[(prot_pair, binding_pair)]:
                        if prot_pair+"."+pdb_chain_pair+"."+binding_pair in added_proteins: continue
                        if end is not None:
                           positions_pair=binding_pair.split("-")
                           if int(positions_pair[0]) >= start and int(positions_pair[1])<= end:
                              added_proteins.append(prot_pair+"."+pdb_chain_pair+"."+binding_pair)
                              for file_pair in [x for x in os.listdir(args["folder"]) if x.startswith(prot_pair+"."+pdb_chain_pair+"."+binding_pair) and x.endswith(".pdb")]:
                                  namemap[alphabet2[index]] = os.path.join(args["folder"], file_pair)
                                  print("\t\t-- Add dimer SELECTED "+file_pair)
                                  filelist.append(os.path.join(args["folder"], file_pair))
                                  index += 1
                                  if index == len(alphabet2): index = 0
                              if prot_pair!=prot:done_hetero=True
                              if prot_pair==prot:done_homo=True
                        else:      
                              for file_pair in [x for x in os.listdir(args["folder"]) if x.startswith(prot_pair+"."+pdb_chain_pair+"."+binding_pair) and x.endswith(".pdb")]:
                                  namemap[alphabet2[index]] = os.path.join(args["folder"], file_pair)
                                  print("\t\t-- Add dimer SELECTED "+file_pair)
                                  filelist.append(os.path.join(args["folder"], file_pair))
                                  index += 1
                                  if index == len(alphabet2): index = 0
                              added_proteins.append(prot_pair+"."+pdb_chain_pair+"."+binding_pair)
                              if prot_pair!=prot:done_hetero=True
                              if prot_pair==prot:done_homo=True
                   except Exception as e:
                     print("\t\t-- Skip dimer check, error: "+str(e))
                 if done_hetero and done_homo: done=True
                 if not done:
                    paired_codes = [x for x in dimers[code_argument] if pdb in x]
                    for code_pair in paired_codes:
                     if done: continue
                     codes_pair    = code_pair.split(".")
                     if len(codes_pair) >= 3:
                        prot_pair=".".join(codes_pair[:-2])
                        pdb_chain_pair = codes_pair[-2] 
                        binding_pair   = codes_pair[-1]
                        pdb_pair,chain_pair   = pdb_chain_pair.split("_")
                     else:
                        continue
                     if prot_pair!=prot and done_hetero:continue
                     if prot_pair==prot and done_homo:continue
                     if prot_pair+"."+pdb_chain_pair+"."+binding_pair in added_proteins: continue
                     if end is not None:
                           positions_pair=binding_pair.split("-")
                           if int(positions_pair[0]) >= start and int(positions_pair[1])<= end:
                              try:
                                added_proteins.append(prot_pair+"."+pdb_chain_pair+"."+binding_pair)
                                file_pair = [x for x in os.listdir(args["folder"]) if x.startswith(prot_pair+"."+pdb_chain_pair+"."+binding_pair) and x.endswith(".pdb") ][0]
                                namemap[alphabet2[index]] = os.path.join(args["folder"], file_pair)
                                print("\t\t-- Add dimer previously NON-SELECTED "+file_pair)
                                filelist.append(os.path.join(args["folder"], file_pair))
                                index += 1
                                if index == len(alphabet2): index = 0
                                if prot_pair!=prot:done_hetero=True
                                if prot_pair==prot:done_homo=True
                              except Exception as e:
                                print("\t\t--skip "+str(code_pair))
                                print("\t\t\t--error "+str(e))
                     else:      
                              try:
                                added_proteins.append(prot_pair+"."+pdb_chain_pair+"."+binding_pair)
                                file_pair = [x for x in os.listdir(args["folder"]) if x.startswith(prot_pair+"."+pdb_chain_pair+"."+binding_pair) and x.endswith(".pdb") ][0]
                                namemap[alphabet2[index]] = os.path.join(args["folder"], file_pair)
                                print("\t\t-- Add dimer previously NON-SELECTED "+file_pair)
                                filelist.append(os.path.join(args["folder"], file_pair))
                                index += 1
                                if index == len(alphabet2): index = 0
                                if prot_pair!=prot:done_hetero=True
                                if prot_pair==prot:done_homo=True
                              except Exception as e:
                                print("\t\t--skip "+str(code_pair))
                                print("\t\t\t--error "+str(e))

      return index,namemap,filelist,added_proteins

def make_monomers_list(index,namemap,filelist,added_proteins,selected,conformers,dimers,args):

      for file in os.listdir(args["folder"]):

        if file.startswith("dna"): continue 
        if not file.endswith(".pdb"): continue
       
        m=re.search("(\S+).(\S+)_(\S+).(\d+)-(\d+):(\d+):(\d+)_(\S+).pdb",file)
        if m:
            start_model = int(m.group(6))
            end_model   = int(m.group(7))
            prot        = m.group(1)
            pdb         = m.group(2)
            chain       = m.group(3)
            binding     = "%d-%d"%(int(m.group(4)),int(m.group(5)))
            pdb_chain   = pdb+"_"+chain
            code_argument = prot+"."+pdb+"_"+chain+"."+binding
            conformers.setdefault((prot,pdb,binding),set()).add(chain)
            modelsize.setdefault((prot,pdb,binding),[]).append(end_model-start_model)
        else:
            start_model = 0
            end_model   = 0
            continue

        #code_argument = file.split(":")[0]
        #codes         = code_argument.split(".")
        if code_argument in dimers.keys(): continue

        filepath = os.path.join(args["folder"], file)

        #if len(codes) == 3:
        #   prot,pdb_chain,binding = codes
        #   pdb,chain   = pdb_chain.split("_")
        #else:
        #   continue
        print("\tAnalyzing file ... "+filepath)
        if (prot, binding) not in selected.keys():continue
        if pdb == selected[(prot, binding)]:
            print("\t\t TEMPLATE "+pdb+ " CHAINS %s"%conformers[(prot,pdb,binding)])
            use=False
            #if file[6] == ".":
            if m:
            #WARNING: DOBLE CHECK the PDB files should always have a UNiprotAccession code of 6 characters followed by a '.'
                if prot+"."+pdb_chain+"."+binding in added_proteins:continue
                if end is not None:
                   fragment=binding
                   positions=fragment.split("-")
                   if int(positions[0]) >= start and int(positions[1])<= end:
                     added_proteins.append(prot+"."+pdb_chain+"."+binding)
                     use=True
                else:      
                     added_proteins.append(prot+"."+pdb_chain+"."+binding)
                     use=True
            else:
                continue
                #if file.split("_")[0] in added_proteins:continue
                #added_proteins.append(file.split("_")[0])
            if use:
              namemap[alphabet2[index]] = filepath # This is for knowing the original file.
              filelist.append(filepath)
              print("\tSelect to use ... "+filepath)
              index += 1
              if index == len(alphabet2): index = 0
              #check dimers
              if code_argument in dimers.keys():
                 paired_codes = [x for x in dimers[code_argument] if pdb in x]
                 print("\t\t-- Dimmerize with %s "%str(paired_codes))
                 done=False
                 done_hetero=False
                 done_homo=False
                 for code_pair in paired_codes:
                     codes_pair    = code_pair.split(".")
                     if len(codes_pair) >= 3:
                        prot_pair=".".join(codes_pair[:-2])
                        pdb_chain_pair = codes_pair[-2] 
                        binding_pair   = codes_pair[-1]
                        pdb_pair,chain_pair   = pdb_chain_pair.split("_")
                     else:
                        continue
                     if prot_pair!=prot and done_hetero:continue
                     if prot_pair==prot and done_homo:continue
                     if (prot_pair, binding_pair) not in selected.keys():continue
                     if pdb_pair ==  selected[(prot_pair, binding_pair)]:
                        if prot_pair+"."+pdb_chain_pair+"."+binding_pair in added_proteins: continue
                        if end is not None:
                           positions_pair=binding_pair.split("-")
                           if int(positions_pair[0]) >= start and int(positions_pair[1])<= end:
                            try:
                              added_proteins.append(prot_pair+"."+pdb_chain_pair+"."+binding_pair)
                              for file_pair in [x for x in os.listdir(args["folder"]) if x.startswith(prot_pair+"."+pdb_chain_pair+"."+binding_pair) and x.endswith(".pdb")]:
                                  namemap[alphabet2[index]] = os.path.join(args["folder"], file_pair)
                                  print("\t\t-- Add dimer SELECTED "+file_pair)
                                  filelist.append(os.path.join(args["folder"], file_pair))
                                  index += 1
                                  if index == len(alphabet2): index = 0
                              if prot_pair!=prot:done_hetero=True
                              if prot_pair==prot:done_homo=True
                            except Exception as e:
                                print("\t\t--skip "+str(code_pair))
                                print("\t\t\t--error "+str(e))
                        else:      
                            try:
                              for file_pair in [x for x in os.listdir(args["folder"]) if x.startswith(prot_pair+"."+pdb_chain_pair+"."+binding_pair) and x.endswith(".pdb")]:
                                  namemap[alphabet2[index]] = os.path.join(args["folder"], file_pair)
                                  print("\t\t-- Add dimer SELECTED "+file_pair)
                                  filelist.append(os.path.join(args["folder"], file_pair))
                                  index += 1
                                  if index == len(alphabet2): index = 0
                              added_proteins.append(prot_pair+"."+pdb_chain_pair+"."+binding_pair)
                              if prot_pair!=prot:done_hetero=True
                              if prot_pair==prot:done_homo=True
                            except Exception as e:
                                print("\t\t--skip "+str(code_pair))
                                print("\t\t\t--error "+str(e))
                 if done_hetero and done_homo: done=True
                 if not done:
                    paired_codes = [x for x in dimers[code_argument] if pdb in x]
                    for code_pair in paired_codes:
                     if done: continue
                     codes_pair    = code_pair.split(".")
                     if len(codes_pair) >= 3:
                        prot_pair=".".join(codes_pair[:-2])
                        pdb_chain_pair = codes_pair[-2] 
                        binding_pair   = codes_pair[-1]
                        pdb_pair,chain_pair   = pdb_chain_pair.split("_")
                     else:
                        continue
                     if prot_pair!=prot and done_hetero:continue
                     if prot_pair==prot and done_homo:continue
                     if prot_pair+"."+pdb_chain_pair+"."+binding_pair in added_proteins: continue
                     if end is not None:
                           positions_pair=binding_pair.split("-")
                           if int(positions_pair[0]) >= start and int(positions_pair[1])<= end:
                            try:
                              added_proteins.append(prot_pair+"."+pdb_chain_pair+"."+binding_pair)
                              file_pair = [x for x in os.listdir(args["folder"]) if x.startswith(prot_pair+"."+pdb_chain_pair+"."+binding_pair) and x.endswith(".pdb") ][0]
                              namemap[alphabet2[index]] = os.path.join(args["folder"], file_pair)
                              print("\t\t-- Add dimer previously NON-SELECTED "+file_pair)
                              filelist.append(os.path.join(args["folder"], file_pair))
                              index += 1
                              if index == len(alphabet2): index = 0
                              if prot_pair!=prot:done_hetero=True
                              if prot_pair==prot:done_homo=True
                            except Exception as e:
                                print("\t\t--skip "+str(code_pair))
                                print("\t\t\t--error "+str(e))
                     else:      
                            try:
                              added_proteins.append(prot_pair+"."+pdb_chain_pair+"."+binding_pair)
                              file_pair = [x for x in os.listdir(args["folder"]) if x.startswith(prot_pair+"."+pdb_chain_pair+"."+binding_pair) and x.endswith(".pdb") ][0]
                              namemap[alphabet2[index]] = os.path.join(args["folder"], file_pair)
                              print("\t\t-- Add dimer previously NON-SELECTED "+file_pair)
                              filelist.append(os.path.join(args["folder"], file_pair))
                              index += 1
                              if index == len(alphabet2): index = 0
                              if prot_pair!=prot:done_hetero=True
                              if prot_pair==prot:done_homo=True
                            except Exception as e:
                                print("\t\t--skip "+str(code_pair))
                                print("\t\t\t--error "+str(e))

      return index,namemap,filelist,added_proteins



def extract_dimers(args):
      potential_dimers={}
      predimers={}
      exact_dimers={}
      orthologs_file = args["orthologs"]
      if os.path.exists(orthologs_file):
          orthologs_list = json.loads(''.join([line for line in parse_file(orthologs_file)]))
      for cluster in orthologs_list:
        for thread in cluster["monomer"]:
            thread_file, score, d_score = thread
            monomers.add(thread_file.lstrip("aux_files/").rstrip(".txt"))
        for thread_dimer in cluster["dimer"]:
            thread_file, score, d_score = thread_dimer[0]
            protein_1=thread_file.lstrip("aux_files/").rstrip(".txt")
            thread_file, score, d_score = thread_dimer[1]
            protein_2=thread_file.lstrip("aux_files/").rstrip(".txt")
            potential_dimers.setdefault(protein_1,set()).add(protein_2)
            potential_dimers.setdefault(protein_2,set()).add(protein_1)
      print("POTENTIAL DIMERS OBTAINED FROM %s "%args["orthologs"])
      for p in potential_dimers.keys():
          print("\t-- %s dimerize with %s "%(p,str(potential_dimers[p])))

      for p in potential_dimers.keys():
        codes = p.split(".")
        if len(codes) >= 3:
              prot=".".join(codes_pair[:-2])
              pdb_chain = codes_pair[-2] 
              binding   = codes_pair[-1]
              pdb,chain   = pdb_chain.split("_")
              start_binding, end_binding = binding.split("-")
              nucleotides_binding = set([int(x) for x in range(int(start_binding),int(end_binding))])
        partners={}
        for partner in potential_dimers[p]:
            codes_pair    = partner.split(".")
            if len(codes_pair) >= 3:
                        prot_pair=".".join(codes_pair[:-2])
                        pdb_chain_pair = codes_pair[-2] 
                        binding_pair   = codes_pair[-1]
                        start_binding_pair, end_binding_pair = binding_pair.split("-")
                        nucleotides_binding_pair = set([int(x) for x in range(int(start_binding_pair),int(end_binding_pair))])
            else:
                continue
            if pdb_pair != pdb: continue
            if chain_pair == chain: continue
            if partners.get(prot_pair) is not None:
               if len(nucleotides_binding_pair.intersection(nucleotides_binding)) > partners[prot_pair][1]:
                  partners[prot_pair] = (partner,len(nucleotides_binding_pair.intersection(nucleotides_binding)))
            else:
               partners.setdefault(prot_pair, (partner,len(nucleotides_binding_pair.intersection(nucleotides_binding))))
        for prot_pair in partners.keys():
               predimers.setdefault(p,set()).add(partners[prot_pair][0])
      #Check true dimers to distinguish from others
      for p in predimers.keys():
          codes = p.split(".")
          prot=".".join(codes_pair[:-2])
          pdb_chain = codes_pair[-2] 
          binding   = codes_pair[-1]
          true_dimer=True
          for pp in predimers[p]:
              pcodes = pp.split(".")
              pprot=".".join(pcodes_pair[:-2])
              ppdb_chain = pcodes_pair[-2] 
              pbinding   = pcodes_pair[-1]
              if pbinding != binding: true_dimer=False
          if true_dimer:
              for ppp in predimers[p]:
                    exact_dimers.setdefault(p,set()).add(ppp)
                    exact_dimers.setdefault(ppp,set()).add(p)
      for p in predimers.keys():
          if p in exact_dimers.keys():
              print("\t-- %s is a True dimer %s"%(p,str(exact_dimers[p])))
              dimers[p]=exact_dimers[p]
          else:
             if len(predimers[p].intersection(set([x for x in exact_dimers.keys()])))<=0:
                print("\t-- %s is a shifted dimer %s"%(p,str(predimers[p])))
                dimers[p]=predimers[p]

      #Print DIMERS
      print("DIMERS OBTAINED FROM %s "%args["orthologs"])
      for p in dimers.keys():
          print("\t-- %s dimerize with %s "%(p,str(dimers[p])))
      return exact_dimers,dimers


def create_dimers(selected,conformers,args):
         potential_dimers={}
         dimers={}
         predimers={}
         exact_dimers={}
         proteins = set([prot for prot,binding in selected.keys()])
         for prot,binding in selected.keys():
             pdb    = selected[(prot,binding)]
             chains = conformers[(prot,pdb,binding)]
             start_binding, end_binding = binding.split("-")
             nucleotides_binding = set([int(x) for x in range(int(start_binding),int(end_binding))])
             for chain in chains:
                 pdb_chain = pdb+"_"+chain
                 protein_1 = prot+"."+pdb_chain+"."+binding
                 for prot_test in proteins:
                     for chain_pair in chains:
                         if chain_pair == chain: continue
                         pdb_chain_pair = pdb+"_"+chain_pair
                         pair_id        = prot_test+"."+pdb_chain_pair+"."+binding
                         accept = [x for x in os.listdir(args["folder"]) if x.startswith(pair_id)]
                         if len(accept)>0:
                            protein_2   = pair_id
                            potential_dimers.setdefault(protein_1,set()).add(protein_2)
                            potential_dimers.setdefault(protein_2,set()).add(protein_1)
                 for file in os.listdir(args["folder"]):
                     if file.startswith("dna"): continue 
                     if not file.endswith(".pdb"): continue
                     for prot_test in proteins:
                         for chain_pair in chains:
                             pdb_chain_pair = pdb+"_"+chain_pair
                             pair_id        = prot_test+"."+pdb_chain_pair+"."
                             if not file.startswith(pair_id):continue
                             if file.startswith(protein_1): continue
                             protein_2 = file.split(":")[0]
                             codes_pair= protein_2.split(".")
                             if len(codes_pair) >= 3:
                                prot_pair=".".join(codes_pair[:-2])
                                pdb_chain_pair = codes_pair[-2] 
                                binding_pair   = codes_pair[-1]
                                pdb_pair,chain_pair                   = pdb_chain_pair.split("_")
                                start_binding_pair, end_binding_pair  = binding_pair.split("-")
                                nucleotides_binding_pair = set([int(x) for x in range(int(start_binding_pair),int(end_binding_pair))])
                                if len(nucleotides_binding_pair.intersection(nucleotides_binding))>0:
                                   potential_dimers.setdefault(protein_1,set()).add(protein_2)
                                   potential_dimers.setdefault(protein_2,set()).add(protein_1)
         print("POTENTIAL DIMERS DERIVED FROM DATA")
         for p in potential_dimers.keys():
             print("\t-- %s dimerize with %s "%(p,str(potential_dimers[p])))
         for p in potential_dimers.keys():
           codes = p.split(".")
           if len(codes) >= 3:
              prot=".".join(codes_pair[:-2])
              pdb_chain = codes_pair[-2] 
              binding   = codes_pair[-1]
              pdb,chain   = pdb_chain.split("_")
              start_binding, end_binding = binding.split("-")
              nucleotides_binding = set([int(x) for x in range(int(start_binding),int(end_binding))])
           partners={}
           for partner in potential_dimers[p]:
              codes_pair    = partner.split(".")
              if len(codes_pair) >= 3:
                        prot_pair=".".join(codes_pair[:-2])
                        pdb_chain_pair = codes_pair[-2] 
                        binding_pair   = codes_pair[-1]
                        pdb_pair,chain_pair   = pdb_chain_pair.split("_")
                        start_binding_pair, end_binding_pair = binding_pair.split("-")
                        nucleotides_binding_pair = set([int(x) for x in range(int(start_binding_pair),int(end_binding_pair))])
              else:
                continue
              if pdb_pair != pdb: continue
              if chain_pair == chain: continue
              if partners.get(prot_pair) is not None:
               if len(nucleotides_binding_pair.intersection(nucleotides_binding)) > partners[prot_pair][1]:
                  partners[prot_pair] = (partner,len(nucleotides_binding_pair.intersection(nucleotides_binding)))
              else:
               partners[prot_pair] = (partner,len(nucleotides_binding_pair.intersection(nucleotides_binding)))
           for prot_pair in partners.keys():
               predimers.setdefault(p,set()).add(partners[prot_pair][0])

         #Check true dimers to distinguish from others
         for p in predimers.keys():
            prot,pdb_chain,binding = p.split(".")
            true_dimer=True
            for pp in predimers[p]:
                pcodes = pp.split(".")
                pprot=".".join(pcodes_pair[:-2])
                ppdb_chain = pcodes_pair[-2] 
                pbinding   = pcodes_pair[-1]
                if pbinding != binding: true_dimer=False
            if true_dimer:
                for ppp in predimers[p]:
                    exact_dimers.setdefault(p,set()).add(ppp)
                    exact_dimers.setdefault(ppp,set()).add(p)
         for p in predimers.keys():
            if p in exact_dimers.keys():
                print("\t-- %s is a True dimer %s"%(p,str(exact_dimers[p])))
                dimers[p]=exact_dimers[p]
            else:
             if len(predimers[p].intersection(set([x for x in exact_dimers.keys()])))<=0:
                print("\t-- %s is a shifted dimer %s"%(p,str(predimers[p])))
                dimers[p]=predimers[p]


         #Print DIMERS
         print("DIMERS DERIVED FROM DATA")
         for p in dimers.keys():
             print("\t-- %s dimerize with %s "%(p,str(dimers[p])))

         return exact_dimers,dimers

def parse_arguments():  # This function passes the arguments given by the user.

    """
    Parsing of arguments.
    -d : Path to the folder with the PDB's Pairwise Interactions that conforms the Complex.
    -output: Path to the folder where generate the output file complex.
    -v: Activate the verbosity mode
    """

    aparser = argparse.ArgumentParser(description='complexbuilder -d FOLDER -o OUTPUT_FOLDER [-maxit MAXIT]')

    aparser.add_argument('-d', action="store", required =True, dest="folder", help="Path to the folder with the PDB's Pairwise Interactions that conforms the Complex.")
    aparser.add_argument('-o', action="store", required=True, dest="output_folder", help="Path to the folder where generate the output file complex.")
    aparser.add_argument('-maxit', action="store", default=5040, dest="maxit", help="Number of iterations to stop the execution. (Default is 7!)")
    aparser.add_argument('-e', action="store", default=2, dest="marginal_error", help="Marginal error to match the DNA sequence position by threading. (Default is 2!)")
    aparser.add_argument('-ortho', action="store", default=None, dest="orthologs", help="JSON file of orthologs list (i.e. with best templates) to identify dimers (default is None)")

    args = vars(aparser.parse_args())
    return args

def parse_file(file_name, gz=False):
    """
    This function parses any file and yields lines one by one.
    
    @input:
    file_name {string}
    @return:
    line {string}

    """
    import gzip
    if os.path.exists(file_name):
        # Initialize #
        f = None
        # Open file handle #
        if gz: f = gzip.open(file_name, "rt")
        else: f = open(file_name, "rt")
        # For each line... #
        for line in f:
            yield line.strip("\n")
        f.close()
    else:
        raise ValueError("Could not open file %s" % file_name)


def Get_fasta(chain): # This function converts an input BioPython chain of a protein into a fasta sequence (str)

    """
    This function converts an input BioPython chain of a protein into a fasta sequence (str)
    """
    sequence = ""
    for residue in chain:
        try:
            sequence += aminoacids[residue.get_resname().strip()]
        except:
            continue
    return sequence

def check_homology(fasta_1, fasta_2): # Check Homology between two fasta chains. Threshold = 55%
   
    """
    This function check the identity percentage between two chains, and returns True if they are homologous
    or False if they are not. The threshold is currently at 70%.

    @ Input - Fasta_1 : String of fasta sequence
              Fasta_2 : String of fasta sequence

    @ Output - Boolean variable indicating if sequence are homologous
    """

    if fasta_1 in fasta_2 or fasta_2 in fasta_1:
        return True

    else:
        return difflib.SequenceMatcher(None, fasta_1, fasta_2).ratio() > 0.55


def return_equal_chain(structure_1, structure_2, positions,marginal):

    testing_positions=[]
    for i in range(-marginal,marginal):
      for j in range(-marginal,marginal):
            testing_positions.append((int(positions[0])+i,int(positions[1])+j))

    for chain in structure_1: 
        for chain_2 in structure_2[0]: 
            if is_dna(chain) and is_dna(chain_2): # If one chain is DNA, use the positions.

                residues_id = [residue.get_id()[1] for residue in chain.get_residues()]

                for a,b in testing_positions:
                
                  try:
                    starting_index = residues_id.index(a)
                  except:
                    continue
            
                  difference = b - a
                  residues = list(chain.get_residues())[starting_index: starting_index + difference + 1]

                  fasta1 = Get_fasta(residues)
                  fasta2 = Get_fasta(chain_2)

                  print("DNA sequence 1:"+str(fasta1))
                  print("DNA sequence 2:"+str(fasta2))

                  if fasta1 in fasta2 or fasta2 in fasta1:

                    s = difflib.SequenceMatcher(None, fasta1, fasta2)
                    matching_blocks = s.get_matching_blocks()
                    if matching_blocks[0].a != 0:
                        i = matching_blocks[0].a
                    else:
                        i = matching_blocks[0].b

                    size = matching_blocks[0].size
                    new_residues = residues[i: i + size]

                    atoms = list()
                    for residue in new_residues:
                        atoms = atoms + list(residue.get_atoms())
                    
                    main_atoms_1 = [atom for atom in atoms if (atom.get_id() == "P" or atom.get_id() == "CA")]
                    main_atoms_2 = [atom for atom in chain_2.get_atoms() if (atom.get_id() == "P" or atom.get_id() == "CA")]

                    if len(main_atoms_1) != len(main_atoms_2):
                        shorter_len = min(len(main_atoms_1), len(main_atoms_2))
                        return (main_atoms_1[:shorter_len],  main_atoms_2[:shorter_len])

                    return main_atoms_1 , main_atoms_2
                
                  else:
                    continue
            
            else: # When both chains are proteins, use homology:

                fasta_1 = Get_fasta(chain)
                fasta_2 = Get_fasta(chain_2)

                if check_homology(Get_fasta(chain), Get_fasta(chain_2)):

                    s = difflib.SequenceMatcher(None, fasta_1, fasta_2)
                    matching_blocks = s.get_matching_blocks()
                    if matching_blocks[0].a != 0:
                        i = matching_blocks[0].a
                    else:
                        i = matching_blocks[0].b

                    size = matching_blocks[0].size
                    residues = list(chain.get_residues())[i: i + size]

                    atoms = list()
                    for residue in residues:
                        atoms = atoms + list(residue.get_atoms())

                    main_atoms_1 = [atom for atom in atoms if atom.get_id() == "CA"]
                    main_atoms_2 = [atom for atom in chain_2.get_atoms() if atom.get_id() == "CA"]


                    if len(main_atoms_1) != len(main_atoms_2):
                        shorter_len = min(len(main_atoms_1), len(main_atoms_2))
                        return (main_atoms_1[:shorter_len],  main_atoms_2[:shorter_len])

                    return main_atoms_1 , main_atoms_2

    return (None, None)

def check_clash(structure_complex, moved_chain,chain_2):

    from Bio.PDB.NeighborSearch import NeighborSearch
    moved_chain.get_parent().get_parent().id="MOBILE"

    #define trace of CA or backbone
    backbone_atoms = ["CA"]
    #backbone_atoms = ["CA", "P"]
    #backbone_atoms = ["CA","CB","C","N","O","P"]
    #backbone_atoms = ["CA","CB","C","P"]
    moving_bb_atoms = [atom for atom in moved_chain.get_atoms() if atom.id in backbone_atoms]

    id_chains=[chain.id for chain in structure_complex]
    contact_chains =[]

    clash_num = 0
    bb_clash_num=0
    structure_residues=0
    bb_clash_percentage_structure=[]
    clash_percentage_structure=[]
    bb_clash_percentage_model=[]
    clash_percentage_model=[]
    bb_clash_percentage_structure.append(0)
    clash_percentage_structure.append(0)
    bb_clash_percentage_model.append(0)
    clash_percentage_model.append(0)

    print("\t\t--CLASHES BY CHAINS:")
    if len(structure_complex) == 2:
        return False
    else:
        for chain in structure_complex:
          moved_contacts=set()
          fixed_contacts=set()
          bb_moved_contacts=set()
          bb_fixed_contacts=set()
          chain.get_parent().get_parent().id="FIX"
          structure_residues += len([x for x in chain.get_residues()])
          structure_bb_atoms = [atom for atom in chain.get_atoms() if atom.id in backbone_atoms]
          list_of_atoms= [x for x in chain.get_atoms()]
          list_of_atoms.extend([x for x in moved_chain.get_atoms()])
          neighbors = NeighborSearch(list_of_atoms,bucket_size=10)
          contacts  = neighbors.search_all(0.8,level='R')
          if chain.id=="X":continue
          if chain.id=="Y":continue
          for res_1,res_2 in contacts:
              rchain_1 = res_1.get_parent()
              rchain_2 = res_2.get_parent()
              if rchain_1.get_parent().get_parent().id==rchain_2.get_parent().get_parent().id:continue
              if (rchain_1.id == chain.id and rchain_1.get_parent().get_parent().id=="FIX") or (rchain_2.id==chain.id and rchain_2.get_parent().get_parent().id=="FIX"):
                if rchain_1.id == chain.id and rchain_1.get_parent().get_parent().id=="FIX":
                    moved_contacts.add(res_2)
                    fixed_contacts.add(res_1)
                if rchain_2.id == chain.id and rchain_2.get_parent().get_parent().id=="FIX":
                    moved_contacts.add(res_1)
                    fixed_contacts.add(res_2)
                contact_chains.append((chain.id,chain_2))
                contact_chains.append((chain_2,chain.id))
                #print("Residue %s %s %s vs Residue %s %s %s "%(str(res_1.id),res_1.get_parent().id,res_1.get_parent().get_parent().get_parent().id,str(res_2.id),res_2.get_parent().id,res_2.get_parent().get_parent().get_parent().id))
          clash_num += len(moved_contacts)
          clash_percentage_structure.append( float(len(fixed_contacts)) / structure_residues * 100)
          clash_percentage_model.append( float(len(moved_contacts)) / len([x for x in moved_chain.get_residues()]) * 100)
          #print("\t\t\t--clashes with chain _%s_ percentage (fix %f mobile %f) total (fix %f mobile %f)"%(chain.id,clash_percentage_structure,clash_percentage_model,float(len(fixed_contacts)),float(len(moved_contacts)) ))
          #clash_percentage_structure = max(float(len(fixed_contacts)) / structure_residues * 100,clash_percentage_structure)
          #clash_percentage_model = max(float(len(moved_contacts)) / len([x for x in moved_chain.get_residues()]) * 100,clash_percentage_model)

          if  len(moved_contacts) > 1 or len(fixed_contacts) > 1 :
          #double check
            for  atom in structure_bb_atoms:
              for atom_2 in moving_bb_atoms:
                distance = (atom_2 - atom)
                if distance < 2.5:
                    bb_moved_contacts.add(atom)
                    bb_fixed_contacts.add(atom_2)
                    contact_chains.append((chain.id,chain_2))
                    contact_chains.append((chain_2,chain.id))

            bb_clash_num += len(bb_moved_contacts)
            bb_clash_percentage_structure.append( float(len(bb_fixed_contacts)) / len(structure_bb_atoms) * 100)
            bb_clash_percentage_model.append(float(len(bb_moved_contacts)) / len(moving_bb_atoms) * 100)
            #print("\t\t\t--clashes in backbone with chain _%s_: percentage (FIX %f MOBIL %f) total (fix %f mobile %f)"%(chain.id,bb_clash_percentage_structure,bb_clash_percentage_model,float(len(fixed_contacts)),float(len(moved_contacts))   ))
            #bb_clash_percentage_structure = max(float(len(bb_fixed_contacts)) / len(structure_bb_atoms) * 100,bb_clash_percentage_structure)
            #bb_clash_percentage_model = max(float(len(bb_moved_contacts)) / len(moving_bb_atoms) * 100,bb_clash_percentage_model)
          else:
            print("\t\t\t--chains %s and %s are not dimers and have no contacts"%(chain.id,chain_2))
    if (chain.id,chain_2) not in paired_chains and (chain_2,chain.id) not in paired_chains:
        print("\t\t--clashes percentage (fix %f mobile %f) total  %f"%(max(clash_percentage_structure),max(clash_percentage_model),float(clash_num) ))
        print("\t\t--clashes in backbone : percentage (FIX %f MOBIL %f) total %f)"%(max(bb_clash_percentage_structure),max(bb_clash_percentage_model),float(bb_clash_num)   ))
        if  max(bb_clash_percentage_structure) > 1 or max(bb_clash_percentage_model) > 1  or  float(bb_clash_num) > 2 :
          if  max(clash_percentage_structure) > 2 or max(clash_percentage_model) > 2 or  float(clash_num) > 10  :
            for chain_1,chain_2  in contact_chains:
               if chain_1 == "X":continue
               if chain_1 == "Y":continue
               incompatible.setdefault(chain_1,set()).add(chain_2)
               incompatible.setdefault(chain_2,set()).add(chain_1)
            return True
    else:
        print("\t\t--paired chains as dimer ")
        print("\t\t--clashes percentage (fix %f mobile %f) total  %f"%(max(clash_percentage_structure),max(clash_percentage_model),float(clash_num) ))
        print("\t\t--clashes in backbone : percentage (FIX %f MOBIL %f) total %f)"%(max(bb_clash_percentage_structure),max(bb_clash_percentage_model),float(bb_clash_num)   ))
        if  max(bb_clash_percentage_structure) > 1 or max(bb_clash_percentage_model) > 1  or  float(bb_clash_num) > 2 :
          if  max(clash_percentage_structure) > 10 or max(clash_percentage_model) > 10 or  float(clash_num) > 50  :
            for chain_1,chain_2  in contact_chains:
               if chain_1 == "X":continue
               if chain_1 == "Y":continue
               incompatible.setdefault(chain_1,set()).add(chain_2)
               incompatible.setdefault(chain_2,set()).add(chain_1)
            return True

    return False

def is_dna(chain):
    residues = chain.get_residues()
    for residue in residues:
        if residue.get_resname().strip() in ["DA", "DT", "DG", "DC", "DU"]:
            return True
        else:
            return False

def build_complex(structure_1, file_2, namemap, offset, marginal): 
    
    #Get the potential new chain 
    letter = [i for i in namemap.keys() if namemap[i] == file_2][0]

    # Based on the positions, it returns the set of equal atoms.
    warnings.filterwarnings("ignore")
    print("\t--BUILDING: check file ... "+file_2)
    structure_2 = parser.get_structure('Complex', file_2)
    try:
        m=re.search("(\S+).(\S+)_(\S+).(\d+)-(\d+):(\d+):(\d+)_(\S+).pdb",file_2)
        if m:
           positions=[]
           positions.append(m.group(4))  
           positions.append(m.group(5))  
        else:
           positions = os.path.basename(file_2).split(".")[2].split(":")[0].split("-")
    except:
        positions = os.path.basename(file_2).split("_")[3].split(":")[0].split("-")
    #Arrange for the fragment
    positions[0]=int(positions[0])-offset
    positions[1]=int(positions[1])-offset

    atoms_fixed, atoms_moving = return_equal_chain(structure_1, structure_2, positions, marginal)

    # If does not exist a set of common atoms, it returns the same structure. 
    if atoms_fixed == None or atoms_moving == None:
               if letter.isalnum() :
                   unreferenced.add(letter)
               else:
                   print("Wrong chain format in superimposing")
               print("\t\t--missing reference atoms for superimposition")
               return True, structure_1

    # Then we apply the rotation matrix. 
    try:
      sup.set_atoms(atoms_fixed, atoms_moving) # This give us the rotation matrix
      sup.apply(list(structure_2.get_atoms())) # We apply this matrix to all the structure  that we are adding.
    except:
      print("Failing to get rotation matrix")
      return True, structure_1

    # Extraction of the part of the structure we want to add (Protein and not the ones used in superposition)
    for chain in structure_2[0].get_chains():
        if chain.id != list(atoms_moving)[0].get_full_id()[2] and not is_dna(chain):
            moved_chain = chain

    # After moving it we check if exists a clash.
    try:
      if check_clash(structure_1, moved_chain,letter):
        print("\t\t--skip due to clashes")
        return True, structure_1
    except:
        if letter.isalnum() :
                   wrong.add(letter)
        else:
                   print("Wrong chain format and missing chains")
        print("\t\t--skip due to missing chains")
        return True, structure_1


    # If not, we extract the index of the dictionary that relates to the file.
    letter = [i for i in namemap.keys() if namemap[i] == file_2][0]
    print("\t--BUILDING: add as chain "+letter)

    # Finally we change the id of the chain, and add it to the structure. 
    try:
      print("\t--BUILDING: add as chain "+letter)
      moved_chain.id = letter
      structure = (structure_1) + [moved_chain]
      return False, structure
    except Exception as e:
      print("\t--BUILDING: skip chain "+letter)
      print("\t\t--skip due to "+str(e))
      return True, structure_1

def recursive_build(structure, start_chain_ids, chainlist, complexname, args, namemap,output_fragment):

    # Check the number of iterations to stop the program
    global itnum

    itnum += 1
    #print("Iteration %d (MAX %d) Size List %d List %s Complexname %s"%(int(itnum),float(args["maxit"]),len(chainlist),str(chainlist),complexname))

    #Check the origin
    fragment=output_fragment.split("_")[1]
    if fragment == "full": 
       offset=0
    else:
       offset=int(fragment.split("-")[0]) - 1


    #Get the marginal error to match the DNA sequence
    marginal=int(args["marginal_error"])

    #Check the limiting number of iterations
    if int(itnum) >= float(args["maxit"]):
        letters = list(complexname[4:-4].split("_"))
        complexname = "dna"+"_".join(sorted(letters)) + ".pdb"
        code_chains=[]
        for c in start_chain_ids:
            code_chains.append(c)
        for c in letters:
            code_chains.append(c)
        if not os.path.exists(os.path.join(args["output_folder"],output_fragment,complexname)):
          print("\t WAIT: Built last complex combination ... "+complexname)
          n_chain=0
          out_file=open( os.path.join(args["output_folder"],output_fragment,complexname) , "wt")
          letters = list(complexname[4:-4].split("_"))
          remark=open( os.path.join(args["output_folder"],output_fragment,"remarks_"+"_".join(sorted(letters)) +".txt") , "wt")
          io=PDBIO()
          for model in list(structure):
                #str_model=model.get_parent()
                #str_model.id=n_chain
                #model.id="A"
                if n_chain >len(alphabet):
                   remark.write("REMARK\tSKIP CHAIN %s"%(os.path.basename(namemap[model.id])))
                   continue
                print("CHAIN %s is now %s"%(model.id,alphabet[n_chain]))
                if model.id in start_chain_ids:
                    if model.id == "X" or model.id == "Y":
                       remark.write("REMARK\t%s = DNA\n"%(alphabet[n_chain]))
                    else:
                       remark.write("REMARK\t%s = ORIGINAL (i.e. HISTONE)\n"%(alphabet[n_chain]))
                else:
                   remark.write("REMARK\t%s = %s\n"%(alphabet[n_chain],os.path.basename(namemap[model.id])))
                #str_model.header["keywords"]=model.header["keywords"]+" "+model.id+"="+alphabet[n_chain]+" "
                mdl=model.copy()
                mdl.id=alphabet[n_chain]
                io.set_structure(mdl)
                io.save(out_file)
                n_chain = n_chain + 1
          remark.close()
        return complexname


    # BASE CASE - When there are no more chains to add in chainlist.
    # When there are no more chains to add, it stores the complex in a file.

    if chainlist == [] :

        letters = list(complexname[4:-4].split("_"))
        complexname = "dna"+"_".join(sorted(letters)) + ".pdb"
        code_chains=[]
        for c in start_chain_ids:
            code_chains.append(c)
        for c in letters:
            code_chains.append(c)
        if not os.path.exists(os.path.join(args["output_folder"],output_fragment,complexname)):
          print("\t WAIT: Intermediate complex combination (not all potential combinations are tested yet)... "+complexname)
          n_chain=0
          io=PDBIO()
          out_file=open( os.path.join(args["output_folder"],output_fragment,complexname) , "wt")
          remark=open( os.path.join(args["output_folder"],output_fragment,"remarks_"+"_".join(sorted(letters)) +".txt") , "wt")
          for model in list(structure):
                #str_model=model.get_parent()
                #str_model.id=n_chain
                mdl=model.copy()
                #model.id="A"
                if n_chain >len(alphabet):
                   remark.write("REMARK\tSKIP CHAIN %s"%(os.path.basename(namemap[model.id])))
                   continue
                print("CHAIN %s is now %s"%(model.id,alphabet[n_chain]))
                if model.id in start_chain_ids:
                    if model.id == "X" or model.id == "Y":
                       remark.write("REMARK\t%s = DNA\n"%(alphabet[n_chain]))
                    else:
                       remark.write("REMARK\t%s = ORIGINAL (i.e. HISTONE)\n"%(alphabet[n_chain]))
                else:
                   remark.write("REMARK\t%s = %s\n"%(alphabet[n_chain],os.path.basename(namemap[model.id])))
                #str_model.header["keywords"]=model.header["keywords"]+" "+model.id+"="+alphabet[n_chain]+" "
                mdl.id=alphabet[n_chain]
                io.set_structure(mdl)
                io.save(out_file)
                n_chain = n_chain + 1
          remark.close()
        #If more iterations are requested re-start the name
        return "Done"

    # RECURSIVE CASE. - When there are more chains to add in chainlist.

    else: 

        onloop=0
        for to_add in chainlist:

            #print("ON LOOP %d Complexname %s Check %s"%(onloop,complexname,to_add))
            onloop= onloop + 1

            letter = [i for i in namemap.keys() if namemap[i] == to_add][0]

            skip = False

            #Check the new chain is OK
            if letter in wrong:skip=True 
            if letter in unreferenced:skip=True
            if letter in list(complexname[4:-4].split("_")): skip = True

            #Check incompatibilities
            for chain in structure:
                chain_id = chain.id
                if incompatible.get(chain_id) is not None:
                   if letter in incompatible.get(chain_id): skip=True
                if incompatible.get(letter) is not None:
                   if chain_id in incompatible.get(letter): skip=True
    
            #Check limit of iterations
            if int(itnum) >= float(args["maxit"]):
               print("Skip file ... "+to_add)
               continue

            #proceed with recursion

            newcomplex = structure[:]
            clashed=False
            if not skip:
               clashed, newcomplex = build_complex(newcomplex, to_add, namemap,offset,marginal)
            newlist = chainlist[:]
            newlist.remove(to_add)

            if clashed or skip:
                dummy_name=recursive_build(newcomplex, start_chain_ids, newlist, complexname, args, namemap,output_fragment)
            else:
                print("\t--RECURSIVE STEP: chains on model "+complexname[4:-4]+" add new "+letter)
                dummy_name=recursive_build(newcomplex, start_chain_ids, newlist, complexname[:-4] + "_"+letter + ".pdb", args, namemap,output_fragment)


        return None

if __name__ == "__main__": # Main part of the program.

    itnum = 0
    args = parse_arguments()

    # Creation of the output folder.
    if not os.path.isdir(args["output_folder"]):

        try:
            os.makedirs(args["output_folder"], exist_ok=True)
        except OSError:
            print("Creation of the directory %s failed. Aborting..." % args["output_folder"])
            sys.exit()

    # Get the DNA fragments
    dna_files=[]
    for file in os.listdir(args["folder"]):
        if file.startswith("dna") and file.endswith("pdb"):
           datafile=file.rstrip(".pdb").split("_")
           if len(datafile)==1:
              dna_files.append(os.path.join(args["folder"],file))
              continue
           if len(datafile)==2:
              fragment=datafile[1]
              positions=fragment.split("-")
              if len(positions)==2:
                 try:
                    start=int(positions[0])
                    end  =int(positions[1])
                    if start>0 and end>0 and start<end:
                       dna_files.append(os.path.join(args["folder"],file))
                 except:
                    continue
                    
              
           
        
    if len(dna_files)<=0:
        print ("Starting DNA file named dna(*).pdb does not exist. Aborting...")
        sys.exit()

    #Read orthologs list if any
    monomers=set()
    dimers={}
    if args["orthologs"] is not None:
      exact_dimers,dimers=  extract_dimers(args)

    for startfile in dna_files:

      #start name of complex
      complexname="dna_.pdb"

      # Mapping of the files inside a dictionary with a value of the chain as index, set up first index at 0
      namemap = dict()
      filelist = list()
      added_proteins =  []
      index = 0

      #clean up the filters of compatibility and accession
      incompatible ={}
      unreferenced =set()
      wrong =set()
      itnum=0

      #Get DNA fragment
      start= 1
      end  = None
      fragment=os.path.basename(startfile).rstrip(".pdb").split("_")
      if len(fragment)==2:
              positions=fragment[1].split("-")
              if len(positions)==2:
                    start=int(positions[0])
                    end  =int(positions[1])
      if end is None:
         output_fragment="fragment_full"
      else:
         output_fragment="fragment_%d-%d"%(start,end)

      #if os.path.exists(os.path.join(args["output_folder"],output_fragment)): 
      #    os.system("rm -rf %s"% os.path.join(args["output_folder"],output_fragment))
      if not os.path.exists(os.path.join(args["output_folder"],output_fragment)):
           os.makedirs(os.path.join(args["output_folder"],output_fragment))
      print ("Model complex of DNA in fragment "+output_fragment)
      modelsize  = {}
      conformers = {}
      selected   = {}
      for file in os.listdir(args["folder"]):
        if file.startswith("dna"): continue
        if not file.endswith(".pdb"): continue
        #WARNING: the PDB files should always have a UNiprotAccession code of 6 characters followed by a '.'
        m=re.search("(\S+).(\S+)_(\S+).(\d+)-(\d+):(\d+):(\d+)_(\S+).pdb",file)
        if m:
            start_model = int(m.group(6))
            end_model   = int(m.group(7))
            prot        = m.group(1)
            pdb         = m.group(2)
            chain       = m.group(3)
            binding     = "%d-%d"%(int(m.group(4)),int(m.group(5)))
            conformers.setdefault((prot,pdb,binding),set()).add(chain)
            modelsize.setdefault((prot,pdb,binding),[]).append(end_model-start_model)
        else:
            start_model = 0
            end_model   = 0
            
        #if file[6] == ".":
        #   code_argument = file.split(":")[0]
        #   codes         = code_argument.split(".")
        #   try:
        #      start_model   = int(file.split(":")[1])
        #   except:
        #      start_model   = 0 
        #   try:
        #      end_model = int(file.split(":")[2].split("_")[0])
        #   except:
        #      end_model = 0
        #   if len(codes) == 3:
        #      prot,pdb_chain,binding = codes
        #      pdb,chain   = pdb_chain.split("_")
        #      conformers.setdefault((prot,pdb,binding),set()).add(chain)
        #      modelsize.setdefault((prot,pdb,binding),[]).append(end_model-start_model)
      
      for   prot,pdb,binding in    conformers.keys():
            print("\t-- TF %s PDB %s REGION %s = Chains %s"%(prot,pdb,binding,conformers[(prot,pdb,binding)]))
            if (prot, binding) in selected.keys():
               if len(conformers[(prot,pdb,binding)]) >= len(conformers[(prot,selected[(prot, binding)],binding)]): 
                  if len(conformers[(prot,pdb,binding)]) > len(conformers[(prot,selected[(prot, binding)],binding)]): 
                         selected[(prot, binding)] = pdb
                  else:
                     if max(modelsize[(prot,pdb,binding)]) > max(modelsize[(prot,selected[(prot, binding)],binding)]): 
                         selected[(prot, binding)] = pdb
            else:
               selected[(prot, binding)] = pdb

      #Create a hypothetical dictionary of dimers if necessary
      if len(dimers.keys())<=0:
         exact_dimers,dimers= create_dimers(selected,conformers,args)

      for prot,binding in selected.keys():
          print("\t-- SELECTED TF %s REGION %s = TEMPLATE %s Chains %s"%(prot,binding,selected[(prot,binding)],conformers[(prot,selected[(prot, binding)],binding)]))

      #Make the list of exact dimers
      index,namemap,filelist,added_proteins = make_dimers_list(index,namemap,filelist,added_proteins,selected,conformers,exact_dimers,args)

      #Make the list of the rest of dimers
      index,namemap,filelist,added_proteins = make_dimers_list(index,namemap,filelist,added_proteins,selected,conformers,dimers,args)

      #Make the list of monomers
      index,namemap,filelist,added_proteins = make_monomers_list(index,namemap,filelist,added_proteins,selected,conformers,dimers,args)

      #Make the list of paired chains
      paired_chains=make_paired_chains(namemap,exact_dimers)
      print("CHAINS OF DIMERS")
      for p,q in paired_chains:
          print("\t-- DIMER "+p+" with "+q)
      # Readme generation.
      with open(os.path.join(args["output_folder"], output_fragment, "readme_"+output_fragment.lstrip("fragment")+".txt"), "w") as readme:
        readme.write("DEFINITIONS.\n======================================================\n\nX and Y are DNA chains. \n")
        for model in namemap.keys():
            readme.write("Chain " + model + " corresponds to file: " + os.path.basename(namemap[model]) + "\n"  )
      ####

      # Start to mount the complex and change the id of DNA to X and Y.
      parser = PDBParser(PERMISSIVE=1)
      complex_pdb = list(parser.get_structure('Complex', startfile)[0].get_chains())
      for chain, letter in zip(complex_pdb, ["X","Y"]):
        chain.id = letter
      start_chain_ids  = [chain.id for chain in complex_pdb]

      dummy_name=recursive_build(complex_pdb, start_chain_ids, filelist, complexname, args, namemap,output_fragment)

      for pdb_file in os.listdir(os.path.join(args["output_folder"],output_fragment)):
        if pdb_file.endswith(".non_opt.pdb"):continue
        if pdb_file.endswith(".pdb") and pdb_file.startswith("dna"):
          pdb_output = "dna_"+output_fragment.lstrip("fragment")+pdb_file.lstrip("dna")
          print("PDB "+pdb_file)
          print("Write PDB "+pdb_output)
          letters = list(pdb_file[4:-4].split("_"))
          print("Chains "+str(letters))
          print("Remarks "+"remarks__"+"_".join(sorted(letters)) +".txt")
          print("Check "+os.path.join(args["output_folder"],output_fragment, "remarks__"+"_".join(sorted(letters)) +".txt"))
          if os.path.exists(os.path.join(args["output_folder"],output_fragment, "remarks__"+"_".join(sorted(letters)) +".txt")):
            print("Add remarks on ... "+pdb_output)
            os.system("cat %s > %s"%(os.path.join(args["output_folder"],output_fragment,"remarks__"+"_".join(sorted(letters)) +".txt"),os.path.join(args["output_folder"],output_fragment,pdb_output)))
            os.system("cat %s >> %s"%(os.path.join(args["output_folder"],output_fragment,pdb_file),os.path.join(args["output_folder"],output_fragment,pdb_output)))
            os.remove(os.path.join(args["output_folder"],output_fragment,pdb_file))
          elif os.path.exists(os.path.join(args["output_folder"],output_fragment, "remarks_"+"_".join(sorted(letters)) +".txt")):
            print("Add remarks on ... "+pdb_output)
            os.system("cat %s > %s"%(os.path.join(args["output_folder"],output_fragment,"remarks_"+"_".join(sorted(letters)) +".txt"),os.path.join(args["output_folder"],output_fragment,pdb_output)))
            os.system("cat %s >> %s"%(os.path.join(args["output_folder"],output_fragment,pdb_file),os.path.join(args["output_folder"],output_fragment,pdb_output)))
            os.remove(os.path.join(args["output_folder"],output_fragment,pdb_file))

      with open(os.path.join(args["output_folder"], output_fragment, "readme_"+output_fragment.lstrip("fragment")+".txt"), "a") as readme:
        readme.write("\nRESTRICTIONS FOUND\n======================================================\n")
        for chain_id in namemap.keys():
          if incompatible.get(chain_id) is not None:
              readme.write("Incompatile chains  %s = %s\n"%(chain_id,str(incompatible.get(chain_id))))
              print("Incompatile chains on DNA %s %s = %s"%(startfile,chain_id,str(incompatible.get(chain_id))))
          if chain_id in unreferenced:
              readme.write("Failed to superimpose chain %s\n"%(chain_id))
              print("Unreferenced chains on DNA %s = %s "%(startfile,chain_id))
          if chain_id in wrong:
              readme.write("Failed to superimpose chain %s\n"%(chain_id))
              print("Wrong chains on DNA %s = %s "%(startfile,chain_id))
      readme.close()
      print("\n======================================================\n")


    print("Done")

