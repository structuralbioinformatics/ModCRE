import os, sys, re
import numpy
import optparse
from Bio import motifs as mm
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import ConfigParser

# Get scripts path (i.e. ".") #
scripts_path = os.path.abspath(os.path.dirname(__file__))

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(scripts_path, "config.ini")
config.read(config_file)

# Imports my functions #
import functions

# Imports jbonet's module #
from SBI.data import aminoacids1to3, aminoacids_polarity_boolean, nitrogenous_bases
from SBI.structure import PDB

#cisbp_sql="Cisbp_2.00.sql"
#tf_info="TF_Information_all_motifs.txt"
#protein="prot_seq.txt"

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python --sql cisbp_sql_file  --tf tf_info_file --ps proteins_file [ -o rootname -v --select TF_selection_file]")

    parser.add_option("--sql", action="store", default=None, type="string", dest="sql_file", help="SQL file (from CIS-BP)", metavar="{file}")
    parser.add_option("--tf", action="store", default=None, type="string", dest="tfs_file", help="TF info of all (from CIS-BP)", metavar="{file}")
    parser.add_option("--ps", action="store", default=None, type="string", dest="prot_file", help="Protein sequences (from CIS-BP)", metavar="{file}")
    parser.add_option("--select", action="store", default=None, type="string", dest="tf_select", help="TF info of selected TFs (from CIS-BP)", metavar="{file}")
    parser.add_option("-o", action="store", type="string",default="CisBP",dest="root", help="CisBP file in SQL format with table of motifs (default CisBP)", metavar="{rootname}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")

    (options, args) = parser.parse_args()

    if options.prot_file is None or options.tfs_file is None  or options.sql_file is None:
         parser.error("missing arguments: type option \"-h\" for help")

    return options

def read_SQL_Tables(sql_file):
    sql=open(sql_file,"r")
    table_name=None
    db_table={}
    start_table=False
    close_table=False
    for line in sql:
        data=line.strip().split()
        block=line.strip().split(",")
        if len(data)<1: continue
        if len(block)<1: continue
        if data[0]=="INSERT" and data[1]=="INTO":
          table_name=None
          start_table=True
          for w in xrange(len(data)):
              if data[w].startswith("(") and table_name is None: 
                 table_name=data[w-1].rstrip("`").lstrip("`")
          keys=[]
          value_data=block[0].split()
          value=value_data[-1].lstrip("(`").rstrip("`")
          keys.append(value)
          for word in block[1:-1]:
            value_data=word.split()
            value=value_data[-1].rstrip("`").lstrip("`")
            keys.append(value)
          value_data=block[-1].split()
          value=value_data[0].lstrip("`").rstrip("`)")
          keys.append(value)
          #print "TABLE",table_name, "KEYS", keys
          continue
        if start_table:
          table=[]
          values=[]
          for word in block:
            value_data=word.split()
            if len(value_data)<1: continue
            value=value_data[-1]
            if value.startswith("("): value=value.lstrip("(")
            if value.endswith(";"):
               close_table=True 
               value=value.rstrip(";")
            if value.endswith(")"): value=value.rstrip(")")
            values.append(value.rstrip("'").lstrip("'"))
          #print "TABLE",table_name, "VALUES",values, "LEN",len( values ),"KEYS",keys,"LEN",len(keys)
          if len( values ) == len(keys) :
            for i in xrange( len(keys) ):
               table.append((keys[i],values[i]))
            db_table.setdefault(table_name,set()).add(tuple(table))
        if close_table:
           start_table=False
           close_table=False
           #db_table.setdefault(table_name,table)
    
    return   db_table

def write_tf_table(rootname,table_name,table_keys,tf_info,prot_info,db_table):

    outfile=open(os.path.abspath(rootname+"."+table_name+".sql"),"w")
    main_key=table_keys[0]
    dbid={}
    for tableSQL,set_of_tuples in db_table.iteritems():
      if tableSQL!="motifs":continue
      for pairedlist in set_of_tuples:
        dict_tuple=dict(pairedlist)
        tf=None
        dom=None
        if dict_tuple.has_key(main_key): tf=dict_tuple.get(main_key)
        if dict_tuple.has_key("DBID"): dom=dict_tuple.get("DBID")
        if tf is not None and dom is not None: dbid.setdefault(tf,dom)
    
    merged_table={}
    for tf_id,set_of_tuples in tf_info.iteritems():
      if not dbid.has_key(tf_id): continue
      for pairedlist in set_of_tuples:  
        merged_tuples=set(list(pairedlist))
        merged_tuples.add(("DBID",dbid[tf_id]))
        main_value=("TF_ID",tf_id)
        dictionary=dict(list(merged_tuples))
        data={}
        for key,value in dictionary.iteritems():
            if key in table_keys:
               data.setdefault(key,value)
        merged_table.setdefault(main_value,[]).append(data)
 
    clean_table={}
    for key,list_of_data in merged_table.iteritems():
        for data in list_of_data:
            clean_table.setdefault(key,set()).add(tuple([(k,v) for k,v in data.iteritems()]))
    #print clean_table
    outfile.write("CREATE TABLE IF NOT EXISTS `%s` (\n"%table_name)
    first=True
    for key in table_keys:
      if first:
        outfile.write("\t`%s` varchar(20) NOT NULL,\n"%key)
        first=False
      else: 
        outfile.write("\t`%s` varchar(20) DEFAULT NULL,\n"%key)
    
    outfile.write("\tPRIMARY KEY (`%s`)\n"%table_keys[0])
    outfile.write(") ENGINE=MyISAM DEFAULT CHARSET=latin1;\n")
    outfile.write("INSERT INTO `%s` ("%table_name)
    for key in table_keys[0:-1]:
        outfile.write("`%s`, "% key)
    outfile.write("`%s`) VALUES\n"% table_keys[-1])
    for key,set_of_data in clean_table.iteritems():
      for data in set_of_data:
        skip=False
        dictionary=dict(data)
        for key_table in table_keys:
            if not dictionary.has_key(key_table): skip=True
        if skip: continue
        outfile.write("(")
        for key_table in table_keys[0:-1]:
            outfile.write("'%s', "%dictionary.get(key_table))
        outfile.write("'%s'),\n"%dictionary.get(table_keys[-1]))
    outfile.write(";\n")
    outfile.close()    

       
        

def write_table(rootname,table_name,table_keys,tf_info,prot_info,db_table):

    outfile=open(os.path.abspath(rootname+"."+table_name+".sql"),"w")
    
    main_key=table_keys[0]
    partial={}
    for tableSQL,set_of_tuples in db_table.iteritems():
        second_key = None
        list_of_tuples=list(set_of_tuples)
        if main_key not in dict(list_of_tuples[0]).iterkeys():
          for check_key in table_keys:
            if check_key in dict(list_of_tuples[0]).iterkeys() and second_key is None:
             second_key=check_key
        else:
          second_key=main_key
        if second_key is None: continue
        #print tableSQL
        for pairedlist in set_of_tuples:
            #print pairedlist
            data=[]
            key_value=None
            for key,value in pairedlist:
                if key in table_keys:
                   data.append((key,value))
                   if key==second_key: key_value=value
            if key_value is not None: 
               #print "Found in", tableSQL
               #print "Key",second_key,key_value
               #print "Add",tuple(data)
               partial.setdefault((second_key,key_value),set()).add(tuple(data))
    #print tf_info
    for tf_id,set_of_tuples in tf_info.iteritems():
        second_key = None
        list_of_tuples=list(set_of_tuples)
        if main_key  not in dict(list_of_tuples[0]).iterkeys():
          for check_key in table_keys:
            if check_key in dict(list_of_tuples[0]).iterkeys() and second_key is None:
             second_key=check_key
        else:
          second_key=main_key
        if second_key is None: continue
        for pairedlist in set_of_tuples:
            data=[]
            for key,value in pairedlist:
                if key in table_keys:
                   data.append((key,value))
                   if key==second_key: key_value=value
            if key_value is not None: 
               #print "Found in", tf_id
               #print "Key",second_key,key_value
               #print "Add",tuple(data)
               partial.setdefault((second_key,key_value),set()).add(tuple(data))
    #print prot_info
    for p_id,set_of_tuples in prot_info.iteritems():
        second_key = None
        list_of_tuples=list(set_of_tuples)
        if main_key  not in dict(list_of_tuples[0]).iterkeys():
          for check_key in table_keys:
            if check_key in dict(list_of_tuples[0]).iterkeys() and second_key is None:
             second_key=check_key
        else:
          second_key=main_key
        if second_key is None: continue
        for pairedlist in set_of_tuples:
            data=[]
            for key,value in pairedlist:
                #print "CHECK PROTEIN",p_id,key,value
                if key in table_keys:
                   data.append((key,value))
                   if key==second_key: key_value=value
            #print "      DATA",p_id,second_key,key_value,data
            if key_value is not None: 
               #print "Found in", p_id
               #print "Key",second_key,key_value
               #print "Add",tuple(data)
               partial.setdefault((second_key,key_value),set()).add(tuple(data))
    merged_table={}
    for double in xrange(2):
      for key_pair,set_of_tuples in partial.iteritems():
        #print "CHECK ROUND",double,key_pair
        if double==1: print "CHECK MAIN", table_name, main_key,key_pair
        #print set_of_tuples
        main_value_set=set()
        for  k,v in merged_table.iteritems():
             if key_pair==k: main_value_set.add(k)
             for d in v:
                 add_values=set([(kk,vv) for kk,vv in d.iteritems()])
                 if key_pair in add_values: main_value_set.add(k)
        if len(main_value_set)>0 and double==1:print "Main value",len(main_value_set),"CHECK MAIN", table_name, main_key,key_pair
        link=False
        for pairedlist in set_of_tuples:
            dictionary=dict(pairedlist)
            keys_in=set([x for x in dictionary.iterkeys()])
            main_value_pass=set()
            if main_key in dictionary.iterkeys():
               main_value=(main_key,dictionary.get(main_key))
               main_value_pass.add(main_value)
               link=True
            if not link:
             for main_value_link in main_value_set:
               #print "CHECK KEYS IN",keys_in,"MAIN VALUE",table_name,main_value_link
               if merged_table.has_key(main_value_link):
                  listdata=merged_table.get(main_value_link)
                  keys_done=set()
                  for i in xrange(len(listdata)):
                      keys_done.update(set([x for x in listdata[i].iterkeys()]))
                  if  len(keys_done.intersection(keys_in))>0: 
                      main_value_pass.add(main_value_link)
                      link=True  
            if not link: continue
            #print "PASSES MAIN",main_key,key_pair,main_value, "PASSES",len(main_value_pass)
            for main_value in main_value_pass:
              #print "PASSES MAIN",main_key,key_pair,main_value
              data={}
              key_data=set()
              for key,value in dictionary.iteritems():
                   if key in table_keys:
                      data.setdefault(key,value)
                      key_data.add(key)
              if merged_table.has_key(main_value):
                listdata=merged_table.get(main_value)
                if len(listdata)>0:
                 for i in xrange(len(listdata)):
                   dictionary_done=listdata[i]
                   #print " Old",table_name,i,main_value,merged_table[main_value][i]
                   keys_done=set([x for x in dictionary_done.iterkeys()])
                   merge=True
                   for key in key_data.intersection(keys_done):
                       if data.get(key)=="." and dictionary_done.get(key) != ".": data[key]=dictionary_done.get(key)
                       if dictionary_done.get(key)=="." and data.get(key) != ".": dictionary_done[key]= data.get(key)
                       if dictionary_done.get(key) != data.get(key):
                          merge=False
                          if  dictionary_done.get(key) in data.get(key) or  data.get(key) in dictionary_done.get(key): merge=True
                   if merge:
                      #print "Merge",table_name,main_value,data,dictionary_done
                      for key,value in dictionary_done.iteritems():
                        merged_table[main_value][i].setdefault(key,value)
                      for key,value in data.iteritems():
                        merged_table[main_value][i].setdefault(key,value)
                      #print "Modify ",table_name,i,main_value,merged_table[main_value][i]
              else:    
                #print "Add", table_name, main_value,    data
                merged_table.setdefault(main_value,[]).append(data)
    #Clean merged list
    #print    merged_table
    clean_table={}
    for key,list_of_data in merged_table.iteritems():
        for data in list_of_data:
            clean_table.setdefault(key,set()).add(tuple([(k,v) for k,v in data.iteritems()]))
    #print clean_table
    outfile.write("CREATE TABLE IF NOT EXISTS `%s` (\n"%table_name)
    first=True
    for key in table_keys:
      if first:
        outfile.write("\t`%s` varchar(20) NOT NULL,\n"%key)
        first=False
      else: 
        outfile.write("\t`%s` varchar(20) DEFAULT NULL,\n"%key)
    
    outfile.write("\tPRIMARY KEY (`%s`)\n"%table_keys[0])
    outfile.write(") ENGINE=MyISAM DEFAULT CHARSET=latin1;\n")
    outfile.write("INSERT INTO `%s` ("%table_name)
    for key in table_keys[0:-1]:
        outfile.write("`%s`, "% key)
    outfile.write("`%s`) VALUES\n"% table_keys[-1])
    for key,set_of_data in clean_table.iteritems():
      for data in set_of_data:
        skip=False
        dictionary=dict(data)
        for key_table in table_keys:
            if not dictionary.has_key(key_table): skip=True
        if skip: continue
        outfile.write("(")
        for key_table in table_keys[0:-1]:
            outfile.write("'%s', "%dictionary.get(key_table))
        outfile.write("'%s'),\n"%dictionary.get(table_keys[-1]))
    outfile.write(";\n")
    outfile.close()    



#-------------#
# Main        #
#-------------#

def main():
    #Parameters required for ModCRE
    cisbp_tables={}
    cisbp_tf_table={}
    cisbp_tables.setdefault("tf_families",["Family_ID", "Family_Name", "DBDs", "DBD_Count", "Cutoff"])
    cisbp_tables.setdefault("motifs",["Motif_ID", "TF_ID", "MSource_ID", "DBID", "Motif_Type", "Motif_Sequence", "IUPAC", "IUPAC_REV"])
    cisbp_tables.setdefault("motif_sources",["MSource_ID", "MSource_Identifier", "MSource_Type", "MSource_Author", "MSource_Year", "PMID", "MSource_Version"])
    cisbp_tf_table.setdefault("tfs",["TF_ID", "Family_ID", "TSource_ID", "DBID", "TF_Name", "TF_Species", "TF_Status"])
    cisbp_tables.setdefault("proteins",["Protein_ID", "TF_ID", "DBID", "TF_Species", "Protein_Sequence"])

    # Arguments & Options #
    options = parse_options()
    db_table= read_SQL_Tables(os.path.abspath(options.sql_file))
    #Read Selection
    use_selection=False
    if options.tf_select is not None:
       use_selection=True
       select=set()
       titles=[]
       for line in functions.parse_file(os.path.abspath(options.tf_select)):
        tf_dict=[]
        if len(titles)<=0:
           titles=line.strip().split()
        else:
           values=line.strip().split()
           tf_id=values[0]
           select.add(tf_id)
    #Read TF_info
    titles=[]
    tf_info={}
    for line in functions.parse_file(os.path.abspath(options.tfs_file)):
        tf_dict=[]
        if len(titles)<=0:
           titles=line.strip().split("\t")
        else:
           values=line.strip().split("\t")
           tf_id=values[0]
           skip=False
           if  use_selection:
               if tf_id in  select: 
                  skip=False
               else: 
                  skip=True
           if skip:continue
           if len(titles) == len(values):
              for i in xrange(len(titles)):
                  tf_dict.append((titles[i],values[i]))
              tf_info.setdefault(tf_id,set()).add(tuple(tf_dict))
    #Read protein_sequences
    titles=[]
    #print "PROTEIN"
    prot_info={}
    #print "Open ",options.prot_file
    for line in functions.parse_file(os.path.abspath(options.prot_file)):
        prot_dict=[]
        if len(titles)<=0:
           titles=line.strip().split("\t")
           #print titles
        else:
           values=line.strip().split("\t")
           p_id=values[0]
           if len(titles) == len(values):
              for i in xrange(len(titles)):
                  prot_dict.append((titles[i],values[i]))
              prot_info.setdefault(p_id,set()).add(tuple(prot_dict))
              #print p_id,prot_dict
    #print prot_info
    #Create tables
    for table_name,table_keys in cisbp_tables.iteritems():
        write_table(options.root,table_name,table_keys,tf_info,prot_info,db_table)
    for table_name,table_keys in cisbp_tf_table.iteritems():
        write_tf_table(options.root,table_name,table_keys,tf_info,prot_info,db_table)


if __name__ == "__main__":
   main()


        
              
