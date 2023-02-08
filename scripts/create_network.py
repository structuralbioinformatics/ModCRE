import argparse
import networkx as nx
import pandas as pd
import sys, os, re


def main():
    """
    Requires:
    TcoF database

    Usage:
    python create_network.py -i <input_file> -d <tcof_data_file> -o <output_file>

    Example of command:
    python /home/quim/PHD/Side-Projects/modcre/scripts/create_network.py -i /home/quim/PHD/Side-Projects/modcre/data/uniprot_IDs.txt -d /home/quim/Databases/TcoF-DB/TcoF_merged.txt -o /home/quim/PHD/Side-Projects/modcre/data/result.json
    python /Users/quim/Dropbox/UPF/PHD/Projects/modcre/scripts/create_network.py -i /Users/quim/Dropbox/UPF/PHD/Projects/modcre/data/uniprot_IDs.txt -d /Users/quim/Documents/Databases/TcoF-DB/TcoF_merged.txt -o /Users/quim/Dropbox/UPF/PHD/Projects/modcre/outputs/result.json
    python /Users/quim/Dropbox/UPF/PHD/Projects/modcre/scripts/create_network.py -i /Users/quim/Dropbox/UPF/PHD/Projects/modcre/data/input_trampa.txt -d /Users/quim/Documents/Databases/TcoF-DB/TcoF_merged.txt -o /Users/quim/Dropbox/UPF/PHD/Projects/modcre/outputs/result_trampa.json
    python /Users/quim/Dropbox/UPF/PHD/Projects/modcre/scripts/create_network.py -i /Users/quim/Dropbox/UPF/PHD/Projects/modcre/data/input_no_matches.txt -d /Users/quim/Documents/Databases/TcoF-DB/TcoF_merged.txt -o /Users/quim/Dropbox/UPF/PHD/Projects/modcre/outputs/result_no_matches.json
    """

    options = parse_user_arguments()
    create_network_by_biana_search(options)

def parse_user_arguments(*args, **kwds):

    parser = argparse.ArgumentParser(
        description = "Generate a protein-protein interaction network",
        epilog      = "@oliva's lab 2017")
    parser.add_argument('-i','--input_file',dest='input_file',action = 'store',
                        help = 'Input file')
    parser.add_argument('-d','--data_file',dest='data_file',action = 'store',default='uniprotaccession',
                        help = 'Data file containing the interacting proteins')
    parser.add_argument('-o','--output_file',dest='output_file',action = 'store',
                        help = 'Output file')
    parser.add_argument('-c','--include_cofactors',dest='cofactors',action = 'store_true',default=False,
                        help = 'Flag to include the interactions with the co-factors (i.e. non DNA-binding proteins)')
    parser.add_argument('-v','--verbose',dest='verbose',action = 'store_true',
                        help = 'Flag to use verbose mode')
    options=parser.parse_args()

    return options


def create_network_by_biana_search(options):
    """
    Generates a protein-protein interaction network extracting information from BIANA.
    """
    # READ INPUT FILE
    binding_site_to_tfs, binding_site_to_info, all_tfs = read_input_file(options.input_file)
    #print(all_tfs)


    # READ DATA FILE AND SUBSET BY INPUT TF
    tcof_df = pd.read_csv(options.data_file, sep='\t', header=0)
    subset_tcof_df = tcof_df[ (tcof_df['uniprot1'].isin(all_tfs)) | (tcof_df['uniprot2'].isin(all_tfs)) ]


    # READ SUBSET DATAFRAME
    network = nx.Graph()
    interactors = set()
    interactors_to_tfs = {}
    for index, row in subset_tcof_df.iterrows():
        id1 = str(row['uniprot1'])
        id2 = str(row['uniprot2'])
        if id1 != '-' and id2 != '-':
            for sub_id1 in id1.split('; '):
                for sub_id2 in id2.split('; '):
                    network.add_edge(sub_id1, sub_id2)
                    # Obtain the interactors and the pairs interactor-tf
                    if sub_id1 not in all_tfs: interactors.add(sub_id1)
                    if sub_id2 not in all_tfs: interactors.add(sub_id2)
                    if sub_id1 in interactors and sub_id2 in all_tfs: interactors_to_tfs.setdefault(sub_id1, set()).add(sub_id2)
                    if sub_id2 in interactors and sub_id1 in all_tfs: interactors_to_tfs.setdefault(sub_id2, set()).add(sub_id1)

    # Add the connections between binding sites (ordering them from 5' to 3')
    bs1 = None
    for bs2 in sorted(binding_site_to_info, key=lambda x: int(x.split()[-1]), reverse = False):
        if bs1:
            network.add_edge(bs1, bs2)
        else:
            #network.add_edge("5'", bs2)
            pass
        bs1 = (bs2 + '.')[:-1] # copy the variable
    #network.add_edge(bs1, "3'")

    # Add the connections between binding sites and transcription factors
    for bs in binding_site_to_tfs:
        for tf in binding_site_to_tfs[bs]:
            network.add_edge(bs, tf)


    # Filter interactors that are interacting with at least 2 transcription factors
    filtered_interactors = set([interactor for interactor in interactors_to_tfs if len(interactors_to_tfs[interactor])>1])

    # Filter interactors that only interact with one transcription factor but interact with another interactor connected with a different transcription factor
    one_link_interactors = set()
    for u,v in network.edges():
        if u in interactors and v in interactors:
            u_tfs = interactors_to_tfs[u]
            v_tfs = interactors_to_tfs[v]
            if u_tfs != v_tfs:
                one_link_interactors.add(u)
                one_link_interactors.add(v)


    new_network = nx.Graph()
    for u,v in network.edges():
        # Skip the interactors that are not filtered
        if u in interactors and u not in filtered_interactors and u not in one_link_interactors:
            continue
        if v in interactors and v not in filtered_interactors and v not in one_link_interactors:
            continue
        if options.cofactors:
           new_network.add_edge(u,v)
        else:
           if u in all_tfs and v in all_tfs: 
               new_network.add_edge(u,v)
               continue
           if u in binding_site_to_tfs and v in binding_site_to_tfs: 
               new_network.add_edge(u,v)
               continue
           if u in all_tfs and v in binding_site_to_tfs:
               new_network.add_edge(u,v)
               continue
           if u in binding_site_to_tfs and v in all_tfs: 
               new_network.add_edge(u,v)
               continue


    # CREATE JSON FILE

    output = []

    # First the nodes
    for bs in binding_site_to_tfs:
        [bs_id, family, start, end] = binding_site_to_info[bs][0]
        output.append('{ "data": { "id": "%s", "label":"%s", "type": "binding site" } }' % (bs, bs_id))
    for tf in all_tfs:
        output.append('{ "data": { "id": "%s", "label":"%s", "type": "transcription factor" } }' % (tf, tf))
    if options.cofactors:
        for interactor in filtered_interactors:
            output.append('{ "data": { "id": "%s", "label":"%s", "type": "interactor" } }' % (interactor, interactor))
    #output.append('{ "data": { "id": "5\'", "label":"5\'", "type": "indication" } }')
    #output.append('{ "data": { "id": "3\'", "label":"3\'", "type": "indication" } }')

    # Then the interactions
    for u,v in new_network.edges():
        if v in binding_site_to_info and u not in binding_site_to_info:
            # Put always the binding site as source if the other molecule is not a binding site
            output.append('{ "data": { "id": "%s-%s", "source": "%s", "target": "%s" }  }' % (v, u, v, u))
        # Put always the tfs before the interactors
        elif u in interactors and v in all_tfs:
            output.append('{ "data": { "id": "%s-%s", "source": "%s", "target": "%s" }  }' % (v, u, v, u))
        else:
            output.append('{ "data": { "id": "%s-%s", "source": "%s", "target": "%s" }  }' % (u, v, u, v))

    # Write them in the output file
    with open(options.output_file, 'w') as output_f:
        output_f.write('[{}]'.format(','.join(output)))

    return


def read_input_file(input_file):
    """
    Reads the input file.
    """
    binding_site_to_tfs = {}
    binding_site_to_info = {}
    all_tfs = set()
    with open(input_file, 'r') as input_f:
        binding_site = None
        for line in input_f:
            fields = line.strip().split('\t')
            if fields[0] == '':
                continue
            elif fields[0].lower().startswith('binding site'):
                #binding_site = fields[0].rstrip(':')
                [binding_site, family, start, end] = fields
                _, num = binding_site.lower().split('binding site ') # Get the number of the binding site
                bs_id = 'BS-'+num # Create a binding site id i.e. BS-1
                try:
                    binding_site_to_info.setdefault(binding_site,[]).append([bs_id, family, int(start), int(end)])
                except:
                    continue
            else:
                if binding_site:
                    binding_site_to_tfs.setdefault(binding_site,set()).update(set(fields))
                    for tf in fields:
                        all_tfs.add(tf.upper())
    return binding_site_to_tfs, binding_site_to_info, all_tfs


def fileExist(file):
    """
    Checks if a file exists AND is a file.
    """
    return os.path.exists(file) and os.path.isfile(file)


def create_directory(directory):
    """
    Checks if a directory exists and if not, creates it.
    """
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
    return


if  __name__ == "__main__":
    main()



