##################
# Configuration  #
##################

  -> Set configuration file parameters (i.e. ./scripts/config.ini).
     Substitute PATH_OF_MODCRE for the address where modcre has been installed (i.e. /sbi/users/modcre)
     Executable binary programs of  DSSP, X3DNA, TMalign, ghostscript, BLAST+ (ncbi-blast-2.2.31+) and HMMER (3) will be installed or linked in folder "src" of modcre
     Modify the addresses in the scripts if you wish to use other  
         program_path = os.path.join(src_path, config.get("Paths", "program_path"))
     to 
         program_path = config.get("Paths", "program_path")
     
     Check the path loacation of these other programs and modify the config file accordingly:
     modppi (from MODPIN)
     Clustal-Omega
     EMBOSS
     MEME
     MMSeqs2
     CD-HIT
     MODELLER
     WEBLOGO (installed in Python)

  -> Necessary commands to submit to queues of a cluster are specified in files
     command_queue  (i.e. command_queues_cluster.txt) for Python 2.7
     command_queue3 (i.e. command_queues_cluster3.txt) for Python 3


##################
# Download files #
##################
-> make folders for the database
mkdir pdb
mkdir pbm
mkdir uniprot

# PDB:
----
-> Go to "http://www.rcsb.org/pdb/search/advSearch.do?search=new".
->  Choose "All/Experimental Type/Molecule Type".
-> "Ignore" Experimental Method and select Molecule Type "Protein/NA complex".
-> Click on "Submit Query".
-> Click on "Download Results".
-> Uncheck all but "PDB".
-> Select Compresion Type "uncompressed".
-> Click on "Launch Download Application" and download "download_rcsb.jnlp".
-> Execute script and download files under "./pdb/all/".
-> Now, uncheck all but "Biological Assemblies".
-> Select Compresion Type "uncompressed".
-> Click on "Launch Download Application" and download "download_rcsb.jnlp".
-> Execute script and download files under "./pdb/biounits/".

# Cis-BP:
-------
-> Go to "http://cisbp.ccbr.utoronto.ca/index.php".
-> Select "PBM" in By Evidence Type" and press "GO" button.
-> Add all TFs to cart and click on "View cart" in the left menu.
-> Click on "Download TFs in cart" and only leave a [check] on "TF info" and "Protein features", then click on "Download TFs in cart!".
-> Go to "http://cisbp.ccbr.utoronto.ca/entireDownload.php".
-> Click on "E-scores" and download "Escores.txt.zip".
-> Unzip downloaded files and you should obtain the following: "TF_Information.txt", "prot_seq.txt" and "Escores.txt".
-> Put the files of CisBP in pbm folder: pbm/CisBP_2019

# UniProt:
--------
-> Go to "http://www.uniprot.org/".
-> Download UniProtKB Swiss-Prot and TrEMBL proteins in FASTA format.
-> Gunzip downloaded files, and concatenate them as: "uniprot_sprot_trembl.fasta" (file is ~40G).
-> Go to "http://www.uniprot.org/docs/speclist".
-> Right click and "Save Page As..." "speclist.txt".
-> On command line, run: grep ">" uniprot_sprot_trembl.fasta > headers.txt -->
-> Get uniprot_sprot and uniref50/90
-> Move to folder uniprot and execute:

curl -O ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
gunzip uniref50.fasta.gz
curl -O ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
gunzip uniref90.fasta.gz
curl -O ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

-> Reformat UniRef in Sprot format
python ../scripts/uniref2SP.py -u uniref50.fasta -o uniref50_sprot.fasta
python ../scripts/uniref2SP.py -u uniref90.fasta -o uniref90_sprot.fasta

-> Get the header we wish to use
grep ">" uniref90_sprot.fasta>header.txt

-> Make the BLAST database
../src/ncbi-blast-2.2.31+/bin/makeblastdb -in uniref90_sprot.fasta -dbtype prot -title uniref90 -out uniref90_sprot.fasta
-> Get idmapping.dat from Uniprot to get cross-references
curl -O ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
gunzip idmapping.dat.gz

##################
# Get TFs        #
##################
-> If SQL files are not available (i.e download full SQL which is incomplete), try getting at least  Cisbp_2.00.motifs.sql and families: 
-> Cisbp_2.00.motifs.sql cisbp_2.0.tf_families.sql
-> get TF information for all TFs plus TF info and proteins for PBM
->  prot_seq_PBM.txt  TF_Information_PBM.txt
-> Modify  prot_seq_PBM.txt: Species must be TF_Species and Protein_seq must be Protein_Sequence
-> Then extract SQL files and TF files from them

python ../../scripts/get_CisBP_Tables.py --sql Cisbp_2.00.motifs.sql --tf TF_Information_all_motifs.txt --ps prot_seq_PBM.txt --select TF_Information_PBM.txt -o CisBP_2.00.PBM &

-> And TF files (you can later select whether you plan to use all TFs downloaded or only those with PBM data extracted as above)

python  scripts/tfinder2.py -o tfinder2019 -p pbm/CisBP_2019/CisBP_2.00.PBM.proteins.sql --pdb=pdb/all -t pbm/CisBP_2019/CisBP_2.00.PBM.tfs.sql  -u uniprot/uniprot_sprot+trembl.fasta -f pbm/CisBP_2019/cisbp_2.0.tf_families.sql --dummy=dummy_tfinder -v
python  scripts/tfinder2.py -o tfinder_all -p pbm/CisBP_2019/CisBP_2.00.all.proteins.sql --pdb=pdb/all -t pbm/CisBP_2019/CisBP_2.00.all.tfs.sql  -u uniprot/uniprot_sprot+trembl.fasta -f pbm/CisBP_2019/CisBP_2.00.all.tf_families.sql --map  uniprot/idmapping.dat --dummy=dummy_tfinder -v 
 

# [warning] Modify the scripts to load the program versions that you have installed, according to the definitions use in the configuration file
##################
# Parse PDB data #
##################

-> All in a run:
-> with all the pdb
  python scripts/pdb.py --dummy=./dummy -o pdb -p pdb -t tfinder/pdb.txt -v
-> with all TFs
  python scripts/pdb.py --dummy=./dummy_tf -o pdb_tf -p pdb -t tfinder2/tfs.txt -v

-> Step by step DB construction
sh pdb_1.sh
sh pdb_2.sh
sh pdb_3.sh
sh pdb_4.sh
sh pdb_5.sh
sh pdb_6.sh
sh pdb_7.sh
sh pdb_7r.sh
sh pdb_8.sh

##################
# Parse PBM data #
##################

-> All in a run to generate data from 2016 (for all data)
-> PBM of  CisBP 2016
  python scripts/pbm.py -e pbm/CisBP_2016/Escores.txt -f pbm/CisBP_2016/cisbp_1.02.tf_families.sql -m pbm/CisBP_2016/cisbp_1.02.motifs.sql \
  -s pbm/CisBP_2016/cisbp_1.02.motif_sources.sql -t pbm/CisBP_2016/cisbp_1.02.tfs.sql -u uniprot/uniprot_sprot+trembl.fasta --dummy=./dummy_tf \
  --pdb=pdb_tf -o pbm_tf --pwm=pbm/CisBP_2016/PWM_2016/pwms/ -v

  sh pbm_1.1.sh

# Step by step DB construction for 2019 data
# Latest PBM CisBP 2019

sh pbm_1.sh
sh pbm_2.sh
sh pbm_3.sh
sh pbm_4.sh
sh pbm_5.sh
sh pbm_6.sh
sh pbm_7.sh
sh pbm_8.sh
sh pbm_9.sh
sh pbm_10.sh
sh pbm_11.sh

#####################################################################################
# Grid Search to get best parameters per family
#####################################################################################
->  Get parameters with script run_grid_search
->  Parameters are inserted in config.ini

sh run_grid_search.sh

##################
# Examples of USE#
##################

#####################################################################################
# Protein to DNA: applications of the model.py, pwm.py and score.py programs
#####################################################################################

1) Model a monomer:
-------------------
-> [example] submit the script "sh model_atoh1.sh" or:
python scripts/model_protein.py -i example/protein2dna/atoh1.fa -l ATOH1 -o example/protein2dna/ATOH1  --pdb=pdb/ -v --dummy=./dummy

2) Model mix of monomer and dimers:
-----------------------------------
python scripts/model_protein.py -i example/protein2dna/esr1.fa -l ESR1 -o example/protein2dna/ESR1  --pdb=pdb/ -v --dummy=./dummy  --full --all --renumerate --chains_fixed --dummy=./dummy --unrestrictive
python scripts/model_protein.py -i example/protein2dna/err1.fa -l ERR1 -o example/protein2dna/ERR1  --pdb=pdb/ -v --dummy=./dummy  --full --all --renumerate --chains_fixed --dummy=./dummy 

3) Model a heterodimer:
-----------------------
-> [warning] if the tfs do not dimerize, it won't model anything 
python scripts/model_protein.py -i example/protein2dna/hetero.fa -l HETERO -o example/protein2dna/HETERO/  --pdb=pdb/ -v --dummy=./dummy


4) Model several protein-dna structures from threading files:
-------------------------------------------------------------
-> inputs a file with a list of threading files for several proteins and make the models for all parallelizing in a cluster
python scripts/model_multiple_proteins.py  -i example/protein2dna/threads.list --threading --pdb=pdb/ -v --parallel --dna -o example/protein2dna/THREAD_LIST --dummy=./dummy


5) Model a DNA sequence:
------------------------
-> [warning] first, calculate DNA interface with interface.py (if not provided, model.py will calculate it)
python scripts/interface.py -i example/protein2dna/ERR1/ERR1:73:153_1by4_B_1.pdb -o example/protein2dna/ERR1/ERR1:73:153_1by4_B_1.interface --dummy=./dummy
python scripts/model_dna.py -p example/protein2dna/ERR1/ERR1:73:153_1by4_B_1.pdb -o example/protein2dna/ERR1/ERR1:73:153_1by4_B_1.AAAAAAAAAAAAAAA \
   --pdb=pdb/ -s AAAAAAAAAAAAAAA -i example/protein2dna/ERR1/ERR1:73:153_1by4_B_1.interface 


6) Predict PWM of a single model:
--------------------------------
-> [warning] not all combinations of statistical potentials exist (should you want them modify pdb.py and pbm.py)
python scripts/pwm_pbm.py -i examples/protein2dna/ATOH1/ATOH1:160:216_1mdy_A_1.pdb -o  examples/protein2dna/ATOH1/ATOH1:160:216_1mdy_A_1 --pbm=pbm/ --pdb=pdb/ -v --auto


7) Predict PWM of a folder with models:
--------------------------------------
-> [warning] not all combinations of statistical potentials exist (should you want them modify pdb.py and pbm.py)
python scripts/pwm_pbm.py -i examples/protein2dna/ATOH1 --pbm=pbm/ --pdb=pdb/ -v --auto --parallel


8) Score an interaction:
------------------------
-> [warning] not all combinations of statistical potentials exist (should you want them modify pdb.py and pbm.py)
python scripts/scorer.py -i example/protein2dna/ATOH1/ATOH1:160:216_1mdy_A_1.pdb  --dummy=dummy -o example/protein2dna/ATOH1/ATOH1:160:216_1mdy_A_1.score  --pbm=pbm/ --pdb=pdb/  --norm -v -a  -s allexample/protein2dna/ATOH1/ATOH1:160:216_1mdy_A_1.score
 
9) Profile an interaction:
------------------------
-> Make a profile of several TF listed in the input
echo "example/protein2dna/ERR1" >  example/protein2dna/TF_profiler_ESR1_ERR1.dat
echo "example/protein2dna/ESR1" >>  example/protein2dna/TF_profiler_ESR1_ERR1.dat
python scripts/xprofiler.py --chains_fixed --dummy=./dummy -d example/protein2dna/ERR1P/dna_large.fa -i example/protein2dna/TF_profiler_ESR1_ERR1.dat  \
  --pbm=./pbm_2.0 --pdb=./pdb  -v  --plot --html --html_types normal,energy,energy_per_nucleotide,energy_best,fimo_log_score,fimo_binding,fimo_score \
  --html_energies s3dc_dd,s3dc,pair   -o  ESR1_ERR1.profile -l thread --auto  --parallel --radius 22.0 --reuse

######################################################################################################
# DNA to protein: characterization of the human enhanceosome at the ifn-beta gene promoter
######################################################################################################

1) Scan a DNA sequence:
-----------------------
-> [warning] execution time grows with the length of the input nucleotide sequence  
python scripts/scanner.py --dummy=./dummy_scan -i example/enhancer/dna.fa -l IFNB_HUMAN  --pbm=./pbm --pdb=./pdb -v \
  -o example/enhancer/IFNB_HUMAN --parallel --reuse --max 100 --ft 0.005 -s 9606 --rank 

2) Build models from threading files:
-------------------------------------
python scripts/model_protein.py -i example/enhancer/IFNB_HUMAN/aux_files/D3DR86.2i9t_B.203-216.txt -l NFKB2.2i9t_B.203-216 \
      -o example/enhancer/IFNB_HUMAN/models/ --pbm=pbm/ --pdb=pdb/ -t -v 
python scripts/model_protein.py -i example/enhancer/IFNB_HUMAN/aux_files/Q04864.2ram_A.208-220.txt -l REL.2ram_A.208-220 \
      -o example/enhancer/IFNB_HUMAN/models/ --pbm=pbm/ --pdb=pdb/ -t -v


######################################################################################################
# Compare the PWM predictions with the nearest-neighbour approach
######################################################################################################

1) Option 1, single run:
----------------------
sh run_nearest_neighbour.sh

2) Option 2, step by step:
------------------------
-> Get data of TF sequences:  input sequences and separated Fasta files
nn_get_data.sh
-> Obtain models
nn_model_tfs.sh
-> Copy sequences in pbm/sequences folder
\cp nearest_neighbour/sequences/* pbm/sequences
-> Obtain theoretical PWMs of each TF
nn_pwm_tfs.sh
-> Get the input file for models of PWMs
nn_get_data_mdl.sh
-> Calculate the nearest neighbours in a cluster
nn_parallel.sh
-> Runs each confdition of NN
-> Compare by ranks
->results for the comparison are in nearest_neighbour/PLOT_NORMAL_RANK
nn_s3N.sh
-> Compare by rank-enrichment
-> Results for the comparison are in nearest_neighbour/PLOT_NORMAL_RANKENRICHED
nn_s5N.sh



