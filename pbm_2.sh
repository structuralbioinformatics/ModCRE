#!/bin/sh
#$ -l h_vmem=200G
#$ --mem=4G

module load EMBOSS/6.6.0-foss-2016b 
module load ClustalW2 
module load Clustal-Omega
module load dssp
module load Modeller
module load MEME
module load ncbi-blast/2.2.31
module load HMMER/3.1b1
module load TMalign/20160521
module load x3dna/2.2
module load hbplus/3.2
module load Python/2.7.12-foss-2016b 

cd [path modcre]

python scripts/pbm.py --dummy=./dummy_pbm -o pbm --pdb=pdb  --parallel -v --start=2 --stop=2   -t pbm/CisBP_2019/CisBP_2.00.all.tfs.sql -e pbm/CisBP_2019/Escores.txt -f pbm/CisBP_2019/CisBP_2.00.all.tf_families.sql -m pbm/CisBP_2019/CisBP_2.00.all.motifs.sql  -s pbm/CisBP_2019/CisBP_2.00.all.motif_sources.sql -t pbm/CisBP_2019/CisBP_2.00.all.tfs.sql -u uniprot/uniprot_sprot+trembl.fasta  --pwm=pbm/CisBP_2019/pwms 

