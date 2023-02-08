#!/bin/sh
#$ -l h_vmem=200G
#$ --mem=128G

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

python scripts/pbm.py -e pbm/CisBP_2016/Escores.txt -f pbm/CisBP_2016/cisbp_1.02.tf_families.sql -m pbm/CisBP_2016/cisbp_1.02.motifs.sql -s pbm/CisBP_2016/cisbp_1.02.motif_sources.sql -t pbm/CisBP_2016/cisbp_1.02.tfs.sql -u uniprot/uniprot_sprot+trembl.fasta --dummy=./dummy_pbm1 --pdb=pdb_tf -o pbm_1.02 --pwm=pbm/CisBP_2016/PWM_2016/pwms/ -v


