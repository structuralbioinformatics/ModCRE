#!/bin/sh
#$ -l h_vmem=20G
#$ --mem=8000

module load dssp
module load CD-HIT
module load EMBOSS
module load MEME
module load Modeller/9.13
module load BLAST+/2.2.31 
module load HMMER/3.1b2 
module load TMalign
module load x3dna/2.3 
module load HBPLUS/3.2
module load Python/2.7.11

cd /sbi/users/interchange/boliva/modcre

python scripts/model_protein.py -i example/protein2dna/atoh1.fa -l ATOH1 -o example/protein2dna/ATOH1  --pdb=/home/boliva/sit_sbi/ModCRE/pdb -v  --full --all  --renumerate


