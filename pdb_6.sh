#!/bin/sh
#$ -l h_vmem=200G
#$ --mem=32G

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

python scripts/pdb.py --dummy=./dummy_pdb -o pdb -p pdb -t tfinder_all/tfs.txt  -v --start=6 --stop=6 --parallel


