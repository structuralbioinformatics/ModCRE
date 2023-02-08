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

#Uncomment next if you plan to use Uniprot+TremBL sequences. This will require a huge computational time
#python scripts/pbm.py --dummy=./dummy_pbm -o pbm --pdb=pdb  --parallel -v --start=9 --stop=9   -t pbm/CisBP_2019/CisBP_2.00.all.tfs.sql -e pbm/CisBP_2019/Escores.txt -f pbm/CisBP_2019/CisBP_2.00.all.tf_families.sql -m pbm/CisBP_2019/CisBP_2.00.all.motifs.sql  -s pbm/CisBP_2019/CisBP_2.00.all.motif_sources.sql -t pbm/CisBP_2019/CisBP_2.00.all.tfs.sql -u uniprot/uniprot_sprot+trembl.fasta  --pwm=pbm/CisBP_2019/pwms 

#Uncomment next if you plan to use all Uniprot sequences. This will require a large computational
#python scripts/pbm.py --dummy=./dummy_pbm -o pbm --pdb=pdb  --parallel -v --start=9 --stop=9   -t pbm/CisBP_2019/CisBP_2.00.all.tfs.sql -e pbm/CisBP_2019/Escores.txt -f pbm/CisBP_2019/CisBP_2.00.all.tf_families.sql -m pbm/CisBP_2019/CisBP_2.00.all.motifs.sql  -s pbm/CisBP_2019/CisBP_2.00.all.motif_sources.sql -t pbm/CisBP_2019/CisBP_2.00.all.tfs.sql -u uniprot/uniprot_sprot.fasta  --pwm=pbm/CisBP_2019/pwms 

#With this option you will use a non-redundant large set of sequences with less than50% sequence i
python scripts/pbm.py --dummy=./dummy_pbm -o pbm --pdb=pdb  --parallel -v --start=9 --stop=9   -t pbm/CisBP_2019/CisBP_2.00.all.tfs.sql -e pbm/CisBP_2019/Escores.txt -f pbm/CisBP_2019/CisBP_2.00.all.tf_families.sql -m pbm/CisBP_2019/CisBP_2.00.all.motifs.sql  -s pbm/CisBP_2019/CisBP_2.00.all.motif_sources.sql -t pbm/CisBP_2019/CisBP_2.00.all.tfs.sql -u uniprot/uniref50_sprot.fasta  --pwm=pbm/CisBP_2019/pwms 

