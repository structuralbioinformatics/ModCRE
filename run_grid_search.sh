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

cd [modcre path]

python scripts/grid_search.py --dummy=dummy_gs/ -o grid_search_2.0/ --pbm=pbm_2.0/ --pdb=pdb/ --pwm=pbm/CisBP_2019/pwms/ --start=1 --stop=1 -v -p >& run_grid_1.2.log

python scripts/grid_search.py --dummy=dummy_gs/ -o grid_search_2.0/ --pbm=pbm_2.0/ --pdb=pdb/ --pwm=pbm/CisBP_2019/pwms/ --start=2 --stop=2 -v -p >& run_grid_2.2.log

python scripts/grid_search.py --dummy=dummy_gs/ -o grid_search_2.0/ --pbm=pbm_2.0/ --pdb=pdb/ --pwm=pbm/CisBP_2019/pwms/ --start=3 --stop=3 -v -p >& run_grid_3.2.log

python scripts/grid_search.py --dummy=dummy_gs/ -o grid_search_2.0/ --pbm=pbm_2.0/ --pdb=pdb/ --pwm=pbm/CisBP_2019/pwms/ --start=4 --stop=4 -v -p >& run_grid_4.2.log

python scripts/grid_search.py --dummy=dummy_gs/ -o grid_search_2.0/ --pbm=pbm_2.0/ --pdb=pdb/ --pwm=pbm/CisBP_2019/pwms/ --start=5 --stop=5 -v -p >& run_grid_5.2.log

python scripts/grid_search.py --dummy=dummy_gs/ -o grid_search_2.0/ --pbm=pbm_2.0/ --pdb=pdb/ --pwm=pbm/CisBP_2019/pwms/ --start=6 --stop=6 -v -p >& run_grid_6.2.log

python scripts/grid_search.py --dummy=dummy_gs/ -o grid_search_2.0/ --pbm=pbm_2.0/ --pdb=pdb/ --pwm=pbm/CisBP_2019/pwms/ --start=7 --stop=7 -v -p >& run_grid_7.2.log

python scripts/grid_search.py --dummy=dummy_gs/ -o grid_search_2.0/ --pbm=pbm_2.0/ --pdb=pdb/ --pwm=pbm/CisBP_2019/pwms/ --start=8 --stop=8 -v -p >& run_grid_8.2.log

#Results are in grid_search_2.0/results/grid_search_parameters.txt
