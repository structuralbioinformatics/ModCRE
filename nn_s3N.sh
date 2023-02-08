#!/bin/sh
#$ -l h_vmem=200G
#$ --mem=50G

module load MEME
module load Python/2.7.15-foss-2018b

cd  [path modcre]

\rm -r nearest_neighbour/TF_MOTIFS_S3N
\rm -r nearest_neighbour/TF_MOTIFS_S3N.dat
\cp -r nearest_neighbour/TF_MOTIFS nearest_neighbour/TF_MOTIFS_S3N

python scripts/nearest_neighbour.py -s 3 --normalize --top=10 --dummy=dummy_nn -i nearest_neighbour/CisBP_2.00_nn_motifs.dat -m nearest_neighbour/CisBP_2.00_nn_models.dat -f pbm/CisBP_2019/CisBP_2.00.all.tf_families.sql -t pbm/CisBP_2019/CisBP_2.00.all.tfs.sql --pwm=pbm/CisBP_2019/pwms --pwm_models=nearest_neighbour/PWM_MODELS --pbm=pbm_2.0 --verbose  -o nearest_neighbour/TF_MOTIFS_S3N

mkdir nearest_neighbour/PLOT_NORMAL_RANK
\mv nearest_neighbour/TF_MOTIFS_S3N_*graph* nearest_neighbour/PLOT_NORMAL_RANK
\mv nearest_neighbour/TF_MOTIFS_S3N_*plot* nearest_neighbour/PLOT_NORMAL_RANK


