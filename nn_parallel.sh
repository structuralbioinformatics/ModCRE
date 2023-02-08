#!/bin/sh
#$BATCH -l h_vmem=8G
#$BATCH --mem=8G


module load Python/2.7.9
module load MEME/5.0.2-foss-2016b

cd [path modcre]

\cp nearest_neighbour/sequences/* pbm/sequences
\rm -r nearest_neighbour/TF_MOTIFS

python scripts/nearest_neighbour.py --dummy=dummy_nn -i nearest_neighbour/CisBP_2.00_nn_motifs.dat -m nearest_neighbour/CisBP_2.00_nn_models.dat -f pbm/CisBP_2019/CisBP_2.00.all.tf_families.sql -t pbm/CisBP_2019/CisBP_2.00.all.tfs.sql --pwm=pbm/CisBP_2019/pwms --pwm_models=nearest_neighbour/PWM_MODELS --pbm=pbm --verbose --parallel --skip -o nearest_neighbour/TF_MOTIFS

