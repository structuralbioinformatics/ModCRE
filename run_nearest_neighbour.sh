#Order of execution
#get data of TF sequences:  input sequences and separated Fasta files
nn_get_data.sh
#obtain models
nn_model_tfs.sh
#copy sequences in pbm/sequences folder
\cp nearest_neighbour/sequences/* pbm/sequences
#obtain theoretical PWMs of each TF
nn_pwm_tfs.sh
#get the input file for models of PWMs
nn_get_data_mdl.sh
#run first step of NN
nn_parallel.sh
#Runs each confdition of NN
nn_s3N.sh
nn_s5N.sh

