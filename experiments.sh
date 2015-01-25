#!/bin/bash -l
source /usr/Modules/init/sh


#SBATCH --partition=cortex
#SBATCH --time=150:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=adijos@berkeley.edu

cd $HOME/eeg_sparse_code/
module load matlab/R2013a
#export MKL_NUM_THREADS=4
pre_post=0

matlab -nodesktop -r "experiment_wrapper($pre_post)'
