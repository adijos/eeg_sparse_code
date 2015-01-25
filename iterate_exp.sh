#!/bin/sh
#Matlab version of ACS looping over various compression factors
#reality checks
JobFile="experiments.sh"
experiment_name='dummy'
for param in 128 96 65 25 16 4
do
	echo eeg-$param.o
	Outfile=$HOME/Logs/eeg-$experiment_name-$param.o
	Errorfile=$HOME/Errors/eeg-$experiment_name-$param.e
	export num_iters
	export SCRATCH
	export overcomp
	export cmprdim
	export beta
	#jobname is ridiculously short to be easily seen by squeue
	sbatch -o $Outfile -e $Errorfile -J eeg_sc $JobFile
done
