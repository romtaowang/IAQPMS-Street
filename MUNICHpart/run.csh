#!/bin/csh -f
## This is a script to run NAQPMS by PBS batch system. #
#########################################################
#PBS -N Street_pre_multi
#PBS -q low
#PBS -l nodes=5:ppn=22
#PBS -e log.$PBS_JOBID.err
#PBS -o log.$PBS_JOBID.out
cd $PBS_O_WORKDIR
source /public/lapc/Wangtao/.cshrc
mpirun -np 110 -machinefile $PBS_NODEFILE python sing_preproc.py sing_preproc.cfg 
