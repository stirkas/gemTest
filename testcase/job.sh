#!/bin/bash -l

#SBATCH -q debug
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -J test
#SBATCH -C haswell
#SBATCH -A mp118

#cd $SLURM_SUBMIT_DIR
mkdir -p matrix
mkdir -p out
mkdir -p dump

#export OMP_PROC_BIND=true
#export OMP_PLACES=threads
export OMP_NUM_THREADS=8
#export OMP_STACKSIZE=1G

#Run locally.
mpirun -n 2 ./gem_main > run.out 2> run.err
#Run on Cori.
#srun -n 2 ./gem_main > run.out 2> run.err

wait
