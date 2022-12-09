#!/bin/sh

#PBS -q F16cpu
#PBS -l select=16:ncpus=128:mpiprocs=128:ompthreads=1
#PBS -l walltime=04:30:00
#PBS -N foam_judge

module load intel intel-mpi
mpiexec -n 2048 ./a.out > foam_judge.log 2>&1
