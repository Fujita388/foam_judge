#!/bin/sh

#PBS -q i2cpu
#PBS -l select=2:ncpus=128:mpiprocs=128:ompthreads=1
#PBS -l walltime=00:30:00
#PBS -N foam_judge

module load intel intel-mpi
mpiexec -n 256 ./a.out > foam_judge.log 2>&1
