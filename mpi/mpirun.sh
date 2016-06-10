#!/bin/sh
################################################################################
# mpirun.sh - script to use 'mpirun' to launch sweep-mpi-omp with parameters
################################################################################
#
# see: run.sh
#
################################################################################

MPINUMPROCS=2

mpirun -np $MPINUMPROCS "$(pwd)/run.sh" # > $(pwd)/mpirun.log
