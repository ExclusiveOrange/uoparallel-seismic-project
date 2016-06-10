#!/bin/sh
################################################################################
# mpirun.sh - script to use 'mpirun' to launch sweep-mpi-omp with parameters
################################################################################
#
# see: run.sh
#
################################################################################

MPINUMPROCS=12
DIR=$HOME/uoparallel-seismic-project/mpi/

module load mpi

cd $DIR

mpirun -np $MPINUMPROCS "$DIR/run.sh" > "$DIR/mpirun.log" 2>&1
