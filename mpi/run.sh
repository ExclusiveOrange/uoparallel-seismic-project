#!/bin/bash
################################################################################
# run.sh - example script for running sweep-mpi-omp with the 241-241-51 data
################################################################################

VELFILE="../data/velocity-241-241-51.vbox"
STARTFILE="../docs/start-1-241-241-51.txt"
FORWARDSTARFILE="../docs/818-FS.txt"
OUTPUTDIRECTORY="../output"
TRAVELTIMEFILENAMEPREFIX="${OUTPUTDIRECTORY}/tt-241-241-51-out"

mkdir $OUTPUTDIRECTORY 2>/dev/null

time ./sweep-mpi-omp $VELFILE $STARTFILE $FORWARDSTARFILE $TRAVELTIMEFILENAMEPREFIX
