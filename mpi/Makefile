CC ?= gcc
MPI_COMPILE_FLAGS = $(shell mpicc --showme:compile) -std=c99 -Winline
MPI_LINK_FLAGS = $(shell mpicc --showme:link)
INCLUDE = -I../include
LIBS = -fopenmp
DOCS=../docs

PROGS = sweep-mpi-omp

all: $(PROGS)

$(PROGS) : % : %.c
	$(CC) $(MPI_COMPILE_FLAGS) $(LIBS) $(INCLUDE) $< $(MPI_LINK_FLAGS) -o $@

job: sweep-mpi-omp
	qsub -q generic mpirun.sh

clean:
	$(RM) mpirun.sh.* output $(PROGS)
