#!/usr/bin/make -f

CC ?= cc
CCFLAGS = -std=c99 -O3 -Wall -Wfatal-errors -Winline
INCLUDE = -I../include
RM ?= rm -f

PROGS = example_point3d
PROGS += example_floatbox
PROGS += example_boxfiler
#PROGS += example_mpihelpers
PROGS += example_intersect

all: $(PROGS)

$(PROGS):%: %.c
	$(CC) $(CCFLAGS) $(INCLUDE) -o $@ $<

clean:
	$(RM) $(PROGS)
