~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 README.txt ~ sweep-mpi-omp
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


 About
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Based on seismic raytracing code from CIS 410 Parallel Computing in Spring 2016,
these codes perform repeated iterations of a sweep routine that eventually
finds the shortest travel time from a single given point to each point in
a regular grid, given potentially unique velocities at each point.

These particular codes use MPI and OpenMP to speed-up the operation through
parallelization over multiple hardware threads and multiple CPU cores.


 Building
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You'll need a C99 compiler that supports OpenMP, and an MPI front-end such as
mpicc.

If you're on a cluster, you may need to do:

    > module load mpi

then run:

    > make

It should spit out:
    
    sweep-mpi-omp


 Running
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sweep-mpi-omp <in:velocity.vbox> <in:startpoints.txt> <in:forwardstar.txt>
    <out:traveltimefilenameprefix>


in: velocity.vbox:
This is a custom-format binary file for storing a regular grid of velocities.
There is a utility for converting to this format from a plaintext velocity file:

    uoparallel-seismic-project/tools/velfileconvert.c

..which can be built with the Makefile in that same directory.


in: startpoints.txt:
Currently only the first start point in the given file is used: this is due to
concerns about the overhead of keeping track of multiple travel time arrays
when multiprocessing.
The work-around is to create a separate startpoint.txt file for each starting
point, and run the program separately for each. The total time should actually
be less this way, since less local memory will be used for each run, leading
to better use of low-level caches.


in: forwardstar.txt:
Should work with any forward-star text file as in the ../docs/ directory; however
the ghost size is hardcoded in sweep-mpi-omp.c as 7, so you'll need to recompile
if you want to use a larger forward star.


out: traveltimefilenameprefix:
After the solution has converged, each MPI node will output a separate text file
containing the travel time solutions for its respective region. The filename it
uses will begin with the given traveltimefilenameprefix, and will be appended
with a description of the region coordinates.
The given prefix may include a path, at your option.


 Scripts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You probably won't want to run sweep-mpi-omp directly, since it won't benefit
from MPI that way. Instead, use:

    ./mpirun.sh

or if you've got access to qsub:

    make job


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
