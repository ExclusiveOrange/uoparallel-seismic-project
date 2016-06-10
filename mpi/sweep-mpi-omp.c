////////////////////////////////////////////////////////////////////////////////
// sweep-mpi-omp.c
////////////////////////////////////////////////////////////////////////////////
//
// TODO:
//   * support multiple start points
//
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// includes
////////////////////////////////////////////////////////////////////////////////

#include "common.h"       // has struct FORWARDSTAR, NEIGHBOR, STATE
#include "sweep.h"        // has some version of do_sweep( STATE* )


////////////////////////////////////////////////////////////////////////////////
// constants
////////////////////////////////////////////////////////////////////////////////


// NOTE: comment-out this next line to sweep until convergence
#define MAXSWEEPS 1


#define	FSRADIUSMAX	7 // maximum radius forward star
#define FSMAX 818     // maximum number of points in a forward star
#define FSDELTA 10.0  // distance / delay multiplier
#define STARTMAX 12   // maximum number of starting points
#define GHOSTDEPTH FSRADIUSMAX


////////////////////////////////////////////////////////////////////////////////
// box file signatures
////////////////////////////////////////////////////////////////////////////////

const char vbox_sig[4] = {'v', 'b', 'o', 'x'};


////////////////////////////////////////////////////////////////////////////////
// prototypes
// note: prefix do_* because of name conflicts with dynamic libraries
////////////////////////////////////////////////////////////////////////////////

void do_freempibuffers( struct STATE *state );
void do_getargs( struct STATE *state, int argc, char *argv[] );
void do_initmpi( struct STATE *state, int argc, char *argv[] );
void do_initstate( struct STATE *state );
void do_loaddatafromfiles( struct STATE *state );
void do_preparempibuffers( struct STATE *state );
void do_preparestar( struct STATE *state );
void do_preparettbox( struct STATE *state );
void do_sharempibuffers( struct STATE *state );
void do_shutdown( struct STATE *state );
void do_workloop( struct STATE *state );
void do_writetttofile( struct STATE *state );


////////////////////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////////////////////

int
main (
  int argc,
  char *argv[]
)
{
  // all the state we care about is in here
  struct STATE state;

  // all MPI ranks do these functions equally
  do_initstate( &state );
  do_initmpi( &state, argc, argv );
  do_getargs( &state, argc, argv );

  // show some OpenMP information
  printf( "%d: openMP max threads: %d\n",
    state.myrank, omp_get_max_threads() );

  // read forward star
  do_preparestar( &state );

  // these functions do different things depending on this rank
  do_loaddatafromfiles( &state );
  do_preparempibuffers( &state );

  // TODO: read start position from somewhere
  state.ttstart = p3d( 152, 20, 1 );

  // ttbox == travel time FLOATBOX: includes ghost regions
  do_preparettbox( &state );

  // will stay in here for a while
  do_workloop( &state );

  // don't need these anymore
  do_freempibuffers( &state );

  // save travel time volume from this node
  do_writetttofile( &state );

  // free buffers and shutdown MPI and such
  do_shutdown( &state );

  // success
  return 0;
}


////////////////////////////////////////////////////////////////////////////////
// other functions
////////////////////////////////////////////////////////////////////////////////

void 
do_freempibuffers (
  struct STATE *state
)
{
  for( int n = 0; n < state->numneighbors; n++ ) {
    struct NEIGHBOR *neighbor = state->neighbors + n;
    boxfree( &neighbor->send );
    boxfree( &neighbor->recv );
  }
  state->numneighbors = 0;
}


void
do_getargs (
  struct STATE *state,
  int argc,
  char *argv[]
)
{
  struct ARGS *args = &state->args;

  // parse program arguments
  int goodargs = 0;
  if( argc == 5 ) {
    args->velocityfilename = argv[1];
    args->startpointsfilename = argv[2];
    args->forwardstarfilename = argv[3];
    args->traveltimefilenameprefix = argv[4];

    goodargs = 1;
  }

  // tell user how to use this program
  if( !goodargs ) {
    if( state->myrank == 0 ) {
      printf(
        "usage: %s"
        " <in:velocity.vbox> <in:startpoints.txt> <in:forwardstar.txt>"
        " <out:traveltimefilenameprefix>\n",
        argv[0]
      );
      fflush( stdout );
    }
    do_shutdown( state );
  }
}


void
do_initmpi (
  struct STATE *state,
  int argc,
  char *argv[]
)
{
  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &state->numranks );
  MPI_Comm_rank( MPI_COMM_WORLD, &state->myrank );

  if( !state->myrank ) {
    printf( "MPI ranks: %d\n", state->numranks );
  }
  
  int bestx = splitsquare_numx( state->numranks );
  int besty = state->numranks / bestx;
  int bestz = 1;

  state->rankdims = p3d( bestx, besty, bestz );
  state->rankcoords = mpifindrankcoordsfromrank( state->rankdims, state->myrank );
}


void
do_initstate (
  struct STATE *state
)
{
  boxinit( &state->vbox );
  boxinit( &state->ttbox );
  state->numneighbors = 0;
}


void
do_loaddatafromfiles (
  struct STATE *state
)
{
  // open file
  struct BOXOPENFILE vboxfile;
  if( !boxfileopenbinary( &vboxfile, state->args.velocityfilename, vbox_sig ) ) {
    fprintf (
      stderr,
      "%d: boxfileopenbinary(.., %s, ..) failed!\n",
      state->myrank, state->args.velocityfilename
    );
    fflush( stdout );
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }

  // store file coordinates
  state->gmin = vboxfile.min;
  state->gmax = vboxfile.max;

  // compute inner coordinates for this rank
  struct POINT3D imin, imax;
  mpifindregionfromrankcoords (
    &imin, &imax,
    state->gmin, state->gmax,
    state->rankdims, state->rankcoords
  );

  // compute outer coordinates from this rank
  struct POINT3D omin, omax;
  p3dgrowinside( &omin, &omax, imin, imax, GHOSTDEPTH, state->gmin, state->gmax );

  // load volume of outer coordinates from velocity file
  if( !boxfileloadbinarysubset( &state->vbox, omin, omax, imin, imax, vboxfile ) ) {
    fprintf (
      stderr,
      "%d: boxfileloadbinarysubset(..) from %s failed!\n",
      state->myrank, state->args.velocityfilename
    );
    fflush( stdout );
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }
  printf(
    "%d: loaded my region from velocity file: (%d, %d, %d) to (%d, %d, %d)\n",
    state->myrank, imin.x, imin.y, imin.z, imax.x, imax.y, imax.z
  );

  boxfileclosebinary( &vboxfile );
}


void
do_sharempibuffers (
  struct STATE *state
)
{
  const struct FLOATBOX ttbox = state->ttbox;

  MPI_Request mpireqs[26 + 26]; // at most 26 neighbors send + recv
  MPI_Status mpistats[26 + 26];
  int reqnumber = 0;

  // copy and send
  for( int n = 0; n < state->numneighbors; n++ ) {
    const struct FLOATBOX send = state->neighbors[n].send;
    boxcopysubset( send, ttbox, send.omin, send.omax );
    MPI_Isend (
      send.flat, send.flat_size, MPI_FLOAT,
      state->neighbors[n].rank, state->myrank,
      MPI_COMM_WORLD, &mpireqs[ reqnumber++ ]
    );
  }

  // receive (don't block)
  for( int n = 0; n < state->numneighbors; n++ ) {
    const struct FLOATBOX recv = state->neighbors[n].recv;
    MPI_Irecv (
      recv.flat, recv.flat_size, MPI_FLOAT,
      state->neighbors[n].rank, state->neighbors[n].rank,
      MPI_COMM_WORLD, &mpireqs[ reqnumber++ ]
    );
  }

  MPI_Waitall( reqnumber, mpireqs, mpistats );

  // copy from recv buffers back to ttbox (ghost regions)
  for( int n = 0; n < state->numneighbors; n++ ) {
    const struct FLOATBOX recv = state->neighbors[n].recv;
    boxcopysubset( ttbox, recv, recv.omin, recv.omax );
  }
}


void
do_preparempibuffers (
  struct STATE *state
)
// allocates send and receive buffers for adjacent neighbors
{
  // start with 0 neighbors, and increment as they are discovered
  state->numneighbors = 0;

  // find all neighbors
  for( int x = -1; x <= 1; x++ ) {
    for( int y = -1; y <= 1; y++ ) {
      for( int z = -1; z <= 1; z++ ) {
        if( x || y || z ) { // (0,0,0) is me

          struct POINT3D relation = p3d( x, y, z );
          struct POINT3D ncoords = p3daddp3d( state->rankcoords, relation );

          int nrank = mpifindrankfromrankcoords( state->rankdims, ncoords );
          if( nrank >= 0 ) {

            struct NEIGHBOR *neighbor = state->neighbors + state->numneighbors++;

            neighbor->relation = relation; 
            neighbor->rank = nrank;

            // neighbor's inner region (no ghost)
            struct POINT3D nimin, nimax;
            mpifindregionfromrankcoords (
              &nimin, &nimax,
              state->gmin, state->gmax,
              state->rankdims, ncoords
            );

            // neighbor's outer region (including ghost)
            struct POINT3D nomin, nomax;
            p3dgrowinside (
              &nomin, &nomax, nimin, nimax,
              GHOSTDEPTH,
              state->gmin, state->gmax
            );

            // send intersection: neighbor's ghost overlapping my inner region
            struct POINT3D smin, smax;
            // recv intersection: my ghost with neighbor's inner region
            struct POINT3D rmin, rmax;

            // compute intersections
            if (
              !intersect3d (
                &smin, &smax, nomin, nomax, state->vbox.imin, state->vbox.imax
              ) ||
              !intersect3d (
                &rmin, &rmax, nimin, nimax, state->vbox.omin, state->vbox.omax
              )
            ) {
              // intersection failed for some reason (this shouldn't happen)
              fprintf (
                stderr,
                "%d: error: no intersection with neighbor %d, but there should be!\n",
                state->myrank, nrank
              );
            }
            else {
              // allocate memory for send and recv buffers for this neighbor
              if (
                !boxalloc( &neighbor->send, smin, smax, smin, smax ) ||
                !boxalloc( &neighbor->recv, rmin, rmax, rmin, rmax )
              ) {
                fprintf (
                  stderr,
                  "%d: error: memory allocation failure when preparing ghost buffers:\n"
                  "%d:        my rank coords: (%d, %d, %d)\n"
                  "%d:        neighbor %d rank coords: (%d, %d, %d)\n",
                  state->myrank,
                  state->myrank,
                  state->rankcoords.x, state->rankcoords.y, state->rankcoords.z,
                  state->myrank, neighbor->rank, ncoords.x, ncoords.y, ncoords.z
                );
                fflush( stdout );
                MPI_Abort( MPI_COMM_WORLD, 1 );
              }
            }
          }
        }
      }
    }
  }
}


void
do_preparestar (
  struct STATE *state
)
{
  FILE *infile;
  if( NULL != (infile = fopen( state->args.forwardstarfilename, "r" )) ) {

    int starsize = 0;
    if( 1 == fscanf( infile, "%d", &starsize ) && starsize > 0 ) {

      struct FORWARDSTAR *star = malloc( starsize * sizeof(struct FORWARDSTAR) );

      int bad = 0;
      for( int i = 0; i < starsize; i++ ) {
        struct POINT3D pos;

        if( 3 == fscanf( infile, "%d %d %d", &pos.x, &pos.y, &pos.z ) ) {
          star[i].pos = pos;
          star[i].halfdistance = FSDELTA * 0.5 * sqrt (
            pos.x * pos.x + pos.y * pos.y + pos.z * pos.z
          );
        }
        else {
          bad = 1;
          break;
        }
      }
      if( !bad ) {
        fclose( infile );
        state->star = star;
        state->numinstar = starsize;
        return;
      }
    }
    // read error: probably wrong file or bad format or something
    fprintf (
      stderr, "%d: error parsing forward star file %d\n",
      state->myrank, state->args.forwardstarfilename
    );
  }
  else { // fopen failed
    fprintf (
      stderr, "%d: failed to open forward star file: %s\n",
      state->myrank, state->args.forwardstarfilename
    );
  }

  MPI_Abort( MPI_COMM_WORLD, 1 );
}


void
do_preparettbox (
  struct STATE *state
)
{
  struct POINT3D omin = state->vbox.omin;
  struct POINT3D omax = state->vbox.omax;
  struct POINT3D imin = state->vbox.imin;
  struct POINT3D imax = state->vbox.imax;

  // allocate memory
  if( !boxalloc( &state->ttbox, omin, omax, imin, imax ) ) {
    fprintf (
      stderr,
      "%d: error: memory allocation failure for travel time volume:\n"
      "%d:        outer: (%d, %d, %d) to (%d, %d, %d)\n"
      "%d:        inner: (%d, %d, %d) to (%d, %d, %d)\n",
      state->myrank,
      state->myrank, omin.x, omin.y, omin.z, omax.x, omax.y, omax.z,
      state->myrank, imin.x, imin.y, imin.z, imax.x, imax.y, imax.z
    );
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }

  // set all values to INFINITY
  boxsetall( state->ttbox, INFINITY );

  // if a start point is in this box, set it to 0
  if (
    !p3disless( state->ttstart, state->ttbox.imin ) &&
    !p3dismore( state->ttstart, state->ttbox.imax )
  ) {
    boxputglobal( state->ttbox, state->ttstart, 0.f );
  }
}


void
do_shutdown (
  struct STATE *state
)
{
  free( state->star );
  do_freempibuffers( state );
  boxfree( &state->ttbox );
  boxfree( &state->vbox );
  printf( "%d: shutting down normally\n", state->myrank );
  fflush( stdout );
  MPI_Barrier( MPI_COMM_WORLD );
  MPI_Finalize();
  exit(0);
}


void
do_workloop (
  struct STATE *state
)
{
  state->numsweeps = 0;

  for(;;) { // infinite loop

#ifdef MAXSWEEPS
    if( state->numsweeps >= MAXSWEEPS ) return;
#endif

    state->numsweeps++;
    printf( "%d: doing sweep %ld...\n", state->myrank, state->numsweeps );
    long sweepchanges = do_sweep( state );
    printf (
      "%d: sweep %ld: number of changes: %ld\n",
      state->myrank, state->numsweeps, sweepchanges
    );

    if( state->numranks == 0 ) {
      // only one rank total, so don't bother with MPI
      if( sweepchanges < 1 ) break;
    }
    else {
      // more than one rank: need to share information

      char anychange = 0;
      if( state->myrank == 0 ) {
        // rank 0 collects and shares global change status

        // start with my own sweep changes
        anychange = sweepchanges > 0 ? 1 : 0;

        // collect others' changes
        for( int r = 1; r < state->numranks; r++ ) {
          char otherchange = 0;
          MPI_Recv( &otherchange, 1, MPI_BYTE, r, r, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
          anychange |= otherchange;
        }
      }
      else {
        // this is NOT rank 0
        char mychange = sweepchanges > 0 ? 1 : 0;
        MPI_Send( &mychange, 1, MPI_BYTE, 0, state->myrank, MPI_COMM_WORLD );
      }

      // share anychange from rank 0 to everyone else
      MPI_Bcast( &anychange, 1, MPI_BYTE, 0, MPI_COMM_WORLD );

      if( !anychange ) return;

      do_sharempibuffers( state );
    }
  }
}


void
do_writetttofile (
  struct STATE *state
)
// each MPI node writes its region to its own file
{
  const struct FLOATBOX ttbox = state->ttbox;
  const struct POINT3D min = ttbox.imin;
  const struct POINT3D max = ttbox.imax;

  char filename[1024];
  sprintf (
    filename,
    "%s_%d-%d-%d_to_%d-%d-%d.txt",
    state->args.traveltimefilenameprefix,
    min.x, min.y, min.z,
    max.x, max.y, max.z
  );

  FILE *outfile = fopen( filename, "w" );
  if( outfile == NULL ) {
    fprintf (
      stderr, "%d: error: can't open output file %s!\n",
      state->myrank, filename
    );
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }

  for( int x = min.x; x <= max.x; x++ ) {
    for( int y = min.y; y <= max.y; y++ ) {
      for( int z = min.z; z <= max.z; z++ ) {
        fprintf (
          outfile, "%d,%d,%d,%g\n",
          x, y, z, boxgetglobal( ttbox, p3d( x, y, z ) )
        );
      }
    }
  }

  fclose( outfile );
}


////////////////////////////////////////////////////////////////////////////////
// vim: set tabstop=2 shiftwidth=2 softtabstop=2 expandtab:
// END
////////////////////////////////////////////////////////////////////////////////
