////////////////////////////////////////////////////////////////////////////////
// sweep.h
////////////////////////////////////////////////////////////////////////////////
//
// Functions:
//   do_sweep( struct STATE *state )
//
//
////////////////////////////////////////////////////////////////////////////////

#pragma once


////////////////////////////////////////////////////////////////////////////////
// includes
////////////////////////////////////////////////////////////////////////////////

#include "common.h"


////////////////////////////////////////////////////////////////////////////////
// compile options
////////////////////////////////////////////////////////////////////////////////

// beware: multiple OMP threads will result in nondeterministic convergence
//         due to races when changing travel times.
#define SWEEP_OMP_MAX_THREADS 1


////////////////////////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////////////////////////

long
do_sweep (
  struct STATE *state  
)
{
  // copy some state into local stack memory for fastness
  const struct FLOATBOX vbox = state->vbox;
  const struct FLOATBOX ttbox = state->ttbox;
  const struct FORWARDSTAR * const star = state->star;
  const int numinstar = state->numinstar;

  // count how many (if any) values we change
  long changes = 0;

  //#pragma omp parallel for default(shared) \
  //  reduction(+:changes) schedule(dynamic) num_threads(SWEEP_OMP_MAX_THREADS)
  for( int x = vbox.imin.x; x <= vbox.imax.x; x++ ) {

    const int minx = vbox.omin.x - x;
    const int maxx = vbox.omax.x - x;

    for( int y = vbox.imin.y; y <= vbox.imax.y; y++ ) {

      const int miny = vbox.omin.y - y;
      const int maxy = vbox.omax.y - y;

      for( int z = vbox.imin.z; z <= vbox.imax.z; z++ ) {

        const int minz = vbox.omin.z - z;
        const int maxz = vbox.omax.z - z;

        const long index_here = boxindex( ttbox, p3d( x, y, z ) ) - ttbox.index_omin;

        const float vel_here = vbox.flat[ index_here ];

        float new_tt_here = ttbox.flat[ index_here ];

        for( int l = 0; l < numinstar; l++ ) {

          // if 'there' is outside the boundaries, then skip
          if(
            star[l].pos.x < minx ||
            star[l].pos.x > maxx ||
            star[l].pos.y < miny ||
            star[l].pos.y > maxy ||
            star[l].pos.z < minz ||
            star[l].pos.z > maxz
          ) {
            continue;
          }

          const long index_there = index_here + star[l].offset;
          
          // compute delay from 'here' to 'there' with endpoint average
          const float delay = star[l].halfdistance *
                              (vel_here + vbox.flat[ index_there ]);

          // current travel time from start point to 'there'
          const float tt_there = ttbox.flat[ index_there ];

          // update (maybe) tt_there with a better time
          if( new_tt_here + delay < tt_there ) {
            changes++;
            ttbox.flat[ index_there ] = new_tt_here + delay;
          }

          // update (maybe) new_tt_here with a potentially better time
          if( tt_there + delay < new_tt_here ) new_tt_here = tt_there + delay;
        }

        // if a faster path was found, use it
        if( new_tt_here < ttbox.flat[ index_here ] ) {
          changes++;
          ttbox.flat[ index_here ] = new_tt_here;
        }
      }
    }
  }

  return changes;
}


////////////////////////////////////////////////////////////////////////////////
// vim: set tabstop=2 shiftwidth=2 softtabstop=2 expandtab:
// END
////////////////////////////////////////////////////////////////////////////////
