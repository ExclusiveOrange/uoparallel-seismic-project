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

  #pragma omp parallel for default(shared) \
    reduction(+:changes) schedule(dynamic) num_threads(16)
  for( int x = vbox.imin.x; x <= vbox.imax.x; x++ ) {
    for( int y = vbox.imin.y; y <= vbox.imax.y; y++ ) {
      for( int z = vbox.imin.z; z <= vbox.imax.z; z++ ) {

        const struct POINT3D here = p3d( x, y, z );

        const long index_here = boxindex( ttbox, here ) - ttbox.index_omin;

        const float vel_here = vbox.flat[ index_here ];
        const float tt_here = ttbox.flat[ index_here ];

        float new_tt_here = tt_here;

        for( int l = 0; l < numinstar; l++ ) {

          // find point in forward star based on offsets
          const struct POINT3D there = p3daddp3d( here, star[l].pos );

          // if 'there' is outside the boundaries, then skip
          if (
            p3disless( there, vbox.omin ) ||
            p3dismore( there, vbox.omax )
          ) {
            continue;
          }

          const long index_there = index_here + star[l].offset;
          
          // compute delay from 'here' to 'there' with endpoint average
          const float vel_there = vbox.flat[ index_there ];
          const float delay = star[l].halfdistance * (vel_here + vel_there);

          // current travel time from start point to 'there'
          const float tt_there = ttbox.flat[ index_there ];

          // update (maybe) tt_there with a better time
          const float maybe_new_tt_there = new_tt_here + delay;
          if( maybe_new_tt_there < tt_there ) {
            changes++;
            ttbox.flat[ index_there ] = maybe_new_tt_there;
          }

          // update (maybe) new_tt_here with a potentially better time
          new_tt_here = fmin( new_tt_here, tt_there + delay );
        }

        // if a faster path was found, use it
        if( new_tt_here < tt_here ) {
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
