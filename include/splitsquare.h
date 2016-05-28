////////////////////////////////////////////////////////////////////////////////
// splitsquare.h - 2015.05.27 - Atlee Brink
////////////////////////////////////////////////////////////////////////////////
//
// Functions:
//   splitsquare_numx
//
//
// Example:
//   int num_nodes = 9;
//   int num_x = splitsquare_numx( num_nodes );
//   int num_y = num_nodes / num_x;
//
//
////////////////////////////////////////////////////////////////////////////////

#pragma once


////////////////////////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////////////////////////

int
splitsquare_numx (
    const int num
)
// returns: nx such that nx*(num/nx) minimizes |(size/nx)-(size/(num/nx))|
//          for a square of 'size'
{
    int last = 1;
    for( int cur = 2; cur*cur <= num; cur++ ) {
        if( cur*cur == num ) return cur;
        if( num % cur == 0 ) last = cur;
    }
    return num/last;
}


////////////////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////////////////
