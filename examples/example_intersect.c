#include "intersect.h"

#include <stdio.h>

int main() {

    int amin = 10, amax = 20;
    
    for( int x = -10; x <= 40; x++ ) {
        int bmin = x, bmax = x;
        int min, max;
        int doesintersect = intersect1d( &min, &max, amin, amax, bmin, bmax );
        printf( "[%d, %d] and [%d, %d]: ", amin, amax, bmin, bmax );
        if( doesintersect ) {
            printf( "[%d, %d]\n", min, max );
        }
        else {
            printf( "null\n" );
        }
    }
    return 0;
}
