#include <math.h>
#include <stdio.h>
#include <stdlib.h>


struct CHARCOORDS {
    char x, y, z;
};


long
cube (
    int side
)
{
    return (long)side*side*side;
}


int
main (
    int argc,
    char *argv[]
)
{
    // get arguments from command line
    int radius;

    if (
        3 != argc ||
        1 != sscanf( argv[1], "%d", &radius ) ||
        radius < 1
    ) {
        printf( "usage: %s <radius> <outfile>\n", argv[0] );
        return 0;
    }

    char *outfilename = argv[2];

    // build spherical star in memory
    struct CHARCOORDS *coords = malloc( sizeof(char) * cube( 2*radius + 1 ) );

    if( coords == NULL ) {
        fprintf( stderr, "couldn't allocate enough memory for radius %d\n", radius );
        return 1;
    }

    long numcoords = 0;
    
    const float maxrsquared = (float)radius * radius;

    for( int x = -radius; x <= radius; x++ ) {
        for( int y = -radius; y <= radius; y++ ) {
            for( int z = -radius; z <= radius; z++ ) {
                const float rsquared = (float)x*x + (float)y*y + (float)z*z;
                if( rsquared <= maxrsquared ) {
                    coords[ numcoords++ ] = (struct CHARCOORDS){ x, y, z };
                }
            }
        }
    }

    // write star to file
    FILE *outfile = fopen( outfilename, "w" );
    if( outfile == NULL ) {
        fprintf( stderr, "couldn't open/create file %s\n", outfilename );
        free( coords );
        return 1;
    }

    fprintf( outfile, "%ld\n", numcoords );

    for( long i = 0; i < numcoords; i++ ) {
        fprintf( outfile, "%d %d %d\n",
            (int)coords[i].x, (int)coords[i].y, (int)coords[i].z );
    }

    fclose( outfile );

    free( coords );
}
