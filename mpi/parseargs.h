////////////////////////////////////////////////////////////////////////////////
// parseargs.h - 2016.05.31 - Atlee Brink
////////////////////////////////////////////////////////////////////////////////

#pragma once


////////////////////////////////////////////////////////////////////////////////
// includes
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>


////////////////////////////////////////////////////////////////////////////////
// structures
////////////////////////////////////////////////////////////////////////////////

struct ARGS {
  const char *velocityfilename;
  const char *startpointsfilename;
  const char *forwardstarfilename;
  const char *traveltimefilename;
};


////////////////////////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////////////////////////

int
parseargs (
  struct ARGS *args,
  int argc,
  char *argv[]
)
// parses command line arguments and stores them in 'args'
// error: returns 0
// success: returns non-0
{
  // check arguments
  if( argc != 5 ) return 0;

  // filenames
  args->velocityfilename = argv[1];
  args->startpointsfilename = argv[2];
  args->forwardstarfilename = argv[3];
  args->traveltimefilename = argv[4];

  // success
  return 1;
}


////////////////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////////////////
