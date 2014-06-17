/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2014 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include <stdio.h>
#include <string.h>

int CDO_optind = 1;
char *CDO_optarg;

int cdo_getopt(int argc, char * const argv[], const char *optstring)
{
  static int optpos = 0;
  int optval = -1, value;
  int opthasarg = 0;
  int optstrlen = strlen(optstring);
  int iargc;

  CDO_optarg = NULL;

  while ( optpos < optstrlen && CDO_optind < argc )
    {
      value = optstring[optpos];
      optpos++;
      if ( optstring[optpos] == ':' )
	{
	  opthasarg = 1;
	  optpos++;
	}
      else
	opthasarg = 0;

      for ( iargc = 1; iargc < argc; iargc++ )
	{
	  if ( *argv[iargc] == '-' && strlen(argv[iargc]) == 2 )
	    {
	      if ( (argv[iargc][1]) == value )
		{
		  optval = value;
		  CDO_optind++;
		  if ( opthasarg )
		    {
		      CDO_optarg = argv[iargc+1];
		      CDO_optind++;
		    }
		  break;
		}
	    }
	}
      if ( iargc < argc ) break;
    }

  if ( opthasarg && CDO_optarg == NULL ) optval = ':';

  return (optval);
}
