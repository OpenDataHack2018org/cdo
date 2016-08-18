/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2016 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

/*
   This module contains the following operators:

*/

#include "cdo.h"
#include "cdo_int.h"
#include "kvlist.h"


static
int read_cmor_table(const char *filename)
{
  void *kvlist = kvlParseFile(filename);
  int nlists = kvlGetNumLists(kvlist);
  printf("# Number of lists: %d\n", nlists);
  for ( int listID = 0; listID < nlists; ++listID )
    {
      const char *listname = kvlGetListName(kvlist, listID);
      int nelements = kvlGetListNumElements(kvlist, listID);
      printf("# list ID: %d;   Number of elements: %d\n", listID, nelements);
      printf("&%s\n", listname);
      for ( int elemID = 0; elemID < nelements; ++elemID )
	{
	  const char *ename  = kvlGetListElementName(kvlist, listID, elemID);
	  const char *evalue = kvlGetListElementValue(kvlist, listID, elemID);
	  printf("  %s = %s\n", ename, evalue);
	}
      printf("/\n");
    }

  kvlDelete(kvlist);

  return 0;
}

static
int conv_cmor_table(const char *filename)
{
  bool hasmissval = false;
  double missval;

  void *kvlist = kvlParseFile(filename);
  int nlists = kvlGetNumLists(kvlist);
  //printf("# Number of lists: %d\n", nlists);
  for ( int listID = 0; listID < nlists; ++listID )
    {
      const char *listname = kvlGetListName(kvlist, listID);
      int nelements = kvlGetListNumElements(kvlist, listID);
      //printf("# list ID: %d;   Number of elements: %d\n", listID, nelements);
      if ( strncmp("global", listname, strlen(listname)) == 0 )
	{
	  for ( int elemID = 0; elemID < nelements; ++elemID )
	    {
	      const char *ename  = kvlGetListElementName(kvlist, listID, elemID);
	      const char *evalue = kvlGetListElementValue(kvlist, listID, elemID);
              size_t len = strlen(ename);

	      if ( strncmp("missing_value", ename, len) == 0 )
		{
		  missval = atof(evalue);
		  hasmissval = true;
		}
	    }
	}
      else if ( strncmp("variable", listname, strlen(listname)) == 0 )
	{
	  int vlen;
	  printf("&%s\n", "parameter");
	  for ( int elemID = 0; elemID < nelements; ++elemID )
	    {
	      const char *ename  = kvlGetListElementName(kvlist, listID, elemID);
	      const char *evalue = kvlGetListElementValue(kvlist, listID, elemID);
	      int len = strlen(ename);
	      vlen = strlen(evalue);

	      if ( vlen > 1 && evalue[0] == '"' && evalue[vlen-1] == '"' ) 
		{
		  vlen -= 2;
		  evalue++;
		}

	      char *ovalue = strdup(evalue);
	      for ( int i = 1; i < vlen; ++i )
		{
		  if ( ovalue[i-1] == '"' && ovalue[i] == '"' )
		    {
		      ovalue [i-1] = '\'';
		      for ( int j = i+1; j < vlen; ++j ) ovalue[j-1] = ovalue[j];
		      vlen -= 1;
		    }
		}

	      if ( strncmp("name", ename, len)            == 0 ||
		   strncmp("standard_name", ename, len)   == 0 ||
		   strncmp("out_name", ename, len)        == 0 ||
		   strncmp("type", ename, len)            == 0 ||
		   strncmp("valid_min", ename, len)       == 0 ||
		   strncmp("valid_max", ename, len)       == 0 ||
		   strncmp("ok_min_mean_abs", ename, len) == 0 ||
		   strncmp("ok_max_mean_abs", ename, len) == 0 )
		printf("  %-15s = %s\n", ename, ovalue);
	      else if ( strncmp("long_name", ename, len)  == 0 ||
		   strncmp("units", ename, len)           == 0 ||
		   strncmp("cell_methods", ename, len)    == 0 ||
		   strncmp("cell_measures", ename, len)   == 0 ||
		   strncmp("comment", ename, len)         == 0 )
		printf("  %-15s = \"%.*s\"\n", ename, vlen, ovalue);

	      Free(ovalue);
	    }
	  if ( hasmissval ) printf("  %-15s = %g\n", "missing_value", missval);
	  printf("/\n");
	}
    }

  kvlDelete(kvlist);

  return 0;
}


void *Kvl(void *argument)
{
  cdoInitialize(argument);

  int READ_CMOR_TABLE = cdoOperatorAdd("read_cmor_table",   0,   0, NULL);
  int CONV_CMOR_TABLE = cdoOperatorAdd("conv_cmor_table",   0,   0, NULL);
  int CONV_PARTAB     = cdoOperatorAdd("conv_partab",   0,   0, NULL);

  int operatorID = cdoOperatorID();

  if ( operatorID == READ_CMOR_TABLE )
    {
      if ( operatorArgc() != 1 ) cdoAbort("Too few arguments!");
      const char *filename = operatorArgv()[0];

      if ( cdoVerbose ) cdoPrint("Parse file: %s", filename);

      read_cmor_table(filename);
    }
  else if ( operatorID == CONV_CMOR_TABLE )
    {
      if ( operatorArgc() != 1 ) cdoAbort("Too few arguments!");
      const char *filename = operatorArgv()[0];

      if ( cdoVerbose ) cdoPrint("Parse file: %s", filename);

      conv_cmor_table(filename);
    }
  else if ( operatorID == CONV_PARTAB )
    {
      if ( operatorArgc() != 1 ) cdoAbort("Too few arguments!");
      const char *filename = operatorArgv()[0];

      if ( cdoVerbose ) cdoPrint("Parse file: %s", filename);

      // conv_partab(filename);
    }

  cdoFinish();

  return 0;
}
