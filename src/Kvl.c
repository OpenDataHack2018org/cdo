/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2012 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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
  void *kvlist;
  int nlists, listID;
  int nelements, elemID;
  const char *listname;
  const char *ename;
  const char *evalue;

  kvlist = kvlParseFile(filename);
  nlists = kvlGetNumLists(kvlist);
  printf("# Number of lists: %d\n", nlists);
  for ( listID = 0; listID < nlists; ++listID )
    {
      listname = kvlGetListName(kvlist, listID);
      nelements = kvlGetListNumElements(kvlist, listID);
      printf("# list ID: %d;   Number of elements: %d\n", listID, nelements);
      printf("&%s\n", listname);
      for ( elemID = 0; elemID < nelements; ++elemID )
	{
	  ename  = kvlGetListElementName(kvlist, listID, elemID);
	  evalue = kvlGetListElementValue(kvlist, listID, elemID);
	  printf("  %s = %s\n", ename, evalue);
	}
      printf("/\n");
    }

  kvlDelete(kvlist);

  return (0);
}

static
int conv_cmor_table(const char *filename)
{
  void *kvlist;
  int nlists, listID;
  int nelements, elemID;
  int len;
  int hasmissval = FALSE;
  double missval;
  const char *listname;
  const char *ename;
  const char *evalue;

  kvlist = kvlParseFile(filename);
  nlists = kvlGetNumLists(kvlist);
  //printf("# Number of lists: %d\n", nlists);
  for ( listID = 0; listID < nlists; ++listID )
    {
      listname = kvlGetListName(kvlist, listID);
      nelements = kvlGetListNumElements(kvlist, listID);
      //printf("# list ID: %d;   Number of elements: %d\n", listID, nelements);
      if ( strncmp("global", listname, strlen(listname)) == 0 )
	{
	  for ( elemID = 0; elemID < nelements; ++elemID )
	    {
	      ename  = kvlGetListElementName(kvlist, listID, elemID);
	      evalue = kvlGetListElementValue(kvlist, listID, elemID);
	      len = strlen(ename);

	      if ( strncmp("missing_value", ename, len) == 0 )
		{
		  missval = atof(evalue);
		  hasmissval = TRUE;
		}
	    }
	}
      else if ( strncmp("variable", listname, strlen(listname)) == 0 )
	{
	  printf("&%s\n", "parameter");
	  for ( elemID = 0; elemID < nelements; ++elemID )
	    {
	      ename  = kvlGetListElementName(kvlist, listID, elemID);
	      evalue = kvlGetListElementValue(kvlist, listID, elemID);
	      len = strlen(ename);

	      if ( strncmp("name", ename, len)            == 0 ||
		   strncmp("standard_name", ename, len)   == 0 ||
		   strncmp("out_name", ename, len)        == 0 ||
		   strncmp("type", ename, len)            == 0 ||
		   strncmp("valid_min", ename, len)       == 0 ||
		   strncmp("valid_max", ename, len)       == 0 ||
		   strncmp("ok_min_mean_abs", ename, len) == 0 ||
		   strncmp("ok_max_mean_abs", ename, len) == 0 )
		printf("  %-15s = %s\n", ename, evalue);
	      else if ( strncmp("long_name", ename, len)  == 0 ||
		   strncmp("units", ename, len)           == 0 ||
		   strncmp("cell_methods", ename, len)    == 0 ||
		   strncmp("cell_measures", ename, len)   == 0 ||
		   strncmp("comment", ename, len)         == 0 )
		printf("  %-15s = \"%s\"\n", ename, evalue);
	    }
	  if ( hasmissval ) printf("  %-15s = %g\n", "missing_value", missval);
	  printf("/\n");
	}
    }

  kvlDelete(kvlist);

  return (0);
}


void *Kvl(void *argument)
{
  int READ_CMOR_TABLE, CONV_CMOR_TABLE;
  int operatorID;
  const char *filename;

  cdoInitialize(argument);

  READ_CMOR_TABLE = cdoOperatorAdd("read_cmor_table",   0,   0, NULL);
  CONV_CMOR_TABLE = cdoOperatorAdd("conv_cmor_table",   0,   0, NULL);

  operatorID = cdoOperatorID();

  if ( operatorID == READ_CMOR_TABLE )
    {
      if ( operatorArgc() != 1 ) cdoAbort("Too few arguments!");
      filename = operatorArgv()[0];

      if ( cdoVerbose ) cdoPrint("Parse file: %s\n", filename);

      read_cmor_table(filename);
    }
  else if ( operatorID == CONV_CMOR_TABLE )
    {
      if ( operatorArgc() != 1 ) cdoAbort("Too few arguments!");
      filename = operatorArgv()[0];

      if ( cdoVerbose ) cdoPrint("Parse file: %s\n", filename);

      conv_cmor_table(filename);
    }

  cdoFinish();

  return (0);
}
