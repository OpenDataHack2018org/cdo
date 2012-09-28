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


/*
int main(int argc, char *argv[])
{
  char *filename;
  void *kvlist;
  int nlists, listID;
  int nelements, elemID;
  const char *listname;
  const char *ename;
  const char *evalue;

  if ( argc != 2 ) 
    {
      fprintf(stderr, "usage: kvlist filename\n");
      return (1);
    }

  filename = argv[1];

  printf("Parse file: %s\n", filename);

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
*/


void *Kvl(void *argument)
{

  cdoInitialize(argument);


  cdoFinish();

  return (0);
}
