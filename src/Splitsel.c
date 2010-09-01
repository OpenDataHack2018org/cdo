/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2010 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Splitsel   splitsel        Split time selection
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


void *Splitsel(void *argument)
{
  static const char *func = "Splitsel";

  /* from Selstat.c */
  int operatorID;
  int operfunc;
  int gridsize;
  int vdate = 0, vtime = 0;
  int nrecs = 0;
  int varID, levelID, recID;
  int tsID, tsID2;
  int nsets;
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nmiss;
/*   int ndates = 0, noffset = 0, nskip = 0, nargc; */
  double ndates, noffset, nskip;
  int i2 = 0;
  int nargc;

  /* from Splittime.c */
  int nchars;
  char *filesuffix;
  char filename[1024];
  int index = 0;
  int lcopy = FALSE;
  double *array = NULL;

  cdoInitialize(argument);

  cdoOperatorAdd("splitsel",  0,  0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  if ( UNCHANGED_RECORD ) lcopy = TRUE;

  /*  operatorInputArg("nsets <noffset <nskip>>"); */

  nargc = operatorArgc();

  if ( nargc < 1 )
    cdoAbort("Too few arguments! Need %d found %d.", 1, nargc);

/*   ndates = atoi(operatorArgv()[0]); */
/*   if ( nargc > 1 ) noffset = atoi(operatorArgv()[1]); */
/*   if ( nargc > 2 ) nskip   = atoi(operatorArgv()[2]); */
/*   printf("%s %s %s\n", operatorArgv()[0],operatorArgv()[1],operatorArgv()[2]); */
  ndates = noffset = nskip = 0.0;
  ndates = atof(operatorArgv()[0]);
  if ( nargc > 1 ) noffset = atof(operatorArgv()[1]);
  if ( nargc > 2 ) nskip   = atof(operatorArgv()[2]);

/*   if ( cdoVerbose ) cdoPrint("nsets = %d, noffset = %d, nskip = %d", ndates, noffset, nskip); */
  if ( cdoVerbose ) cdoPrint("nsets = %f, noffset = %f, nskip = %f", ndates, noffset, nskip);

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
/*   taxisID2 = taxisCreate(TAXIS_ABSOLUTE); */
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  strcpy(filename, cdoStreamName(1));
  nchars = strlen(filename);
  filesuffix = streamFilesuffix(cdoDefaultFileType);

  if ( ! lcopy )
    {
      gridsize = vlistGridsizeMax(vlistID1);
      if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
      array = (double *) malloc(gridsize*sizeof(double));
    }


  /* offset */
  for ( tsID = 0; tsID < noffset; tsID++ )
    {
      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) { /* printf ("break\n");*/ break; }
    }
  if ( tsID < noffset )
    {
      cdoWarning("noffset larger than number of timesteps!");
      goto LABEL_END;
    }

  index = 0;
  nsets = 0;
  while ( TRUE )
    {
      sprintf(filename+nchars, "%06d", index);
      sprintf(filename+nchars+6, "%s", filesuffix);
	  
      if ( cdoVerbose ) cdoPrint("create file %s", filename);
      streamID2 = streamOpenWrite(filename, cdoFiletype());
      if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", filename);

      streamDefVlist(streamID2, vlistID2);

      tsID2 = 0;

/*       printf("comp: %f %d\n",ndates*(index+1),(int)(ndates*(index+1))); */
      for ( ; nsets < (int)(ndates*(index+1)); nsets++ ) 
	{

	  nrecs = streamInqTimestep(streamID1, tsID);
	  if ( nrecs == 0 ) break;

	  vdate = taxisInqVdate(taxisID1);
	  vtime = taxisInqVtime(taxisID1);

	  /* printf("vdate: %d vtime: %d\n", vdate, vtime); */

	  taxisCopyTimestep(taxisID2, taxisID1);
	  streamDefTimestep(streamID2, tsID2);

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      
	      streamInqRecord(streamID1, &varID, &levelID);
	      streamDefRecord(streamID2,  varID,  levelID);
	      if ( lcopy )
		{
		  streamCopyRecord(streamID2, streamID1);
		}
	      else
		{
		  streamReadRecord(streamID1, array, &nmiss);
		  streamWriteRecord(streamID2, array, nmiss);
		}
	    }
	  
	  tsID++;
	  tsID2++;	  
	}
      
      streamClose(streamID2);
      if ( nrecs == 0 ) break;

      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;

/*       for ( i = 0; i < nskip; i++ ) */
      for ( ; i2 < (int)(nskip*(index+1)); i2++ )
	{
	  nrecs = streamInqTimestep(streamID1, tsID);
	  if ( nrecs == 0 ) break;
	  tsID++;
	}

      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;

      index++;
    }

 LABEL_END:

  streamClose(streamID1);
 
  if ( ! lcopy )
    if ( array ) free(array);

  cdoFinish();

  return (0);
}
