/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2006 Uwe Schulzweida, schulzweida@dkrz.de
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

#include <string.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "dtypes.h"

  /* from Splittime.c */
//#define  MAX_STREAMS 32
#define  MAX_STREAMS 1000

void *Splitsel(void *argument)
{
  static char func[] = "Splitsel";

  /* from Selstat.c */
  int operatorID;
  int operfunc;
  int gridsize;
  int vdate = 0, vtime = 0;
  int vdate0 = 0, vtime0 = 0;
  int nrecs = 0, nrecords;
  int gridID, varID, levelID, recID;
  int tsID;
  int otsID;
  int nsets;
  int i;
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nmiss;
  int nvars, nlevel;
/*   int ndates = 0, noffset = 0, nskip = 0, nargc; */
  double ndates, noffset, nskip;
  int i2 = 0;
  int nargc;
  int *recVarID, *recLevelID;
  FIELD field;

  /* from Splittime.c */
  int nchars;
  int  streamIDs[MAX_STREAMS], tsIDs[MAX_STREAMS];
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

  for ( i = 0; i < MAX_STREAMS; i++ ) streamIDs[i] = -1;
  for ( i = 0; i < MAX_STREAMS; i++ ) tsIDs[i] = 0;

  //  operatorInputArg("nsets <noffset <nskip>>");

  nargc = operatorArgc();

  if (nargc<1)
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

  otsID = 0;
  index = 0;
  nsets = 0;
  while ( TRUE )
    {
/*       for ( nsets = 0; nsets < ndates; nsets++ ) */
/* 	{ */
/*       printf("index: %d tsID: %d nrecs: %d\n",index,tsID,nrecs); */
/*       nrecs = streamInqTimestep(streamID1, tsID); */
      //      if ( nrecs == 0 ) break;

      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

/*       printf("vdate: %d vtime: %d\n",vdate,vtime); */
      if ( index < 0 || index >= MAX_STREAMS )
	cdoAbort("Index out of range!");

      streamID2 = streamIDs[index];
      if ( streamID2 < 0 )
	{
	  sprintf(filename+nchars, "%03d", index);
	  sprintf(filename+nchars+3, "%s", filesuffix);
	  
	  if ( cdoVerbose ) cdoPrint("create file %s", filename);
	  streamID2 = streamOpenWrite(filename, cdoFiletype());
	  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", filename);

	  streamDefVlist(streamID2, vlistID2);

	  streamIDs[index] = streamID2;
	}


/*       printf("comp: %f %d\n",ndates*(index+1),(int)(ndates*(index+1))); */
      for ( ; nsets < (int)(ndates*(index+1)); nsets++ ) 
	{

	  nrecs = streamInqTimestep(streamID1, tsID);
/* 	  printf( "nsets: %d mod: %d\n",nsets, nsets % (int)(ndates) ); */

	  for ( recID = 0; recID < nrecs; recID++ )
	    {
	      
	      taxisCopyTimestep(taxisID2, taxisID1);
	      
/* 	      streamDefTimestep(streamID2, nsets % (int)ndates); */
	      streamDefTimestep(streamID2,tsIDs[index]);
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
	  tsIDs[index]++;	  
	}
      


      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;

/*       for ( i = 0; i < nskip; i++ ) */
/* 	{ */
/* 	  nrecs = streamInqTimestep(streamID1, tsID); */
/* 	  if ( nrecs == 0 ) break; */
/* 	  tsID++; */
/* 	} */
      for ( ; i2 < (int)(nskip*(index+1)); i2++ )
	{
	  nrecs = streamInqTimestep(streamID1, tsID);
	  if ( nrecs == 0 ) break;
	  tsID++;
	}

      nrecs = streamInqTimestep(streamID1, tsID);
      if ( nrecs == 0 ) break;

      index++;
    } // while(1)

 LABEL_END:


  streamClose(streamID1);

  for ( index = 0; index < MAX_STREAMS; index++ )
    {
      streamID2 = streamIDs[index];
      if ( streamID2 >= 0 )  streamClose(streamID2);
    }
 
  if ( ! lcopy )
    if ( array ) free(array);

  cdoFinish();

  return (0);

  return (0);
}
