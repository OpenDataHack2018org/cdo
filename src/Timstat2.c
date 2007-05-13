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

        Timcor        timcor      correlates two data files on the same grid
*/


#include <string.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

void *Timstat2(void *argument)
{
  static char func[] = "func_cor";
  int operatorID;
  int operfunc;
  int intvdat;
  int indate1, indate2 = 0;
  int streamID1, streamID2, streamID3;
  int gridsize;
  int vdate = 0, vtime = 0;
  int vdate0 = 0, vtime0 = 0;
  int nrecs, nrecs2, nrecs3, nvars, nlevs ;
  int i;
  int tsID, otsID;
  int varID, recID, levelID, gridID;
  int nsets, nmiss1, nmiss2;
  int *recVarID, *recLevelID;
  int offset;
  int vlistID1, vlistID2, vlistID3;
  int taxisID1, taxisID2, taxisID3;
  int  cor0;
  double missval, missval1, missval2;
  FIELD **vars1 = NULL, **samp1 = NULL;
  FIELD **vars2 = NULL, ***vars3 = NULL;
  FIELD ***temp;
  int ***nofvals;
  FIELD field1, field2;


  cdoInitialize(argument);
  cdoOperatorAdd("timcor", 2, 1, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);
  vlistID3 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, func_sft);
 
  gridsize = vlistGridsizeMax(vlistID1);
  nvars = vlistNvars(vlistID1);
  nrecs = vlistNrecs(vlistID1);
  nrecs3 = nrecs;
  recVarID   = (int *) malloc(nrecs*sizeof(int));
  recLevelID = (int *) malloc(nrecs*sizeof(int));
  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = vlistInqTaxis(vlistID2);
  taxisID3 = taxisDuplicate(taxisID1);
 
  vlistDefTaxis(vlistID3, taxisID3);
  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());
  if ( streamID3 < 0 ) cdiError(streamID3, "Open failed on %s", cdoStreamName(2));

  streamDefVlist(streamID3, vlistID3);
 
  field1.ptr = (double *) malloc(gridsize*sizeof(double));
  field2.ptr = (double *) malloc(gridsize*sizeof(double));
				 

  vars1 = (FIELD **)  malloc(nvars*sizeof(FIELD *));
  vars2 = (FIELD **)  malloc(nvars*sizeof(FIELD *));
  vars3 = (FIELD ***) malloc(nvars*sizeof(FIELD **));
  nofvals = (int ***)  malloc(nvars*sizeof(int **));
  temp  = (FIELD ***) malloc(nvars*sizeof(FIELD **));
  samp1 = (FIELD **)  malloc(nvars*sizeof(FIELD *));
 
  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, 0);
      gridsize = gridInqSize(gridID);
      nlevs    = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      missval = missval1 = vlistInqVarMissval(vlistID1, varID);
      missval2 = vlistInqVarMissval(vlistID1, varID); 

      vars1[varID] = (FIELD *)   malloc(nlevs*sizeof(FIELD));
      samp1[varID] = (FIELD *)   malloc(nlevs*sizeof(FIELD));
      vars2[varID] = (FIELD *)   malloc(nlevs*sizeof(FIELD));
      vars3[varID] = (FIELD **)  malloc(nlevs*sizeof(FIELD*));
      temp[varID]  = (FIELD **)  malloc(nlevs*sizeof(FIELD*));
      nofvals[varID] = (int **)  malloc(nlevs*sizeof(int*));
      for ( i = 0; i < nlevs; i++ )
	{
	  /* XXX */
	  nofvals[varID][i] = (int *) malloc(gridsize*sizeof(int));
	  vars3[varID][i] = (FIELD *) malloc(7*sizeof(FIELD));
	  temp[varID][i] =  (FIELD *) malloc(7*sizeof(FIELD));
	}

      for ( levelID = 0; levelID < nlevs; levelID++ )
	{
	  vars1[varID][levelID].grid    = gridID;
	  vars1[varID][levelID].nmiss   = 0;
	  vars1[varID][levelID].missval = missval1;
	  vars1[varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
	  samp1[varID][levelID].grid    = gridID;
	  samp1[varID][levelID].nmiss   = 0;
	  samp1[varID][levelID].missval = missval1;
	  samp1[varID][levelID].ptr     = NULL;

	  vars2[varID][levelID].grid    = gridID;
	  vars2[varID][levelID].nmiss   = 0;
	  vars2[varID][levelID].missval = missval2;
	  vars2[varID][levelID].ptr     = (double *) malloc(gridsize*sizeof(double));
	  /* XXX */
	  memset(nofvals[varID][levelID], 0, gridsize*sizeof(int));

	  for ( i = 0; i < 7; i++ )
	    {
	      vars3[varID][levelID][i].grid    = gridID;
	      vars3[varID][levelID][i].nmiss   = 0;
	      vars3[varID][levelID][i].missval = missval;
	      vars3[varID][levelID][i].ptr     = (double *) malloc(gridsize*sizeof(double));
	      memset(vars3[varID][levelID][i].ptr, 0, gridsize*sizeof(double));

	      temp[varID][levelID][i].grid    = gridID;
	      temp[varID][levelID][i].nmiss   = 0;
	      temp[varID][levelID][i].missval = missval;
	      temp[varID][levelID][i].ptr     = (double *) malloc(gridsize*sizeof(double));
	      memset(temp[varID][levelID][i].ptr, 0, gridsize*sizeof(double));
	    }
	}
    }

  nsets = 0;
  tsID=0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      nrecs2 = streamInqTimestep(streamID2, tsID);

      for ( recID = 0; recID<nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamInqRecord(streamID2, &varID, &levelID);

	  if ( tsID == 0 )
	    {
	      recVarID[recID] = varID;
	      recLevelID[recID] = levelID;
	    }
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  if ( nsets == 0 )
	    {
	      streamReadRecord(streamID1, field1.ptr, &nmiss1);
	      streamReadRecord(streamID2, field2.ptr, &nmiss2);
	      vars1[varID][levelID].nmiss = nmiss1;
	      vars2[varID][levelID].nmiss = nmiss2;

	      if ( nmiss1 > 0 || nmiss2 > 0  || samp1[varID][levelID].ptr )
		{
		  if ( samp1[varID][levelID].ptr == NULL )
		    samp1[varID][levelID].ptr = (double *) malloc(gridsize*sizeof(double));
		  for ( i = 0; i < gridsize; i++ )
		    if ( DBL_IS_EQUAL(field1.ptr[i], vars1[varID][levelID].missval) )
		      samp1[varID][levelID].ptr[i] = -1;
		    else if (DBL_IS_EQUAL(field2.ptr[i], vars2[varID][levelID].missval) )
		      samp1[varID][levelID].ptr[i] = -1;
		    else
		      samp1[varID][levelID].ptr[i] = 0;
		}
	      for ( i = 0; i < gridsize; i++)
		{
                  if ( ( ! DBL_IS_EQUAL(field1.ptr[i], vars1[varID][levelID].missval ) ) && 
	       	     ( ! DBL_IS_EQUAL(field2.ptr[i], vars2[varID][levelID].missval ) ) )
		    {
		      vars3[varID][levelID][0].ptr[i] += field1.ptr[i];
		      vars3[varID][levelID][1].ptr[i] += field2.ptr[i];
		      vars3[varID][levelID][2].ptr[i] += ( field1.ptr[i]*field1.ptr[i] );
		      vars3[varID][levelID][3].ptr[i] += ( field2.ptr[i]*field2.ptr[i] );
		      vars3[varID][levelID][4].ptr[i] += ( field1.ptr[i]*field2.ptr[i] );
		      nofvals[varID][levelID][i]++;
		    }
		}
	    }
	  else
	    {
	      streamReadRecord(streamID1, field1.ptr, &field1.nmiss);
	      field1.grid    = vars1[varID][levelID].grid;
	      field1.missval = vars1[varID][levelID].missval;
	    
	      streamReadRecord(streamID2, field2.ptr, &field2.nmiss);
	      field2.grid    = vars2[varID][levelID].grid;
	      field2.missval = vars2[varID][levelID].missval;
	      gridsize = gridInqSize(gridID);

	      if ( field1.nmiss > 0 || field2.nmiss > 0 || samp1[varID][levelID].ptr )
		{
		  if ( samp1[varID][levelID].ptr == NULL )
		    {
		      samp1[varID][levelID].ptr = (double *) malloc(gridsize*sizeof(double));
		      for ( i = 0; i < gridsize; i++ )
			samp1[varID][levelID].ptr[i] = (nsets-1);
		    }
		  for ( i = 0; i < gridsize; i ++ )
			if ( ( ! DBL_IS_EQUAL(field1.ptr[i], vars1[varID][levelID].missval ) ) &&
			     ( ! DBL_IS_EQUAL(field2.ptr[i], vars2[varID][levelID].missval ) ) )
			  samp1[varID][levelID].ptr[i]++;
		}
	      for ( i = 0; i < gridsize; i++)
		{
                  if ( ( ! DBL_IS_EQUAL(field1.ptr[i], vars1[varID][levelID].missval ) ) && 
	       	     ( ! DBL_IS_EQUAL(field2.ptr[i], vars2[varID][levelID].missval ) ) )
		    {
		      vars3[varID][levelID][0].ptr[i] += field1.ptr[i];
		      vars3[varID][levelID][1].ptr[i] += field2.ptr[i];
		      vars3[varID][levelID][2].ptr[i] += ( field1.ptr[i]*field1.ptr[i] );
		      vars3[varID][levelID][3].ptr[i] += ( field2.ptr[i]*field2.ptr[i] );
		      vars3[varID][levelID][4].ptr[i] += ( field1.ptr[i]*field2.ptr[i] );
		      nofvals[varID][levelID][i]++;
		    }
		}
	    } 
	}
      /* XXX */
      tsID++;
      nsets++;
    }

  printf("nsets %d %d\n", tsID, nsets );
  tsID = 0;
  taxisDefVdate(taxisID3, vdate);
  taxisDefVtime(taxisID3, vtime);
  streamDefTimestep(streamID3, tsID);

  for ( recID = 0;  recID<nrecs3; recID++ )
    {
      varID =	recVarID[recID];
      levelID = recLevelID[recID];
      missval1 = vlistInqVarMissval(vlistID1, varID);
      missval2 = vlistInqVarMissval(vlistID1, varID);
      for ( i = 0; i < gridsize; i++ )
	{
	  if ( samp1[varID][levelID].ptr )
	    nofvals[varID][levelID][i] = samp1[varID][levelID].ptr[i];


	  if (nofvals[varID][levelID][i] > 0 )
	    {
	      /* XXX */
	      if ( i == 0 ) printf("nofvals %d\n", nofvals[varID][levelID][i]);
	      temp[varID][levelID][0].ptr[i] = MUL (vars3[varID][levelID][0].ptr[i], vars3[varID][levelID][1].ptr[i]);
	      temp[varID][levelID][1].ptr[i] = SUB (vars3[varID][levelID][4].ptr[i], DIV(temp[varID][levelID][0].ptr[i], nofvals[varID][levelID][i]));
	      temp[varID][levelID][2].ptr[i] = MUL (vars3[varID][levelID][0].ptr[i], vars3[varID][levelID][0].ptr[i]);
	      temp[varID][levelID][3].ptr[i] = MUL (vars3[varID][levelID][1].ptr[i], vars3[varID][levelID][1].ptr[i]);
	      temp[varID][levelID][4].ptr[i] = SUB (vars3[varID][levelID][2].ptr[i], DIV(temp[varID][levelID][2].ptr[i], nofvals[varID][levelID][i]));
	      temp[varID][levelID][5].ptr[i] = SUB (vars3[varID][levelID][3].ptr[i], DIV(temp[varID][levelID][3].ptr[i], nofvals[varID][levelID][i]));
	      temp[varID][levelID][6].ptr[i] = MUL (temp[varID][levelID][4].ptr[i], temp[varID][levelID][5].ptr[i]);
	      vars3[varID][levelID][0].ptr[i] = DIV (temp[varID][levelID][1].ptr[i], sqrt(temp[varID][levelID][6].ptr[i]));
	      if ( vars3[varID][levelID][0].ptr[i] < -1)
		vars3[varID][levelID][0].ptr[i] = -1;
	      else if ( vars3[varID][levelID][0].ptr[i] > 1)
	      vars3[varID][levelID][0].ptr[i] = 1;
	    }
	  if ( samp1[varID][levelID].ptr && samp1[varID][levelID].ptr[i] <= 20 )
	    {
	      vars3[varID][levelID][0].nmiss++;
	      vars3[varID][levelID][0].ptr[i] = missval;
	    }
	}
      /* XXX */ 									       
      /* if ( otsID == 1 || vlistInqVarTime(vlistID1, varID == TIME_VARIABLE ) ) */
      if ( /* otsID == 1 || */ vlistInqVarTime(vlistID1, varID == TIME_VARIABLE ) )
	{
	  streamDefRecord(streamID3, varID, levelID);
	  streamWriteRecord(streamID3, vars3[varID][levelID][0].ptr, vars3[varID][levelID][0].nmiss);
	}
    }
  for ( varID = 0; varID < nvars; varID++ )
    {
      nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for ( levelID = 0; levelID < nlevs; levelID++ )
	{
	  free(vars1[varID][levelID].ptr);
	  free(vars2[varID][levelID].ptr);
	  for ( i = 0; i < 7; i++ )
	    free(vars3[varID][levelID][i].ptr);
	  if ( samp1[varID][levelID].ptr) free(samp1[varID][levelID].ptr);
	}
    
      free(vars1[varID]);
      free(vars2[varID]);
      free(vars3[varID]);
      free(temp[varID]);
      free(samp1[varID]);
    }


  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);
    
  free(vars1);
  free(vars2);
  free(vars3);
  free(temp);
  free(samp1);
  if ( field1.ptr ) free(field1.ptr);
  if ( field2.ptr ) free(field2.ptr);
  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);
    
  cdoFinish();
    
  return (0);
}
