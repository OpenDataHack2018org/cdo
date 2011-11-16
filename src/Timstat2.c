/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2011 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

        Timstat2        timcor      correlates two data files on the same grid
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


static
int correlation(long gridsize, double missval1, double missval2, int *nofvals, 
		double *work0, double *work1, double *work2, double *work3, double *work4)
{
  long i;
  int nvals, nmiss = 0;
  double temp0, temp1, temp2, temp3, temp4, temp5, temp6;
  double cor;

  for ( i = 0; i < gridsize; i++ )
    {	  
      nvals = nofvals[i];

      if ( nvals > 0 )
	{
	  temp0 = MUL(work0[i], work1[i]);
	  temp1 = SUB(work4[i], DIV(temp0, nvals));
	  temp2 = MUL(work0[i], work0[i]);
	  temp3 = MUL(work1[i], work1[i]);
	  temp4 = SUB(work2[i], DIV(temp2, nvals));
	  temp5 = SUB(work3[i], DIV(temp3, nvals));
	  temp6 = MUL(temp4, temp5);

	  cor = DIV(temp1, SQRT(temp6));
	  /*
	    if      ( cor < -1)  cor = -1;
	    else if ( cor >  1)  cor =  1;
	  */
	  if ( DBL_IS_EQUAL(cor, missval1) ) nmiss++;
	}
      else
	{
	  nmiss++;
	  cor = missval1;
	}

      work0[i] = cor;
    }

  return nmiss;
}


void *Timstat2(void *argument)
{
  int operatorID;
  int operfunc;
  int nwork = 0;
  int streamID1, streamID2, streamID3;
  int vdate = 0, vtime = 0;
  int nrecs, nrecs2, nrecs3, nvars, nlevs;
  long i, gridsize;
  int tsID;
  int varID, recID, levelID, gridID;
  int nmiss;
  int *recVarID, *recLevelID;
  int vlistID1, vlistID2, vlistID3;
  int taxisID1, taxisID2, taxisID3;
  double missval1, missval2;
  double ****work = NULL;
  int ***nofvals = NULL;
  double *array1 = NULL, *array2 = NULL;

  cdoInitialize(argument);

  cdoOperatorAdd("timcor",   func_cor,   0, NULL);
  cdoOperatorAdd("timcovar", func_covar, 0, NULL);

  operatorID = cdoOperatorID();
  operfunc   = cdoOperatorF1(operatorID);

  if      ( operfunc == func_cor   ) nwork = 5;
  else if ( operfunc == func_covar ) nwork = 3;

  streamID1 = streamOpenRead(cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);
  vlistID3 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_ALL);
 
  nvars  = vlistNvars(vlistID1);
  nrecs  = vlistNrecs(vlistID1);
  nrecs3 = nrecs;
  recVarID   = (int *) malloc(nrecs*sizeof(int));
  recLevelID = (int *) malloc(nrecs*sizeof(int));
  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = vlistInqTaxis(vlistID2);
  taxisID3 = taxisDuplicate(taxisID1);
 
  vlistDefTaxis(vlistID3, taxisID3);
  streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());

  streamDefVlist(streamID3, vlistID3);
 
  gridsize = vlistGridsizeMax(vlistID1);

  array1  = (double *) malloc(gridsize*sizeof(double));
  array2  = (double *) malloc(gridsize*sizeof(double));
  				 
  work    = (double ****) malloc(nvars*sizeof(double ***));
  nofvals = (int ***)     malloc(nvars*sizeof(int **));

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, 0);  
      gridsize = gridInqSize(gridID);
      nlevs    = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));

      work[varID]    = (double ***)  malloc(nlevs*sizeof(double **));
      nofvals[varID] = (int **)      malloc(nlevs*sizeof(int *));  

      for ( levelID = 0; levelID < nlevs; levelID++ )
	{
	  nofvals[varID][levelID] = (int *) malloc(gridsize*sizeof(int));
	  memset(nofvals[varID][levelID], 0, gridsize*sizeof(int));
      
	  work[varID][levelID] = (double **) malloc(nwork*sizeof(double *));
	  for ( i = 0; i < nwork; i++ )
	    {
	      work[varID][levelID][i] = (double *) malloc(gridsize*sizeof(double));
	      memset(work[varID][levelID][i], 0, gridsize*sizeof(double));
	    }
	}
    }
 
  tsID=0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      vdate = taxisInqVdate(taxisID1);
      vtime = taxisInqVtime(taxisID1);

      nrecs2 = streamInqTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamInqRecord(streamID2, &varID, &levelID);

	  if ( tsID == 0 )
	    {
	      recVarID[recID]   = varID;
	      recLevelID[recID] = levelID;	     	     
	    }	 

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

	  missval1 = vlistInqVarMissval(vlistID1, varID);
	  missval2 = vlistInqVarMissval(vlistID2, varID);

	  streamReadRecord(streamID1, array1, &nmiss);
	  streamReadRecord(streamID2, array2, &nmiss);
	      	     
	  for ( i = 0; i < gridsize; i++)
	    {
	      if ( ( ! DBL_IS_EQUAL(array1[i], missval1) ) && 
		   ( ! DBL_IS_EQUAL(array2[i], missval2) ) )
		{
		  work[varID][levelID][0][i] += array1[i];
		  work[varID][levelID][1][i] += array2[i];
		  work[varID][levelID][2][i] += ( array1[i]*array1[i] );
		  work[varID][levelID][3][i] += ( array2[i]*array2[i] );
		  work[varID][levelID][4][i] += ( array1[i]*array2[i] );
		  nofvals[varID][levelID][i]++;
		}
	    }	 
	}

      tsID++;
    }

  tsID = 0;
  taxisDefVdate(taxisID3, vdate);
  taxisDefVtime(taxisID3, vtime);
  streamDefTimestep(streamID3, tsID);

  for ( recID = 0; recID < nrecs3; recID++ )
    {
      varID    = recVarID[recID];
      levelID  = recLevelID[recID];
   
      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

      missval1 = vlistInqVarMissval(vlistID1, varID);
      missval2 = vlistInqVarMissval(vlistID2, varID);

      nmiss = correlation(gridsize, missval1, missval2, nofvals[varID][levelID],
			  work[varID][levelID][0], work[varID][levelID][1],
			  work[varID][levelID][2], work[varID][levelID][3], 
			  work[varID][levelID][4]);

      streamDefRecord(streamID3, varID, levelID);
      streamWriteRecord(streamID3, work[varID][levelID][0], nmiss);
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for ( levelID = 0; levelID < nlevs; levelID++ )
	{
	  free(nofvals[varID][levelID]);
	  for ( i = 0; i < nwork; i++ )
	    free(work[varID][levelID][i]);
	  free(work[varID][levelID]);
	}
    
      free(nofvals[varID]);
      free(work[varID]);
    }
    
  free(nofvals);
  free(work);

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  if ( array1 )     free(array1);
  if ( array2 )     free(array2);
  if ( recVarID   ) free(recVarID);
  if ( recLevelID ) free(recLevelID);
    
  cdoFinish();   
 
  return (0);
}
