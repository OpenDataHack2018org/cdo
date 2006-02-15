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

#include <string.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

/*
@BeginDoc

@BeginModule

@Name      = Cond2
@Title     = 
@Section   = Conditions
@Class     = Comparison
@Arguments = ifile1 ifile2 ifile3 ofile
@Operators = ifthenelse

@BeginDesciption
@EndDesciption

@EndModule


@BeginOperator_ifthenelse

@Title     = If then else

@BeginDesciption
@IfMan
          / i_2(t,x) if i_1(t,x) NE 0  AND  i_1(t,x) NE miss
o(t,x) = <  i_3(t,x) if i_1(t,x) EQ 0  AND  i_1(t,x) NE miss
          \ miss     if i_1(t,x) EQ miss
@EndifMan
@IfDoc
@BeginMath
o(t,x) = \left\{
\begin{array}{cll}
   i_2(t,x) & \mbox{if} \;\;  i_1(t,x) \neq 0 & \wedge \;\; i_1(t,x) \neq \mbox{miss} \\
   i_3(t,x) & \mbox{if} \;\;  i_1(t,x)    = 0 & \wedge \;\; i_1(t,x) \neq \mbox{miss} \\
\mbox{miss} & \mbox{if} \;\;  i_1(t,x)    = \mbox{miss} \\
\end{array}   \right.
@EndMath
@EndifDoc
@EndDesciption

@EndOperator

@EndDoc
*/


void *Cond2(void *argument)
{
  static char func[] = "Cond2";
  int IFTHENELSE;
  int operatorID;
  int streamID1, streamID2, streamID3, streamID4;
  int gridsize;
  int nrecs, nrecs2, nvars = 0, nlev, recID;
  int tsID;
  int varID, levelID;
  int offset;
  int vlistID1, vlistID2, vlistID3, vlistID4;
  int nmiss1, nmiss2, nmiss3, nmiss4;
  int i;
  int lptype = 0;
  double missval1 = -9.E33;
  double *array1, *array2, *array3, *array4;
  int **varnmiss1 = NULL;
  double **vardata1 = NULL;
  int taxisID2, taxisID3, taxisID4;

  cdoInitialize(argument);

  IFTHENELSE = cdoOperatorAdd("ifthenelse",    0, 0, NULL);

  operatorID = cdoOperatorID();

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));
  streamID3 = streamOpenRead(cdoStreamName(2));
  if ( streamID3 < 0 ) cdiError(streamID3, "Open failed on %s", cdoStreamName(2));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);
  vlistID3 = streamInqVlist(streamID3);
  vlistID4 = vlistDuplicate(vlistID2);

  taxisID2 = vlistInqTaxis(vlistID2);
  taxisID3 = vlistInqTaxis(vlistID3);
  taxisID4 = taxisDuplicate(taxisID2);
  vlistDefTaxis(vlistID4, taxisID4);

  if ( vlistNrecs(vlistID2) == 1 && vlistNrecs(vlistID1) != 1 )
    lptype = 2;
  else
    vlistCompare(vlistID1, vlistID2, func_sft);

  nospec(vlistID1);

  streamID4 = streamOpenWrite(cdoStreamName(3), cdoFiletype());
  if ( streamID4 < 0 ) cdiError(streamID4, "Open failed on %s", cdoStreamName(3));

  streamDefVlist(streamID4, vlistID4);

  gridsize = vlistGridsizeMax(vlistID1);

  array1 = (double *) malloc(gridsize*sizeof(double));
  array2 = (double *) malloc(gridsize*sizeof(double));
  array3 = (double *) malloc(gridsize*sizeof(double));
  array4 = (double *) malloc(gridsize*sizeof(double));

  if ( cdoVerbose )
    cdoPrint("Number of timesteps: file1 %d, file2 %d, file3 %d",
	     vlistNtsteps(vlistID1), vlistNtsteps(vlistID2), vlistNtsteps(vlistID3));

  if ( vlistNtsteps(vlistID1) == 1 &&
       vlistNtsteps(vlistID2) != 1 && vlistNtsteps(vlistID3) != 1 && lptype == 0 )
    {
      lptype = 1;
      nvars  = vlistNvars(vlistID1);
      vardata1  = (double **) malloc(nvars*sizeof(double *));
      varnmiss1 = (int **) malloc(nvars*sizeof(int *));
      for ( varID = 0; varID < nvars; varID++ )
	{
	  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	  nlev     = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  vardata1[varID]  = (double *) malloc(nlev*gridsize*sizeof(double));
	  varnmiss1[varID] = (int *) malloc(nlev*sizeof(int));
	}

      cdoPrint("Filling up stream >%s< by copying the first timestep.", cdoStreamName(0));
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID2, tsID)) )
    {
      nrecs = streamInqTimestep(streamID3, tsID);
      if ( nrecs == 0 ) cdoAbort("Input streams have different number of timesteps!");

      if ( tsID == 0 || lptype == 0 )
	{
	  nrecs2 = streamInqTimestep(streamID1, tsID);
	  if ( nrecs2 == 0 )
	    cdoAbort("Input streams have different number of timesteps!");
	}

      taxisCopyTimestep(taxisID4, taxisID2);

      streamDefTimestep(streamID4, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID2, &varID, &levelID);
	  streamReadRecord(streamID2, array2, &nmiss2);

	  streamInqRecord(streamID3, &varID, &levelID);
	  streamReadRecord(streamID3, array3, &nmiss3);

	  if ( tsID == 0 || lptype == 0 )
	    {
	      if ( recID == 0 || lptype != 2 )
		{
		  streamInqRecord(streamID1, &varID, &levelID);
		  streamReadRecord(streamID1, array1, &nmiss1);
		}

	      if ( lptype == 1 )
		{
		  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
		  offset   = gridsize*levelID;
		  memcpy(vardata1[varID]+offset, array1, gridsize*sizeof(double));
		  varnmiss1[varID][levelID] = nmiss1;
		}
	    }
	  else if ( lptype == 1 )
	    {
	      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
	      offset   = gridsize*levelID;
	      memcpy(array1, vardata1[varID]+offset, gridsize*sizeof(double));
	      nmiss1 = varnmiss1[varID][levelID];
	    }

	  gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID));

	  if ( recID == 0 || lptype != 2 )
	    {
	      /* gridsize = gridInqSize(vlistInqVarGrid(vlistID2, varID)); */
	      missval1  = vlistInqVarMissval(vlistID1, varID);
	    }

	  if ( operatorID == IFTHENELSE )
	    {
	      for ( i = 0; i < gridsize; i++ )
		array4[i] = DBL_IS_EQUAL(array1[i], missval1) ?
		  missval1 : !DBL_IS_EQUAL(array1[i], 0) ? array2[i] : array3[i];
	    }
	  else
	    {
	      cdoAbort("Operator not implemented!");
	    }

	  nmiss4 = 0;
	  for ( i = 0; i < gridsize; i++ )
	    if ( DBL_IS_EQUAL(array4[i], missval1) ) nmiss4++;

	  streamDefRecord(streamID4, varID, levelID);
	  streamWriteRecord(streamID4, array4, nmiss4);
	}

      tsID++;
    }

  streamClose(streamID4);
  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  if ( vardata1 )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  free(vardata1[varID]);
	  free(varnmiss1[varID]);
	}

      free(vardata1);
      free(varnmiss1);
    }

  if ( array4 ) free(array4);
  if ( array3 ) free(array3);
  if ( array2 ) free(array2);
  if ( array1 ) free(array1);

  cdoFinish();

  return (0);
}
