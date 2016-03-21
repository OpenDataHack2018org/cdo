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

   Pack      reduce
*/

#if defined(_OPENMP)
#  include <omp.h>
#endif

#include <limits.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"

/* read only the first data variable from input filename into a given double
 * pointer */
/* DON'T MOVE IT! is necessary to have the pstream.h file included AFTER this
 * function definition */
void read_first_record(char *filename, int gridSize, double *field)
{
  int nmiss,varID,levelID;
  int streamID = streamOpenRead(filename);
  int nrecs = streamInqTimestep(streamID,0);
  streamInqRecord(streamID,&varID,&levelID);
  streamReadRecord(streamID, field, &nmiss);
  streamClose(streamID);
  nrecs = nrecs;
}

#include "pstream.h"

/*
 * count the number of locations, for which the mask is TRUE
 * */
int countMask(double *maskField, int gridSize, double falseVal)
{
  int counter;

  counter = 0;

  for (int i = 0; i < gridSize; i++)
  {
    if (!DBL_IS_EQUAL(maskField[i],falseVal)) counter += 1;
  }
  return counter;
}

/* 
 * the operators argument has to be a single horizontal field,
 * non-zero values are used to mark the relevant locations
 */
void *MapReduce(void *argument)
{
  int nrecs;
  int tsID;
  int varID, levelID, recID;
  int nmiss;
  /*double missval1, missval2; */

  cdoInitialize(argument);

  /* check input grid type and size - this will be used for selecting relevant
   * variables from the input file*/
  int inputGridID   = cdoDefineGrid(operatorArgv()[0]);
  int inputGridSize = gridInqSize(inputGridID);
  int inputGridType = gridInqType(inputGridID);
  if ( cdoVerbose ) cdoPrint("MapReduce: input gridSize:%d", inputGridSize);

  /* creata an index list of the relevant locations  {{{ */
  tsID = 0;
  double *inputMaskField = (double*) Malloc(inputGridSize*sizeof(double));
  read_first_record(operatorArgv()[0],inputGridSize, inputMaskField);

  /* non-zero values mark the relevant points */
  int maskSize = countMask(inputMaskField, inputGridSize, 0.0);
  cdoPrint("MapReduce: maskSize = %d",maskSize);

  int *maskIndexList = (int *) Malloc(maskSize*sizeof(int));
  for (int m = 0; m < maskSize; m++) maskIndexList[m] = -1;

  int k = 0;
  for (int i = 0; i < inputGridSize; i++)
  {
    if (!DBL_IS_EQUAL(inputMaskField[i],0.0))
    {
      maskIndexList[k] = i;
      k += 1;
    }
  }
  /* }}} */

  /* create unstructured output grid including bounds*/
  int outputGridID = gridToUnstructuredSelecton(inputGridID, maskSize, maskIndexList);

  /* create output vlist: Only variabes which have the same gridtype and
   * gridsize as the input mask should be proessed. Everything else is ignoreds
   * {{{ */
  int streamID1 = streamOpenRead(cdoStreamName(0));
  int vlistID1  = streamInqVlist(streamID1);
  int nvars     = vlistNvars(vlistID1);
  int *vars     = (int*) Malloc(nvars*sizeof(int));

  /* use vlist flags for marking the corresponding variables */
  vlistClearFlag(vlistID1);
  for ( varID = 0; varID < nvars; varID++ )
  {
    vars[varID] = FALSE;

    int gridID = vlistInqVarGrid(vlistID1, varID);
    if (inputGridType == gridInqType(gridID) && inputGridSize == gridInqSize(gridID))
    {
      vars[varID] = TRUE;
      int zaxisID  = vlistInqVarZaxis(vlistID1, varID);
      int nlevs    = zaxisInqSize(zaxisID);
      for ( int levID = 0; levID < nlevs; levID++ )
      {
        vlistDefFlag(vlistID1, varID, levID, TRUE);
      }
    }
  }
  int vlistID2 = vlistCreate();
  vlistCopyFlag(vlistID2, vlistID1);
  /* }}} */

  int taxisID1  = vlistInqTaxis(vlistID1);
  int taxisID2  = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  /* use the new selection grid for all output variables */
  int ngrids = vlistNgrids(vlistID2);
  for ( int index = 0; index < ngrids; index++ ) vlistChangeGridIndex(vlistID2, index, outputGridID);

  /* loop over input fields and mask the data values {{{ */
  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  streamDefVlist(streamID2, vlistID2);

  double *arrayIn  = (double *)Malloc(inputGridSize*sizeof(double));
  double *arrayOut = (double *)Malloc(maskSize*sizeof(double));

  tsID = 0; 
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
  {
    taxisCopyTimestep(taxisID2, taxisID1);
    streamDefTimestep(streamID2, tsID);

    for ( recID = 0; recID < nrecs; recID++ )
    {
      streamInqRecord(streamID1, &varID, &levelID);
      if (TRUE == vars[varID])
      {
        int varID2   = vlistFindVar(vlistID2, varID);
        int levelID2 = vlistFindLevel(vlistID2, varID, levelID);

        streamReadRecord(streamID1, arrayIn, &nmiss);

        for (int i = 0; i < maskSize;  i++)
          arrayOut[i] = arrayIn[maskIndexList[i]];


        streamDefRecord(streamID2, varID2, levelID2);
        streamWriteRecord(streamID2, arrayOut, 0);

      }
    }
    tsID++;
  }
  /* }}} */


  streamClose(streamID2);
  streamClose(streamID1);

  Free(inputMaskField);
  Free(maskIndexList);

  cdoFinish();

  return 0;
}
