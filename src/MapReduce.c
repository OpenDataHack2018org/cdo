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
   Pack    unreduce
   */

#if defined(_OPENMP)
#  include <omp.h>
#endif

#include <limits.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"

int matrix2vector(int i, int j, int Ni, int Nj)
{
  return i*Nj + j;
}
int *vector2matrix(int k, int Ni, int Nj)
{
  int *retval;
  retval = (int *) malloc(2*sizeof(int));

  retval[0] = k/Nj; retval[1] = k%Nj;
  return retval;
}

/* count the number of locations, for which the mask is TRUE, i.e. has a
 * non-zero value */
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
 * collect the indices of relevant location out of the given mask
 */
void collectLocations(double *maskField, int gridSize, double falseVal, int *maskIndexList)
{
  int k = 0;
  for (int i = 0; i < gridSize; i++)
  {
    if (!DBL_IS_EQUAL(maskField[i],falseVal))
    {
      printf("found at:%d -",i);
      k += 1;
    }
  }
      printf(" k =%d ",k);
}

/*
 * the operators argument has to be a single horizontal field,
 * non-zero values are used to mark the relevant locations
 */
void *MapReduce(void *argument)
{
  int gridsize;
  int nrecs;
  int tsID;
  int gridID, varID, levelID, recID;
  int i;
  int nts;
  int nalloc = 0;
  int nmiss;
  int nlevel;
  int datatype = DATATYPE_INT16;
  dtlist_type *dtlist = dtlist_new();
  double missval1, missval2;
  field_t ***vars = NULL;

  cdoInitialize(argument);

  int streamID1 = streamOpenRead(cdoStreamName(0));
  int vlistID1  = streamInqVlist(streamID1);
  int nvars     = vlistNvars(vlistID1);
  int vlistID2  = vlistDuplicate(vlistID1);

  /* check input grid type and size */
  int inputGridID = cdoDefineGrid(operatorArgv()[0]);
  int inputGridSize = gridInqSize(inputGridID);
  if ( cdoVerbose ) cdoPrint("input gridSize:%d", inputGridSize);

  /* search the number of relevant locations */
  tsID = 0; nrecs = 0;
  double *inputMaskField = (double*) Malloc(inputGridSize*sizeof(double));
  streamReadRecord(streamID1, inputMaskField, &nmiss);
  minmaxval(inputGridSize, inputMaskField, NULL,&missval1, &missval2);
  cdoPrint("min: %g | max: %g ",missval1, missval2);
  /* count points {{{*/
  int maskSize = countMask(inputMaskField, inputGridSize, 0.0);
  if (cdoDebug) cdoPrint("maskSize = %d",maskSize); /* }}} */

  /* collect the original coordinates */
   int *maskIndexList = (int *) Malloc(maskSize*sizeof(int));
   for (int m = 0; m < maskSize; m++) maskIndexList[m] = -1;

   /* create an index list of relevant points {{{ */
   int k = 0;
   for (int i = 0; i < inputGridSize; i++)
   {
     if (!DBL_IS_EQUAL(inputMaskField[i],0.0))
     {
       if (cdoDebug) printf("found at:%d -",i);
       maskIndexList[k] = i;
       k += 1;
     }
   }

   if (cdoDebug) for (int l = 0; l < maskSize; l++) cdoPrint("maskIndexList[%d] = %d",l,maskIndexList[l]);
   /* }}} */
 
  /* create unstructured output grid */
  int outputGridID = gridToUnstructuredSelecton(inputGridID,1, maskSize, maskIndexList);

  int inputGridType = gridInqType(inputGridID);

  /* copy time axis */
  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);


  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
  {
    tsID++;
  }

  vlistChangeGridIndex(vlistID2, 0, outputGridID);
  int tsteptype = TSTEP_CONSTANT;
  int zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);

  varID = vlistDefVar(vlistID2, outputGridID, zaxisID, tsteptype);
  vlistDefVarName(vlistID2    , varID , "mask");
  vlistDefVarStdname(vlistID2 , varID , "grid_mask");
  vlistDefVarLongname(vlistID2, varID , "mask");
  vlistDefVarUnits(vlistID2   , varID , "");

  double *values = (double *)Malloc(maskSize*sizeof(double));
  for (int i = 0; i < maskSize;  i++) values[i] = 1.0;

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  streamDefVlist(streamID2, vlistID2);

  streamDefRecord(streamID2, 0, 0);
  streamWriteRecord(streamID2, values, 0);

  streamClose(streamID2);
  streamClose(streamID1);

  Free(inputMaskField);
  Free(maskIndexList);

  cdoFinish();

  return 0;
}
