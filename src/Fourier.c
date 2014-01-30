/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2014 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#if defined(_OPENMP)
#  include <omp.h>
#endif

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "statistic.h"


#define  NALLOC_INC  1024


void *Fourier(void *argument)
{
  int ompthID;
  int bit, sign;
  int gridsize;
  int nrecs;
  int gridID, varID, levelID, recID;
  int tsID;
  int i;
  int nts;
  int nalloc = 0;
  int streamID1, streamID2;
  int vlistID1, vlistID2, taxisID1, taxisID2;
  int nmiss;
  int nvars, nlevel;
  int *vdate = NULL, *vtime = NULL;
  double missval;
  field_t ***vars = NULL;
  typedef struct
  {
    double *real;
    double *imag;
    double *work_r;
    double *work_i;
  } memory_t;
  memory_t *mem = NULL;


  cdoInitialize(argument);

  operatorInputArg("the sign of the exponent (-1 for normal or 1 for reverse transformation)!");
  sign = atoi(operatorArgv()[0]);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  nvars = vlistNvars(vlistID1);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      if ( tsID >= nalloc )
	{
	  nalloc += NALLOC_INC;
	  vdate = (int *) realloc(vdate, nalloc*sizeof(int));
	  vtime = (int *) realloc(vtime, nalloc*sizeof(int));
	  vars  = (field_t ***) realloc(vars, nalloc*sizeof(field_t **));
	}

      vdate[tsID] = taxisInqVdate(taxisID1);
      vtime[tsID] = taxisInqVtime(taxisID1);

      vars[tsID] = field_malloc(vlistID1, FIELD_NONE);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  gridID   = vlistInqVarGrid(vlistID1, varID);
	  gridsize = gridInqSize(gridID);
	  vars[tsID][varID][levelID].ptr = malloc(2*gridsize*sizeof(double));
	  streamReadRecord(streamID1, vars[tsID][varID][levelID].ptr, &nmiss);
	  vars[tsID][varID][levelID].nmiss = nmiss;
	}

      tsID++;
    }

  nts = tsID;

  for ( bit = nts; !(bit & 1); bit >>= 1 );

  mem = malloc(ompNumThreads*sizeof(memory_t));
  for ( i = 0; i < ompNumThreads; i++ )
    {
      mem[i].real = malloc(nts*sizeof(double));
      mem[i].imag = malloc(nts*sizeof(double));
      if ( bit != 1 )
	{
	  mem[i].work_r = malloc(nts*sizeof(double));
	  mem[i].work_i = malloc(nts*sizeof(double));
	}
    }

  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID   = vlistInqVarGrid(vlistID1, varID);
      missval  = vlistInqVarMissval(vlistID1, varID);
      gridsize = gridInqSize(gridID);
      nlevel   = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for ( levelID = 0; levelID < nlevel; levelID++ )
	{
#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(i, ompthID, tsID)
#endif
	  for ( i = 0; i < gridsize; i++ )
	    {
	      int lmiss = 0;
#if defined(_OPENMP)
              ompthID = omp_get_thread_num();
#else
              ompthID = 0;
#endif
	      for ( tsID = 0; tsID < nts; tsID++ )
		{
		  mem[ompthID].real[tsID] = vars[tsID][varID][levelID].ptr[2*i];
		  mem[ompthID].imag[tsID] = vars[tsID][varID][levelID].ptr[2*i+1];
		  if ( DBL_IS_EQUAL(mem[ompthID].real[tsID], missval) ||
		       DBL_IS_EQUAL(mem[ompthID].imag[tsID], missval) ) lmiss = 1;
		}

	      if ( lmiss == 0 )
		{
		  if ( bit == 1 )	/* nts is a power of 2 */
		    fft(mem[ompthID].real, mem[ompthID].imag, nts, sign);
		  else
		    ft_r(mem[ompthID].real, mem[ompthID].imag, nts, sign, mem[ompthID].work_r, mem[ompthID].work_i);

		  for ( tsID = 0; tsID < nts; tsID++ )
		    {
		      vars[tsID][varID][levelID].ptr[2*i]   = mem[ompthID].real[tsID];
		      vars[tsID][varID][levelID].ptr[2*i+1] = mem[ompthID].imag[tsID];
		    }
		}
	      else
		{
		  for ( tsID = 0; tsID < nts; tsID++ )
		    {
		      vars[tsID][varID][levelID].ptr[2*i]   = missval;
		      vars[tsID][varID][levelID].ptr[2*i+1] = missval;
		    }
		}
	    }
	}
    }

  for ( i = 0; i < ompNumThreads; i++ )
    {
      free(mem[i].real);
      free(mem[i].imag);
      if ( bit != 1 )
	{
	  free(mem[i].work_r);
	  free(mem[i].work_i);
	}
    }
  free(mem);

  for ( tsID = 0; tsID < nts; tsID++ )
    {
      taxisDefVdate(taxisID2, vdate[tsID]);
      taxisDefVtime(taxisID2, vtime[tsID]);
      streamDefTimestep(streamID2, tsID);

      for ( varID = 0; varID < nvars; varID++ )
	{
	  nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
	  for ( levelID = 0; levelID < nlevel; levelID++ )
	    {
	      if ( vars[tsID][varID][levelID].ptr )
		{
		  nmiss = vars[tsID][varID][levelID].nmiss;
		  streamDefRecord(streamID2, varID, levelID);
		  streamWriteRecord(streamID2, vars[tsID][varID][levelID].ptr, nmiss);
		  free(vars[tsID][varID][levelID].ptr);
		  vars[tsID][varID][levelID].ptr = NULL;
		}
	    }
	}

      field_free(vars[tsID], vlistID1);
    }

  if ( vars  ) free(vars);
  if ( vdate ) free(vdate);
  if ( vtime ) free(vtime);

  streamClose(streamID2);
  streamClose(streamID1);

  cdoFinish();

  return (NULL);
}
