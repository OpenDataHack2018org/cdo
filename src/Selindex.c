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

      Selindex     selindex    Select grid indices
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "listarray.h"
#include "pstream.h"


int gengridcell(int gridID1, int gridsize2, int *cellidx);

static
int genindexgrid(int gridID1, int gridsize2, int *cellidx)
{
  int gridID0 = gridID1;
  int gridtype1 = gridInqType(gridID1);

  if ( gridtype1 == GRID_LONLAT || gridtype1 == GRID_GAUSSIAN || gridtype1 == GRID_PROJECTION )
    {
      gridID1 = gridToCurvilinear(gridID1, 0);
      gridtype1 = GRID_CURVILINEAR;
    }

  int gridID2 = -1;
  if ( gridtype1 == GRID_UNSTRUCTURED || gridtype1 == GRID_CURVILINEAR )
    gridID2 = gengridcell(gridID1, gridsize2, cellidx);

  if ( gridID0 != gridID1 ) gridDestroy(gridID1);

  return gridID2;
}

static
void sel_index(double *array1, double *array2, int nind, int *indarr)
{
  for ( int i = 0; i < nind; ++i )
    {
      array2[i] = array1[indarr[i]];
    }
}


void *Selindex(void *argument)
{
  int nrecs;
  int varID, levelID;
  int gridID1 = -1, gridID2;
  int index, gridtype = -1;
  int nmiss;
  double missval;
  typedef struct {
    int gridID1, gridID2;
  } sindex_t;
  lista_t *ilista = lista_new(INT_LISTA);

  cdoInitialize(argument);

  cdoOperatorAdd("selindex", 0, 0, "grid cell indices (1-N)");

  operatorInputArg(cdoOperatorEnter(0));

  int nind = args2int_lista(operatorArgc(), operatorArgv(), ilista);
  int *indarr = (int*) lista_dataptr(ilista);

  int indmin = indarr[0];
  int indmax = indarr[0];
  for ( int i = 1; i < nind; i++ )
    {
      if ( indmax < indarr[i] ) indmax = indarr[i];
      if ( indmin > indarr[i] ) indmin = indarr[i];
    }

  if ( cdoVerbose )
    for ( int i = 0; i < nind; i++ )
      cdoPrint("int %d = %d", i+1, indarr[i]);

  if ( indmin < 1 ) cdoAbort("Index < 1 not allowed!");

  for ( int i = 0; i < nind; i++ ) indarr[i] -= 1;


  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int nvars = vlistNvars(vlistID1);
  bool *vars  = (bool*) Malloc(nvars*sizeof(bool));
  for ( varID = 0; varID < nvars; varID++ ) vars[varID] = false;

  int ngrids = vlistNgrids(vlistID1);
  sindex_t *sindex = (sindex_t *) Malloc(ngrids*sizeof(sindex_t));

  for ( index = 0; index < ngrids; index++ )
    {
      gridID1  = vlistGrid(vlistID1, index);
      gridtype = gridInqType(gridID1);

      int gridsize = gridInqSize(gridID1);
      if ( gridsize == 1 ) continue;
      if ( indmax > gridsize )
        {
          cdoWarning("Max grid index is greater than grid size, skipped grid %d!", index+1);
          continue;
        }

      gridID2 = genindexgrid(gridID1, nind, indarr);

      if ( gridID2 == -1 )
        {
          cdoWarning("Unsupported grid type >%s<, skipped grid %d!", gridNamePtr(gridtype), index+1);
          continue;
        }
	  
      sindex[index].gridID1 = gridID1;
      sindex[index].gridID2 = gridID2;

      vlistChangeGridIndex(vlistID2, index, gridID2);

      for ( varID = 0; varID < nvars; varID++ )
        if ( gridID1 == vlistInqVarGrid(vlistID1, varID) )
          vars[varID] = true;
    }

  for ( varID = 0; varID < nvars; varID++ )
    if ( vars[varID] ) break;

  if ( varID >= nvars ) cdoAbort("No variables selected!");

  
  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
  double *array1 = (double*) Malloc(gridsize*sizeof(double));

  int gridsize2 = vlistGridsizeMax(vlistID2);
  if ( vlistNumber(vlistID2) != CDI_REAL ) gridsize2 *= 2;
  double *array2 = (double*) Malloc(gridsize2*sizeof(double));

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss);

	  streamDefRecord(streamID2, varID, levelID);

	  if ( vars[varID] )
	    {	      
	      gridID1 = vlistInqVarGrid(vlistID1, varID);

	      for ( index = 0; index < ngrids; index++ )
		if ( gridID1 == sindex[index].gridID1 ) break;

	      if ( index == ngrids ) cdoAbort("Internal problem, grid not found!");

	      gridsize2 = gridInqSize(sindex[index].gridID2);

              sel_index(array1, array2, nind, indarr); 

	      if ( nmiss )
		{
		  nmiss = 0;
		  missval = vlistInqVarMissval(vlistID2, varID);
		  for ( int i = 0; i < gridsize2; i++ )
		    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;
		}

	      streamWriteRecord(streamID2, array2, nmiss);
	    }
	  else
	    {
	      streamWriteRecord(streamID2, array1, nmiss);
	    }
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  vlistDestroy(vlistID2);

  if ( vars   ) Free(vars);
  if ( array2 ) Free(array2);
  if ( array1 ) Free(array1);
  if ( sindex ) Free(sindex);

  lista_destroy(ilista);

  cdoFinish();

  return 0;
}
