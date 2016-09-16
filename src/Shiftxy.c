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

      Selbox     sellonlatbox    Select lon/lat box
      Selbox     selindexbox     Select index box
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "pstream.h"


static
void shiftx(bool lcyclic, int nshift, int gridID, double *array1, double *array2, double missval)
{
  int nx = gridInqXsize(gridID);
  int ny = gridInqYsize(gridID);

  for ( int i = 0; i < nx; i++ )
    {
      bool is_cyclic = false;
      int ins = i + nshift%nx;
      while ( ins >= nx ) { ins -= nx; is_cyclic = true; }
      while ( ins <   0 ) { ins += nx; is_cyclic = true; }
        
      if ( !lcyclic && is_cyclic )
        {
          for ( int j = 0; j < ny; j++ )
            array2[IX2D(j,ins,nx)] = missval;
        }
      else
        {
          for ( int j = 0; j < ny; j++ )
            array2[IX2D(j,ins,nx)] = array1[IX2D(j,i,nx)];
        }
    }
}

static
void shifty(bool lcyclic, int nshift, int gridID, double *array1, double *array2, double missval)
{
  int nx = gridInqXsize(gridID);
  int ny = gridInqYsize(gridID);

  for ( int j = 0; j < ny; j++ )
    {
      bool is_cyclic = false;
      int jns = j + nshift%ny;

      while ( jns >= ny ) { jns -= ny; is_cyclic = true; }
      while ( jns <   0 ) { jns += ny; is_cyclic = true; }

      if ( !lcyclic && is_cyclic )
        {
          for ( int i = 0; i < nx; i++ )
            array2[IX2D(jns,i,nx)] = missval;
        }
      else
        {
          for ( int i = 0; i < nx; i++ )
            array2[IX2D(jns,i,nx)] = array1[IX2D(j,i,nx)];
        }
    }
}


void *Shiftxy(void *argument)
{
  int nrecs;
  int varID, levelID;
  int nmiss;

  cdoInitialize(argument);

  int SHIFTX = cdoOperatorAdd("shiftx", 0, 0, NULL);
  int SHIFTY = cdoOperatorAdd("shifty", 0, 0, NULL);

  int operatorID = cdoOperatorID();

  int nshift = 1;
  if ( operatorArgc() > 0 )
    {
      operatorCheckArgc(1);
      nshift = parameter2int(operatorArgv()[0]);
    }

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

  for ( int index = 0; index < ngrids; index++ )
    {
      int gridID1  = vlistGrid(vlistID1, index);
      int gridtype = gridInqType(gridID1);

      if ( gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN || gridtype == GRID_CURVILINEAR ||
           (gridtype == GRID_PROJECTION && gridInqProjType(gridID1) == CDI_PROJ_RLL) ||
	   (gridtype == GRID_GENERIC && gridInqXsize(gridID1) > 0 && gridInqYsize(gridID1) > 0) )
	{
          if ( operatorID == SHIFTX )
            {
              
            }
          /*
	  vlistChangeGridIndex(vlistID2, index, gridID2);
          */

	  for ( varID = 0; varID < nvars; varID++ )
	    if ( gridID1 == vlistInqVarGrid(vlistID1, varID) )
	      vars[varID] = true;
	}
      else if ( gridtype == GRID_GENERIC && gridInqXsize(gridID1) <= 1 && gridInqYsize(gridID1) <=1 )
	{
	}
      else
	{
	  cdoPrint("Unsupported grid type: %s", gridNamePtr(gridtype));
	  if ( gridtype == GRID_GAUSSIAN_REDUCED )
	    cdoPrint("Use option -R to convert Gaussian reduced grid to a regular grid!");
	  cdoAbort("Unsupported grid type!");
	}
    }

  for ( varID = 0; varID < nvars; varID++ )
    if ( vars[varID] ) break;

  if ( varID >= nvars ) cdoWarning("No variables selected!");

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);
  double *array1 = (double*) Malloc(gridsize*sizeof(double));
  double *array2 = (double*) Malloc(gridsize*sizeof(double));

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
	      int gridID1 = vlistInqVarGrid(vlistID1, varID);
	      int gridsize = gridInqSize(gridID1);
              double missval = vlistInqVarMissval(vlistID2, varID);
                  
              if      ( operatorID == SHIFTX ) shiftx(false, nshift, gridID1, array1, array2, missval);
              else if ( operatorID == SHIFTY ) shifty(false, nshift, gridID1, array1, array2, missval);

              nmiss = 0;
              for ( int i = 0; i < gridsize; i++ )
                if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;

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

  cdoFinish();

  return 0;
}
