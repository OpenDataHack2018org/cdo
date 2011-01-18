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

      Setgrid    setgrid         Set grid
      Setgrid    setgridtype     Set grid type
      Setgrid    setgridarea     Set grid area
      Setgrid    setgridmask     Set grid mask
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


void *Setgrid(void *argument)
{
  int SETGRID, SETGRIDTYPE, SETGRIDAREA, SETGRIDMASK, UNSETGRIDMASK;
  int operatorID;
  int streamID1, streamID2 = CDI_UNDEFID;
  int nrecs;
  int tsID, recID, varID, levelID;
  int vlistID1, vlistID2;
  int taxisID1, taxisID2;
  int gridID1, gridID2 = -1;
  int ngrids, index;
  int gridtype = -1;
  int nmiss;
  int found;
  long i, gridsize;
  long areasize = 0;
  long masksize = 0;
  int lregular = 0;
  char *gridname = NULL;
  double *gridmask = NULL;
  double *areaweight = NULL;
  double *array = NULL;

  cdoInitialize(argument);

  SETGRID       = cdoOperatorAdd("setgrid",       0, 0, "grid description file or name");
  SETGRIDTYPE   = cdoOperatorAdd("setgridtype",   0, 0, "grid type");
  SETGRIDAREA   = cdoOperatorAdd("setgridarea",   0, 0, "filename with area weights");
  SETGRIDMASK   = cdoOperatorAdd("setgridmask",   0, 0, "filename with grid mask");
  UNSETGRIDMASK = cdoOperatorAdd("unsetgridmask", 0, 0, NULL);

  operatorID = cdoOperatorID();

  if ( operatorID != UNSETGRIDMASK )
    operatorInputArg(cdoOperatorEnter(operatorID));  

  if ( operatorID == SETGRID )
    {
      gridID2 = cdoDefineGrid(operatorArgv()[0]);
    }
  else if ( operatorID == SETGRIDTYPE )
    {
      gridname = operatorArgv()[0];

      if      ( strcmp(gridname, "curvilinear") == 0 ) gridtype = GRID_CURVILINEAR;
      else if ( strcmp(gridname, "cell") == 0 )        gridtype = GRID_CELL;
      else if ( strcmp(gridname, "lonlat") == 0 )      gridtype = GRID_LONLAT;
      else if ( strcmp(gridname, "gaussian") == 0 )    gridtype = GRID_GAUSSIAN;
      else if ( strcmp(gridname, "regular") == 0 )    {gridtype = GRID_GAUSSIAN; lregular = 1;}
      else cdoAbort("Unsupported grid name: %s", gridname);
    }
  else if ( operatorID == SETGRIDAREA )
    {
      int streamID, vlistID, gridID;
      char *areafile;

      areafile = operatorArgv()[0];

      streamID = streamOpenRead(areafile);

      vlistID = streamInqVlist(streamID);

      nrecs = streamInqTimestep(streamID, 0);
      streamInqRecord(streamID, &varID, &levelID);

      gridID = vlistInqVarGrid(vlistID, varID);
      areasize = gridInqSize(gridID);
      areaweight = (double *) malloc(areasize*sizeof(double));
  
      streamReadRecord(streamID, areaweight, &nmiss);

      streamClose(streamID);

      if ( cdoVerbose )
	{
	  double arrmean, arrmin, arrmax;

	  arrmean = areaweight[0];
	  arrmin  = areaweight[0];
	  arrmax  = areaweight[0];
	  for ( i = 1; i < areasize; i++ )
	    {
	      if ( areaweight[i] < arrmin ) arrmin = areaweight[i];
	      if ( areaweight[i] > arrmax ) arrmax = areaweight[i];
	      arrmean += areaweight[i];
	    }
	  arrmean = arrmean/areasize;

	  cdoPrint("areaweights: %d %#12.5g%#12.5g%#12.5g", areasize, arrmin, arrmean, arrmax);
	}
    }
  else if ( operatorID == SETGRIDMASK )
    {
      int streamID, vlistID, gridID;
      char *maskfile;
      double missval;

      maskfile = operatorArgv()[0];
      streamID = streamOpenRead(maskfile);

      vlistID = streamInqVlist(streamID);

      nrecs = streamInqTimestep(streamID, 0);
      streamInqRecord(streamID, &varID, &levelID);

      missval  = vlistInqVarMissval(vlistID, varID);
      gridID   = vlistInqVarGrid(vlistID, varID);
      masksize = gridInqSize(gridID);
      gridmask = (double *) malloc(masksize*sizeof(double));
  
      streamReadRecord(streamID, gridmask, &nmiss);

      streamClose(streamID);

      for ( i = 0; i < masksize; i++ )
	if ( DBL_IS_EQUAL(gridmask[i], missval) ) gridmask[i] = 0;
    }

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  if ( operatorID == SETGRID )
    {
      found = 0;
      ngrids = vlistNgrids(vlistID1);
      for ( index = 0; index < ngrids; index++ )
	{
	  gridID1 = vlistGrid(vlistID1, index);

	  if ( gridInqSize(gridID1) == gridInqSize(gridID2) )
	    {
	      vlistChangeGridIndex(vlistID2, index, gridID2);
	      found++;
	    }
	}
      if ( ! found ) cdoWarning("No grid with %d points found!", gridInqSize(gridID2));
    }
  else if ( operatorID == SETGRIDTYPE )
    {
      ngrids = vlistNgrids(vlistID1);
      for ( index = 0; index < ngrids; index++ )
	{
	  gridID1 = vlistGrid(vlistID1, index);

	  if ( lregular )
	    {
	      if ( gridInqType(gridID1) == GRID_GAUSSIAN_REDUCED )
		{
		  gridID2 = gridToRegular(gridID1);
		}
	    }
	  else
	    {
	      if      ( gridtype == GRID_CURVILINEAR ) gridID2 = gridToCurvilinear(gridID1);
	      else if ( gridtype == GRID_CELL )        gridID2 = gridToCell(gridID1);
	      else cdoAbort("Unsupported grid name: %s", gridname);
	    }

	  /*	  gridCompress(gridID2); */
	  vlistChangeGridIndex(vlistID2, index, gridID2);
	}
    }
  else if ( operatorID == SETGRIDAREA )
    {
      ngrids = vlistNgrids(vlistID1);
      for ( index = 0; index < ngrids; index++ )
	{
	  gridID1  = vlistGrid(vlistID1, index);
	  gridsize = gridInqSize(gridID1);
	  if ( gridsize == areasize )
	    {
	      gridID2 = gridDuplicate(gridID1);
	      gridDefArea(gridID2, areaweight);
	      vlistChangeGridIndex(vlistID2, index, gridID2);
	    }
	}
    }
  else if ( operatorID == SETGRIDMASK )
    {
      ngrids = vlistNgrids(vlistID1);
      for ( index = 0; index < ngrids; index++ )
	{
	  gridID1  = vlistGrid(vlistID1, index);
	  gridsize = gridInqSize(gridID1);
	  if ( gridsize == masksize )
	    {
	      int *mask;
	      mask = (int *) malloc(masksize*sizeof(int));
	      for ( i = 0; i < masksize; i++ )
		{
		  if ( gridmask[i] < 0 || gridmask[i] > 255 )
		    mask[i] = 0;
		  else
		    mask[i] = NINT(gridmask[i]);
		}
	      gridID2 = gridDuplicate(gridID1);
	      gridDefMask(gridID2, mask);
	      vlistChangeGridIndex(vlistID2, index, gridID2);
	      free(mask);
	    }
	}
    }
  else if ( operatorID == UNSETGRIDMASK )
    {
      ngrids = vlistNgrids(vlistID1);
      for ( index = 0; index < ngrids; index++ )
	{
	  gridID1  = vlistGrid(vlistID1, index);
	  gridID2 = gridDuplicate(gridID1);
	  gridDefMask(gridID2, NULL);
	  vlistChangeGridIndex(vlistID2, index, gridID2);
	}
    }

  streamDefVlist(streamID2, vlistID2);
  //vlistPrint(vlistID2);

  if ( lregular )
    gridsize = vlistGridsizeMax(vlistID2);
  else
    gridsize = vlistGridsizeMax(vlistID1);

  array = (double *) malloc(gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamDefRecord(streamID2,  varID,  levelID);
	  
	  streamReadRecord(streamID1, array, &nmiss);
	  if ( lregular )
	    {
	      gridID1 = vlistInqVarGrid(vlistID1, varID);
	      gridID2 = vlistInqVarGrid(vlistID2, varID);
	      if ( gridInqType(gridID1) == GRID_GAUSSIAN_REDUCED )
		{
		  double missval = vlistInqVarMissval(vlistID1, varID);
		  field2regular(gridID1, gridID2, missval, array, nmiss);
		}
	    }

	  streamWriteRecord(streamID2, array, nmiss);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( gridmask ) free(gridmask);
  if ( areaweight ) free(areaweight);
  if ( array ) free(array);

  cdoFinish();

  return (0);
}
