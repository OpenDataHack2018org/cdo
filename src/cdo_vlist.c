/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2007 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "error.h"
#include "functs.h"

void vlistCompare(int vlistID1, int vlistID2, int function)
{
  static char func[] = "vlistCompare";
  int varID, nvars;

  if ( vlistNvars(vlistID1) != vlistNvars(vlistID2) )
    cdoAbort("Input streams have different number of variables per timestep!");

  if ( vlistNrecs(vlistID1) != vlistNrecs(vlistID2) )
    cdoAbort("Input streams have different number of records per timestep!");

  nvars = vlistNvars(vlistID1);

  if ( function == func_hrd )
    for ( varID = 0; varID < nvars; varID++ )
      {
	if ( vlistInqVarCode(vlistID1, varID) != vlistInqVarCode(vlistID2, varID) )
	  cdoAbort("Input streams have different structure!");

	if ( gridInqSize(vlistInqVarGrid(vlistID1, varID)) !=
	     gridInqSize(vlistInqVarGrid(vlistID2, varID)) )
	  cdoAbort("Grid size of the input fields do not match!");
      }
  else if ( function == func_code )
    for ( varID = 0; varID < nvars; varID++ )
      {
	if ( vlistInqVarCode(vlistID1, varID) != vlistInqVarCode(vlistID2, varID) )
	  cdoAbort("Input streams have different structure!");

	if ( zaxisInqSize(vlistInqVarZaxis(vlistID1, varID)) !=
	     zaxisInqSize(vlistInqVarZaxis(vlistID2, varID)) )
	  cdoAbort("Number of levels of the input fields do not match!");
      }
  else if ( function == func_sft )
    {
      for ( varID = 0; varID < nvars; varID++ )
	{
	  if ( gridInqSize(vlistInqVarGrid(vlistID1, varID)) !=
	       gridInqSize(vlistInqVarGrid(vlistID2, varID)) )
	    cdoAbort("Grid size of the input fields do not match!");

	  if ( zaxisInqSize(vlistInqVarZaxis(vlistID1, varID)) !=
	       zaxisInqSize(vlistInqVarZaxis(vlistID2, varID)) )
	    cdoAbort("Number of levels of the input fields do not match!");
	}
      
      /* compare grids of first variable */
      {
	int gridID1, gridID2;
	int xsize, ysize;

	gridID1 = vlistInqVarGrid(vlistID1, 0);
	gridID2 = vlistInqVarGrid(vlistID2, 0);

	if ( gridInqType(gridID1) == GRID_GAUSSIAN || gridInqType(gridID1) == GRID_LONLAT )
	  {
	    xsize = gridInqXsize(gridID1);
	    ysize = gridInqYsize(gridID1);

	    if ( ysize == gridInqYsize(gridID2) )
	      {
		if ( ysize > 1 )
		  {
		    double *yvals1, *yvals2;

		    yvals1 = (double *) malloc(ysize*sizeof(double));
		    yvals2 = (double *) malloc(ysize*sizeof(double));

		    gridInqYvals(gridID1, yvals1);
		    gridInqYvals(gridID2, yvals2);
		
		    if ( DBL_IS_EQUAL(yvals1[0], yvals2[ysize-1]) &&
			 DBL_IS_EQUAL(yvals1[ysize-1], yvals2[0]) )
		      {
			if ( yvals1[0] > yvals2[0] )
			  cdoWarning("Grid orientation differ! First grid: N->S; second grid: S->N");
			else
			  cdoWarning("Grid orientation differ! First grid: S->N; second grid: N->S");
		      }

		    free(yvals1);
		    free(yvals2);
		  }
	      }
	    else
	      cdoWarning("ysize of input grids differ!");
	  }
      }
    }
  else
    cdoAbort("Internal problem! Invalid function %d", function);
}


int vlistIsSzipped(int vlistID)
{
  int lszip = FALSE;
  int nvars, varID;

  nvars = vlistNvars(vlistID);

  for ( varID = 0; varID < nvars; varID++ )
    {
      if ( vlistInqVarSzip(vlistID, varID) )
	{
	  lszip = TRUE;
	  break;
	}
    }      

  return (lszip);
}


