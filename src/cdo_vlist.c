/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2005 Uwe Schulzweida, schulzweida@dkrz.de
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
#include "error.h"
#include "functs.h"

void vlistCompare(int vlistID1, int vlistID2, int function)
{
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
	  cdoAbort("Input files have different structure!");

	if ( gridInqSize(vlistInqVarGrid(vlistID1, varID)) !=
	     gridInqSize(vlistInqVarGrid(vlistID2, varID)) )
	  cdoAbort("Grid size of the input fields do not match!");
      }
  else if ( function == func_code )
    for ( varID = 0; varID < nvars; varID++ )
      {
	if ( vlistInqVarCode(vlistID1, varID) != vlistInqVarCode(vlistID2, varID) )
	  cdoAbort("Input files have different structure!");

	if ( zaxisInqSize(vlistInqVarZaxis(vlistID1, varID)) !=
	     zaxisInqSize(vlistInqVarZaxis(vlistID2, varID)) )
	  cdoAbort("Number of level of the input fields do not match!");
      }
  else if ( function == func_sft )
    for ( varID = 0; varID < nvars; varID++ )
      {
	if ( gridInqSize(vlistInqVarGrid(vlistID1, varID)) !=
	     gridInqSize(vlistInqVarGrid(vlistID2, varID)) )
	  cdoAbort("Grid size of the input fields do not match!");

	if ( zaxisInqSize(vlistInqVarZaxis(vlistID1, varID)) !=
	     zaxisInqSize(vlistInqVarZaxis(vlistID2, varID)) )
	  cdoAbort("Number of level of the input fields do not match!");
      }
  else
    cdoAbort("Internal problem! Invalid function %d", function);
}
