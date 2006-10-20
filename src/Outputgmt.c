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

/*
   This module contains the following operators:

*/
/*
  Output center or bounderies for GMT plotting

  Plotting example:

    - outputcenter
    - outputbounds:
    - outputboundscpt:
*/

#if  defined  (HAVE_CONFIG_H)
#  include "config.h" /* VERSION */
#endif

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "color.h"


void *Outputgmt(void *argument)
{
  static char func[] = "Outputgmt";
  int OUTPUTCENTER, OUTPUTBOUNDS, OUTPUTBOUNDSCPT;
  int operatorID;
  int i;
  int varID0, varID, recID;
  int gridsize = 0;
  int gridID, code;
  int nrecs;
  int levelID;
  int tsID;
  int streamID = 0;
  int vlistID;
  int nmiss;
  int nlon, nlat, nalloc;
  int nlev, lzon = FALSE, lmer = FALSE;
  int gridcorners = 0, ic;
  int status;
  int lgrid_gen_bounds = FALSE, luse_grid_corner = FALSE;
  int zaxisID, taxisID;
  int vdate, vtime;
  int year, month, day, hour, minute;
  char varname[256];
  double level;
  double missval;
  double *array = NULL;
  double *grid_center_lat = NULL, *grid_center_lon = NULL;
  double *grid_corner_lat = NULL, *grid_corner_lon = NULL;
  FILE *cpt_fp;
  CPT cpt;

  cdoInitialize(argument);

  OUTPUTCENTER    = cdoOperatorAdd("outputcenter",    0, 0, NULL);
  OUTPUTBOUNDS    = cdoOperatorAdd("outputbounds",    0, 0, NULL);
  OUTPUTBOUNDSCPT = cdoOperatorAdd("outputboundscpt", 0, 0, NULL);

  operatorID = cdoOperatorID();

  if ( operatorID == OUTPUTBOUNDS || operatorID == OUTPUTBOUNDSCPT )
    luse_grid_corner = TRUE;

  if ( operatorID == OUTPUTBOUNDSCPT )
    {
      char *cpt_file;

      cpt_file = operatorArgv()[0];

      if ( (cpt_fp = fopen (cpt_file, "r")) == NULL )
	cdoAbort("Open failed on color palette table %s", cpt_file);

      status = cptRead(cpt_fp, &cpt);
      if ( status != 0 )
	cdoAbort("Error during read of color palette table %s", cpt_file);
      
      if ( cdoVerbose ) cptWrite(stderr, cpt);
    }

  streamID = streamOpenRead(cdoStreamName(0));
  if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(0));

  vlistID = streamInqVlist(streamID);
  taxisID = vlistInqTaxis(vlistID);

  varID = 0;
  vlistInqVarName(vlistID, varID, varname);
  code    = vlistInqVarCode(vlistID, varID);
  gridID  = vlistInqVarGrid(vlistID, varID);
  zaxisID = vlistInqVarZaxis(vlistID, varID);
  missval = vlistInqVarMissval(vlistID, varID);

  if ( gridInqType(gridID) != GRID_LONLAT      &&
       gridInqType(gridID) != GRID_GAUSSIAN    &&
       gridInqType(gridID) != GRID_GME         &&
       gridInqType(gridID) != GRID_CURVILINEAR &&
       gridInqType(gridID) != GRID_CELL )
    cdoAbort("Output of %s data failed!", gridNamePtr(gridInqType(gridID)));
  
  if ( gridInqType(gridID) != GRID_CELL && gridInqType(gridID) != GRID_CURVILINEAR )
    {
      if ( gridInqType(gridID) == GRID_GME )
	{
	  gridID = gridToCell(gridID);
	}
      else
	{
	  gridID = gridToCurvilinear(gridID);
	  lgrid_gen_bounds = TRUE;
	}
    }

  gridsize = gridInqSize(gridID);
  nlon     = gridInqXsize(gridID);
  nlat     = gridInqYsize(gridID);
  nlev     = zaxisInqSize(zaxisID);

  if ( gridInqType(gridID) != GRID_CELL )
    {
      if ( nlon == 1 && nlat  > 1 && nlev > 1 ) lzon = TRUE;
      if ( nlon  > 1 && nlat == 1 && nlev > 1 ) lmer = TRUE;
    }

  if ( cdoVerbose && lzon ) cdoPrint("Process zonal data");
  if ( cdoVerbose && lmer ) cdoPrint("Process meridional data");

  if ( lzon || lmer ) 
    {
      if ( operatorID == OUTPUTBOUNDS || operatorID == OUTPUTBOUNDSCPT )
	cdoAbort("Bounds not available for zonal/meridional data!");
    }

  if ( gridInqType(gridID) == GRID_CELL )
    gridcorners = gridInqNvertex(gridID);
  else
    gridcorners = 4;

  grid_center_lat = (double *) realloc(grid_center_lat, gridsize*sizeof(double));
  grid_center_lon = (double *) realloc(grid_center_lon, gridsize*sizeof(double));

  gridInqYvals(gridID, grid_center_lat);
  gridInqXvals(gridID, grid_center_lon);

  if ( luse_grid_corner )
    {
      if ( gridcorners == 0 ) cdoAbort("grid corner missing!");
      nalloc = gridcorners*gridsize;
      grid_corner_lat = (double *) realloc(grid_corner_lat, nalloc*sizeof(double));
      grid_corner_lon = (double *) realloc(grid_corner_lon, nalloc*sizeof(double));

      if ( gridInqYbounds(gridID, NULL) && gridInqXbounds(gridID, NULL) )
	{
	  gridInqYbounds(gridID, grid_corner_lat);
	  gridInqXbounds(gridID, grid_corner_lon);
	}
      else
	{
	  if ( lgrid_gen_bounds )
	    {
	      genXbounds(nlon, nlat, grid_center_lon, grid_corner_lon);
	      genYbounds(nlon, nlat, grid_center_lat, grid_corner_lat);
	    }
	  else
	    cdoAbort("grid corner missing!");
	}
    }

  array = (double *) malloc(gridsize*sizeof(double));

  tsID = 0;
  nrecs = streamInqTimestep(streamID, tsID);
  /*  while ( (nrecs = streamInqTimestep(streamID, tsID)) ) */
    {
      vdate = taxisInqVdate(taxisID);
      vtime = taxisInqVtime(taxisID);

      decode_date(vdate, &year, &month, &day);
      decode_time(vtime, &hour, &minute);

      fprintf(stdout, "# Generated by CDO version %s\n", VERSION);
      fprintf(stdout, "#\n");
      fprintf(stdout, "# File  = %s\n", cdoStreamName(0));
      fprintf(stdout, "# Date  = %4.4d-%2.2d-%2.2d\n", year, month, day);
      fprintf(stdout, "# Time  = %2.2d:%2.2d\n", hour, minute);
      fprintf(stdout, "# Name  = %s\n", varname);
      fprintf(stdout, "# Code  = %d\n", code);

      varID0 = varID;

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID, &varID, &levelID);

	  if ( varID != varID0 ) continue;
	  if ( recID > 0 && !lzon && !lmer ) continue;

	  level   = zaxisInqLevel(zaxisID, levelID);

	  fprintf(stdout, "# Level = %g\n", level);
	  fprintf(stdout, "#\n");

	  streamReadRecord(streamID, array, &nmiss);

	  if ( operatorID == OUTPUTCENTER )
	    {
	      for ( i = 0; i < gridsize; i++ )
		{
		  /* if ( !DBL_IS_EQUAL(array[i], missval) ) */
		    {
		      if ( lzon )
			fprintf(stdout, " %g  %g  %g\n", grid_center_lat[i], level, array[i]);
		      else if ( lmer )
			fprintf(stdout, " %g  %g  %g\n", grid_center_lon[i], level, array[i]);
		      else
			fprintf(stdout, " %g  %g  %g\n", grid_center_lon[i], grid_center_lat[i], array[i]);
		    }
		}
	      fprintf(stdout, "#\n");
	    }
	  else if ( operatorID == OUTPUTBOUNDS || operatorID == OUTPUTBOUNDSCPT )
	    {
	      for ( i = 0; i < gridsize; i++ )
		{
		  if ( !DBL_IS_EQUAL(array[i], missval) )
		    fprintf(stdout, "> -Z%g", array[i]);
		  else
		    fprintf(stdout, ">");

		  if ( operatorID == OUTPUTBOUNDSCPT )
		    {
		      int r = 0, g = 0, b = 0, n;

		      if ( !DBL_IS_EQUAL(array[i], missval) )
			{
			  for ( n = 0; n < cpt.ncolors; n++ )
			    if ( array[i] > cpt.lut[n].z_low && array[i] <= cpt.lut[n].z_high ) break;

			  if ( n == cpt.ncolors )
			    {
			      r = cpt.bfn[0].rgb[0];  g = cpt.bfn[0].rgb[1];  b = cpt.bfn[0].rgb[2];
			    }
			  else
			    {
			      r = cpt.lut[n].rgb_high[0];  g = cpt.lut[n].rgb_high[1];  b = cpt.lut[n].rgb_high[2];
			    }
			}
		      else
			{
			  r = cpt.bfn[2].rgb[0];  g = cpt.bfn[2].rgb[1];  b = cpt.bfn[2].rgb[2]; 
			}

		      fprintf(stdout, " -G%d/%d/%d", r, g, b);
		    }

		  fprintf(stdout, "\n");

		  for ( ic = 0; ic < gridcorners; ic++ )
		    fprintf(stdout, "   %g  %g\n",
			    grid_corner_lon[i*gridcorners+ic], grid_corner_lat[i*gridcorners+ic]);
		  fprintf(stdout, "   %g  %g\n",
			  grid_corner_lon[i*gridcorners+0], grid_corner_lat[i*gridcorners+0]);
		}
	      fprintf(stdout, "\n");
	    }
	}
      tsID++;
    }

  streamClose(streamID);

  if ( array ) free(array);
  if ( grid_center_lon ) free(grid_center_lon);
  if ( grid_center_lat ) free(grid_center_lat);
  if ( grid_corner_lon ) free(grid_corner_lon);
  if ( grid_corner_lat ) free(grid_corner_lat);

  cdoFinish();

  return (0);
}
