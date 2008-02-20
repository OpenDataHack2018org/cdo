/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2008 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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
    - outputbounds
    - outputboundscpt
    - outputvector
*/

#if  defined  (HAVE_CONFIG_H)
#  include "config.h" /* VERSION */
#endif

#include <string.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "color.h"



int pnpoly(int npol, double *xp, double *yp, double x, double y)
{
  int i, j, c = 0;
  for (i = 0, j = npol-1; i < npol; j = i++) {
    if ((((yp[i]<=y) && (y<yp[j])) ||
	 ((yp[j]<=y) && (y<yp[i]))) &&
	(x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]))
      
      c = !c;
  }
  return c;
}



void verify_grid(int gridtype, int gridsize, int xsize, int ysize, int ncorner,
		double *grid_center_lon, double *grid_center_lat,
		double *grid_corner_lon, double *grid_corner_lat)
{
  int i, k;
  int nout = 0;
  int isinside;
  double lon, lat;
  double lon_bounds[9], lat_bounds[9];

  /* check that all centers are inside the bounds */

  for ( i = 0; i < gridsize; ++i )
    {
      lon = grid_center_lon[i];
      lat = grid_center_lat[i];

      for ( k = 0; k < ncorner; ++k )
	{
	  lon_bounds[k] = grid_corner_lon[i*ncorner+k];
	  lat_bounds[k] = grid_corner_lat[i*ncorner+k];
	}

      for ( k = 0; k < ncorner; ++k )
	{
	  if ( (grid_center_lon[i] - lon_bounds[k]) > 270 ) lon_bounds[k] += 360;
	  if ( (lon_bounds[k] - grid_center_lon[i]) > 270 ) lon_bounds[k] -= 360;
	}

      lon_bounds[ncorner] = lon_bounds[0];
      lat_bounds[ncorner] = lat_bounds[0];

      isinside = pnpoly(ncorner+1, lon_bounds, lat_bounds, lon, lat);

      if ( !isinside ) nout++;

      if ( !isinside && cdoVerbose )
	printf("%d %d %g %g %g %g %g %g %g %g %g %g\n", nout, i, lon, lat, lon_bounds[0], lat_bounds[0],
	       lon_bounds[1], lat_bounds[1], lon_bounds[2], lat_bounds[2], lon_bounds[3], lat_bounds[3]);
    }

  if ( nout > 0 )
    cdoWarning("%d of %d points out of bounds!\n", nout, gridsize);
}


void make_cyclic(double *array1, double *array2, int nlon, int nlat)
{
  int i, j;
  int ij1, ij2;

  for ( j = 0; j < nlat; ++j )
    {
      for ( i = 0; i < nlon; ++i )
	{
	  ij1 = j*nlon+i;
	  ij2 = j*(nlon+1)+i;
	  array2[ij2] = array1[ij1];
	}
    }

  for ( j = 0; j < nlat; ++j )
    {
      ij2 = j*(nlon+1);
      array2[ij2+nlon] = array2[ij2];
    }
}


void *Outputgmt(void *argument)
{
  static char func[] = "Outputgmt";
  int GRIDVERIFY, OUTPUTCENTER, OUTPUTCENTER2, OUTPUTCENTERCPT, OUTPUTBOUNDS;
  int OUTPUTBOUNDSCPT, OUTPUTVECTOR, OUTPUTTRI;
  int operatorID;
  int process_data = TRUE;
  int i, j;
  int varID0, varID, recID;
  int nvals;
  int gridsize = 0;
  int gridsize2 = 0;
  int gridID, code;
  int nrecs;
  int levelID;
  int tsID;
  int streamID = 0;
  int vlistID;
  int nmiss;
  int nlon, nlat, nalloc;
  int nlev, lzon = FALSE, lmer = FALSE, lhov = FALSE;
  int gridcorners = 0, ic;
  int status;
  int lgrid_gen_bounds = FALSE, luse_grid_corner = FALSE;
  int zaxisID, taxisID;
  int ninc = 1;
  int vdate, vtime;
  int year, month, day, hour, minute;
  char varname[256];
  double level;
  double missval;
  double *array = NULL;
  double *array2 = NULL;
  double *parray;
  double *uf = NULL, *vf = NULL, *alpha = NULL, *auv = NULL;
  double *grid_center_lat = NULL, *grid_center_lon = NULL;
  double *grid_center_lat2 = NULL, *grid_center_lon2 = NULL;
  double *grid_corner_lat = NULL, *grid_corner_lon = NULL;
  double *plat, *plon;
  double *zaxis_center_lev, *zaxis_lower_lev, *zaxis_upper_lev;
  FILE *cpt_fp;
  CPT cpt;
  int grid_is_circular;

  cdoInitialize(argument);

  GRIDVERIFY      = cdoOperatorAdd("gridverify",      0, 0, NULL);
  OUTPUTCENTER    = cdoOperatorAdd("outputcenter",    0, 0, NULL);
  OUTPUTCENTER2   = cdoOperatorAdd("outputcenter2",   0, 0, NULL);
  OUTPUTCENTERCPT = cdoOperatorAdd("outputcentercpt", 0, 0, NULL);
  OUTPUTBOUNDS    = cdoOperatorAdd("outputbounds",    0, 0, NULL);
  OUTPUTBOUNDSCPT = cdoOperatorAdd("outputboundscpt", 0, 0, NULL);
  OUTPUTVECTOR    = cdoOperatorAdd("outputvector",    0, 0, NULL);
  OUTPUTTRI       = cdoOperatorAdd("outputtri",       0, 0, NULL);

  operatorID = cdoOperatorID();

  if ( operatorID == OUTPUTVECTOR )
    {
      operatorInputArg("increment");
      operatorCheckArgc(1);
      ninc  = atoi(operatorArgv()[0]);
      if ( ninc < 1 ) cdoAbort("Increment must be greater than 0!");
    }

  if ( operatorID == GRIDVERIFY  )
    {
      process_data = FALSE;
      luse_grid_corner = TRUE;
    }

  if ( operatorID == OUTPUTBOUNDS || operatorID == OUTPUTBOUNDSCPT )
    luse_grid_corner = TRUE;

  if ( operatorID == OUTPUTCENTERCPT || operatorID == OUTPUTBOUNDSCPT )
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
      if ( nlon == 1 && nlat  > 1 && nlev == 1 ) lhov = TRUE;
      if ( nlon == 1 && nlat  > 1 && nlev  > 1 ) lzon = TRUE;
      if ( nlon  > 1 && nlat == 1 && nlev  > 1 ) lmer = TRUE;
    }
  else
    {
      nlat = 1;
    }

  if ( cdoVerbose && lhov ) cdoPrint("Process hovmoeller data");
  if ( cdoVerbose && lzon ) cdoPrint("Process zonal data");
  if ( cdoVerbose && lmer ) cdoPrint("Process meridional data");
  /*
  if ( lzon || lmer ) 
    {
      if ( operatorID == OUTPUTBOUNDS || operatorID == OUTPUTBOUNDSCPT )
	cdoAbort("Bounds not available for zonal/meridional data!");
    }
  */
  if ( lhov ) 
    {
      if ( operatorID == OUTPUTBOUNDS || operatorID == OUTPUTBOUNDSCPT )
	cdoAbort("Bounds not available hovmoeller data!");
    }

  if ( gridInqType(gridID) == GRID_CELL )
    gridcorners = gridInqNvertex(gridID);
  else
    gridcorners = 4;

  grid_is_circular = gridIsCircular(gridID);

  grid_center_lat = (double *) malloc(gridsize*sizeof(double));
  grid_center_lon = (double *) malloc(gridsize*sizeof(double));

  gridInqYvals(gridID, grid_center_lat);
  gridInqXvals(gridID, grid_center_lon);

  nvals = gridsize;
  plon = grid_center_lon;
  plat = grid_center_lat;

  if ( operatorID == OUTPUTCENTER2 && grid_is_circular )
    {
      int ij2;

      gridsize2 = nlat*(nlon+1);

      grid_center_lat2 = (double *) malloc(gridsize2*sizeof(double));
      grid_center_lon2 = (double *) malloc(gridsize2*sizeof(double));

      make_cyclic(grid_center_lat, grid_center_lat2, nlon, nlat);
      make_cyclic(grid_center_lon, grid_center_lon2, nlon, nlat);

      for ( j = 0; j < nlat; ++j )
	{
	  ij2 = j*(nlon+1);
	  grid_center_lon2[ij2+nlon] += 360;
	}

      nvals = gridsize2;
      plon = grid_center_lon2;
      plat = grid_center_lat2;
    }

  zaxis_center_lev = (double *) malloc(nlev*sizeof(double));
  zaxis_lower_lev  = (double *) malloc(nlev*sizeof(double));
  zaxis_upper_lev  = (double *) malloc(nlev*sizeof(double));

  zaxisInqLevels(zaxisID, zaxis_center_lev);

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
	      if ( ! (lzon || lmer) )
		genXbounds(nlon, nlat, grid_center_lon, grid_corner_lon);
	      genYbounds(nlon, nlat, grid_center_lat, grid_corner_lat);
	    }
	  else
	    cdoAbort("grid corner missing!");
	}

      if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
	{
	  zaxisInqLbounds(zaxisID, zaxis_lower_lev);
	  zaxisInqUbounds(zaxisID, zaxis_upper_lev);
	}
      else
	{
	  zaxis_lower_lev[0] = zaxis_center_lev[0];
	  for ( i = 1; i < nlev; ++i )
	    zaxis_lower_lev[i] = 0.5*(zaxis_center_lev[i] + zaxis_center_lev[i-1]);

	  zaxis_upper_lev[nlev-1] = zaxis_center_lev[nlev-1];
	  for ( i = 0; i < nlev-1; ++i )
	    zaxis_upper_lev[i] = zaxis_lower_lev[i+1];

	  if ( cdoVerbose )
	    for ( i = 0; i < nlev; ++i )
	      printf("level: %d %g %g %g\n",
		     i+1, zaxis_lower_lev[i], zaxis_center_lev[i], zaxis_upper_lev[i]);
	}
    }

  array = (double *) malloc(gridsize*sizeof(double));
  parray = array;
						
  if ( operatorID == OUTPUTCENTER2 && grid_is_circular )
    {
      array2 = (double *) malloc(nlat*(nlon+1)*sizeof(double));
      parray = array2;
    }

  if ( operatorID == OUTPUTVECTOR )
    {
      uf    = (double *) malloc(gridsize*sizeof(double));
      vf    = (double *) malloc(gridsize*sizeof(double));
      alpha = (double *) malloc(gridsize*sizeof(double));
      auv   = (double *) malloc(gridsize*sizeof(double));
    }

  if ( operatorID == GRIDVERIFY )
    verify_grid(gridInqType(gridID), gridsize, nlon, nlat, gridcorners,
		grid_center_lon, grid_center_lat,
		grid_corner_lon, grid_corner_lat);

  tsID = 0;
  if ( process_data )
  while ( (nrecs = streamInqTimestep(streamID, tsID)) )
    {
      vdate = taxisInqVdate(taxisID);
      vtime = taxisInqVtime(taxisID);

      decode_date(vdate, &year, &month, &day);
      decode_time(vtime, &hour, &minute);

      if ( tsID == 0 && operatorID != OUTPUTTRI )
	{
#if defined (VERSION)
	  fprintf(stdout, "# Generated by CDO version %s\n", VERSION);
	  fprintf(stdout, "#\n");
#endif
	  fprintf(stdout, "# Operator = %s\n", cdoOperatorName(operatorID));
	  if      ( lhov )  fprintf(stdout, "# Mode     = hovmoeller\n");
	  else if ( lzon )  fprintf(stdout, "# Mode     = zonal\n");
	  else if ( lmer )  fprintf(stdout, "# Mode     = meridional\n");
	  else              fprintf(stdout, "# Mode     = horizonal\n");

	  if ( operatorID == OUTPUTVECTOR )
	    fprintf(stdout, "# Increment = %d\n", ninc);
	  fprintf(stdout, "#\n");
	  fprintf(stdout, "# File  = %s\n", cdoStreamName(0));
	  fprintf(stdout, "# Date  = %4.4d-%2.2d-%2.2d\n", year, month, day);
	  fprintf(stdout, "# Time  = %2.2d:%2.2d\n", hour, minute);
	  fprintf(stdout, "# Name  = %s\n", varname);
	  fprintf(stdout, "# Code  = %d\n", code);
	}

      varID0 = varID;

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID, &varID, &levelID);

	  if ( varID != varID0 ) continue;
	  if ( recID > 0 && !lzon && !lmer ) continue;

	  streamReadRecord(streamID, array, &nmiss);

	  if ( operatorID == OUTPUTCENTER2 && grid_is_circular )
	    make_cyclic(array, array2, nlon, nlat);

	  level = zaxis_center_lev[levelID];

	  if ( (tsID == 0 || lzon || lmer) && operatorID != OUTPUTTRI )
	    fprintf(stdout, "# Level = %g\n", level);
	  if ( lhov )
	    fprintf(stdout, "# Timestep = %d\n", tsID+1);

	  if ( operatorID != OUTPUTTRI ) fprintf(stdout, "#\n");

	  if ( operatorID == OUTPUTCENTER || operatorID == OUTPUTCENTER2 || operatorID == OUTPUTCENTERCPT )
	    {
	      for ( i = 0; i < nvals; i++ )
		{
		  if ( operatorID == OUTPUTCENTERCPT )
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
		    }

		  if ( operatorID == OUTPUTCENTER )
		    {
		      if ( lzon )
			fprintf(stdout, " %g  %g  %g\n", grid_center_lat[i], level, array[i]);
		      else if ( lmer )
			fprintf(stdout, " %g  %g  %g\n", grid_center_lon[i], level, array[i]);
		      else if ( lhov )
			fprintf(stdout, " %d  %g  %g\n", tsID+1, grid_center_lat[i], array[i]);
		      else
			fprintf(stdout, " %g  %g  %g\n", grid_center_lon[i], grid_center_lat[i], array[i]);
		    }
		  else if ( operatorID == OUTPUTCENTER2 )
		    {
		      fprintf(stdout, " %g  %g  %g\n", plon[i], plat[i], parray[i]);
		    }
		  else
		    {
		      if ( lzon )
			fprintf(stdout, " %g  %g  %g\n", grid_center_lat[i], level, array[i]);
		      else if ( lmer )
			fprintf(stdout, " %g  %g  %g\n", grid_center_lon[i], level, array[i]);
		      else
			fprintf(stdout, " %g  %g  %g  %g\n", grid_center_lon[i], grid_center_lat[i], array[i], array[i]);
		    }
		}
	      fprintf(stdout, "#\n");
	    }
	  else if ( operatorID == OUTPUTTRI )
	    {
	      int ij, c1, c2, c3;
	      int mlon, ip1;
	      if ( gridInqType(gridID) != GRID_CURVILINEAR ) cdoAbort("Unsupported grid!");

	      mlon = nlon-1;
	      /* if ( gridIsCircular(gridID) ) mlon = nlon; */
	      for ( j = 0; j < nlat-1; ++j )
		{
		  for ( i = 0; i < mlon; ++i )
		    {
		      ip1 = i+1;
		      if ( i == nlon-1 ) ip1 = 0;
		      ij = j*nlon+i;
		      c1 = (j)*nlon+ip1;
		      c2 = (j)*nlon+i;
		      c3 = (j+1)*nlon+i;
		      fprintf(stdout, "%d   %d   %d\n", c1, c2, c3);
		      c1 = (j)*nlon+i+1;
		      c2 = (j+1)*nlon+i;
		      c3 = (j+1)*nlon+ip1;
		      fprintf(stdout, "%d   %d   %d\n", c1, c2, c3);
		    }
		}
	    }
	  else if ( operatorID == OUTPUTVECTOR )
	    {
	      if ( nrecs < 2 ) cdoAbort("No enough fields!");

	      memcpy(uf, array, gridsize*sizeof(double));
	      streamInqRecord(streamID, &varID, &levelID);
	      streamReadRecord(streamID, vf, &nmiss);

	      for ( j = 0; j < nlat; j += ninc )
		for ( i = 0; i < nlon; i += ninc )
		  {
		    /* compute length of velocity vector */
		    auv[IX2D(j,i,nlon)] = sqrt(uf[IX2D(j,i,nlon)]*uf[IX2D(j,i,nlon)] + 
					       vf[IX2D(j,i,nlon)]*vf[IX2D(j,i,nlon)]);

		    alpha[IX2D(j,i,nlon)] = atan2(vf[IX2D(j,i,nlon)],uf[IX2D(j,i,nlon)]);
		    alpha[IX2D(j,i,nlon)] = 90. - alpha[IX2D(j,i,nlon)]*RAD2DEG;

		    if ( alpha[IX2D(j,i,nlon)] <   0 ) alpha[IX2D(j,i,nlon)] += 360;
		    if ( alpha[IX2D(j,i,nlon)] > 360 ) alpha[IX2D(j,i,nlon)] -= 360;

		    if ( fabs(auv[IX2D(j,i,nlon)]) > 0 )
		      fprintf(stdout, " %g  %g  %g  %g\n",
			      grid_center_lon[IX2D(j,i,nlon)], grid_center_lat[IX2D(j,i,nlon)],
			      alpha[IX2D(j,i,nlon)], auv[IX2D(j,i,nlon)]);
		  }

	      fprintf(stdout, "#\n");
	      break;
	    }
	  else if ( operatorID == OUTPUTBOUNDS || operatorID == OUTPUTBOUNDSCPT )
	    {
	      for ( i = 0; i < gridsize; i++ )
		{
		  if ( !DBL_IS_EQUAL(array[i], missval) )
		    fprintf(stdout, "> -Z%g", array[i]);
		  else
		    fprintf(stdout, "> -ZNaN");

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

		  if ( lzon )
		    {
		      double xlev[4];
		      xlev[0] = zaxis_lower_lev[levelID];
		      xlev[1] = zaxis_upper_lev[levelID];
		      xlev[2] = zaxis_upper_lev[levelID];
		      xlev[3] = zaxis_lower_lev[levelID];
		      for ( ic = 0; ic < 4; ic++ )
			fprintf(stdout, "   %g  %g\n",
				grid_corner_lat[i*4+ic], xlev[ic]);
		      fprintf(stdout, "   %g  %g\n",
			      grid_corner_lat[i*4], xlev[0]);
		    }
		  else if ( lmer )
		    {
		      cdoAbort("Implementation for meridional data missing!\n");
		    }
		  else if ( lhov )
		    {
		      cdoAbort("Implementation for hovmoeller data missing!\n");
		    }
		  else
		    {
		      for ( ic = 0; ic < gridcorners; ic++ )
			fprintf(stdout, "   %g  %g\n",
				grid_corner_lon[i*gridcorners+ic], grid_corner_lat[i*gridcorners+ic]);
		      fprintf(stdout, "   %g  %g\n",
			      grid_corner_lon[i*gridcorners], grid_corner_lat[i*gridcorners]);
		    }
		}
	      fprintf(stdout, "\n");
	    }
	}

      if ( ! lhov ) break;

      tsID++;
    }

  streamClose(streamID);

  if ( array  ) free(array);
  if ( array2 ) free(array2);
  if ( grid_center_lon ) free(grid_center_lon);
  if ( grid_center_lat ) free(grid_center_lat);
  if ( grid_center_lon2 ) free(grid_center_lon2);
  if ( grid_center_lat2 ) free(grid_center_lat2);
  if ( grid_corner_lon ) free(grid_corner_lon);
  if ( grid_corner_lat ) free(grid_corner_lat);

  free(zaxis_center_lev);
  free(zaxis_lower_lev);
  free(zaxis_upper_lev);

  cdoFinish();

  return (0);
}
