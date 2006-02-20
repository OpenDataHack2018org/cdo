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

#include <string.h>
#include <math.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "remap.h"

/*
@BeginDoc

@BeginModule

@Name      = Remap
@Title     = SCRIP interpolation
@Section   = Interpolation
@Class     = Interpolation
@Arguments = ifile ofile
@Operators = remapcon remapbil remapbic remapdis gencon remap

@EndModule


@BeginOperator_remapcon

@Title     = Conservative remapping
@Parameter = grid

@BeginDesciption
SCRIP first order conservative remapping.
@EndDesciption

@BeginParameter
@Item = grid
STRING  Grid description file or name of the target grid
@EndParameter

@BeginEnvironment
@Item = NORMALIZE_OPT
This variable is used to choose the normalization
of the remapping. By default, NORMALIZE_OPT is set to be
'fracarea' and will include the destination area fraction in
the output weights; other options are 'none' and 'destarea'
(for more information see \cite{SCRIP}).
@EndEnvironment

@EndOperator


@BeginOperator_remapbil

@Title     = Bilinear interpolation
@Parameter = grid

@BeginDesciption
SCRIP bilinear interpolation (only rectangular grids).
@EndDesciption

@BeginParameter
@Item = grid
STRING  Grid description file or name of the target grid
@EndParameter

@EndOperator


@BeginOperator_remapbic

@Title     = Bicubic interpolation
@Parameter = grid

@BeginDesciption
SCRIP bicubic interpolation (only rectangular grids).
@EndDesciption

@BeginParameter
@Item = grid
STRING  Grid description file or name of the target grid
@EndParameter

@EndOperator


@BeginOperator_remapdis

@Title     = Distance-weighted averaging
@Parameter = grid

@BeginDesciption
SCRIP distance-weighted average of the four nearest neighbor values.
@EndDesciption

@BeginParameter
@Item = grid
STRING  Grid description file or name of the target grid
@EndParameter

@EndOperator


@BeginOperator_gencon

@Title     = Generate conservative interpolation weights
@Parameter = grid

@BeginDesciption
Generate SCRIP first order conservative interpolation weights and write the result to
a file.
@EndDesciption

@BeginParameter
@Item = grid
STRING  Grid description file or name of the target grid
@EndParameter

@BeginEnvironment
@Item = NORMALIZE_OPT
This variable is used to choose the normalization
of the remapping. By default, NORMALIZE_OPT is set to be
'fracarea' and will include the destination area fraction in
the output weights; other options are 'none' and 'destarea'
(for more information see \cite{SCRIP}).
@EndEnvironment

@EndOperator


@BeginOperator_remap

@Title     = Remapping
@Parameter = grid weights

@BeginDesciption
Remapping with the interpolation weights from a netCDF file.
The netCDF file must follow the SCRIP convention.
@EndDesciption

@BeginParameter weights
@Item = grid
STRING  Grid description file or name of the target grid
@Item = weights
STRING  Interpolation weights (SCRIP netCDF file)
@EndParameter

@EndOperator

@EndDoc
*/

void gme_grid_restore(double *p, int ni, int nd);

void *Remap(void *argument)
{
  static char func[] = "Remap";
  enum {REMAPCON, REMAPBIL, REMAPBIC, REMAPDIS, REMAPDIS1, REMAPCON1, 
        REMAPCONF, REMAPBILF, REMAPBICF, REMAPDISF, REMAPXXX};
  int operatorID;
  int operfunc;
  int streamID1, streamID2 = -1;
  int nrecs, ngrids;
  int nzaxis, zaxisID, zaxissize;
  int index;
  int tsID, recID, varID, levelID;
  int gridsize, gridsize2;
  int vlistID1, vlistID2;
  int taxisID1, taxisID2;
  int gridID1 = -1, gridID2;
  int nmiss1, nmiss2, i, j, r;
  int *imask = NULL;
  int nremaps = 0;
  int norm_opt = NORM_OPT_NONE;
  int map_type = -1;
  int max_remaps = 0;
  int lgridboxinfo = TRUE;
  char varname[128];
  double missval;
  double *array1 = NULL, *array2 = NULL;
  double *grad1_lat = NULL, *grad1_lon = NULL, *grad1_latlon = NULL;
  REMAP *remaps;
  char *envstring;
  char *remap_file = NULL;
  int lwrite_remap;

  cdoInitialize(argument);

  cdoOperatorAdd("remapcon",    REMAPCON,    0, NULL);
  cdoOperatorAdd("remapbil",    REMAPBIL,    0, NULL);
  cdoOperatorAdd("remapbic",    REMAPBIC,    0, NULL);
  cdoOperatorAdd("remapdis",    REMAPDIS,    0, NULL);
  cdoOperatorAdd("remapdis1",   REMAPDIS1,   0, NULL);
  cdoOperatorAdd("remapcon1",   REMAPCON1,   0, NULL);
  cdoOperatorAdd("gencon",      REMAPCONF,   1, NULL);
  cdoOperatorAdd("genbil",      REMAPBILF,   1, NULL);
  cdoOperatorAdd("genbic",      REMAPBICF,   1, NULL);
  cdoOperatorAdd("gendis",      REMAPDISF,   1, NULL);
  cdoOperatorAdd("remap",       REMAPXXX,    0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);
  lwrite_remap = cdoOperatorIntval(operatorID);

  envstring = getenv("MAX_REMAPS");
  if ( envstring )
    {
      int ival;
      ival = atoi(envstring);
      if ( ival > 0 )
	{
	  max_remaps = ival;
	  if ( cdoVerbose )
	    cdoPrint("Set MAX_REMAPS to %d", max_remaps);
	}
    }

  if ( operfunc == REMAPCON || operfunc == REMAPCON1 || operfunc == REMAPCONF )
    {
      norm_opt = NORM_OPT_FRACAREA;

      envstring = getenv("NORMALIZE_OPT");

      if ( envstring )
        {
	  if ( strcmp(envstring, "fracarea") == 0 )
	    norm_opt = NORM_OPT_FRACAREA;
	  else if ( strcmp(envstring, "destarea") == 0 )
	    norm_opt = NORM_OPT_DESTAREA;
	  else if ( strcmp(envstring, "none") == 0 )
	    norm_opt = NORM_OPT_NONE;
	  else
	    cdoWarning("NORMALIZE_OPT=%s unsupported!\n", envstring);
	}

      if ( cdoVerbose )
        {
	  if ( norm_opt == NORM_OPT_FRACAREA )
	    cdoPrint("Normalization option: fracarea");
	  else if ( norm_opt == NORM_OPT_DESTAREA )
	    cdoPrint("Normalization option: destarea");
	  else
	    cdoPrint("Normalization option: none");
	}
    }

  if ( operfunc == REMAPXXX )
    {
      operatorInputArg("grid description file or name, remap file (SCRIP netCDF)");
      operatorCheckArgc(2);
      gridID2 = cdoDefineGrid(operatorArgv()[0]);
      remap_file = operatorArgv()[1];
    }
  else
    {
      operatorInputArg("grid description file or name");
      operatorCheckArgc(1);
      gridID2 = cdoDefineGrid(operatorArgv()[0]);
    }

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  ngrids = vlistNgrids(vlistID1);
  for ( index = 0; index < ngrids; index++ )
    {
      gridID1 = vlistGrid(vlistID1, index);

      if ( gridInqType(gridID1) != GRID_LONLAT      &&
	   gridInqType(gridID1) != GRID_GAUSSIAN    &&
	   gridInqType(gridID1) != GRID_GME         &&
	   gridInqType(gridID1) != GRID_CURVILINEAR &&
	   gridInqType(gridID1) != GRID_CELL )
	cdoAbort("Remapping of %s data failed!", gridNamePtr(gridInqType(gridID1)) );

      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

  if ( max_remaps == 0 )
    {
      nzaxis = vlistNzaxis(vlistID1);
      for ( index = 0; index < nzaxis; index++ )
        {
	  zaxisID = vlistZaxis(vlistID1, index);
	  zaxissize = zaxisInqSize(zaxisID);
          if ( zaxissize > max_remaps ) max_remaps = zaxissize;
	}

      max_remaps++;

      if ( cdoVerbose )
        cdoPrint("Set max_remaps to %d", max_remaps);
    }

  remaps = (REMAP *) malloc(max_remaps*sizeof(REMAP));
  for ( r = 0; r < max_remaps; r++ )
    {
      remaps[r].gridID   = -1;
      remaps[r].gridsize = 0;
      remaps[r].nmiss    = 0;
    }

  if ( operfunc == REMAPXXX )
    {
      read_remap_scrip(remap_file, gridID1, gridID2, &map_type, &remaps[0].grid, &remaps[0].vars);
      nremaps = 1;
      gridsize = remaps[0].grid.grid1_size;
      remaps[0].gridID = gridID1;
      remaps[0].gridsize = gridInqSize(gridID1);
      remaps[0].nmiss = 0;
      if ( gridInqType(gridID1) == GRID_GME ) gridsize = remaps[0].grid.grid1_nvgp;
      if ( gridsize != remaps[0].gridsize )
	cdoAbort("Size of data grid and weights from %s differ!", remap_file);

      if ( gridInqType(gridID1) == GRID_GME ) gridsize = remaps[0].grid.grid1_size;

      for ( i = 0; i < gridsize; i++ )
        if ( remaps[0].grid.grid1_mask[i] == FALSE )
          remaps[0].nmiss++;

      if ( gridInqType(gridID2) == GRID_GME )
	{
	  int gridID2_gme;
	  remaps[0].grid.grid2_nvgp = gridInqSize(gridID2);
	  remaps[0].grid.grid2_vgpm = (int *) realloc(remaps[0].grid.grid2_vgpm,
						      gridInqSize(gridID2)*sizeof(int));
	  gridID2_gme = gridToCell(gridID2);
	  gridInqMask(gridID2_gme, remaps[0].grid.grid2_vgpm);
	}

      if ( map_type == MAP_TYPE_CONSERV )
        {
	  operfunc = REMAPCON;
	  cdoPrint("Using remapcon");
	}
      else if ( map_type == MAP_TYPE_BILINEAR )
        {
	  operfunc = REMAPBIL;
	  cdoPrint("Using remapbil");
	}
      else if ( map_type == MAP_TYPE_BICUBIC )
        {
	  operfunc = REMAPBIC;
	  cdoPrint("Using remapbic");
	}
      else if ( map_type == MAP_TYPE_DISTWGT )
        {
	  operfunc = REMAPDIS;
	  cdoPrint("Using remapdis");
	}
      else
	cdoAbort("unsupported mapping method (map_type = %d)", map_type);
    }

  switch ( operfunc )
    {
    case REMAPCON:
    case REMAPCON1:
    case REMAPCONF:
      map_type = MAP_TYPE_CONSERV;
      break;
    case REMAPBIL:
    case REMAPBILF:
      map_type = MAP_TYPE_BILINEAR;
      break;
    case REMAPBIC:
    case REMAPBICF:
      map_type = MAP_TYPE_BICUBIC;
      break;
    case REMAPDIS:
    case REMAPDISF:
      map_type = MAP_TYPE_DISTWGT;
      break;
    case REMAPDIS1:
      map_type = MAP_TYPE_DISTWGT1;
      break;
    default:
      cdoAbort("unknown mapping method");
    }

  gridsize = vlistGridsizeMax(vlistID1);

  if ( map_type == MAP_TYPE_BICUBIC )
    {
      grad1_lat    = (double *) malloc(gridsize*sizeof(double));
      grad1_lon    = (double *) malloc(gridsize*sizeof(double));
      grad1_latlon = (double *) malloc(gridsize*sizeof(double));
    }

  array1 = (double *) malloc(gridsize*sizeof(double));
  imask  = (int *) malloc(gridsize*sizeof(int));

  gridsize = gridInqSize(gridID2);
  array2   = (double *) malloc(gridsize*sizeof(double));

  if ( ! lwrite_remap )
    {
      streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
      if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

      streamDefVlist(streamID2, vlistID2);
    }

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      if ( ! lwrite_remap ) 
	streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss1);

	  gridID1 = vlistInqVarGrid(vlistID1, varID);

	  if ( map_type != MAP_TYPE_CONSERV && 
	       gridInqType(gridID1) == GRID_GME && gridInqType(gridID2) == GRID_GME )
	    cdoAbort("Only conservative remapping is available to remap between GME grids!");
	  /*
	  if ( gridIsRotated(gridID1) && map_type != MAP_TYPE_CONSERV )
	    cdoAbort("Only conservative remapping is available for rotated grids!");
	  */
	  missval = vlistInqVarMissval(vlistID1, varID);
	  gridsize = gridInqSize(gridID1);

	  for ( i = 0; i < gridsize; i++ )
	    if ( DBL_IS_EQUAL(array1[i], missval) )
	      imask[i] = FALSE;
	    else
	      imask[i] = TRUE;

	  for ( r = nremaps-1; r >= 0; r-- )
	    {
	      if ( gridID1 == remaps[r].gridID && nmiss1 == remaps[r].nmiss )
		{
		  if ( memcmp(imask, remaps[r].grid.grid1_mask, remaps[r].grid.grid1_size*sizeof(int)) == 0 )
		    break;
		}	      
	    }

	  if ( cdoVerbose && r >= 0 ) cdoPrint("Using remap %d\n", r);

	  if ( r < 0 )
	    {
	      if ( nremaps < max_remaps )
		{
		  r = nremaps;
		  nremaps++;
		}
	      else
		{
		  r = nremaps - 1;
		}

	      if ( remaps[r].gridID != gridID1 )
		{
		  /*
		    remaps[r].grid.luse_grid1_area = FALSE;
		    remaps[r].grid.luse_grid2_area = FALSE;
		  */
		  remaps[r].grid.restrict_type = RESTRICT_LATITUDE;
		  /* remaps[r].grid.restrict_type = RESTRICT_LATLON; */
		  remaps[r].grid.num_srch_bins = 180;
		  remaps[r].grid.pinit = FALSE;

		  remaps[r].vars.norm_opt = norm_opt;
		  remaps[r].vars.pinit = FALSE;
		  
		  /* initialize grid information for both grids */
		  remapGridInit(map_type, gridID1, gridID2, &remaps[r].grid);
		}

	      remaps[r].gridID = gridID1;
	      remaps[r].nmiss  = nmiss1;

	      if ( gridInqType(gridID1) == GRID_GME )
		{
		  j = 0;
		  for ( i = 0; i < gridsize; i++ )
		    if ( remaps[r].grid.grid1_vgpm[i] ) imask[j++] = imask[i];
		}

	      memcpy(remaps[r].grid.grid1_mask, imask, remaps[r].grid.grid1_size*sizeof(int));
	      /*
	      for ( i = 0; i < gridsize; i++ )
		if ( remaps[r].grid.grid1_mask[i] ) remaps[r].grid.grid1_mask[i] = imask[i];
	      */

	      if ( map_type == MAP_TYPE_CONSERV )
		{
		  memset(remaps[r].grid.grid1_area, 0, remaps[r].grid.grid1_size*sizeof(double));
		  memset(remaps[r].grid.grid1_frac, 0, remaps[r].grid.grid1_size*sizeof(double));
		  memset(remaps[r].grid.grid2_area, 0, remaps[r].grid.grid2_size*sizeof(double));
		}
	      memset(remaps[r].grid.grid2_frac, 0, remaps[r].grid.grid2_size*sizeof(double));

	      /* initialize some remapping variables */
	      remapVarsInit(map_type, &remaps[r].grid, &remaps[r].vars);

	      if      ( map_type == MAP_TYPE_CONSERV  ) remap_conserv(&remaps[r].grid, &remaps[r].vars);
	      else if ( map_type == MAP_TYPE_BILINEAR ) remap_bilin(&remaps[r].grid, &remaps[r].vars);
	      else if ( map_type == MAP_TYPE_BICUBIC  ) remap_bicub(&remaps[r].grid, &remaps[r].vars);
	      else if ( map_type == MAP_TYPE_DISTWGT  ) remap_distwgt(&remaps[r].grid, &remaps[r].vars);
	      else if ( map_type == MAP_TYPE_DISTWGT1 ) remap_distwgt1(&remaps[r].grid, &remaps[r].vars);

	      if ( remaps[r].vars.num_links != remaps[r].vars.max_links )
		resize_remap_vars(&remaps[r].vars, remaps[r].vars.num_links-remaps[r].vars.max_links);

	      sort_add(remaps[r].vars.num_links, remaps[r].vars.num_wts,
		       remaps[r].vars.grid2_add, remaps[r].vars.grid1_add, remaps[r].vars.wts);
	      	      
	      if ( lwrite_remap ) goto WRITE_REMAP;
	    }

	  if ( gridInqType(gridID1) == GRID_GME )
	    {
	      j = 0;
	      for ( i = 0; i < gridsize; i++ )
		if ( remaps[r].grid.grid1_vgpm[i] ) array1[j++] = array1[i];
	    }
	  
	  if ( map_type == MAP_TYPE_BICUBIC )
	    remap_gradients(remaps[r].grid, array1, grad1_lat, grad1_lon, grad1_latlon);

	  if ( operfunc == REMAPCON1 )
	    remap_con1(array2, missval, gridInqSize(gridID2), remaps[r].vars.num_links, remaps[r].vars.wts,
		  remaps[r].vars.grid2_add, remaps[r].vars.grid1_add, array1);
	  else
	    remap(array2, missval, gridInqSize(gridID2), remaps[r].vars.num_links, remaps[r].vars.wts,
		  remaps[r].vars.num_wts, remaps[r].vars.grid2_add, remaps[r].vars.grid1_add,
		  array1, grad1_lat, grad1_lon, grad1_latlon);

	  gridsize2 = gridInqSize(gridID2);

	  /*
	  if ( operfunc == REMAPCON )
	    {
	      double grid2_err;

	      if ( remaps[r].vars.norm_opt == NORM_OPT_NONE )
		{
		  for ( i = 0; i < gridsize2; i++ )
		    {
		      grid2_err = remaps[r].grid.grid2_frac[i]*remaps[r].grid.grid2_area[i];
		      if ( grid2_err != 0 )
			array2[i] = array2[i]/grid2_err;
		      else
			array2[i] = missval;
		    }
		}
	      else if ( remaps[r].vars.norm_opt == NORM_OPT_DESTAREA )
		{
		  for ( i = 0; i < gridsize2; i++ )
		    {
		      if ( remaps[r].grid.grid2_frac[i] != 0 )
			array2[i] = array2[i]/remaps[r].grid.grid2_frac[i];
		      else
			array2[i] = missval;
		    }
		}
	    }
	  */

	  vlistInqVarName(vlistID1, varID, varname);
	  if ( operfunc == REMAPCON )
	    if ( strcmp(varname, "gridbox_area") == 0 )
	      {
		double array1sum = 0;
		double array2sum = 0;

		for ( i = 0; i < gridsize; i++ )
		  array1sum += array1[i];

		for ( i = 0; i < gridsize2; i++ )
		  array2sum += remaps[r].grid.grid2_area[i];

		for ( i = 0; i < gridsize2; i++ )
		  array2[i] = remaps[r].grid.grid2_area[i]/array2sum*array1sum;

		if ( lgridboxinfo )
		  {
		    cdoPrint("%s replaced and scaled to %g", varname, array1sum);
		    lgridboxinfo = FALSE;
		  }
	      }

	  /* calculate some statistics */
	  if ( cdoVerbose ) remap_stat(remaps[r].grid, remaps[r].vars, array1, array2, missval);

	  if ( gridInqType(gridID2) == GRID_GME )
	    {
	      int ni, nd;
	      ni = gridInqGMEni(gridID2);
	      nd = gridInqGMEnd(gridID2);
	      j = remaps[r].grid.grid2_size;

	      for ( i = gridsize2-1; i >=0 ; i-- )
		if ( remaps[r].grid.grid2_vgpm[i] ) array2[i] = array2[--j];

	      gme_grid_restore(array2, ni, nd);
	    }

	  nmiss2 = 0;
	  for ( i = 0; i < gridsize2; i++ )
	    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss2++;

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, array2, nmiss2);
	}
      tsID++;
    }

  streamClose(streamID2);

  WRITE_REMAP:
 
  if ( lwrite_remap ) 
    write_remap_scrip(cdoStreamName(1), map_type, remaps[r].grid, remaps[r].vars);

  streamClose(streamID1);

  if ( imask )  free(imask);
  if ( array2 ) free(array2);
  if ( array1 ) free(array1);

  if ( grad1_latlon ) free(grad1_latlon);
  if ( grad1_lon ) free(grad1_lon);
  if ( grad1_lat ) free(grad1_lat);

  remapVarsFree(&remaps[r].vars);
  remapGridFree(&remaps[r].grid);

  if ( remaps ) free(remaps);

  cdoFinish();

  return (0);
}
