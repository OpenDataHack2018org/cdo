/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2012 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "interpol.h"

#if defined (HAVE_LIBYAC)
#include "points.h"
#include "grid_reg2d.h"
#include "event.h"
#include "search.h"
#endif

static
void gen_xbounds(int nx, double *xvals, double *xbounds)
{
  int i;

  for ( i = 0; i < nx-1; i++ )
    {
      xbounds[i+1]   = 0.5*(xvals[i] + xvals[i+1]);
    }

  xbounds[0]  = 2*xvals[0] - xbounds[1];
  xbounds[nx] = 2*xvals[nx-1] - xbounds[nx-1];
}

static
void gen_ybounds(int ny, double *yvals, double *ybounds)
{
  int i;

  for ( i = 0; i < ny-1; i++ )
    {
      ybounds[i+1]   = 0.5*(yvals[i] + yvals[i+1]);
    }

  ybounds[0]  = 2*yvals[0] - ybounds[1];
  ybounds[ny] = 2*yvals[ny-1] - ybounds[ny-1];

  if ( yvals[0] > yvals[ny-1] )
    {
      if ( ybounds[0]  >  88 ) ybounds[0]  =  90;
      if ( ybounds[ny] < -88 ) ybounds[ny] = -90;
    }
  else
    {
      if ( ybounds[0]  < -88 ) ybounds[0]  = -90;
      if ( ybounds[ny] >  88 ) ybounds[ny] =  90;
    }
}


void set_source_data(double * source_data, double init_value,
                     unsigned size_x, unsigned size_y) {

   for (unsigned i = 0; i < size_x; ++i)
      for (unsigned j = 0; j < size_y; ++j)
         source_data[i + j * size_x] = init_value;
}


int grid_search( int *restrict src_add, double *restrict src_lats, 
		double *restrict src_lons, double plat, double plon, const int *restrict src_grid_dims)
{
  return 0;
}

void testint_p(field_t *field1, field_t *field2)
{
  int nlonIn, nlatIn;
  int nlonOut, nlatOut;
  int ilat, ilon;
  int gridIDin, gridIDout;
  int i, nmiss;
  int gridsize1, gridsize2;
  double *lonIn, *latIn;
  double *lonOut, *latOut;
  double **fieldIn;
  double **field;
  double *array = NULL;
  double *arrayIn, *arrayOut;
  double missval;
  int testit = 1;
  /* static int index = 0; */

  gridIDin  = field1->grid;
  gridIDout = field2->grid;
  arrayIn   = field1->ptr;
  arrayOut  = field2->ptr;
  missval   = field1->missval;

  if ( ! (gridInqXvals(gridIDin, NULL) && gridInqYvals(gridIDin, NULL)) )
    cdoAbort("Source grid has no values");

  nlonIn = gridInqXsize(gridIDin);
  nlatIn = gridInqYsize(gridIDin);
  gridsize1 = gridInqSize(gridIDin);
  lonIn = (double *) malloc(nlonIn*sizeof(double));
  latIn = (double *) malloc(nlatIn*sizeof(double));
  gridInqXvals(gridIDin, lonIn);
  gridInqYvals(gridIDin, latIn);

  if ( ! (gridInqXvals(gridIDout, NULL) && gridInqYvals(gridIDout, NULL)) )
    cdoAbort("Target grid has no values");

  nlonOut = gridInqXsize(gridIDout);
  nlatOut = gridInqYsize(gridIDout);
  gridsize2 = gridInqSize(gridIDout);
  lonOut = (double *) malloc(nlonOut*sizeof(double));
  latOut = (double *) malloc(nlatOut*sizeof(double));
  gridInqXvals(gridIDout, lonOut);
  gridInqYvals(gridIDout, latOut);

#if defined (HAVE_LIBYAC)

  //--------------------------------------------
  // define a grid
  //--------------------------------------------
  unsigned num_source_cells[2];
  unsigned num_target_cells[2];
  num_source_cells[0] = nlonIn;
  num_source_cells[1] = nlatIn;
  num_target_cells[0] = nlonOut;
  num_target_cells[1] = nlatOut;

  unsigned cyclic[2] = {0,0};
  struct grid source_grid, target_grid;

  init_reg2d_grid(&source_grid, NULL, NULL, num_source_cells, cyclic);
  init_reg2d_grid(&target_grid, NULL, NULL, num_target_cells, cyclic);

  struct points source_points, target_points;

  //--------------------------------------------
  // define points
  //--------------------------------------------
  init_points(&source_points, &source_grid, CELL, lonIn, latIn);
  init_points(&target_points, &target_grid, CELL, lonOut, latOut);

  //--------------------------------------------
  // initialise interpolation
  //--------------------------------------------

  struct dep_list tgt_to_src_cell;
  unsigned search_id;
  //struct interpolation interpolation;

  // seg fault:  printf("src num_grid_corners %d\n", get_num_grid_corners(*get_point_grid(&source_grid)));
  // seg fault:  printf("tgt num_grid_corners %d\n", get_num_grid_corners(*get_point_grid(&target_grid)));
  printf("src num_grid_corners %d\n", get_num_grid_corners(*get_point_grid(&source_points)));
  printf("tgt num_grid_corners %d\n", get_num_grid_corners(*get_point_grid(&target_points)));

  search_id = search_init(get_point_grid(&source_points));
 
  do_point_search_p(*get_point_grid(&target_points), search_id, &tgt_to_src_cell);

  printf("total_num_dependencies: %d\n", get_total_num_dependencies(tgt_to_src_cell));

  unsigned const * curr_src_corners;
  for ( int i = 0; i < 10; ++i )
    {
      printf("num_deps_per_element %d %d\n", i, tgt_to_src_cell.num_deps_per_element[i]);
      curr_src_corners = get_dependencies_of_element(tgt_to_src_cell,i);
      for ( int k = 0; k < tgt_to_src_cell.num_deps_per_element[i]; ++k )
	printf("  curr_src_corners: %d %d\n", k, curr_src_corners[k]);
    }

  for ( int i = 0; i < 10; ++i )
    {
      double lon = lonOut[i];
      double lat = latOut[i];
      printf("num_deps_per_element %d %d\n", i, tgt_to_src_cell.num_deps_per_element[i]);
      curr_src_corners = get_dependencies_of_element(tgt_to_src_cell,i);
      for ( int k = 0; k < tgt_to_src_cell.num_deps_per_element[i]; ++k )
	printf("  curr_src_corners: %d %d\n", k, curr_src_corners[k]);
    }

  for ( int i = 0; i < gridsize2; ++i )
    {
      double lon = lonOut[i];
      double lat = latOut[i];
      curr_src_corners = get_dependencies_of_element(tgt_to_src_cell,i);
      for ( int k = 0; k < tgt_to_src_cell.num_deps_per_element[i]; ++k )
	{
	}
    }

  /*
  for ( int j = 0; j < 10; ++j )
    {
      for ( int i = 0; i < 10; ++i )
	printf("%g ", arrayOut[j*10+i]);
      printf("\n");
    }
  */
  // free (search_result);
  // if (array) free(array);
  //free(lonIn);
  //free(latIn);
  //free(lonOut);
  //free(latOut);
  //free(fieldIn);
#endif
}


void testint_c(field_t *field1, field_t *field2)
{
  int nlonIn, nlatIn;
  int nlonOut, nlatOut;
  int ilat, ilon;
  int gridIDin, gridIDout;
  int i, nmiss;
  int gridsize1, gridsize2;
  double *lonIn, *latIn;
  double *lonOut, *latOut;
  double **fieldIn;
  double **field;
  double *array = NULL;
  double *arrayIn, *arrayOut;
  double missval;
  double dxIn, dxOut;
  /* static int index = 0; */

  gridIDin  = field1->grid;
  gridIDout = field2->grid;
  arrayIn   = field1->ptr;
  arrayOut  = field2->ptr;
  missval   = field1->missval;

  if ( ! (gridInqXvals(gridIDin, NULL) && gridInqYvals(gridIDin, NULL)) )
    cdoAbort("Source grid has no values");

  nlonIn = gridInqXsize(gridIDin);
  nlatIn = gridInqYsize(gridIDin);
  gridsize1 = gridInqSize(gridIDin);
  lonIn = (double *) malloc(nlonIn*sizeof(double));
  latIn = (double *) malloc(nlatIn*sizeof(double));
  gridInqXvals(gridIDin, lonIn);
  gridInqYvals(gridIDin, latIn);
  dxIn = lonIn[1] - lonIn[0];

  if ( ! (gridInqXvals(gridIDout, NULL) && gridInqYvals(gridIDout, NULL)) )
    cdoAbort("Target grid has no values");

  nlonOut = gridInqXsize(gridIDout);
  nlatOut = gridInqYsize(gridIDout);
  gridsize2 = gridInqSize(gridIDout);
  lonOut = (double *) malloc(nlonOut*sizeof(double));
  latOut = (double *) malloc(nlatOut*sizeof(double));
  gridInqXvals(gridIDout, lonOut);
  gridInqYvals(gridIDout, latOut);
  dxOut = lonOut[1] - lonOut[0];

  printf("dxIn: %g   dxOut: %g\n", dxIn, dxOut);
#if defined (HAVE_LIBYAC)

  //--------------------------------------------
  // define a grid
  //--------------------------------------------
  unsigned num_source_cells[2];
  unsigned num_target_cells[2];
  num_source_cells[0] = nlonIn;
  num_source_cells[1] = nlatIn;
  num_target_cells[0] = nlonOut;
  num_target_cells[1] = nlatOut;

  unsigned cyclic[2] = {0,0};
  struct grid source_grid, target_grid;

  init_reg2d_grid(&source_grid, NULL, NULL, num_source_cells, cyclic);
  init_reg2d_grid(&target_grid, NULL, NULL, num_target_cells, cyclic);

  struct points source_points, target_points;

  //--------------------------------------------
  // define points
  //--------------------------------------------
  init_points(&source_points, &source_grid, CELL, lonIn, latIn);
  init_points(&target_points, &target_grid, CELL, lonOut, latOut);

  //--------------------------------------------
  // initialise interpolation
  //--------------------------------------------

  struct dep_list tgt_to_src_cell;
  unsigned search_id;
  //struct interpolation interpolation;

  // seg fault:  printf("src num_grid_corners %d\n", get_num_grid_corners(*get_point_grid(&source_grid)));
  // seg fault:  printf("tgt num_grid_corners %d\n", get_num_grid_corners(*get_point_grid(&target_grid)));
  printf("src num_grid_corners %d\n", get_num_grid_corners(*get_point_grid(&source_points)));
  printf("tgt num_grid_corners %d\n", get_num_grid_corners(*get_point_grid(&target_points)));

  search_id = search_init(get_point_grid(&source_points));
 
  do_cell_search(*get_point_grid(&target_points), search_id, &tgt_to_src_cell);

  printf("total_num_dependencies: %d\n", get_total_num_dependencies(tgt_to_src_cell));


  enum edge_type quad_type[] = {GREAT_CIRCLE, GREAT_CIRCLE, GREAT_CIRCLE, GREAT_CIRCLE};

  double *weight;
  weight = (double *) malloc(gridsize1*sizeof(double));

  struct grid_cell *SourceCell;
  SourceCell = malloc (gridsize1  * sizeof(*SourceCell) );

  for ( int n = 0; n <  gridsize1; n++ ) {
    SourceCell[n].num_corners   = 4;
    SourceCell[n].edge_type     = quad_type;
    SourceCell[n].coordinates_x = malloc ( 4 * sizeof(SourceCell[n].coordinates_x[0]) );
    SourceCell[n].coordinates_y = malloc ( 4 * sizeof(SourceCell[n].coordinates_y[0]) );
  }

  struct grid_cell  TargetCell;

  TargetCell.num_corners   = 4;
  TargetCell.edge_type     = quad_type;

  TargetCell.coordinates_x = malloc ( 4 * sizeof(*TargetCell.coordinates_x) );
  TargetCell.coordinates_y = malloc ( 4 * sizeof(*TargetCell.coordinates_y) );

  unsigned const * curr_deps;

  for ( int i = 0; i < gridsize2; ++i )
    {
      int index2 = i;
      int ilat2 = index2/nlonOut;
      int ilon2 = index2 - ilat2*nlonOut;

      TargetCell.coordinates_x[0] =  lonOut[ilon2]-dxOut/2;
      TargetCell.coordinates_y[0] =  latOut[ilat2]-dxOut/2;
      TargetCell.coordinates_x[1] =  lonOut[ilon2]+dxOut/2;
      TargetCell.coordinates_y[1] =  latOut[ilat2]-dxOut/2;
      TargetCell.coordinates_x[2] =  lonOut[ilon2]+dxOut/2;
      TargetCell.coordinates_y[2] =  latOut[ilat2]+dxOut/2;
      TargetCell.coordinates_x[3] =  lonOut[ilon2]-dxOut/2;
      TargetCell.coordinates_y[3] =  latOut[ilat2]-dxOut/2;

      if ( i < 10 )
	printf("num_deps_per_element %d %d\n", i, tgt_to_src_cell.num_deps_per_element[i]);
      int num_deps = tgt_to_src_cell.num_deps_per_element[i];
      if ( num_deps > 0 ) curr_deps = get_dependencies_of_element(tgt_to_src_cell, i);
      for ( int k = 0; k < num_deps; ++k )
	{
	  int index1 = curr_deps[k];
	  int ilat1 = index1/nlonIn;
	  int ilon1 = index1 - ilat1*nlonIn;
	  if ( i < 10 )
	    printf("  dep: %d %d %d %d %d %d\n", k, nlonOut, nlatOut, index1, ilon1, ilat1);
	
	  SourceCell[k].coordinates_x[0] =  lonIn[ilon1]-dxIn/2;
	  SourceCell[k].coordinates_y[0] =  latIn[ilat1]-dxIn/2;
	  SourceCell[k].coordinates_x[1] =  lonIn[ilon1]+dxIn/2;
	  SourceCell[k].coordinates_y[1] =  latIn[ilat1]-dxIn/2;
	  SourceCell[k].coordinates_x[2] =  lonIn[ilon1]+dxIn/2;
	  SourceCell[k].coordinates_y[2] =  latIn[ilat1]+dxIn/2;
	  SourceCell[k].coordinates_x[3] =  lonIn[ilon1]-dxIn/2;
	  SourceCell[k].coordinates_y[3] =  latIn[ilat1]-dxIn/2;
	}
      
      polygon_partial_weights (num_deps, SourceCell, TargetCell, weight );
      for ( int k = 0; k < num_deps; ++k )
	{
	  int index1 = curr_deps[k];
	  int ilat1 = index1/nlonIn;
	  int ilon1 = index1 - ilat1*nlonIn;
	  // if ( i < 10 )
	    printf("  dep: %d %d %d %d %d %d  %g\n", k, nlonOut, nlatOut, index1, ilon1, ilat1, weight[k]);
	}
      // correct_weights ( nSourceCells, weight );
    }

  /* xxxxxx1
  init_interpolation_g(&interpolation, *get_point_grid(&source_points), *get_point_grid(&target_points),
		       *(struct const_dep_list*)&tgt_to_src_cell,
		       AVERAGE);
  */
  //--------------------------------------------
  // interpolate data
  //--------------------------------------------
  /* xxxxx1
  do_interpolation(interpolation, &arrayIn, arrayOut);
  */
  /*
  for ( int j = 0; j < 10; ++j )
    {
      for ( int i = 0; i < 10; ++i )
	printf("%g ", arrayOut[j*10+i]);
      printf("\n");
    }
  */
  // free (search_result);
  // if (array) free(array);
  //free(lonIn);
  //free(latIn);
  //free(lonOut);
  //free(latOut);
  //free(fieldIn);
#endif
}


void *YAR(void *argument)
{
  int YARBIL, YARCON;
  int operatorID;
  int streamID1, streamID2;
  int nrecs, ngrids;
  int index;
  int tsID, recID, varID, levelID;
  int gridsize;
  int vlistID1, vlistID2;
  int gridID1 = -1, gridID2 = -1;
  int nmiss;
  int xinc = 0, yinc = 0;
  double missval;
  double slon, slat;
  double *array1 = NULL, *array2 = NULL;
  field_t field1, field2;
  int taxisID1, taxisID2;

  cdoInitialize(argument);

  YARBIL = cdoOperatorAdd("yarbil",  0, 0, NULL);
  YARCON = cdoOperatorAdd("yarcon",  0, 0, NULL);

  operatorID = cdoOperatorID();

  operatorInputArg("grid description file or name");
  gridID2 = cdoDefineGrid(operatorArgv()[0]);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  ngrids = vlistNgrids(vlistID1);
  for ( index = 0; index < ngrids; index++ )
    {
      gridID1 = vlistGrid(vlistID1, index);

      if ( gridInqType(gridID1) != GRID_LONLAT && gridInqType(gridID1) != GRID_GAUSSIAN )
	cdoAbort("Interpolation of %s data unsupported!", gridNamePtr(gridInqType(gridID1)) );

      if ( gridIsRotated(gridID1) )
	cdoAbort("Rotated grids not supported!");

      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  gridsize = vlistGridsizeMax(vlistID1);
  array1   = (double *) malloc(gridsize*sizeof(double));

  gridsize = gridInqSize(gridID2);
  array2   = (double *) malloc(gridsize*sizeof(double));

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);
	       
      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss);

	  gridID1 = vlistInqVarGrid(vlistID1, varID);
	  missval = vlistInqVarMissval(vlistID1, varID);

	  field1.grid    = gridID1;
	  field1.nmiss   = nmiss;
	  field1.missval = missval;
	  field1.ptr     = array1;
	  field2.grid    = gridID2;
	  field2.ptr     = array2;
	  field2.nmiss   = 0;

	  if ( operatorID == YARBIL )
	    testint_p(&field1, &field2);
	  else if ( operatorID == YARCON )
	    testint_c(&field1, &field2);
	  else
	    cdoAbort("Not implemented!");

	  nmiss = field2.nmiss;

	  streamDefRecord(streamID2, varID, levelID);
	  streamWriteRecord(streamID2, array2, nmiss);
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  if ( array2 ) free(array2);
  if ( array1 ) free(array1);

  cdoFinish();

  return (0);
}
