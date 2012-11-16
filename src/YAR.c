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

void testint(field_t *field1, field_t *field2)
{
  int nlonIn, nlatIn;
  int nlonOut, nlatOut;
  int ilat, ilon;
  int gridIDin, gridIDout;
  int i, nmiss;
  double *lonIn, *latIn;
  double *lonOut, *latOut;
  double **fieldIn;
  double **field;
  double *array = NULL;
  double *arrayIn, *arrayOut;
  double missval;
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
  lonIn = (double *) malloc(nlonIn*sizeof(double));
  latIn = (double *) malloc(nlatIn*sizeof(double));
  gridInqXvals(gridIDin, lonIn);
  gridInqYvals(gridIDin, latIn);

  if ( ! (gridInqXvals(gridIDout, NULL) && gridInqYvals(gridIDout, NULL)) )
    cdoAbort("Target grid has no values");

  nlonOut = gridInqXsize(gridIDout);
  nlatOut = gridInqYsize(gridIDout);
  lonOut = (double *) malloc(nlonOut*sizeof(double));
  latOut = (double *) malloc(nlatOut*sizeof(double));
  gridInqXvals(gridIDout, lonOut);
  gridInqYvals(gridIDout, latOut);

#if defined (HAVE_LIBYAC)
  //--------------------------------------------
  // define a grid
  //--------------------------------------------
  /*
  unsigned num_source_cells[2] = {11,11};
  unsigned num_target_cells[2] = {10,10};
  unsigned cyclic[2] = {0,0};

  struct grid source_grid, target_grid;

  init_reg2d_grid(&source_grid, NULL, NULL, num_source_cells, cyclic);
  init_reg2d_grid(&target_grid, NULL, NULL, num_target_cells, cyclic);
  */
  //--------------------------------------------
  // define points
  //--------------------------------------------
  /*
  double source_point_coordinates_x[11] = {0,1,2,3,4,5,6,7,8,9,10};
  double source_point_coordinates_y[11] = {0,1,2,3,4,5,6,7,8,9,10};
  double target_point_coordinates_x[10] = {0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5};
  double target_point_coordinates_y[10] = {0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5};

  struct points source_points, target_points;

  init_points(&source_points, &source_grid, CELL, source_point_coordinates_x, source_point_coordinates_y);
  init_points(&target_points, &target_grid, CELL, target_point_coordinates_x, target_point_coordinates_y);
  */
  //--------------------------------------------
  // initialise interpolation
  //--------------------------------------------
  /*
  struct dep_list tgt_to_src_cell;
  unsigned search_id;
  struct interpolation interpolation;

  search_id = search_init(*get_point_grid(&source_points));

  do_point_search_p(*get_point_grid(&target_points), search_id, &tgt_to_src_cell);

  init_interpolation_g(&interpolation, *get_point_grid(&source_points), *get_point_grid(&target_points),
		       *(struct const_dep_list*)&tgt_to_src_cell,
		       AVERAGE);
  */
  //--------------------------------------------
  // interpolate data
  //--------------------------------------------
  /*
  double source_data[11][11];
  double * p_source_data;
  double target_data[10][10];

  set_source_data((double*)source_data, 0.2, 11, 11);

  p_source_data = &(source_data[0][0]);
  do_interpolation(interpolation, &p_source_data, &(target_data[0][0]));

  for ( int j = 0; j < 10; ++j )
    {
      for ( int i = 0; i < 10; ++i )
	printf("%g ", target_data[j][i]);
      printf("\n");
    }

   free_search(search_id);
   free_interpolation(&interpolation);
   free_dep_list(&tgt_to_src_cell);
   free_points(&source_points);
   free_points(&target_points);
   free_grid(&source_grid);
   free_grid(&target_grid);
  */

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
  /*
  double *source_point_coordinates_x;
  double *source_point_coordinates_y;
  double *target_point_coordinates_x;
  double *target_point_coordinates_y;

  source_point_coordinates_x = (double *) malloc((nlonIn+1)*sizeof(double));
  source_point_coordinates_y = (double *) malloc((nlatIn+1)*sizeof(double));
  target_point_coordinates_x = (double *) malloc((nlonOut+1)*sizeof(double));
  target_point_coordinates_y = (double *) malloc((nlatOut+1)*sizeof(double));

  gen_xbounds(nlonIn, lonIn, source_point_coordinates_x);
  gen_ybounds(nlatIn, latIn, source_point_coordinates_y);
  gen_xbounds(nlonOut, lonOut, target_point_coordinates_x);
  gen_ybounds(nlatOut, latOut, target_point_coordinates_y);

  for ( i = 0; i < nlonIn+1; ++i ) printf("lon %d %g\n", i, source_point_coordinates_x[i]);
  for ( i = 0; i < nlatIn+1; ++i ) printf("lat %d %g\n", i, source_point_coordinates_y[i]);
  */
  /*
  for ( i = 0; i < nlonIn; ++i ) printf("%g,", lonIn[i]); printf("\n");
  for ( i = 0; i < nlatIn; ++i ) printf("%g,", latIn[i]); printf("\n");
  for ( i = 0; i < nlonOut; ++i ) printf("%g,", lonOut[i]); printf("\n");
  for ( i = 0; i < nlatOut; ++i ) printf("%g,", latOut[i]); printf("\n");
  */
  init_points(&source_points, &source_grid, CELL, lonIn, latIn);
  init_points(&target_points, &target_grid, CELL, lonOut, latOut);
  /*
  double source_point_coordinates_x[11] = {0,32,65,98,130,163,196,229,261,294,327};
  double source_point_coordinates_y[11] = {-89,-72,-54,-36,-18,0,18,36,54,72,89};
  double target_point_coordinates_x[10] = {0,36,72,108,144,180,216,252,288,324};
  double target_point_coordinates_y[10] = {-81,-63,-45,-27,-9,9,27,45,63,81};
  */
 /*
  double source_point_coordinates_x[11] = {0,1,2,3,4,5,6,7,8,9,10};
  double source_point_coordinates_y[11] = {0,1,2,3,4,5,6,7,8,9,10};
  double target_point_coordinates_x[10] = {0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5};
  double target_point_coordinates_y[10] = {0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5};
  */
  /*
  init_points(&source_points, &source_grid, CELL, source_point_coordinates_x, source_point_coordinates_y);
  init_points(&target_points, &target_grid, CELL, target_point_coordinates_x, target_point_coordinates_y);
  */
  // printf("num_grid_corners source: %d\n", get_num_grid_corners(source_grid));
  // printf("num_grid_corners target: %d\n", get_num_grid_corners(target_grid));

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

  for ( int i = 0; i < 10; ++i )
    {
      printf("num_deps_per_element %d %d\n", i, tgt_to_src_cell.num_deps_per_element[i]);
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

  cdoOperatorAdd("yarbil",  0, 0, NULL);

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

	  testint(&field1, &field2);

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
