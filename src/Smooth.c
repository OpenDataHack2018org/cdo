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

      Smoothstat       smooth9             running 9-point-average
*/
#include <time.h> // clock()

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"
#include "pmlist.h"

#include "grid_search.h"

void grid_search_nbr(struct gridsearch *gs, int num_neighbors, int *restrict nbr_add, double *restrict nbr_dist, double plon, double plat);

double intlin(double x, double y1, double x1, double y2, double x2);

double smooth_nbr_compute_weights(unsigned num_neighbors, const int *restrict src_grid_mask, int *restrict nbr_mask, const int *restrict nbr_add, double *restrict nbr_dist, double search_radius, double weight0, double weightInf)
{
  // Compute weights based on inverse distance if mask is false, eliminate those points

  double dist_tot = 0.; // sum of neighbor distances (for normalizing)

  for ( unsigned n = 0; n < num_neighbors; ++n )
    {
      nbr_mask[n] = FALSE;
      if ( nbr_add[n] >= 0 && src_grid_mask[nbr_add[n]] )
        {
          nbr_dist[n] = intlin(nbr_dist[n], weight0, 0, weightInf, search_radius);
          dist_tot += nbr_dist[n];
          nbr_mask[n] = TRUE;
        }
    }

  return dist_tot;
}


unsigned smooth_nbr_normalize_weights(unsigned num_neighbors, double dist_tot, const int *restrict nbr_mask, int *restrict nbr_add, double *restrict nbr_dist)
{
  // Normalize weights and store the link

  unsigned nadds = 0;

  for ( unsigned n = 0; n < num_neighbors; ++n )
    {
      if ( nbr_mask[n] )
        {
          nbr_dist[nadds] = nbr_dist[n]/dist_tot;
          nbr_add[nadds]  = nbr_add[n];
          nadds++;
        }
    }

  return nadds;
}


static
void smooth(int gridID, double missval, const double *restrict array1, double *restrict array2, int *nmiss,
            int num_neighbors, double search_radius, double weight0, double weightInf)
{
  *nmiss = 0;
  int gridID0 = gridID;
  unsigned gridsize = gridInqSize(gridID);

  int *mask = (int*) Malloc(gridsize*sizeof(double));
  for ( unsigned i = 0; i < gridsize; ++i )
    if ( DBL_IS_EQUAL(array1[i], missval) )
      mask[i] = 0;
    else
      mask[i] = 1;
  
  double *xvals = (double*) Malloc(gridsize*sizeof(double));
  double *yvals = (double*) Malloc(gridsize*sizeof(double));

  if ( gridInqType(gridID) == GRID_GME ) gridID = gridToUnstructured(gridID, 0);

  if ( gridInqType(gridID) != GRID_UNSTRUCTURED && gridInqType(gridID) != GRID_CURVILINEAR )
    gridID = gridToCurvilinear(gridID, 0);

  gridInqXvals(gridID, xvals);
  gridInqYvals(gridID, yvals);

  /* Convert lat/lon units if required */
  char units[CDI_MAX_NAME];
  gridInqXunits(gridID, units);
  grid_to_radian(units, gridsize, xvals, "grid center lon");
  gridInqYunits(gridID, units);
  grid_to_radian(units, gridsize, yvals, "grid center lat");

  int nbr_mask[num_neighbors];    /* mask at nearest neighbors                   */
  int nbr_add[num_neighbors];     /* source address at nearest neighbors         */
  double nbr_dist[num_neighbors]; /* angular distance four nearest neighbors     */

  clock_t start, finish;
  start = clock();

  struct gridsearch *gs = NULL;

  if ( num_neighbors == 1 )
    gs = gridsearch_create_nn(gridsize, xvals, yvals);
  else
    gs = gridsearch_create(gridsize, xvals, yvals);

  gs->search_radius = search_radius;
  
  finish = clock();

  if ( cdoVerbose ) printf("gridsearch created: %.2f seconds\n", ((double)(finish-start))/CLOCKS_PER_SEC);

  progressInit();

  start = clock();

  double findex = 0;

  /*
#pragma omp parallel for default(none) shared(findex, array1, array2, xvals, yvals, gs, gridsize, num_neighbors) \
                                      private(nbr_mask, nbr_add, nbr_dist)
  */
  for ( unsigned i = 0; i < gridsize; ++i )
    {
      /*
#if defined(_OPENMP)
#include "pragma_omp_atomic_update.h"
#endif
      */
      findex++;
      if ( cdo_omp_get_thread_num() == 0 ) progressStatus(0, 1, findex/gridsize);

      grid_search_nbr(gs, num_neighbors, nbr_add, nbr_dist, xvals[i], yvals[i]);

      /* Compute weights based on inverse distance if mask is false, eliminate those points */
      double dist_tot = smooth_nbr_compute_weights(num_neighbors, mask, nbr_mask, nbr_add, nbr_dist,
                                                   search_radius, weight0, weightInf);

      /* Normalize weights and store the link */
      unsigned nadds = smooth_nbr_normalize_weights(num_neighbors, dist_tot, nbr_mask, nbr_add, nbr_dist);
      if ( nadds )
        {
          /*
          printf("n %u %d nadds %u dis %g\n", i, nbr_add[0], nadds, nbr_dist[0]);
          for ( unsigned n = 0; n < nadds; ++n )
            printf("   n %u add %d dis %g\n", n, nbr_add[n], nbr_dist[n]);
          */
          double result = 0;
          for ( unsigned n = 0; n < nadds; ++n ) result += array1[nbr_add[n]]*nbr_dist[n];
          array2[i] = result;
        }
      else
        {
          (*nmiss)++;
          array2[i] = missval;
        }
    }

  finish = clock();

  if ( cdoVerbose ) printf("gridsearch nearest: %.2f seconds\n", ((double)(finish-start))/CLOCKS_PER_SEC);

  if ( gs ) gridsearch_delete(gs);

  if ( gridID0 != gridID ) gridDestroy(gridID);

  Free(mask);
  Free(xvals);
  Free(yvals);
}

static inline
void smooth9_sum(size_t ij, short *mask, double sfac, const double *restrict array, double *avg, double *divavg)
{
  if ( mask[ij] ) { *avg += sfac*array[ij]; *divavg += sfac; }
}

static
void smooth9(int gridID, double missval, const double *restrict array1, double *restrict array2, int *nmiss)
{
  double avg,divavg;
  size_t i, ij , j;
  size_t gridsize = gridInqSize(gridID);
  size_t nlon = gridInqXsize(gridID);	 
  size_t nlat = gridInqYsize(gridID);
  int grid_is_cyclic = gridIsCircular(gridID);

  short *mask = (short*) Malloc(gridsize*sizeof(short));

  for ( size_t i = 0; i < gridsize; i++) 
    {		
      if ( DBL_IS_EQUAL(missval, array1[i]) ) mask[i] = 0;
      else mask[i] = 1;
    }
 
  *nmiss = 0;
  for ( i = 0; i < nlat; i++ )
    {
      for ( j = 0; j < nlon; j++ )
        {		      
          avg = 0; divavg = 0; 	  		     			

          if ( (i == 0) || (j == 0) || (i == (nlat-1)) || (j == (nlon-1)) )
            {
              ij = j+nlon*i;
              if ( mask[ij] )
                {
                  avg += array1[ij];  divavg+= 1;					     		       
                  /* upper left corner */
                  if ( (i != 0) && (j != 0) ) 
                    smooth9_sum(((i-1)*nlon)+j-1, mask, 0.3, array1, &avg, &divavg);
                  else if ( i != 0 && grid_is_cyclic ) 
                    smooth9_sum((i-1)*nlon+j-1+nlon, mask, 0.3, array1, &avg, &divavg);
			      
                  /* upper cell */
                  if ( i != 0 ) 
                    smooth9_sum(((i-1)*nlon)+j, mask, 0.5, array1, &avg, &divavg);
                  
                  /* upper right corner */
                  if ( (i != 0) && (j != (nlon-1)) ) 
                    smooth9_sum(((i-1)*nlon)+j+1, mask, 0.3, array1, &avg, &divavg);
                  else if ( (i !=0 ) && grid_is_cyclic )
                    smooth9_sum((i-1)*nlon+j+1-nlon, mask, 0.3, array1, &avg, &divavg);
                  
                  /* left cell */
                  if  ( j != 0 ) 
                    smooth9_sum(((i)*nlon)+j-1, mask, 0.5, array1, &avg, &divavg);
                  else if ( grid_is_cyclic )
                    smooth9_sum(i*nlon-1+nlon, mask, 0.5, array1, &avg, &divavg);
                  
                  /* right cell */
                  if ( j!=(nlon-1) ) 
                    smooth9_sum((i*nlon)+j+1, mask, 0.5, array1, &avg, &divavg);
                  else if ( grid_is_cyclic )
                    smooth9_sum(i*nlon+j+1-nlon, mask, 0.5, array1, &avg, &divavg);
                  
                  /* lower left corner */
                  if ( mask[ij] &&  ( (i!=(nlat-1))&& (j!=0) ) )
                    smooth9_sum(((i+1)*nlon+j-1), mask, 0.3, array1, &avg, &divavg);
                  else if ( (i != (nlat-1)) && grid_is_cyclic ) 
                    smooth9_sum((i+1)*nlon-1+nlon, mask, 0.3, array1, &avg, &divavg);
                  
                  /* lower cell */
                  if  ( i != (nlat-1) ) 
                    smooth9_sum(((i+1)*nlon)+j, mask, 0.5, array1, &avg, &divavg);
                  
                  /* lower right corner */
                  if ( (i != (nlat-1)) && (j != (nlon-1)) )
                    smooth9_sum(((i+1)*nlon)+j+1, mask, 0.3, array1, &avg, &divavg);
                  else if ( (i != (nlat-1)) && grid_is_cyclic )
                    smooth9_sum(((i+1)*nlon)+j+1-nlon, mask, 0.3, array1, &avg, &divavg);
                }
            }
          else if ( mask[j+nlon*i] )
            {			 
              avg += array1[j+nlon*i]; divavg+= 1;
			    
              smooth9_sum(((i-1)*nlon)+j-1, mask, 0.3, array1, &avg, &divavg);
              smooth9_sum(((i-1)*nlon)+j,   mask, 0.5, array1, &avg, &divavg);
              smooth9_sum(((i-1)*nlon)+j+1, mask, 0.3, array1, &avg, &divavg);
              smooth9_sum(((i)*nlon)+j-1,   mask, 0.5, array1, &avg, &divavg);
              smooth9_sum((i*nlon)+j+1,     mask, 0.5, array1, &avg, &divavg);
              smooth9_sum(((i+1)*nlon+j-1), mask, 0.3, array1, &avg, &divavg);
              smooth9_sum(((i+1)*nlon)+j,   mask, 0.5, array1, &avg, &divavg);
              smooth9_sum(((i+1)*nlon)+j+1, mask, 0.3, array1, &avg, &divavg);
            }

          if ( fabs(divavg) > 0 )
            {
              array2[i*nlon+j] = avg/divavg;			
            }
          else 
            {
              array2[i*nlon+j] = missval;					
              (*nmiss)++;
            }
        }			    	     
    }

  Free(mask);
}


void *Smooth(void *argument)
{
  int gridID;
  int nrecs;
  int varID, levelID;
  int nmiss;
  int gridtype;
  int xnsmooth = 1;
  int xmax_points = 5;
  double xsearch_radius = 90, xweight0 = 0.25, xweightInf = 0.25;

  cdoInitialize(argument);

  int SMOOTHP = cdoOperatorAdd("smoothpoint",  0,   0, NULL);
  int SMOOTH9 = cdoOperatorAdd("smooth9",      0,   0, NULL);
 
  int operatorID = cdoOperatorID();

  if ( operatorID == SMOOTHP )
    {      
      int pargc = operatorArgc();

      if ( pargc )
        {
          char **pargv = operatorArgv();
          pml_t *pml = pml_create("SMOOTH");

          PML_ADD_INT(pml, nsmooth,           1, "Number of smooth iterations");
          PML_ADD_INT(pml, max_points,        1, "Maximum number of points");
          PML_ADD_FLT(pml, search_radius,     1, "Search radius");
          PML_ADD_FLT(pml, weight0,           1, "weight0");
          PML_ADD_FLT(pml, weightInf,         1, "weightInf");
      
          pml_read(pml, pargc, pargv);
          if ( cdoVerbose ) pml_print(pml);
      
          if ( PML_NOCC(pml, nsmooth) )       xnsmooth = par_nsmooth[0];
          if ( PML_NOCC(pml, max_points) )    xmax_points = par_max_points[0];
          if ( PML_NOCC(pml, search_radius) ) xsearch_radius = par_search_radius[0];
          if ( PML_NOCC(pml, weight0) )       xweight0 = par_weight0[0];
          if ( PML_NOCC(pml, weightInf) )     xweightInf = par_weightInf[0];

          UNUSED(nsmooth);
          UNUSED(max_points);
          UNUSED(search_radius);
          UNUSED(weight0);
          UNUSED(weightInf);

          pml_destroy(pml);
        }
      
      if ( cdoVerbose )
        cdoPrint("nsmooth = %d, max_points = %d, search_radius = %g, weight0 = %g, weightInf = %g",
                 xnsmooth, xmax_points, xsearch_radius, xweight0, xweightInf);
      
    }

  xsearch_radius *= DEG2RAD;

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int nvars = vlistNvars(vlistID1);
  int *varIDs = (int*) Malloc(nvars*sizeof(int)); 

  for ( varID = 0; varID < nvars; ++varID )
    {
      gridID = vlistInqVarGrid(vlistID1, varID);
      gridtype = gridInqType(gridID);
      if ( gridtype == GRID_GAUSSIAN ||
           gridtype == GRID_LONLAT   ||
           gridtype == GRID_CURVILINEAR )
	{
	  varIDs[varID] = 1;
	}
      else
	{
	  varIDs[varID] = 0;
	  cdoWarning("Unsupported grid for varID %d", varID);
	}
    }

  size_t gridsize = vlistGridsizeMax(vlistID1);
  double *array1 = (double*) Malloc(gridsize*sizeof(double));
  double *array2 = (double*) Malloc(gridsize*sizeof(double));
 
  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      streamDefTimestep(streamID2, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss);
	
	  if ( varIDs[varID] )
	    {	    
	      double missval = vlistInqVarMissval(vlistID1, varID);
	      gridID = vlistInqVarGrid(vlistID1, varID);

              for ( int i = 0; i < xnsmooth; ++i )
                {
                  if ( operatorID == SMOOTHP )
                    smooth(gridID, missval, array1, array2, &nmiss, xmax_points, xsearch_radius, xweight0, xweightInf);
                  else if ( operatorID == SMOOTH9 )
                    smooth9(gridID, missval, array1, array2, &nmiss);
                  
                  memcpy(array1, array2, gridsize*sizeof(double));
                }
          
	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, array2, nmiss);		
	    }     	   
	  else 
	    {
	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, array1, nmiss);
	    }
	}

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  Free(varIDs);
  if ( array2 ) Free(array2);
  if ( array1 ) Free(array1);

  cdoFinish();

  return 0;
}
