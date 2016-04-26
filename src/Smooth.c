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
#include "constants.h" // planet radius
#include "pmlist.h"

#include "grid_search.h"

enum {FORM_LINEAR};
static const char *Form[] = {"linear"};

typedef struct {
  int maxpoints;
  int form;
  double radius;
  double weight0;
  double weightR;
} smoothpoint_t;


double intlin(double x, double y1, double x1, double y2, double x2);

double smooth_nbr_compute_weights(unsigned num_neighbors, const int *restrict src_grid_mask, int *restrict nbr_mask, const int *restrict nbr_add, double *restrict nbr_dist, double search_radius, double weight0, double weightR)
{
  // Compute weights based on inverse distance if mask is false, eliminate those points
  double dist_tot = 0.; // sum of neighbor distances (for normalizing)

  for ( unsigned n = 0; n < num_neighbors; ++n )
    {
      nbr_mask[n] = FALSE;
      if ( nbr_add[n] >= 0 && src_grid_mask[nbr_add[n]] )
        {
          nbr_dist[n] = intlin(nbr_dist[n], weight0, 0, weightR, search_radius);
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
void smooth(int gridID, double missval, const double *restrict array1, double *restrict array2, int *nmiss, smoothpoint_t spoint)
{
  *nmiss = 0;
  int gridID0 = gridID;
  unsigned gridsize = gridInqSize(gridID);
  unsigned num_neighbors = spoint.maxpoints;
  if ( num_neighbors > gridsize ) num_neighbors = gridsize;

  int *mask = (int*) Malloc(gridsize*sizeof(int));
  for ( unsigned i = 0; i < gridsize; ++i )
    mask[i] = DBL_IS_EQUAL(array1[i], missval) ? 0 : 1;
  
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
  
  int *nbr_mask = (int*) Malloc(num_neighbors*sizeof(int));          /* mask at nearest neighbors                */

  struct gsknn *knn = gridsearch_knn_new(num_neighbors);

  clock_t start, finish;

  start = clock();

  struct gridsearch *gs = NULL;

  if ( num_neighbors == 1 )
    gs = gridsearch_create_nn(gridsize, xvals, yvals);
  else
    gs = gridsearch_create(gridsize, xvals, yvals);

  gs->search_radius = spoint.radius;

  finish = clock();

  if ( cdoVerbose ) printf("gridsearch created: %.2f seconds\n", ((double)(finish-start))/CLOCKS_PER_SEC);

  if ( cdoVerbose ) progressInit();

  start = clock();

  double findex = 0;

  /*
#pragma omp parallel for default(none) shared(findex, array1, array2, xvals, yvals, gs, gridsize) \
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
      if ( cdoVerbose && cdo_omp_get_thread_num() == 0 ) progressStatus(0, 1, findex/gridsize);

      int *nbr_add = knn->add;
      double *nbr_dist = knn->dist;
      
      unsigned nadds = gridsearch_knn(gs, knn, xvals[i], yvals[i]);

      /* Compute weights based on inverse distance if mask is false, eliminate those points */
      double dist_tot = smooth_nbr_compute_weights(nadds, mask, nbr_mask, nbr_add, nbr_dist,
                                                   spoint.radius, spoint.weight0, spoint.weightR);

      /* Normalize weights and store the link */
      nadds = smooth_nbr_normalize_weights(nadds, dist_tot, nbr_mask, nbr_add, nbr_dist);
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

  gridsearch_knn_delete(knn);

  if ( gridID0 != gridID ) gridDestroy(gridID);

  Free(nbr_mask);
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
  size_t gridsize = gridInqSize(gridID);
  size_t nlon = gridInqXsize(gridID);	 
  size_t nlat = gridInqYsize(gridID);
  int grid_is_cyclic = gridIsCircular(gridID);

  short *mask = (short*) Malloc(gridsize*sizeof(short));

  for ( size_t i = 0; i < gridsize; ++i ) 
    {		
      if ( DBL_IS_EQUAL(missval, array1[i]) ) mask[i] = 0;
      else mask[i] = 1;
    }
 
  *nmiss = 0;
  for ( size_t i = 0; i < nlat; i++ )
    {
      for ( size_t j = 0; j < nlon; j++ )
        {		      
          double avg = 0;
          double divavg = 0; 	  		     			

          if ( (i == 0) || (j == 0) || (i == (nlat-1)) || (j == (nlon-1)) )
            {
              size_t ij = j+nlon*i;
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

static
double convert_radius(const char *string)
{
  char *endptr = NULL;
  double radius = strtod(string, &endptr);

  if ( *endptr != 0 )
    {
      if ( strcmp(endptr, "km") == 0 )
        radius = 360*((radius*1000)/(2*PlanetRadius*M_PI));
      else if ( strncmp(endptr, "m", 1) == 0 )
        radius = 360*((radius)/(2*PlanetRadius*M_PI));
      else if ( strncmp(endptr, "deg", 3) == 0 )
        ;
      else if ( strncmp(endptr, "rad", 3) == 0 )
        radius *= RAD2DEG;
      else
        cdoAbort("Float parameter >%s< contains invalid character at position %d!",
                 string, (int)(endptr-string+1));
    }

  if ( radius > 180. ) radius = 180.;

  return radius;
}

static
int convert_form(const char *formstr)
{
  int form = FORM_LINEAR;

  if ( strcmp(formstr, "linear") == 0 ) form = FORM_LINEAR;
  else cdoAbort("form=%s unsupported!", formstr);

  return form;
}


void *Smooth(void *argument)
{
  int nrecs;
  int varID, levelID;
  int nmiss;
  int xnsmooth = 1;
  smoothpoint_t spoint;
  spoint.maxpoints = INT_MAX;
  spoint.radius    = 1;
  spoint.form      = FORM_LINEAR;
  spoint.weight0   = 0.25;
  spoint.weightR   = 0.25;

  cdoInitialize(argument);

  int SMOOTH  = cdoOperatorAdd("smooth",   0,   0, NULL);
  int SMOOTH9 = cdoOperatorAdd("smooth9",  0,   0, NULL);
 
  int operatorID = cdoOperatorID();

  if ( operatorID == SMOOTH )
    {
      int pargc = operatorArgc();

      if ( pargc )
        {
          char **pargv = operatorArgv();
          pml_t *pml = pml_create("SMOOTH");

          PML_ADD_INT(pml, nsmooth,   1, "Number of times to smooth");
          PML_ADD_INT(pml, maxpoints, 1, "Maximum number of points");
          PML_ADD_FLT(pml, weight0,   1, "Weight at distance 0");
          PML_ADD_FLT(pml, weightR,   1, "Weight at the search radius");
          PML_ADD_WORD(pml, radius,   1, "Search radius");
          PML_ADD_WORD(pml, form,     1, "Form of the curve (linear, exponential, gauss)");
      
          pml_read(pml, pargc, pargv);
          if ( cdoVerbose ) pml_print(pml);
      
          if ( PML_NOCC(pml, nsmooth) )   xnsmooth         = par_nsmooth[0];
          if ( PML_NOCC(pml, maxpoints) ) spoint.maxpoints = par_maxpoints[0];
          if ( PML_NOCC(pml, weight0) )   spoint.weight0   = par_weight0[0];
          if ( PML_NOCC(pml, weightR) )   spoint.weightR   = par_weightR[0];
          if ( PML_NOCC(pml, radius) )    spoint.radius    = convert_radius(par_radius[0]);
          if ( PML_NOCC(pml, form) )      spoint.form      = convert_form(par_form[0]);

          UNUSED(nsmooth);
          UNUSED(maxpoints);
          UNUSED(radius);
          UNUSED(form);
          UNUSED(weight0);
          UNUSED(weightR);

          pml_destroy(pml);
        }
      
      if ( cdoVerbose )
        cdoPrint("nsmooth = %d, maxpoints = %d, radius = %gdegree, form = %s, weight0 = %g, weightR = %g",
                 xnsmooth, spoint.maxpoints, spoint.radius, Form[spoint.form], spoint.weight0, spoint.weightR);
    }

  spoint.radius *= DEG2RAD;

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
      int gridID = vlistInqVarGrid(vlistID1, varID);
      int gridtype = gridInqType(gridID);
      if ( gridtype == GRID_GAUSSIAN ||
           gridtype == GRID_LONLAT   ||
           gridtype == GRID_CURVILINEAR )
	{
	  varIDs[varID] = 1;
	}
      else if ( gridtype == GRID_UNSTRUCTURED && operatorID == SMOOTH )
        {
	  varIDs[varID] = 1;
        }
      else
	{
          char varname[CDI_MAX_NAME];
          vlistInqVarName(vlistID1, varID, varname);
	  varIDs[varID] = 0;
	  cdoWarning("Unsupported grid for variable %s", varname);
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
	      int gridID = vlistInqVarGrid(vlistID1, varID);

              for ( int i = 0; i < xnsmooth; ++i )
                {
                  if ( operatorID == SMOOTH )
                    smooth(gridID, missval, array1, array2, &nmiss, spoint);
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
