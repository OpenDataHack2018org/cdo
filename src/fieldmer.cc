/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2018 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/


#include "cdo_int.h"
#include <cdi.h>
#include "percentiles.h"


void merfun(field_type field1, field_type *field2, int function)
{
  // clang-format off
  switch (function)
    {
    case func_min:   mermin(field1, field2);    break;
    case func_max:   mermax(field1, field2);    break;
    case func_range: merrange(field1, field2);  break;
    case func_sum:   mersum(field1, field2);    break;
    case func_meanw: mermeanw(field1, field2);  break;
    case func_avgw:  meravgw(field1, field2);   break;
    case func_stdw:  merstdw(field1, field2);   break;
    case func_std1w: merstd1w(field1, field2);  break;
    case func_varw:  mervarw(field1, field2);   break;
    case func_var1w: mervar1w(field1, field2);  break;
    default: cdoAbort("function %d not implemented!", function);
    }
  // clang-format on
}


void mermin(field_type field1, field_type *field2)
{
  size_t rnmiss = 0;
  int    grid    = field1.grid;
  size_t nmiss   = field1.nmiss;
  double missval = field1.missval;
  double *array  = field1.ptr;
  double rmin = 0;

  size_t nx = gridInqXsize(grid);
  size_t ny = gridInqYsize(grid);

  for ( size_t i = 0; i < nx; i++ )
    {
      if ( nmiss > 0 )
	{
	  rmin = DBL_MAX;
	  for ( size_t j = 0; j < ny; j++ )
	    if ( !DBL_IS_EQUAL(array[j*nx+i], missval) )
	      if ( array[j*nx+i] < rmin ) rmin = array[j*nx+i];

	  if ( IS_EQUAL(rmin, DBL_MAX) )
	    {
	      rnmiss++;
	      rmin = missval;
	    }
	}
      else
	{
	  rmin = DBL_MAX;
	  for ( size_t j = 0; j < ny; j++ )
	    if ( array[j*nx+i] < rmin )  rmin = array[j*nx+i];
	}

      field2->ptr[i] = rmin;
    }

  field2->nmiss  = rnmiss;
}


void mermax(field_type field1, field_type *field2)
{
  size_t rnmiss = 0;
  int    grid    = field1.grid;
  size_t nmiss   = field1.nmiss;
  double missval = field1.missval;
  double *array  = field1.ptr;
  double rmax = 0;

  size_t nx = gridInqXsize(grid);
  size_t ny = gridInqYsize(grid);

  for ( size_t i = 0; i < nx; i++ )
    {
      if ( nmiss > 0 )
	{
	  rmax = -DBL_MAX;
	  for ( size_t j = 0; j < ny; j++ )
	    if ( !DBL_IS_EQUAL(array[j*nx+i], missval) )
	      if ( array[j*nx+i] > rmax ) rmax = array[j*nx+i];

	  if ( IS_EQUAL(rmax, -DBL_MAX) )
	    {
	      rnmiss++;
	      rmax = missval;
	    }
	}
      else
	{
	  rmax = DBL_MIN;
	  for ( size_t j = 0; j < ny; j++ )
	    if ( array[j*nx+i] > rmax )  rmax = array[j*nx+i];
	}

      field2->ptr[i] = rmax;
    }

  field2->nmiss  = rnmiss;
}


void merrange(field_type field1, field_type *field2)
{
  size_t rnmiss = 0;
  int    grid    = field1.grid;
  size_t nmiss   = field1.nmiss;
  double missval = field1.missval;
  double *array  = field1.ptr;
  double rmin = 0;
  double rmax = 0;
  double rrange = 0;

  size_t nx = gridInqXsize(grid);
  size_t ny = gridInqYsize(grid);

  for ( size_t i = 0; i < nx; i++ )
    {
      if ( nmiss > 0 )
	{
	  rmin =  DBL_MAX;
	  rmax = -DBL_MAX;
	  for ( size_t j = 0; j < ny; j++ )
	    if ( !DBL_IS_EQUAL(array[j*nx+i], missval) )
              {
		if      ( array[j*nx+i] < rmin ) rmin = array[j*nx+i];
                else if ( array[j*nx+i] > rmax ) rmax = array[j*nx+i];
              }

	  if ( IS_EQUAL(rmin, DBL_MAX) || IS_EQUAL(rmax, -DBL_MAX) )
	    {
	      rnmiss++;
	      rrange = missval;
	    }
	  else
	    {
	      rrange = rmax - rmin;
	    }
	}
      else
	{
	  rmin = DBL_MAX;
	  rmax = DBL_MIN;
	  for ( size_t j = 0; j < ny; j++ )
	    {
	      if      ( array[j*nx+i] < rmin )  rmin = array[j*nx+i];
	      else if ( array[j*nx+i] > rmax )  rmax = array[j*nx+i];
	    }

          rrange = rmax - rmin;
	}

      field2->ptr[i] = rrange;
    }

  field2->nmiss  = rnmiss;
}


void mersum(field_type field1, field_type *field2)
{
  size_t rnmiss = 0;
  double rsum = 0;

  size_t nx = gridInqXsize(field1.grid);
  size_t ny = gridInqYsize(field1.grid);

  std::vector<double> v(ny);

  for ( size_t i = 0; i < nx; ++i )
    {
      for ( size_t j = 0; j < ny; j++ ) v[j] = field1.ptr[j*nx+i];

      if ( field1.nmiss > 0 )
	{
          rsum = arraySumMV(ny, &v[0], field1.missval);
	  if ( DBL_IS_EQUAL(rsum, field1.missval) ) rnmiss++;
	}
      else
	{
          rsum = arraySum(ny, &v[0]);
	}

      field2->ptr[i] = rsum;
    }

  field2->nmiss  = rnmiss;
}


void mermeanw(field_type field1, field_type *field2)
{
  size_t rnmiss = 0;
  double rmean = 0;

  size_t nx = gridInqXsize(field1.grid);
  size_t ny = gridInqYsize(field1.grid);

  std::vector<double> v(ny);
  std::vector<double> w(ny);

  for ( size_t i = 0; i < nx; ++i )
    {
      for ( size_t j = 0; j < ny; j++ ) v[j] = field1.ptr[j*nx+i];
      for ( size_t j = 0; j < ny; j++ ) w[j] = field1.weight[j*nx+i];

      if ( field1.nmiss )
	{
          rmean = arrayWeightedMeanMV(ny, &v[0], &w[0], field1.missval);
	}
      else
	{
          rmean = arrayWeightedMean(ny, &v[0], &w[0], field1.missval);
	}

      if ( DBL_IS_EQUAL(rmean, field1.missval) ) rnmiss++;

      field2->ptr[i] = rmean;
    }

  field2->nmiss  = rnmiss;
}


void meravgw(field_type field1, field_type *field2)
{
  size_t rnmiss = 0;
  double ravg = 0;

  size_t nx = gridInqXsize(field1.grid);
  size_t ny = gridInqYsize(field1.grid);

  std::vector<double> v(ny);
  std::vector<double> w(ny);

  for ( size_t i = 0; i < nx; ++i )
    {
      for ( size_t j = 0; j < ny; j++ ) v[j] = field1.ptr[j*nx+i];
      for ( size_t j = 0; j < ny; j++ ) w[j] = field1.weight[j*nx+i];

      if ( field1.nmiss )
	{
          ravg = arrayWeightedAvgMV(ny, &v[0], &w[0], field1.missval);
	}
      else
	{
          ravg = arrayWeightedMean(ny, &v[0], &w[0], field1.missval);
	}

      if ( DBL_IS_EQUAL(ravg, field1.missval) ) rnmiss++;

      field2->ptr[i] = ravg;
    }

  field2->nmiss  = rnmiss;
}

static
void prevarsum_merw(const double *restrict array, const double *restrict w, size_t  nx, size_t  ny, size_t nmiss, 
                    double missval, double *restrict rsum, double *restrict rsumw, double *restrict rsumq, double *restrict rsumwq)
{ 
  *rsum   = 0;
  *rsumq  = 0;
  *rsumw  = 0;
  *rsumwq = 0;

  if ( nmiss > 0 )
    {
      for ( size_t  j = 0; j < ny; j++ )
        if ( !DBL_IS_EQUAL(array[j*nx], missval) &&
             !DBL_IS_EQUAL(w[j*nx], missval) )
          {
            *rsum   += w[j*nx] * array[j*nx];
            *rsumq  += w[j*nx] * array[j*nx] * array[j*nx];
            *rsumw  += w[j*nx];
            *rsumwq += w[j*nx] * w[j*nx];
          }
    }
  else
    {
      for ( size_t  j = 0; j < ny; j++ )
        {
          *rsum   += w[j*nx] * array[j*nx];
          *rsumq  += w[j*nx] * array[j*nx] * array[j*nx];
          *rsumw  += w[j*nx];
          *rsumwq += w[j*nx] * w[j*nx];
        }
    }
}


void mervarw(field_type field1, field_type *field2)
{
  size_t rnmiss = 0;
  int    grid    = field1.grid;
  size_t nmiss   = field1.nmiss;
  double missval = field1.missval;
  double *array  = field1.ptr;
  double *w      = field1.weight;
  double rsum = 0, rsumw = 0, rvar = 0;
  double rsumq = 0, rsumwq = 0;

  size_t nx = gridInqXsize(grid);
  int ny = gridInqYsize(grid);

  for ( size_t i = 0; i < nx; i++ )
    {
      prevarsum_merw(array+i, w+i, nx, ny, nmiss, missval, &rsum, &rsumw, &rsumq, &rsumwq);

      rvar = IS_NOT_EQUAL(rsumw, 0) ? (rsumq*rsumw - rsum*rsum) / (rsumw*rsumw) : missval;
      if ( rvar < 0 && rvar > -1.e-5 ) rvar = 0;

      if ( DBL_IS_EQUAL(rvar, missval) ) rnmiss++;

      field2->ptr[i] = rvar;
    }

  field2->nmiss  = rnmiss;
}


void mervar1w(field_type field1, field_type *field2)
{
  size_t rnmiss = 0;
  int    grid    = field1.grid;
  size_t nmiss   = field1.nmiss;
  double missval = field1.missval;
  double *array  = field1.ptr;
  double *w      = field1.weight;
  double rsum = 0, rsumw = 0, rvar = 0;
  double rsumq = 0, rsumwq = 0;

  size_t nx = gridInqXsize(grid);
  int ny = gridInqYsize(grid);

  for ( size_t i = 0; i < nx; i++ )
    {
      prevarsum_merw(array+i, w+i, nx, ny, nmiss, missval, &rsum, &rsumw, &rsumq, &rsumwq);

      rvar = (rsumw*rsumw > rsumwq) ? (rsumq*rsumw - rsum*rsum) / (rsumw*rsumw - rsumwq) : missval;
      if ( rvar < 0 && rvar > -1.e-5 ) rvar = 0;

      if ( DBL_IS_EQUAL(rvar, missval) ) rnmiss++;

      field2->ptr[i] = rvar;
    }

  field2->nmiss  = rnmiss;
}


void merstdw(field_type field1, field_type *field2)
{
  size_t rnmiss = 0;
  int    grid    = field1.grid;
  double missval = field1.missval;
  double rstd;

  size_t nx = gridInqXsize(grid);

  mervarw(field1, field2);

  for ( size_t i = 0; i < nx; i++ )
    {
      rstd = var_to_std(field2->ptr[i], missval);

      if ( DBL_IS_EQUAL(rstd, missval) ) rnmiss++;

      field2->ptr[i] = rstd;
    }

  field2->nmiss  = rnmiss;
}


void merstd1w(field_type field1, field_type *field2)
{
  size_t rnmiss = 0;
  int    grid    = field1.grid;
  double missval = field1.missval;
  double rstd;

  size_t nx = gridInqXsize(grid);

  mervar1w(field1, field2);

  for ( size_t i = 0; i < nx; i++ )
    {
      rstd = var_to_std(field2->ptr[i], missval);

      if ( DBL_IS_EQUAL(rstd, missval) ) rnmiss++;

      field2->ptr[i] = rstd;
    }

  field2->nmiss  = rnmiss;
}

/* RQ */
void merpctl(field_type field1, field_type *field2, int p)
{
  size_t rnmiss = 0;
  int    grid    = field1.grid;
  size_t nmiss   = field1.nmiss;
  double missval = field1.missval;
  double *array  = field1.ptr;

  size_t nx = gridInqXsize(grid);
  size_t ny = gridInqYsize(grid);
  
  double *array2 = (double*) Malloc(nx*sizeof(double));
  
  if ( nmiss > 0 )
    {
      for ( size_t i = 0; i < nx; i++ )
        {
          size_t l = 0;
          for ( size_t j = 0; j < ny; j++ )
	    if ( !DBL_IS_EQUAL(array[j*nx+i], missval) )
	      array2[l++] = array[j*nx+i];
	    
          if ( l > 0 )
            {
              field2->ptr[i] = percentile(array2, l, p);
            }
          else
            {
              field2->ptr[i] = missval;
              rnmiss++;
            }
        }
    }
  else
    {
      for ( size_t i = 0; i < nx; i++ )
      	{
          if ( ny > 0 )
            {
              for ( size_t j = 0; j < ny; j++ )
                array2[j] = array[j*nx+i];
              field2->ptr[i] = percentile(array2, ny, p);
            }
          else
            {
              field2->ptr[i] = missval;
              rnmiss++;
            }
      	}
    }

  Free(array2);

  field2->nmiss = rnmiss;
}
/* QR */
