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
#include "merge_sort2.h"
#include "cdoOptions.h"
#include "array.h"


double fldfun(field_type field, int function)
{
  double rval = 0;

  // clang-format off
  switch (function)
    {
    case func_range:  rval = fldrange(field);  break;
    case func_min:    rval = fldmin(field);    break;
    case func_max:    rval = fldmax(field);    break;
    case func_sum:    rval = fldsum(field);    break;
    case func_mean:   rval = fldmean(field);   break;
    case func_avg:    rval = fldavg(field);    break;
    case func_std:    rval = fldstd(field);    break;
    case func_std1:   rval = fldstd1(field);   break;
    case func_var:    rval = fldvar(field);    break;
    case func_var1:   rval = fldvar1(field);   break;
    case func_meanw:  rval = fldmeanw(field);  break;
    case func_avgw:   rval = fldavgw(field);   break;
    case func_stdw:   rval = fldstdw(field);   break;
    case func_std1w:  rval = fldstd1w(field);  break;
    case func_varw:   rval = fldvarw(field);   break;
    case func_var1w:  rval = fldvar1w(field);  break;
    case func_brs:    rval = fldbrs(field);    break;
    case func_rank:   rval = fldrank(field);   break;
    case func_roc:    rval = fldroc(field);    break;
    case func_skew:   rval = fldskew(field);   break;
    case func_kurt:   rval = fldkurt(field);   break;
    default: cdoAbort("%s: function %d not implemented!", __func__, function);
    }
  // clang-format on

  return rval;
}


double fldrange(field_type field)
{
  double range;

  if ( field.nmiss )
    {
      range = arrayRangeMV(field.size, field.ptr, field.missval);
    }
  else
    {
      range = arrayRange(field.size, field.ptr);
    }

  return range;
}


double fldmin(field_type field)
{
  double rmin;

  if ( field.nmiss )
    {
      rmin = arrayMinMV(field.size, field.ptr, field.missval);
    }
  else
    {
      rmin = arrayMin(field.size, field.ptr);
    }

  return rmin;
}


double fldmax(field_type field)
{
  double rmax;

  if ( field.nmiss )
    {
      rmax = arrayMaxMV(field.size, field.ptr, field.missval);
    }
  else
    {
      rmax = arrayMax(field.size, field.ptr);
    }

  return rmax;
}


double fldsum(field_type field)
{
  double rsum = 0;

  if ( field.nmiss )
    {
      rsum = arraySumMV(field.size, field.ptr, field.missval);
    }
  else
    {
      rsum = arraySum(field.size, field.ptr);
    }

  return rsum;
}


double fldmean(field_type field)
{
  double rmean = 0;

  if ( field.nmiss )
    {
      rmean = arrayMeanMV(field.size, field.ptr, field.missval);
    }
  else
    {
      rmean = arrayMean(field.size, field.ptr);
    }

  return rmean;
}


double fldmeanw(field_type field)
{
  double rmean = 0;

  if ( field.nmiss )
    {
      rmean = arrayWeightedMeanMV(field.size, field.ptr, field.weight, field.missval);
    }
  else
    {
      rmean = arrayWeightedMean(field.size, field.ptr, field.weight, field.missval);
    }

  return rmean;
}


double fldavg(field_type field)
{
  double ravg = 0;

  if ( field.nmiss )
    {
      ravg = arrayMeanMV(field.size, field.ptr, field.missval);
    }
  else
    {
      ravg = arrayMean(field.size, field.ptr);
    }

  return ravg;
}


double fldavgw(field_type field)
{
  double ravg = 0;

  if ( field.nmiss )
    {
      ravg = arrayWeightedAvgMV(field.size, field.ptr, field.weight, field.missval);
    }
  else
    {
      ravg = arrayWeightedMean(field.size, field.ptr, field.weight, field.missval);
    }

  return ravg;
}

static
void prevarsum(const double *restrict array, size_t len, size_t nmiss,
               double missval, double *rsum, double *rsumw, double *rsumq, double *rsumwq)
{
  assert(array!=NULL);

  double xsum = 0, xsumw = 0;
  double xsumq = 0, xsumwq = 0;

  if ( nmiss )
    {
      for ( size_t i = 0; i < len; ++i )
        if ( !DBL_IS_EQUAL(array[i], missval) )
          {
            xsum   += array[i];
            xsumq  += array[i] * array[i];
            xsumw  += 1;
            xsumwq += 1;
          }
    }
  else
    {
      for ( size_t i = 0; i < len; ++i )
        {
          xsum   += array[i];
          xsumq  += array[i] * array[i];
        }
      xsumw = len;
      xsumwq = len;
    }

  *rsum   = xsum;
  *rsumq  = xsumq;
  *rsumw  = xsumw;
  *rsumwq = xsumwq;
}

static
void preskewsum(const double *restrict array, size_t len, size_t nmiss, const double mean,
               double missval, double *rsum3w, double *rsum4w,
               double *rsum3diff, double *rsum2diff)
{
  assert(array!=NULL);

  double xsum3w = 0, xsum3diff = 0;
  double xsum4w = 0, xsum2diff = 0;

  if ( nmiss )
    {
      for ( size_t i = 0; i < len; ++i )
        if ( !DBL_IS_EQUAL(array[i], missval) )
          {
            xsum3diff += (array[i]-mean) * (array[i]-mean) * (array[i]-mean);
            xsum2diff += (array[i]-mean) * (array[i]-mean);
            xsum3w    += 1;
            xsum4w    += 1;
          }
    }
  else
    {
      for ( size_t i = 0; i < len; ++i )
        {
          xsum3diff += (array[i]-mean) * (array[i]-mean) * (array[i]-mean);
          xsum2diff += (array[i]-mean) * (array[i]-mean);
        }
      xsum3w = len;
      xsum4w = len;
    }

  *rsum3diff = xsum3diff;
  *rsum2diff = xsum2diff;
  *rsum3w    = xsum3w;
  *rsum4w    = xsum4w;
}

static
void prekurtsum(const double *restrict array, size_t len, size_t nmiss, const double mean,
               double missval, double *rsum3w, double *rsum4w,
               double *rsum2diff, double *rsum4diff)
{
  assert(array!=NULL);

  double xsum3w = 0, xsum4diff = 0;
  double xsum4w = 0, xsum2diff = 0;

  if ( nmiss )
    {
      for ( size_t i = 0; i < len; ++i )
        if ( !DBL_IS_EQUAL(array[i], missval) )
          {
            xsum2diff += (array[i]-mean) * (array[i]-mean);
            xsum4diff += (array[i]-mean) * (array[i]-mean) * (array[i]-mean) * (array[i]-mean);
            xsum3w    += 1;
            xsum4w    += 1;
          }
    }
  else
    {
      for ( size_t i = 0; i < len; ++i )
        {
          xsum2diff += (array[i]-mean) * (array[i]-mean);
          xsum4diff += (array[i]-mean) * (array[i]-mean) * (array[i]-mean) * (array[i]-mean);
        }
      xsum3w = len;
      xsum4w = len;
    }

  *rsum4diff = xsum4diff;
  *rsum2diff = xsum2diff;
  *rsum3w    = xsum3w;
  *rsum4w    = xsum4w;
}

double fldvar(field_type field)
{
  const size_t nmiss   = field.nmiss > 0;
  const size_t len     = field.size;
  const double missval = field.missval;
  double rsum, rsumw;
  double rsumq, rsumwq;

  prevarsum(field.ptr, len, nmiss, missval, &rsum, &rsumw, &rsumq, &rsumwq);

  double rvar = IS_NOT_EQUAL(rsumw, 0) ? (rsumq*rsumw - rsum*rsum) / (rsumw*rsumw) : missval;
  if ( rvar < 0 && rvar > -1.e-5 ) rvar = 0;

  return rvar;
}


double fldvar1(field_type field)
{
  const size_t nmiss   = field.nmiss > 0;
  const size_t len     = field.size;
  const double missval = field.missval;
  double rsum, rsumw;
  double rsumq, rsumwq;

  prevarsum(field.ptr, len, nmiss, missval, &rsum, &rsumw, &rsumq, &rsumwq);

  double rvar = (rsumw*rsumw > rsumwq) ? (rsumq*rsumw - rsum*rsum) / (rsumw*rsumw - rsumwq) : missval;
  if ( rvar < 0 && rvar > -1.e-5 ) rvar = 0;

  return rvar;
}
double fldkurt(field_type field)
{
  const size_t nmiss   = field.nmiss > 0;
  const size_t len     = field.size;
  const double missval = field.missval;
  double rsum, rsumw;
  double rsumq, rsumwq;
  double rsum3w;  /* 3rd moment variables */
  double rsum4w;  /* 4th moment variables */
  double rsum2diff, rsum4diff;


  prevarsum(field.ptr, len, nmiss, missval, &rsum, &rsumw, &rsumq, &rsumwq);
  prekurtsum(field.ptr, len, nmiss, (rsum/rsumw), missval, &rsum3w, &rsum4w, &rsum2diff, &rsum4diff);

  if(IS_EQUAL(rsum3w,0.0)|| IS_EQUAL(rsum2diff,0.0)) return missval;
  double rvar = ((rsum4diff/rsum3w)/ pow((rsum2diff)/(rsum3w),2)) - 3.0;

  if ( rvar < 0 && rvar > -1.e-5 ) rvar = 0;

  return rvar;
}

double fldskew(field_type field)
{
  const size_t nmiss   = field.nmiss > 0;
  const size_t len     = field.size;
  const double missval = field.missval;
  double rsum, rsumw;
  double rsumq, rsumwq;
  double rsum3w;  /* 3rd moment variables */
  double rsum4w;  /* 4th moment variables */
  double rsum3diff, rsum2diff;

  prevarsum(field.ptr, len, nmiss, missval, &rsum, &rsumw, &rsumq, &rsumwq);
  preskewsum(field.ptr, len, nmiss, (rsum/rsumw), missval, &rsum3w, &rsum4w, &rsum3diff, &rsum2diff);

  if(IS_EQUAL(rsum3w, 0.0) || IS_EQUAL(rsum3w, 1.0) || IS_EQUAL(rsum2diff,0.0)) return missval; 
  double rvar = (rsum3diff/rsum3w)/pow((rsum2diff)/(rsum3w-1.0),1.5);

  if ( rvar < 0 && rvar > -1.e-5 ) rvar = 0;

  return rvar;
}

static
void prevarsumw(const double *restrict array, const double *restrict w, size_t len, size_t nmiss,
                double missval, double *rsum, double *rsumw, double *rsumq, double *rsumwq)
{
  assert(array!=NULL);
  assert(w!=NULL);

  double xsum = 0, xsumw = 0;
  double xsumq = 0, xsumwq = 0;

  if ( nmiss )
    {
      for ( size_t i = 0; i < len; ++i )
        if ( !DBL_IS_EQUAL(array[i], missval) && !DBL_IS_EQUAL(w[i], missval) )
          {
            xsum   += w[i] * array[i];
            xsumq  += w[i] * array[i] * array[i];
            xsumw  += w[i];
            xsumwq += w[i] * w[i];
          }
    }
  else
    {
      for ( size_t i = 0; i < len; ++i )
        {
          xsum   += w[i] * array[i];
          xsumq  += w[i] * array[i] * array[i];
          xsumw  += w[i];
          xsumwq += w[i] * w[i];
        }
    }

  *rsum   = xsum;
  *rsumq  = xsumq;
  *rsumw  = xsumw;
  *rsumwq = xsumwq;
}


double fldvarw(field_type field)
{
  const size_t nmiss   = field.nmiss > 0;
  const size_t len     = field.size;
  const double missval = field.missval;
  double rsum, rsumw;
  double rsumq, rsumwq;

  prevarsumw(field.ptr, field.weight, len, nmiss, missval, &rsum, &rsumw, &rsumq, &rsumwq);

  double rvar = IS_NOT_EQUAL(rsumw, 0) ? (rsumq*rsumw - rsum*rsum) / (rsumw*rsumw) : missval;
  if ( rvar < 0 && rvar > -1.e-5 ) rvar = 0;

  return rvar;
}


double fldvar1w(field_type field)
{
  const size_t nmiss   = field.nmiss > 0;
  const size_t len     = field.size;
  const double missval = field.missval;
  double rsum, rsumw;
  double rsumq, rsumwq;

  prevarsumw(field.ptr, field.weight, len, nmiss, missval, &rsum, &rsumw, &rsumq, &rsumwq);

  double rvar = (rsumw*rsumw > rsumwq) ? (rsumq*rsumw - rsum*rsum) / (rsumw*rsumw - rsumwq) : missval;
  if ( rvar < 0 && rvar > -1.e-5 ) rvar = 0;

  return rvar;
}


double varToStd(double rvar, double missval)
{
  double rstd;

  if ( DBL_IS_EQUAL(rvar, missval) || rvar < 0 )
    {
      rstd = missval;
    }
  else
    {
      rstd = IS_NOT_EQUAL(rvar, 0) ? sqrt(rvar) : 0;
    }

  return rstd;
}


double fldstd(field_type field)
{
  return varToStd(fldvar(field), field.missval);
}


double fldstd1(field_type field)
{
  return varToStd(fldvar1(field), field.missval);
}


double fldstdw(field_type field)
{
  return varToStd(fldvarw(field), field.missval);
}


double fldstd1w(field_type field)
{
  return varToStd(fldvar1w(field), field.missval);
}


void fldrms(field_type field, field_type field2, field_type *field3)
{
  size_t   i;
  size_t len;
  size_t rnmiss = 0;
  int    grid1    = field.grid;
  //  size_t nmiss1   = field.nmiss;
  double *array1  = field.ptr;
  int    grid2    = field2.grid;
  //  size_t nmiss2   = field2.nmiss;
  double *array2  = field2.ptr;
  const double missval1 = field.missval;
  const double missval2 = field2.missval;
  double *w       = field.weight;
  double rsum = 0, rsumw = 0, ravg = 0;

  len    = gridInqSize(grid1);
  if ( len != (size_t) gridInqSize(grid2) )
    cdoAbort("fields have different size!");

  /*
  if ( nmiss1 > 0 )
  */
    {
      for ( i = 0; i < len; i++ )
	if ( !DBL_IS_EQUAL(w[i], missval1) )
	  {
	    rsum  = ADDMN(rsum, MULMN(w[i], MULMN( SUBMN(array2[i], array1[i]),
                                            SUBMN(array2[i], array1[i]))));
	    rsumw = ADDMN(rsumw, w[i]);
	  }
    }
    /*
  else
    {
      for ( i = 0; i < len; i++ )
	{
	  rsum  += w[i] * array1[i];
	  rsumw += w[i];
	}
    }
    */

  ravg = SQRTMN( DIVMN(rsum, rsumw));

  if ( DBL_IS_EQUAL(ravg, missval1) ) rnmiss++;

  field3->ptr[0] = ravg;
  field3->nmiss  = rnmiss;
}


void varrms(field_type field, field_type field2, field_type *field3)
{
  size_t   i, k, nlev, len;
  size_t rnmiss = 0;
  int    zaxis    = field.zaxis;
  int    grid1    = field.grid;
  //  size_t nmiss1   = field.nmiss;
  double *array1  = field.ptr;
  int    grid2    = field2.grid;
  //  size_t nmiss2   = field2.nmiss;
  double *array2  = field2.ptr;
  const double missval1 = field.missval;
  const double missval2 = field2.missval;
  double *w       = field.weight;
  double rsum = 0, rsumw = 0, ravg = 0;

  nlev   = zaxisInqSize(zaxis);
  len    = gridInqSize(grid1);
  if ( len != (size_t) gridInqSize(grid2) )
    cdoAbort("fields have different size!");

  /*
  if ( nmiss1 > 0 )
  */
    {
      for ( k = 0; k < nlev; k++ )
	for ( i = 0; i < len; i++ )
	  /*	  if ( !DBL_IS_EQUAL(w[i], missval1) ) */
	    {
	      rsum  = ADDMN(rsum, MULMN(w[i], MULMN( SUBMN(array2[k*len+i], array1[k*len+i]),
                                              SUBMN(array2[k*len+i], array1[k*len+i]))));
	      rsumw = ADDMN(rsumw, w[i]);
	    }
    }
    /*
  else
    {
      for ( i = 0; i < len; i++ )
	{
	  rsum  += w[i] * array1[i];
	  rsumw += w[i];
	}
    }
    */

  ravg = SQRTMN( DIVMN(rsum, rsumw));

  if ( DBL_IS_EQUAL(ravg, missval1) ) rnmiss++;

  field3->ptr[0] = ravg;
  field3->nmiss  = rnmiss;
}

/* RQ */
double fldpctl(field_type field, const double pn)
{
  double *array = field.ptr;
  double pctl = field.missval;

  if ( field.size - field.nmiss > 0 )
    {
      if ( field.nmiss )
        {
          double *array2 = (double*) Malloc((field.size - field.nmiss)*sizeof(double));

          size_t j = 0;
          for ( size_t i = 0; i < field.size; i++ )
            if ( !DBL_IS_EQUAL(array[i], field.missval) )
              array2[j++] = array[i];

          pctl = percentile(array2, j, pn);

          Free(array2);
        }
      else
        {
          pctl = percentile(array, field.size, pn);
        }
    }

  return pctl;
}
/* QR */

double fldrank(field_type field)
{
  double res = 0;
  // Using first value as reference (observation)
  double *array  =  &(field.ptr[1]);
  double val     = array[-1];
  const double missval = field.missval;
  size_t nmiss      = field.nmiss;
  const size_t len       = field.size-1;
  size_t j;

  if ( nmiss ) return missval;

  sort_iter_single(len,array, 1);

  if ( val > array[len-1] )
    res=(double)len;
  else
    for ( j=0; j<len; j++ )
      if ( array[j] >= val ) {
	res=(double)j;
	break;
      }

  return res;
}


double fldroc(field_type field)
{
  return field.missval;
}


double fldbrs(field_type field)
{
  const size_t  nmiss   = field.nmiss;
  const size_t    len   = field.size;
  double *array   = field.ptr;
  const double missval  = field.missval;

  double brs = 0;
  size_t i, count=0;

  // Using first value as reference
  if ( nmiss == 0 )
    {
      for ( i=1; i<len; i++ )
	brs += (array[i] - array[0]) * (array[i] - array[0]);
      count = i-1;
    }
  else
    {
      if ( DBL_IS_EQUAL(array[0], missval) ) return missval;

      for ( i=1; i<len; i++ )
	if ( !DBL_IS_EQUAL(array[i], missval) )
	  {
	    brs += (array[i] - array[0]) * (array[i] - array[0]);
	    count ++;
	  }
    }

  return brs/count;
}

/*  field_type UTILITIES */
/*  update the number non missing values */
void fldunm(field_type *field)
{
  field->nmiss = arrayNumMV(field->size, field->ptr, field->missval);
}

/*  check for non missval values */
int fldhvs(field_type *fieldPtr, const size_t nlevels)
{
  field_type field;

  for ( size_t level = 0; level < nlevels; level++ )
    {
      field = fieldPtr[level];
      if ( field.nmiss != field.size ) return TRUE;
    }

  return FALSE;
}
