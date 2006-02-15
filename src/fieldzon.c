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

#include "cdo.h"
#include "cdo_int.h"
#include "cdi.h"


void zonfun(FIELD field1, FIELD *field2, int function)
{
  if      ( function == func_min )  zonmin(field1, field2);
  else if ( function == func_max )  zonmax(field1, field2);  
  else if ( function == func_sum )  zonsum(field1, field2);  
  else if ( function == func_mean ) zonmean(field1, field2);  
  else if ( function == func_avg )  zonavg(field1, field2);  
  else if ( function == func_std )  zonstd(field1, field2);  
  else if ( function == func_var )  zonvar(field1, field2);
  else cdoAbort("function %d not implemented!", function);
}


void zonmin(FIELD field1, FIELD *field2)
{
  int i, j, nx, ny, rnmiss = 0;
  int    grid    = field1.grid;
  int    nmiss   = field1.nmiss;
  double missval = field1.missval;
  double *array  = field1.ptr;
  double rmin = 0;

  nx    = gridInqXsize(grid);
  ny    = gridInqYsize(grid);

  for ( j = 0; j < ny; j++ )
    {
      if ( nmiss > 0 )
	{
	  rmin = DBL_MAX;
	  for ( i = 0; i < nx; i++ )
	    if ( !DBL_IS_EQUAL(array[j*nx+i], missval) )
	      if ( array[j*nx+i] < rmin ) rmin = array[j*nx+i];

	  if ( DBL_IS_EQUAL(rmin, DBL_MAX) )
	    {
	      rnmiss++;
	      rmin = missval;
	    }
	}
      else
	{
	  rmin = array[j*nx];
	  for ( i = 1; i < nx; i++ )
	    if ( array[j*nx+i] < rmin )  rmin = array[j*nx+i];
	}

      field2->ptr[j] = rmin;
    }

  field2->nmiss  = rnmiss;
}


void zonmax(FIELD field1, FIELD *field2)
{
  int i, j, nx, ny, rnmiss = 0;
  int    grid    = field1.grid;
  int    nmiss   = field1.nmiss;
  double missval = field1.missval;
  double *array  = field1.ptr;
  double rmax = 0;

  nx    = gridInqXsize(grid);
  ny    = gridInqYsize(grid);

  for ( j = 0; j < ny; j++ )
    {
      if ( nmiss > 0 )
	{
	  rmax = DBL_MIN;
	  for ( i = 0; i < nx; i++ )
	    if ( !DBL_IS_EQUAL(array[j*nx+i], missval) )
	      if ( array[j*nx+i] > rmax ) rmax = array[j*nx+i];

	  if ( DBL_IS_EQUAL(rmax, DBL_MIN) )
	    {
	      rnmiss++;
	      rmax = missval;
	    }
	}
      else
	{
	  rmax = array[j*nx];
	  for ( i = 1; i < nx; i++ ) 
	    if ( array[j*nx+i] > rmax )  rmax = array[j*nx+i];
	}

      field2->ptr[j] = rmax;
    }

  field2->nmiss  = rnmiss;
}


void zonsum(FIELD field1, FIELD *field2)
{
  int i, j, nx, ny, rnmiss = 0;
  int    grid    = field1.grid;
  int    nmiss   = field1.nmiss;
  double missval = field1.missval;
  double *array  = field1.ptr;
  double rsum = 0;

  nx    = gridInqXsize(grid);
  ny    = gridInqYsize(grid);

  for ( j = 0; j < ny; j++ )
    {
      if ( nmiss > 0 )
	{
	  rsum = 0;
	  for ( i = 0; i < nx; i++ )
	    if ( !DBL_IS_EQUAL(array[j*nx+i], missval) )
	      rsum += array[j*nx+i];
	}
      else
	{
	  rsum = 0;
	  for ( i = 0; i < nx; i++ )
	    rsum += array[j*nx+i];
	}

      field2->ptr[j] = rsum;
    }

  field2->nmiss  = rnmiss;
}


void zonmean(FIELD field1, FIELD *field2)
{
  int i, j, nx, ny, rnmiss = 0;
  int    grid     = field1.grid;
  int    nmiss    = field1.nmiss;
  double missval1 = field1.missval;
  double missval2 = field1.missval;
  double *array   = field1.ptr;
  double rsum = 0, rsumw = 0, ravg = 0;

  nx    = gridInqXsize(grid);
  ny    = gridInqYsize(grid);

  for ( j = 0; j < ny; j++ )
    {
      rsum  = 0;
      rsumw = 0;
      if ( nmiss > 0 )
	{
	  for ( i = 0; i < nx; i++ )
	    if ( !DBL_IS_EQUAL(array[j*nx+i], missval1) )
	      {
		rsum  += array[j*nx+i];
		rsumw += 1;
	      }
	}
      else
	{
	  for ( i = 0; i < nx; i++ )
	    {
	      rsum  += array[j*nx+i];
	      rsumw += 1;
	    }
	}

      ravg = DIV(rsum, rsumw);

      if ( DBL_IS_EQUAL(ravg, missval1) ) rnmiss++;

      field2->ptr[j] = ravg;
    }

  field2->nmiss  = rnmiss;
}


void zonavg(FIELD field1, FIELD *field2)
{
  int i, j, nx, ny, rnmiss = 0;
  int    grid     = field1.grid;
  int    nmiss    = field1.nmiss;
  double missval1 = field1.missval;
  double missval2 = field1.missval;
  double *array   = field1.ptr;
  double rsum = 0, rsumw = 0, ravg = 0;

  nx    = gridInqXsize(grid);
  ny    = gridInqYsize(grid);

  for ( j = 0; j < ny; j++ )
    {
      rsum  = 0;
      rsumw = 0;
      if ( nmiss > 0 )
	{
	  for ( i = 0; i < nx; i++ )
	    {
	      rsum   = ADD(rsum, array[j*nx+i]);
	      rsumw += 1;
	    }
	}
      else
	{
	  for ( i = 0; i < nx; i++ )
	    {
	      rsum  += array[j*nx+i];
	      rsumw += 1;
	    }
	}

      ravg = DIV(rsum, rsumw);

      if ( DBL_IS_EQUAL(ravg, missval1) ) rnmiss++;

      field2->ptr[j] = ravg;
    }

  field2->nmiss  = rnmiss;
}


void zonvar(FIELD field1, FIELD *field2)
{
  int i, j, nx, ny, rnmiss = 0;
  int    grid     = field1.grid;
  int    nmiss    = field1.nmiss;
  double missval1 = field1.missval;
  double *array   = field1.ptr;
  double rsum = 0, rsumw = 0, rvar = 0;
  double rsumq = 0, rsumwq = 0;

  nx    = gridInqXsize(grid);
  ny    = gridInqYsize(grid);

  for ( j = 0; j < ny; j++ )
    {
      rsum   = 0;
      rsumq  = 0;
      rsumw  = 0;
      rsumwq = 0;
      if ( nmiss > 0 )
	{
	  for ( i = 0; i < nx; i++ )
	    if ( !DBL_IS_EQUAL(array[j*nx+i], missval1) )
	      {
		rsum   += array[j*nx+i];
		rsumq  += array[j*nx+i] * array[j*nx+i];
		rsumw  += 1;
		rsumwq += 1;
	      }
	}
      else
	{
	  for ( i = 0; i < nx; i++ )
	    {
	      rsum   += array[j*nx+i];
	      rsumq  += array[j*nx+i] * array[j*nx+i];
	      rsumw  += 1;
	      rsumwq += 1;
	    }
	}

      rvar = !DBL_IS_EQUAL(rsumw, 0) ? (rsumq*rsumw - rsum*rsum) / (rsumw*rsumw) : missval1;

      if ( DBL_IS_EQUAL(rvar, missval1) ) rnmiss++;

      field2->ptr[j] = rvar;
    }

  field2->nmiss  = rnmiss;
}


void zonstd(FIELD field1, FIELD *field2)
{
  int j, ny, rnmiss = 0;
  int    grid    = field1.grid;
  double missval = field1.missval;
  double rvar, rstd;

  ny    = gridInqYsize(grid);

  zonvar(field1, field2);

  for ( j = 0; j < ny; j++ )
    {
      rvar = field2->ptr[j];
      rstd = (!DBL_IS_EQUAL(rvar, 0) && !DBL_IS_EQUAL(rvar, missval)) ? sqrt(rvar) : missval;

      if ( DBL_IS_EQUAL(rvar, missval) ) rnmiss++;

      field2->ptr[j] = rstd;
    }

  field2->nmiss  = rnmiss;
}
