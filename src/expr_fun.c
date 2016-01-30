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

#include "cdo_int.h"
#include "grid.h"


void fld_field_init(field_t *field, size_t nmiss, double missval, size_t ngp, double *array, double *w)
{
  field_init(field);

  field->size    = ngp;
  field->nmiss   = nmiss;
  field->missval = missval;
  field->ptr     = array;
  field->weight  = w;
}


double *fld_weights(int gridID, size_t ngp)
{
  static bool lwarn = true;
  double *weights = (double*) Malloc(ngp*sizeof(double));
  for ( size_t i = 0; i < ngp; ++i ) weights[i] = 1;

  if ( ngp > 1 )
    {
      int wstatus = gridWeights(gridID, weights);
      if ( wstatus != 0 && lwarn )
        {
          lwarn = false;
          cdoWarning("Grid cell bounds not available, using constant grid cell area weights!");
        }
    }

  return weights;
}

/*
double fun_fldmin(size_t nmiss, double missval, size_t ngp, double *array, double *weights)
{
  field_t field;
  fld_field_init(&field, nmiss, missval, ngp, array, weights);
  
  return fldmin(field);
}


double fun_fldmax(size_t nmiss, double missval, size_t ngp, double *array, double *weights)
{
  field_t field;
  fld_field_init(&field, nmiss, missval, ngp, array, weights);
  
  return fldmax(field);
}


double fun_fldsum(size_t nmiss, double missval, size_t ngp, double *array, double *weights)
{
  field_t field;
  fld_field_init(&field, nmiss, missval, ngp, array, weights);
  
  return fldsum(field);
}


double fun_fldmean(size_t nmiss, double missval, size_t ngp, double *array, double *weights)
{
  field_t field;
  fld_field_init(&field, nmiss, missval, ngp, array, weights);
  
  return fldmean(field);
}


double fun_fldavg(size_t nmiss, double missval, size_t ngp, double *array, double *weights)
{
  field_t field;
  fld_field_init(&field, nmiss, missval, ngp, array, weights);
  
  return fldavg(field);
}
*/
