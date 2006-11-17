/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include <assert.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "ecautil.h"


/**
 * Counts the number of nonmissing values. The result of the operation
 * is computed according to the following rules:
 * 
 * field1  field2  mode  result
 * a       b       0     a + 1
 * a       miss    0     a
 * miss    b       0     1
 * miss    miss    0     0
 * 
 * a       b       1     a + 1
 * a       miss    1     0
 * miss    b       1     1
 * miss    miss    1     0
 * 
 * a       b       n     a + 1 if b > n; a + n if b = n; a if b < n
 * a       miss    n     a
 * miss    b       n     b if b > n; n if b = n; 0 if b < n
 * miss    miss    n     0    
 * 
 * @param field1 the 1st input field, also holds the result
 * @param field2 the 2nd input field
 * @param mode   the counting mode, must be an exact mathematical
 *               integer
 */  
static void num(FIELD *field1, const FIELD *field2, double mode)
{
  static const char func[] = "num";
  int   i, len;
  const int     grid1    = field1->grid;
  const int     nmiss1   = field1->nmiss;
  const double  missval1 = field1->missval;
  double       *array1   = field1->ptr;
  const int     grid2    = field2->grid;
  const int     nmiss2   = field2->nmiss;
  const double  missval2 = field2->missval;
  const double *array2   = field2->ptr;
  
  len = gridInqSize(grid1);

  if ( len != gridInqSize(grid2) )
    cdoAbort("Fields have different gridsize (%s)", func);

  if ( nmiss1 > 0 )
    {
      for ( i = 0; i < len; i++ )
        {
          if ( DBL_IS_EQUAL(array2[i], missval2) )
            {
              if ( mode == 1.0 || DBL_IS_EQUAL(array1[i], missval1) )
                array1[i] = 0.0; 
              continue;
            }

          if ( !DBL_IS_EQUAL(array1[i], missval1) )
            {
              if ( mode == 0.0 || mode == 1.0 || array2[i] > mode )
                array1[i] += 1.0;
              else if ( array2[i] == mode )
                array1[i] += mode;
            } 
          else
            {
              if ( mode == 0.0 || mode == 1.0 )
                array1[i] = 1.0;
              else if ( array2[i] < mode )
                array1[i] = 0.0;
              else
                array1[i] = array2[i];
            } 
        }

      field1->nmiss = 0;
      for ( i = 0; i < len; i++ )
        if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
    }
  else 
    {
      for ( i = 0; i < len; i++ )
        {
          if ( DBL_IS_EQUAL(array2[i], missval2) )
            {
              if ( mode == 1.0 )
                array1[i] = 0.0; 
              continue;
            }

          if ( mode == 0.0 || mode == 1.0 || array2[i] > mode )
            array1[i] += 1.0;
          else if ( array2[i] == mode )
            array1[i] += mode;
        }
    }
}


/**
 * Selects all field elements that compare to the corresponding
 * element of a reference field. The result of the operation is
 * computed according to the following rules:
 * 
 * field1  field2  result
 * a       b       comp(a, b) is true; miss otherwise
 * a       miss    miss
 * miss    b       miss
 * miss    miss    miss    
 * 
 * @param field1  the input field, also holds the result
 * @param field2  the reference field
 * @param compare the comparator
 */  
static void selcomp(FIELD *field1, const FIELD *field2, int (*compare)(double, double))
{
  static const char func[] = "selcomp";
  int   i, len;
  const int     grid1    = field1->grid;
  const int     nmiss1   = field1->nmiss;
  const double  missval1 = field1->missval;
  double       *array1   = field1->ptr;
  const int     grid2    = field2->grid;
  const int     nmiss2   = field2->nmiss;
  const double  missval2 = field2->missval;
  const double *array2   = field2->ptr;
  
  len = gridInqSize(grid1);

  if ( len != gridInqSize(grid2) )
    cdoAbort("Fields have different gridsize (%s)", func);

  if ( nmiss1 > 0 || nmiss2 > 0 )
    {
      for ( i = 0; i < len; i++ )
        if ( DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2) || !compare(array1[i], array2[i]) ) 
          array1[i] = missval1;
    }
  else
    {
      for ( i = 0; i < len; i++ )
        if ( !compare(array1[i], array2[i]) ) 
          array1[i] = missval1;
    }
      
  field1->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array1[i], missval1) ) field1->nmiss++;
}


/**
 * Selects all field elements that compare to a certain reference
 * value. The result of the operation is computed according to the
 * following rules:
 * 
 * field  c      result
 * a      c      comp(a, c) is true; miss otherwise
 * a      miss   miss
 * miss   c      miss
 * miss   miss   miss    
 * 
 * @param field   the input field, also holds the result
 * @param c       the refence value
 * @param compare the comparator
 */  
static void selcompc(FIELD *field, double c, int (*compare)(double, double))
{
  int   i, len;
  const int     grid    = field->grid;
  const int     nmiss   = field->nmiss;
  const double  missval = field->missval;
  double       *array   = field->ptr;
  
  len = gridInqSize(grid);

  if ( DBL_IS_EQUAL(c, missval) )
    {
      for ( i = 0; i < len; i++ )
        array[i] = missval;
    }
  else if ( nmiss > 0 )
    {
      for ( i = 0; i < len; i++ )
        if ( DBL_IS_EQUAL(array[i], missval) || !compare(array[i], c) ) 
          array[i] = missval;
    }
  else
    {
      for ( i = 0; i < len; i++ )
        if ( !compare(array[i], c) ) 
          array[i] = missval;
    }
      
  field->nmiss = 0;
  for ( i = 0; i < len; i++ )
    if ( DBL_IS_EQUAL(array[i], missval) ) field->nmiss++;
}


static int and(double a, double b)
{
  return a <= b;
}


static int le(double a, double b)
{
  return a <= b;
}


static int lt(double a, double b)
{
  return a < b;
}


static int ge(double a, double b)
{
  return a >= b;
}


static int gt(double a, double b)
{
  return a > b;
}


static int eq(double a, double b)
{
  return DBL_IS_EQUAL(a, b);
}


static int ne(double a, double b)
{
  return !DBL_IS_EQUAL(a, b);
}


void farnum(FIELD *field1, FIELD field2)
{
  return num(field1, &field2, 0.0);
}


void farnum2(FIELD *field1, FIELD field2)
{
  return num(field1, &field2, 1.0);
}


void farnum3(FIELD *field1, FIELD field2, double n)
{
  return num(field1, &field2, n);
}


void farselle(FIELD *field1, FIELD field2)
{
  return selcomp(field1, &field2, le);
}


void farsellt(FIELD *field1, FIELD field2)
{
  return selcomp(field1, &field2, lt);
}


void farselge(FIELD *field1, FIELD field2)
{
  return selcomp(field1, &field2, ge);
}


void farselgt(FIELD *field1, FIELD field2)
{
  return selcomp(field1, &field2, gt);
}


void farseleq(FIELD *field1, FIELD field2)
{
  return selcomp(field1, &field2, eq);
}


void farselne(FIELD *field1, FIELD field2)
{
  return selcomp(field1, &field2, ne);
}


void farsellec(FIELD *field, double c)
{
  return selcompc(field, c, le);
}


void farselltc(FIELD *field, double c)
{
  return selcompc(field, c, lt);
}


void farselgec(FIELD *field, double c)
{
  return selcompc(field, c, ge);
}


void farseleqc(FIELD *field, double c)
{
  return selcompc(field, c, eq);
}


void farselnec(FIELD *field, double c)
{
  return selcompc(field, c, ne);
}


void farselgtc(FIELD *field, double c)
{
  return selcompc(field, c, gt);
}
