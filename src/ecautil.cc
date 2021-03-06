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

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "ecautil.h"

/**
 * Convert a Gregorian/Julian date to a Julian day number.
 *
 * The Gregorian calendar was adopted midday, October 15, 1582.
 */
static unsigned long
gregdate_to_julday(int year,  /* Gregorian year */
                   int month, /* Gregorian month (1-12) */
                   int day    /* Gregorian day (1-31) */
                   )
{
#if INT_MAX <= 0X7FFF
  long igreg = 15 + 31 * (10 + (12 * 1582));
  long iy; /* signed, origin 0 year */
  long ja; /* Julian century */
  long jm; /* Julian month */
  long jy; /* Julian year */
#else
  int igreg = 15 + 31 * (10 + (12 * 1582));
  int iy; /* signed, origin 0 year */
  int ja; /* Julian century */
  int jm; /* Julian month */
  int jy; /* Julian year */
#endif
  unsigned long julday; /* returned Julian day number */

  /*
   * Because there is no 0 BC or 0 AD, assume the user wants the start of
   * the common era if they specify year 0.
   */
  if (year == 0) year = 1;

  iy = year;
  if (year < 0) iy++;
  if (month > 2)
    {
      jy = iy;
      jm = month + 1;
    }
  else
    {
      jy = iy - 1;
      jm = month + 13;
    }

  /*
   *  Note: SLIGHTLY STRANGE CONSTRUCTIONS REQUIRED TO AVOID PROBLEMS WITH
   *        OPTIMISATION OR GENERAL ERRORS UNDER VMS!
   */
  julday = day + (int) (30.6001 * jm);
  if (jy >= 0)
    {
      julday += 365 * jy;
      julday += (unsigned long) (0.25 * jy);
    }
  else
    {
      double xi = 365.25 * jy;

      if ((int) xi != xi) xi -= 1;
      julday += (int) xi;
    }
  julday += 1720995;

  if (day + (31 * (month + (12 * iy))) >= igreg)
    {
      ja = jy / 100;
      julday -= ja;
      julday += 2;
      julday += ja / 4;
    }

  return julday;
}

/**
 * Computes the day-of-year correspnding a given Gregorian date.
 *
 * @param date a Gregorian date in the form YYYYMMDD
 *
 * @return the day-of-year
 */
unsigned long
day_of_year(int64_t date)
{
  const int year = date / 10000;
  const int month = (date - year * 10000) / 100;
  const int day = date - year * 10000 - month * 100;

  return gregdate_to_julday(year, month, day) - gregdate_to_julday(year, 1, 1) + 1;
}

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
 * a       b       n     b < n ? a : b > n ? a + 1 : a + n
 * a       miss    n     a
 * miss    b       n     b < n ? 0 : b
 * miss    miss    n     0
 *
 * @param field1 the 1st input field, also holds the result
 * @param field2 the 2nd input field
 * @param mode   the counting mode, must be an exact mathematical
 *               integer
 */
static void
count(Field *field1, const Field *field2, double mode)
{
  size_t i;
  const int grid1 = field1->grid;
  const size_t nmiss1 = field1->nmiss;
  const double missval1 = field1->missval;
  double *array1 = field1->ptr;
  const int grid2 = field2->grid;
  const double missval2 = field2->missval;
  const double *array2 = field2->ptr;

  size_t len = gridInqSize(grid1);

  if (len != gridInqSize(grid2)) cdoAbort("Fields have different gridsize (%s)", __func__);

  if (nmiss1 > 0)
    {
      for (i = 0; i < len; i++)
        {
          if (DBL_IS_EQUAL(array2[i], missval2))
            {
              if (IS_EQUAL(mode, 1.0) || DBL_IS_EQUAL(array1[i], missval1)) array1[i] = 0.0;
              continue;
            }

          if (!DBL_IS_EQUAL(array1[i], missval1))
            {
              if (IS_EQUAL(mode, 0.0) || IS_EQUAL(mode, 1.0) || array2[i] > mode)
                array1[i] += 1.0;
              else if (DBL_IS_EQUAL(array2[i], mode))
                array1[i] += mode;
            }
          else
            {
              if (IS_EQUAL(mode, 0.0) || IS_EQUAL(mode, 1.0))
                array1[i] = 1.0;
              else if (array2[i] < mode)
                array1[i] = 0.0;
              else
                array1[i] = array2[i];
            }
        }

      field1->nmiss = arrayNumMV(len, array1, missval1);
    }
  else
    {
      for (i = 0; i < len; i++)
        {
          if (DBL_IS_EQUAL(array2[i], missval2))
            {
              if (IS_EQUAL(mode, 1.0)) array1[i] = 0.0;
              continue;
            }

          if (IS_EQUAL(mode, 0.0) || IS_EQUAL(mode, 1.0) || array2[i] > mode)
            array1[i] += 1.0;
          else if (DBL_IS_EQUAL(array2[i], mode))
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
 * a       b       comp(a, b) ? a : miss
 * a       miss    miss
 * miss    b       miss
 * miss    miss    miss
 *
 * @param field1  the input field, also holds the result
 * @param field2  the reference field
 * @param compare the comparator
 */
static void
selcomp(Field *field1, const Field *field2, int (*compare)(double, double))
{
  size_t i;
  const int grid1 = field1->grid;
  const size_t nmiss1 = field1->nmiss;
  const double missval1 = field1->missval;
  double *array1 = field1->ptr;
  const int grid2 = field2->grid;
  const size_t nmiss2 = field2->nmiss;
  const double missval2 = field2->missval;
  const double *array2 = field2->ptr;

  size_t len = gridInqSize(grid1);

  if (len != gridInqSize(grid2)) cdoAbort("Fields have different gridsize (%s)", __func__);

  if (nmiss1 > 0 || nmiss2 > 0)
    {
      for (i = 0; i < len; i++)
        if (DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2) || !compare(array1[i], array2[i]))
          array1[i] = missval1;
    }
  else
    {
      for (i = 0; i < len; i++)
        if (!compare(array1[i], array2[i])) array1[i] = missval1;
    }

  field1->nmiss = arrayNumMV(len, array1, missval1);
}

/**
 * Selects all field elements that compare to a certain reference
 * value. The result of the operation is computed according to the
 * following rules:
 *
 * field  c      result
 * a      c      comp(a, c) ? a : miss
 * a      miss   miss
 * miss   c      miss
 * miss   miss   miss
 *
 * @param field   the input field, also holds the result
 * @param c       the refence value
 * @param compare the comparator
 */
static void
selcompc(Field *field, double c, int (*compare)(double, double))
{
  size_t i;
  const int grid = field->grid;
  const size_t nmiss = field->nmiss;
  const double missval = field->missval;
  double *array = field->ptr;

  size_t len = gridInqSize(grid);

  if (DBL_IS_EQUAL(c, missval))
    {
      for (i = 0; i < len; i++) array[i] = missval;
    }
  else if (nmiss > 0)
    {
      for (i = 0; i < len; i++)
        if (DBL_IS_EQUAL(array[i], missval) || !compare(array[i], c)) array[i] = missval;
    }
  else
    {
      for (i = 0; i < len; i++)
        if (!compare(array[i], c)) array[i] = missval;
    }

  field->nmiss = arrayNumMV(len, array, missval);
}

static int
le(double a, double b)
{
  return a <= b;
}

static int
lt(double a, double b)
{
  return a < b;
}

static int
ge(double a, double b)
{
  return a >= b;
}

static int
gt(double a, double b)
{
  return a > b;
}

static int
eq(double a, double b)
{
  return DBL_IS_EQUAL(a, b);
}

static int
ne(double a, double b)
{
  return !DBL_IS_EQUAL(a, b);
}

void
farnum(Field *field1, Field field2)
{
  count(field1, &field2, 0.0);
}

void
farnum2(Field *field1, Field field2)
{
  count(field1, &field2, 1.0);
}

void
farnum3(Field *field1, Field field2, double n)
{
  count(field1, &field2, n);
}

void
farsel(Field *field1, Field field2)
{
  size_t i;
  const int grid1 = field1->grid;
  const double missval1 = field1->missval;
  double *array1 = field1->ptr;
  const int grid2 = field2.grid;
  const size_t nmiss2 = field2.nmiss;
  const double missval2 = field2.missval;
  const double *array2 = field2.ptr;

  size_t len = gridInqSize(grid1);

  if (len != gridInqSize(grid2)) cdoAbort("Fields have different gridsize (%s)", __func__);

  if (nmiss2 > 0)
    {
      for (i = 0; i < len; i++)
        if (DBL_IS_EQUAL(array2[i], missval2) || DBL_IS_EQUAL(array2[i], 0.0)) array1[i] = missval1;
    }
  else
    {
      for (i = 0; i < len; i++)
        if (IS_EQUAL(array2[i], 0.0)) array1[i] = missval1;
    }

  field1->nmiss = arrayNumMV(len, array1, missval1);
}

void
farselle(Field *field1, Field field2)
{
  selcomp(field1, &field2, le);
}

void
farsellt(Field *field1, Field field2)
{
  selcomp(field1, &field2, lt);
}

void
farselge(Field *field1, Field field2)
{
  selcomp(field1, &field2, ge);
}

void
farselgt(Field *field1, Field field2)
{
  selcomp(field1, &field2, gt);
}

void
farseleq(Field *field1, Field field2)
{
  selcomp(field1, &field2, eq);
}

void
farselne(Field *field1, Field field2)
{
  selcomp(field1, &field2, ne);
}

void
farsellec(Field *field, double c)
{
  selcompc(field, c, le);
}

void
farselltc(Field *field, double c)
{
  selcompc(field, c, lt);
}

void
farselgec(Field *field, double c)
{
  selcompc(field, c, ge);
}

void
farseleqc(Field *field, double c)
{
  selcompc(field, c, eq);
}

void
farselnec(Field *field, double c)
{
  selcompc(field, c, ne);
}

void
farselgtc(Field *field, double c)
{
  selcompc(field, c, gt);
}

void
updateHist(Field *field[2], int nlevels, size_t gridsize, double *yvals, int onlyNorth)
{
  int levelID;

  for (levelID = 0; levelID < nlevels; levelID++)
    for (size_t i = 0; i < gridsize; i++)
      if (onlyNorth)
        {
          if (yvals[i] >= 0.0) field[1][levelID].ptr[i] = field[0][levelID].ptr[i];
        }
      else
        field[1][levelID].ptr[i] = field[0][levelID].ptr[i];
}

void
adjustEndDate(int nlevels, size_t gridsize, double *yvals, double missval, int64_t ovdate, Field *startDateWithHist[2],
              Field *endDateWithHist[2])
{
  int64_t ovdateSouth = MIN(cdiEncodeDate(ovdate / 10000, 6, 30), ovdate);

  for (int levelID = 0; levelID < nlevels; levelID++)
    {
      for (size_t i = 0; i < gridsize; i++)
        {
          /* start with southern sphere */
          if (yvals[i] < 0)
            {
              if (DBL_IS_EQUAL(startDateWithHist[1][levelID].ptr[i], missval))
                {
                  endDateWithHist[0][levelID].ptr[i] = missval;
                  continue;
                }
              if (DBL_IS_EQUAL(endDateWithHist[0][levelID].ptr[i], missval))
                {
                  endDateWithHist[0][levelID].ptr[i] = ovdateSouth;
                }
            }
          else
            {
              if (DBL_IS_EQUAL(startDateWithHist[0][levelID].ptr[i], missval))
                {
                  endDateWithHist[0][levelID].ptr[i] = missval;
                  continue;
                }

              if (DBL_IS_EQUAL(endDateWithHist[0][levelID].ptr[i], missval))
                {
                  endDateWithHist[0][levelID].ptr[i] = ovdate;
                }
            }
        }
    }
}

void
computeGsl(int nlevels, size_t gridsize, double *yvals, double missval, Field *startDateWithHist[2], Field *endDateWithHist[2],
           Field *gslDuration, Field *gslFirstDay, int useCurrentYear)
{
  int levelID;
  double firstDay, duration;

  if (!useCurrentYear)
    {
      for (levelID = 0; levelID < nlevels; levelID++)
        {
          for (size_t i = 0; i < gridsize; i++)
            {
              /* start with southern sphere */
              if (yvals[i] < 0.0)
                {
                  duration = (double) (date_to_julday(CALENDAR_PROLEPTIC, (int64_t) endDateWithHist[0][levelID].ptr[i])
                                       - date_to_julday(CALENDAR_PROLEPTIC, (int64_t) startDateWithHist[1][levelID].ptr[i]));
                }
              else
                {
                  duration = (double) (date_to_julday(CALENDAR_PROLEPTIC, (int64_t) endDateWithHist[1][levelID].ptr[i])
                                       - date_to_julday(CALENDAR_PROLEPTIC, (int64_t) startDateWithHist[1][levelID].ptr[i]));
                }

              if (DBL_IS_EQUAL(startDateWithHist[1][levelID].ptr[i], missval))
                firstDay = missval;
              else
                firstDay = (double) day_of_year((int64_t) startDateWithHist[1][levelID].ptr[i]);

              gslDuration[levelID].ptr[i] = duration;
              gslFirstDay[levelID].ptr[i] = firstDay;
            }
        }
    }
  else
    {
      /* the current year can only have values for the northern hemisphere */
      for (levelID = 0; levelID < nlevels; levelID++)
        {
          for (size_t i = 0; i < gridsize; i++)
            {
              /* start with southern sphere */
              if (yvals[i] < 0.0)
                {
                  gslDuration[levelID].ptr[i] = missval;
                  gslFirstDay[levelID].ptr[i] = missval;
                }
              else
                {
                  duration = (double) (date_to_julday(CALENDAR_PROLEPTIC, (int64_t) endDateWithHist[0][levelID].ptr[i])
                                       - date_to_julday(CALENDAR_PROLEPTIC, (int64_t) startDateWithHist[0][levelID].ptr[i]));

                  if (DBL_IS_EQUAL(startDateWithHist[0][levelID].ptr[i], missval))
                    firstDay = missval;
                  else
                    firstDay = (double) day_of_year((int64_t) startDateWithHist[0][levelID].ptr[i]);

                  gslDuration[levelID].ptr[i] = duration;
                  gslFirstDay[levelID].ptr[i] = firstDay;
                }
            }
        }
    }

  for (levelID = 0; levelID < nlevels; levelID++)
    {
      gslDuration[levelID].nmiss = arrayNumMV(gridsize, gslDuration[levelID].ptr, missval);
      gslFirstDay[levelID].nmiss = arrayNumMV(gridsize, gslFirstDay[levelID].ptr, missval);
    }
}

void
writeGslStream(int ostreamID, int otaxisID, int otsID, int ovarID1, int ovarID2, int ivlistID1, int first_var_id,
               Field *gslDuration, Field *gslFirstDay, int64_t vdate, int vtime, int nlevels)
{
  (void) ivlistID1;
  (void) first_var_id;

  taxisDefVdate(otaxisID, vdate);
  taxisDefVtime(otaxisID, vtime);
  pstreamDefTimestep(ostreamID, otsID);

  for (int levelID = 0; levelID < nlevels; levelID++)
    {
      pstreamDefRecord(ostreamID, ovarID1, levelID);
      pstreamWriteRecord(ostreamID, gslDuration[levelID].ptr, gslDuration[levelID].nmiss);
    }
  for (int levelID = 0; levelID < nlevels; levelID++)
    {
      pstreamDefRecord(ostreamID, ovarID2, levelID);
      pstreamWriteRecord(ostreamID, gslFirstDay[levelID].ptr, gslFirstDay[levelID].nmiss);
    }
}
