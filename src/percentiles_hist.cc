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
#include <stdlib.h>
#include <math.h>

#include <cdi.h>

#include "cdo_int.h"
#include "percentiles.h"
#include "percentiles_hist.h"

#define NBINS_DEFAULT 101
#define NBINS_MINIMUM 11

#define DBL_CAPACITY(n) ((int) (((n) * sizeof(int)) / sizeof(double)))
#define DBL_PTR(p) ((double *) (p))
#define INT_PTR(p) ((int *) (p))

static int
histGetEnvNBins()
{
  const char *str = getenv("CDO_PCTL_NBINS");

  return str != NULL ? MAX(atoi(str), NBINS_MINIMUM) : NBINS_DEFAULT;
}

static void
histDefBounds(HISTOGRAM *hist, double a, double b)
{
  int i;

  assert(hist != NULL);
  assert(hist->nbins > 0);

  hist->min = MIN(a, b);
  hist->max = MAX(a, b);
  hist->step = (hist->max - hist->min) / hist->nbins;
  hist->nsamp = 0;

  for (i = 0; i < hist->nbins; i++) INT_PTR(hist->ptr)[i] = 0;
}

static void
histBinValue(HISTOGRAM *hist, double value)
{
  int bin;

  assert(hist->step > 0);

  bin = MIN((int) ((value - hist->min) / hist->step), hist->nbins - 1);
  if (bin >= 0 && bin < hist->nbins) INT_PTR(hist->ptr)[bin]++;
}

static void
histBin(HISTOGRAM *hist)
{
  int i;
  double *values;

  assert(hist->nsamp == DBL_CAPACITY(hist->nbins));

  values = (double *) Malloc(hist->nsamp * sizeof(double));

  for (i = 0; i < hist->nsamp; i++) values[i] = DBL_PTR(hist->ptr)[i];
  for (i = 0; i < hist->nbins; i++) INT_PTR(hist->ptr)[i] = 0;
  for (i = 0; i < hist->nsamp; i++) histBinValue(hist, values[i]);

  Free(values);
}

static int
histAddValue(HISTOGRAM *hist, double value)
{
  assert(hist != NULL);
  assert(hist->nbins > 0);

  /* 2011-08-01 Uwe Schulzweida: added check for rounding errors */
  if (value < hist->min && (hist->min - value) < 1e5) value = hist->min;
  if (value > hist->max && (value - hist->max) < 1e5) value = hist->max;

  if (IS_EQUAL(hist->min, hist->max)) return 0;
  if (value < hist->min || value > hist->max) return 1;

  if (hist->nsamp < DBL_CAPACITY(hist->nbins))
    {
      DBL_PTR(hist->ptr)[hist->nsamp] = value;
      hist->nsamp++;
    }
  else if (hist->nsamp > DBL_CAPACITY(hist->nbins))
    {
      histBinValue(hist, value);
      hist->nsamp++;
    }
  else
    {
      histBin(hist);
      histBinValue(hist, value);
      hist->nsamp++;
    }

  return 0;
}

static double
histGetPercentile(const HISTOGRAM *hist, double p)
{
  int i = 0, count = 0;

  assert(hist != NULL);
  assert(hist->nsamp > 0);
  assert(hist->nbins > 0);
  assert(p >= 0);
  assert(p <= 100);

  if (hist->nsamp > DBL_CAPACITY(hist->nbins))
    {
      double s = hist->nsamp * (p / 100.0);

      do
        count += INT_PTR(hist->ptr)[i++];
      while (count < s);

      assert(i > 0);
      assert(i - 1 < hist->nbins);
      assert(INT_PTR(hist->ptr)[i - 1] > 0);
      assert(hist->step > 0.0);

      double t = (count - s) / INT_PTR(hist->ptr)[i - 1];
      return hist->min + (i - t) * hist->step;
    }
  else
    {
      return percentile(DBL_PTR(hist->ptr), hist->nsamp, p);
    }
}

HISTOGRAM_SET *
hsetCreate(int nvars)
{
  int varID;
  HISTOGRAM_SET *hset;

  assert(nvars > 0);

  hset = (HISTOGRAM_SET *) Malloc(sizeof(HISTOGRAM_SET));
  if (hset == NULL) cdoAbort("Not enough memory (%s)", __func__);

  hset->nvars = nvars;
  hset->nlevels = (int *) Malloc(nvars * sizeof(int));
  hset->grids = (int *) Malloc(nvars * sizeof(int));
  hset->histograms = (HISTOGRAM ***) Malloc(nvars * sizeof(HISTOGRAM **));
  if (hset->histograms == NULL) cdoAbort("Not enough memory (%s)", __func__);

  for (varID = 0; varID < nvars; varID++)
    {
      hset->nlevels[varID] = 0;
      hset->grids[varID] = 0;
      hset->histograms[varID] = NULL;
    }

  return hset;
}

void
hsetCreateVarLevels(HISTOGRAM_SET *hset, int varID, int nlevels, int grid)
{
  int nvars, nhists, nbins, levelID, histID;
  HISTOGRAM *hists;

  nbins = histGetEnvNBins();

  assert(hset != NULL);
  assert(nlevels > 0);
  assert(nbins > 0);

  nvars = hset->nvars;

  assert(nvars > 0);

  if (varID < 0 || varID >= nvars) cdoAbort("Illegal argument: varID %d is undefined (%s)", varID, __func__);

  nhists = gridInqSize(grid);

  assert(nhists > 0);

  hset->nlevels[varID] = nlevels;
  hset->grids[varID] = grid;

  hset->histograms[varID] = (HISTOGRAM **) Malloc(nlevels * sizeof(HISTOGRAM *));
  if (hset->histograms[varID] == NULL) cdoAbort("Not enough memory (%s)", __func__);

  for (levelID = 0; levelID < nlevels; levelID++)
    {
      hists = hset->histograms[varID][levelID] = (HISTOGRAM *) Malloc(nhists * sizeof(HISTOGRAM));
      if (hists == NULL) cdoAbort("Not enough memory (%s)", __func__);

      for (histID = 0; histID < nhists; histID++)
        {
          hists[histID].min = 0.0;
          hists[histID].max = 0.0;
          hists[histID].step = 0.0;
          hists[histID].nbins = nbins;
          hists[histID].nsamp = 0;

          hists[histID].ptr = (int *) Malloc(nbins * sizeof(int));
          if (hists[histID].ptr == NULL) cdoAbort("Not enough memory (%s)", __func__);
        }
    }
}

void
hsetDestroy(HISTOGRAM_SET *hset)
{
  int varID, levelID, histID, nhists;

  if (hset != NULL)
    {
      for (varID = hset->nvars; varID-- > 0;)
        {
          nhists = gridInqSize(hset->grids[varID]);
          for (levelID = hset->nlevels[varID]; levelID-- > 0;)
            {
              for (histID = nhists; histID-- > 0;) Free(hset->histograms[varID][levelID][histID].ptr);
              Free(hset->histograms[varID][levelID]);
            }
          Free(hset->histograms[varID]);
        }

      Free(hset->histograms);
      Free(hset->grids);
      Free(hset->nlevels);
      Free(hset);

      hset = NULL;
    }
}

void
hsetDefVarLevelBounds(HISTOGRAM_SET *hset, int varID, int levelID, const Field *field1, const Field *field2)
{
  const double *array1 = field1->ptr;
  const double *array2 = field2->ptr;

  int i, nvars, nlevels, nhists, grid;
  double a, b;
  HISTOGRAM *hists;

  assert(hset != NULL);
  assert(field1 != NULL);
  assert(field2 != NULL);
  assert(array1 != NULL);
  assert(array2 != NULL);

  nvars = hset->nvars;

  assert(nvars > 0);

  if (varID < 0 || varID >= nvars) cdoAbort("Illegal argument: varID %d is undefined (%s)", varID, __func__);

  nlevels = hset->nlevels[varID];

  assert(nlevels > 0);

  if (levelID < 0 || levelID >= nlevels) cdoAbort("Illegal argument: levelID %d is undefined (%s)", levelID, __func__);

  grid = hset->grids[varID];

  if (grid != field1->grid || grid != field2->grid) cdoAbort("Grids are different (%s)", __func__);

  hists = hset->histograms[varID][levelID];

  assert(hists != NULL);

  nhists = gridInqSize(grid);

  assert(nhists > 0);

  for (i = 0; i < nhists; i++)
    {
      a = array1[i];
      b = array2[i];

      if (DBL_IS_EQUAL(a, field1->missval) || DBL_IS_EQUAL(b, field2->missval) || DBL_IS_EQUAL(a, b))
        histDefBounds(&hists[i], 0.0, 0.0);
      else
        histDefBounds(&hists[i], a, b);
    }
}

void
hsetAddVarLevelValues(HISTOGRAM_SET *hset, int varID, int levelID, const Field *field)
{
  const double *array = field->ptr;
  int i, grid, nvars, nlevels, nhists, nign = 0;
  HISTOGRAM *hists;

  assert(hset != NULL);
  assert(field != NULL);
  assert(array != NULL);

  nvars = hset->nvars;

  assert(nvars > 0);

  if (varID < 0 || varID >= nvars) cdoAbort("Illegal argument: varID %d is undefined (%s)", varID, __func__);

  nlevels = hset->nlevels[varID];

  assert(nlevels > 0);

  if (levelID < 0 || levelID >= nlevels) cdoAbort("Illegal argument: levelID %d is undefined (%s)", levelID, __func__);

  grid = hset->grids[varID];
  if (grid != field->grid) cdoAbort("Grids are different (%s)", __func__);

  hists = hset->histograms[varID][levelID];

  assert(hists != NULL);

  nhists = gridInqSize(grid);

  assert(nhists > 0);

  if (field->nmiss)
    {
      for (i = 0; i < nhists; i++)
        if (!DBL_IS_EQUAL(array[i], field->missval)) nign += histAddValue(&hists[i], array[i]);
    }
  else
    {
      for (i = 0; i < nhists; i++) nign += histAddValue(&hists[i], array[i]);
    }

  if (nign) cdoWarning("%d out of %d grid values are out of bounds and have been ignored (%s)", nign, nhists, __func__);
}

void
hsetGetVarLevelPercentiles(Field *field, const HISTOGRAM_SET *hset, int varID, int levelID, double p)
{
  double *array = field->ptr;
  int i, nvars, nlevels, nhists, grid;
  HISTOGRAM *hists;

  assert(hset != NULL);
  assert(field != NULL);
  assert(array != NULL);

  nvars = hset->nvars;

  assert(nvars > 0);

  if (varID < 0 || varID >= nvars) cdoAbort("Illegal argument: varID %d is undefined (%s)", varID, __func__);

  nlevels = hset->nlevels[varID];

  assert(nlevels > 0);

  if (levelID < 0 || levelID >= nlevels) cdoAbort("Illegal argument: levelID %d is undefined (%s)", levelID, __func__);

  grid = hset->grids[varID];

  if (grid != field->grid) cdoAbort("Grids are different (%s)", __func__);

  hists = hset->histograms[varID][levelID];

  assert(hists != NULL);

  nhists = gridInqSize(grid);

  assert(nhists > 0);

  field->nmiss = 0;
  for (i = 0; i < nhists; i++)
    {
      if (hists[i].nsamp)
        {
          array[i] = histGetPercentile(&hists[i], p);
        }
      else
        {
          array[i] = field->missval;
          field->nmiss++;
        }
    }
}
