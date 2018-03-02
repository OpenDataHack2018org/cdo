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

/*
   This module contains the following operators:

      Smoothstat       smooth9             running 9-point-average
*/
#include <time.h>  // clock()

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "grid.h"
#include "constants.h"  // planet radius
#include "pmlist.h"
#include "cdoOptions.h"

#include "grid_search.h"

enum
{
  FORM_LINEAR
};
static const char *Form[] = { "linear" };

typedef struct
{
  size_t maxpoints;
  int form;
  double radius;
  double weight0;
  double weightR;
} smoothpoint_t;

static void
smooth(int gridID, double missval, const double *restrict array1, double *restrict array2, size_t *nmiss,
       smoothpoint_t spoint)
{
  int gridID0 = gridID;
  size_t gridsize = gridInqSize(gridID);
  size_t numNeighbors = spoint.maxpoints;
  if (numNeighbors > gridsize) numNeighbors = gridsize;

  uint8_t *mask = (uint8_t *) Malloc(gridsize * sizeof(uint8_t));
  for (size_t i = 0; i < gridsize; ++i) mask[i] = !DBL_IS_EQUAL(array1[i], missval);

  std::vector<double> xvals(gridsize);
  std::vector<double> yvals(gridsize);

  if (gridInqType(gridID) == GRID_GME) gridID = gridToUnstructured(gridID, 0);

  if (gridInqType(gridID) != GRID_UNSTRUCTURED && gridInqType(gridID) != GRID_CURVILINEAR)
    gridID = gridToCurvilinear(gridID, 0);

  gridInqXvals(gridID, &xvals[0]);
  gridInqYvals(gridID, &yvals[0]);

  /* Convert lat/lon units if required */
  char units[CDI_MAX_NAME];
  gridInqXunits(gridID, units);
  grid_to_radian(units, gridsize, &xvals[0], "grid center lon");
  gridInqYunits(gridID, units);
  grid_to_radian(units, gridsize, &yvals[0], "grid center lat");

  std::vector<knnWeightsType> knnWeights;
  for (int i = 0; i < Threading::ompNumThreads; ++i) knnWeights.push_back(knnWeightsType(numNeighbors));

  clock_t start, finish;

  start = clock();

  bool xIsCyclic = false;
  size_t dims[2] = { gridsize, 0 };
  struct gridsearch *gs = gridsearch_create(xIsCyclic, dims, gridsize, &xvals[0], &yvals[0]);

  gs->search_radius = spoint.radius;

  finish = clock();

  if (cdoVerbose) printf("gridsearch created: %.2f seconds\n", ((double) (finish - start)) / CLOCKS_PER_SEC);

  if (cdoVerbose) progressInit();

  start = clock();

  size_t nmissx = 0;
  double findex = 0;

#ifdef HAVE_OPENMP4
#pragma omp parallel for schedule(dynamic) default(none) reduction(+ : findex) reduction(+ : nmissx) shared( \
    cdoVerbose, knnWeights, spoint, mask, array1, array2, xvals, yvals, gs, gridsize, missval)
#endif
  for (size_t i = 0; i < gridsize; ++i)
    {
      int ompthID = cdo_omp_get_thread_num();

      findex++;
      if (cdoVerbose && cdo_omp_get_thread_num() == 0) progressStatus(0, 1, findex / gridsize);

      grid_search_nbr(gs, knnWeights[ompthID], xvals[i], yvals[i]);

      // Compute weights based on inverse distance if mask is false, eliminate those points
      size_t nadds = knnWeights[ompthID].compute_weights(mask, spoint.radius, spoint.weight0, spoint.weightR);
      if (nadds)
        {
          array2[i] = knnWeights[ompthID].array_weights_sum(array1);
        }
      else
        {
          nmissx++;
          array2[i] = missval;
        }
    }

  *nmiss = nmissx;

  finish = clock();

  if (cdoVerbose) printf("gridsearch nearest: %.2f seconds\n", ((double) (finish - start)) / CLOCKS_PER_SEC);

  if (gs) gridsearch_delete(gs);

  if (gridID0 != gridID) gridDestroy(gridID);

  Free(mask);
}

static inline void
smooth9_sum(size_t ij, uint8_t *mask, double sfac, const double *restrict array, double *avg, double *divavg)
{
  if (mask[ij])
    {
      *avg += sfac * array[ij];
      *divavg += sfac;
    }
}

static void
smooth9(int gridID, double missval, const double *restrict array1, double *restrict array2, size_t *nmiss)
{
  size_t gridsize = gridInqSize(gridID);
  size_t nlon = gridInqXsize(gridID);
  size_t nlat = gridInqYsize(gridID);
  int grid_is_cyclic = gridIsCircular(gridID);

  uint8_t *mask = (uint8_t *) Malloc(gridsize * sizeof(uint8_t));

  for (size_t i = 0; i < gridsize; ++i) mask[i] = !DBL_IS_EQUAL(missval, array1[i]);

  *nmiss = 0;
  for (size_t i = 0; i < nlat; i++)
    {
      for (size_t j = 0; j < nlon; j++)
        {
          double avg = 0;
          double divavg = 0;

          if ((i == 0) || (j == 0) || (i == (nlat - 1)) || (j == (nlon - 1)))
            {
              size_t ij = j + nlon * i;
              if (mask[ij])
                {
                  avg += array1[ij];
                  divavg += 1;
                  /* upper left corner */
                  if ((i != 0) && (j != 0))
                    smooth9_sum(((i - 1) * nlon) + j - 1, mask, 0.3, array1, &avg, &divavg);
                  else if (i != 0 && grid_is_cyclic)
                    smooth9_sum((i - 1) * nlon + j - 1 + nlon, mask, 0.3, array1, &avg, &divavg);

                  /* upper cell */
                  if (i != 0) smooth9_sum(((i - 1) * nlon) + j, mask, 0.5, array1, &avg, &divavg);

                  /* upper right corner */
                  if ((i != 0) && (j != (nlon - 1)))
                    smooth9_sum(((i - 1) * nlon) + j + 1, mask, 0.3, array1, &avg, &divavg);
                  else if ((i != 0) && grid_is_cyclic)
                    smooth9_sum((i - 1) * nlon + j + 1 - nlon, mask, 0.3, array1, &avg, &divavg);

                  /* left cell */
                  if (j != 0)
                    smooth9_sum(((i) *nlon) + j - 1, mask, 0.5, array1, &avg, &divavg);
                  else if (grid_is_cyclic)
                    smooth9_sum(i * nlon - 1 + nlon, mask, 0.5, array1, &avg, &divavg);

                  /* right cell */
                  if (j != (nlon - 1))
                    smooth9_sum((i * nlon) + j + 1, mask, 0.5, array1, &avg, &divavg);
                  else if (grid_is_cyclic)
                    smooth9_sum(i * nlon + j + 1 - nlon, mask, 0.5, array1, &avg, &divavg);

                  /* lower left corner */
                  if (mask[ij] && ((i != (nlat - 1)) && (j != 0)))
                    smooth9_sum(((i + 1) * nlon + j - 1), mask, 0.3, array1, &avg, &divavg);
                  else if ((i != (nlat - 1)) && grid_is_cyclic)
                    smooth9_sum((i + 1) * nlon - 1 + nlon, mask, 0.3, array1, &avg, &divavg);

                  /* lower cell */
                  if (i != (nlat - 1)) smooth9_sum(((i + 1) * nlon) + j, mask, 0.5, array1, &avg, &divavg);

                  /* lower right corner */
                  if ((i != (nlat - 1)) && (j != (nlon - 1)))
                    smooth9_sum(((i + 1) * nlon) + j + 1, mask, 0.3, array1, &avg, &divavg);
                  else if ((i != (nlat - 1)) && grid_is_cyclic)
                    smooth9_sum(((i + 1) * nlon) + j + 1 - nlon, mask, 0.3, array1, &avg, &divavg);
                }
            }
          else if (mask[j + nlon * i])
            {
              avg += array1[j + nlon * i];
              divavg += 1;

              smooth9_sum(((i - 1) * nlon) + j - 1, mask, 0.3, array1, &avg, &divavg);
              smooth9_sum(((i - 1) * nlon) + j, mask, 0.5, array1, &avg, &divavg);
              smooth9_sum(((i - 1) * nlon) + j + 1, mask, 0.3, array1, &avg, &divavg);
              smooth9_sum(((i) *nlon) + j - 1, mask, 0.5, array1, &avg, &divavg);
              smooth9_sum((i * nlon) + j + 1, mask, 0.5, array1, &avg, &divavg);
              smooth9_sum(((i + 1) * nlon + j - 1), mask, 0.3, array1, &avg, &divavg);
              smooth9_sum(((i + 1) * nlon) + j, mask, 0.5, array1, &avg, &divavg);
              smooth9_sum(((i + 1) * nlon) + j + 1, mask, 0.3, array1, &avg, &divavg);
            }

          if (fabs(divavg) > 0)
            {
              array2[i * nlon + j] = avg / divavg;
            }
          else
            {
              array2[i * nlon + j] = missval;
              (*nmiss)++;
            }
        }
    }

  Free(mask);
}

double
radius_str_to_deg(const char *string)
{
  char *endptr = NULL;
  double radius = strtod(string, &endptr);

  if (*endptr != 0)
    {
      if (strncmp(endptr, "km", 2) == 0)
        radius = 360 * ((radius * 1000) / (2 * PlanetRadius * M_PI));
      else if (strncmp(endptr, "m", 1) == 0)
        radius = 360 * ((radius) / (2 * PlanetRadius * M_PI));
      else if (strncmp(endptr, "deg", 3) == 0)
        ;
      else if (strncmp(endptr, "rad", 3) == 0)
        radius *= RAD2DEG;
      else
        cdoAbort("Float parameter >%s< contains invalid character at position %d!", string,
                 (int) (endptr - string + 1));
    }

  if (radius > 180.) radius = 180.;

  return radius;
}

static int
convert_form(const char *formstr)
{
  int form = FORM_LINEAR;

  if (strcmp(formstr, "linear") == 0)
    form = FORM_LINEAR;
  else
    cdoAbort("form=%s unsupported!", formstr);

  return form;
}

static void
smoothGetParameter(int *xnsmooth, smoothpoint_t *spoint)
{
  int pargc = operatorArgc();
  if (pargc)
    {
      char **pargv = operatorArgv();

      list_t *kvlist = list_new(sizeof(keyValues_t *), free_keyval, "SMOOTH");
      if (kvlist_parse_cmdline(kvlist, pargc, pargv) != 0) cdoAbort("Parse error!");
      if (cdoVerbose) kvlist_print(kvlist);

      for (listNode_t *kvnode = kvlist->head; kvnode; kvnode = kvnode->next)
        {
          keyValues_t *kv = *(keyValues_t **) kvnode->data;
          const char *key = kv->key;
          if (kv->nvalues > 1) cdoAbort("Too many values for parameter key >%s<!", key);
          if (kv->nvalues < 1) cdoAbort("Missing value for parameter key >%s<!", key);
          const char *value = kv->values[0];

          if (STR_IS_EQ(key, "nsmooth"))
            *xnsmooth = parameter2int(value);
          else if (STR_IS_EQ(key, "maxpoints"))
            spoint->maxpoints = parameter2sizet(value);
          else if (STR_IS_EQ(key, "weight0"))
            spoint->weight0 = parameter2double(value);
          else if (STR_IS_EQ(key, "weightR"))
            spoint->weightR = parameter2double(value);
          else if (STR_IS_EQ(key, "radius"))
            spoint->radius = radius_str_to_deg(value);
          else if (STR_IS_EQ(key, "form"))
            spoint->form = convert_form(value);
          else
            cdoAbort("Invalid parameter key >%s<!", key);
        }

      list_destroy(kvlist);
    }
}

void *
Smooth(void *process)
{
  int nrecs;
  int xnsmooth = 1;
  smoothpoint_t spoint;
  spoint.maxpoints = SIZE_MAX;
  spoint.radius = 1;
  spoint.form = FORM_LINEAR;
  spoint.weight0 = 0.25;
  spoint.weightR = 0.25;

  cdoInitialize(process);

  // clang-format off
  int SMOOTH  = cdoOperatorAdd("smooth",   0,   0, NULL);
  int SMOOTH9 = cdoOperatorAdd("smooth9",  0,   0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  if (operatorID == SMOOTH) smoothGetParameter(&xnsmooth, &spoint);

  if (spoint.radius < 0 || spoint.radius > 180) cdoAbort("%s=%g out of bounds (0-180 deg)!", "radius", spoint.radius);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int nvars = vlistNvars(vlistID1);
  int *varIDs = (int *) Malloc(nvars * sizeof(int));

  for (int varID = 0; varID < nvars; ++varID)
    {
      int gridID = vlistInqVarGrid(vlistID1, varID);
      int gridtype = gridInqType(gridID);
      if (gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT || gridtype == GRID_CURVILINEAR)
        {
          varIDs[varID] = 1;
        }
      else if (gridtype == GRID_UNSTRUCTURED && operatorID == SMOOTH)
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

  size_t gridsizemax = vlistGridsizeMax(vlistID1);
  if (gridsizemax < spoint.maxpoints) spoint.maxpoints = gridsizemax;
  if (cdoVerbose)
    cdoPrint("nsmooth = %d, maxpoints = %zu, radius = %gdeg, form = %s, weight0 = %g, weightR = %g", xnsmooth,
             spoint.maxpoints, spoint.radius, Form[spoint.form], spoint.weight0, spoint.weightR);

  spoint.radius *= DEG2RAD;

  double *array1 = (double *) Malloc(gridsizemax * sizeof(double));
  double *array2 = (double *) Malloc(gridsizemax * sizeof(double));

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          int varID, levelID;
          size_t nmiss;

          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, array1, &nmiss);

          if (varIDs[varID])
            {
              double missval = vlistInqVarMissval(vlistID1, varID);
              int gridID = vlistInqVarGrid(vlistID1, varID);

              for (int i = 0; i < xnsmooth; ++i)
                {
                  if (operatorID == SMOOTH)
                    smooth(gridID, missval, array1, array2, &nmiss, spoint);
                  else if (operatorID == SMOOTH9)
                    smooth9(gridID, missval, array1, array2, &nmiss);

                  arrayCopy(gridsizemax, array2, array1);
                }

              pstreamDefRecord(streamID2, varID, levelID);
              pstreamWriteRecord(streamID2, array2, nmiss);
            }
          else
            {
              pstreamDefRecord(streamID2, varID, levelID);
              pstreamWriteRecord(streamID2, array1, nmiss);
            }
        }

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  Free(varIDs);
  if (array2) Free(array2);
  if (array1) Free(array1);

  cdoFinish();

  return 0;
}
