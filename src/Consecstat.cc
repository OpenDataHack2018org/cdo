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

      Consectstep  consecsum  For each timestep, the current number of
                              onsecutive timsteps is counted
      Consectstep  consects   For each period of consecutive timesteps, only
                              count its lenght + last contributing timesteps

   =============================================================================
   Created:  04/08/2010 11:58:01 AM
    Author:  Ralf Mueller (ram), ralf.mueller@mpimet.mpg.de
   Company:  Max-Planck-Institute for Meteorology
   =============================================================================
 */

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

enum
{
  CONSECSUM,
  CONSECTS
};
#define SWITCHWARN "Hit default case!This should never happen (%s).\n"

static void
selEndOfPeriod(Field *periods, Field history, Field current, int isLastTimestep)
{
  size_t i;
  double pmissval = periods->missval;
  double *parray = periods->ptr;
  double hmissval = history.missval;
  double *harray = history.ptr;
  double cmissval = current.missval;
  double *carray = current.ptr;

  size_t len = gridInqSize(periods->grid);
  if (len != gridInqSize(current.grid) || (gridInqSize(current.grid) != gridInqSize(history.grid)))
    cdoAbort("Fields have different gridsize (%s)", __func__);

  if (!isLastTimestep)
    {
      if (current.nmiss > 0 || history.nmiss > 0)
        {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
          for (i = 0; i < len; i++)
            {
              if (!DBL_IS_EQUAL(harray[i], hmissval))
                {
                  if (!DBL_IS_EQUAL(carray[i], cmissval))
                    {
                      parray[i] = (DBL_IS_EQUAL(carray[i], 0.0) && IS_NOT_EQUAL(harray[i], 0.0)) ? harray[i] : pmissval;
                    }
                  else /* DBL_IS_EQUAL(carray[i], cmissval) */
                    {
                      parray[i] = (IS_NOT_EQUAL(harray[i], 0.0)) ? harray[i] : pmissval;
                    }
                }
              else /* DBL_IS_EQUAL(harray[i], hmissval) */
                {
                  parray[i] = pmissval;
                }
            }
        }
      else
        {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
          for (i = 0; i < len; i++)
            parray[i] = (DBL_IS_EQUAL(carray[i], 0.0) && IS_NOT_EQUAL(harray[i], 0.0)) ? harray[i] : pmissval;
        }
    }
  else
    {
      if (current.nmiss > 0)
        {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
          for (i = 0; i < len; i++)
            if (!DBL_IS_EQUAL(carray[i], cmissval))
              {
                parray[i] = (DBL_IS_EQUAL(carray[i], 0.0)) ? pmissval : carray[i];
              }
            else /* DBL_IS_EQUAL(carray[i], cmissval) */
              {
                parray[i] = pmissval;
              }
        }
      else
        {
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
          for (i = 0; i < len; i++) parray[i] = DBL_IS_EQUAL(carray[i], 0.0) ? pmissval : carray[i];
        }
    }

  periods->nmiss = arrayNumMV(len, parray, pmissval);
}

void *
Consecstat(void *process)
{
  int64_t vdate = 0, histvdate = 0;
  int vtime = 0, histvtime = 0;
  int nrecs;
  int varID;
  int levelID, nlevels;
  size_t nmiss;
  double refval = 0.0;

  cdoInitialize(process);
  cdoOperatorAdd("consecsum", CONSECSUM, 0, "refval");
  cdoOperatorAdd("consects", CONSECTS, 0, NULL);
  int operatorID = cdoOperatorID();

  if (operatorID == CONSECSUM)
    if (operatorArgc() > 0) refval = parameter2double(operatorArgv()[0]);

  int istreamID = cdoStreamOpenRead(cdoStreamName(0));

  int ivlistID = cdoStreamInqVlist(istreamID);
  int itaxisID = vlistInqTaxis(ivlistID);
  int ovlistID = vlistDuplicate(ivlistID);
  int otaxisID = taxisDuplicate(itaxisID);
  vlistDefTaxis(ovlistID, otaxisID);

  Field field;
  field_init(&field);
  field.ptr = (double *) Malloc(vlistGridsizeMax(ovlistID) * sizeof(double));

  int nvars = vlistNvars(ivlistID);
  Field **vars = field_calloc(ivlistID, FIELD_PTR);
  Field **hist = NULL, **periods = NULL;
  if (operatorID == CONSECTS)
    {
      hist = field_malloc(ivlistID, FIELD_PTR);
      periods = field_malloc(ivlistID, FIELD_PTR);
    }

  for (varID = 0; varID < nvars; varID++)
    {
      vlistDefVarUnits(ovlistID, varID, "steps"); /* TODO */
    }

  int ostreamID = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(ostreamID, ovlistID);

  int itsID = 0;
  int otsID = 0;
  while ((nrecs = cdoStreamInqTimestep(istreamID, itsID)))
    {
      vdate = taxisInqVdate(itaxisID);
      vtime = taxisInqVtime(itaxisID);
      switch (operatorID)
        {
        case CONSECSUM:
          taxisDefVdate(otaxisID, vdate);
          taxisDefVtime(otaxisID, vtime);
          pstreamDefTimestep(ostreamID, otsID);
          break;
        case CONSECTS:
          if (itsID != 0)
            {
              taxisDefVdate(otaxisID, histvdate);
              taxisDefVtime(otaxisID, histvtime);
              pstreamDefTimestep(ostreamID, otsID - 1);
            }
          break;
        default: printf(SWITCHWARN, __func__); break;
        }

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(istreamID, &varID, &levelID);
          pstreamReadRecord(istreamID, field.ptr, &nmiss);
          field.nmiss = nmiss;
          field.grid = vlistInqVarGrid(ovlistID, varID);
          field.missval = vlistInqVarMissval(ovlistID, varID);

          farsumtr(&vars[varID][levelID], field, refval);

          switch (operatorID)
            {
            case CONSECSUM:
              pstreamDefRecord(ostreamID, varID, levelID);
              pstreamWriteRecord(ostreamID, vars[varID][levelID].ptr, vars[varID][levelID].nmiss);
              break;
            case CONSECTS:
              if (itsID != 0)
                {
                  selEndOfPeriod(&periods[varID][levelID], hist[varID][levelID], vars[varID][levelID], FALSE);
                  pstreamDefRecord(ostreamID, varID, levelID);
                  pstreamWriteRecord(ostreamID, periods[varID][levelID].ptr, periods[varID][levelID].nmiss);
                }
#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
              for (size_t i = 0; i < gridInqSize(vars[varID][levelID].grid); i++)
                hist[varID][levelID].ptr[i] = vars[varID][levelID].ptr[i];
#else
              memcpy(hist[varID][levelID].ptr, vars[varID][levelID].ptr, gridInqSize(vars[varID][levelID].grid) * sizeof(double));
#endif
              break;
            default: printf(SWITCHWARN, __func__); break;
            }
        }
      histvdate = vdate;
      histvtime = vtime;
      itsID++;
      otsID++;
    }

  if (operatorID == CONSECTS) /* Save the last timestep */
    {
      taxisDefVdate(otaxisID, vdate);
      taxisDefVtime(otaxisID, vtime);
      pstreamDefTimestep(ostreamID, otsID - 1);
      for (varID = 0; varID < nvars; varID++)
        {
          nlevels = zaxisInqSize(vlistInqVarZaxis(ovlistID, varID));
          for (levelID = 0; levelID < nlevels; levelID++)
            {
              selEndOfPeriod(&periods[varID][levelID], hist[varID][levelID], vars[varID][levelID], TRUE);
              pstreamDefRecord(ostreamID, varID, levelID);
              pstreamWriteRecord(ostreamID, periods[varID][levelID].ptr, periods[varID][levelID].nmiss);
            }
        }
    }

  if (vars) field_free(vars, ivlistID);
  if (hist) field_free(hist, ivlistID);
  if (periods) field_free(periods, ivlistID);

  pstreamClose(istreamID);
  pstreamClose(ostreamID);

  cdoFinish();

  return 0;
}
