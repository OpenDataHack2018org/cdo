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
#include <stdio.h>
#include <string.h>

#include <cdi.h>

#include <cdo_int.h>
#include "dmemory.h"
#include "field.h"
#include "util.h"

void
field_init(Field *field)
{
  memset(field, 0, sizeof(Field));
}

Field **
field_allocate(int vlistID, int ptype, int init)
{
  int nvars = vlistNvars(vlistID);

  Field **field = (Field **) Malloc(nvars * sizeof(Field *));

  for (int varID = 0; varID < nvars; ++varID)
    {
      int nwpv = vlistInqNWPV(vlistID, varID);  // number of words per value; real:1  complex:2
      int gridID = vlistInqVarGrid(vlistID, varID);
      size_t gridsize = gridInqSize(gridID);
      int zaxisID = vlistInqVarZaxis(vlistID, varID);
      int nlevel = zaxisInqSize(zaxisID);
      double missval = vlistInqVarMissval(vlistID, varID);

      field[varID] = (Field *) Malloc(nlevel * sizeof(Field));

      for (int levelID = 0; levelID < nlevel; ++levelID)
        {
          field_init(&field[varID][levelID]);

          field[varID][levelID].nwpv = nwpv;
          field[varID][levelID].grid = gridID;
          field[varID][levelID].size = gridsize;
          field[varID][levelID].nsamp = 0;
          field[varID][levelID].nmiss = 0;
          field[varID][levelID].nmiss2 = 0;
          if (ptype & FIELD_FLT) field[varID][levelID].memtype = MEMTYPE_FLOAT;
          field[varID][levelID].missval = missval;
          field[varID][levelID].ptr = NULL;
          field[varID][levelID].ptrf = NULL;
          field[varID][levelID].ptr2 = NULL;
          field[varID][levelID].ptr2f = NULL;
          field[varID][levelID].weight = NULL;

          if (ptype & FIELD_PTR)
            {
              if (ptype & FIELD_FLT)
                {
                  field[varID][levelID].ptrf = (float *) Malloc(nwpv * gridsize * sizeof(float));
                  if (init) arrayFill(nwpv * gridsize, field[varID][levelID].ptrf, 0.0f);
                }
              else
                {
                  field[varID][levelID].ptr = (double *) Malloc(nwpv * gridsize * sizeof(double));
                  if (init) arrayFill(nwpv * gridsize, field[varID][levelID].ptr, 0.0);
                }
            }

          if (ptype & FIELD_PTR2)
            {
              if (ptype & FIELD_FLT)
                {
                  field[varID][levelID].ptr2f = (float *) Malloc(nwpv * gridsize * sizeof(float));
                  if (init) arrayFill(nwpv * gridsize, (float *) field[varID][levelID].ptr2f, 0.0f);
                }
              else
                {
                  field[varID][levelID].ptr2 = (double *) Malloc(nwpv * gridsize * sizeof(double));
                  if (init) arrayFill(nwpv * gridsize, (double *) field[varID][levelID].ptr2, 0.0);
                }
            }

          if (ptype & FIELD_WGT)
            {
              field[varID][levelID].weight = (double *) Malloc(nwpv * gridsize * sizeof(double));
              if (init) arrayFill(nwpv * gridsize, field[varID][levelID].weight, 0.0);
            }
        }
    }

  return field;
}

Field **
field_malloc(int vlistID, int ptype)
{
  return field_allocate(vlistID, ptype, 0);
}

Field **
field_calloc(int vlistID, int ptype)
{
  return field_allocate(vlistID, ptype, 1);
}

void
field_free(Field **field, int vlistID)
{
  int nvars = vlistNvars(vlistID);
  for (int varID = 0; varID < nvars; ++varID)
    {
      int nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
      for (int levelID = 0; levelID < nlevel; ++levelID)
        {
          if (field[varID][levelID].ptr) Free(field[varID][levelID].ptr);
          if (field[varID][levelID].ptrf) Free(field[varID][levelID].ptrf);
          if (field[varID][levelID].ptr2) Free(field[varID][levelID].ptr2);
          if (field[varID][levelID].ptr2f) Free(field[varID][levelID].ptr2f);
          if (field[varID][levelID].weight) Free(field[varID][levelID].weight);
        }

      Free(field[varID]);
    }

  Free(field);
}
