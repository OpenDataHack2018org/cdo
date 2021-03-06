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

      Setmiss    setmissval      Set a new missing value
      Setmiss    setctomiss      Set constant to missing value
      Setmiss    setmisstoc      Set missing value to constant
      Setmiss    setrtomiss      Set range to missing value
      Setmiss    setvrange       Set range of valid value
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cmath>
#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

void *
Setmiss(void *process)
{
  int nrecs;
  int varID, levelID;
  size_t nmiss;
  double missval2 = 0;
  double rconst = 0, rmin = 0, rmax = 0;

  cdoInitialize(process);

  // clang-format off
  int SETMISSVAL = cdoOperatorAdd("setmissval", 0, 0, "missing value");
  int SETCTOMISS = cdoOperatorAdd("setctomiss", 0, 0, "constant");
  int SETMISSTOC = cdoOperatorAdd("setmisstoc", 0, 0, "constant");
  int SETRTOMISS = cdoOperatorAdd("setrtomiss", 0, 0, "range (min, max)");
  int SETVRANGE  = cdoOperatorAdd("setvrange",  0, 0, "range (min, max)");
  // clang-format on

  int operatorID = cdoOperatorID();

  if (operatorID == SETMISSVAL)
    {
      operatorCheckArgc(1);
      missval2 = parameter2double(operatorArgv()[0]);
    }
  else if (operatorID == SETCTOMISS || operatorID == SETMISSTOC)
    {
      operatorCheckArgc(1);
      /*
      if ( operatorArgv()[0][0] == 'n' || operatorArgv()[0][0] == 'N' )
        {
#if ! defined(HAVE_ISNAN)
          cdoWarning("Function >isnan< not available!");
#endif
          rconst = 0.0/0.0;
        }
      else
      */
      rconst = parameter2double(operatorArgv()[0]);
    }
  else
    {
      operatorCheckArgc(2);
      rmin = parameter2double(operatorArgv()[0]);
      rmax = parameter2double(operatorArgv()[1]);
    }

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if (operatorID == SETMISSVAL)
    {
      int nvars = vlistNvars(vlistID2);
      for (varID = 0; varID < nvars; varID++) vlistDefVarMissval(vlistID2, varID, missval2);
    }
  else if (operatorID == SETMISSTOC)
    {
      int nvars = vlistNvars(vlistID2);
      for (varID = 0; varID < nvars; varID++)
        {
          double missval = vlistInqVarMissval(vlistID2, varID);
          if (DBL_IS_EQUAL(rconst, missval))
            {
              cdoWarning("Missing value and constant have the same value!");
              break;
            }
        }
    }

  /*
  if ( operatorID == SETVRANGE )
    {
      double range[2];
      range[0] = rmin;
      range[1] = rmax;

      nvars = vlistNvars(vlistID2);
      for ( varID = 0; varID < nvars; varID++ )
        cdiDefAttFlt(vlistID2, varID, "valid_range", CDI_DATATYPE_FLT64, 2,
  range);
    }
  */
  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  size_t gridsizemax = vlistGridsizeMax(vlistID1);
  double *array = (double *) Malloc(gridsizemax * sizeof(double));

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, array, &nmiss);

          size_t gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
          double missval = vlistInqVarMissval(vlistID1, varID);

          if (operatorID == SETMISSVAL)
            {
              nmiss = 0;
              for (size_t i = 0; i < gridsize; i++)
                if (DBL_IS_EQUAL(array[i], missval) || DBL_IS_EQUAL(array[i], (float) missval) || DBL_IS_EQUAL(array[i], missval2)
                    || DBL_IS_EQUAL(array[i], (float) missval2))
                  {
                    array[i] = missval2;
                    nmiss++;
                  }
            }
          else if (operatorID == SETCTOMISS)
            {
              if (DBL_IS_NAN(rconst))
                {
                  for (size_t i = 0; i < gridsize; i++)
                    if (DBL_IS_NAN(array[i]))
                      {
                        array[i] = missval;
                        nmiss++;
                      }
                }
              else
                {
                  for (size_t i = 0; i < gridsize; i++)
                    if (DBL_IS_EQUAL(array[i], rconst) || DBL_IS_EQUAL(array[i], (float) rconst))
                      {
                        array[i] = missval;
                        nmiss++;
                      }
                }
            }
          else if (operatorID == SETMISSTOC)
            {
              nmiss = 0;
              for (size_t i = 0; i < gridsize; i++)
                if (DBL_IS_EQUAL(array[i], missval) || DBL_IS_EQUAL(array[i], (float) missval))
                  {
                    array[i] = rconst;
                  }
            }
          else if (operatorID == SETRTOMISS)
            {
              for (size_t i = 0; i < gridsize; i++)
                if (array[i] >= rmin && array[i] <= rmax)
                  {
                    array[i] = missval;
                    nmiss++;
                  }
            }
          else if (operatorID == SETVRANGE)
            {
              for (size_t i = 0; i < gridsize; i++)
                if (array[i] < rmin || array[i] > rmax) array[i] = missval;

              nmiss = arrayNumMV(gridsize, array, missval);
            }

          pstreamDefRecord(streamID2, varID, levelID);
          pstreamWriteRecord(streamID2, array, nmiss);
        }

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if (array) Free(array);

  cdoFinish();

  return 0;
}
