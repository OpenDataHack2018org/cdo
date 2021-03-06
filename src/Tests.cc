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

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "statistic.h"

void *
Tests(void *process)
{
  int nrecs;
  int varID, levelID;
  size_t nmiss;
  double degree_of_freedom = 0, p = 0, q = 0, n = 0, d = 0;
  double missval;

  cdoInitialize(process);

  // clang-format off
  int NORMAL    = cdoOperatorAdd("normal",    0, 0, NULL);
  int STUDENTT  = cdoOperatorAdd("studentt",  0, 0, "degree of freedom");
  int CHISQUARE = cdoOperatorAdd("chisquare", 0, 0, "degree of freedom");
  int BETA      = cdoOperatorAdd("beta",      0, 0, "p and q");
  int FISHER    = cdoOperatorAdd("fisher",    0, 0, "degree of freedom of nominator and of denominator");
  // clang-format on

  int operatorID = cdoOperatorID();

  if (operatorID == STUDENTT || operatorID == CHISQUARE)
    {
      operatorInputArg(cdoOperatorEnter(operatorID));

      operatorCheckArgc(1);

      degree_of_freedom = parameter2double(operatorArgv()[0]);

      if (degree_of_freedom <= 0) cdoAbort("degree of freedom must be positive!");
    }
  else if (operatorID == BETA)
    {
      operatorInputArg(cdoOperatorEnter(operatorID));

      operatorCheckArgc(2);

      p = parameter2double(operatorArgv()[0]);
      q = parameter2double(operatorArgv()[1]);

      if (p <= 0 || q <= 0) cdoAbort("p and q must be positive!");
    }
  else if (operatorID == FISHER)
    {
      operatorInputArg(cdoOperatorEnter(operatorID));

      operatorCheckArgc(2);

      n = parameter2double(operatorArgv()[0]);
      d = parameter2double(operatorArgv()[1]);

      if (n <= 0 || d <= 0) cdoAbort("both degrees must be positive!");
    }

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  size_t gridsize = vlistGridsizeMax(vlistID1);
  double *array1 = (double *) Malloc(gridsize * sizeof(double));
  double *array2 = (double *) Malloc(gridsize * sizeof(double));

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, array1, &nmiss);

          gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
          missval = vlistInqVarMissval(vlistID1, varID);

          if (operatorID == NORMAL)
            {
              for (size_t i = 0; i < gridsize; i++)
                array2[i] = DBL_IS_EQUAL(array1[i], missval) ? missval : normal(array1[i], processInqPrompt());
            }
          else if (operatorID == STUDENTT)
            {
              for (size_t i = 0; i < gridsize; i++)
                array2[i]
                    = DBL_IS_EQUAL(array1[i], missval) ? missval : student_t(degree_of_freedom, array1[i], processInqPrompt());
            }
          else if (operatorID == CHISQUARE)
            {
              for (size_t i = 0; i < gridsize; i++)
                array2[i]
                    = DBL_IS_EQUAL(array1[i], missval) ? missval : chi_square(degree_of_freedom, array1[i], processInqPrompt());
            }
          else if (operatorID == BETA)
            {
              for (size_t i = 0; i < gridsize; i++)
                {
                  if (array1[i] < 0 || array1[i] > 1) cdoAbort("Value out of range (0-1)!");

                  array2[i] = DBL_IS_EQUAL(array1[i], missval) ? missval : beta_distr(p, q, array1[i], processInqPrompt());
                }
            }
          else if (operatorID == FISHER)
            {
              for (size_t i = 0; i < gridsize; i++)
                array2[i] = DBL_IS_EQUAL(array1[i], missval) ? missval : fisher(n, d, array1[i], processInqPrompt());
            }
          else
            {
              cdoAbort("Internal problem, operator not implemented!");
            }

          pstreamDefRecord(streamID2, varID, levelID);
          pstreamWriteRecord(streamID2, array2, nmiss);
        }

      tsID++;
    }

  pstreamClose(streamID1);
  pstreamClose(streamID2);

  vlistDestroy(vlistID2);

  if (array1) Free(array1);
  if (array2) Free(array2);

  cdoFinish();

  return 0;
}
