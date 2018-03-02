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

      Comp       eq              Equal
      Comp       ne              Not equal
      Comp       le              Less equal
      Comp       lt              Less than
      Comp       ge              Greater equal
      Comp       gt              Greater than
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

void *
Comp(void *process)
{
  enum
  {
    FILL_NONE,
    FILL_TS,
    FILL_REC
  };
  int filltype = FILL_NONE;
  int nrecs, nrecs2, nvars = 0, nlev;
  int varID, levelID;
  double missval1, missval2 = 0;
  std::vector<std::vector<double>> vardata;

  cdoInitialize(process);

  int EQ = cdoOperatorAdd("eq", 0, 0, NULL);
  int NE = cdoOperatorAdd("ne", 0, 0, NULL);
  int LE = cdoOperatorAdd("le", 0, 0, NULL);
  int LT = cdoOperatorAdd("lt", 0, 0, NULL);
  int GE = cdoOperatorAdd("ge", 0, 0, NULL);
  int GT = cdoOperatorAdd("gt", 0, 0, NULL);

  int operatorID = cdoOperatorID();

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));
  int streamID2 = cdoStreamOpenRead(cdoStreamName(1));

  int streamIDx1 = streamID1;
  int streamIDx2 = streamID2;

  double *missvalx1 = &missval1;
  double *missvalx2 = &missval2;

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = cdoStreamInqVlist(streamID2);
  int vlistIDx1 = vlistID1;
  int vlistIDx2 = vlistID2;

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = vlistInqTaxis(vlistID2);
  int taxisIDx1 = taxisID1;

  int ntsteps1 = vlistNtsteps(vlistID1);
  int ntsteps2 = vlistNtsteps(vlistID2);
  if (ntsteps1 == 0) ntsteps1 = 1;
  if (ntsteps2 == 0) ntsteps2 = 1;

  bool fillstream1 = false;

  if (vlistNrecs(vlistID1) != 1 && vlistNrecs(vlistID2) == 1)
    {
      filltype = FILL_REC;
      cdoPrint("Filling up stream2 >%s< by copying the first record.", cdoGetStreamName(1).c_str());
      if (ntsteps2 != 1) cdoAbort("stream2 has more than 1 timestep!");
    }
  else if (vlistNrecs(vlistID1) == 1 && vlistNrecs(vlistID2) != 1)
    {
      filltype = FILL_REC;
      cdoPrint("Filling up stream1 >%s< by copying the first record.", cdoGetStreamName(0).c_str());
      if (ntsteps1 != 1) cdoAbort("stream1 has more than 1 timestep!");
      fillstream1 = true;
      streamIDx1 = streamID2;
      streamIDx2 = streamID1;
      vlistIDx1 = vlistID2;
      vlistIDx2 = vlistID1;
      taxisIDx1 = taxisID2;
    }

  if (filltype == FILL_NONE) vlistCompare(vlistID1, vlistID2, CMP_ALL);

  nospec(vlistID1);
  nospec(vlistID2);

  size_t gridsizemax = vlistGridsizeMax(vlistIDx1);

  double *array1 = (double *) Malloc(gridsizemax * sizeof(double));
  double *array2 = (double *) Malloc(gridsizemax * sizeof(double));
  double *array3 = (double *) Malloc(gridsizemax * sizeof(double));

  double *arrayx1 = array1;
  double *arrayx2 = array2;

  if (cdoVerbose) cdoPrint("Number of timesteps: file1 %d, file2 %d", ntsteps1, ntsteps2);

  if (filltype == FILL_NONE)
    {
      if (ntsteps1 != 1 && ntsteps2 == 1)
        {
          filltype = FILL_TS;
          cdoPrint("Filling up stream2 >%s< by copying the first timestep.", cdoGetStreamName(1).c_str());
        }
      else if (ntsteps1 == 1 && ntsteps2 != 1)
        {
          filltype = FILL_TS;
          cdoPrint("Filling up stream1 >%s< by copying the first timestep.", cdoGetStreamName(0).c_str());
          fillstream1 = true;
          streamIDx1 = streamID2;
          streamIDx2 = streamID1;
          vlistIDx1 = vlistID2;
          vlistIDx2 = vlistID1;
          taxisIDx1 = taxisID2;
        }

      if (filltype == FILL_TS)
        {
          nvars = vlistNvars(vlistIDx2);
          vardata.resize(nvars);
          for (varID = 0; varID < nvars; varID++)
            {
              size_t gridsize = gridInqSize(vlistInqVarGrid(vlistIDx2, varID));
              nlev = zaxisInqSize(vlistInqVarZaxis(vlistIDx2, varID));
              vardata[varID].resize(nlev * gridsize);
            }
        }
    }

  if (fillstream1)
    {
      arrayx1 = array2;
      arrayx2 = array1;
      missvalx1 = &missval2;
      missvalx2 = &missval1;
    }

  int vlistID3 = vlistDuplicate(vlistIDx1);

  int taxisID3 = taxisDuplicate(taxisIDx1);
  vlistDefTaxis(vlistID3, taxisID3);

  int streamID3 = cdoStreamOpenWrite(cdoStreamName(2), cdoFiletype());
  pstreamDefVlist(streamID3, vlistID3);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamIDx1, tsID)))
    {
      if (tsID == 0 || filltype == FILL_NONE)
        {
          nrecs2 = cdoStreamInqTimestep(streamIDx2, tsID);
          if (nrecs2 == 0) cdoAbort("Input streams have different number of timesteps!");
        }

      taxisCopyTimestep(taxisID3, taxisIDx1);

      pstreamDefTimestep(streamID3, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          size_t nmiss1;
          pstreamInqRecord(streamIDx1, &varID, &levelID);
          pstreamReadRecord(streamIDx1, arrayx1, &nmiss1);

          if (tsID == 0 || filltype == FILL_NONE)
            {
              if (recID == 0 || filltype != FILL_REC)
                {
                  size_t nmiss2;
                  pstreamInqRecord(streamIDx2, &varID, &levelID);
                  pstreamReadRecord(streamIDx2, arrayx2, &nmiss2);
                }

              if (filltype == FILL_TS)
                {
                  size_t gridsize = gridInqSize(vlistInqVarGrid(vlistIDx2, varID));
                  size_t offset = gridsize * levelID;
                  arrayCopy(gridsize, arrayx2, &vardata[varID][offset]);
                }
            }
          else if (filltype == FILL_TS)
            {
              size_t gridsize = gridInqSize(vlistInqVarGrid(vlistIDx2, varID));
              size_t offset = gridsize * levelID;
              arrayCopy(gridsize, &vardata[varID][offset], arrayx2);
            }

          int datatype1 = vlistInqVarDatatype(vlistIDx1, varID);
          int datatype2 = CDI_UNDEFID;
          size_t gridsize1 = gridInqSize(vlistInqVarGrid(vlistIDx1, varID));
          size_t gridsize2;
          *missvalx1 = vlistInqVarMissval(vlistIDx1, varID);

          if (filltype == FILL_REC)
            {
              datatype2 = vlistInqVarDatatype(vlistIDx2, 0);
              gridsize2 = gridInqSize(vlistInqVarGrid(vlistIDx2, 0));
              *missvalx2 = vlistInqVarMissval(vlistIDx2, 0);
            }
          else
            {
              datatype2 = vlistInqVarDatatype(vlistIDx2, varID);
              gridsize2 = gridInqSize(vlistInqVarGrid(vlistIDx2, varID));
              *missvalx2 = vlistInqVarMissval(vlistIDx2, varID);
            }

          if (gridsize1 != gridsize2)
            cdoAbort("Streams have different gridsize (gridsize1 = %zu; gridsize2 = %zu)!", gridsize1, gridsize2);

          size_t gridsize = gridsize1;

          if (datatype1 != datatype2)
            {
              if (datatype1 == CDI_DATATYPE_FLT32 && datatype2 == CDI_DATATYPE_FLT64)
                {
                  missval2 = (float) missval2;
                  for (size_t i = 0; i < gridsize; i++) array2[i] = (float) array2[i];
                }
              else if (datatype1 == CDI_DATATYPE_FLT64 && datatype2 == CDI_DATATYPE_FLT32)
                {
                  missval1 = (float) missval1;
                  for (size_t i = 0; i < gridsize; i++) array1[i] = (float) array1[i];
                }
            }

          if (operatorID == EQ)
            {
              for (size_t i = 0; i < gridsize; i++)
                array3[i] = (DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2)
                                 ? missval1
                                 : IS_EQUAL(array1[i], array2[i]));
            }
          else if (operatorID == NE)
            {
              for (size_t i = 0; i < gridsize; i++)
                array3[i] = (DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2)
                                 ? missval1
                                 : IS_NOT_EQUAL(array1[i], array2[i]));
            }
          else if (operatorID == LE)
            {
              for (size_t i = 0; i < gridsize; i++)
                array3[i]
                    = (DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2) ? missval1
                                                                                              : array1[i] <= array2[i]);
            }
          else if (operatorID == LT)
            {
              for (size_t i = 0; i < gridsize; i++)
                array3[i]
                    = (DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2) ? missval1
                                                                                              : array1[i] < array2[i]);
            }
          else if (operatorID == GE)
            {
              for (size_t i = 0; i < gridsize; i++)
                array3[i]
                    = (DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2) ? missval1
                                                                                              : array1[i] >= array2[i]);
            }
          else if (operatorID == GT)
            {
              for (size_t i = 0; i < gridsize; i++)
                array3[i]
                    = (DBL_IS_EQUAL(array1[i], missval1) || DBL_IS_EQUAL(array2[i], missval2) ? missval1
                                                                                              : array1[i] > array2[i]);
            }
          else
            {
              cdoAbort("Operator not implemented!");
            }

          size_t nmiss3 = arrayNumMV(gridsize, array3, missval1);
          pstreamDefRecord(streamID3, varID, levelID);
          pstreamWriteRecord(streamID3, array3, nmiss3);
        }

      tsID++;
    }

  pstreamClose(streamID3);
  pstreamClose(streamID2);
  pstreamClose(streamID1);

  if (array3) Free(array3);
  if (array2) Free(array2);
  if (array1) Free(array1);

  cdoFinish();

  return 0;
}
