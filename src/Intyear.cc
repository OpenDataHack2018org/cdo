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

      Intyear    intyear         Year interpolation
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "listarray.h"

void *
Intyear(void *process)
{
  int nrecs;
  int varID, levelID;
  size_t nmiss1, nmiss2, nmiss3;
  char filesuffix[32];
  char filename[8192];

  cdoInitialize(process);

  operatorInputArg("years");

  ListArray<int> listArrayInt;
  int nyears = listArrayInt.argvToInt(operatorArgc(), operatorArgv());
  int *iyears = listArrayInt.data();

  std::vector<int> streamIDs(nyears);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));
  int streamID2 = cdoStreamOpenRead(cdoStreamName(1));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = cdoStreamInqVlist(streamID2);
  int vlistID3 = vlistDuplicate(vlistID1);

  vlistCompare(vlistID1, vlistID2, CMP_ALL);

  size_t gridsize = vlistGridsizeMax(vlistID1);
  std::vector<double> array1(gridsize);
  std::vector<double> array2(gridsize);
  std::vector<double> array3(gridsize);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = vlistInqTaxis(vlistID2);
  int taxisID3 = taxisDuplicate(taxisID1);
  if (taxisHasBounds(taxisID3)) taxisDeleteBounds(taxisID3);
  vlistDefTaxis(vlistID3, taxisID3);

  strcpy(filename, cdoGetStreamName(2).c_str());
  int nchars = strlen(filename);

  const char *refname = cdoGetObase();
  filesuffix[0] = 0;
  cdoGenFileSuffix(filesuffix, sizeof(filesuffix), pstreamInqFiletype(streamID1), vlistID1, refname);

  for (int iy = 0; iy < nyears; iy++)
    {
      sprintf(filename + nchars, "%04d", iyears[iy]);
      if (filesuffix[0]) sprintf(filename + nchars + 4, "%s", filesuffix);

      streamIDs[iy] = cdoStreamOpenWrite(filename, cdoFiletype());
      pstreamDefVlist(streamIDs[iy], vlistID3);
    }

  int tsID = 0;
  while (TRUE)
    {
      nrecs = cdoStreamInqTimestep(streamID1, tsID);
      if (nrecs == 0) break;
      nrecs = cdoStreamInqTimestep(streamID2, tsID);
      if (nrecs == 0) cdoAbort("Too few timesteps in second inputfile!");

      int vtime = taxisInqVtime(taxisID1);
      int vdate1 = taxisInqVdate(taxisID1);
      int year1 = vdate1 / 10000;
      int vdate2 = taxisInqVdate(taxisID2);
      int year2 = vdate2 / 10000;

      for (int iy = 0; iy < nyears; iy++)
        {
          if (iyears[iy] < year1 || iyears[iy] > year2)
            cdoAbort("Year %d out of bounds (first year %d; last year %d)!", iyears[iy], year1, year2);
          int vdate3 = vdate1 - year1 * 10000 + iyears[iy] * 10000;
          taxisDefVdate(taxisID3, vdate3);
          taxisDefVtime(taxisID3, vtime);
          pstreamDefTimestep(streamIDs[iy], tsID);
        }

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamInqRecord(streamID2, &varID, &levelID);

          pstreamReadRecord(streamID1, array1.data(), &nmiss1);
          pstreamReadRecord(streamID2, array2.data(), &nmiss2);

          gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));

          for (int iy = 0; iy < nyears; iy++)
            {
              double fac1 = ((double) year2 - iyears[iy]) / (year2 - year1);
              double fac2 = ((double) iyears[iy] - year1) / (year2 - year1);

              nmiss3 = 0;

              if (nmiss1 > 0 || nmiss2 > 0)
                {
                  double missval1 = vlistInqVarMissval(vlistID1, varID);
                  double missval2 = vlistInqVarMissval(vlistID2, varID);

                  for (size_t i = 0; i < gridsize; i++)
                    {
                      if (!DBL_IS_EQUAL(array1[i], missval1) && !DBL_IS_EQUAL(array2[i], missval2))
                        array3[i] = array1[i] * fac1 + array2[i] * fac2;
                      /* 2010-04-19 Uwe Schulzweida: removed
                      else if ( DBL_IS_EQUAL(array1[i], missval1) &&
                                !DBL_IS_EQUAL(array2[i], missval2) && fac2 >=
                      0.5 ) array3[i] = array2[i]; else if (
                      DBL_IS_EQUAL(array2[i], missval2) &&
                                !DBL_IS_EQUAL(array1[i], missval1) && fac1 >=
                      0.5 ) array3[i] = array1[i];
                      */
                      else
                        {
                          array3[i] = missval1;
                          nmiss3++;
                        }
                    }
                }
              else
                {
                  for (size_t i = 0; i < gridsize; i++) array3[i] = array1[i] * fac1 + array2[i] * fac2;
                }

              pstreamDefRecord(streamIDs[iy], varID, levelID);
              pstreamWriteRecord(streamIDs[iy], array3.data(), nmiss3);
            }
        }

      tsID++;
    }

  for (int iy = 0; iy < nyears; iy++) pstreamClose(streamIDs[iy]);

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}
