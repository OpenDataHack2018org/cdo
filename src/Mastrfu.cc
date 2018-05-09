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

      Mastrfu    mastrfu         Mass stream function
*/

#include <cdi.h>

#include "cdo_int.h"
#include "grid.h"
#include "pstream_int.h"


static void
mastrfu(int gridID, int zaxisID, std::vector<std::vector<double>> &field1, std::vector<std::vector<double>> &field2, size_t nmiss, double missval)
{
  size_t ilat;
  int ilev, n;
  double fact = 4 * atan(1.0) * 6371000 / 9.81;
  char units[CDI_MAX_NAME];

  size_t nlat = gridInqSize(gridID);
  int nlev = zaxisInqSize(zaxisID);
  std::vector<double> phi(nlat);
  std::vector<double> cosphi(nlat);
  std::vector<double> plevel(nlev);

  cdoZaxisInqLevels(zaxisID, plevel.data());

  gridInqYvals(gridID, phi.data());
  gridInqYunits(gridID, units);

  if (memcmp(units, "degree", 6) == 0)
    for (ilat = 0; ilat < nlat; ilat++) phi[ilat] *= DEG2RAD;

  for (ilat = 0; ilat < nlat; ilat++) phi[ilat] = sin(phi[ilat]);

  for (ilat = 0; ilat < nlat; ilat++) cosphi[ilat] = sqrt(1.0 - phi[ilat] * phi[ilat]);

  for (ilev = 0; ilev < nlev; ilev++)
    for (ilat = 0; ilat < nlat; ilat++) field2[ilev][ilat] = 0.0;

  if (nmiss == 0)
    {
      for (ilev = nlev - 1; ilev >= 0; ilev--)
        for (n = ilev; n < nlev - 1; n++)
          for (ilat = 0; ilat < nlat; ilat++)
            {
              field2[ilev][ilat] += fact * (field1[n][ilat] + field1[n + 1][ilat]) * cosphi[ilat] * (plevel[n] - plevel[n + 1]);
            }
    }
  else
    {
      for (ilat = 0; ilat < nlat; ilat++)
        for (ilev = nlev - 1; ilev >= 0; ilev--)
          for (n = ilev; n < nlev - 1; n++)
            {
              if (DBL_IS_EQUAL(field1[n][ilat], missval))
                {
                  field2[ilev][ilat] = missval;
                  break;
                }
              else
                field2[ilev][ilat]
                    += fact * (field1[n][ilat] + field1[n + 1][ilat]) * cosphi[ilat] * (plevel[n] - plevel[n + 1]);
            }
    }
}

void *
Mastrfu(void *process)
{
  int nrecs;
  int varID, levelID;
  size_t nmiss, nmiss1;

  cdoInitialize(process);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);

  int nvars = vlistNvars(vlistID1);
  if (nvars != 1) cdoAbort("This operator works only with one variable!");

  int code = vlistInqVarCode(vlistID1, 0);
  if (code > 0 && code != 132) cdoWarning("Unexpected code %d!", code);

  double missval = vlistInqVarMissval(vlistID1, 0);

  int zaxisID = vlistInqVarZaxis(vlistID1, 0);
  if (zaxisInqType(zaxisID) != ZAXIS_PRESSURE && zaxisInqType(zaxisID) != ZAXIS_GENERIC)
    {
      char longname[CDI_MAX_NAME];
      zaxisInqLongname(zaxisID, longname);
      cdoWarning("Unexpected vertical grid %s!", longname);
    }

  int gridID = vlistInqVarGrid(vlistID1, 0);
  if (gridInqXsize(gridID) > 1) cdoAbort("Grid must be a zonal mean!");

  size_t gridsize = gridInqSize(gridID);
  int nlev = zaxisInqSize(zaxisID);

  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  vlistDefVarCode(vlistID2, 0, 272);
  vlistDefVarName(vlistID2, 0, "mastrfu");
  vlistDefVarLongname(vlistID2, 0, "mass stream function");
  vlistDefVarUnits(vlistID2, 0, "kg/s");
  vlistDefVarDatatype(vlistID2, 0, CDI_DATATYPE_FLT32);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  VECTOR_2D(double, array1, nlev, gridsize);
  VECTOR_2D(double, array2, nlev, gridsize);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      nmiss = 0;
      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamReadRecord(streamID1, array1[levelID].data(), &nmiss1);
          nmiss += nmiss1;
        }

      mastrfu(gridID, zaxisID, array1, array2, nmiss, missval);

      for (int recID = 0; recID < nrecs; recID++)
        {
          varID = 0;
          levelID = recID;
          pstreamDefRecord(streamID2, varID, levelID);
          pstreamWriteRecord(streamID2, array2[levelID].data(), nmiss);
        }

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}
