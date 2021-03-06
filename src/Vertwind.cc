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

      Vertwind    vertwind      Convert the vertical velocity to [m/s]
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "after_vertint.h"
#include "util_string.h"

#define R 287.07  /* spezielle Gaskonstante fuer Luft */
#define G 9.80665 /* Erdbeschleunigung */

void *
Vertwind(void *process)
{
  int nrecs;
  int varID, levelID;
  int nvct = 0;
  size_t nmiss;
  int tempID = -1, sqID = -1, psID = -1, omegaID = -1;
  char varname[CDI_MAX_NAME];
  std::vector<double> vct;
  std::vector<double> hpress, ps_prog;

  cdoInitialize(process);

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);

  vlist_check_gridsize(vlistID1);

  int temp_code = 130;
  int sq_code = 133;
  int ps_code = 134;
  int omega_code = 135;

  int nvars = vlistNvars(vlistID1);
  for (varID = 0; varID < nvars; ++varID)
    {
      int code = vlistInqVarCode(vlistID1, varID);

      if (code <= 0)
        {
          vlistInqVarName(vlistID1, varID, varname);
          strtolower(varname);

          // clang-format off
          if      (strcmp(varname, "st") == 0) code = temp_code;
          else if (strcmp(varname, "sq") == 0) code = sq_code;
          else if (strcmp(varname, "aps") == 0) code = ps_code;
          else if (strcmp(varname, "omega") == 0) code = omega_code;
          // clang-format on
        }

      // clang-format off
      if      (code == temp_code)  tempID = varID;
      else if (code == sq_code)    sqID = varID;
      else if (code == ps_code)    psID = varID;
      else if (code == omega_code) omegaID = varID;
      // clang-format on
    }

  if (tempID == -1 || sqID == -1 || omegaID == -1)
    {
      if (tempID == -1) cdoWarning("Temperature (code 130) not found!");
      if (sqID == -1) cdoWarning("Specific humidity (code 133) not found!");
      if (omegaID == -1) cdoWarning("Vertical velocity (code 135) not found!");
      cdoAbort("Parameter not found!");
    }

  /* Get missing values */
  double missval_t = vlistInqVarMissval(vlistID1, tempID);
  double missval_sq = vlistInqVarMissval(vlistID1, sqID);
  double missval_wap = vlistInqVarMissval(vlistID1, omegaID);
  double missval_out = missval_wap;

  int gridID = vlistInqVarGrid(vlistID1, omegaID);
  int zaxisID = vlistInqVarZaxis(vlistID1, omegaID);

  if (psID == -1 && zaxisInqType(zaxisID) == ZAXIS_HYBRID) cdoAbort("Surface pressure (code 134) not found!");

  size_t gridsize = gridInqSize(gridID);
  int nlevel = zaxisInqSize(zaxisID);
  std::vector<double> level(nlevel);
  cdoZaxisInqLevels(zaxisID, level.data());

  std::vector<double> temp(gridsize * nlevel);
  std::vector<double> sq(gridsize * nlevel);
  std::vector<double> omega(gridsize * nlevel);
  std::vector<double> wms(gridsize * nlevel);
  std::vector<double> fpress(gridsize * nlevel);

  if (zaxisInqType(zaxisID) == ZAXIS_PRESSURE)
    {
      for (levelID = 0; levelID < nlevel; ++levelID)
        {
          size_t offset = (size_t) levelID * gridsize;
          for (size_t i = 0; i < gridsize; ++i) fpress[offset + i] = level[levelID];
        }
    }
  else if (zaxisInqType(zaxisID) == ZAXIS_HYBRID)
    {
      ps_prog.resize(gridsize);
      hpress.resize(gridsize * (nlevel + 1));

      nvct = zaxisInqVctSize(zaxisID);
      if (nlevel == (nvct / 2 - 1))
        {
          vct.resize(nvct);
          zaxisInqVct(zaxisID, vct.data());
        }
      else
        cdoAbort("Unsupported vertical coordinate table format!");
    }
  else
    cdoAbort("Unsupported Z-Axis type!");

  vlistClearFlag(vlistID1);
  for (levelID = 0; levelID < nlevel; ++levelID) vlistDefFlag(vlistID1, omegaID, levelID, TRUE);

  int vlistID2 = vlistCreate();
  cdoVlistCopyFlag(vlistID2, vlistID1);
  vlistDefVarCode(vlistID2, 0, 40);
  vlistDefVarName(vlistID2, 0, "W");
  vlistDefVarLongname(vlistID2, 0, "Vertical velocity");
  vlistDefVarUnits(vlistID2, 0, "m/s");
  vlistDefVarMissval(vlistID2, 0, missval_out);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID2, vlistID2);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);

          size_t offset = (size_t) levelID * gridsize;

          if (varID == tempID)
            pstreamReadRecord(streamID1, &temp[offset], &nmiss);
          else if (varID == sqID)
            pstreamReadRecord(streamID1, &sq[offset], &nmiss);
          else if (varID == omegaID)
            pstreamReadRecord(streamID1, &omega[offset], &nmiss);
          else if (varID == psID && zaxisInqType(zaxisID) == ZAXIS_HYBRID)
            pstreamReadRecord(streamID1, ps_prog.data(), &nmiss);
        }

      if (zaxisInqType(zaxisID) == ZAXIS_HYBRID) presh(fpress.data(), hpress.data(), vct.data(), ps_prog.data(), nlevel, gridsize);

      for (levelID = 0; levelID < nlevel; ++levelID)
        {
          size_t offset = (size_t) levelID * gridsize;

          for (size_t i = 0; i < gridsize; ++i)
            {
              if (DBL_IS_EQUAL(temp[offset + i], missval_t) || DBL_IS_EQUAL(omega[offset + i], missval_wap)
                  || DBL_IS_EQUAL(sq[offset + i], missval_sq))
                {
                  wms[offset + i] = missval_out;
                }
              else
                {
                  // Virtuelle Temperatur bringt die Feuchteabhaengigkeit hinein
                  double tv = temp[offset + i] * (1. + 0.608 * sq[offset + i]);

                  // Die Dichte erhaelt man nun mit der Gasgleichung rho=p/(R*tv) Level in Pa!
                  double rho = fpress[offset + i] / (R * tv);
                  /*
                    Nun daraus die Vertikalgeschwindigkeit im m/s, indem man die Vertikalgeschwindigkeit
                    in Pa/s durch die Erdbeschleunigung und die Dichte teilt
                  */
                  wms[offset + i] = omega[offset + i] / (G * rho);
                }
            }
        }

      for (levelID = 0; levelID < nlevel; ++levelID)
        {
          size_t offset = (size_t) levelID * gridsize;

          size_t nmiss_out = 0;
          for (size_t i = 0; i < gridsize; i++)
            if (DBL_IS_EQUAL(wms[offset + i], missval_out)) nmiss_out++;

          pstreamDefRecord(streamID2, 0, levelID);
          pstreamWriteRecord(streamID2, &wms[offset], nmiss_out);
        }

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  vlistDestroy(vlistID2);

  cdoFinish();

  return 0;
}
