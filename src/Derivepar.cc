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

      Derivepar     gheight          geopotential height
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "after_vertint.h"
#include "stdnametable.h"
#include "util_string.h"

void MakeGeopotHeight(double *geop, double *gt, double *gq, double *ph, int nhor, int nlev);

double *vlist_hybrid_vct(int vlistID, int *rzaxisIDh, int *rnvct, int *rnhlevf);

void *
Derivepar(void *process)
{
  ModelMode mode(ModelMode::UNDEF);
  int nrecs;
  int i, offset;
  int varID, levelID;
  int zaxisID;
  int nlevel;
  int surfaceID = -1;
  int sgeopotID = -1, geopotID = -1, tempID = -1, humID = -1, psID = -1, lnpsID = -1, presID = -1, gheightID = -1;
  // int clwcID = -1, ciwcID = -1;
  int code, param;
  int pnum, pcat, pdis;
  char paramstr[32];
  char varname[CDI_MAX_NAME], stdname[CDI_MAX_NAME];
  double *single2;
  // double *lwater = NULL, *iwater = NULL;
  size_t nmiss, nmissout = 0;
  double minval, maxval;
  int instNum, tableNum;
  gribcode_t gribcodes;
  memset(&gribcodes, 0, sizeof(gribcode_t));

  cdoInitialize(process);

  // clang-format off
  int GHEIGHT          = cdoOperatorAdd("gheight",            0, 0, NULL);
  int SEALEVELPRESSURE = cdoOperatorAdd("sealevelpressure",   0, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);

  int gridID = vlistGrid(vlistID1, 0);
  if (gridInqType(gridID) == GRID_SPECTRAL) cdoAbort("Spectral data unsupported!");

  size_t gridsize = vlist_check_gridsize(vlistID1);

  int zaxisIDh = -1;
  int nvct = 0;
  int nhlevf = 0;
  double *vct = vlist_hybrid_vct(vlistID1, &zaxisIDh, &nvct, &nhlevf);

  if (cdoVerbose)
    for (i = 0; i < nvct / 2; ++i) cdoPrint("vct: %5d %25.17f %25.17f", i, vct[i], vct[nvct / 2 + i]);

  if (zaxisIDh == -1) cdoAbort("No 3D variable with hybrid sigma pressure coordinate found!");

  int nvars = vlistNvars(vlistID1);

  bool useTable = false;
  for (varID = 0; varID < nvars; varID++)
    {
      tableNum = tableInqNum(vlistInqVarTable(vlistID1, varID));

      if (tableNum > 0 && tableNum != 255)
        {
          useTable = true;
          break;
        }
    }

  if (cdoVerbose && useTable) cdoPrint("Using code tables!");

  for (varID = 0; varID < nvars; varID++)
    {
      gridID = vlistInqVarGrid(vlistID1, varID);
      zaxisID = vlistInqVarZaxis(vlistID1, varID);
      nlevel = zaxisInqSize(zaxisID);
      instNum = institutInqCenter(vlistInqVarInstitut(vlistID1, varID));
      tableNum = tableInqNum(vlistInqVarTable(vlistID1, varID));

      code = vlistInqVarCode(vlistID1, varID);
      param = vlistInqVarParam(vlistID1, varID);

      cdiParamToString(param, paramstr, sizeof(paramstr));
      cdiDecodeParam(param, &pnum, &pcat, &pdis);
      if (pdis >= 0 && pdis < 255) code = -1;

      if (useTable)
        {
          if (tableNum == 2)
            {
              mode = ModelMode::WMO;
              wmo_gribcodes(&gribcodes);
            }
          else if (tableNum == 128 || tableNum == 0)
            {
              mode = ModelMode::ECHAM;
              echam_gribcodes(&gribcodes);
            }
          //  KNMI: HIRLAM model version 7.2 uses tableNum=1    (LAMH_D11*)
          //  KNMI: HARMONIE model version 36 uses tableNum=1   (grib*)
          //  (opreational NWP version) KNMI: HARMONIE model version 38 uses
          //  tableNum=253 (grib,grib_md) and tableNum=1 (grib_sfx) (research
          //  version)
          else if (tableNum == 1 || tableNum == 253)
            {
              mode = ModelMode::HIRLAM;
              hirlam_harmonie_gribcodes(&gribcodes);
            }
        }
      else
        {
          mode = ModelMode::ECHAM;
          echam_gribcodes(&gribcodes);
        }

      if (cdoVerbose) cdoPrint("Mode = %d  Center = %d  Param = %s", static_cast<int>(mode), instNum, paramstr);

      if (code <= 0 || code == 255)
        {
          vlistInqVarName(vlistID1, varID, varname);
          strtolower(varname);

          vlistInqVarStdname(vlistID1, varID, stdname);
          strtolower(stdname);

          code = echamcode_from_stdname(stdname);

          if (code < 0)
            {
              // clang-format off
              if      (sgeopotID == -1 && strcmp(varname, "geosp") == 0) code = gribcodes.geopot;
              else if (psID == -1 && strcmp(varname, "aps") == 0) code = gribcodes.ps;
              else if (psID == -1 && strcmp(varname, "ps") == 0) code = gribcodes.ps;
              else if (lnpsID == -1 && strcmp(varname, "lsp") == 0) code = gribcodes.lsp;
              else if (tempID == -1 && strcmp(varname, "t") == 0) code = gribcodes.temp;
              else if (humID == -1 && strcmp(varname, "q") == 0) code = gribcodes.hum;
              // else if ( geopotID  == -1 && strcmp(stdname, "geopotential_full") == 0 ) code = gribcodes.geopot;
              // else if ( strcmp(varname, "clwc") == 0 ) code = 246;
              // else if ( strcmp(varname, "ciwc") == 0 ) code = 247;
              // clang-format on
            }
        }

      // clang-format off
      if      (code == gribcodes.geopot && nlevel == 1) sgeopotID = varID;
      else if (code == gribcodes.geopot && nlevel == nhlevf) geopotID = varID;
      else if (code == gribcodes.temp && nlevel == nhlevf) tempID = varID;
      else if (code == gribcodes.hum && nlevel == nhlevf) humID = varID;
      else if (code == gribcodes.ps && nlevel == 1) psID = varID;
      else if (code == gribcodes.lsp && nlevel == 1) lnpsID = varID;
      else if (code == gribcodes.gheight && nlevel == nhlevf) gheightID = varID;
      // else if ( code == 246 ) clwcID    = varID;
      // else if ( code == 247 ) ciwcID    = varID;
      // clang-format on

      if (operatorID == SEALEVELPRESSURE) humID = -1;

      if (gridInqType(gridID) == GRID_SPECTRAL && zaxisInqType(zaxisID) == ZAXIS_HYBRID)
        cdoAbort("Spectral data on model level unsupported!");

      if (gridInqType(gridID) == GRID_SPECTRAL) cdoAbort("Spectral data unsupported!");
    }

  if (cdoVerbose)
    {
      cdoPrint("Found:");
      if (tempID != -1) cdoPrint("  %s", var_stdname(air_temperature));
      if (psID != -1) cdoPrint("  %s", var_stdname(surface_air_pressure));
      if (lnpsID != -1) cdoPrint("  LOG(%s)", var_stdname(surface_air_pressure));
      if (sgeopotID != -1) cdoPrint("  %s", var_stdname(surface_geopotential));
      if (geopotID != -1) cdoPrint("  %s", var_stdname(geopotential));
      if (gheightID != -1) cdoPrint("  %s", var_stdname(geopotential_height));
    }

  if (tempID == -1) cdoAbort("%s not found!", var_stdname(air_temperature));

  std::vector<double> array(gridsize);
  std::vector<double> sgeopot(gridsize);
  std::vector<double> ps(gridsize);
  std::vector<double> temp(gridsize * nhlevf);

  std::vector<double> half_press(gridsize * (nhlevf + 1));

  std::vector<double> hum;
  std::vector<double> gheight;
  if (operatorID == GHEIGHT)
    {
      if (humID == -1)
        cdoWarning("%s not found - using algorithm without %s!", var_stdname(specific_humidity), var_stdname(specific_humidity));
      else
        hum.resize(gridsize * nhlevf);

      gheight.resize(gridsize * (nhlevf + 1));
    }

  std::vector<double> full_press;
  std::vector<double> sealevelpressure;
  if (operatorID == SEALEVELPRESSURE)
    {
      full_press.resize(gridsize * nhlevf);

      surfaceID = zaxisFromName("surface");
      sealevelpressure.resize(gridsize);
    }

  if (zaxisIDh != -1 && sgeopotID == -1)
    {
      if (geopotID == -1)
        cdoWarning("%s not found - set to zero!", var_stdname(surface_geopotential));
      else
        cdoPrint("%s not found - using bottom layer of %s!", var_stdname(surface_geopotential), var_stdname(geopotential));

      arrayFill(gridsize, sgeopot.data(), 0.0);
    }

  presID = lnpsID;
  if (zaxisIDh != -1 && lnpsID == -1)
    {
      if (psID == -1)
        cdoAbort("%s not found!", var_stdname(surface_air_pressure));
      else
        presID = psID;
    }

  if (cdoVerbose)
    {
      if (presID == lnpsID)
        cdoPrint("using LOG(%s)", var_stdname(surface_air_pressure));
      else
        cdoPrint("using %s", var_stdname(surface_air_pressure));
    }

  int vlistID2 = vlistCreate();

  int var_id = -1;

  if (operatorID == GHEIGHT)
    {
      var_id = geopotential_height;
      varID = vlistDefVar(vlistID2, gridID, zaxisIDh, TIME_VARYING);
    }
  else if (operatorID == SEALEVELPRESSURE)
    {
      var_id = air_pressure_at_sea_level;
      varID = vlistDefVar(vlistID2, gridID, surfaceID, TIME_VARYING);
    }
  else
    cdoAbort("Internal problem, invalid operatorID: %d!", operatorID);

  vlistDefVarParam(vlistID2, varID, cdiEncodeParam(var_echamcode(var_id), 128, 255));
  vlistDefVarName(vlistID2, varID, var_name(var_id));
  vlistDefVarStdname(vlistID2, varID, var_stdname(var_id));
  vlistDefVarUnits(vlistID2, varID, var_units(var_id));

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
          zaxisID = vlistInqVarZaxis(vlistID1, varID);
          nlevel = zaxisInqSize(zaxisID);
          offset = gridsize * levelID;
          pstreamReadRecord(streamID1, array.data(), &nmiss);

          if (zaxisIDh != -1)
            {
              if (varID == sgeopotID)
                {
                  arrayCopy(gridsize, array.data(), sgeopot.data());
                }
              else if (varID == geopotID && sgeopotID == -1 && (levelID + 1) == nhlevf)
                {
                  arrayCopy(gridsize, array.data(), sgeopot.data());
                }
              else if (varID == presID)
                {
                  if (lnpsID != -1)
                    for (size_t i = 0; i < gridsize; ++i) ps[i] = exp(array[i]);
                  else if (psID != -1)
                    arrayCopy(gridsize, array.data(), ps.data());
                }
              else if (varID == tempID)
                arrayCopy(gridsize, array.data(), &temp[offset]);
              else if (varID == humID)
                arrayCopy(gridsize, array.data(), &hum[offset]);
            }
        }

      if (zaxisIDh != -1)
        {
          /* check range of ps_prog */
          arrayMinMaxMask(gridsize, ps.data(), NULL, &minval, &maxval);
          if (minval < MIN_PS || maxval > MAX_PS) cdoWarning("Surface pressure out of range (min=%g max=%g)!", minval, maxval);

          /* check range of surface geopot */
          arrayMinMaxMask(gridsize, sgeopot.data(), NULL, &minval, &maxval);
          if (minval < MIN_FIS || maxval > MAX_FIS) cdoWarning("Orography out of range (min=%g max=%g)!", minval, maxval);
        }

      varID = tempID;
      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      for (levelID = 0; levelID < nlevel; levelID++)
        {
          offset = gridsize * levelID;
          single2 = &temp[offset];

          arrayMinMaxMask(gridsize, single2, NULL, &minval, &maxval);
          if (minval < MIN_T || maxval > MAX_T)
            cdoWarning("Input temperature at level %d out of range (min=%g max=%g)!", levelID + 1, minval, maxval);
        }

      if (humID != -1)
        {
          varID = humID;
          nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
          for (levelID = 0; levelID < nlevel; levelID++)
            {
              offset = gridsize * levelID;
              single2 = &hum[offset];

              // corr_hum(gridsize, single2, MIN_Q);

              arrayMinMaxMask(gridsize, single2, NULL, &minval, &maxval);
              if (minval < -0.1 || maxval > MAX_Q)
                cdoWarning("Input humidity at level %d out of range (min=%g max=%g)!", levelID + 1, minval, maxval);
            }
        }

      if (operatorID == GHEIGHT)
        {
          presh(NULL, half_press.data(), vct, ps.data(), nhlevf, gridsize);
          arrayCopy(gridsize, sgeopot.data(), gheight.data() + gridsize * nhlevf);
          MakeGeopotHeight(gheight.data(), temp.data(), hum.data(), half_press.data(), gridsize, nhlevf);

          nmissout = 0;
          varID = 0;
          nlevel = nhlevf;
          for (levelID = 0; levelID < nlevel; levelID++)
            {
              pstreamDefRecord(streamID2, varID, levelID);
              pstreamWriteRecord(streamID2, gheight.data() + levelID * gridsize, nmissout);
            }
        }
      else if (operatorID == SEALEVELPRESSURE)
        {
          presh(full_press.data(), half_press.data(), vct, ps.data(), nhlevf, gridsize);

          extra_P(sealevelpressure.data(), &half_press[gridsize * (nhlevf)], &full_press[gridsize * (nhlevf - 1)], sgeopot.data(),
                  &temp[gridsize * (nhlevf - 1)], gridsize);

          pstreamDefRecord(streamID2, 0, 0);
          pstreamWriteRecord(streamID2, sealevelpressure.data(), 0);
        }
      else
        cdoAbort("Internal error");

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  vlistDestroy(vlistID2);

  if (vct) Free(vct);

  cdoFinish();

  return 0;
}
