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

      Setzaxis   setzaxis        Set zaxis
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

void genLayerBounds(int nlev, double *levels, double *lbounds, double *ubounds);

int
getkeyval_dp(const char *keyval, const char *key, double *val)
{
  int status = 0;
  size_t keylen = strlen(key);

  if (strncmp(keyval, key, keylen) == 0)
    {
      const char *pkv = keyval + keylen;
      if (pkv[0] == '=' && pkv[1] != 0)
        {
          *val = parameter2double(&pkv[1]);
          status = 1;
        }
      else
        {
          cdoAbort("Syntax error for parameter %s!", keyval);
        }
    }

  return status;
}

void *
Setzaxis(void *process)
{
  int nrecs;
  int varID, levelID;
  int zaxisID1, zaxisID2 = -1;
  int nzaxis, index;
  size_t nmiss;
  bool lztop = false, lzbot = false;
  double ztop = 0, zbot = 0;

  cdoInitialize(process);

  // clang-format off
  int SETZAXIS       = cdoOperatorAdd("setzaxis",        0, 0, "zaxis description file");
  int GENLEVELBOUNDS = cdoOperatorAdd("genlevelbounds",  0, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  if (operatorID == SETZAXIS)
    {
      operatorInputArg(cdoOperatorEnter(operatorID));
      operatorCheckArgc(1);
      zaxisID2 = cdoDefineZaxis(operatorArgv()[0]);
    }
  else if (operatorID == GENLEVELBOUNDS)
    {
      unsigned npar = operatorArgc();
      char **parnames = operatorArgv();

      for (unsigned i = 0; i < npar; i++)
        {
          if (cdoVerbose) cdoPrint("keyval[%d]: %s", i + 1, parnames[i]);

          if (!lzbot && getkeyval_dp(parnames[i], "zbot", &zbot))
            lzbot = true;
          else if (!lztop && getkeyval_dp(parnames[i], "ztop", &ztop))
            lztop = true;
          else
            cdoAbort("Parameter >%s< unsupported! Supported parameter are: zbot, ztop", parnames[i]);
        }
    }

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());

  if (operatorID == SETZAXIS)
    {
      int found = 0;
      nzaxis = vlistNzaxis(vlistID1);
      for (index = 0; index < nzaxis; index++)
        {
          zaxisID1 = vlistZaxis(vlistID1, index);

          if (zaxisInqSize(zaxisID1) == zaxisInqSize(zaxisID2))
            {
              vlistChangeZaxisIndex(vlistID2, index, zaxisID2);
              found++;
            }
        }
      if (!found) cdoWarning("No zaxis with %d levels found!", zaxisInqSize(zaxisID2));
    }
  else if (operatorID == GENLEVELBOUNDS)
    {
      nzaxis = vlistNzaxis(vlistID1);
      for (index = 0; index < nzaxis; index++)
        {
          zaxisID1 = vlistZaxis(vlistID1, index);
          int nlev = zaxisInqSize(zaxisID1);
          if (nlev > 1)
            {
              std::vector<double> levels(nlev);
              std::vector<double> lbounds(nlev);
              std::vector<double> ubounds(nlev);

              cdoZaxisInqLevels(zaxisID1, levels.data());
              zaxisID2 = zaxisDuplicate(zaxisID1);
              if (!zaxisInqLevels(zaxisID1, NULL)) zaxisDefLevels(zaxisID2, levels.data());

              genLayerBounds(nlev, levels.data(), lbounds.data(), ubounds.data());

              if (lzbot) lbounds[0] = zbot;
              if (lztop) ubounds[nlev - 1] = ztop;
              zaxisDefLbounds(zaxisID2, lbounds.data());
              zaxisDefUbounds(zaxisID2, ubounds.data());
              vlistChangeZaxisIndex(vlistID2, index, zaxisID2);
            }
        }
    }

  pstreamDefVlist(streamID2, vlistID2);

  size_t gridsize = vlistGridsizeMax(vlistID1);
  if (vlistNumber(vlistID1) != CDI_REAL) gridsize *= 2;
  std::vector<double> array(gridsize);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamDefRecord(streamID2, varID, levelID);

          pstreamReadRecord(streamID1, array.data(), &nmiss);
          pstreamWriteRecord(streamID2, array.data(), nmiss);
        }

      tsID++;
    }

  pstreamClose(streamID1);
  pstreamClose(streamID2);

  cdoFinish();

  return 0;
}
