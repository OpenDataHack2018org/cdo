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

      Change     chcode          Change code number
      Change     chtabnum        Change GRIB1 parameter table number
      Change     chparam         Change parameter identifier
      Change     chname          Change variable name
      Change     chlevel         Change level
      Change     chlevelc        Change level of one code
      Change     chlevelv        Change level of one variable
      Change     chltype         Change GRIB level type
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

int stringToParam(const char *paramstr);

#define MAXARG 16384

void *
Change(void *process)
{
  int nrecs;
  int varID = 0, levelID;
  int chints[MAXARG];
  char *chnames[MAXARG];
  char varname[CDI_MAX_NAME];
  char *chname = NULL;
  int chcode = 0;
  int i;
  size_t nmiss;
  int zaxisID1, zaxisID2, k, index;
  double chlevels[MAXARG];
  int chltypes[MAXARG];

  cdoInitialize(process);

  // clang-format off
  int CHCODE   = cdoOperatorAdd("chcode",   0, 0, "pairs of old and new code numbers");
  int CHTABNUM = cdoOperatorAdd("chtabnum", 0, 0, "pairs of old and new GRIB1 table numbers");
  int CHPARAM  = cdoOperatorAdd("chparam",  0, 0, "pairs of old and new parameter identifiers");
  int CHNAME   = cdoOperatorAdd("chname",   0, 0, "pairs of old and new variable names");
  int CHUNIT   = cdoOperatorAdd("chunit",   0, 0, "pairs of old and new variable units");
  int CHLEVEL  = cdoOperatorAdd("chlevel",  0, 0, "pairs of old and new levels");
  int CHLEVELC = cdoOperatorAdd("chlevelc", 0, 0, "code number, old and new level");
  int CHLEVELV = cdoOperatorAdd("chlevelv", 0, 0, "variable name, old and new level");
  int CHLTYPE  = cdoOperatorAdd("chltype",  0, 0, "pairs of old and new type");
  // clang-format on

  int operatorID = cdoOperatorID();

  operatorInputArg(cdoOperatorEnter(operatorID));

  int nch = operatorArgc();

  if (operatorID == CHCODE || operatorID == CHTABNUM)
    {
      if (nch % 2) cdoAbort("Odd number of input arguments!");
      for (i = 0; i < nch; i++) chints[i] = parameter2int(operatorArgv()[i]);
    }
  else if (operatorID == CHPARAM || operatorID == CHNAME || operatorID == CHUNIT)
    {
      if (nch % 2) cdoAbort("Odd number of input arguments!");
      for (i = 0; i < nch; i++) chnames[i] = operatorArgv()[i];
    }
  else if (operatorID == CHLEVEL)
    {
      if (nch % 2) cdoAbort("Odd number of input arguments!");
      for (i = 0; i < nch; i++) chlevels[i] = parameter2double(operatorArgv()[i]);
    }
  else if (operatorID == CHLEVELC)
    {
      operatorCheckArgc(3);

      chcode = parameter2int(operatorArgv()[0]);
      chlevels[0] = parameter2double(operatorArgv()[1]);
      chlevels[1] = parameter2double(operatorArgv()[2]);
    }
  else if (operatorID == CHLEVELV)
    {
      operatorCheckArgc(3);

      chname = operatorArgv()[0];
      chlevels[0] = parameter2double(operatorArgv()[1]);
      chlevels[1] = parameter2double(operatorArgv()[2]);
    }
  else if (operatorID == CHLTYPE)
    {
      if (nch % 2) cdoAbort("Odd number of input arguments!");
      for (i = 0; i < nch; i++) chltypes[i] = parameter2int(operatorArgv()[i]);
    }

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  if (operatorID == CHCODE)
    {
      int nvars = vlistNvars(vlistID2);
      for (varID = 0; varID < nvars; varID++)
        {
          int code = vlistInqVarCode(vlistID2, varID);
          for (i = 0; i < nch; i += 2)
            if (code == chints[i]) vlistDefVarCode(vlistID2, varID, chints[i + 1]);
        }
    }
  else if (operatorID == CHTABNUM)
    {
      int tableID;
      int nvars = vlistNvars(vlistID2);
      for (varID = 0; varID < nvars; varID++)
        {
          int tabnum = tableInqNum(vlistInqVarTable(vlistID2, varID));
          for (i = 0; i < nch; i += 2)
            if (tabnum == chints[i])
              {
                tableID = tableDef(-1, chints[i + 1], NULL);
                vlistDefVarTable(vlistID2, varID, tableID);
              }
        }
    }
  else if (operatorID == CHPARAM)
    {
      int nvars = vlistNvars(vlistID2);
      for (varID = 0; varID < nvars; varID++)
        {
          int param = vlistInqVarParam(vlistID2, varID);
          if (cdoVerbose)
            {
              int pnum, pcat, pdis;
              cdiDecodeParam(param, &pnum, &pcat, &pdis);
              cdoPrint("pnum, pcat, pdis: %d.%d.%d", pnum, pcat, pdis);
            }
          for (i = 0; i < nch; i += 2)
            if (param == stringToParam(chnames[i])) vlistDefVarParam(vlistID2, varID, stringToParam(chnames[i + 1]));
        }
    }
  else if (operatorID == CHNAME)
    {
      int nvars = vlistNvars(vlistID2);
      for (varID = 0; varID < nvars; varID++)
        {
          vlistInqVarName(vlistID2, varID, varname);
          for (i = 0; i < nch; i += 2)
            if (strcmp(varname, chnames[i]) == 0) vlistDefVarName(vlistID2, varID, chnames[i + 1]);
        }
    }
  else if (operatorID == CHUNIT)
    {
      int nvars = vlistNvars(vlistID2);
      for (varID = 0; varID < nvars; varID++)
        {

          vlistInqVarUnits(vlistID2, varID, varname);
          for (i = 0; i < nch; i += 2)
            if (strcmp(varname, chnames[i]) == 0) vlistDefVarUnits(vlistID2, varID, chnames[i + 1]);
        }
    }
  else if (operatorID == CHLEVEL)
    {
      int nzaxis = vlistNzaxis(vlistID2);
      for (index = 0; index < nzaxis; index++)
        {
          zaxisID1 = vlistZaxis(vlistID2, index);
          if (zaxisInqLevels(zaxisID1, NULL))
            {
              int nlevs = zaxisInqSize(zaxisID1);
              std::vector<double> levels(nlevs);
              std::vector<double> newlevels(nlevs);
              zaxisInqLevels(zaxisID1, &levels[0]);

              for (k = 0; k < nlevs; k++) newlevels[k] = levels[k];

              int nfound = 0;
              for (i = 0; i < nch; i += 2)
                for (k = 0; k < nlevs; k++)
                  if (fabs(levels[k] - chlevels[i]) < 0.0001) nfound++;

              if (nfound)
                {
                  int zaxisID2 = zaxisDuplicate(zaxisID1);
                  for (i = 0; i < nch; i += 2)
                    for (k = 0; k < nlevs; k++)
                      if (fabs(levels[k] - chlevels[i]) < 0.001) newlevels[k] = chlevels[i + 1];

                  zaxisDefLevels(zaxisID2, &newlevels[0]);
                  vlistChangeZaxis(vlistID2, zaxisID1, zaxisID2);
                }
            }
        }
    }
  else if (operatorID == CHLEVELC || operatorID == CHLEVELV)
    {
      int nvars = vlistNvars(vlistID2);
      if (operatorID == CHLEVELC)
        {
          for (varID = 0; varID < nvars; varID++)
            {
              int code = vlistInqVarCode(vlistID2, varID);
              if (code == chcode) break;
            }
          if (varID == nvars) cdoAbort("Code %d not found!", chcode);
        }
      else
        {
          for (varID = 0; varID < nvars; varID++)
            {
              vlistInqVarName(vlistID2, varID, varname);
              if (strcmp(varname, chname) == 0) break;
            }
          if (varID == nvars) cdoAbort("Variable name %s not found!", chname);
        }

      int zaxisID1 = vlistInqVarZaxis(vlistID2, varID);
      if (zaxisInqLevels(zaxisID1, NULL))
        {
          int nlevs = zaxisInqSize(zaxisID1);
          std::vector<double> levels(nlevs);
          zaxisInqLevels(zaxisID1, &levels[0]);
          int nfound = 0;
          for (k = 0; k < nlevs; k++)
            if (fabs(levels[k] - chlevels[0]) < 0.0001) nfound++;

          if (nfound)
            {
              int zaxisID2 = zaxisDuplicate(zaxisID1);
              for (k = 0; k < nlevs; k++)
                if (fabs(levels[k] - chlevels[0]) < 0.001) levels[k] = chlevels[1];

              zaxisDefLevels(zaxisID2, &levels[0]);
              vlistChangeVarZaxis(vlistID2, varID, zaxisID2);
            }
          else
            cdoAbort("Level %g not found!", chlevels[0]);
        }
    }
  else if (operatorID == CHLTYPE)
    {
      int ltype, ltype1, ltype2;
      int nzaxis = vlistNzaxis(vlistID2);
      for (index = 0; index < nzaxis; index++)
        {
          zaxisID1 = vlistZaxis(vlistID2, index);
          zaxisID2 = zaxisDuplicate(zaxisID1);

          ltype = zaxisInqLtype(zaxisID1);

          for (i = 0; i < nch; i += 2)
            {
              ltype1 = chltypes[i];
              ltype2 = chltypes[i + 1];

              if (ltype1 == ltype)
                {
                  zaxisChangeType(zaxisID2, ZAXIS_GENERIC);
                  zaxisDefLtype(zaxisID2, ltype2);
                  vlistChangeZaxis(vlistID2, zaxisID1, zaxisID2);
                }
            }
        }
    }

  size_t gridsizemax = vlistGridsizeMax(vlistID2);
  std::vector<double> array(gridsizemax);

  int streamID2 = CDI_UNDEFID;

  int tsID1 = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID1)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      if (streamID2 == CDI_UNDEFID)
        {
          streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
          pstreamDefVlist(streamID2, vlistID2);
        }
      pstreamDefTimestep(streamID2, tsID1);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);
          pstreamDefRecord(streamID2, varID, levelID);

          pstreamReadRecord(streamID1, &array[0], &nmiss);
          pstreamWriteRecord(streamID2, &array[0], nmiss);
        }
      tsID1++;
    }

  pstreamClose(streamID1);
  pstreamClose(streamID2);

  cdoFinish();

  return 0;
}
