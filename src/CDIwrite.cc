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
#include "grid.h"
#include "timer.h"

const char *filetypestr(int filetype);
const char *datatypestr(int datatype);

static void
print_stat(const char *sinfo, int memtype, int datatype, int filetype, off_t nvalues, double data_size, double file_size,
           double tw)
{
  nvalues /= 1000000;
  data_size /= 1024. * 1024. * 1024.;

  double rout = 0;
  if (tw > 0) rout = nvalues / tw;

  if (memtype == MEMTYPE_FLOAT)
    cdoPrint("%s Wrote %.1f GB of 32 bit floats to %s %s, %.1f MVal/s", sinfo, data_size, datatypestr(datatype),
             filetypestr(filetype), rout);
  else
    cdoPrint("%s Wrote %.1f GB of 64 bit floats to %s %s, %.1f MVal/s", sinfo, data_size, datatypestr(datatype),
             filetypestr(filetype), rout);

  file_size /= 1024. * 1024. * 1024.;

  rout = 0;
  if (tw > 0) rout = 1024 * file_size / tw;

  cdoPrint("%s Wrote %.1f GB in %.1f seconds, total %.1f MB/s", sinfo, file_size, tw, rout);
}

void *
CDIwrite(void *process)
{
  int memtype = CDO_Memtype;
  int nvars = 10, nlevs = 0, ntimesteps = 30;
  const char *defaultgrid = "global_.2";
  int tsID, varID, levelID;
  int i;
  int zaxisID;
  int filetype = -1, datatype = -1;
  int irun, nruns = 1;
  unsigned int seed = 1;
  char sinfo[64];
  off_t nvalues = 0;
  double file_size = 0, data_size = 0;
  double tw, tw0, t0, twsum = 0;

  srand(seed);
  sinfo[0] = 0;

  cdoInitialize(process);

  if (cdoVerbose) cdoPrint("parameter: <nruns, <grid, <nlevs, <ntimesteps, <nvars>>>>>");

  if (operatorArgc() > 5) cdoAbort("Too many arguments!");

  const char *gridfile = defaultgrid;
  if (operatorArgc() >= 1) nruns = parameter2int(operatorArgv()[0]);
  if (operatorArgc() >= 2) gridfile = operatorArgv()[1];
  if (operatorArgc() >= 3) nlevs = parameter2int(operatorArgv()[2]);
  if (operatorArgc() >= 4) ntimesteps = parameter2int(operatorArgv()[3]);
  if (operatorArgc() >= 5) nvars = parameter2int(operatorArgv()[4]);

  if (nruns < 0) nruns = 0;
  if (nruns > 9999) nruns = 9999;

  if (nlevs <= 0) nlevs = 1;
  if (nlevs > 255) nlevs = 255;
  if (ntimesteps <= 0) ntimesteps = 1;
  if (nvars <= 0) nvars = 1;

  int gridID = cdoDefineGrid(gridfile);
  size_t gridsize = gridInqSize(gridID);

  if (nlevs == 1)
    zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);
  else
    {
      std::vector<double> levels(nlevs);
      for (i = 0; i < nlevs; ++i) levels[i] = 100 * i;
      zaxisID = zaxisCreate(ZAXIS_HEIGHT, nlevs);
      zaxisDefLevels(zaxisID, &levels[0]);
    }

  if (cdoVerbose)
    {
      cdoPrint("nruns      : %d", nruns);
      cdoPrint("gridsize   : %zu", gridsize);
      cdoPrint("nlevs      : %d", nlevs);
      cdoPrint("ntimesteps : %d", ntimesteps);
      cdoPrint("nvars      : %d", nvars);
    }

  std::vector<double> array(gridsize);
  std::vector<double> xvals(gridsize);
  std::vector<double> yvals(gridsize);

  int gridID2 = gridID;
  if (gridInqType(gridID) == GRID_GME) gridID2 = gridToUnstructured(gridID, 0);

  if (gridInqType(gridID) != GRID_UNSTRUCTURED && gridInqType(gridID) != GRID_CURVILINEAR)
    gridID2 = gridToCurvilinear(gridID, 0);

  gridInqXvals(gridID2, &xvals[0]);
  gridInqYvals(gridID2, &yvals[0]);

  /* Convert lat/lon units if required */
  char units[CDI_MAX_NAME];
  gridInqXunits(gridID2, units);
  grid_to_radian(units, gridsize, &xvals[0], "grid center lon");
  gridInqYunits(gridID2, units);
  grid_to_radian(units, gridsize, &yvals[0], "grid center lat");

  for (size_t i = 0; i < gridsize; i++) array[i] = 2 - cos(acos(cos(xvals[i]) * cos(yvals[i])) / 1.2);

  std::vector<std::vector<std::vector<double>>> vars(nvars);
  for (varID = 0; varID < nvars; varID++)
    {
      vars[varID].resize(nlevs);
      for (levelID = 0; levelID < nlevs; levelID++)
        {
          vars[varID][levelID].resize(gridsize);
          for (size_t i = 0; i < gridsize; ++i) vars[varID][levelID][i] = varID + array[i] * (levelID + 1);
        }
    }

  std::vector<float> farray;
  if (memtype == MEMTYPE_FLOAT) farray.resize(gridsize);

  int vlistID = vlistCreate();

  for (i = 0; i < nvars; ++i)
    {
      varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
      vlistDefVarParam(vlistID, varID, cdiEncodeParam(varID + 1, 255, 255));
      //    vlistDefVarName(vlistID, varID, );
    }

  int taxisID = taxisCreate(TAXIS_RELATIVE);
  vlistDefTaxis(vlistID, taxisID);

  // vlistDefNtsteps(vlistID, 1);

  for (irun = 0; irun < nruns; ++irun)
    {
      tw0 = timer_val(timer_write);
      data_size = 0;
      nvalues = 0;

      int streamID = cdoStreamOpenWrite(cdoStreamName(0), cdoFiletype());

      pstreamDefVlist(streamID, vlistID);

      filetype = pstreamInqFiletype(streamID);
      datatype = vlistInqVarDatatype(vlistID, 0);
      if (datatype == CDI_UNDEFID) datatype = CDI_DATATYPE_FLT32;

      int julday = date_to_julday(CALENDAR_PROLEPTIC, 19870101);

      t0 = timer_val(timer_write);

      for (tsID = 0; tsID < ntimesteps; tsID++)
        {
          int vdate = julday_to_date(CALENDAR_PROLEPTIC, julday + tsID);
          int vtime = 0;
          taxisDefVdate(taxisID, vdate);
          taxisDefVtime(taxisID, vtime);
          pstreamDefTimestep(streamID, tsID);

          for (varID = 0; varID < nvars; varID++)
            {
              for (levelID = 0; levelID < nlevs; levelID++)
                {
                  nvalues += gridsize;
                  pstreamDefRecord(streamID, varID, levelID);
                  if (memtype == MEMTYPE_FLOAT)
                    {
                      for (size_t i = 0; i < gridsize; ++i) farray[i] = vars[varID][levelID][i];
                      pstreamWriteRecordF(streamID, &farray[0], 0);
                      data_size += gridsize * 4;
                    }
                  else
                    {
                      pstreamWriteRecord(streamID, &vars[varID][levelID][0], 0);
                      data_size += gridsize * 8;
                    }
                }
            }

          if (cdoVerbose)
            {
              tw = timer_val(timer_write) - t0;
              t0 = timer_val(timer_write);
              cdoPrint("Timestep %d: %.2f seconds", tsID + 1, tw);
            }
        }

      pstreamClose(streamID);

      tw = timer_val(timer_write) - tw0;
      twsum += tw;

      file_size = (double) fileSize(cdoGetStreamName(0).c_str());

      if (nruns > 1) snprintf(sinfo, sizeof(sinfo), "(run %d)", irun + 1);

      print_stat(sinfo, memtype, datatype, filetype, nvalues, data_size, file_size, tw);
    }

  if (nruns > 1) print_stat("(mean)", memtype, datatype, filetype, nvalues, data_size, file_size, twsum / nruns);

  vlistDestroy(vlistID);

  cdoFinish();

  return 0;
}
