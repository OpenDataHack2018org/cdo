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

      Output     output          ASCII output
      Output     outputf         Formatted output
      Output     outputint       Integer output
      Output     outputsrv       SERVICE output
      Output     outputext       EXTRA output
      Output     outputtab       Table output
*/

#include <cdi.h>

#include "cdo_int.h"
#include "grid.h"
#include "pstream_int.h"

static
void outputarr(size_t gridsize, std::vector<double> &array)
{
  for (size_t i = 0; i < gridsize; i++)
    {
      fprintf(stdout, "  arr[%zu] = %12.6g;\n", i, array[i]);
    }
}

static
void outputsp(size_t gridsize, std::vector<double> &array, long ntr)
{
  double minval = array[0];
  double maxval = array[0];
  arrayMinMax(gridsize, array.data(), &minval, &maxval);
  if ( /* T11 */ minval >= -1 && maxval <= 12 )
    {
      double *spc = array.data();
      for (long m = 0; m <= ntr; m++)
        {
          for (long n = m; n <= ntr; n++)
            {
              fprintf(stdout, "%3d", (int) *spc++);
              fprintf(stdout, "%3d", (int) *spc++);
            }
          fprintf(stdout, "\n");
        }
    }
}

static
void output(size_t gridsize, std::vector<double> &array)
{
  int nout = 0;
  for (size_t i = 0; i < gridsize; i++)
    {
      if (nout == 6)
        {
          nout = 0;
          fprintf(stdout, "\n");
        }
      fprintf(stdout, " %12.6g", array[i]);
      nout++;
    }
  fprintf(stdout, "\n");
}

static
void outputxyz(size_t gridsize, std::vector<double> &array, double missval, size_t nlon, size_t nlat, std::vector<double> &lon, std::vector<double> &lat)
{
  double fmin = 0;
  double x, y, z;
  for (size_t i = 0; i < gridsize; i++)
    if (!DBL_IS_EQUAL(array[i], missval))
      {
        if (array[i] < fmin) fmin = array[i];
        fprintf(stdout, "%g\t%g\t%g\t%g\n", lon[i], lat[i], array[i], array[i]);
      }
  const char *fname = "frontplane.xyz";
  FILE *fp = fopen(fname, "w");
  if (fp == NULL) cdoAbort("Open failed on %s", fname);
  // first front plane
  double dx = (lon[1] - lon[0]);
  double x0 = lon[0] - dx / 2;
  double y0 = lat[0] - dx / 2;
  double z0 = fmin;
  fprintf(fp, ">\n");
  for (size_t i = 0; i < nlon; ++i)
    {
      x = x0;
      y = y0;
      z = z0;
      fprintf(fp, "%g %g %g\n", x, y, z);
      x = x0;
      y = y0;
      z = array[i];
      fprintf(fp, "%g %g %g\n", x, y, z);
      x = x0 + dx;
      y = y0;
      fprintf(fp, "%g %g %g\n", x, y, z);
      x0 = x; /*y0 = y0;*/
      z0 = z;
    }
  x = x0;
  y = y0;
  z = fmin;
  fprintf(fp, "%g %g %g\n", x, y, z);
  x = lon[0] - dx / 2;
  fprintf(fp, "%g %g %g\n", x, y, z);

  // second front plane
  x0 = lon[0] - dx / 2;
  y0 = lat[0] - dx / 2;
  z0 = fmin;
  fprintf(fp, ">\n");
  for (size_t i = 0; i < nlat; ++i)
    {
      x = x0;
      y = y0;
      z = z0;
      fprintf(fp, "%g %g %g\n", x, y, z);
      x = x0;
      y = y0;
      z = array[i * nlon];
      fprintf(fp, "%g %g %g\n", x, y, z);
      x = x0;
      y = y0 + dx;
      fprintf(fp, "%g %g %g\n", x, y, z);
      /*x0 = x0;*/ y0 = y;
      z0 = z;
    }
  x = x0;
  y = y0;
  z = fmin;
  fprintf(fp, "%g %g %g\n", x, y, z);
  y = lat[0] - dx / 2;
  fprintf(fp, "%g %g %g\n", x, y, z);

  fclose(fp);
}

void *
Output(void *process)
{
  int varID;
  int gridID;
  int nrecs;
  int levelID;
  size_t nmiss;
  int nelem = 1;
  int len;
  int index;
  const char *format = NULL;
  char paramstr[32];
  char vdatestr[32], vtimestr[32];
  std::vector<double> grid_center_lon, grid_center_lat;
  char name[CDI_MAX_NAME];
  int year, month, day;
  std::vector<int> keys;
  int nkeys = 0, k;
  int nKeys;

  // clang-format off
  int Keylen[]           = {      0,        8,      11,      4,      8,     6,     6,     6,     6,      4,      4,          6,     10,      8,      5,       2,     2 };
  enum                     {knohead,   kvalue,  kparam,  kcode,  kname,  klon,  klat,  klev,  kbin,  kxind,  kyind,  ktimestep,  kdate,  ktime,  kyear,  kmonth,  kday };
  const char *Keynames[] = {"nohead",  "value", "param", "code", "name", "lon", "lat", "lev", "bin", "xind", "yind", "timestep", "date", "time", "year", "month", "day"};


  cdoInitialize(process);

  int OUTPUT    = cdoOperatorAdd("output",    0, 1, NULL);
  int OUTPUTINT = cdoOperatorAdd("outputint", 0, 0, NULL);
  int OUTPUTSRV = cdoOperatorAdd("outputsrv", 0, 0, NULL);
  int OUTPUTEXT = cdoOperatorAdd("outputext", 0, 0, NULL);
  int OUTPUTF   = cdoOperatorAdd("outputf",   0, 0, NULL);
  int OUTPUTTS  = cdoOperatorAdd("outputts",  0, 0, NULL);
  int OUTPUTFLD = cdoOperatorAdd("outputfld", 0, 0, NULL);
  int OUTPUTARR = cdoOperatorAdd("outputarr", 0, 0, NULL);
  int OUTPUTXYZ = cdoOperatorAdd("outputxyz", 0, 0, NULL);
  int OUTPUTTAB = cdoOperatorAdd("outputtab", 0, 0, NULL);
  // clang-format on

  UNUSED(OUTPUT);

  int operatorID = cdoOperatorID();
  bool opercplx = cdoOperatorF2(operatorID);

  if (operatorID == OUTPUTF)
    {
      operatorInputArg("format and number of elements [optional]");

      if (operatorArgc() < 1) cdoAbort("Too few arguments!");

      format = operatorArgv()[0];
      if (operatorArgc() == 2) nelem = parameter2int(operatorArgv()[1]);
    }
  else if (operatorID == OUTPUTTAB)
    {
      bool lhead = true;

      operatorInputArg("keys to print");

      int npar = operatorArgc();
      char **parnames = operatorArgv();

      if (cdoVerbose)
        for (int i = 0; i < npar; i++) cdoPrint("key %d = %s", i + 1, parnames[i]);

      keys.resize(npar);
      nkeys = 0;
      nKeys = sizeof(Keynames) / sizeof(char *);
      for (int i = 0; i < npar; i++)
        {
          for (k = 0; k < nKeys; ++k)
            {
              //	      len = strlen(parnames[i]);
              len = strlen(Keynames[k]);
              if (len < 3) len = 3;
              if (strncmp(parnames[i], Keynames[k], len) == 0)
                {
                  int len2 = strlen(parnames[i]);
                  if (len2 > len && parnames[i][len] != ':')
                    cdoAbort("Key parameter >%s< contains invalid character at position %d!", parnames[i], len + 1);

                  if (k == knohead)
                    lhead = false;
                  else
                    {
                      keys[nkeys++] = k;
                      if (len2 > len && parnames[i][len] == ':' && isdigit(parnames[i][len + 1]))
                        Keylen[k] = atoi(&parnames[i][len + 1]);
                    }
                  break;
                }
            }

          if (k == nKeys) cdoAbort("Key %s unsupported!", parnames[i]);
        }

      if (cdoVerbose)
        for (k = 0; k < nkeys; ++k)
          cdoPrint("keynr = %d  keyid = %d  keylen = %d  keyname = %s", k, keys[k], Keylen[keys[k]], Keynames[keys[k]]);

      if (lhead)
        {
          fprintf(stdout, "#");
          for (k = 0; k < nkeys; ++k)
            {
              len = Keylen[keys[k]];
              //   if ( k == 0 ) len -= 1;
              fprintf(stdout, "%*s ", len, Keynames[keys[k]]);
            }
          fprintf(stdout, "\n");
        }
    }

  for (int indf = 0; indf < cdoStreamCnt(); indf++)
    {
      int streamID = cdoStreamOpenRead(cdoStreamName(indf));

      int vlistID = cdoStreamInqVlist(streamID);

      int ngrids = vlistNgrids(vlistID);
      int ndiffgrids = 0;
      for (index = 1; index < ngrids; index++)
        if (vlistGrid(vlistID, 0) != vlistGrid(vlistID, index)) ndiffgrids++;

      if (ndiffgrids > 0) cdoAbort("Too many different grids!");

      gridID = vlistGrid(vlistID, 0);
      int gridtype = gridInqType(gridID);
      size_t gridsize = gridInqSize(gridID);
      size_t nwpv = (vlistNumber(vlistID) == CDI_COMP) ? 2 : 1;
      if (nwpv == 2 && opercplx == false) cdoAbort("Fields with complex numbers are not supported by this operator!");
      size_t gridsizemax = nwpv*gridInqSize(gridID);
      std::vector<double> array(gridsizemax);

      if (operatorID == OUTPUTFLD || operatorID == OUTPUTXYZ || operatorID == OUTPUTTAB)
        {
          if (gridInqType(gridID) == GRID_GME) gridID = gridToUnstructured(gridID, 0);

          if (gridInqType(gridID) != GRID_UNSTRUCTURED && gridInqType(gridID) != GRID_CURVILINEAR)
            gridID = gridToCurvilinear(gridID, 0);

          gridtype = gridInqType(gridID);

          grid_center_lon.resize(gridsize);
          grid_center_lat.resize(gridsize);
          gridInqXvals(gridID, grid_center_lon.data());
          gridInqYvals(gridID, grid_center_lat.data());

          /* Convert lat/lon units if required */
          {
            char units[CDI_MAX_NAME];
            gridInqXunits(gridID, units);
            grid_to_degree(units, gridsize, grid_center_lon.data(), "grid center lon");
            gridInqYunits(gridID, units);
            grid_to_degree(units, gridsize, grid_center_lat.data(), "grid center lat");
          }
        }

      int tsID = 0;
      int taxisID = vlistInqTaxis(vlistID);
      while ((nrecs = cdoStreamInqTimestep(streamID, tsID)))
        {
          int64_t vdate = taxisInqVdate(taxisID);
          int vtime = taxisInqVtime(taxisID);
          date2str(vdate, vdatestr, sizeof(vdatestr));
          time2str(vtime, vtimestr, sizeof(vtimestr));

          cdiDecodeDate(vdate, &year, &month, &day);

          for (int recID = 0; recID < nrecs; recID++)
            {
              pstreamInqRecord(streamID, &varID, &levelID);

              vlistInqVarName(vlistID, varID, name);
              int param = vlistInqVarParam(vlistID, varID);
              int code = vlistInqVarCode(vlistID, varID);
              int gridID = vlistInqVarGrid(vlistID, varID);
              int zaxisID = vlistInqVarZaxis(vlistID, varID);
              size_t nwpv = (vlistInqVarNumber(vlistID, varID) == CDI_COMP) ? 2 : 1;
              size_t gridsize = nwpv*gridInqSize(gridID);
              size_t nlon = gridInqXsize(gridID);
              size_t nlat = gridInqYsize(gridID);
              double level = cdoZaxisInqLevel(zaxisID, levelID);
              double missval = vlistInqVarMissval(vlistID, varID);

              cdiParamToString(param, paramstr, sizeof(paramstr));

              if (nlon * nlat != gridsize)
                {
                  nlon = gridsize;
                  nlat = 1;
                }

              pstreamReadRecord(streamID, array.data(), &nmiss);

              if (operatorID == OUTPUTSRV)
                fprintf(stdout, "%4d %8g %8lld %4d %8zu %8zu %d %d\n", code, level, vdate, vtime, nlon, nlat, 0, 0);

              if (operatorID == OUTPUTEXT) fprintf(stdout, "%8lld %4d %8g %8zu\n", vdate, code, level, gridsize);

              if (operatorID == OUTPUTINT)
                {
                  int nout = 0;
                  for (size_t i = 0; i < gridsize; i++)
                    {
                      if (nout == 8)
                        {
                          nout = 0;
                          fprintf(stdout, "\n");
                        }
                      fprintf(stdout, " %8d", (int) array[i]);
                      nout++;
                    }
                  fprintf(stdout, "\n");
                }
              else if (operatorID == OUTPUTF)
                {
                  int nout = 0;
                  for (size_t i = 0; i < gridsize; i++)
                    {
                      if (nout == nelem)
                        {
                          nout = 0;
                          fprintf(stdout, "\n");
                        }
                      fprintf(stdout, format, array[i]);
                      nout++;
                    }
                  fprintf(stdout, "\n");
                }
              else if (operatorID == OUTPUTTS)
                {
                  char vdatestr[32], vtimestr[32];

                  if (gridsize > 1) cdoAbort("operator works only with one gridpoint!");

                  date2str(vdate, vdatestr, sizeof(vdatestr));
                  time2str(vtime, vtimestr, sizeof(vtimestr));

                  fprintf(stdout, "%s %s %12.12g\n", vdatestr, vtimestr, array[0]);
                }
              else if (operatorID == OUTPUTFLD)
                {
                  int hour, minute, second;
                  cdiDecodeTime(vtime, &hour, &minute, &second);
                  double xdate = vdate - (vdate / 100) * 100 + (hour * 3600 + minute * 60 + second) / 86400.;
                  for (size_t i = 0; i < gridsize; i++)
                    if (!DBL_IS_EQUAL(array[i], missval))
                      fprintf(stdout, "%g\t%g\t%g\t%g\n", xdate, grid_center_lat[i], grid_center_lon[i], array[i]);
                }
              else if (operatorID == OUTPUTTAB)
                {
                  bool l2d = false;
                  int xsize = gridInqXsize(gridID);
                  // int ysize = gridInqYsize(gridID);
                  if (gridtype == GRID_CURVILINEAR) l2d = true;

                  for (size_t i = 0; i < gridsize; i++)
                    {
                      int yind = i;
                      int xind = i;
                      if (l2d)
                        {
                          yind /= xsize;
                          xind -= yind * xsize;
                        }
                      double lon = grid_center_lon[i];
                      double lat = grid_center_lat[i];

                      for (k = 0; k < nkeys; ++k)
                        {
                          len = Keylen[keys[k]];
                          switch (keys[k])
                            {
                            case kvalue: fprintf(stdout, "%*g ", len, array[i]); break;
                            case kparam: fprintf(stdout, "%*s ", len, paramstr); break;
                            case kcode: fprintf(stdout, "%*d ", len, code); break;
                            case kname: fprintf(stdout, "%*s ", len, name); break;
                            case klon: fprintf(stdout, "%*g ", len, lon); break;
                            case klat: fprintf(stdout, "%*g ", len, lat); break;
                            case klev: fprintf(stdout, "%*g ", len, level); break;
                            case kbin: fprintf(stdout, "%*g ", len, level); break;
                            case kxind: fprintf(stdout, "%*d ", len, xind + 1); break;
                            case kyind: fprintf(stdout, "%*d ", len, yind + 1); break;
                            case ktimestep: fprintf(stdout, "%*d ", len, tsID + 1); break;
                            case kdate: fprintf(stdout, "%*s ", len, vdatestr); break;
                            case ktime: fprintf(stdout, "%*s ", len, vtimestr); break;
                            case kyear: fprintf(stdout, "%*d ", len, year); break;
                            case kmonth: fprintf(stdout, "%*d ", len, month); break;
                            case kday: fprintf(stdout, "%*d ", len, day); break;
                            }
                        }
                      fprintf(stdout, "\n");
                    }
                }
              else if (operatorID == OUTPUTXYZ)
                {
                  if (tsID == 0 && recID == 0)
                    {
                      outputxyz(gridsize, array, missval, nlon, nlat, grid_center_lon, grid_center_lat);
                    }
                }
              else if (operatorID == OUTPUTARR)
                {
                  outputarr(gridsize, array);
                }
              else
                {
                  if (gridInqType(gridID) == GRID_SPECTRAL && gridsize <= 156)
                    {
                      outputsp(gridsize, array, gridInqTrunc(gridID));
                    }
                  else
                    {
                      output(gridsize, array);
                    }
                }
            }

          tsID++;
        }

      pstreamClose(streamID);
    }

  cdoFinish();

  return 0;
}
