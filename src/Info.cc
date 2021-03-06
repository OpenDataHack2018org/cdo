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

      Info       info            Dataset information
      Info       map             Dataset information and simple map
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "text.h"
#include "array.h"
#include "datetime.h"

static void
printMap(int nlon, int nlat, double *array, double missval, double min, double max)
{
  /* source code from PINGO */
  int ilon, ilat, i;
  double x, a, b;
  double level[10];
  int min_n, max_n;
  int bmin = 1, bmax = 1;
  unsigned char c;

  double step = (max - min) / 10;

  if (IS_NOT_EQUAL(step, 0))
    {
      a = pow(10, floor(log(step) / M_LN10));
      b = step / a;

      if (b > 5)
        b = 0.5 * ceil(b / 0.5);
      else if (b > 2)
        b = 0.2 * ceil(b / 0.2);
      else if (b > 1)
        b = 0.1 * ceil(b / 0.1);
      else
        b = 1;

      step = b * a;

      if (min < 0 && max > 0)
        {
          min_n = (int) floor(10 * (-min) / (max - min) - 0.5);
          max_n = (int) ceil(10 * (-min) / (max - min) - 0.5);
          level[min_n] = 0;
          for (i = min_n - 1; i >= 0; i--) level[i] = level[i + 1] - step;
          for (i = max_n; i < 9; i++) level[i] = level[i - 1] + step;
        }
      else
        {
          level[0] = step * ceil(min / step + 0.5);
          for (i = 1; i < 9; i++) level[i] = level[i - 1] + step;
        }
    }
  else
    for (i = 0; i < 9; i++) level[i] = min;

  fputc('\n', stdout);
  fflush(stdout);

  if (nlon >= 1000)
    {
      printf("%.*s", nlat < 10 ? 2 : nlat < 100 ? 3 : nlat < 1000 ? 4 : 5, "     ");
      for (ilon = 0; ilon < nlon; ilon++) printf("%d", ((ilon + 1) / 1000) % 10);
      putchar('\n');
      fflush(stdout);
    }

  if (nlon >= 100)
    {
      printf("%.*s", nlat < 10 ? 2 : nlat < 100 ? 3 : nlat < 1000 ? 4 : 5, "     ");
      for (ilon = 0; ilon < nlon; ilon++) printf("%d", ((ilon + 1) / 100) % 10);
      putchar('\n');
      fflush(stdout);
    }

  if (nlon >= 10)
    {
      printf("%.*s", nlat < 10 ? 2 : nlat < 100 ? 3 : nlat < 1000 ? 4 : 5, "     ");
      for (ilon = 0; ilon < nlon; ilon++) printf("%d", ((ilon + 1) / 10) % 10);
      putchar('\n');
      fflush(stdout);
    }

  printf("%.*s", nlat < 10 ? 2 : nlat < 100 ? 3 : nlat < 1000 ? 4 : 5, "     ");

  for (ilon = 0; ilon < nlon; ilon++) printf("%d", (ilon + 1) % 10);
  putchar('\n');
  fflush(stdout);
  putchar('\n');
  fflush(stdout);

  for (ilat = 0; ilat < nlat; ilat++)
    {
      printf("%0*d ", nlat < 10 ? 1 : nlat < 100 ? 2 : nlat < 1000 ? 3 : 4, ilat + 1);
      for (ilon = 0; ilon < nlon; ilon++)
        {
          x = array[ilat * nlon + ilon];
          if (DBL_IS_EQUAL(x, missval))
            c = '.';
          else if (DBL_IS_EQUAL(x, min) && !DBL_IS_EQUAL(min, max))
            c = 'm';
          else if (DBL_IS_EQUAL(x, max) && !DBL_IS_EQUAL(min, max))
            c = 'M';
          else if (DBL_IS_EQUAL(x, 0.))
            c = '*';
          else if (x < 0)
            {
              c = '9';
              for (i = 0; i < 9; i++)
                if (level[i] > x)
                  {
                    c = i + '0';
                    break;
                  }
            }
          else
            {
              c = '0';
              for (i = 8; i >= 0; i--)
                if (level[i] < x)
                  {
                    c = i + 1 + '0';
                    break;
                  }
            }

          if (c == '0')
            {
              set_text_color(stdout, BRIGHT, BLUE);
            }
          else if (c == '1')
            {
              set_text_color(stdout, RESET, BLUE);
            }
          else if (c == '2')
            {
              set_text_color(stdout, BRIGHT, CYAN);
            }
          else if (c == '3')
            {
              set_text_color(stdout, RESET, CYAN);
            }
          else if (c == '4')
            {
              set_text_color(stdout, RESET, GREEN);
            }
          else if (c == '5')
            {
              set_text_color(stdout, RESET, YELLOW);
            }
          else if (c == '6')
            {
              set_text_color(stdout, RESET, RED);
            }
          else if (c == '7')
            {
              set_text_color(stdout, BRIGHT, RED);
            }
          else if (c == '8')
            {
              set_text_color(stdout, RESET, MAGENTA);
            }
          else if (c == '9')
            {
              set_text_color(stdout, BRIGHT, MAGENTA);
            }
          else if (c == 'M')
            {
              if (bmax)
                {
                  bmax = 0;
                  set_text_color(stdout, BLINK, BLACK);
                }
              else
                {
                  set_text_color(stdout, RESET, BLACK);
                }
            }
          else if (c == 'm')
            {
              if (bmin)
                {
                  bmin = 0;
                  set_text_color(stdout, BLINK, BLACK);
                }
              else
                {
                  set_text_color(stdout, RESET, BLACK);
                }
            }
          putchar(c);
          reset_text_color(stdout);
        }
      printf(" %0*d\n", nlat < 10 ? 1 : nlat < 100 ? 2 : nlat < 1000 ? 3 : 4, ilat + 1);
      fflush(stdout);
    }
  putchar('\n');
  fflush(stdout);

  if (nlon >= 1000)
    {
      printf("%.*s", nlat < 10 ? 2 : nlat < 100 ? 3 : nlat < 1000 ? 4 : 5, "     ");
      for (ilon = 0; ilon < nlon; ilon++) printf("%d", ((ilon + 1) / 1000) % 10);
      putchar('\n');
      fflush(stdout);
    }

  if (nlon >= 100)
    {
      printf("%.*s", nlat < 10 ? 2 : nlat < 100 ? 3 : nlat < 1000 ? 4 : 5, "     ");
      for (ilon = 0; ilon < nlon; ilon++) printf("%d", ((ilon + 1) / 100) % 10);
      putchar('\n');
      fflush(stdout);
    }

  if (nlon >= 10)
    {
      printf("%.*s", nlat < 10 ? 2 : nlat < 100 ? 3 : nlat < 1000 ? 4 : 5, "     ");
      for (ilon = 0; ilon < nlon; ilon++) printf("%d", ((ilon + 1) / 10) % 10);
      putchar('\n');
      fflush(stdout);
    }

  printf("%.*s", nlat < 10 ? 2 : nlat < 100 ? 3 : nlat < 1000 ? 4 : 5, "     ");
  for (ilon = 0; ilon < nlon; ilon++) printf("%d", (ilon + 1) % 10);
  putchar('\n');
  fflush(stdout);
  putchar('\n');
  fflush(stdout);

  for (i = 0; i < 10; i++)
    {
      printf("%d=%c%+9.3e,%+9.3e%c%s", (int) i, i == 0 || level[i - 1] >= 0 ? '[' : '[', i == 0 ? min : level[i - 1],
             i == 9 ? max : level[i], i == 9 || level[i] <= 0 ? ']' : ']', i != 2 && i != 5 && i != 8 ? "  " : "");

      if (i == 2 || i == 5 || i == 8)
        {
          fputc('\n', stdout);
          fflush(stdout);
        }
    }

  printf("*=0  .=miss  m=min=%+9.3e  M=max=%+9.3e\n", min, max);
  fflush(stdout);
  putchar('\n');
  fflush(stdout);
}

struct Infostat
{
  double min, max, sum, sumi;
  size_t nvals, nmiss, nlevs;
};

static void
infostatInit(Infostat &infostat)
{
  infostat.nvals = 0;
  infostat.nmiss = 0;
  infostat.nlevs = 0;
  infostat.min = DBL_MAX;
  infostat.max = -DBL_MAX;
  infostat.sum = 0;
  infostat.sumi = 0;
}

void *
Info(void *process)
{
  enum
  {
    E_NAME,
    E_CODE,
    E_PARAM
  };
  int fpeRaised = 0;
  int varID, levelID;
  int nrecs;
  size_t nmiss;
  size_t imiss = 0;
  char varname[CDI_MAX_NAME];
  char paramstr[32];
  char vdatestr[32], vtimestr[32];

  cdoInitialize(process);

  // clang-format off
  int INFO   = cdoOperatorAdd("info",   E_PARAM,  0, NULL);
  int INFOP  = cdoOperatorAdd("infop",  E_PARAM,  0, NULL);
  int INFON  = cdoOperatorAdd("infon",  E_NAME,   0, NULL);
  int INFOC  = cdoOperatorAdd("infoc",  E_CODE,   0, NULL);
  int XINFON = cdoOperatorAdd("xinfon", E_NAME,   0, NULL);
  int MAP    = cdoOperatorAdd("map",    E_PARAM,  0, NULL);
  // clang-format on

  UNUSED(INFO);
  UNUSED(INFOP);
  UNUSED(INFON);
  UNUSED(INFOC);
  UNUSED(XINFON);

  int operatorID = cdoOperatorID();

  int operfunc = cdoOperatorF1(operatorID);

  DateTimeList dtlist;

  for (int indf = 0; indf < cdoStreamCnt(); indf++)
    {
      int streamID = cdoStreamOpenRead(indf);

      int vlistID = cdoStreamInqVlist(streamID);
      int taxisID = vlistInqTaxis(vlistID);

      int nvars = vlistNvars(vlistID);
      if (nvars == 0) continue;

      std::vector<Infostat> infostat(nvars);

      size_t gridsizemax = vlistGridsizeMax(vlistID);
      if (vlistNumber(vlistID) != CDI_REAL) gridsizemax *= 2;

      std::vector<double> array(gridsizemax);

      int indg = 0;
      int tsID = 0;
      while ((nrecs = cdoStreamInqTimestep(streamID, tsID)))
        {
          dtlist.taxisInqTimestep(taxisID, 0);
          int64_t vdate = dtlist.getVdate(0);
          int vtime = dtlist.getVtime(0);

          date2str(vdate, vdatestr, sizeof(vdatestr));
          time2str(vtime, vtimestr, sizeof(vtimestr));

          for (varID = 0; varID < nvars; ++varID) infostatInit(infostat[varID]);

          for (int recID = 0; recID < nrecs; ++recID)
            {
              if ((tsID == 0 && recID == 0) || operatorID == MAP)
                {
                  set_text_color(stdout, BRIGHT, BLACK);
                  fprintf(stdout,
                          "%6d :       Date     Time   %s Gridsize    Miss :     Minimum        Mean     Maximum : ", -(indf + 1),
                          operatorID == XINFON ? "Nlevs" : "Level");

                  if (operfunc == E_NAME)
                    fprintf(stdout, "Parameter name");
                  else if (operfunc == E_CODE)
                    fprintf(stdout, "Code number");
                  else
                    fprintf(stdout, "Parameter ID");

                  if (cdoVerbose) fprintf(stdout, " : Extra");
                  reset_text_color(stdout);
                  fprintf(stdout, "\n");
                }

              pstreamInqRecord(streamID, &varID, &levelID);
              pstreamReadRecord(streamID, array.data(), &nmiss);

              indg = (operatorID == XINFON) ? varID + 1 : indg + 1;

              int param = vlistInqVarParam(vlistID, varID);
              int code = vlistInqVarCode(vlistID, varID);
              int gridID = vlistInqVarGrid(vlistID, varID);
              int zaxisID = vlistInqVarZaxis(vlistID, varID);
              int number = vlistInqVarNumber(vlistID, varID);
              size_t gridsize = gridInqSize(gridID);
              size_t nlevs = zaxisInqSize(zaxisID);
              double level = cdoZaxisInqLevel(zaxisID, levelID);
              double missval = vlistInqVarMissval(vlistID, varID);

              bool loutput = (operatorID != XINFON);

              if (loutput) infostatInit(infostat[varID]);

              Infostat *infostatp = &infostat[varID];
              infostatp->nlevs += 1;
              infostatp->nmiss += nmiss;

              if (nlevs == infostatp->nlevs) loutput = true;

              if (loutput)
                {
                  cdiParamToString(param, paramstr, sizeof(paramstr));

                  if (operfunc == E_NAME) vlistInqVarName(vlistID, varID, varname);

                  set_text_color(stdout, BRIGHT, BLACK);
                  fprintf(stdout, "%6d ", indg);
                  reset_text_color(stdout);
                  set_text_color(stdout, RESET, BLACK);
                  fprintf(stdout, ":");
                  reset_text_color(stdout);

                  set_text_color(stdout, RESET, MAGENTA);
                  fprintf(stdout, "%s %s ", vdatestr, vtimestr);
                  reset_text_color(stdout);

                  set_text_color(stdout, RESET, GREEN);
                  if (operatorID == XINFON)
                    fprintf(stdout, "%7zu ", nlevs);
                  else
                    fprintf(stdout, "%7g ", level);

                  fprintf(stdout, "%8zu %7zu ", gridsize, infostatp->nmiss);

                  set_text_color(stdout, RESET, BLACK);
                  fprintf(stdout, ":");
                  reset_text_color(stdout);

                  set_text_color(stdout, RESET, BLUE);
                }

              if (number == CDI_REAL)
                {
                  fpeRaised = 0;

                  if (infostatp->nmiss > 0)
                    {
                      size_t nvals
                          = arrayMinMaxSumMV(gridsize, array.data(), missval, &infostatp->min, &infostatp->max, &infostatp->sum);
                      imiss = gridsize - nvals;
                      infostatp->nvals += nvals;
                    }
                  else if (gridsize == 1)
                    {
                      if (infostatp->nvals == 0)
                        infostatp->sum = array[0];
                      else
                        infostatp->sum += array[0];
                      infostatp->nvals += 1;
                    }
                  else
                    {
                      arrayMinMaxSum(gridsize, array.data(), &infostatp->min, &infostatp->max, &infostatp->sum);
                      infostatp->nvals += gridsize;
                    }

                  if (loutput)
                    {
                      if (infostatp->nvals)
                        {
                          if (infostatp->nvals == 1)
                            {
                              fprintf(stdout, "            %#12.5g            ", infostatp->sum);
                            }
                          else
                            {
                              double mean = infostatp->sum / (double) infostatp->nvals;
                              fprintf(stdout, "%#12.5g%#12.5g%#12.5g", infostatp->min, mean, infostatp->max);
                            }
                        }
                      else
                        {
                          fprintf(stdout, "                     nan            ");
                        }
                    }
                }
              else
                {
                  size_t nvals = 0;
                  for (size_t i = 0; i < gridsize; i++)
                    {
                      if (!DBL_IS_EQUAL(array[i * 2], missval) && !DBL_IS_EQUAL(array[i * 2 + 1], missval))
                        {
                          infostatp->sum += array[i * 2];
                          infostatp->sumi += array[i * 2 + 1];
                          nvals++;
                        }
                    }
                  fpeRaised = 0;

                  imiss = gridsize - nvals;
                  infostatp->nvals += nvals;

                  if (loutput)
                    {
                      double arrmean_r = 0, arrmean_i = 0;
                      if (infostatp->nvals > 0) arrmean_r = infostatp->sum / infostatp->nvals;
                      if (infostatp->nvals > 0) arrmean_i = infostatp->sumi / infostatp->nvals;
                      fprintf(stdout, "   -  (%#12.5g,%#12.5g)  -", arrmean_r, arrmean_i);
                    }
                }

              if (loutput)
                {
                  reset_text_color(stdout);

                  set_text_color(stdout, RESET, BLACK);
                  fprintf(stdout, " : ");
                  reset_text_color(stdout);

                  set_text_color(stdout, BRIGHT, GREEN);
                  if (operfunc == E_NAME)
                    fprintf(stdout, "%-14s", varname);
                  else if (operfunc == E_CODE)
                    fprintf(stdout, "%4d   ", code);
                  else
                    fprintf(stdout, "%-14s", paramstr);
                  reset_text_color(stdout);

                  if (cdoVerbose)
                    {
                      char varextra[CDI_MAX_NAME];
                      vlistInqVarExtra(vlistID, varID, varextra);
                      fprintf(stdout, " : %s", varextra);
                    }

                  fprintf(stdout, "\n");
                }

              if (imiss != nmiss && nmiss > 0) cdoPrint("Found %zu of %zu missing values!", imiss, nmiss);

              if (fpeRaised > 0) cdoWarning("floating-point exception reported: %s!", fpe_errstr(fpeRaised));

              if (operatorID == MAP)
                {
                  size_t nlon = gridInqXsize(gridID);
                  size_t nlat = gridInqYsize(gridID);

                  if (gridInqType(gridID) == GRID_GAUSSIAN || gridInqType(gridID) == GRID_LONLAT
                      || gridInqType(gridID) == GRID_CURVILINEAR
                      || (gridInqType(gridID) == GRID_GENERIC && nlon * nlat == gridInqSize(gridID) && nlon < 1024))
                    {
                      printMap(nlon, nlat, array.data(), missval, infostatp->min, infostatp->max);
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
