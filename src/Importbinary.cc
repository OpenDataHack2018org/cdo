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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "grid.h"

#include "gradsdeslib.h"

static void
get_dim_vals(dsets_t *pfi, double *vals, int dimlen, int dim)
{
  gadouble (*conv)(gadouble *, gadouble);
  gadouble *cvals;
  int i;

  assert(dimlen == pfi->dnum[dim]);

  if (pfi->linear[dim] == 0)
    {
      for (i = 0; i < dimlen; ++i)
        {
          vals[i] = pfi->grvals[dim][i + 1];
          /* printf("%d %g\n", i, vals[i]); */
        }
    }
  else if (pfi->linear[dim] == 1)
    {
      /*
      for ( i = 0; i < 3; ++i )
        printf("%d %g %g\n", i, pfi->grvals[dim][i] , pfi->abvals[dim][i]);
      */
      conv = pfi->gr2ab[dim];
      cvals = pfi->grvals[dim];
      for (i = 0; i < dimlen; ++i)
        {
          vals[i] = conv(cvals, i + 1);
          /* printf("%d %g\n", i, vals[i]); */
        }
    }
}

static void
rev_vals(double *vals, int n)
{
  double dum;

  for (int i = 0; i < n / 2; ++i)
    {
      dum = vals[i];
      vals[i] = vals[n - 1 - i];
      vals[n - 1 - i] = dum;
    }
}

static int
y_is_gauss(double *gridyvals, int ysize)
{
  int lgauss = FALSE;
  int i;

  if (ysize > 2)
    {
      double *yvals, *yw;
      yvals = (double *) Malloc(ysize * sizeof(double));
      yw = (double *) Malloc(ysize * sizeof(double));
      gaussaw(yvals, yw, (size_t) ysize);
      Free(yw);
      for (i = 0; i < (int) ysize; i++) yvals[i] = asin(yvals[i]) / M_PI * 180.0;

      for (i = 0; i < (int) ysize; i++)
        if (fabs(yvals[i] - gridyvals[i]) > ((yvals[0] - yvals[1]) / 500)) break;

      if (i == (int) ysize) lgauss = TRUE;

      /* check S->N */
      if (lgauss == FALSE)
        {
          for (i = 0; i < (int) ysize; i++)
            if (fabs(yvals[i] - gridyvals[ysize - i - 1]) > ((yvals[0] - yvals[1]) / 500)) break;

          if (i == (int) ysize) lgauss = TRUE;
        }

      Free(yvals);
    }

  return lgauss;
}

static int
define_grid(dsets_t *pfi)
{
  int gridID, gridtype;
  int lgauss = FALSE;

  int nx = pfi->dnum[0];
  int ny = pfi->dnum[1];

  double *xvals = (double *) Malloc(nx * sizeof(double));
  double *yvals = (double *) Malloc(ny * sizeof(double));

  get_dim_vals(pfi, xvals, nx, 0);
  get_dim_vals(pfi, yvals, ny, 1);

  if (pfi->yrflg) rev_vals(yvals, ny);

  if (pfi->linear[1] == 0) lgauss = y_is_gauss(yvals, (size_t) ny);

  if (lgauss)
    gridtype = GRID_GAUSSIAN;
  else
    gridtype = GRID_LONLAT;

  gridID = gridCreate(gridtype, nx * ny);
  gridDefXsize(gridID, nx);
  gridDefYsize(gridID, ny);

  gridDefXvals(gridID, xvals);
  gridDefYvals(gridID, yvals);

  Free(xvals);
  Free(yvals);

  return gridID;
}

static int
define_level(dsets_t *pfi, int nlev)
{
  int zaxisID = -1;
  int nz = pfi->dnum[2];

  if (nz)
    {
      double *zvals = (double *) Malloc(nz * sizeof(double));

      get_dim_vals(pfi, zvals, nz, 2);

      if (nz == 1 && IS_EQUAL(zvals[0], 0))
        zaxisID = zaxisCreate(ZAXIS_SURFACE, nz);
      else
        {
          if (nlev > 0 && nlev < nz) nz = nlev;
          if (pfi->zrflg) rev_vals(zvals, nz);
          zaxisID = zaxisCreate(ZAXIS_GENERIC, nz);
        }
      zaxisDefLevels(zaxisID, zvals);

      Free(zvals);
    }
  else
    {
      double level = 0;
      nz = 1;

      zaxisID = zaxisCreate(ZAXIS_SURFACE, nz);
      zaxisDefLevels(zaxisID, &level);
    }

  return zaxisID;
}

void *
Importbinary(void *process)
{
  size_t i;
  size_t nmiss = 0, n_nan;
  int ivar;
  int varID = -1, levelID, tsID;
  size_t gridsize;
  int datatype;
  dsets_t pfi;
  int64_t vdate;
  int vtime;
  int tcur, told, fnum;
  int tmin = 0, tmax = 0;
  char *ch = NULL;
  int nlevels;
  int e, flag;
  size_t rc, recsize;
  int recoffset;
  char *rec = NULL;
  struct dt dtim, dtimi;
  double missval;
  double fmin, fmax;
  double *array;
  double sfclevel = 0;
  char vdatestr[32], vtimestr[32];

  cdoInitialize(process);

  dsets_init(&pfi);

  int status = read_gradsdes((char *) cdoGetStreamName(0).c_str(), &pfi);
  if (cdoVerbose) fprintf(stderr, "status %d\n", status);
  // if ( status ) cdoAbort("Open failed on %s!", pfi.name);
  if (status) cdoAbort("Open failed!");

  int nrecs = pfi.trecs;
  int nvars = pfi.vnum;
  struct gavar *pvar = pfi.pvar1;

  if (nvars == 0) cdoAbort("No variables found!");

  int gridID = define_grid(&pfi);
  int zaxisID = define_level(&pfi, 0);

  int zaxisIDsfc = zaxisCreate(ZAXIS_SURFACE, 1);
  zaxisDefLevels(zaxisIDsfc, &sfclevel);

  int vlistID = vlistCreate();

  int *var_zaxisID = (int *) Malloc(nvars * sizeof(int));
  int *var_dfrm = (int *) Malloc(nrecs * sizeof(int));
  std::vector<RecordInfo> recinfo(nrecs);

  int recID = 0;
  for (ivar = 0; ivar < nvars; ++ivar)
    {
      /*
      if ( cdoVerbose )
        fprintf(stderr, "1:%s 2:%s %d %d %d %d 3:%s %d \n",
                pvar->abbrv, pvar->longnm, pvar->offset, pvar->recoff,
      pvar->levels, pvar->nvardims, pvar->varnm, pvar->var_t);
      */
      nlevels = pvar->levels;

      if (nlevels == 0)
        {
          nlevels = 1;
          varID = vlistDefVar(vlistID, gridID, zaxisIDsfc, TIME_VARYING);
        }
      else
        {
          if (nlevels > zaxisInqSize(zaxisID))
            cdoAbort("Variable %s has too many number of levels!", pvar->abbrv);
          else if (nlevels < zaxisInqSize(zaxisID))
            {
              int vid, zid = -1, nlev;
              for (vid = 0; vid < ivar; ++vid)
                {
                  zid = var_zaxisID[vid];
                  nlev = zaxisInqSize(zid);
                  if (nlev == nlevels) break;
                }

              if (vid == ivar) zid = define_level(&pfi, nlevels);
              varID = vlistDefVar(vlistID, gridID, zid, TIME_VARYING);
            }
          else
            varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARYING);
        }

      var_zaxisID[varID] = vlistInqVarZaxis(vlistID, varID);

      vlistDefVarName(vlistID, varID, pvar->abbrv);
      {
        char *longname = pvar->varnm;
        size_t len = strlen(longname);
        if (longname[0] == '\'' && longname[len - 1] == '\'')
          {
            longname[len - 1] = 0;
            longname++;
          }
        if (longname[0] == '\t') longname++;
        vlistDefVarLongname(vlistID, varID, longname);
      }

      missval = pfi.undef;
      datatype = CDI_DATATYPE_FLT32;

      if (pvar->dfrm == 1)
        {
          datatype = CDI_DATATYPE_UINT8;
          if (missval < 0 || missval > 255) missval = 255;
        }
      else if (pvar->dfrm == 2)
        {
          datatype = CDI_DATATYPE_UINT16;
          if (missval < 0 || missval > 65535) missval = 65535;
        }
      else if (pvar->dfrm == -2)
        {
          datatype = CDI_DATATYPE_INT16;
          if (missval < -32768 || missval > 32767) missval = -32768;
        }
      else if (pvar->dfrm == 4)
        {
          datatype = CDI_DATATYPE_INT32;
          if (missval < -2147483648 || missval > 2147483647) missval = -2147483646;
        }
      else if (pfi.flt64)
        datatype = CDI_DATATYPE_FLT64;

      vlistDefVarDatatype(vlistID, varID, datatype);
      vlistDefVarMissval(vlistID, varID, missval);

      for (levelID = 0; levelID < nlevels; ++levelID)
        {
          if (recID >= nrecs) cdoAbort("Internal problem with number of records!");
          recinfo[recID].varID = varID;
          recinfo[recID].levelID = levelID;
          var_dfrm[recID] = pvar->dfrm;
          recID++;
        }

      pvar++;
    }

  int calendar = CALENDAR_STANDARD;

  int taxisID = taxisCreate(TAXIS_RELATIVE);

  taxisDefCalendar(taxisID, calendar);

  vlistDefTaxis(vlistID, taxisID);

  int streamID = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());

  pstreamDefVlist(streamID, vlistID);

  gridsize = pfi.dnum[0] * pfi.dnum[1];
  if (pfi.flt64)
    recoffset = pfi.xyhdr * 8;
  else
    recoffset = pfi.xyhdr * 4;

  if (pfi.seqflg) recoffset += 4;

  // recsize = pfi.gsiz*4;
  recsize = pfi.gsiz * 8;
  rec = (char *) Malloc(recsize);

  array = (double *) Malloc(gridsize * sizeof(double));

  /*
  if (pfi.tmplat)
    for ( i = 0; i <  pfi.dnum[3]; ++i )
      printf("%d %d\n", i, pfi.fnums[i]);
  */

  pfi.infile = NULL;
  tcur = 0;
  e = 1;
  while (1)
    { /* loop over all times for this ensemble */
      if (pfi.tmplat)
        {
          /* make sure no file is open */
          if (pfi.infile != NULL)
            {
              fclose(pfi.infile);
              pfi.infile = NULL;
            }
          /* advance to first valid time step for this ensemble */
          if (tcur == 0)
            {
              told = 0;
              tcur = 1;
              while (pfi.fnums[tcur - 1] == -1) tcur++;
            }
          else
            { /* tcur!=0 */
              told = pfi.fnums[tcur - 1];
              /* increment time step until fnums changes */
              while (told == pfi.fnums[tcur - 1] && tcur <= pfi.dnum[3])
                {
                  tcur++;
                  if (tcur > pfi.dnum[3]) break;
                }
            }

          /* make sure we haven't advanced past end of time axis */
          if (tcur > pfi.dnum[3]) break;

          /* check if we're past all valid time steps for this ensemble */
          if ((told != -1) && (pfi.fnums[tcur - 1] == -1)) break;

          /* Find the range of t indexes that have the same fnums value.
             These are the times that are contained in this particular file */
          tmin = tcur;
          tmax = tcur - 1;
          fnum = pfi.fnums[tcur - 1];
          if (fnum != -1)
            {
              while (fnum == pfi.fnums[tmax])
                {
                  tmax++;
                  if (tmax == pfi.dnum[3]) break;
                }
              gr2t(pfi.grvals[3], (gadouble) tcur, &dtim);
              gr2t(pfi.grvals[3], (gadouble) 1, &dtimi);
              ch = gafndt(pfi.name, &dtim, &dtimi, pfi.abvals[3], pfi.pchsub1, NULL, tcur, e, &flag);
              if (ch == NULL) cdoAbort("Couldn't determine data file name for e=%d t=%d!", e, tcur);
            }
        }
      else
        {
          /* Data set is not templated */
          ch = pfi.name;
          tmin = 1;
          tmax = pfi.dnum[3];
        }

      /* Open this file and position to start of first record */
      if (cdoVerbose) cdoPrint("Opening file: %s", ch);
      pfi.infile = fopen(ch, "rb");
      if (pfi.infile == NULL)
        {
          if (pfi.tmplat)
            {
              cdoWarning("Could not open file: %s", ch);
              break;
            }
          else
            {
              cdoAbort("Could not open file: %s", ch);
            }
        }
      if (pfi.tmplat) gree(ch, "312");

      /* file header */
      if (pfi.fhdr > 0) fseeko(pfi.infile, pfi.fhdr, SEEK_SET);

      /* Get file size */
      /*
      fseeko(pfi.infile,0L,2);
      flen = ftello(pfi.infile);

      printf("flen %d tsiz %d\n", flen, pfi.tsiz);

      fseeko (pfi.infile,0,0);
      */
      for (tsID = tmin - 1; tsID < tmax; ++tsID)
        {
          gr2t(pfi.grvals[3], (gadouble)(tsID + 1), &dtim);
          vdate = cdiEncodeDate(dtim.yr, dtim.mo, dtim.dy);
          vtime = cdiEncodeTime(dtim.hr, dtim.mn, 0);

          date2str(vdate, vdatestr, sizeof(vdatestr));
          time2str(vtime, vtimestr, sizeof(vtimestr));

          if (cdoVerbose) cdoPrint(" Reading timestep: %3d %s %s", tsID + 1, vdatestr, vtimestr);

          taxisDefVdate(taxisID, vdate);
          taxisDefVtime(taxisID, vtime);
          pstreamDefTimestep(streamID, tsID);

          for (int recID = 0; recID < nrecs; ++recID)
            {
              /* record size depends on data type */
              if (var_dfrm[recID] == 1)
                {
                  recsize = pfi.gsiz;
                }
              else if ((var_dfrm[recID] == 2) || (var_dfrm[recID] == -2))
                {
                  recsize = pfi.gsiz * 2;
                }
              else
                {
                  recsize = pfi.flt64 ? pfi.gsiz * 8 : pfi.gsiz * 4;
                }

              rc = fread(rec, 1, recsize, pfi.infile);
              if (rc < recsize) cdoAbort("I/O error reading record=%d of timestep=%d!", recID + 1, tsID + 1);

              /* convert */
              if (var_dfrm[recID] == 1)
                {
                  unsigned char *carray = (unsigned char *) (rec + recoffset);
                  for (i = 0; i < gridsize; ++i) array[i] = (double) carray[i];
                }
              else if (var_dfrm[recID] == 2)
                {
                  unsigned short *sarray = (unsigned short *) (rec + recoffset);
                  if (pfi.bswap) gabswp2(sarray, gridsize);
                  for (i = 0; i < gridsize; ++i) array[i] = (double) sarray[i];
                }
              else if (var_dfrm[recID] == -2)
                {
                  short *sarray = (short *) (rec + recoffset);
                  if (pfi.bswap) gabswp2(sarray, gridsize);
                  for (i = 0; i < gridsize; ++i) array[i] = (double) sarray[i];
                }
              else if (var_dfrm[recID] == 4)
                {
                  int *iarray = (int *) (rec + recoffset);
                  if (pfi.bswap) gabswp(iarray, gridsize);
                  for (i = 0; i < gridsize; ++i) array[i] = (double) iarray[i];
                }
              else
                {
                  if (pfi.flt64)
                    {
                      double *darray = (double *) (rec + recoffset);
                      if (pfi.bswap) gabswp(darray, gridsize);
                      for (i = 0; i < gridsize; ++i) array[i] = darray[i];
                    }
                  else
                    {
                      float *farray = (float *) (rec + recoffset);
                      if (pfi.bswap) gabswp(farray, gridsize);
                      for (i = 0; i < gridsize; ++i) array[i] = (double) farray[i];
                    }
                }

              fmin = 1.e99;
              fmax = -1.e99;
              nmiss = 0;
              n_nan = 0;
              for (i = 0; i < gridsize; ++i)
                {
                  if (array[i] > pfi.ulow && array[i] < pfi.uhi)
                    {
                      array[i] = pfi.undef;
                      nmiss++;
                    }
                  else if (DBL_IS_NAN(array[i]))
                    {
                      array[i] = pfi.undef;
                      nmiss++;
                      n_nan++;
                    }
                  else
                    {
                      if (array[i] < fmin) fmin = array[i];
                      if (array[i] > fmax) fmax = array[i];
                    }
                }
              /*
              if ( cdoVerbose )
                printf("%3d %4d %3d %6zu %6zu %12.5g %12.5g\n", tsID, recID,
              recoffset, nmiss, n_nan, fmin, fmax);
              */
              varID = recinfo[recID].varID;
              levelID = recinfo[recID].levelID;
              pstreamDefRecord(streamID, varID, levelID);
              pstreamWriteRecord(streamID, array, nmiss);
            }
        }

      /* break out if not templating */
      if (!pfi.tmplat) break;

    } /* end of while (1) loop */

  processDefVarNum(vlistNvars(vlistID));

  pstreamClose(streamID);

  vlistDestroy(vlistID);
  gridDestroy(gridID);
  zaxisDestroy(zaxisID);
  taxisDestroy(taxisID);

  Free(array);
  Free(rec);

  if (var_zaxisID) Free(var_zaxisID);
  if (var_dfrm) Free(var_dfrm);

  cdoFinish();

  return 0;
}
