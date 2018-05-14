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

      Sort sortcode  Sort by code number
*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"

typedef struct
{
  int levelID;
  size_t nmiss;
  double level;
} levinfo_t;

typedef struct
{
  int varID;
  int nlevs;
  int code;
  char param[CDI_MAX_NAME];
  char name[CDI_MAX_NAME];
  levinfo_t *levInfo;
} varinfo_t;

static int
cmpvarcode(const void *a, const void *b)
{
  const int x = ((const varinfo_t *) a)->code;
  const int y = ((const varinfo_t *) b)->code;
  return ((x > y) - (x < y)) * 2 + (x > y) - (x < y);
}

static int
cmpvarparam(const void *a, const void *b)
{
  const char *x = ((const varinfo_t *) a)->param;
  const char *y = ((const varinfo_t *) b)->param;
  return strcmp(x, y);
}

static int
cmpvarname(const void *a, const void *b)
{
  const char *x = ((const varinfo_t *) a)->name;
  const char *y = ((const varinfo_t *) b)->name;
  return strcmp(x, y);
}

static int
cmpvarlevel(const void *a, const void *b)
{
  const double x = ((const levinfo_t *) a)->level;
  const double y = ((const levinfo_t *) b)->level;
  return ((x > y) - (x < y)) * 2 + (x > y) - (x < y);
}

static int
cmpvarlevelrev(const void *a, const void *b)
{
  const double x = ((const levinfo_t *) a)->level;
  const double y = ((const levinfo_t *) b)->level;
  return -1 * (((x > y) - (x < y)) * 2 + (x > y) - (x < y));
}

static void
setNmiss(int varID, int levelID, int nvars, varinfo_t *varInfo, size_t nmiss)
{
  int vindex, lindex;

  for (vindex = 0; vindex < nvars; vindex++)
    if (varInfo[vindex].varID == varID) break;

  if (vindex == nvars) cdoAbort("Internal problem; varID not found!");

  int nlevs = varInfo[vindex].nlevs;
  for (lindex = 0; lindex < nlevs; lindex++)
    if (varInfo[vindex].levInfo[lindex].levelID == levelID) break;

  if (lindex == nlevs) cdoAbort("Internal problem; levelID not found!");

  varInfo[vindex].levInfo[lindex].nmiss = nmiss;
}

void
paramToStringLong(int param, char *paramstr, int maxlen)
{
  int len;

  int dis, cat, num;
  cdiDecodeParam(param, &num, &cat, &dis);

  size_t umaxlen = maxlen >= 0 ? (unsigned) maxlen : 0U;
  if (dis == 255 && (cat == 255 || cat == 0))
    len = snprintf(paramstr, umaxlen, "%03d", num);
  else if (dis == 255)
    len = snprintf(paramstr, umaxlen, "%03d.%03d", num, cat);
  else
    len = snprintf(paramstr, umaxlen, "%03d.%03d.%03d", num, cat, dis);

  if (len >= maxlen || len < 0) fprintf(stderr, "Internal problem (%s): size of input string is too small!\n", __func__);
}

void *
Sort(void *process)
{
  int varID, levelID, zaxisID;
  int vindex, lindex;
  int nrecs, nlevs, offset;
  size_t gridsize;
  size_t nmiss;
  double *single;
  int (*cmpvarlev)(const void *, const void *) = cmpvarlevel;

  cdoInitialize(process);

  // clang-format off
  int SORTCODE  = cdoOperatorAdd("sortcode",  0, 0, NULL);
  int SORTPARAM = cdoOperatorAdd("sortparam", 0, 0, NULL);
  int SORTNAME  = cdoOperatorAdd("sortname",  0, 0, NULL);
  int SORTLEVEL = cdoOperatorAdd("sortlevel", 0, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  if (operatorArgc() > 1) cdoAbort("Too many arguments!");

  if (operatorID == SORTLEVEL && operatorArgc() == 1)
    {
      int iarg = parameter2int(operatorArgv()[0]);
      if (iarg < 0) cmpvarlev = cmpvarlevelrev;
    }

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID1 = cdoStreamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);
  /*
  if ( operatorID == SORTCODE )
      vlistSortCode(vlistID2);
   else if ( operatorID == SORTNAME )
      ;
   else if ( operatorID == SORTLEVEL )
      ;
  */

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  int nvars = vlistNvars(vlistID1);

  varinfo_t *varInfo = (varinfo_t *) Malloc(nvars * sizeof(varinfo_t));
  for (varID = 0; varID < nvars; ++varID)
    {
      nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      varInfo[varID].nlevs = nlevs;
      varInfo[varID].levInfo = (levinfo_t *) Malloc(nlevs * sizeof(levinfo_t));
    }

  double **vardata = (double **) Malloc(nvars * sizeof(double *));

  for (varID = 0; varID < nvars; varID++)
    {
      gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
      nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, varID));
      vardata[varID] = (double *) Malloc(gridsize * nlevs * sizeof(double));
    }

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);

          if (tsID == 0)
            {
              varInfo[varID].varID = varID;
              varInfo[varID].code = vlistInqVarCode(vlistID1, varID);
              int iparam = vlistInqVarParam(vlistID1, varID);
              paramToStringLong(iparam, varInfo[varID].param, sizeof(varInfo[varID].param));
              vlistInqVarName(vlistID1, varID, varInfo[varID].name);
              zaxisID = vlistInqVarZaxis(vlistID1, varID);
              varInfo[varID].levInfo[levelID].levelID = levelID;
              varInfo[varID].levInfo[levelID].level = cdoZaxisInqLevel(zaxisID, levelID);
            }

          gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
          offset = gridsize * levelID;
          single = vardata[varID] + offset;

          pstreamReadRecord(streamID1, single, &nmiss);

          setNmiss(varID, levelID, nvars, varInfo, nmiss);
          // varInfo[varID].levInfo[levelID].nmiss = nmiss;
        }

      if (tsID == 0)
        {
          if (cdoVerbose)
            for (vindex = 0; vindex < nvars; vindex++)
              {
                nlevs = varInfo[vindex].nlevs;
                for (lindex = 0; lindex < nlevs; ++lindex)
                  printf("sort in: %d %s %d %d %g\n", vindex, varInfo[vindex].name, varInfo[vindex].code, varInfo[vindex].nlevs,
                         varInfo[vindex].levInfo[lindex].level);
              }

          if (operatorID == SORTCODE)
            qsort(varInfo, nvars, sizeof(varinfo_t), cmpvarcode);
          else if (operatorID == SORTPARAM)
            qsort(varInfo, nvars, sizeof(varinfo_t), cmpvarparam);
          else if (operatorID == SORTNAME)
            qsort(varInfo, nvars, sizeof(varinfo_t), cmpvarname);
          else if (operatorID == SORTLEVEL)
            {
              for (vindex = 0; vindex < nvars; vindex++)
                {
                  nlevs = varInfo[vindex].nlevs;
                  qsort(varInfo[vindex].levInfo, nlevs, sizeof(levinfo_t), cmpvarlev);
                }
            }

          if (cdoVerbose)
            for (vindex = 0; vindex < nvars; vindex++)
              {
                nlevs = varInfo[vindex].nlevs;
                for (lindex = 0; lindex < nlevs; ++lindex)
                  printf("sort out: %d %s %d %d %g\n", vindex, varInfo[vindex].name, varInfo[vindex].code, varInfo[vindex].nlevs,
                         varInfo[vindex].levInfo[lindex].level);
              }
        }

      for (vindex = 0; vindex < nvars; vindex++)
        {
          varID = varInfo[vindex].varID;
          nlevs = varInfo[vindex].nlevs;
          for (lindex = 0; lindex < nlevs; ++lindex)
            {
              levelID = varInfo[vindex].levInfo[lindex].levelID;
              nmiss = varInfo[vindex].levInfo[lindex].nmiss;

              if (tsID == 0 || vlistInqVarTimetype(vlistID1, varID) != TIME_CONSTANT)
                {
                  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));
                  offset = gridsize * levelID;
                  single = vardata[varID] + offset;

                  pstreamDefRecord(streamID2, varID, levelID);
                  pstreamWriteRecord(streamID2, single, nmiss);
                }
            }
        }

      tsID++;
    }

  pstreamClose(streamID1);
  pstreamClose(streamID2);

  for (varID = 0; varID < nvars; varID++) Free(vardata[varID]);
  Free(vardata);

  for (vindex = 0; vindex < nvars; vindex++) Free(varInfo[vindex].levInfo);
  Free(varInfo);

  cdoFinish();

  return 0;
}
