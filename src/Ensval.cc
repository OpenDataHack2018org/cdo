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
   Ensval       enscrps          Ensemble cumulative ranked probability score &
   decomposition Ensval       ensbrs           Ensemble Brier score &
   decomposition

   The implementation of the decomposition and score calculation
   as carried out in this routine follows the paper
     Hans Hersbach (2000): Decomposition of the Continuous Ranked Probability
     Score for Ensemble Prediction Systems, in: Weather and Forecasting (15)
     pp. 559-570
*/

#include <cmath>
#include <cdi.h>

#include "statistic.h"
#include "cdo_int.h"
#include "pstream_int.h"
#include "merge_sort2.h"
#include "util_files.h"

enum OPERTYPE
{
  CRPS,
  BRS
};

enum RESTYPE_BRS
{
  BRS_RES,
  BRS_RELI,
  BRS_RESOL,
  BRS_UNCTY
};
enum RESTYPE_CRPS
{
  CRPS_RES,
  CRPS_RELI,
  CRPS_POT
};

void *
Ensval(void *process)
{
  int i, k;
  int nrecs = 0, nrecs0, nostreams = 0, ngrids;
  size_t nmiss;
  int levelID, varID;
  size_t gridsize = 0;
  int vlistID;
  int gridID = -1;
  int have_miss = 0;
  int stream, streamID = 0;
  // int xsize,ysize;
  double missval = 0;
  double *alpha, *beta, *alpha_weights, *beta_weights;
  double *brs_g, *brs_o, *brs_g_weights, *brs_o_weights;
  double xval = 0;
  double yval = 0;
  double xa, *x;
  double *weights, sum_weights = 0;
  double crps_reli, crps_pot, crps;
  double heavyside0, heavysideN;
  double brs_reli, brs_resol, brs_uncty, brs_thresh = 0;
  double g, o, p;

  int fileID;
  char file_suffix[32];
  char type_suffix[10];

  struct ens_file_t
  {
    int streamID;
    int vlistID;
    double *array;
  };

  // INITIALIZE POINTERS
  alpha = NULL;
  beta = NULL;
  alpha_weights = NULL;
  beta_weights = NULL;
  brs_g = NULL;
  brs_o = NULL;
  brs_g_weights = NULL;
  brs_o_weights = NULL;
  x = NULL;
  weights = NULL;
  // int vlistCheck, gridsizeCheck;

  cdoInitialize(process);

  cdoOperatorAdd("enscrps", CRPS, 0, NULL);
  cdoOperatorAdd("ensbrs", BRS, 0, NULL);

  int operatorID = cdoOperatorID();
  int operfunc = cdoOperatorF1(operatorID);

  int nfiles = cdoStreamCnt() - 1;
  int nens = nfiles - 1;

  if (operfunc == CRPS)
    {
      nostreams = 3;
    }
  else if (operfunc == BRS)
    {
      operatorInputArg("Threshold for Brier score?");
      operatorCheckArgc(1);
      brs_thresh = parameter2double(operatorArgv()[0]);
      nostreams = 4;

      fprintf(stderr, "brs_thres %10.6f\n", brs_thresh);
    }

  // allocate array to hold results
  double *r = (double *) Malloc(nostreams * sizeof(double));

  // one stream for each value of the decomposition
  int *streamID2 = (int *) Malloc(nostreams * sizeof(int));
  int *vlistID2 = (int *) Malloc(nostreams * sizeof(int));
  int *taxisID2 = (int *) Malloc(nostreams * sizeof(int));
  int *zaxisID2 = (int *) Malloc(nostreams * sizeof(int));

  double *val = (double *) Calloc(nfiles, sizeof(double));

  if (operfunc == CRPS)
    {
      alpha = (double *) Calloc(nens + 1, sizeof(double));
      beta = (double *) Calloc(nens + 1, sizeof(double));
      alpha_weights = (double *) Calloc(nens + 1, sizeof(double));
      beta_weights = (double *) Calloc(nens + 1, sizeof(double));
    }
  else if (operfunc == BRS)
    {
      brs_g = (double *) Calloc(nens + 1, sizeof(double));
      brs_o = (double *) Calloc(nens + 1, sizeof(double));
    }
  if (cdoVerbose) cdoPrint("Ensemble over %d files (Ensstat5).", nfiles - 1);

  ens_file_t *ef = (ens_file_t *) Malloc(nfiles * sizeof(ens_file_t));

  for (fileID = 0; fileID < nfiles; fileID++)
    {
      streamID = cdoStreamOpenRead(cdoStreamName(fileID));
      vlistID = cdoStreamInqVlist(streamID);

      ef[fileID].streamID = streamID;
      ef[fileID].vlistID = vlistID;
      ef[fileID].array = NULL;
    }

  int streamID1 = ef[0].streamID;

  if (cdoVerbose) cdoPrint("Opened %i Input Files for Ensemble Operator", nfiles);

  /* check for identical contents of all ensemble members */
  int nvars = vlistNvars(ef[0].vlistID);
  if (cdoVerbose) cdoPrint("nvars %i", nvars);

  for (fileID = 1; fileID < nfiles; fileID++) vlistCompare(ef[0].vlistID, ef[fileID].vlistID, CMP_ALL);

  int vlistID1 = ef[0].vlistID;
  int taxisID1 = vlistInqTaxis(vlistID1);
  int zaxisID1 = vlistInqVarZaxis(vlistID1, 0);

  int gridID2 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID2, 1);
  gridDefYsize(gridID2, 1);
  gridDefXvals(gridID2, &xval);
  gridDefYvals(gridID2, &yval);

  const char *ofilebase = cdoGetStreamName(nfiles).c_str();

  const char *refname = cdoGetObase();
  memset(file_suffix, 0, sizeof(file_suffix));
  cdoGenFileSuffix(file_suffix, sizeof(file_suffix), pstreamInqFiletype(streamID1), vlistID1, refname);

  for (stream = 0; stream < nostreams; stream++)
    {
      int namelen = strlen(ofilebase) + 9 /*type_suffix*/
                    + 32                  /*file_suffix*/
                    + 3 /*separating dots and EOS*/;

      switch (operfunc)
        {
        case CRPS:
          switch (stream)
            {
            case 0: sprintf(type_suffix, "crps"); break;
            case 1: sprintf(type_suffix, "crps_reli"); break;
            case 2: sprintf(type_suffix, "crps_pot"); break;
            }
          break;
        case BRS:
          switch (stream)
            {
            case 0: sprintf(type_suffix, "brs"); break;
            case 1: sprintf(type_suffix, "brs_reli"); break;
            case 2: sprintf(type_suffix, "brs_reso"); break;
            case 3: sprintf(type_suffix, "brs_unct"); break;
            }
          break;
        }

      char *ofilename = (char *) Calloc(namelen, sizeof(char));

      sprintf(ofilename, "%s.%s%s", ofilebase, type_suffix, file_suffix);
      // fprintf(stderr, "StreamID %i: %s\n", stream, ofilename);

      if (!cdoOverwriteMode && fileExists(ofilename) && !userFileOverwrite(ofilename))
        cdoAbort("Outputfile %s already exists!", ofilename);

      streamID2[stream] = cdoStreamOpenWrite(ofilename, cdoFiletype());

      Free(ofilename);

      zaxisID2[stream] = zaxisDuplicate(zaxisID1);
      taxisID2[stream] = taxisDuplicate(taxisID1);
      vlistID2[stream] = vlistDuplicate(vlistID1);

      ngrids = vlistNgrids(vlistID2[stream]);
      // fprintf(stderr,"ngrids %i\n",ngrids);
      for (i = 0; i < ngrids; i++) vlistChangeGridIndex(vlistID2[stream], i, gridID2);

      vlistDefTaxis(vlistID2[stream], taxisID2[stream]);
      pstreamDefVlist(streamID2[stream], vlistID2[stream]);

      // vlistCheck = cdoStreamInqVlist(streamID2[stream]);
      // gridsizeCheck = vlistGridsizeMax(vlistCheck);

      // fprintf(stderr,"stream %i vlist %3i gridsize
      // %4i\n",stream,vlistCheck,gridsizeCheck);
    }

  if (cdoVerbose) cdoPrint(" sum_weights %10.6f", sum_weights);

  int tsID = 0;
  do
    {
      nrecs0 = cdoStreamInqTimestep(ef[0].streamID, tsID);
      for (fileID = 1; fileID < nfiles; fileID++)
        {
          streamID = ef[fileID].streamID;
          nrecs = cdoStreamInqTimestep(streamID, tsID);
          if (nrecs != nrecs0)
            cdoAbort("Number of records at time step %d of %s and %s differ!", tsID + 1, cdoGetStreamName(0).c_str(),
                     cdoGetStreamName(fileID).c_str());
        }

      for (stream = 0; stream < nostreams; stream++)
        {
          taxisCopyTimestep(taxisID2[stream], taxisID1);
          if (nrecs0 > 0) pstreamDefTimestep(streamID2[stream], tsID);
        }

      for (int recID = 0; recID < nrecs0; recID++)
        {
          for (fileID = 0; fileID < nfiles; fileID++)
            {
              pstreamInqRecord(fileID, &varID, &levelID);

              if (fileID == 0)
                {
                  gridID = vlistInqVarGrid(vlistID1, varID);
                  gridsize = gridInqSize(gridID);
                  missval = vlistInqVarMissval(vlistID1, varID);
                  gridsize = gridInqSize(vlistInqVarGrid(vlistID1, varID));  // vlistGridsizeMax(vlistID1);
                  if (weights) Free(weights);
                  weights = (double *) Malloc(gridsize * sizeof(double));
                }

              if (ef[fileID].array) Free(ef[fileID].array);
              ef[fileID].array = (double *) Malloc(gridsize * sizeof(double));

              streamID = ef[fileID].streamID;
              pstreamReadRecord(streamID, ef[fileID].array, &nmiss);
            }

          // xsize = gridInqXsize(gridID);
          // ysize = gridInqYsize(gridID);

          /*	  if ( xsize > 1 && ysize > 1 )  {
            gridWeights(gridID, weights);
            sum_weights=0;
            for ( i=0; i<gridsize; i++ )
              sum_weights += weights[i];
          }
          else*/
          {
            for (size_t i = 0; i < gridsize; i++) weights[i] = 1. / gridsize;
            sum_weights = 1.;
          }

          nmiss = 0;
          heavyside0 = 0;
          heavysideN = 0;

          for (size_t i = 0; i < gridsize; i++)
            {
              have_miss = 0;
              for (fileID = 0; fileID < nfiles; fileID++)
                {
                  val[fileID] = ef[fileID].array[i];
                  if (DBL_IS_EQUAL(val[fileID], missval))
                    {
                      have_miss = 1;
                      break;
                    }
                }

              xa = val[0];                  /* 1st file contains reference  */
              x = &val[1];                  /* Ensembles start at 2nd   file*/
              sort_iter_single(nens, x, 1); /* Sort The Ensemble Array      */
                                            /* to ascending order           */
              // only process if no missing value in ensemble

              if (!have_miss && operfunc == CRPS)
                {

                  if (xa < x[0])
                    { /* Consider outliers            */
                      beta[0] += (x[0] - xa) * weights[i];
                      heavyside0 += 1.;
                    }
                  if (xa > x[nens - 1])
                    {
                      alpha[nens] += (xa - x[nens - 1]) * weights[i];
                      alpha_weights[nens] += weights[i];
                      heavysideN += 1.;
                    }

                  /* Loop start at zero ==> 1st ensemble (c-indexing) */
                  for (k = 0; k < nens - 1; k++)
                    {                                                   /* Cumulate alpha and beta      */
                      if (xa > x[k + 1])                                /* left of heavyside            */
                        alpha[k + 1] += (x[k + 1] - x[k]) * weights[i]; /*                              */
                      else if (xa < x[k])                               /* right of heavyside           */
                        beta[k + 1] += (x[k + 1] - x[k]) * weights[i];  /*                              */
                      else if (x[k + 1] >= xa && xa >= x[k])
                        {                                              /* hitting jump pf heavyside    */
                          alpha[k + 1] += (xa - x[k]) * weights[i];    /* (occurs exactly once!) */
                          beta[k + 1] += (x[k + 1] - xa) * weights[i]; /* ****************************
                                                                        */
                        }
                    }
                }
              else if (operfunc == BRS)
                {
                  //  int occ = xa > brs_thresh? 1 : 0;

                  // brs_g[i] - number of enemble members with rank i that
                  // forecast event
                  //          - event: value > brs_thresh
                  //
                  if (x[0] > brs_thresh)
                    brs_g[0] += weights[i];
                  else if (x[nens - 1] < brs_thresh)
                    brs_g[nens] += weights[i];
                  else
                    for (k = 0; k < nens - 1; k++)
                      {
                        if (x[k + 1] >= brs_thresh && brs_thresh >= x[k])
                          {
                            brs_g[k + 1] += weights[i];
                            break;
                          }
                      }

                  // brs_o[i] - number of times that the obs is between
                  //            Ensemble i-1 and i
                  if (1)
                    {
                      if (x[0] > xa)
                        brs_o[0] += weights[i];
                      else if (x[nens - 1] < xa)
                        brs_o[nens] += weights[i];
                      else
                        for (k = 0; k < nens - 1; k++)
                          {
                            if (x[k + 1] >= xa && xa >= x[k])
                              {
                                brs_o[k + 1] += weights[i];
                                break;
                              }
                          }
                    }
                }
            }  // for ( i=0; i<gridsize; i++ )

          if (operfunc == CRPS)
            {

              // First Bin
              p = 0.;
              g = 0.;
              o = heavyside0 / gridsize;
              if (o > 0.)
                {
                  g = beta[0] / o;
                }
              crps_reli = g * (o - p) * (o - p);
              crps_pot = g * o * (1. - o);
              crps = g * ((1. - o) * p * p + o * (1. - p) * (1. - p));

              // Middle Bins
              for (k = 1; k < nens; k++)
                {
                  p = (double) k / (double) nens;

                  if (!DBL_IS_EQUAL(sum_weights, 1.))
                    {
                      alpha[k] /= sum_weights;
                      beta[k] /= sum_weights;
                    }

                  g = alpha[k] + beta[k];
                  o = beta[k] / (alpha[k] + beta[k]);

                  crps_reli += g * (o - p) * (o - p);
                  crps_pot += g * o * (1. - o);
                  crps += g * ((1. - o) * p * p + o * (1. - p) * (1. - p));
                }

              // Last Bin
              p = 1.;
              g = 0.;
              o = 1. - heavysideN / gridsize;
              if (IS_NOT_EQUAL(o, 1.))
                {
                  g = alpha[nens] / (1 - o);

                  crps_reli += g * (o - p) * (o - p);
                  crps_pot += g * o * (1 - o);
                  crps += g * ((1 - o) * p * p + o * (1 - p) * (1 - p));
                }
              r[CRPS_RES] = crps;
              r[CRPS_RELI] = crps_reli;
              r[CRPS_POT] = crps_pot;
            }
          else if (operfunc == BRS)
            {
              double gsum = 0;
              double obar = 0;
              double osum = 0;
              double o, g, p;

              brs_reli = 0;
              brs_resol = 0;
              brs_uncty = 0;

              for (k = 0; k <= nens; k++)
                {
                  obar += brs_g[k] * brs_o[k];
                  gsum += brs_g[k];
                  osum += brs_o[k];
                }

              if (fabs(osum - 1) > 1.e-06 || fabs(gsum - 1) > 1.e-06)
                {
                  cdoAbort("Internal error - normalization constraint of problem not fulfilled");
                  cdoAbort("This is likely due to missing values");
                }
              o = 0;
              p = 0;
              g = 0;
              brs_uncty = obar * (1 - obar);

              for (k = 0; k <= nens; k++)
                {

                  g = brs_g[k];
                  o = brs_o[k];
                  p = 1. - k / (float) nens;
                  // need p = 1 - k/nens here as k=0 if all members forecast
                  // event and k=nens if none does so.

                  brs_reli += g * (o - p) * (o - p);
                  brs_resol += g * (o - obar) * (o - obar);
                }

              r[BRS_RES] = brs_reli - brs_resol + brs_uncty;
              r[BRS_RELI] = brs_reli;
              r[BRS_RESOL] = brs_resol;
              r[BRS_UNCTY] = brs_uncty;

              if (cdoVerbose)
                {
                  //	      cdoPrint("Brier score for var %i level %i
                  // calculated",varID, levelID);
                  cdoPrint("BRS: obar %12.6g brs  %12.6g reli %12.6g resol %12.6g u %12.6g",
                           obar, brs_reli - brs_resol + brs_uncty, brs_reli, brs_resol, brs_uncty);
                }
            }

          if (cdoVerbose && operfunc == CRPS)
            cdoPrint("CRPS:%12.6g reli:%12.6g crps_pot:%12.6g crps:%12.6g", crps, crps_reli, crps_pot, crps_reli + crps_pot);

          for (stream = 0; stream < nostreams; stream++)
            {
              pstreamDefRecord(streamID2[stream], varID, levelID);
              if (DBL_IS_NAN(r[stream]))
                {
                  r[stream] = missval;
                  have_miss = 1;
                }
              pstreamWriteRecord(streamID2[stream], &r[stream], have_miss);
            }

          switch (operfunc)
            {
            case (CRPS):
              arrayFill(nens + 1, alpha, 0.0);
              arrayFill(nens + 1, beta, 0.0);
              heavyside0 = 0;
              heavysideN = 0;
              break;
            case (BRS):
              arrayFill(nens + 1, brs_o, 0.0);
              arrayFill(nens + 1, brs_g, 0.0);
              break;
            }
        }  // for ( int recID = 0; recID < nrecs; recID++ )
      tsID++;
    }
  while (nrecs);

  for (fileID = 0; fileID < nfiles; fileID++)
    {
      streamID = ef[fileID].streamID;
      pstreamClose(streamID);
    }

  for (stream = 0; stream < nostreams; stream++) pstreamClose(streamID2[stream]);

  for (fileID = 0; fileID < nfiles; fileID++)
    if (ef[fileID].array) Free(ef[fileID].array);

  for (stream = 0; stream < nostreams; stream++)
    {
      vlistDestroy(vlistID2[stream]);
      taxisDestroy(taxisID2[stream]);
      zaxisDestroy(zaxisID2[stream]);
    }

  //  vlistDestroy(vlistID1);
  //  taxisDestroy(taxisID1);
  //  zaxisDestroy(zaxisID1);

  gridDestroy(gridID);
  gridDestroy(gridID2);

  if (ef) Free(ef);
  if (weights) Free(weights);
  if (r) Free(r);
  if (alpha) Free(alpha);
  if (beta) Free(beta);
  if (alpha_weights) Free(alpha_weights);
  if (beta_weights) Free(beta_weights);
  if (brs_g) Free(brs_g);
  if (brs_o) Free(brs_o);
  if (brs_g_weights) Free(brs_g_weights);
  if (brs_o_weights) Free(brs_o_weights);
  if (val) Free(val);
  if (vlistID2) Free(vlistID2);
  if (streamID2) Free(streamID2);
  if (zaxisID2) Free(zaxisID2);
  if (taxisID2) Free(taxisID2);
  cdoFinish();

  return 0;
}
