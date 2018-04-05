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

*/

#include <cdi.h>

#include "cdo_int.h"
#include "pstream_int.h"
#include "grid.h"

static int
gentpngrid(int gridID1)
{
  size_t ilat, ilon, ilonr, k, kr;

  size_t nlon1 = gridInqXsize(gridID1);
  size_t nlat1 = gridInqYsize(gridID1);

  size_t nlon2 = nlon1;
  size_t nlat2 = nlat1 + 2;

  int gridtype = gridInqType(gridID1);
  int prec = gridInqDatatype(gridID1);

  int gridID2 = gridCreate(gridtype, nlon2 * nlat2);
  gridDefXsize(gridID2, nlon2);
  gridDefYsize(gridID2, nlat2);

  gridDefDatatype(gridID2, prec);

  grid_copy_attributes(gridID1, gridID2);

  char xunits[CDI_MAX_NAME];
  xunits[0] = 0;
  char yunits[CDI_MAX_NAME];
  yunits[0] = 0;
  cdiGridInqKeyStr(gridID1, CDI_KEY_XUNITS, CDI_MAX_NAME, xunits);
  cdiGridInqKeyStr(gridID1, CDI_KEY_YUNITS, CDI_MAX_NAME, yunits);

  if (gridInqXvals(gridID1, NULL) && gridInqYvals(gridID1, NULL))
    {
      if (gridtype == GRID_CURVILINEAR)
        {
          std::vector<double> xvals1(nlon1 * nlat1);
          std::vector<double> yvals1(nlon1 * nlat1);
          std::vector<double> xvals2(nlon2 * nlat2);
          std::vector<double> yvals2(nlon2 * nlat2);

          gridInqXvals(gridID1, xvals1.data());
          gridInqYvals(gridID1, yvals1.data());

          for (ilat = 0; ilat < nlat1; ilat++)
            {
              for (ilon = 0; ilon < nlon1; ilon++)
                {
                  xvals2[(ilat + 2) * nlon1 + ilon] = xvals1[ilat * nlon1 + ilon];
                  yvals2[(ilat + 2) * nlon1 + ilon] = yvals1[ilat * nlon1 + ilon];
                }
            }

          for (ilon = 0; ilon < nlon1; ilon++)
            {
              ilonr = nlon1 - ilon - 1;
              xvals2[1 * nlon1 + ilon] = xvals2[2 * nlon1 + ilonr]; /* syncronise line 2 with line 3 */
              xvals2[0 * nlon1 + ilon] = xvals2[3 * nlon1 + ilonr]; /* syncronise line 1 with line 4 */
              yvals2[1 * nlon1 + ilon] = yvals2[2 * nlon1 + ilonr]; /* syncronise line 2 with line 3 */
              yvals2[0 * nlon1 + ilon] = yvals2[3 * nlon1 + ilonr]; /* syncronise line 1 with line 4 */
            }

          gridDefXvals(gridID2, xvals2.data());
          gridDefYvals(gridID2, yvals2.data());
        }
    }

  if (gridInqXbounds(gridID1, NULL) && gridInqYbounds(gridID1, NULL))
    {
      if (gridtype == GRID_CURVILINEAR)
        {
          std::vector<double> xbounds1(4 * nlon1 * nlat1);
          std::vector<double> ybounds1(4 * nlon1 * nlat1);
          std::vector<double> xbounds2(4 * nlon2 * nlat2);
          std::vector<double> ybounds2(4 * nlon2 * nlat2);

          gridInqXbounds(gridID1, xbounds1.data());
          gridInqYbounds(gridID1, ybounds1.data());

          if (gridtype == GRID_CURVILINEAR)
            {
              gridDefNvertex(gridID2, 4);

              for (ilat = 0; ilat < nlat1; ilat++)
                {
                  for (ilon = 0; ilon < 4 * nlon1; ilon++)
                    {
                      xbounds2[4 * (ilat + 2) * nlon1 + ilon] = xbounds1[4 * ilat * nlon1 + ilon];
                      ybounds2[4 * (ilat + 2) * nlon1 + ilon] = ybounds1[4 * ilat * nlon1 + ilon];
                    }
                }

              for (ilon = 0; ilon < nlon1; ilon++)
                {
                  ilonr = nlon1 - ilon - 1;
                  for (k = 0; k < 4; ++k)
                    {
                      kr = 3 - k;
                      xbounds2[4 * 1 * nlon1 + 4 * ilon + k] = xbounds2[4 * 2 * nlon1 + 4 * ilonr + kr];
                      xbounds2[4 * 0 * nlon1 + 4 * ilon + k] = xbounds2[4 * 3 * nlon1 + 4 * ilonr + kr];
                      ybounds2[4 * 1 * nlon1 + 4 * ilon + k] = ybounds2[4 * 2 * nlon1 + 4 * ilonr + kr];
                      ybounds2[4 * 0 * nlon1 + 4 * ilon + k] = ybounds2[4 * 3 * nlon1 + 4 * ilonr + kr];
                    }
                }
              /*
              for ( ilon = 0; ilon < 4*nlon1; ilon++ )
                {
                  ilonr = 4*nlon1 - ilon - 1;
                  xbounds2[4*1*nlon1 + ilon] = xbounds2[4*2*nlon1 + ilonr]; xbounds2[4*0*nlon1 + ilon] = xbounds2[4*3*nlon1 + ilonr];
                  ybounds2[4*1*nlon1 + ilon] = ybounds2[4*2*nlon1 + ilonr]; ybounds2[4*0*nlon1 + ilon] = ybounds2[4*3*nlon1 + ilonr];
                }
              */
            }

          gridDefXbounds(gridID2, xbounds2.data());
          gridDefYbounds(gridID2, ybounds2.data());
        }
    }

  return gridID2;
}

static int
gengrid(int gridID1, int lhalo, int rhalo)
{
  long i, ilat, ilon;
  double cpi2 = M_PI * 2;

  long nlon1 = gridInqXsize(gridID1);
  long nlat1 = gridInqYsize(gridID1);

  long nlon2 = nlon1 + lhalo + rhalo;
  long nlat2 = nlat1;

  long nmin = 0;
  long nmax = nlon1;
  if (lhalo < 0) nmin = -lhalo;
  if (rhalo < 0) nmax += rhalo;
  /*
  printf("nlon1=%d, nlon2=%d, lhalo=%d, rhalo=%d nmin=%d, nmax=%d\n",
         nlon1, nlon2, lhalo, rhalo, nmin, nmax);
  */
  int gridtype = gridInqType(gridID1);
  int prec = gridInqDatatype(gridID1);

  int gridID2 = gridCreate(gridtype, nlon2 * nlat2);
  gridDefXsize(gridID2, nlon2);
  gridDefYsize(gridID2, nlat2);

  gridDefDatatype(gridID2, prec);

  grid_copy_attributes(gridID1, gridID2);

  char xunits[CDI_MAX_NAME];
  xunits[0] = 0;
  char yunits[CDI_MAX_NAME];
  yunits[0] = 0;
  cdiGridInqKeyStr(gridID1, CDI_KEY_XUNITS, CDI_MAX_NAME, xunits);
  cdiGridInqKeyStr(gridID1, CDI_KEY_YUNITS, CDI_MAX_NAME, yunits);

  if (memcmp(xunits, "degree", 6) == 0) cpi2 *= RAD2DEG;

  if (gridInqXvals(gridID1, NULL) && gridInqYvals(gridID1, NULL))
    {
      std::vector<double> xvals1, yvals1, xvals2, yvals2;
      if (gridtype == GRID_CURVILINEAR)
        {
          xvals1.resize(nlon1 * nlat1);
          yvals1.resize(nlon1 * nlat1);
          xvals2.resize(nlon2 * nlat2);
          yvals2.resize(nlon2 * nlat2);
        }
      else
        {
          xvals1.resize(nlon1);
          yvals1.resize(nlat1);
          xvals2.resize(nlon2);
          yvals2.resize(nlat2);
        }

      double *pxvals2 = xvals2.data();
      double *pyvals2 = yvals2.data();

      gridInqXvals(gridID1, xvals1.data());
      gridInqYvals(gridID1, yvals1.data());

      if (gridtype == GRID_CURVILINEAR)
        {
          for (ilat = 0; ilat < nlat2; ilat++)
            {
              for (ilon = nlon1 - lhalo; ilon < nlon1; ilon++)
                {
                  *pxvals2++ = xvals1[ilat * nlon1 + ilon];
                  *pyvals2++ = yvals1[ilat * nlon1 + ilon];
                }

              for (ilon = nmin; ilon < nmax; ilon++)
                {
                  *pxvals2++ = xvals1[ilat * nlon1 + ilon];
                  *pyvals2++ = yvals1[ilat * nlon1 + ilon];
                }

              for (ilon = 0; ilon < rhalo; ilon++)
                {
                  *pxvals2++ = xvals1[ilat * nlon1 + ilon];
                  *pyvals2++ = yvals1[ilat * nlon1 + ilon];
                }
            }
        }
      else
        {
          for (i = nlon1 - lhalo; i < nlon1; i++) *pxvals2++ = xvals1[i] - cpi2;
          for (i = nmin; i < nmax; i++) *pxvals2++ = xvals1[i];
          for (i = 0; i < rhalo; i++) *pxvals2++ = xvals1[i] + cpi2;

          for (i = 0; i < nlat1; i++) yvals2[i] = yvals1[i];
        }
      /*
      for ( i = 0; i < nlat2; i++ ) printf("lat : %d %g\n", i+1, yvals2[i]);
      for ( i = 0; i < nlon2; i++ ) printf("lon : %d %g\n", i+1, xvals2[i]);
      */
      gridDefXvals(gridID2, xvals2.data());
      gridDefYvals(gridID2, yvals2.data());
    }

  if (gridInqXbounds(gridID1, NULL) && gridInqYbounds(gridID1, NULL))
    {
      std::vector<double> xbounds1, ybounds1, xbounds2, ybounds2;
      if (gridtype == GRID_CURVILINEAR)
        {
          xbounds1.resize(4 * nlon1 * nlat1);
          ybounds1.resize(4 * nlon1 * nlat1);
          xbounds2.resize(4 * nlon2 * nlat2);
          ybounds2.resize(4 * nlon2 * nlat2);
        }
      else
        {
          xbounds1.resize(2 * nlon1);
          ybounds1.resize(2 * nlat1);
          xbounds2.resize(2 * nlon2);
          ybounds2.resize(2 * nlat2);
        }

      double *pxbounds2 = xbounds2.data();
      double *pybounds2 = ybounds2.data();

      gridInqXbounds(gridID1, xbounds1.data());
      gridInqYbounds(gridID1, ybounds1.data());

      if (gridtype == GRID_CURVILINEAR)
        {
          gridDefNvertex(gridID2, 4);
          for (ilat = 0; ilat < nlat1; ilat++)
            {
              for (ilon = 4 * (nlon1 - lhalo); ilon < 4 * nlon1; ilon++)
                {
                  *pxbounds2++ = xbounds1[4 * ilat * nlon1 + ilon];
                  *pybounds2++ = ybounds1[4 * ilat * nlon1 + ilon];
                }

              for (ilon = 4 * nmin; ilon < 4 * nmax; ilon++)
                {
                  *pxbounds2++ = xbounds1[4 * ilat * nlon1 + ilon];
                  *pybounds2++ = ybounds1[4 * ilat * nlon1 + ilon];
                }

              for (ilon = 0; ilon < 4 * rhalo; ilon++)
                {
                  *pxbounds2++ = xbounds1[4 * ilat * nlon1 + ilon];
                  *pybounds2++ = ybounds1[4 * ilat * nlon1 + ilon];
                }
            }
        }
      else
        {
          gridDefNvertex(gridID2, 2);
          for (i = 2 * (nlon1 - lhalo); i < 2 * nlon1; i++) *pxbounds2++ = xbounds1[i] - cpi2;
          for (i = 2 * nmin; i < 2 * nmax; i++) *pxbounds2++ = xbounds1[i];
          for (i = 0; i < 2 * rhalo; i++) *pxbounds2++ = xbounds1[i] + cpi2;

          for (i = 0; i < 2 * nlat2; i++) ybounds2[i] = ybounds1[i];
        }

      gridDefXbounds(gridID2, xbounds2.data());
      gridDefYbounds(gridID2, ybounds2.data());
    }

  return gridID2;
}

static int
genindexgrid(int gridID1, int *lhalo, int *rhalo)
{
  operatorCheckArgc(2);

  *lhalo = parameter2int(operatorArgv()[0]);
  *rhalo = parameter2int(operatorArgv()[1]);

  int nlon1 = gridInqXsize(gridID1);

  if (*lhalo > nlon1)
    {
      *lhalo = nlon1;
      cdoWarning("left halo out of range. Set to %d.", *lhalo);
    }

  if (*lhalo < 0 && -(*lhalo) > nlon1 / 2)
    {
      *lhalo = -nlon1 / 2;
      cdoWarning("left halo out of range. Set to %d.", *lhalo);
    }

  if (*rhalo > nlon1)
    {
      *rhalo = nlon1;
      cdoWarning("right halo out of range. Set to %d.", *rhalo);
    }

  if (*rhalo < 0 && -(*rhalo) > nlon1 / 2)
    {
      *rhalo = -nlon1 / 2;
      cdoWarning("right halo out of range. Set to %d.", *rhalo);
    }

  int gridID2 = gengrid(gridID1, *lhalo, *rhalo);

  return gridID2;
}

static void
halo(double *array1, int gridID1, double *array2, int lhalo, int rhalo)
{
  int nlon1 = gridInqXsize(gridID1);
  int nlat = gridInqYsize(gridID1);

  int nmin = 0;
  int nmax = nlon1;
  if (lhalo < 0) nmin = -lhalo;
  if (rhalo < 0) nmax += rhalo;

  for (int ilat = 0; ilat < nlat; ilat++)
    {
      for (int ilon = nlon1 - lhalo; ilon < nlon1; ilon++) *array2++ = array1[ilat * nlon1 + ilon];

      for (int ilon = nmin; ilon < nmax; ilon++) *array2++ = array1[ilat * nlon1 + ilon];

      for (int ilon = 0; ilon < rhalo; ilon++) *array2++ = array1[ilat * nlon1 + ilon];
    }
}

static void
tpnhalo(double *array1, int gridID1, double *array2)
{
  size_t nlon = gridInqXsize(gridID1);
  size_t nlat = gridInqYsize(gridID1);

  for (size_t ilat = 0; ilat < nlat; ilat++)
    for (size_t ilon = 0; ilon < nlon; ilon++) array2[(ilat + 2) * nlon + ilon] = array1[ilat * nlon + ilon];

  for (size_t ilon = 0; ilon < nlon; ilon++)
    {
      size_t ilonr = nlon - ilon - 1;
      array2[1 * nlon + ilon] = array2[2 * nlon + ilonr]; /* syncronise line 2 with line 3 */
      array2[0 * nlon + ilon] = array2[3 * nlon + ilonr]; /* syncronise line 1 with line 4 */
    }
}

void *
Sethalo(void *process)
{
  int nrecs;
  int varID, levelID;
  int gridID1 = -1, gridID2;
  int index;
  size_t nmiss;
  int lhalo = 0, rhalo = 0;

  cdoInitialize(process);

  // clang-format off
  int SETHALO = cdoOperatorAdd("sethalo", 0, 0, NULL);
                cdoOperatorAdd("tpnhalo", 0, 0, NULL);
  // clang-format on

  int operatorID = cdoOperatorID();

  int streamID1 = cdoStreamOpenRead(cdoStreamName(0));
  int vlistID1 = cdoStreamInqVlist(streamID1);

  int ngrids = vlistNgrids(vlistID1);
  int ndiffgrids = 0;
  for (index = 1; index < ngrids; index++)
    if (vlistGrid(vlistID1, 0) != vlistGrid(vlistID1, index)) ndiffgrids++;

  for (index = 0; index < ngrids; index++)
    {
      gridID1 = vlistGrid(vlistID1, index);
      int gridtype = gridInqType(gridID1);
      if (gridtype == GRID_LONLAT || gridtype == GRID_GAUSSIAN) break;
      if (gridtype == GRID_CURVILINEAR) break;
      if (gridtype == GRID_GENERIC && gridInqXsize(gridID1) > 0 && gridInqYsize(gridID1) > 0) break;
    }

  if (gridInqType(gridID1) == GRID_GAUSSIAN_REDUCED)
    cdoAbort("Gaussian reduced grid found. Use option -R to convert it to a regular grid!");

  if (index == ngrids) cdoAbort("No regular grid found!");
  if (ndiffgrids > 0) cdoAbort("Too many different grids!");

  if (operatorID == SETHALO)
    {
      operatorInputArg("left and right halo");
      gridID2 = genindexgrid(gridID1, &lhalo, &rhalo);
    }
  else
    {
      gridID2 = gentpngrid(gridID1);
    }

  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  ngrids = vlistNgrids(vlistID1);
  for (index = 0; index < ngrids; index++)
    {
      if (gridID1 == vlistGrid(vlistID1, index))
        {
          vlistChangeGridIndex(vlistID2, index, gridID2);
          break;
        }
    }

  int nvars = vlistNvars(vlistID1);
  std::vector<bool> vars(nvars);
  for (varID = 0; varID < nvars; varID++) vars[varID] = (gridID1 == vlistInqVarGrid(vlistID1, varID));

  int streamID2 = cdoStreamOpenWrite(cdoStreamName(1), cdoFiletype());
  pstreamDefVlist(streamID2, vlistID2);

  size_t gridsize = gridInqSize(gridID1);
  std::vector<double> array1(gridsize);

  size_t gridsize2 = gridInqSize(gridID2);
  std::vector<double> array2(gridsize2);

  int tsID = 0;
  while ((nrecs = cdoStreamInqTimestep(streamID1, tsID)))
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      pstreamDefTimestep(streamID2, tsID);

      for (int recID = 0; recID < nrecs; recID++)
        {
          pstreamInqRecord(streamID1, &varID, &levelID);

          if (vars[varID])
            {
              pstreamReadRecord(streamID1, array1.data(), &nmiss);

              if (operatorID == SETHALO)
                halo(array1.data(), gridID1, array2.data(), lhalo, rhalo);
              else
                tpnhalo(array1.data(), gridID1, array2.data());

              if (nmiss)
                {
                  double missval = vlistInqVarMissval(vlistID1, varID);
                  nmiss = arrayNumMV(gridsize2, array2.data(), missval);
                }

              pstreamDefRecord(streamID2, varID, levelID);
              pstreamWriteRecord(streamID2, array2.data(), nmiss);
            }
        }

      tsID++;
    }

  pstreamClose(streamID2);
  pstreamClose(streamID1);

  cdoFinish();

  return 0;
}
