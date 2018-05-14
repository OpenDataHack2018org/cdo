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
#include <stdio.h>
#include <stdlib.h> /* exit */
#include <string.h>
#include <math.h>
#include <assert.h>

#include "compare.h"
#include "array.h"
#include "constants.h"
#include "after_vertint.h"

#define SCALEHEIGHT (-7000.)
#define SCALESLP (101325.0)

int Mars = 0;

void
height2pressure(double *restrict phlev, const double *restrict hlev, long nphlev)
{
  double exp_arg;
  double height;

  for (long k = 0; k < nphlev; k++)
    {
      height = hlev[k];
      /*
        unitsel == 1 : hlev[k] is given in meters
        unitsel == 2 : hlev[k] is given in kilometers
        height2pressure needs meters (MKSC-standard)
      */

      exp_arg = height / SCALEHEIGHT;

      phlev[k] = SCALESLP * exp(exp_arg);
    }
}

void
pressure2height(double *restrict hlev, const double *restrict plev, long nphlev)
{
  for (long k = 0; k < nphlev; k++)
    {
      hlev[k] = log(plev[k] / SCALESLP) * SCALEHEIGHT;
    }
}

void
presh(double *restrict fullp, double *halfp, const double *restrict vct, const double *restrict ps, long nhlev, long ngp)
{
  if (ps == NULL)
    {
      fprintf(stderr, "ps undefined!\n");
      exit(EXIT_FAILURE);
    }

  double *halfpres = halfp;
  for (long lh = 0; lh < nhlev; lh++)
    {
      double zp = vct[lh];
      double ze = vct[lh + nhlev + 1];
      for (long i = 0; i < ngp; i++) halfpres[i] = zp + ze * ps[i];
      halfpres += ngp;
    }
  arrayCopy(ngp, ps, halfpres);

  if (fullp)
    {
      halfpres = halfp;
      for (long i = 0; i < ngp * nhlev; i++) fullp[i] = 0.5 * (halfpres[i] + halfpres[i + ngp]);
    }
}

void
genind(int *nx, const double *restrict plev, const double *restrict fullp, long ngp, long nplev, long nhlev)
{
  arrayFill(ngp * nplev, nx, 0);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(nx, plev, fullp, ngp, nplev, nhlev)
#endif
  for (long lp = 0; lp < nplev; lp++)
    {
      const double pres = plev[lp];
      int *restrict nxl = nx + lp * ngp;
      for (long lh = 0; lh < nhlev; lh++)
        {
          const double *restrict fullpx = fullp + lh * ngp;
          for (long i = 0; i < ngp; i++)
            {
              if (pres > fullpx[i]) nxl[i] = lh;
            }
        }
    }
}

void
genindmiss(int *nx, const double *restrict plev, int ngp, int nplev, const double *restrict ps_prog, size_t *restrict pnmiss)
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(nx, plev, ngp, nplev, ps_prog, pnmiss)
#endif
  for (long lp = 0; lp < nplev; lp++)
    {
      pnmiss[lp] = 0;
      const double pres = plev[lp];
      int *restrict nxl = nx + lp * ngp;
      for (long i = 0; i < ngp; i++)
        {
          if (pres > ps_prog[i])
            {
              nxl[i] = -1;
              pnmiss[lp]++;
            }
        }
    }

} /* genindmiss */

void
extra_P(double *restrict slp, const double *restrict halfp, const double *restrict fullp, const double *restrict geop,
        const double *restrict temp, long ngp)
{
  double alpha, tstar, tmsl, zprt, zprtal;
  const double zlapse = 0.0065;
  const double zrg = 1.0 / PlanetGrav;

  for (long j = 0; j < ngp; ++j)
    {
      if (geop[j] < 0.0001 && geop[j] > -0.0001)
        slp[j] = halfp[j];
      else
        {
          alpha = PlanetRD * zlapse * zrg;
          tstar = (1.0 + alpha * (halfp[j] / fullp[j] - 1.0)) * temp[j];

          if (tstar < 255.0) tstar = 0.5 * (255.0 + tstar);

          tmsl = tstar + zlapse * zrg * geop[j];
          if (tmsl > 290.5 && tstar > 290.5)
            {
              tstar = 0.5 * (290.5 + tstar);
              tmsl = tstar;
            }

          if (tmsl - tstar < 0.000001 && tstar - tmsl < 0.000001)
            alpha = 0.0;
          else if (geop[j] > 0.0001 || geop[j] < -0.0001)
            alpha = PlanetRD * (tmsl - tstar) / geop[j];

          zprt = geop[j] / (PlanetRD * tstar);
          zprtal = zprt * alpha;
          slp[j] = halfp[j] * exp(zprt * (1.0 - zprtal * (0.5 - zprtal / 3.0)));
        }
    }
} /* extrap */

static double
extra_T(double pres, double halfp, double fullp, double geop, double temp)
{
  double peval;
  const double zlapse = 0.0065;
  const double zrg = 1.0 / PlanetGrav;
  double tstar = (1.0 + zlapse * PlanetRD * zrg * (halfp / fullp - 1.0)) * temp;
  const double ztsz = tstar;
  const double z1 = tstar + zlapse * zrg * geop;

  if (tstar < 255.0) tstar = 0.5 * (255.0 + tstar);

  double ztmsl = tstar + zlapse * zrg * geop;

  if (ztmsl > 290.5 && tstar > 290.5)
    {
      tstar = 0.5 * (290.5 + tstar);
      ztmsl = tstar;
    }

  if (ztmsl > 290.5 && tstar <= 290.5) ztmsl = 290.5;

  if (pres <= halfp)
    peval = ((halfp - pres) * temp + (pres - fullp) * tstar) / (halfp - fullp);
  else
    {
      double ztmsl = z1;
      tstar = ztsz;
      const double zhts = geop * zrg;

      if (zhts > 2000. && z1 > 298.)
        {
          ztmsl = 298.;
          if (zhts < 2500.) ztmsl = 0.002 * ((2500. - zhts) * z1 + (zhts - 2000.) * ztmsl);
        }

      double zalph;
      if ((ztmsl - tstar) < 0.000001)
        zalph = 0.;
      else if (geop > 0.0001 || geop < -0.0001)
        zalph = PlanetRD * (ztmsl - tstar) / geop;
      else
        zalph = PlanetRD * zlapse * zrg;

      double zalp = zalph * log(pres / halfp);
      peval = tstar * (1.0 + zalp * (1.0 + zalp * (0.5 + 0.16666666667 * zalp)));
    }

  return peval;
} /* extra_T */

static double
extra_Z(double pres, double halfp, double fullp, double geop, double temp)
{
  const double zlapse = 0.0065;
  const double ztlim = 290.5;
  const double zrg = 1.0 / PlanetGrav;
  double alpha = PlanetRD * zlapse * zrg;
  double tstar = (1.0 + alpha * (halfp / fullp - 1.0)) * temp;

  if (tstar < 255.0) tstar = 0.5 * (255.0 + tstar);

  double tmsl = tstar + zlapse * zrg * geop;

  if (tmsl > ztlim && tstar > ztlim)
    {
      tstar = 0.5 * (ztlim + tstar);
      tmsl = tstar;
    }

  if (tmsl > ztlim && tstar <= ztlim) tmsl = ztlim;

  if (tmsl - tstar < 0.000001 && tstar - tmsl < 0.000001)
    alpha = 0.0;
  else if (geop > 0.0001 || geop < -0.0001)
    alpha = PlanetRD * (tmsl - tstar) / geop;

  const double zalp = log(pres / halfp);
  const double zalpal = zalp * alpha;

  return (geop - PlanetRD * tstar * zalp * (1.0 + zalpal * (0.5 + zalpal / 6.0))) * zrg;
} /* extra_Z */

void
interp_X(const double *restrict gt, double *pt, const double *restrict hyb_press, const int *nx, const double *restrict plev,
         long nplev, long ngp, long nhlev, double missval)
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(gt, pt, hyb_press, nx, plev, nplev, ngp, nhlev, missval)
#endif
  for (long lp = 0; lp < nplev; lp++)
    {
      long nl, nh;
      const double pres = plev[lp];
      const int *restrict nxl = nx + lp * ngp;
      double *restrict ptl = pt + lp * ngp;
      for (long i = 0; i < ngp; i++)
        {
          if (nxl[i] == -1)
            ptl[i] = missval;
          else
            {
              nl = nxl[i] * ngp + i;
              nh = nl + ngp;
              ptl[i] = (nh >= ngp * nhlev) ? gt[nl]
                                           : gt[nl] + (pres - hyb_press[nl]) * (gt[nh] - gt[nl]) / (hyb_press[nh] - hyb_press[nl]);
            }
        }
    }
} /* interp_X */

void
interp_T(const double *restrict geop, const double *restrict gt, double *pt, const double *restrict fullp,
         const double *restrict halfp, const int *nx, const double *restrict plev, long nplev, long ngp, long nhlev, double missval)
{
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(geop, gt, pt, fullp, halfp, nx, plev, nplev, ngp, nhlev, missval, Mars)
#endif
  for (long lp = 0; lp < nplev; lp++)
    {
      long nl, nh;
      const double pres = plev[lp];
      const int *restrict nxl = nx + lp * ngp;
      double *restrict ptl = pt + lp * ngp;
#if defined(CRAY)
#pragma _CRI inline extra_T
#endif
      for (long i = 0; i < ngp; i++)
        {
          nl = nxl[i];
          if (nl < 0)
            ptl[i] = missval;
          else
            {
              if (nl > nhlev - 2)
                {
                  if (Mars)
                    ptl[i] = gt[(nhlev - 1) * ngp + i];
                  else
#if defined(SX)
#pragma cdir inline
#endif
                    ptl[i]
                        = extra_T(pres, halfp[nhlev * ngp + i], fullp[(nhlev - 1) * ngp + i], geop[i], gt[(nhlev - 1) * ngp + i]);
                }
              else
                {
                  nh = nl + 1;
                  ptl[i] = gt[nl * ngp + i]
                           + (pres - fullp[nl * ngp + i]) * (gt[nh * ngp + i] - gt[nl * ngp + i])
                                 / (fullp[nh * ngp + i] - fullp[nl * ngp + i]);
                }
            }
        }
    }
} /* interp_T */

void
interp_Z(const double *restrict geop, const double *restrict gz, double *pz, const double *restrict fullp,
         const double *restrict halfp, const int *nx, const double *restrict gt, const double *restrict plev, long nplev, long ngp,
         long nhlev, double missval)
{
  assert(geop != NULL);
  assert(gz != NULL);
  assert(pz != NULL);
  assert(fullp != NULL);
  assert(halfp != NULL);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(geop, gz, pz, fullp, halfp, nx, gt, plev, nplev, ngp, nhlev, missval, Mars)
#endif
  for (long lp = 0; lp < nplev; lp++)
    {
      long nl, nh;
      const double pres = plev[lp];
      const int *restrict nxl = nx + lp * ngp;
      double *restrict pzl = pz + lp * ngp;
#if defined(CRAY)
#pragma _CRI inline extra_Z
#endif
      for (long i = 0; i < ngp; i++)
        {
          nl = nxl[i];
          if (nl < 0)
            pzl[i] = missval;
          else
            {
              if (pres > halfp[(nl + 1) * ngp + i]) nl++;

              if (nl > nhlev - 1)
                {
                  if (Mars)
                    pzl[i] = gt[(nhlev - 1) * ngp + i];
                  else
#if defined(SX)
#pragma cdir inline
#endif
                    pzl[i]
                        = extra_Z(pres, halfp[nhlev * ngp + i], fullp[(nhlev - 1) * ngp + i], geop[i], gt[(nhlev - 1) * ngp + i]);
                }
              else
                {
                  nh = nl + 1;
                  pzl[i] = gz[nl * ngp + i]
                           + (pres - halfp[nl * ngp + i]) * (gz[nh * ngp + i] - gz[nl * ngp + i])
                                 / (halfp[nh * ngp + i] - halfp[nl * ngp + i]);
                }
            }
        }
    }
} /* interp_Z */
