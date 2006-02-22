#include <stdio.h>
#include <math.h>
#include <string.h>

#include "cdo.h"
#include "cdi.h"
#ifndef _DMEMORY_H
#  include "dmemory.h"
#endif

#ifndef _ERROR_H
#  include "error.h"
#endif
#include "specspace.h"


#define  C_EARTH_RADIUS  (6371000.0)
double PlanetRadius = C_EARTH_RADIUS;

void geninx(int nt, double *f, double *g)
{
  int m2,n2;
  int m, n ;

  for ( m = 0; m <= nt; m++ )
    {
      m2 = m * m;
      for ( n = m; n <= nt; n++ )
	{
	  n2 = n * n;
	  if ( n )
	    {
	      *g++ = -PlanetRadius / n * sqrt((double)(n2-m2)/(double)(4*n2-1));
	      *f++ = -PlanetRadius * m / (double)(n2+n);
	    }
	  else
	    {
	      *g++ = 0.0;
	      *f++ = 0.0;
	    }
	}
    }
}


void legini(int trunc, int nlat, double *poli, double *pold, double *rcoslat)
{
  static char func[] = "legini";
  int waves, dimsp, dimpnm;
  int jgl, jm, jn, is;
  int isp, latn, lats, pdim;
  double *gmu, *gwt, *pnm, *work;

  waves  = trunc + 1;
  dimsp  = (trunc + 1)*(trunc + 2);
  dimpnm = (trunc + 1)*(trunc + 4)/2;
  pdim   = dimsp / 2 * nlat;

  gmu  = (double *) malloc(nlat * sizeof(double));
  gwt  = (double *) malloc(nlat * sizeof(double));
  pnm  = (double *) malloc(dimpnm * sizeof(double));
  work = (double *) malloc(3*waves * sizeof(double));

  gaussaw(gmu, gwt, nlat);
  for ( jgl = 0; jgl < nlat; jgl++ ) gwt[jgl] *= 0.5;

  for ( jgl = 0; jgl < nlat; jgl++ )
    rcoslat[jgl] = 1.0 / sqrt(1.0 - gmu[jgl]*gmu[jgl]);

  for ( jgl = 0; jgl < nlat/2; jgl++ )
    {
      jspleg1(pnm, gmu[jgl], trunc, work);

      latn = jgl;
      isp = 0;
      for ( jm = 0; jm < waves; jm++ )
	{
#if defined (SX)
#pragma vdir nodep
#endif
	  for ( jn = 0; jn < waves - jm; jn++ )
	    {
	      is = (jn+1)%2 * 2 - 1;
	      lats = latn - jgl + nlat - jgl - 1;
	      poli[latn] = pnm[isp];
	      pold[latn] = pnm[isp] * gwt[jgl];
	      poli[lats] = pnm[isp] * is;
	      pold[lats] = pnm[isp] * gwt[jgl] * is;
	      latn += nlat;
	      isp++;
	    }
	  isp += 1;
	}
    }

  free(work);
  free(pnm);
  free(gwt);
  free(gmu);
}


void grid2spec(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut)
{
  static char func[] = "grid2spec";
  int trunc, nlat, nlon, nfc;
  int nlev = 1;
  int waves;
  double *fpwork;
    
  trunc = gridInqTrunc(gridIDout);
  nlon  = gridInqXsize(gridIDin);
  nlat  = gridInqYsize(gridIDin);

  waves = trunc + 1;
  nfc   = waves * 2;

  fpwork = (double *) malloc(nlat*nfc*nlev*sizeof(double));

  gp2fc(sptrans->trig, sptrans->ifax, arrayIn, fpwork, nlat, nlon, nlev, nfc);
  fc2sp(fpwork, arrayOut, sptrans->pold, nlev, nlat, nfc, trunc);

  free(fpwork);
}
	   
   
void spec2grid(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut)
{
  static char func[] = "spec2grid";
  int trunc, nlat, nlon, nfc;
  int nlev = 1;
  int waves;
  double *fpwork;
    
  trunc = gridInqTrunc(gridIDin);
  nlon  = gridInqXsize(gridIDout);
  nlat  = gridInqYsize(gridIDout);

  waves = trunc + 1;
  nfc   = waves * 2;

  fpwork = (double *) malloc(nlat*nfc*nlev*sizeof(double));

  sp2fc(arrayIn, fpwork, sptrans->poli, nlev, nlat, nfc, trunc);
  fc2gp(sptrans->trig, sptrans->ifax, fpwork, arrayOut, nlat, nlon, nlev, nfc);

  free(fpwork);
}


void spec2spec(int gridIDin, double *arrayIn, int gridIDout, double *arrayOut)
{
  int truncIn, truncOut;

  truncIn  = gridInqTrunc(gridIDin);
  truncOut = gridInqTrunc(gridIDout);

  sp2sp(arrayIn, truncIn, arrayOut, truncOut);
}


void speccut(int gridIDin, double *arrayIn, double *arrayOut, int waves[])
{
  int trunc;

  trunc = gridInqTrunc(gridIDin);

  spcut(arrayIn, arrayOut, trunc, waves);
}


SPTRANS *sptrans_new(int nlon, int nlat, int trunc)
{
  static char func[] = "sptrans_new";
  SPTRANS *sptrans;
  int dimsp;

  sptrans = (SPTRANS *) malloc(sizeof(SPTRANS));

  sptrans->nlon   = nlon;
  sptrans->nlat   = nlat;
  sptrans->trunc  = trunc;

  dimsp = (trunc + 1)*(trunc + 2);
  sptrans->poldim = dimsp / 2 * nlat;

  sptrans->trig = (double *) malloc(nlon * sizeof(double));
  fft_set(sptrans->trig, sptrans->ifax, nlon);

  sptrans->poli = (double *) malloc(sptrans->poldim * sizeof(double));
  sptrans->pold = (double *) malloc(sptrans->poldim * sizeof(double));

  sptrans->rcoslat = (double *) malloc(nlat * sizeof(double));

  legini(trunc, nlat, sptrans->poli, sptrans->pold, sptrans->rcoslat);

  return (sptrans);
}


void sptrans_delete(SPTRANS *sptrans)
{
  static char func[] = "sptrans_delete";

  if ( sptrans )
    {
      if ( sptrans->trig ) { free(sptrans->trig);  sptrans->trig = NULL; }
      if ( sptrans->poli ) { free(sptrans->poli);  sptrans->poli = NULL; }
      if ( sptrans->pold ) { free(sptrans->pold);  sptrans->pold = NULL; }

      free(sptrans); sptrans = NULL;
    }
}


DVTRANS *dvtrans_new(int ntr)
{
  static char func[] = "dvtrans_new";
  DVTRANS *dvtrans;
  int dimsp;

  dvtrans = (DVTRANS *) malloc(sizeof(DVTRANS));

  dvtrans->trunc  = ntr;

  dimsp = (ntr + 1)*(ntr + 2);
  dvtrans->fdim = dimsp / 2;

  dvtrans->f1 = (double *) malloc(dvtrans->fdim * sizeof(double));
  dvtrans->f2 = (double *) malloc(dvtrans->fdim * sizeof(double));

  geninx(ntr, dvtrans->f1, dvtrans->f2);

  return (dvtrans);
}


void dvtrans_delete(DVTRANS *dvtrans)
{
  static char func[] = "dvtrans_delete";

  if ( dvtrans )
    {
      if ( dvtrans->f1 ) { free(dvtrans->f1);  dvtrans->f1 = NULL; }
      if ( dvtrans->f2 ) { free(dvtrans->f2);  dvtrans->f2 = NULL; }

      free(dvtrans); dvtrans = NULL;
    }
}


void dv2uv(double *d, double *o, double *u, double *v, double *f, double *g,
	   int nt, int nsp, int nlev)
{
  /* d(nsp,nlev), o(nsp,nlev)     ! divergence, vorticity        */
  /* u(nsp,nlev), v(nsp,nlev)     ! zonal wind, meridional wind  */
  /* f(nsp/2)   , g(nsp/2)        ! factor tables                */

  int l, m, n;
  int i;

  for ( l = 0; l < nlev; l++ )
    {
      i = 0;

      for ( m = 0; m < nt-1; m++ )
	{
	  /*********/
	  /* n = m */
	  /*********/

	  if ( m == 0 )
	    {
	      *u++ = -g[i+1] * o[2*(i+1)  ];
	      *u++ = -g[i+1] * o[2*(i+1)+1];
	      *v++ =  g[i+1] * d[2*(i+1)  ];
	      *v++ =  g[i+1] * d[2*(i+1)+1];
	    }
	  else
	    {
	      *u++ = -f[i] * d[2*i+1] - g[i+1] * o[2*(i+1)  ];
	      *u++ =  f[i] * d[2*i  ] - g[i+1] * o[2*(i+1)+1];
	      *v++ = -f[i] * o[2*i+1] + g[i+1] * d[2*(i+1)  ];
	      *v++ =  f[i] * o[2*i  ] + g[i+1] * d[2*(i+1)+1];
	    }
	  ++i;

	  /****************/
	  /* m < n < nt-1 */
	  /****************/

	  for ( n = m+1; n < nt-1; n++ )
	    {
	      *u++ =  g[i] * o[2*(i-1)  ] - f[i] * d[2*i+1] - g[i+1] * o[2*(i+1)  ];
	      *u++ =  g[i] * o[2*(i-1)+1] + f[i] * d[2*i  ] - g[i+1] * o[2*(i+1)+1];
	      *v++ = -g[i] * d[2*(i-1)  ] - f[i] * o[2*i+1] + g[i+1] * d[2*(i+1)  ];
	      *v++ = -g[i] * d[2*(i-1)+1] + f[i] * o[2*i  ] + g[i+1] * d[2*(i+1)+1];
	      ++i;
	    }

	  /************/
	  /* n = nt-1 */
	  /************/

	  *u++ =  g[i] * o[2*(i-1)  ] - f[i] * d[2*i+1];
	  *u++ =  g[i] * o[2*(i-1)+1] + f[i] * d[2*i  ];
	  *v++ = -g[i] * d[2*(i-1)  ] - f[i] * o[2*i+1];
	  *v++ = -g[i] * d[2*(i-1)+1] + f[i] * o[2*i  ];
	  ++i;

	  /**********/
	  /* n = nt */
	  /**********/

	  *u++ =  g[i] * o[2*(i-1)  ];
	  *u++ =  g[i] * o[2*(i-1)+1];
	  *v++ = -g[i] * d[2*(i-1)  ];
	  *v++ = -g[i] * d[2*(i-1)+1];
	  ++i;
	}

      /***************************/
      /* m = nt-1  and  n = nt-1 */
      /***************************/

      *u++ = -f[i] * d[2*i+1];
      *u++ =  f[i] * d[2*i  ];
      *v++ = -f[i] * o[2*i+1];
      *v++ =  f[i] * o[2*i  ];
      ++i;

      /*************************/
      /* m = nt-1  and  n = nt */
      /*************************/

      *u++ =  g[i] * o[2*(i-1)  ];
      *u++ =  g[i] * o[2*(i-1)+1];
      *v++ = -g[i] * d[2*(i-1)  ];
      *v++ = -g[i] * d[2*(i-1)+1];
      ++i;

      /***********************/
      /* m = nt  and  n = nt */
      /***********************/

      *u++ = 0.0;
      *u++ = 0.0;
      *v++ = 0.0;
      *v++ = 0.0;

      d += nsp;
      o += nsp;
    }
}


void scaluv(double *fu, double *rclat, int nlat, int lot)
{
  int l,lat;

  for (l = 0; l < lot; l++)
    for (lat = 0; lat < nlat; lat++)
      {
        *fu *= rclat[lat];
        fu++;
      }
}


void trans_dv2uv(SPTRANS *sptrans, DVTRANS *dvtrans, int nlev,
		 int gridID1, double *sd, double *svo,
		 int gridID2, double *gu, double *gv)
{
  static char func[] = "spec2grid";
  int ntr, nlat, nlon, nfc;
  int waves;
  int dimsp;
  double *fpwork;
  double *su, *sv;
    
  ntr  = gridInqTrunc(gridID1);
  nlon = gridInqXsize(gridID2);
  nlat = gridInqYsize(gridID2);

  waves = ntr + 1;
  nfc   = waves * 2;

  dimsp = (ntr + 1)*(ntr + 2);

  su = gu;
  sv = gv;

  dv2uv(sd, svo, su, sv, dvtrans->f1, dvtrans->f2, ntr, dimsp, nlev);

  fpwork = (double *) malloc(nlat*nfc*nlev*sizeof(double));

  sp2fc(su, fpwork, sptrans->poli, nlev, nlat, nfc, ntr);
  scaluv(fpwork, sptrans->rcoslat, nlat, nfc*nlev);
  fc2gp(sptrans->trig, sptrans->ifax, fpwork, gu, nlat, nlon, nlev, nfc);

  sp2fc(sv, fpwork, sptrans->poli, nlev, nlat, nfc, ntr);
  scaluv(fpwork, sptrans->rcoslat, nlat, nfc*nlev);
  fc2gp(sptrans->trig, sptrans->ifax, fpwork, gv, nlat, nlon, nlev, nfc);

  free(fpwork);
}

