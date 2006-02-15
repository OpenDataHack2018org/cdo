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



void legini(int trunc, int nlat, double *poli, double *pold)
{
  static char func[] = "legini";
  int waves, dimsp, dimpnm;
  int jgl, jm, jn, is;
  int isp, jsp, latn, lats, pdim;
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

  for ( jgl = 0; jgl < nlat/2; jgl++ )
    {
      jspleg1(pnm, gmu[jgl], trunc, work);

      latn = jgl;
      isp = 0;
      jsp = 0;
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
	      jsp++;
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

  fpwork = (double *) malloc(nlat*nfc*sizeof(double));

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

  fpwork = (double *) malloc(nlat*nfc*sizeof(double));

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

  legini(trunc, nlat, sptrans->poli, sptrans->pold);

  return (sptrans);
}


void sptrans_delete(SPTRANS *sptrans)
{
  static char func[] = "sptrans_delete";

  if ( sptrans )
    {
      if ( sptrans->trig  ) { free(sptrans->trig);  sptrans->trig = NULL; }
      if ( sptrans->poli )  { free(sptrans->poli);  sptrans->poli = NULL; }
      if ( sptrans->pold )  { free(sptrans->pold);  sptrans->pold = NULL; }

      free(sptrans); sptrans = NULL;
    }
}
