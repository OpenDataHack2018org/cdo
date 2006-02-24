/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2006 Uwe Schulzweida, schulzweida@dkrz.de
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
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

/*
@BeginDoc

@BeginModule

@Name      = Specinfo
@Title     = Sectral information
@Section   = Information
@Class     = Information
@Arguments = 
@Operators = specinfo

@EndModule


@BeginOperator_specinfo

@Title     = Spectral information

@BeginDesciption
@EndDesciption

@EndOperator

@EndDoc
*/

#define NSP2NTR(nsp) ((int) ((((sqrt((double)(4*nsp+1)))-3)/2)))


static void fac(int nlonin, int *nlonout, int *ierr)
{
  int n2, n3, n5;
  int m;
  
  n2 = 0;
  n3 = 0;
  n5 = 0;

  m = nlonin;

  while (m%2 == 0)
    {
      m = m/2;
      n2++;
    }
  while (m%3 == 0)
    {
      m = m/3;
      n3++;
    }
  while (m%5 == 0)
    {
      m = m/5;
      n5++;
    }

  if (m == 1) {
    *nlonout = nlonin;
    *ierr = 0;
  } else {
    *nlonout =  nlonin+1;
    *ierr = 1;
  }

  return;
}


static int nlat2nlon(int nlat)
{
  int nlon, m, ierr;

  if ( nlat == 0 )
    cdoAbort("nlat = 0!");

  nlon = 2*nlat;

  fac(nlon, &m, &ierr);
  /* adjust till fft is possible */
  while (ierr != 0) 
    {
      nlon = m;
      /* correct here nlon so that nlat keeps always even */
      while (nlon%4 != 0) nlon++;
      fac(nlon, &m, &ierr);
    }

  return (nlon);
}


void *Specinfo(void *argument)
{
  char arg[128], *parg;
  int len, i, nout1 = 0, nout2 = 0;
  int ntr1 = 0, nsp1 = 0, nlat1 = 0, nlon1 = 0, ngp1 = 0, ni1 = 0, ngp_gme1 = 0;
  int ntr2 = 0, nsp2 = 0, nlat2 = 0, nlon2 = 0, ngp2 = 0, ni2 = 0, ngp_gme2 = 0;

  cdoInitialize(argument);

  operatorInputArg("TX, TLX, NLON=X, NLAT=X or NIX");

  len = strlen(operatorArgv()[0]);

  if ( (len+1) >= 128 ) cdoAbort("Parameter string to large!");

  for ( i = 0; i < len; i++ ) arg[i] = toupper(operatorArgv()[0][i]);
  arg[len] = 0;

  if ( arg[0] == 'T' && arg[1] == 'L' )
    {
      parg = &arg[2];
      if ( *parg == '=' ) *parg++;
      if ( ! isdigit((int) *parg) ) cdoAbort("Wrong parameter: %s", arg);
      ntr2   = atoi(parg);
      nlat2  = ntr2nlat_linear(ntr2);
      nlon2  = compNlon(nlat2);
      ngp2   = nlon2*nlat2;
      ni2    = (int) sqrt((double)ngp2/10) - 1;
      nout1  = FALSE;
      nout2  = TRUE;
    }
  else if ( arg[0] == 'T' )
    {
      parg = &arg[1];
      if ( *parg == '=' ) *parg++;
      if ( ! isdigit((int) *parg) ) cdoAbort("Wrong parameter: %s", arg);
      ntr1   = atoi(parg);
      nlat1  = ntr2nlat(ntr1);
      nlon1  = compNlon(nlat1);
      ngp1   = nlon1*nlat1;
      ni1    = (int) sqrt((double)ngp1/10) - 1;
      nout1  = TRUE;
      nout2  = FALSE;
    }
  else if ( arg[0] == 'N' && arg[1] == 'I' )
    {
      parg = &arg[2];
      if ( *parg == '=' ) *parg++;
      if ( ! isdigit((int) *parg) ) cdoAbort("Wrong parameter: %s", arg);
      ni1    = atoi(parg);
      ni2    = ni1;
      ngp_gme1 = (ni1+1)*(ni1+1)*10;
      ngp_gme2 = (ni2+1)*(ni2+1)*10;
      nsp1   = ngp_gme1 / 2;
      nsp2   = ngp_gme2 / 2;
      ntr1   = NSP2NTR(nsp1);
      ntr2   = NSP2NTR(nsp2);
      nlat1  = ntr1 + ntr1%2;
      nlat2  = ntr2 + ntr2%2;
      nlon1  = nlat2nlon(nlat1);
      nlon2  = nlat2nlon(nlat2);
      nlat1  = nlon1 / 2;
      nlat2  = nlon2 / 2;
      ntr1   = (nlat1*2-1)/3;
      ntr2   = (nlat2*2-1)/2;
      nout1  = TRUE;
      nout2  = TRUE;
    }
  else if ( arg[0] == 'N' && arg[1] == 'L' && arg[2] == 'A' && arg[3] == 'T' )
    {
      parg = &arg[4];
      if ( *parg == '=' ) *parg++;
      if ( ! isdigit((int) *parg) ) cdoAbort("Wrong parameter: %s", arg);
      nlat1  = atoi(parg);
      nlat2  = nlat1;
      nlon1  = nlat2nlon(nlat1);
      nlon2  = nlat2nlon(nlat2);
      nlat1  = nlon1 / 2;
      nlat2  = nlon2 / 2;
      ntr1   = (nlat1*2-1)/3;
      ntr2   = (nlat2*2-1)/2;
      ngp1   = nlon1*nlat1;
      ngp2   = nlon2*nlat2;
      ni1    = (int) sqrt((double)ngp1/10) - 1;
      ni2    = (int) sqrt((double)ngp2/10) - 1;
      nout1  = TRUE;
      nout2  = TRUE;
    }
  else if ( arg[0] == 'N' && arg[1] == 'L' && arg[2] == 'O' && arg[3] == 'N' )
    {
      parg = &arg[4];
      if ( *parg == '=' ) *parg++;
      if ( ! isdigit((int) *parg) ) cdoAbort("Wrong parameter: %s", arg);
      nlon1  = atoi(parg);
      nlon2  = nlon1;
      nlat1  = nlon1 / 2;
      nlat2  = nlon2 / 2;
      nlon1  = nlat2nlon(nlat1);
      nlon2  = nlat2nlon(nlat2);
      nlat1  = nlon1 / 2;
      nlat2  = nlon2 / 2;
      ntr1   = (nlat1*2-1)/3;
      ntr2   = (nlat2*2-1)/2;
      ngp1   = nlon1*nlat1;
      ngp2   = nlon2*nlat2;
      ni1    = (int) sqrt((double)ngp1/10) - 1;
      ni2    = (int) sqrt((double)ngp2/10) - 1;
      nout1  = TRUE;
      nout2  = TRUE;
    }
  else
    cdoAbort("Unsupported parameter: %s", arg);

  nsp1     = (ntr1+1)*(ntr1+2);
  nsp2     = (ntr2+1)*(ntr2+2);
  ngp1     = nlon1*nlat1;
  ngp2     = nlon2*nlat2;
  ngp_gme1 = (ni1+1)*(ni1+1)*10;
  ngp_gme2 = (ni2+1)*(ni2+1)*10;

  fprintf(stdout, "truncation     nsp  nlat  nlon      ngp   ni  ngp_gme\n");

  if ( nout1 ) fprintf(stdout, "   T%-4d  %8d %5d %5d %8d %4d %8d\n",
		       ntr1, nsp1, nlat1, nlon1, ngp1, ni1, ngp_gme1);

  if ( nout2 ) fprintf(stdout, "   TL%-4d %8d %5d %5d %8d %4d %8d\n",
		       ntr2, nsp2, nlat2, nlon2, ngp2, ni2, ngp_gme2);

  cdoFinish();

  return (0);
}
