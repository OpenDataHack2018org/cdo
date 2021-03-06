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
#include <math.h>

#include "grid.h"

double
lamrot_to_lam(double phirot, double lamrot, double polphi, double pollam, double polgam)
{
  /*
    Name of the original Fortran function: PHTOPHS

    This function converts lambda from one rotated system to lambda in another
    system. If the optional argument polgam is present, the other system can
    also be a rotated one, where polgam is the angle between the two north
    poles. If polgam is not present, the other system is the real geographical
    system.

    phirot : latitude in the rotated system
    lamrot : longitude in the rotated system (E>0)
    polphi : latitude of the rotated north pole
    pollam : longitude of the rotated north pole

    result : longitude in the geographical system
  */
  double zarg1, zarg2;

  double zsinpol = sin(DEG2RAD * polphi);
  double zcospol = cos(DEG2RAD * polphi);

  double zlampol = DEG2RAD * pollam;
  double zphirot = DEG2RAD * phirot;
  if (lamrot > 180.0) lamrot -= 360.0;
  double zlamrot = DEG2RAD * lamrot;

  if (polgam > 0)
    {
      double zgam = DEG2RAD * polgam;
      zarg1 = sin(zlampol)
                  * (-zsinpol * cos(zphirot) * (cos(zlamrot) * cos(zgam) - sin(zlamrot) * sin(zgam)) + zcospol * sin(zphirot))
              - cos(zlampol) * cos(zphirot) * (sin(zlamrot) * cos(zgam) + cos(zlamrot) * sin(zgam));

      zarg2 = cos(zlampol)
                  * (-zsinpol * cos(zphirot) * (cos(zlamrot) * cos(zgam) - sin(zlamrot) * sin(zgam)) + zcospol * sin(zphirot))
              + sin(zlampol) * cos(zphirot) * (sin(zlamrot) * cos(zgam) + cos(zlamrot) * sin(zgam));
    }
  else
    {
      zarg1 = sin(zlampol) * (-zsinpol * cos(zlamrot) * cos(zphirot) + zcospol * sin(zphirot))
              - cos(zlampol) * sin(zlamrot) * cos(zphirot);
      zarg2 = cos(zlampol) * (-zsinpol * cos(zlamrot) * cos(zphirot) + zcospol * sin(zphirot))
              + sin(zlampol) * sin(zlamrot) * cos(zphirot);
    }

  double result = 0;
  if (fabs(zarg2) > 0) result = RAD2DEG * atan2(zarg1, zarg2);
  if (fabs(result) < 9.e-14) result = 0;

  return result;
}

double
phirot_to_phi(double phirot, double lamrot, double polphi, double polgam)
{
  /*
    Name of the original Fortran function: PHSTOPH

    This function converts phi from one rotated system to phi in another
    system. If the optional argument polgam is present, the other system
    can also be a rotated one, where polgam is the angle between the two
    north poles.
    If polgam is not present, the other system is the real geographical
    system.

    phirot : latitude in the rotated system
    lamrot : longitude in the rotated system (E>0)
    polphi : latitude of the rotated north pole
    polgam : angle between the north poles of the systems

    result : latitude in the geographical system
  */
  double zarg;

  double zsinpol = sin(DEG2RAD * polphi);
  double zcospol = cos(DEG2RAD * polphi);

  double zphirot = DEG2RAD * phirot;
  if (lamrot > 180.0) lamrot -= 360.0;
  double zlamrot = DEG2RAD * lamrot;

  if (polgam > 0)
    {
      double zgam = DEG2RAD * polgam;
      zarg = zsinpol * sin(zphirot) + zcospol * cos(zphirot) * (cos(zlamrot) * cos(zgam) - sin(zgam) * sin(zlamrot));
    }
  else
    zarg = zcospol * cos(zphirot) * cos(zlamrot) + zsinpol * sin(zphirot);

  return RAD2DEG * asin(zarg);
}

static double
lam_to_lamrot(double phi, double rla, double polphi, double pollam)
{
  /*
    Name of the original Fortran function: RLSTORL

    Umrechnung von rla (geo. System) auf rlas (rot. System)

    phi    : Breite im geographischen System (N>0)
    rla    : Laenge im geographischen System (E>0)
    polphi : Geographische Breite des Nordpols des rot. Systems
    pollam : Geographische Laenge des Nordpols des rot. Systems

    result : Rotierte Laenge
  */
  double zsinpol = sin(DEG2RAD * polphi);
  double zcospol = cos(DEG2RAD * polphi);
  double zlampol = DEG2RAD * pollam;

  if (rla > 180.0) rla -= 360.0;

  double zrla = DEG2RAD * rla;
  double zphi = DEG2RAD * phi;

  double zarg1 = -sin(zrla - zlampol) * cos(zphi);
  double zarg2 = -zsinpol * cos(zphi) * cos(zrla - zlampol) + zcospol * sin(zphi);

  if (fabs(zarg2) < 1.0e-20) zarg2 = 1.0e-20;

  return RAD2DEG * atan2(zarg1, zarg2);
}

#ifdef TEST_GRID_ROT
static double
phi_to_phirot(double phi, double lam, double polphi, double pollam)
{
  /*
    Name of the original Fortran function: PHTOPHS

    Umrechnung von phi (geo. System) auf phis (rot. System)

    phi    : Breite im geographischen System (N>0)
    lam    : Laenge im geographischen System (E>0)
    polphi : Geographische Breite des Nordpols des rot. Systems
    pollam : Geographische Laenge des Nordpols des rot. Systems

    result : Rotierte Breite
  */
  double zsinpol = sin(DEG2RAD * polphi);
  double zcospol = cos(DEG2RAD * polphi);
  double zlampol = DEG2RAD * pollam;

  double zphi = DEG2RAD * phi;
  if (lam > 180.0) lam -= 360.0;
  double zlam = DEG2RAD * lam;

  double zarg = zcospol * cos(zphi) * cos(zlam - zlampol) + zsinpol * sin(zphi);

  return RAD2DEG * asin(zarg);
}
#endif

void
usvs_to_uv(double us, double vs, double phi, double rla, double polphi, double pollam, double *u, double *v)
{
  /*
    Umrechnen der windkomponenten us, vs im rotierten sphaerischen
    system in die windkomponenten u, v, im geographischen system

    us     : 'zonaler wind im rotierten system
    vs     : 'merid. wind im rotierten  system
    phi    : Breite im geographischen system (N>0)
    rla    : Laenge im geographischen system (E>0)
    polphi : Geographische breite des Nordpols des rot. Systems
    pollam : Geographische laenge des Nordpols des rot. Systems

    u      : zonaler wind im geographischen system
    v      : merid. wind im geographischen system
  */
  /* umrechnung von grad in bogenmass */
  double zpolphi = polphi * DEG2RAD;
  double zpollam = pollam * DEG2RAD;
  // Added by Uwe Schulzweida (17/11/2017)
  if (pollam < 0 && rla < pollam) rla += 360.0;
  // if ( pollam < 0 && rla < 0 ) rla += 360.0;
  double zrla = rla * DEG2RAD;
  double pollamd = pollam;
  if (pollamd < 0.0) pollamd += 360.0;

  // laenge im rotierten system berechnen
  double zrlas = lam_to_lamrot(phi, rla, polphi, pollam) * DEG2RAD;

  // winkel zbeta berechen (schnittwinkel der breitenkreise)
  double zarg = -sin(zpolphi) * sin(zrla - zpollam) * sin(zrlas) - cos(zrla - zpollam) * cos(zrlas);
  if (zarg > 1.0) zarg = 1.0;
  if (zarg < -1.0) zarg = -1.0;
  /*
  zbeta = acos(zarg);
  zbeta = sign(zbeta, -(rla - (pollamd-180.0)));
  */
  double zbeta = fabs(acos(zarg));
  // if ( -(rla - (pollamd-180.0)) < 0 ) zbeta = -zbeta;
  if ((-(rla - (pollamd - 180.0)) < 0) && (-(rla - (pollamd - 180.0)) >= -180)) zbeta = -zbeta;

  // us - wind transformieren
  *u = us * cos(zbeta) - vs * sin(zbeta);

  // vs - wind transformieren
  *v = us * sin(zbeta) + vs * cos(zbeta);
}

#ifdef TEST_GRID_ROT
int
main(void)
{
  double x0, y0, x1, y1, x2, y2;
  double polphi, pollam;
  double angle = 0;

  polphi = 90.0;
  pollam = 0.0;

  polphi = 32.5;
  pollam = -170.0;

  x0 = -20.0;
  y0 = 0.0;

  for (int i = 0; i < 10; i++)
    {
      x0 = i * 20.0;
      printf("rot in: %g %g\n", x0, y0);

      x1 = lamrot_to_lam(y0, x0, polphi, pollam, angle);
      y1 = phirot_to_phi(y0, x0, polphi, angle);
      printf("geo: %g %g\n", x1, y1);

      x2 = lam_to_lamrot(y1, x1, polphi, pollam);
      y2 = phi_to_phirot(y1, x1, polphi, pollam);
      printf("rot out:%g %g\n", x2, y2);
    }

  usvs_to_uv(30.0, 20.0, 30.0, 0.0, polphi, pollam, &x1, &x2);
  printf("usvs_to_uv: %g %g %g %g\n", polphi, pollam, x1, x2);
  printf("usvs_to_uv: 32.5 -170 26.3124 24.6507 <-- reference\n");

  return 0;
}
#endif
