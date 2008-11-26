#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI		3.14159265358979323846	/* pi */
#endif

#ifndef  rad2deg
#define  rad2deg  (180./M_PI)   /* conversion for rad to deg */
#endif

#ifndef  deg2rad
#define  deg2rad  (M_PI/180.)   /* conversion for deg to rad */
#endif


double rls_to_rl(double phis, double rlas, double polphi, double pollam)
{
  /*
    Umrechnung von rlas (rot. System) auf rla (geo. System)

    phis   : Breite im rotierten System (N>0)
    rlas   : Laenge im rotierten System (E>0)
    polphi : Geographische Breite des Nordpols des rot. Systems
    pollam : Geographische Laenge des Nordpols des rot. Systems

    result : Geographische Laenge
  */
  double zsinpol, zcospol, zlampol;
  double zphis, zrlas, zarg1, zarg2;

  zsinpol = sin(deg2rad*polphi);
  zcospol = cos(deg2rad*polphi);

  zlampol = deg2rad*pollam;
  zphis   = deg2rad*phis;
  if ( rlas > 180.0 ) rlas -= 360.0;
  zrlas   = deg2rad*rlas;

  zarg1   = sin(zlampol)*(- zsinpol*cos(zrlas)*cos(zphis)  +
                            zcospol*           sin(zphis)) -
            cos(zlampol)*           sin(zrlas)*cos(zphis);
  zarg2   = cos(zlampol)*(- zsinpol*cos(zrlas)*cos(zphis)  +
                            zcospol*           sin(zphis)) +
            sin(zlampol)*           sin(zrlas)*cos(zphis);

  if ( fabs(zarg2) < 1.0e-20 ) zarg2 = 1.0e-20;

  return (rad2deg*atan2(zarg1, zarg2));
}


double phs_to_ph(double phis, double rlas, double polphi)
{
  /*
    Umrechnung von phis (rot. System) auf phi (geo. System)

    phis   : Breite im rotierten System (N>0)
    rlas   : Laenge im rotierten System (E>0)
    polphi : Geographische Breite des Nordpols des rot. Systems

    result : Geographische Breite
  */
  double zsinpol, zcospol;
  double zphis, zrlas, zarg;

  zsinpol = sin(deg2rad*polphi);
  zcospol = cos(deg2rad*polphi);

  zphis   = deg2rad*phis;
  if ( rlas > 180.0 ) rlas -= 360.0;
  zrlas   = deg2rad*rlas;

  zarg    = zcospol*cos(zphis)*cos(zrlas) + zsinpol*sin(zphis);

  return (rad2deg*asin(zarg));
}


double rl_to_rls(double phi, double rla, double polphi, double pollam)
{
  /*
    Umrechnung von rla (geo. System) auf rlas (rot. System)

    phi    : Breite im geographischen System (N>0)
    rla    : Laenge im geographischen System (E>0)
    polphi : Geographische Breite des Nordpols des rot. Systems
    pollam : Geographische Laenge des Nordpols des rot. Systems

    result : Rotierte Laenge
  */
  double zsinpol, zcospol, zlampol;
  double zphi, zrla, zarg1, zarg2;

  zsinpol = sin(deg2rad*polphi);
  zcospol = cos(deg2rad*polphi);
  zlampol =     deg2rad*pollam;

  if ( rla > 180.0 ) rla -= 360.0;

  zrla = deg2rad*rla;
  zphi = deg2rad*phi;

  zarg1  = - sin(zrla-zlampol)*cos(zphi);
  zarg2  = - zsinpol*cos(zphi)*cos(zrla-zlampol)+zcospol*sin(zphi);

  if ( fabs(zarg2) < 1.0e-20 ) zarg2 = 1.0e-20;

  return (rad2deg*atan2(zarg1,zarg2));
}


double ph_to_phs(double phi, double rla, double polphi, double pollam)
{
  /*
    Umrechnung von phi (geo. System) auf phis (rot. System)

    phi    : Breite im geographischen System (N>0)
    rla    : Laenge im geographischen System (E>0)
    polphi : Geographische Breite des Nordpols des rot. Systems
    pollam : Geographische Laenge des Nordpols des rot. Systems

    result : Rotierte Breite
  */
  double zsinpol, zcospol, zlampol;
  double zphi, zrla, zarg;

  zsinpol = sin(deg2rad*polphi);
  zcospol = cos(deg2rad*polphi);
  zlampol =     deg2rad*pollam;

  zphi = deg2rad*phi;
  if ( rla > 180.0 ) rla -= 360.0;
  zrla = deg2rad*rla;

  zarg = zcospol*cos(zphi)*cos(zrla-zlampol) + zsinpol*sin(zphi);

  return (rad2deg*asin(zarg));
}


void usvs_to_uv(double us, double vs, double phi, double rla,
		double polphi, double pollam, double *u, double *v)
{
  /*
    Umrechnen der windkomponente us, vs im rotierten sphaerischen
    system in die windkomponenten u, v, im geographischen system

    us     : 'zonaler wind im rotierten system
    vs     : 'merid. wind im rotierten  system
    phi    : breite im geographischen system (n>0)
    rla    : laenge im geographischen system (e>0)
    polphi : geographische breite des n-pols des rot. sys.
    pollam : geographische laenge des n-pols des rot. sys.
 
    u      : zonaler wind im geographischen system
    v      : merid. wind im geographischen system
  */
  double zpolphi, zpollam, zrla, zphi, pollamd, zrlas, zarg, zbeta;

  /* umrechnung von grad in bogenmass */
  zpolphi = polphi*deg2rad;
  zpollam = pollam*deg2rad;
  zrla    = rla   *deg2rad;
  zphi    = phi   *deg2rad;
  pollamd = pollam;
  if ( pollamd < 0.0 ) pollamd += 360.0;

  /* laenge im rotierten system berechnen */
  zrlas = rl_to_rls(phi, rla, polphi, pollam)*deg2rad;

  /* winkel zbeta berechen (schnittwinkel der breitenkreise) */
  zarg = - sin(zpolphi)*sin(zrla-zpollam)*sin(zrlas) - cos(zrla-zpollam)*cos(zrlas);
  if ( zarg >  1.0 ) zarg =  1.0;
  if ( zarg < -1.0 ) zarg = -1.0;
  /*
  zbeta = acos(zarg);
  zbeta = sign(zbeta, -(rla - (pollamd-180.0)));
  */
  zbeta = fabs(acos(zarg));
  /*  if ( -(rla - (pollamd-180.0)) < 0 ) zbeta = -zbeta;*/
  if ( (-(rla - (pollamd-180.0)) < 0) && (-(rla - (pollamd-180.0)) >= -180) ) zbeta = -zbeta;

  /* us - wind transformieren */
  *u = us*cos(zbeta) - vs*sin(zbeta);
  
  /* vs - wind transformieren */
  *v = us*sin(zbeta) + vs*cos(zbeta);
}

/*
int main(void)
{
  double polphi, pollam;
  double x0, y0, x1, y1, x2, y2;
  int i;

  polphi = 90.0;
  pollam = 0.0;

  polphi = 32.5;
  pollam = -170.0;

  x0 = -20.0;
  y0 = 0.0;

  for ( i = 0; i < 10; i++ )
    {
      x0 = i *20.0;
      printf("%g %g\n", x0, y0);

      x1 = rls_to_rl(y0, x0, polphi, pollam);
      y1 = phs_to_ph(y0, x0, polphi);

      printf("%g %g\n", x1, y1);

      x2 = rl_to_rls(y1, x1, polphi, pollam);
      y2 = ph_to_phs(y1, x1, polphi, pollam);

      printf("%g %g\n", x2, y2);
    }

  usvs_to_uv(30.0, 20.0, 30.0, 0.0, polphi, pollam, &x1, &x2);
  printf("usvs_to_uv %g %g %g %g\n", polphi, pollam, x1, x2);

  return (0);
}
*/
