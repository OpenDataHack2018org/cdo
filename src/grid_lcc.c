#include <stdio.h>
#include <math.h>

#ifndef  M_PI
#define  M_PI		3.14159265358979323846	/* pi */
#endif

int W3FB12(double xi, double xj, double alat1, double elon1, double dx,
	   double elonv, double alatan, double *alat, double *elon)
{
  /*
    SUBPROGRAM  DOCUMENTATION  BLOCK

    SUBPROGRAM:  W3FB12        LAMBERT(I,J) TO LAT/LON FOR GRIB

    PRGMMR: STACKPOLE        ORG: NMC42       DATE:88-11-28

    ABSTRACT:
      CONVERTS THE COORDINATES OF A LOCATION ON EARTH GIVEN IN A
      GRID COORDINATE SYSTEM OVERLAID ON A LAMBERT CONFORMAL TANGENT
      CONE PROJECTION TRUE AT A GIVEN N OR S LATITUDE TO THE
      NATURAL COORDINATE SYSTEM OF LATITUDE/LONGITUDE
      W3FB12 IS THE REVERSE OF W3FB11.
      USES GRIB SPECIFICATION OF THE LOCATION OF THE GRID

    PROGRAM HISTORY LOG:
      1988-11-25  ORIGINAL AUTHOR:   Stackpole, W/NMC42
      2007-12-24  CONVERT TO ANSI C: Uwe Schulzweida, MPIMET

    ierr = W3FB12(xi,xj,alat1,elon1,dx,elonv,alatan,&alat,&elon)

    INPUT ARGUMENT LIST:
      xi       - I COORDINATE OF THE POINT
      xj       - J COORDINATE OF THE POINT
      alat1    - LATITUDE  OF LOWER LEFT POINT OF GRID (POINT 1,1)
                 LATITUDE <0 FOR SOUTHERN HEMISPHERE;
      elon1    - LONGITUDE OF LOWER LEFT POINT OF GRID (POINT 1,1)
                   EAST LONGITUDE USED THROUGHOUT;
      dx       - MESH LENGTH OF GRID IN METERS AT TANGENT LATITUDE
      elonv    - THE ORIENTATION OF THE GRID.  I.E.,
                 THE EAST LONGITUDE VALUE OF THE VERTICAL MERIDIAN
                 WHICH IS PARALLEL TO THE Y-AXIS (OR COLUMNS OF
                 THE GRID) ALONG WHICH LATITUDE INCREASES AS
                 THE Y-COORDINATE INCREASES.
                 THIS IS ALSO THE MERIDIAN (ON THE OTHER SIDE OF THE
                 TANGENT CONE) ALONG WHICH THE CUT IS MADE TO LAY
                 THE CONE FLAT.
      alatan   - THE LATITUDE AT WHICH THE LAMBERT CONE IS TANGENT TO
                 (TOUCHES OR OSCULATES) THE SPHERICAL EARTH.
                  SET NEGATIVE TO INDICATE A
                  SOUTHERN HEMISPHERE PROJECTION;

    OUTPUT ARGUMENT LIST:
      alat     - LATITUDE IN DEGREES (NEGATIVE IN SOUTHERN HEMI.)
      elon     - EAST LONGITUDE IN DEGREES
      ierr     - .EQ. 0   IF NO PROBLEM
                 .GE. 1   IF THE REQUESTED XI,XJ POINT IS IN THE
                          FORBIDDEN ZONE, I.E. OFF THE LAMBERT MAP
                          IN THE OPEN SPACE WHERE THE CONE IS CUT.
                   IF IERR.GE.1 THEN ALAT=999. AND ELON=999.

    REMARKS:
      FORMULAE AND NOTATION LOOSELY BASED ON HOKE, HAYES,
      AND RENNINGER'S "MAP PROJECTIONS AND GRID SYSTEMS...", MARCH 1981
      AFGWC/TN-79/003
  */

  int ierr = 0;
  const double REARTH = 6371200.0;
  /* const double REARTH = 6367470.0; *//* GRIB !!! */
  double h, piby2, radpd, degprd, rebydx, alatn1, an, cosltn;
  double elon1l, ala1, rmll, elo1, arg, polei, polej;
  double xx, yy, r2, theta, beta, elonvr, aninv, aninv2, thing;

  /*
      PRELIMINARY VARIABLES AND REDIFINITIONS

      H = 1 FOR NORTHERN HEMISPHERE; = -1 FOR SOUTHERN
  */
  if ( alatan > 0 ) h =  1.;
  else              h = -1.;

  piby2  = M_PI/2.;
  radpd  = M_PI/180.0;
  degprd = 1./radpd;
  rebydx = REARTH/dx;
  alatn1 = alatan * radpd;
  an     = h * sin(alatn1);
  cosltn = cos(alatn1);
  /*
    MAKE SURE THAT INPUT LONGITUDE DOES NOT PASS THROUGH
    THE CUT ZONE (FORBIDDEN TERRITORY) OF THE FLAT MAP
    AS MEASURED FROM THE VERTICAL (REFERENCE) LONGITUDE
  */
  elon1l = elon1;
  if ( (elon1-elonv) >  180. ) elon1l = elon1 - 360.;
  if ( (elon1-elonv) < -180. ) elon1l = elon1 + 360.;

  elonvr = elonv * radpd;
  /*
    RADIUS TO LOWER LEFT HAND (LL) CORNER
  */
  ala1 = alat1 * radpd;
  rmll = rebydx * (pow(cosltn,(1.-an))*pow((1.+an),an)) *
         pow(((cos(ala1))/(1.+h*sin(ala1))),an)/an;

  /* USE LL POINT INFO TO LOCATE POLE POINT */
  elo1 = elon1l * radpd;
  arg = an * (elo1-elonvr);
  polei = 1. - h * rmll * sin(arg);
  polej = 1. + rmll * cos(arg);

  /*
    RADIUS TO THE I,J POINT (IN GRID UNITS)
    YY REVERSED SO POSITIVE IS DOWN
  */
  xx = xi - polei;
  yy = polej - xj;
  r2 = xx*xx + yy*yy;
  /*
    CHECK THAT THE REQUESTED I,J IS NOT IN THE FORBIDDEN ZONE
    YY MUST BE POSITIVE UP FOR THIS TEST
  */
  theta = M_PI*(1.-an);
  beta = fabs(atan2(xx,-yy));
  if ( beta <= theta )
    {
      ierr = 1;
      *alat = 999.;
      *elon = 999.;
      /*  return (ierr); */
    }
  /*
    NOW THE MAGIC FORMULAE
  */
  if ( !(fabs(r2) > 0.) )
    {
      *alat = h * 90.;
      *elon = elonv;
    }
  else
    {
      /*
	FIRST THE LONGITUDE
      */
      *elon = elonv + degprd * atan2(h*xx,yy)/an;
      *elon = fmod(*elon+360., 360.);
      /*
	NOW THE LATITUDE
	RECALCULATE THE THING ONLY IF MAP IS NEW SINCE LAST TIME
      */
      aninv  = 1./an;
      aninv2 = aninv/2.;
      thing  = pow((an/rebydx),aninv)/ (pow(cosltn,((1.-an)*aninv))*(1.+ an));

      *alat  = h*(piby2 - 2.*atan(thing*pow(r2,aninv2)))*degprd;
    }
  /*
    FOLLOWING TO ASSURE ERROR VALUES IF FIRST TIME THRU
    IS OFF THE MAP
  */
  if ( ierr != 0 )
    {
      *alat = 999.;
      *elon = 999.;
      ierr = 2;
    }

  return(ierr);
}

/*
int main(void)
{
  int    status;
  int    nlon = 245;
  int    nlat = 277;
  double xi, xj;
  double lat_ll_p   =  47.806;
  double lon_ll_p   = -10.063;
  double lat_tan_p  =  59.2;
  double dx_p       =  11000.0;
  double lon_xx_p   = -10.0;
  double zlat, zlon;

  xi = 1;
  xj = 1;
  status =  W3FB12 (xi, xj,lat_ll_p,360.+lon_ll_p,dx_p, 360.+lon_xx_p,lat_tan_p,&zlat,&zlon);
  printf("1 1 47.806 349.937 0\n");
  printf("%g %g %g %g %d\n", xj, xi, zlat, zlon, status);

  xi = nlon;
  xj = nlat;
  status =  W3FB12 (xi, xj,lat_ll_p,360.+lon_ll_p,dx_p, 360.+lon_xx_p,lat_tan_p,&zlat,&zlon);
  printf("277 245 63.086 51.4192 0\n");
  printf("%g %g %g %g %d\n", xj, xi, zlat, zlon, status);

  {
    int i, j;
    for ( j = 1; j <= nlat; j++ )
      for ( i = 1; i <= nlon; i++ )
	{
	  xi = i;
	  xj = j;
	  status =  W3FB12 (xi, xj,lat_ll_p,360.+lon_ll_p,dx_p, 360.+lon_xx_p,lat_tan_p,&zlat,&zlon);
	}
  }

  return (0);
}
*/
