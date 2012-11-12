#include <stdio.h>
#include <math.h>

#define  deg2rad  (M_PI/180.)

struct cart {
  double x[3];
};

struct geo {
  double lon;
  double lat;
};

struct cart gc2cc(struct geo *position)
{
  double cln;
  double sln;
  double clt;
  double slt;
  
  struct cart x;

  sln = sin(position->lon);
  cln = cos(position->lon);
  slt = sin(position->lat);
  clt = cos(position->lat);

  x.x[0] = cln*clt;
  x.x[1] = sln*clt;
  x.x[2] = slt;

  return (x);
}

static
double areas(struct cart *dv1, struct cart *dv2, struct cart *dv3)
{
  double a1, a2, a3;
  double ca1, ca2, ca3;
  double s12, s23, s31;

  struct cart u12, u23, u31;

  double areas;

  /* compute cross products Uij = Vi X Vj */

  u12.x[0] = dv1->x[1]*dv2->x[2] - dv1->x[2]*dv2->x[1];
  u12.x[1] = dv1->x[2]*dv2->x[0] - dv1->x[0]*dv2->x[2];
  u12.x[2] = dv1->x[0]*dv2->x[1] - dv1->x[1]*dv2->x[0];
  
  u23.x[0] = dv2->x[1]*dv3->x[2] - dv2->x[2]*dv3->x[1];
  u23.x[1] = dv2->x[2]*dv3->x[0] - dv2->x[0]*dv3->x[2];
  u23.x[2] = dv2->x[0]*dv3->x[1] - dv2->x[1]*dv3->x[0];
  
  u31.x[0] = dv3->x[1]*dv1->x[2] - dv3->x[2]*dv1->x[1];
  u31.x[1] = dv3->x[2]*dv1->x[0] - dv3->x[0]*dv1->x[2];
  u31.x[2] = dv3->x[0]*dv1->x[1] - dv3->x[1]*dv1->x[0];
  
  /* normalize Uij to unit vectors */
  
  s12 = u12.x[0]*u12.x[0]+u12.x[1]*u12.x[1]+u12.x[2]*u12.x[2];
  s23 = u23.x[0]*u23.x[0]+u23.x[1]*u23.x[1]+u23.x[2]*u23.x[2];
  s31 = u31.x[0]*u31.x[0]+u31.x[1]*u31.x[1]+u31.x[2]*u31.x[2];

  /* test for a degenerate triangle associated with collinear vertices */
  
  if ( !(fabs(s12) > 0.0) || !(fabs(s23) > 0.0) || !(fabs(s31) > 0.0) ) {
    areas = 0.0;
    return areas;
  }

  s12 = sqrt(s12);
  s23 = sqrt(s23);
  s31 = sqrt(s31);
  
  u12.x[0] = u12.x[0]/s12; u12.x[1] = u12.x[1]/s12; u12.x[2] = u12.x[2]/s12;
  u23.x[0] = u23.x[0]/s23; u23.x[1] = u23.x[1]/s23; u23.x[2] = u23.x[2]/s23;
  u31.x[0] = u31.x[0]/s31; u31.x[1] = u31.x[1]/s31; u31.x[2] = u31.x[2]/s31;
  
  /*
   *  Compute interior angles Ai as the dihedral angles between planes:
   *  CA1 = cos(A1) = -<U12,U31>
   *  CA2 = cos(A2) = -<U23,U12>
   *  CA3 = cos(A3) = -<U31,U23>
   */

  ca1 = -( u12.x[0]*u31.x[0]+u12.x[1]*u31.x[1]+u12.x[2]*u31.x[2] );
  ca2 = -( u23.x[0]*u12.x[0]+u23.x[1]*u12.x[1]+u23.x[2]*u12.x[2] );
  ca3 = -( u31.x[0]*u23.x[0]+u31.x[1]*u23.x[1]+u31.x[2]*u23.x[2] );

#if ! defined (FMAX)
#define  FMAX(a,b)  ((a) > (b) ? (a) : (b))
#endif
#if ! defined (FMIN)
#define  FMIN(a,b)  ((a) < (b) ? (a) : (b))
#endif

  ca1 = FMAX(ca1, -1.0);
  ca1 = FMIN(ca1, +1.0);
  ca2 = FMAX(ca2, -1.0);
  ca2 = FMIN(ca2, +1.0);
  ca3 = FMAX(ca3, -1.0);
  ca3 = FMIN(ca3, +1.0);
  
  a1 = acos(ca1);
  a2 = acos(ca2);
  a3 = acos(ca3);
  
  /* compute AREAS = A1 + A2 + A3 - PI */
  
  areas = a1 + a2 + a3 - M_PI;

  if ( areas < 0.0 ) {
    areas = 0.0;
  }

  return areas;
}

static
double cell_area(long i, long nv, double *grid_center_lon, double *grid_center_lat,
		 double *grid_corner_lon, double *grid_corner_lat, int *status)
{
  long k;
  double xa;
  double area;
  struct geo p1, p2, p3;
  struct cart c1, c2, c3;

  area = 0;
      
  p3.lon = grid_center_lon[i]; 
  p3.lat = grid_center_lat[i];
  c3 = gc2cc(&p3);
      
  for ( k = 1; k < nv; ++k )
    {
      p1.lon = grid_corner_lon[i*nv+k-1]; 
      p1.lat = grid_corner_lat[i*nv+k-1];
      c1 = gc2cc(&p1);
      p2.lon = grid_corner_lon[i*nv+k]; 
      p2.lat = grid_corner_lat[i*nv+k];
      c2 = gc2cc(&p2);

      xa = areas(&c1, &c2, &c3);
      area += xa;
    }

  p1.lon = grid_corner_lon[i*nv+0]; 
  p1.lat = grid_corner_lat[i*nv+0];
  c1 = gc2cc(&p1);
  p2.lon = grid_corner_lon[i*nv+nv-1]; 
  p2.lat = grid_corner_lat[i*nv+nv-1];
  c2 = gc2cc(&p2);

  xa = areas(&c1, &c2, &c3);
  area += xa;

  return (area);
}

static
double cellarea(int nv, double grid_center_lon, double grid_center_lat, double *grid_corner_lon, double *grid_corner_lat)
{
  int status, i;
  double area;
  double PlanetRadius = 6371229;
  double center_lon, center_lat;
  double corner_lon[nv], corner_lat[nv]; 

  center_lon = grid_center_lon * deg2rad;
  center_lat = grid_center_lat * deg2rad;

  for ( i = 0; i < nv; ++i )
    {
      corner_lon[i] = grid_corner_lon[i] * deg2rad;
      corner_lat[i] = grid_corner_lat[i] * deg2rad;
    }

  area = cell_area(0, nv, &center_lon, &center_lat, corner_lon, corner_lat, &status);

  return area*PlanetRadius*PlanetRadius*1.e-6;
}

int main(void)
{
  double area;
  double grid_center_lon, grid_center_lat;
  double grid_corner_lon[4], grid_corner_lat[4]; 

  /* set the test data over the Equator

             0.0 (lon)
             0.5 (lat)
             / \
            /   \
  -0.5     /     \  0.5
   0.0     \     /  0.0
            \   /
             \ /
             0.0
            -0.5         
  */

  grid_center_lon    =    0;
  grid_center_lat    =    0;
  grid_corner_lon[0] =  0.5;
  grid_corner_lat[0] =  0.0;
  grid_corner_lon[1] =  0.0;
  grid_corner_lat[1] =  0.5;
  grid_corner_lon[2] = -0.5;
  grid_corner_lat[2] =  0.0;
  grid_corner_lon[3] =  0.0;
  grid_corner_lat[3] = -0.5;

  area = cellarea(4, grid_center_lon, grid_center_lat, grid_corner_lon, grid_corner_lat);

  printf ( "CDO    area of quad     over Equator    is %f sqr km\n", area );

  /* Triangle, half of the above Quad */

  grid_center_lon    =    0;
  grid_center_lat    =    0.25;

  area = cellarea(3, grid_center_lon, grid_center_lat, grid_corner_lon, grid_corner_lat);

  printf ( "CDO    area of triangle over Equator    is %f sqr km\n", area );

  /* set the test data over the Equator

  -0.5               0.5 (lon)
  0.5  -----------  0.5 (lat)
       |         |
       |         |
       |         |
       |         |
       |         |
  -0.5  -----------  0.5
  -0.5              -0.5

  */

  grid_center_lon    =    0;
  grid_center_lat    =    0;
  grid_corner_lon[0] = -0.5;
  grid_corner_lat[0] = -0.5;
  grid_corner_lon[1] =  0.5;
  grid_corner_lat[1] = -0.5;
  grid_corner_lon[2] =  0.5;
  grid_corner_lat[2] =  0.5;
  grid_corner_lon[3] = -0.5;
  grid_corner_lat[3] =  0.5;

  area = cellarea(4, grid_center_lon, grid_center_lat, grid_corner_lon, grid_corner_lat);

  printf ( "CDO    area of quad     over Equator    is %f sqr km\n", area );

  /* Half of the above Quad */

  grid_center_lon    =    0.1;
  grid_center_lat    =   -0.1;

  area = cellarea(3, grid_center_lon, grid_center_lat, grid_corner_lon, grid_corner_lat);

  printf ( "CDO    area of triangle over Equator    is %f sqr km\n", area );

  /* set the test data over the North Pole

  90.0                0.0 (lon)
  89.5  -----------  89.5 (lat)
        |         |
        |         |
        |    x    |
        |         |
        |         |
  180.0  ----------- 270.0
  89.5               89.5

  */

  grid_center_lon    =    0;
  grid_center_lat    =   90;
  grid_corner_lon[0] =   0.0;
  grid_corner_lat[0] =  89.5;
  grid_corner_lon[1] =  90.0;
  grid_corner_lat[1] =  89.5;
  grid_corner_lon[2] =-180.0;
  grid_corner_lat[2] =  89.5;
  grid_corner_lon[3] = -90.0;
  grid_corner_lat[3] =  89.5;

  area = cellarea(4, grid_center_lon, grid_center_lat, grid_corner_lon, grid_corner_lat);

  printf ( "CDO    area of quad     over North Pole is %f sqr km\n", area );

  /* set the test data over the globe

  -180               180 (lon)
    90  -----------   90 (lat)
        |         |
        |         |
        |    x    |
        |         |
        |         |
  -180  ----------- 180
   -90              -90

  */

  grid_center_lon    =    0;
  grid_center_lat    =    0;
  grid_corner_lon[0] =-180.0;
  grid_corner_lat[0] = -90.0;
  grid_corner_lon[1] = 180.0;
  grid_corner_lat[1] = -90.0;
  grid_corner_lon[2] = 180;
  grid_corner_lat[2] =  90.0;
  grid_corner_lon[3] =-180.0;
  grid_corner_lat[3] =  90.0;

  area = cellarea(4, grid_center_lon, grid_center_lat, grid_corner_lon, grid_corner_lat);

  printf ( "CDO    area of quad     over globe is %f sqr km\n", area );

  return 0;
}
