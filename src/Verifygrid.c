/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2016 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
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

#if defined(HAVE_CONFIG_H)
#  include "config.h" /* VERSION */
#endif

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "pstream.h"
#include "clipping/geometry.h"
#include "clipping/clipping.c"
#include "math.h"

static double euclidean_norm (double a[]) {

  /* Computes the Euclidean norm of a vector given in Cartesian coordinates. */

  return sqrt(pow(a[0],2) + pow(a[1],2) + pow(a[2],2));
}

static void divide_by_scalar (double (*a)[3], double scalar) {

  /* Component-wise scalar division of a three dimensional vector given in Cartesian coordinates. */

  (*a)[0] = (*a)[0]/scalar;
  (*a)[1] = (*a)[1]/scalar;
  (*a)[2] = (*a)[2]/scalar;
}

static void normalize_vector(double (*a)[3]){

  /* Normalizes a vector by dividing it though its magnitude. */

  divide_by_scalar(a, euclidean_norm(*a));
}

static int is_point_left_of_edge(double point_on_line_1[], double point_on_line_2[], double point[]){

  /* 
     Computes whether a point is left of the line through point_on_line_1 and point_on_line_2. This is part of the solution to the point in polygon problem.
     Returns 0 if the point is on the line, > 0 if the point is left of the line, and < 0 if the point is right of the line.
     This algorithm is by Dan Sunday (geomalgorithms.com) and is completely free for use and modification.
  */

  return ((point_on_line_2[0] - point_on_line_1[0]) * (point[1] - point_on_line_1[1]) - (point[0] - point_on_line_1[0]) * (point_on_line_2[1] - point_on_line_1[1]));

}

static int winding_numbers_algorithm(double (*cell_corners)[][2], int number_corners, double point[]){
  
  /* 
     Computes whether a point is inside the bounds of a cell. This is the solution to the point in polygon problem.
     Returns 0 if the point is outside, returns 1 if the point is inside the cell.
     Based on an algorithm by Dan Sunday (geomalgorithms.com). His algorithm is completely free for use and modification.
   */
  
  int winding_number;
  
  for (int i = 0;  i < number_corners; i++){
    if ((*cell_corners)[i][1] <= point[1]){
      if ((*cell_corners)[i + 1][1] > point[1])
	if (is_point_left_of_edge((*cell_corners)[i], (*cell_corners)[i + 1], point) > 0)
	  ++winding_number;
    }
    else { if ((*cell_corners)[i + 1][1] <= point[1])
	if (is_point_left_of_edge((*cell_corners)[i], (*cell_corners)[i + 1], point) < 0)
	  --winding_number;
    }
  }

  return winding_number;

}



double intlin(double x, double y1, double x1, double y2, double x2);

static
int pnpoly(int npol, double *xp, double *yp, double x, double y)
{
  int i, j, c = 0;

  for (i = 0, j = npol-1; i < npol; j = i++) {
    if ((((yp[i]<=y) && (y<yp[j])) ||
	 ((yp[j]<=y) && (y<yp[i]))) &&
	(x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]))
      
      c = !c;
  }
  return c;
}


static
double PolygonArea_old(int np, double *xp, double *yp)
{
  int i, j;
  double area = 0;

  for ( i = 0; i < np; i++ )
    {
      j = (i + 1) % np;
      area += xp[i] * yp[j];
      area -= yp[i] * xp[j];
    }

  area /= 2;
  /* return(area < 0 ? -area : area); */
  return (area);
}


static
double PolygonArea(int np, double *xp, double *yp, double yc)
{
  int i, j;
  double area = 0.;

  /* Process area in Radians */
   
  for ( i = 0; i < np; i++ )
    {
      j = (i + 1) % np;
      area += DEG2RAD*xp[i] * DEG2RAD*yp[j];
      area -= DEG2RAD*yp[i] * DEG2RAD*xp[j];
    }
  area *= 0.5 * cos(DEG2RAD*yc);
  return (area);
}

static
int ccw(double p0x, double p0y, double p1x, double p1y, double p2x, double p2y)
{
  /*
    This function says wether the point are orientated clockwise
    +1 positive orientation
    -1 negative orientation
     0 points are on a line --> no orientation
    
    This is done by a comparision of the gradient of
    dy1/dx1 = p1 - p0 vs.
    dy2/dx2 = p2 - p0
    To avoid singularities at dx1=0 OR dx2 = 0 we multiply with dx1*dx2
  */
  double dx1, dx2, dy1, dy2;

  dx1 = p1x - p0x; dy1 = p1y - p0y;
  dx2 = p2x - p0x; dy2 = p2y - p0y;
  if ( dx1*dy2 > dy1*dx2 ) return +1;
  if ( dx1*dy2 < dy1*dx2 ) return -1;
  if ( (dx1*dx2 < 0 ) || (dy1*dy2 < 0)) return -1;
  if ( (dx1*dx1 + dy1*dy1) < (dx2*dx2 + dy2*dy2)) return +1;

  return 0;
}

static
int intersect(double pix, double piy, double pjx, double pjy,
              double pkx, double pky, double plx, double ply)
{
  /*This function returns if there is an intersection between the lines 
    line1 between pi and pj and
    line2 between pk and pl,
    whereas pi = (pix, piy).
      
    This can done by means of ccw since the product of ccw(pi,pj,pk)*ccw(pi,pj,pl)
    shows if pk and pl are on different or the same side(s) of the line1 (They must
    have different signums to be on different sides).
      
    Consequently if and ONLY IF pk as well as pl are on different sides of line1
    AND pi as well as pj are on different sides of line2 there HAS TO be an intersection.
  */
    
  return ( ( ccw(pix, piy, pjx, pjy, pkx, pky) *
	     ccw(pix, piy, pjx, pjy, plx, ply) <= 0 ) &&
	   ( ccw(pkx, pky, plx, ply, pix, piy) *
	     ccw(pkx, pky, plx, ply, pjx, pjy) <= 0 ) );
}

static
int check_ncorner(int ncorner, const double *lon_bounds, const double *lat_bounds)
{
  int ncorner_new = ncorner;
  int k;

  for ( k=ncorner-1; k>0; --k )
    if ( IS_NOT_EQUAL(lon_bounds[k], lon_bounds[k-1]) ||
	 IS_NOT_EQUAL(lat_bounds[k], lat_bounds[k-1]) ) break;

  if ( k < ncorner-1 ) ncorner_new = k+1;

  return ncorner_new;
}

static
void verify_grid(int gridsize, int ncorner,
		double *grid_center_lon, double *grid_center_lat,
		double *grid_corner_lon, double *grid_corner_lat)
{
  int i0, i, j, k, l;
  int l0;
  int nout;
  int isinside, convex, alone, isnegative;
  const int mnv = ncorner+1;
  int cuts[mnv][mnv];  
  int *alone_cell;          
  int check_corners;
  double lon, lat = 0;
  double lon_bounds[mnv], lat_bounds[mnv];
  double area, sumarea;

  alone_cell = (int*) Malloc(gridsize*ncorner*sizeof(int));

  check_corners = 0; /* don't execute corner checking (last loop) */
  nout = 0;
  sumarea = 0;
  /*
  for ( i = 0; i < gridsize; ++i )
    {
      lon = grid_center_lon[i];
      lat = grid_center_lat[i];
      for ( k = 0; k < ncorner; ++k )
        {
          lon_bounds[k] = grid_corner_lon[i*ncorner+k];
          lat_bounds[k] = grid_corner_lat[i*ncorner+k];
          if ( (lon - lon_bounds[k]) > 270 ) lon_bounds[k] += 360;
          if ( (lon_bounds[k] - lon) > 270 ) lon_bounds[k] -= 360;
        }      
      lon_bounds[ncorner] = lon_bounds[0];
      lat_bounds[ncorner] = lat_bounds[0];
      fprintf(stdout, " %6i %6i %9.4f %9.4f :",  nout, i+1, lon, lat);
      for ( k = 0; k < ncorner; k++ )
	fprintf(stdout, " %9.4f %9.4f : ", lon_bounds[k], lat_bounds[k]);
      fprintf(stdout, "\n");
    }
  */

  /* Check if center is inside bounds of cell */
  for ( i = 0; i < gridsize; ++i )
    {
      lon = grid_center_lon[i];
      lat = grid_center_lat[i];

      for ( k = 0; k < ncorner; ++k )
        {
          lon_bounds[k] = grid_corner_lon[i*ncorner+k];
          lat_bounds[k] = grid_corner_lat[i*ncorner+k];
          if ( (lon - lon_bounds[k]) > 270 ) lon_bounds[k] += 360;
          if ( (lon_bounds[k] - lon) > 270 ) lon_bounds[k] -= 360;
        }      
      lon_bounds[ncorner] = lon_bounds[0];
      lat_bounds[ncorner] = lat_bounds[0];
      
      isinside = pnpoly(ncorner+1, lon_bounds, lat_bounds, lon, lat);

      if ( !isinside ) nout++;
      if ( !isinside && cdoVerbose )
        {
          if ( nout == 1 )
            {
              fprintf(stdout,"\n CENTER IS OUT OF BOUNDS");
              fprintf(stdout,"\n                                               :");
              for ( k = 0; k < ncorner; k++ )
                fprintf(stdout, "          Corner %2i : ", k+1);
              fprintf(stdout,"\n Number  Index center_lon center_lat area*10^6 :");
              for ( k = 0; k < ncorner; k++ )
                fprintf(stdout, "   lon_%2.2i    lat_%2.2i : ", k+1, k+1);
              fprintf(stdout, "\n");
            }
          area = PolygonArea(ncorner+1, lon_bounds, lat_bounds,lat);
          fprintf(stdout, " %6i %6i  %9.4f  %9.4f %9.5f :", 
		  nout, i+1, lon, lat, area*pow(10,6));

	  int ncorner_new = check_ncorner(ncorner, lon_bounds, lat_bounds);

          for ( k = 0; k < ncorner_new; k++ )
	    fprintf(stdout, "%9.4f %9.4f : ", lon_bounds[k], lat_bounds[k]);
           for ( k = ncorner_new; k < ncorner; k++ )
	     fprintf(stdout, "     ----      ---- : ");
          fprintf(stdout, "\n");
        }
    }

  if ( nout )
    cdoWarning("%d of %d points out of bounds!", nout, gridsize);
  
  /* check that all cell bounds have the same orientation */
  
  nout = 0;
  for ( i = 0; i < gridsize; ++i )
    {
      lon = grid_center_lon[i];
      lat = grid_center_lat[i];
      
      for ( k = 0; k < ncorner; ++k )
	{
          lon_bounds[k] = grid_corner_lon[i*ncorner+k];
          lat_bounds[k] = grid_corner_lat[i*ncorner+k];
          if ( (grid_center_lon[i] - lon_bounds[k]) > 270 ) lon_bounds[k] += 360;
          if ( (lon_bounds[k] - grid_center_lon[i]) > 270 ) lon_bounds[k] -= 360;
	}
      lon_bounds[ncorner] = lon_bounds[0];
      lat_bounds[ncorner] = lat_bounds[0];
      
      area = PolygonArea(ncorner+1, lon_bounds, lat_bounds, lat);
      
      isnegative = area < 0 ? 1 : 0;
      sumarea += area < 0 ? -area : area;
      
      if ( isnegative ) nout++;
      
      if ( isnegative && cdoVerbose )
        {
          if ( nout == 1 )
            {
              fprintf(stdout,"\n                                     :");
              for ( k = 0; k < ncorner; k++ )
                fprintf(stdout, "          Corner %2i : ", k+1);
              fprintf(stdout,"\n Number  Index center_lon center_lat :");
              for ( k = 0; k < ncorner; k++ )
                fprintf(stdout, "   lon_%2.2i    lat_%2.2i : ", k+1, k+1);
              fprintf(stdout, "\n");
            }
          fprintf(stdout, " %6i %6i  %9.4f  %9.4f :", nout, i+1, lon, lat);

	  int ncorner_new = check_ncorner(ncorner, lon_bounds, lat_bounds);

          for ( k = 0; k < ncorner_new; k++ )
	    fprintf(stdout, "%9.4f %9.4f : ", lon_bounds[k], lat_bounds[k]);
           for ( k = ncorner_new; k < ncorner; k++ )
	     fprintf(stdout, "     ----      ---- : ");

          fprintf(stdout, "\n");
        }
    }

  if ( nout )
    cdoWarning("%d of %d grid cells have wrong orientation!", nout, gridsize);

  if ( cdoVerbose ) 
    fprintf(stdout, "area-error: %9.5f%%\n", 100.*(sumarea - 4.*M_PI)/4.*M_PI );

  if ( fabs(100.*(sumarea - 4.*M_PI)/4.*M_PI) > 0.1)
    cdoWarning("area-error: %9.5f%%", 100.*(sumarea - 4.*M_PI)/4.*M_PI );
  
  /* check that all cells are convex */
  
  nout = 0;
  for ( i0 = 0; i0 < gridsize; i0++ )
    {
      lon = grid_center_lon[i0];
      lat = grid_center_lat[i0];

      for ( k = 0; k < ncorner; k++ )
	{
	  lon_bounds[k] = grid_corner_lon[i0*ncorner+k];
	  lat_bounds[k] = grid_corner_lat[i0*ncorner+k];
	  /* Find cells that cover left and right border of the grid and adjust
	     coordinates --> they become closed polygons on theta-phi plane! */
	  if ( (lon - lon_bounds[k]) > 270 ) lon_bounds[k] += 360; 
	  if ( (lon_bounds[k] - lon) > 270 ) lon_bounds[k] -= 360;
	}
      
      /* Reset found cuts for the current cell before starting the search */
      for ( i = 0; i < ncorner; i++ )
	for ( j = 0; j < ncorner; j++ )
	  cuts[i][j] = 0;
      
      /* Loops cover all combinations between inner lines of the Polygon
	 Check whether each inner line is cut by an other (inner) one at least once. 
	 - Only if there is a cut every inner line the Polygon is convex
	 - We assume: Points are in either cyclic or anticyclic order
      */
      for ( i = 0; i < ncorner-1; i++ )
	{
          /* j = i+2 excludes lines from one corner to an other (j=i+1) and
	     from one point to itself (j=i)*/
          for ( j = i+2 ; j < ncorner; j++ )
	    {
              /* Exclude the line between the last and first corner */
              if ( i == 0 && j == ncorner-1 ) continue;

	      /* k = i+1: if starting point is in common lines to different corners
		 do not intersect */
              for ( k = i+1; k < ncorner - 1; k++ )
		{                  
                  if ( i == k ) l0 = j+1;
                  else          l0 = k+2;

                  for ( l = l0; l < ncorner; l++ )
		    {
                      if ( cuts[k][l] && cuts[i][j] ) continue;
		      /* Exlude the line between the last and first corner 
			 Exlude the line itself (l!=i, k!=j)
			 Check if line ij and kl intersect each other.
			 If so increment respective counters for intersections. 
			 It is not relevant by which line a line is intersected - 
			 it is only relevant if they is itersected! */
                      if ( ! ( k==0 && l == ncorner-1 ) && ( l != j ) && ( k != j )  )
			{
                          if ( intersect(lon_bounds[i], lat_bounds[i], lon_bounds[j], lat_bounds[j],
                                         lon_bounds[k], lat_bounds[k], lon_bounds[l], lat_bounds[l]) )
			    {
			      cuts[i][j]++; cuts[k][l]++; cuts[j][i]++; cuts[l][k]++;
			    }
			}
		    }
		}                  
	    }
	}

      convex = 1;
      /* The following loop covers all inner lines of the Polygon 
	 (The assumption applies that the points are in cyclic order) */
      for ( i = 0; i < ncorner-1; i++ )
	for ( j = i+2; j < ncorner; j++)
	  {
	    if ( i == 0 && j == ncorner-1 ) continue;	   
	    if ( ! cuts[i][j] ) convex = 0;
	  }
      if ( !convex ) nout++;        
      if ( cdoVerbose && ( !convex ) )
	{
          if ( nout == 1 )
	    {
              fprintf(stdout,"\n NO CONVEX POLYGON");
              fprintf(stdout,"\n                                       :");
              for ( k = 0; k < ncorner; k++ )
		fprintf(stdout, "            Corner %2i : ", k);
              fprintf(stdout,"\n Number  Index  center_lon  center_lat :");
              for ( k = 0; k < ncorner; k++ )
		fprintf(stdout, "    lon_%2.2i     lat_%2.2i : ", k, k);
              fprintf(stdout, "\n");
	    }
          
          fprintf(stdout, " %6i %6i   %9.4f   %9.4f :", nout, i0+1, lon, lat);
          for ( k = 0; k < ncorner; k++ )
	    fprintf(stdout, "  %9.4f %9.4f : ", lon_bounds[k], lat_bounds[k]);
          fprintf(stdout, "\n");         
	}     
    }

  if ( nout )
    cdoWarning("%d of %d cells are not Convex!", nout, gridsize);

  if ( check_corners )
    {
      /* 
	 Check if there is a corner at the same point of 
	 an other cell foreach corner of each cell 
      */
      nout = 0;
      for ( i = 0; i < gridsize*ncorner; i++ )
	alone_cell[i] = 1;
      
      for ( i = 0; i < gridsize*ncorner; i++ )
	{
	  if ( ! alone_cell[i] ) continue;
	  alone = 1;
	  lon = grid_corner_lon[i];
	  lat = grid_corner_lat[i];			
	  for ( j = 0; j < gridsize*ncorner; j++ )
	    if ( j != i && 
		 IS_EQUAL(grid_corner_lat[j], lat) && 
		 IS_EQUAL(grid_corner_lon[j], lon) )
	      { alone = 0; alone_cell[i] = alone_cell[j] = 1; break; }
	  if ( alone )
	    {
	      if      ( lon >= 180. ) lon -= 360.;
	      else if ( lon  < 180. ) lon += 360.;
	      for ( j = i+1; j < gridsize*ncorner; j++ )
		if (j != i  && 
		    IS_EQUAL(grid_corner_lat[j], lat) && 
		    IS_EQUAL(grid_corner_lon[j], lon) )
		  { alone = 0; alone_cell[i] = alone_cell[j] = 0; break; }
	    }
	  if ( alone )
	    { 
	      nout++;
	      if ( cdoVerbose )
		{
		  if ( nout == 1 )
		    {
		      fprintf(stdout,"\n VERTEX ALONE ON GRID\n");
		      fprintf(stdout," number cell-Index  Vert-Index :        lon        lat\n");
		    }							
		  fprintf(stdout, " %6i     %6i      %6i : %10.4f %10.4f\n", 
			  nout, i/ncorner, i, grid_corner_lon[i], grid_corner_lat[i]);
		}					
	    }
	}

      if ( nout )
	cdoWarning("%d of %d corners are lonely on the grid!", nout, gridsize*ncorner);
    }

  Free(alone_cell);
}


void verify_grid_old(int gridsize, int ncorner,
		double *grid_center_lon, double *grid_center_lat,
		double *grid_corner_lon, double *grid_corner_lat)
{
  int i, k;
  int nout;
  int isinside;
  int isnegative;
  double area;
  double lon, lat;
  double lon_bounds[ncorner], lat_bounds[ncorner];

  /* check that all centers are inside the bounds */

  nout = 0;
  for ( i = 0; i < gridsize; ++i )
    {
      lon = grid_center_lon[i];
      lat = grid_center_lat[i];

      for ( k = 0; k < ncorner; ++k )
	{
	  lon_bounds[k] = grid_corner_lon[i*ncorner+k];
	  lat_bounds[k] = grid_corner_lat[i*ncorner+k];
	}

      for ( k = 0; k < ncorner; ++k )
	{
	  if ( (lon - lon_bounds[k]) > 270 ) lon_bounds[k] += 360;
	  if ( (lon_bounds[k] - lon) > 270 ) lon_bounds[k] -= 360;
	}

      lon_bounds[ncorner] = lon_bounds[0];
      lat_bounds[ncorner] = lat_bounds[0];

      isinside = pnpoly(ncorner+1, lon_bounds, lat_bounds, lon, lat);

      if ( !isinside ) nout++;

      if ( !isinside && cdoVerbose )
	printf("center: %d %d %g %g %g %g %g %g %g %g %g %g\n", nout, i, lon, lat, lon_bounds[0], lat_bounds[0],
	       lon_bounds[1], lat_bounds[1], lon_bounds[2], lat_bounds[2], lon_bounds[3], lat_bounds[3]);
    }

  if ( nout > 0 )
    cdoWarning("%d of %d points out of bounds!", nout, gridsize);


  /* check that all cell bounds have the same orientation */

  nout = 0;
  for ( i = 0; i < gridsize; ++i )
    {
      lon = grid_center_lon[i];
      lat = grid_center_lat[i];

      for ( k = 0; k < ncorner; ++k )
	{
	  lon_bounds[k] = grid_corner_lon[i*ncorner+k];
	  lat_bounds[k] = grid_corner_lat[i*ncorner+k];
	}

      for ( k = 0; k < ncorner; ++k )
	{
	  if ( (grid_center_lon[i] - lon_bounds[k]) > 270 ) lon_bounds[k] += 360;
	  if ( (lon_bounds[k] - grid_center_lon[i]) > 270 ) lon_bounds[k] -= 360;
	}

      lon_bounds[ncorner] = lon_bounds[0];
      lat_bounds[ncorner] = lat_bounds[0];

      area = PolygonArea_old(ncorner+1, lon_bounds, lat_bounds);

      if ( area < 0 ) isnegative = 1;
      else            isnegative = 0;

      if ( isnegative ) nout++;


      if ( isnegative && cdoVerbose )
	printf("bounds: %d %d %g %g %g %g %g %g %g %g %g %g\n", nout, i, lon, lat, lon_bounds[0], lat_bounds[0],
	       lon_bounds[1], lat_bounds[1], lon_bounds[2], lat_bounds[2], lon_bounds[3], lat_bounds[3]);
    }

  if ( nout > 0 )
    cdoWarning("%d of %d grid cells have wrong orientation!", nout, gridsize);
}


static void verify_grid_test(int gridsize, int ncorner, double *grid_center_lon, double *grid_center_lat, double *grid_corner_lon, double *grid_corner_lat){
  
  /* 
     This function performs three tests on each cell of a given grid:

     1) it tests whether the center point is within the bounds of the cell
     2) it tests whether the corners all cell bounds have the same orientation, i.e. the corners of the cell are in clockwise or counterclockwise order
     3) it tests whether the cell is convex

     It performs these tests after longitude and latitude on the unit circle have been converted first to Cartesian coordinates in Euclidean space and subsequently to two dimensional coordinates on the plane each cell occupies.  
  */
  
  double center_point_in_Euclidean_space[3];
  double cell_corners_in_Euclidean_space[ncorner][3];
  double center_point_on_cell_plane[2];
  double cell_corners_on_cell_plane[ncorner+1][2];

  int cell_no;
  int corner_no;

  for (cell_no = 0; cell_no < gridsize; ++cell_no)
    {

       printf("Cell number %d:\n\n", cell_no+1);

       /* Latitude and longitude are spherical coordinates on a unit circle. Each such coordinate tuple is transformed into a triple of cartesian coordinates in Euclidean space. LLtoXYZ is located in clipping/geometry.h */

       /* This is first done for the presumed center point of the cell. LLtoXYZ is defined in clipping/geometry.h */

      LLtoXYZ(grid_center_lon[cell_no], grid_center_lat[cell_no], center_point_in_Euclidean_space);

      /* Then the transformation is done for all the corners of the cell. */

      for (corner_no = 0; corner_no < ncorner; ++corner_no)
	{
	  LLtoXYZ(grid_corner_lon[cell_no*ncorner+corner_no], grid_corner_lat[cell_no*ncorner+corner_no], cell_corners_in_Euclidean_space[corner_no]);
	}

      /*
	Each cell corresponds to a two-dimensional polygon now in unknown orientation in three-dimensional space. Each cell and its center point are coplanar. THIS IS A GIVEN.
	In order to solve the two-dimensional point-in-polygon problem for each cell, the three-dimensional coordinates of the polygon and its center point are projected onto the two-dimensional plane they form.
        
	This is done in the following steps:
      
	1) Compute the normal of the plane using three corners. The normal is the new z-axis. AT LEAST THREE CORNERS ARE A GIVEN.
	2) Compute the new y-axis by computing the cross product of the new z-axis and the old x-axis.
	3) Compute the new x-axis by computing the cross product of the new z-axis and the new y-axis.
	4) Divide all three axis vectors by their length cutting them to unit length.
	5) Project every corner point onto the new x- and y-axes by using the dot product.

	The result is a xy tuple for each corner and the presumend center point  which is a projection onto the plane the xyz corner points form.
      */

      double new_x_axis[3];
      double new_y_axis[3];
      double new_z_axis[3];
      double old_x_axis[3] = {1, 0, 0};
      double old_y_axis[3] = {0, 1, 0};
      double old_z_axis[3] = {0, 0, 1};

      /* The surface normal is the result of the cross product of two edges of a cell polygon. The edges must not be parallel, i.e their cross product must not be zero. */

      /* The two edges A and B are calculated from the first three corners of the cell. The first edge A is corner no. 2 - corner no. 1. */

      double A[3] = {cell_corners_in_Euclidean_space[1][0] - cell_corners_in_Euclidean_space[0][0], cell_corners_in_Euclidean_space[1][1] - cell_corners_in_Euclidean_space[0][1], cell_corners_in_Euclidean_space[1][2] - cell_corners_in_Euclidean_space[0][2]};
      
      /* The second edge B is corner no. 3 - corner no. 1. */

      double B[3] = {cell_corners_in_Euclidean_space[2][0] - cell_corners_in_Euclidean_space[0][0], cell_corners_in_Euclidean_space[2][1] - cell_corners_in_Euclidean_space[0][1], cell_corners_in_Euclidean_space[2][2] - cell_corners_in_Euclidean_space[0][2]};

      /* The cross product of the two edges A and B is the surface normal and the new z-axis. crossproduct_d is defined in clipping/geometry.h */

      crossproduct_d(A, B, new_z_axis);

      printf("The cross product of vector A (%f, %f, %f) and vector B (%f, %f, %f) is the normal vector (%f, %f, %f), the new z-axis.\n", A[0], A[1], A[2], B[0], B[1], B[2], new_z_axis[0], new_z_axis[1], new_z_axis[2] );

      /* Then the new y-axis is the result of the cross product of the new z-axis and the old x-axis, and the new x-axis is the result of the cross product of the new z-axis and the new y-axis. */

      crossproduct_d(old_x_axis, new_z_axis, new_y_axis);

      printf("The cross product of the old x-axis (%f, %f, %f) and the new z-axis (%f, %f, %f) is the new y-axis (%f, %f, %f).\n", old_x_axis[0], old_x_axis[1], old_x_axis[2], new_z_axis[0], new_z_axis[1], new_z_axis[2], new_y_axis[0], new_y_axis[1], new_y_axis[2]);

      crossproduct_d(new_z_axis, new_y_axis, new_x_axis);

      printf("The cross product of the new z-axis (%f, %f, %f) and the new y-axis (%f, %f, %f) is the new x-axis (%f, %f, %f).\n\n", new_z_axis[0], new_z_axis[1], new_z_axis[2], new_y_axis[0], new_y_axis[1], new_y_axis[2], new_x_axis[0], new_x_axis[1], new_x_axis[2]);
      	
      /* The axis vectors need to be normalized in order to function as basis for a coordinate system. */
      
      printf("The z-axis coordinates before normalization are: (%f, %f, %f)\n", new_z_axis[0], new_z_axis[1], new_z_axis[2]);
      printf("The y-axis coordinates before normalization are: (%f, %f, %f)\n", new_y_axis[0], new_y_axis[1], new_y_axis[2]);
      printf("The x-axis coordinates before normalization are: (%f, %f, %f)\n\n", new_x_axis[0], new_x_axis[1], new_x_axis[2]);
        
      normalize_vector(&new_x_axis);
      normalize_vector(&new_y_axis);
      normalize_vector(&new_z_axis);

      printf("The z-axis coordinates after normalization are: (%f, %f, %f)\n", new_z_axis[0], new_z_axis[1], new_z_axis[2]);
      printf("The y-axis coordinates after normalization are: (%f, %f, %f)\n", new_y_axis[0], new_y_axis[1], new_y_axis[2]);
      printf("The x-axis coordinates after normalization are: (%f, %f, %f)\n\n", new_x_axis[0], new_x_axis[1], new_x_axis[2]);
      
      printf("Checking orthogonality of coordinate axes after normalization...\n");
      printf("Scalar product of x- and z-axes: %f\n", dotproduct(new_x_axis, new_z_axis));
      printf("Scalar product of y- and z-axes: %f\n", dotproduct(new_y_axis, new_z_axis));
      printf("Scalar product of x- and y-axes: %f\n\n", dotproduct(new_x_axis, new_y_axis));


      /* All corner points are projected onto the new x- and y-axes. dotproduct is defined in clipping/clipping.c */

      for (corner_no = 0; corner_no < ncorner; ++corner_no)
	{
	  cell_corners_on_cell_plane[corner_no][0] = dotproduct(new_x_axis, cell_corners_in_Euclidean_space[corner_no]);
	  cell_corners_on_cell_plane[corner_no][1] = dotproduct(new_y_axis, cell_corners_in_Euclidean_space[corner_no]);  
	}
      
      cell_corners_on_cell_plane[ncorner][0] = cell_corners_on_cell_plane[0][0];
      cell_corners_on_cell_plane[ncorner][1] = cell_corners_on_cell_plane[0][1];
       


      /* As well as the center point. */

      center_point_on_cell_plane[0] = dotproduct(new_x_axis, center_point_in_Euclidean_space);
      center_point_on_cell_plane[1] = dotproduct(new_y_axis, center_point_in_Euclidean_space);
      
      printf("The center xyz coordinates in respect to the original reference frame are: %f, %f, %f\n", center_point_in_Euclidean_space[0], center_point_in_Euclidean_space[1], center_point_in_Euclidean_space[2]);
      printf("The center xy coordinates in respect to the reference frame contructed from the surface normal of the cell as its z-axis are: %f, %f\n\n", center_point_on_cell_plane[0], center_point_on_cell_plane[1]);
     
      /* The winding numbers algorithm is used to test whether the presumed center point is within the bounds of the cell. */
        
      int winding_number = winding_numbers_algorithm(&cell_corners_on_cell_plane, ncorner, center_point_on_cell_plane);

      printf("If the winding number is not 0 the presumed center point is within the bounds of the cell. The winding number is: %u\n\n", winding_number);

}

 
  cdoAbort("implementation missing");
}

void *Verifygrid(void *argument)
{
  bool lgrid_gen_bounds = false, luse_grid_corner = true;
  double *grid_corner_lat = NULL, *grid_corner_lon = NULL;
  char units[CDI_MAX_NAME];

  cdoInitialize(argument);

  int VERIFYGRID     = cdoOperatorAdd("verifygrid",  0,   0, NULL);
  int VERIFYGRIDTEST = cdoOperatorAdd("verifygridtest",  0,   0, NULL);

  int operatorID = cdoOperatorID();

  int streamID = streamOpenRead(cdoStreamName(0));

  int vlistID = streamInqVlist(streamID);

  int gridID  = vlistInqVarGrid(vlistID, 0);

  if ( gridInqType(gridID) == GRID_GME ) gridID = gridToUnstructured(gridID, 1);

  if ( gridInqType(gridID) != GRID_UNSTRUCTURED && gridInqType(gridID) != GRID_CURVILINEAR )
    {
      gridID = gridToCurvilinear(gridID, 1);
      lgrid_gen_bounds = TRUE;
    }

  int gridsize = gridInqSize(gridID);
  /*
  if ( gridInqMaskGME(gridID, NULL) )
    {
      int *grid_mask = (int*) Malloc(gridsize*sizeof(int));
      gridInqMaskGME(gridID, grid_mask);
      free(grid_mask);
    }
  */
  int ncorner = 4;
  if ( gridInqType(gridID) == GRID_UNSTRUCTURED )
    ncorner = gridInqNvertex(gridID);

  double *grid_center_lat = (double*) Malloc(gridsize*sizeof(double));
  double *grid_center_lon = (double*) Malloc(gridsize*sizeof(double));

  gridInqYvals(gridID, grid_center_lat);
  gridInqXvals(gridID, grid_center_lon);

  /* Convert lat/lon units if required */
  gridInqXunits(gridID, units);
  grid_to_degree(units, gridsize, grid_center_lon, "grid center lon");
  gridInqYunits(gridID, units);
  grid_to_degree(units, gridsize, grid_center_lat, "grid center lat");

  if ( luse_grid_corner )
    {
      if ( ncorner == 0 ) cdoAbort("grid corner missing!");
      int nalloc = ncorner*gridsize;
      grid_corner_lat = (double*) Realloc(grid_corner_lat, nalloc*sizeof(double));
      grid_corner_lon = (double*) Realloc(grid_corner_lon, nalloc*sizeof(double));

      if ( gridInqYbounds(gridID, NULL) && gridInqXbounds(gridID, NULL) )
	{
	  gridInqYbounds(gridID, grid_corner_lat);
	  gridInqXbounds(gridID, grid_corner_lon);
	}
      else
	{
	  if ( lgrid_gen_bounds )
	    {
	      char xunitstr[CDI_MAX_NAME];
	      char yunitstr[CDI_MAX_NAME];
	      gridInqXunits(gridID, xunitstr);
	      gridInqYunits(gridID, yunitstr);
	    }
	  else
	    cdoAbort("Grid corner missing!");
	}


      /* Note: using units from latitude instead from bounds */
      grid_to_degree(units, ncorner*gridsize, grid_corner_lon, "grid corner lon");
      grid_to_degree(units, ncorner*gridsize, grid_corner_lat, "grid corner lat");

    }

  streamClose(streamID);

  if ( operatorID == VERIFYGRID )
    verify_grid(gridsize, ncorner, grid_center_lon, grid_center_lat, grid_corner_lon, grid_corner_lat);
  else
    verify_grid_test(gridsize, ncorner, grid_center_lon, grid_center_lat, grid_corner_lon, grid_corner_lat);

  if ( grid_center_lon ) Free(grid_center_lon);
  if ( grid_center_lat ) Free(grid_center_lat);
  if ( grid_corner_lon ) Free(grid_corner_lon);
  if ( grid_corner_lat ) Free(grid_corner_lat);

  cdoFinish();

  return 0;
}
