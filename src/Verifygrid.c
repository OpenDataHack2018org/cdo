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
#include "stdio.h"


static void get_next_point_and_calculate_delta(double *p_delta, double *p_previous_point, double *p_current_point, double cell_corners[], int *p_corner_to_read){
  
  fprintf(stderr,"The previous_point is: (%f, %f)\n", p_previous_point[0], p_previous_point[1]);
  
  p_current_point[0] = cell_corners[*p_corner_to_read * 2 + 0];		
  p_current_point[1] = cell_corners[*p_corner_to_read * 2 + 1];

  fprintf(stderr,"The current_point is: (%f, %f)\n", p_current_point[0], p_current_point[1]);

  *p_corner_to_read = *p_corner_to_read + 1;
  p_delta[0] = p_current_point[0] - p_previous_point[0];			
  p_delta[1] = p_current_point[1] - p_previous_point[1];

  fprintf(stderr,"Resulting delta is: (%f, %f)\n\n", p_delta[0], p_delta[1]);
			 
}		

static int calculate_direction(double *p_delta){
 
  if(p_delta[0] > 0){
    return -1;}
  
  if(p_delta[0] < 0){
    return 1;}
   
  if(p_delta[1] > 0){ 
    return -1;}
  
  if(p_delta[1] < 0){
    return 1;}
  
  return 0;
}

static double x_calculate_direction(double *p, double *q){
  if (p[0] < q[0]) return -1;    
  if (p[0] > q[0]) return  1;	
  if (p[1] < q[1]) return -1;	
  if (p[1] > q[1]) return  1;  
  return 0;	
}

static double cross(double *p, double *q){ 
  return p[0] * q[1] - p[1] * q[0];
}

static int x_ConvexCheckTriple(int *p_current_direction, int *p_this_direction, int *p_direction_changes, double *p_current_delta, double *p_previous_delta, int *p_angle_sign, double *p_second_corner, double *p_third_corner){

  if( ((*p_this_direction) = calculate_direction(p_current_delta)) == -(*p_current_direction)){ 
    
    fprintf(stderr,"A direction change occured!\n");
    fprintf(stderr,"No of direction changes before this one: %d\n", *p_direction_changes);
    
    (*p_direction_changes) = (*p_direction_changes) + 1;						
    
    if(*p_direction_changes > 2){
      fprintf(stderr,"FIRST RETURN!\n\n");
      return 1;
    }				
  } else { 
    fprintf(stderr,"No direction change occured!\n"); 
    fprintf(stderr,"No of direction changes so far: %d\n", *p_direction_changes);
  }									
  
  fprintf(stderr,"The old and new direction are: %d, %d\n", *p_current_direction, *p_this_direction);

  *p_current_direction = *p_this_direction;

  fprintf(stderr,"The new 'old' direction is: %d\n\n", *p_current_direction);

  double cross_product = cross(p_previous_delta, p_current_delta);		

  if(cross_product > 0){ 
    if (*p_angle_sign == -1){ 
      fprintf(stderr,"SECOND RETURN!\n\n");
      return 1;
    }	
    *p_angle_sign = 1;							
  }									
  else if (cross_product < 0){ 
    if (*p_angle_sign == 1){ 
      fprintf(stderr,"THIRD RETURN!\n\n");
      return 1;
    }	
    *p_angle_sign = -1;							
  }									
  p_second_corner[0] = p_third_corner[0];
  p_second_corner[1] = p_third_corner[1];
  p_previous_delta[0] = p_current_delta[0];		
  p_previous_delta[1] = p_current_delta[1];	
}


static int x_classify_polygon(int number_corners, double cell_corners[]){
 
  int current_direction= 0, this_direction = 0, direction_changes = 0, angle_sign = 0, corner_to_read = 1;
  int *p_current_direction, *p_this_direction, *p_direction_changes, *p_angle_sign, *p_corner_to_read;

  p_current_direction = &current_direction;
  p_this_direction = &this_direction;
  p_direction_changes = &direction_changes;
  p_angle_sign = &angle_sign;
  p_corner_to_read = &corner_to_read;
  
  double second_corner[2], third_corner[2], first_corner[2], save_second_corner[2], previous_delta[2], current_delta[2];
  double *p_second_corner, *p_third_corner, *p_first_corner, *p_save_second_corner, *p_previous_delta, *p_current_delta;

  p_second_corner = &second_corner[0];
  p_third_corner = &third_corner[0];
  p_first_corner = &first_corner[0];
  p_save_second_corner = &save_second_corner[0];
  p_previous_delta = &previous_delta[0];
  p_current_delta = &current_delta[0];
    
  if (number_corners < 3){ 
    return 3;
  }
  
  first_corner[0] = cell_corners[0];
  first_corner[1] = cell_corners[1];

  fprintf(stderr,"The first corner set in x_classify_polygon is: (%f, %f)\n\n", p_first_corner[0], p_first_corner[1]);

  while(1){
    
    fprintf(stderr,"The (first) call of x_ConvexGetPointDelta\n\n");

    get_next_point_and_calculate_delta(p_previous_delta, p_first_corner, p_second_corner, cell_corners, p_corner_to_read);

    if (previous_delta[0] || previous_delta[1]){
      break;
    }

    /* Check if there are any corners left. */
    
    if (corner_to_read >= number_corners){
      return 3;
    }
  }

  p_save_second_corner = p_second_corner;
  
  /* This finds the initial direction. */

  current_direction = calculate_direction(p_previous_delta);	

  fprintf(stderr,"The initial direction is: %d\n\n", current_direction);

  while (corner_to_read < number_corners) {
    
    /* Get different point, break if no more points */
    
    fprintf(stderr,"A call of x_ConvexGetPointDelta\n\n");
    
    get_next_point_and_calculate_delta(p_current_delta, p_second_corner, p_third_corner, cell_corners, p_corner_to_read);

    if ( current_delta[0] == 0.0  && current_delta[1] == 0.0 ){
      continue;
    }
    
    fprintf(stderr,"Call of x_ConvexCheckTriple.\n\n");

    x_ConvexCheckTriple(p_current_direction, p_this_direction, p_direction_changes, p_current_delta, p_previous_delta, p_angle_sign, p_second_corner, p_third_corner);
  }

  /* Must check for direction changes from last vertex back to first */
  third_corner[0] = cell_corners[0];
  third_corner[1] = cell_corners[1];
  
  /* Prepare for 'ConvexCheckTriple' */
  current_delta[0] = third_corner[0] - second_corner[0];
  current_delta[1] = third_corner[1] - second_corner[1];
  if (calculate_direction(p_current_delta)){
    x_ConvexCheckTriple(p_current_direction, p_this_direction, p_direction_changes, p_current_delta, p_previous_delta, p_angle_sign, p_second_corner, p_third_corner);
  }
  
  /* and check for direction changes back to second vertex */
  current_delta[0] = save_second_corner[0] - second_corner[0];
  current_delta[1] = save_second_corner[1] - second_corner[1];
  x_ConvexCheckTriple(p_current_direction, p_this_direction, p_direction_changes, p_current_delta, p_previous_delta, p_angle_sign, p_second_corner, p_third_corner);
  
  /* Decide on polygon type given accumulated status */
  if (direction_changes > 2)
    return angle_sign ? 1 : 2;
  
  if (angle_sign > 0) return 4;
  if (angle_sign < 0) return 5;
  return 3;
  
}

static double determinant(double matrix[3][3]){
  return matrix[0][0] * matrix[1][1] * matrix[2][2] + matrix[0][1] * matrix[1][2] * matrix[2][0] + matrix[0][2] * matrix[1][0] * matrix[2][1] - matrix[0][2] * matrix[1][1] * matrix[2][0] - matrix[0][1] * matrix[1][0] * matrix[2][2] - matrix[0][0] * matrix[1][2] * matrix[2][1];
}

static void find_unit_normal(double a[3], double b[3], double c[3], double * normal){
  
  double matrix_for_x[3][3] = {{1, a[1], a[2]},
			       {1, b[1], b[2]},
			       {1, c[1], c[2]}			 
  };

  double x = determinant(matrix_for_x);

  double matrix_for_y[3][3] = {{a[0], 1, a[2]},
			       {b[0], 1, b[2]},
			       {c[0], 1, c[2]}
  };
  
  double y = determinant(matrix_for_y);
  
  double matrix_for_z[3][3] = {{a[0], a[1], 1},
			       {b[0], b[1], 1},
			       {c[0], c[1], 1}
  };

  double z = determinant(matrix_for_z);

  double magnitude = sqrt(x * x + y * y + z * z);

  normal[0] = x / magnitude;
  normal[1] = y / magnitude;
  normal[2] = z / magnitude;

}

static void find_unit_normal_without_det(double a[3], double b[3], double c[3], double * normal){

  double edge_one[3] = {b[0] - a[0],
			b[1] - a[1],
			b[2] - a[2]};
  
  double edge_two[3] = {c[0] - a[0],
			c[1] - a[1],
			c[2] - a[2]};

  fprintf(stderr,"Edge one is: (%f, %f, %f)\n", edge_one[0], edge_one[1], edge_one[2]);
  fprintf(stderr,"Edge two is: (%f, %f, %f)\n", edge_two[0], edge_two[1], edge_two[2]);

  double surface_normal[3];

  crossproduct_ld(edge_one, edge_two, surface_normal);

  fprintf(stderr,"Surface normal is: (%f, %f, %f)\n", surface_normal[0], surface_normal[1], surface_normal[2]);

  double magnitude = sqrt(surface_normal[0] * surface_normal[0] + surface_normal[1] * surface_normal[1] + surface_normal[2] * surface_normal[2]);

  fprintf(stderr, "Magnitude is: %f\n", magnitude);

  normal[0] = surface_normal[0] / magnitude;
  normal[1] = surface_normal[1] / magnitude;
  normal[2] = surface_normal[2] / magnitude;
}

static int no_of_duplicates_in_this_list_of_vertices(double cell_corners[], int array_length){
  
  /* Ensure that the lenght of the array is a multiple of 3. */

  if ((array_length % 3) != 0){
    return -1;
  }

  /* A brute force search for duplicate Cartesian coordinates. */

  int no_duplicates = 0;

  for (int i = 0; i < array_length; i = i + 3){
    for (int j = i + 3; j < array_length; j = j + 3 ){
      if (cell_corners[i + 0] == cell_corners[j]){
	if (cell_corners[i + 1] == cell_corners[j + 1]){
	  if (cell_corners[i + 2] == cell_corners[j + 2]){
	    no_duplicates += 1;
	  }
	}
      }
    }
  }
  return no_duplicates;
}

static double is_point_left_of_edge(double point_on_line_1[2], double point_on_line_2[2], double point[2]){

  /* 
     Computes whether a point is left of the line through point_on_line_1 and point_on_line_2. This is part of the solution to the point in polygon problem.
     Returns 0 if the point is on the line, > 0 if the point is left of the line, and < 0 if the point is right of the line.
     This algorithm is by Dan Sunday (geomalgorithms.com) and is completely free for use and modification.
  */
  
  /*
  printf("Testing whether point (%f, %f) is left of the the edge from point (%f, %f) and point (%f, %f)... ", point[0], point[1], point_on_line_1[0], point_on_line_1[1], point_on_line_2[0], point_on_line_2[1]);
  */
  double answer = ((point_on_line_2[0] - point_on_line_1[0]) * (point[1] - point_on_line_1[1]) 
		- (point[0] - point_on_line_1[0]) * (point_on_line_2[1] - point_on_line_1[1]));

  if (answer == 0){
    printf("the point lies on the edge.\n");
  }

  if (answer > 0){
    printf("the point lies left of the edge.\n");
  }

  if (answer < 0){
    printf("the point lies to the right of the edge.\n");
  }

  return answer;
}

static double is_point_left_of_edge_without_printfs(double point_on_line_1[2], double point_on_line_2[2], double point[2]){

  /* 
     Computes whether a point is left of the line through point_on_line_1 and point_on_line_2. This is part of the solution to the point in polygon problem.
     Returns 0 if the point is on the line, > 0 if the point is left of the line, and < 0 if the point is right of the line.
     This algorithm is by Dan Sunday (geomalgorithms.com) and is completely free for use and modification.
  */
  
  /*
  printf("Testing whether point (%f, %f) is left of the the edge from point (%f, %f) and point (%f, %f)... ", point[0], point[1], point_on_line_1[0], point_on_line_1[1], point_on_line_2[0], point_on_line_2[1]);
  */
  double answer = ((point_on_line_2[0] - point_on_line_1[0]) * (point[1] - point_on_line_1[1]) 
		- (point[0] - point_on_line_1[0]) * (point_on_line_2[1] - point_on_line_1[1]));
  return answer;
}

static int winding_numbers_algorithm_without_printfs(double cell_corners[], int number_corners, double point[]){
  
  /* 
     Computes whether a point is inside the bounds of a cell. This is the solution to the point in polygon problem.
     Returns 0 if the point is outside, returns 1 if the point is inside the cell.
     Based on an algorithm by Dan Sunday (geomalgorithms.com). His algorithm is completely free for use and modification.
  */
  
  int winding_number = 0;
  
  for (int i = 0;  i < number_corners; i++){
    if (cell_corners[i * 2 + 1] <= point[1]){
      if (cell_corners[(i + 1) * 2 + 1] > point[1]){
	
	double point_on_edge_1[2] = {cell_corners[i * 2 + 0], cell_corners[i * 2 + 1]};
	double point_on_edge_2[2] = {cell_corners[(i + 1) * 2 + 0], cell_corners[(i + 1) * 2 + 1]};

	if (is_point_left_of_edge_without_printfs(point_on_edge_1, point_on_edge_2, point) > 0){
	  winding_number++;
	}
      }       
    }
    else { 
      if (cell_corners[(i + 1) * 2 + 1] <= point[1]){
	
	double point_on_edge_1[2] = {cell_corners[i * 2 + 0], cell_corners[i * 2 + 1]};
	double point_on_edge_2[2] = {cell_corners[(i + 1) * 2 + 0], cell_corners[(i + 1) * 2 + 1]};

	if (is_point_left_of_edge_without_printfs(point_on_edge_1, point_on_edge_2, point) < 0){
	  winding_number--;
	}
      }
    }
  }
  return winding_number;
}

static int winding_numbers_algorithm_with_printfs(double cell_corners[], int number_corners, double point[]){
  
  /* 
     Computes whether a point is inside the bounds of a cell. This is the solution to the point in polygon problem.
     Returns 0 if the point is outside, returns 1 if the point is inside the cell.
     Based on an algorithm by Dan Sunday (geomalgorithms.com). His algorithm is completely free for use and modification.
  */

  printf("Checking if the presumed center point lies within the bounds of the cell ... \n\n");
  
  int winding_number = 0;
  
  for (int i = 0;  i < number_corners; i++){
    printf("Edge number %u from vertex (%f, %f) to vertex (%f, %f):\n\n", i+1, cell_corners[i * 2 + 0], cell_corners[i * 2 + 1], cell_corners[(i + 1) * 2 + 0], cell_corners[(i + 1) * 2 + 1]);
    if (cell_corners[i * 2 + 1] <= point[1]){
      if (cell_corners[(i + 1) * 2 + 1] > point[1]){
	printf("There is an upward crossing of the ray y = %f.\n", point[1]);

	double point_on_edge_1[2] = {cell_corners[i * 2 + 0], cell_corners[i * 2 + 1]};
	double point_on_edge_2[2] = {cell_corners[(i + 1) * 2 + 0], cell_corners[(i + 1) * 2 + 1]};

	if (is_point_left_of_edge(point_on_edge_1, point_on_edge_2, point) > 0){
	  winding_number++;
	  printf("The upward edge crossing happened on the right side of the presumed center point. The new winding number is increased to  %u.\n", winding_number);
	}
	else {
	  printf("The downward edge crossing happened on the left side of the presumed center point. The winding number remains unchanged.\n");
	}
      } 
      else {
	printf("There is NO crossing of the ray y = %f.\n", point[1]);
      }
    }
    else { 
      if (cell_corners[(i + 1) * 2 + 1] <= point[1]){
	printf("There is a downward crossing of the ray y = %f.\n", point[1]);

	double point_on_edge_1[2] = {cell_corners[i * 2 + 0], cell_corners[i * 2 + 1]};
	double point_on_edge_2[2] = {cell_corners[(i + 1) * 2 + 0], cell_corners[(i + 1) * 2 + 1]};

	if (is_point_left_of_edge(point_on_edge_1, point_on_edge_2, point) < 0){
	  winding_number--;
	  printf("The downward edge crossing happened on the right side of the presumed center point. The new winding number is decreased to %u.\n", winding_number);
	}
	else {
	  printf("The downward edge crossing happened on the left side of the presumed center point. The winding number remains unchanged.\n");
	}
      }
      else {
	printf("There is NO crossing of the ray y = %f.\n", point[1]);
      }
    }
    printf("\n");
  }
  return winding_number;
}

static double perp_dot_product(double vector_one[2], double vector_two[2]){

  /* 
     The perp-dot product is used for testing if a simple polygon is convex. From Hill, F. S. Jr. "The Pleasures of 'Perp Dot' Products." Chapter II.5 in Graphics Gems IV, Academic Press, 1994
     It uses a vector that is perpendicular to vector_one (rotated 90 degrees counterclockwise) to calculate the dot product with vector_two.
     It is positive if vector_one is less than 90 degrees away from vector_two indicating a left turn between vector_one and vector_two, and is negative otherwise.
     If all edges of a simple polygon wind the same way, it is convex.
  */
  
  return (vector_one[0] * vector_two[1]) - (vector_one[1] - vector_two[0]);

}

static double sign(double x){

  /* Is +1 if x is positive, -1 if x is negative and 0 if x is zero.*/

  return (x > 0) -  (x < 0);
}

static int is_simple_polygon_convex(double cell_corners[], int number_corners){

  /* Uses the perp-dot product to tell if a simple polygon, a cell, is convex. */

  double direction = 0;

  for (int i = 0; i < number_corners - 2; i++){

    /* Tests in which direction the polygon winds when walking along its edges. Does so for all edges of the polygon. */

    double turns_to = (cell_corners[i * 2 + 0] - cell_corners[(i + 1) * 2 + 0]) * (cell_corners[(i + 1) * 2 + 1] - cell_corners[(i + 2) * 2 + 1]) - (cell_corners[i * 2 + 1] - cell_corners[(i + 1) * 2 + 1]) * (cell_corners[(i + 1) * 2 + 0] - cell_corners[(i + 2) * 2 + 0]); 

    /* In the first iteration the direction of winding of the entire polygon is set. Better not be 0.*/

    if (i == 1){
      direction = turns_to;
    }

    if (sign(direction) != sign(turns_to)){
      if (direction != 0){
	return 0;
      }
    }
    else{
      direction = turns_to;
    }      
  }
  return 1;
}

static double calculate_the_polygon_area(double cell_corners[], int number_corners){

  /* This algorithm is based on the calculation from Wolfram Mathworld Polygon Area. It results in the area of planar non-self-intersecting polygon. */
  
  double twice_the_polygon_area = 0;

  for (int i = 0; i < number_corners - 1; i++){
    
    twice_the_polygon_area += (cell_corners[i * 2 + 0] * cell_corners[(i + 1) * 2 + 1]) - (cell_corners[(i + 1) * 2 + 0] * cell_corners[i * 2 + 1]);
       
  }
  return twice_the_polygon_area / 2;
}

static int are_polygon_vertices_arranged_in_clockwise_order(double cell_area){

  /* A negative area indicates a clockwise arrangement of vertices, a positive area a counterclockwise arrangement. There should be an area to begin with. */
  
  if (cell_area > 0){
    return 0;
  }
  if (cell_area < 0){
    return 1;
  }
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
    This function says whether the point are orientated clockwise
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
     This function performs the following tests on each cell of a given grid:

     1) it tests whether all cells are convex and  all cell bounds have the same orientation, i.e. the corners of the cell are in clockwise or counterclockwise order
     2) it tests whether the center point is within the bounds of the cell
  */
  
  double center_point_in_Euclidean_space[3];
  double cell_corners_in_Euclidean_space_open_cell[ncorner * 3];
  double no_cells_with_a_specific_no_of_corners[ncorner];
  double no_cells_with_classification_no[6];
  double corner_coordinates[3];
  double second_corner_coordinates[3];
  double third_corner_coordinates[3];
  double surface_normal_of_the_cell[3];

  double no_of_coordinates_ignored[4];
  
  int cell_no = 0;
  int cell_class = 0;
  int corner_no = 0;
  int actual_number_of_corners = 0;
  int no_of_cells_with_duplicates = 0;
  int no_convex_cells = 0;
  int no_clockwise_cells = 0;
  int no_counterclockwise_cells = 0;
  int winding_number = 0;
  int no_of_cells_with_center_points_within_their_bounds = 0;  
  int coordinate_to_ignore = 0;
  
  int no_inverted_results = 0;
  int clockwise_match = 0;
  int counterclockwise_match = 0;

  double abs_x = 0; 
  double abs_y = 0; 
  double abs_z = 0;
  
  double * p_surface_normal_of_the_cell;
  p_surface_normal_of_the_cell = &surface_normal_of_the_cell[0];


  /* 
     Latitude and longitude are spherical coordinates on a unit circle. Each such coordinate tuple is transformed into a triple of Cartesian coordinates in Euclidean space. 
     This is first done for the presumed center point of the cell and then for all the corners of the cell. LLtoXYZ is defined in clipping/geometry.h 
  */
  

  for (cell_no = 0; cell_no < 4; cell_no++)
    {
      printf("\n\n");
      printf("Cell Number %d:\n\n", cell_no + 1);
      printf("Euclidean coordinates of the still open cell are:\n\n");
         
      /* Conversion of center point spherical coordinates to Cartesian coordinates. */

      LLtoXYZ_deg(grid_center_lon[cell_no], grid_center_lat[cell_no], center_point_in_Euclidean_space);
      
      for (corner_no = 0; corner_no < ncorner; corner_no++)
	{	  
	  /* Conversion of corner spherical coordinates to Cartesian coordinates. */

	  LLtoXYZ_deg(grid_corner_lon[cell_no * ncorner + corner_no], grid_corner_lat[cell_no * ncorner + corner_no], corner_coordinates);
	  
	  /* The components of the result vector are appended to the list of cell corner coordinates. */
	  
	  cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + 0] = corner_coordinates[0];	  
	  cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + 1] = corner_coordinates[1];	  
	  cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + 2] = corner_coordinates[2];	  
	  

	  printf("(%f, %f, %f)\n", cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + 0], cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + 1], cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + 2]);
	}
      
      printf("\n");
      printf("Euclidean coordinates of the center point are:\n\n");
      printf("(%f, %f, %f)\n\n", center_point_in_Euclidean_space[0], center_point_in_Euclidean_space[1], center_point_in_Euclidean_space[2]);

      printf("\n");
      
      /* 
	 Not all cells have the same number of corners. The array, however, has ncorner * 3  values for each cell, where ncorner is the maximum number of corners. Unused values have been filled with the values of the final cell.
	 The following identifies the surplus corners and gives the correct length of the cell.
      */
      
      actual_number_of_corners = ncorner;

      for (corner_no = ncorner - 1; corner_no > 0; corner_no--){
	if (cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + 0] == cell_corners_in_Euclidean_space_open_cell[(corner_no - 1) * 3 + 0]){
	  if (cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + 1] == cell_corners_in_Euclidean_space_open_cell[(corner_no - 1) * 3 + 1]){
	    if (cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + 2] == cell_corners_in_Euclidean_space_open_cell[(corner_no - 1) * 3 + 2]){
	      actual_number_of_corners = actual_number_of_corners - 1;
	    }
	  }
	} else {
	  break;
	}	
      }            
      
      printf("The actual number of corners is: %u\n\n", actual_number_of_corners);

      no_cells_with_a_specific_no_of_corners[actual_number_of_corners - 1] += 1;
      
      /* If there are less than three corners in the cell, it is unusable and considered degenerate. No area can be computed. */
      
      if (actual_number_of_corners < 3){
	printf("Less than three corners found. Skipping!\n\n");
	continue;
      }
      
      /* Checks if there are any duplicate vertices in the list of corners. Note that the last (additional) corner has not been set yet. */

      if (no_of_duplicates_in_this_list_of_vertices(cell_corners_in_Euclidean_space_open_cell, actual_number_of_corners * 3) > 0){
	no_of_cells_with_duplicates += 1;
	printf("Duplicate corner found. Skipping!\n\n");
	continue;
      }

      /* We are creating a closed polygon/cell by setting the additional last corner to be the same as the first one. */

      double cell_corners_in_Euclidean_space[(actual_number_of_corners + 1) * 3];

      for (corner_no = 0; corner_no < actual_number_of_corners; corner_no++){
	cell_corners_in_Euclidean_space[corner_no * 3 + 0] = cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + 0];
	cell_corners_in_Euclidean_space[corner_no * 3 + 1] = cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + 1];
	cell_corners_in_Euclidean_space[corner_no * 3 + 2] = cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + 2];
      }

      cell_corners_in_Euclidean_space[actual_number_of_corners * 3 + 0] = cell_corners_in_Euclidean_space[0];
      cell_corners_in_Euclidean_space[actual_number_of_corners * 3 + 1] = cell_corners_in_Euclidean_space[1];
      cell_corners_in_Euclidean_space[actual_number_of_corners * 3 + 2] = cell_corners_in_Euclidean_space[2];

      /* Takes the first three corners/vertices of the cell and calculates the unit normal via determinants. */
      
      corner_coordinates[0] = cell_corners_in_Euclidean_space[0];
      corner_coordinates[1] = cell_corners_in_Euclidean_space[1];
      corner_coordinates[2] = cell_corners_in_Euclidean_space[2];

      second_corner_coordinates[0] = cell_corners_in_Euclidean_space[3 + 0];
      second_corner_coordinates[1] = cell_corners_in_Euclidean_space[3 + 1];
      second_corner_coordinates[2] = cell_corners_in_Euclidean_space[3 + 2];

      third_corner_coordinates[0] = cell_corners_in_Euclidean_space[6 + 0];
      third_corner_coordinates[1] = cell_corners_in_Euclidean_space[6 + 1];
      third_corner_coordinates[2] = cell_corners_in_Euclidean_space[6 + 2];

      
      find_unit_normal(corner_coordinates, second_corner_coordinates, third_corner_coordinates, p_surface_normal_of_the_cell);

      /* The surface normal is used to choose the coordinate to ignore. */
      
      if (surface_normal_of_the_cell[0] > 0){
	abs_x = surface_normal_of_the_cell[0];
      } else {
	abs_x = surface_normal_of_the_cell[0] * (-1);
      }
      
      if (surface_normal_of_the_cell[1] > 0){
	abs_y = surface_normal_of_the_cell[1];
      } else {
	abs_y = surface_normal_of_the_cell[1] * (-1);
      }
      
      if (surface_normal_of_the_cell[2] > 0){
	abs_z = surface_normal_of_the_cell[2];
      } else {
	abs_z = surface_normal_of_the_cell[2] * (-1);
      }
      
      coordinate_to_ignore = 3;
      
      if (abs_x > abs_y){
	if (abs_x > abs_z){
	  coordinate_to_ignore = 1;
	}
      } else {
	if (abs_y > abs_z){
	  coordinate_to_ignore = 2;
	}
      }

      no_of_coordinates_ignored[coordinate_to_ignore] += 1;

      /* The remaining two-dimensional coordinates are extracted into one array for all the cell's corners and into one array for the center point. */

      printf("Coordinate to ignore is: %u\n\n",coordinate_to_ignore);

      double cell_corners_plane_projection[(actual_number_of_corners +1) * 2];
      double center_point_plane_projection[2];

      
      /* The following projection on the plane that two coordinate axes lie on changes the arrangement of the polygon vertices if the coordinate to be ignored along the third axis is smaller than 0.
	 In this case, the result of the computation of the orientation of vertices needs to be inverted. Clockwise becomes counterclockwise and vice versa. */

      int invert_result = 0;
  
      if (cell_corners_in_Euclidean_space[coordinate_to_ignore - 1] < 0){
	invert_result = 1;
	no_inverted_results += 1;
      }
      
      switch(coordinate_to_ignore){
      case 1:
	for(corner_no = 0; corner_no <= actual_number_of_corners; corner_no++){
	  cell_corners_plane_projection[corner_no * 2 + 0] = cell_corners_in_Euclidean_space[corner_no * 3 + 1];
	  cell_corners_plane_projection[corner_no * 2 + 1] = cell_corners_in_Euclidean_space[corner_no * 3 + 2];
	}
	center_point_plane_projection[0] = center_point_in_Euclidean_space[1];
	center_point_plane_projection[1] = center_point_in_Euclidean_space[2];		
	break;
      case 2:
	for(int corner_no = 0; corner_no <= actual_number_of_corners; corner_no++){
	  cell_corners_plane_projection[corner_no * 2 + 0] = cell_corners_in_Euclidean_space[corner_no * 3 + 0];
	  cell_corners_plane_projection[corner_no * 2 + 1] = cell_corners_in_Euclidean_space[corner_no * 3 + 2];
	}
	center_point_plane_projection[0] = center_point_in_Euclidean_space[0];
	center_point_plane_projection[1] = center_point_in_Euclidean_space[2];	
	break;
      case 3:
	for(int corner_no = 0; corner_no <= actual_number_of_corners; corner_no++){
	  cell_corners_plane_projection[corner_no * 2 + 0] = cell_corners_in_Euclidean_space[corner_no * 3 + 0];
	  cell_corners_plane_projection[corner_no * 2 + 1] = cell_corners_in_Euclidean_space[corner_no * 3 + 1];
	}
	center_point_plane_projection[0] = center_point_in_Euclidean_space[0];
	center_point_plane_projection[1] = center_point_in_Euclidean_space[1];	
	break;
      }

      printf("Cell coordinates of the closed cell in 2D are:\n\n");

      for(corner_no = 0; corner_no <= actual_number_of_corners; corner_no++){
	printf("(%f, %f)\n", cell_corners_plane_projection[corner_no * 2 + 0], cell_corners_plane_projection[corner_no * 2 + 1]);
      }
      
      printf("\n");

      printf("Center point coordinates in 2D are:\n\n");
      printf("(%f, %f)", center_point_plane_projection[0], center_point_plane_projection[1]);
      printf("\n\n");

      /* Checking for convexity of the cell. */

      if(is_simple_polygon_convex(cell_corners_plane_projection, actual_number_of_corners +1)){
	no_convex_cells += 1;
      }

      

      /* Checking the arrangement of cell vertices. */

      double polygon_area = calculate_the_polygon_area(cell_corners_plane_projection, actual_number_of_corners);

      fprintf(stderr,"The polygon area is: %f\n\n", polygon_area);

      int is_clockwise = are_polygon_vertices_arranged_in_clockwise_order(polygon_area);
      
      if((is_clockwise == 1) && (invert_result == 1)){
	clockwise_match += 1;
      }

      if((is_clockwise == 0) && (invert_result == 1)){
	counterclockwise_match += 1;
      }
      
      fprintf(stderr,"Is_clockwise: %u\n\n", is_clockwise);
      fprintf(stderr,"Invert result?: %u\n\n", invert_result);

      if(invert_result == 1){
	if(is_clockwise == 1){
	  is_clockwise = 0;
	} 
	
	if(is_clockwise == 0){
	  is_clockwise = 1;
	}
      }
      
      fprintf(stderr,"Is_clockwise: %u\n\n", is_clockwise);
      
      if(is_clockwise){
	no_clockwise_cells += 1;
      }
      else{
	no_counterclockwise_cells +=1;
      }
      
     

      continue;

      /* The winding numbers algorithm is used to test whether the presumed center point is within the bounds of the cell. */
        
      winding_number = winding_numbers_algorithm_with_printfs(cell_corners_plane_projection, actual_number_of_corners + 1, center_point_plane_projection);

      if (winding_number == 0){
	printf("The presumed center point lies OUTSIDE the bounds of the cell.\n\n\n\n");
      } else {
	printf("The presumed center point lies INSIDE the bounds of the cell.\n\n\n\n");
	no_of_cells_with_center_points_within_their_bounds += 1;
      }

      printf("Cell coordinates of the closed cell in 2D are:\n\n");

      for(corner_no = 0; corner_no <= actual_number_of_corners; corner_no++){
	printf("(%f, %f)\n", cell_corners_plane_projection[corner_no * 2 + 0], cell_corners_plane_projection[corner_no * 2 + 1]);
      }
      
      printf("\n");

      printf("Center point coordinates in 2D are:\n\n");
      printf("(%f, %f)", center_point_plane_projection[0], center_point_plane_projection[1]);
      printf("\n\n");
     
    

      

      /* 
	 Determining the class of a polygon cell.
	 The following polygon classifications are possible:
	 
	 NotConvex = 1 
	 NotConvexDegenerate = 2
	 ConvexDegenerate = 3
	 ConvexCCW = 4
	 ConvexCW = 5	       
      */

      cell_class = 0;

      fprintf(stderr,"Calling x_classify_polygon for cell no: %d\n\n", cell_no + 1);
      cell_class = x_classify_polygon(actual_number_of_corners, cell_corners_plane_projection);
      
      no_cells_with_classification_no[cell_class] +=  1;
      
      if (cell_class == 1) printf("Cell classification resulted in polygon class: Not Convex.\n\n");
      if (cell_class == 2) printf("Cell classification resulted in polygon class: Not Convex Degenerate.\n\n");
      if (cell_class == 3) printf("Cell classification resulted in polygon class: Convex Degenerate.\n\n");
      if (cell_class == 4) printf("Cell classification resulted in polygon class: Convex Counterclockwise.\n\n");
      if (cell_class == 5) printf("Cell classification resulted in polygon class: Convex Clockwise.\n\n");
    

    }

  printf("\n\n");
  printf("There are %d cells in all.", gridsize );

  printf("\n\n");
  printf("%d cells are convex.", no_convex_cells);
  
  printf("\n\n");
  fprintf(stderr,"%d cells have vertices arranged in a clockwise manner.", no_clockwise_cells);
  
  printf("\n\n");
  fprintf(stderr,"%d cells have vertices arranged in a counterclockwise manner.", no_counterclockwise_cells);

  printf("\n\n");
  printf("Number of results that should be inverted: %d", no_inverted_results);

  printf("\n\n");
  printf("Number of clockwise results that should be inverted: %d", clockwise_match);

  printf("\n\n");
  printf("Number of counterclockwise results that should be inverted: %d", counterclockwise_match);
  
  printf("\n\n");
  printf("Number of cells with center points within their bounds: %d", no_of_cells_with_center_points_within_their_bounds);
  

  printf("\n\n");
  printf("Number of cells with a certain number of corners:\n\n");
  for(int i = 0; i < ncorner; i++){
    printf("%d corners in %f cells.\n", i+1, no_cells_with_a_specific_no_of_corners[i]);
  }

  printf("\n");
  printf("Number of cells with duplicate vertices: %u", no_of_cells_with_duplicates);

  /*

  printf("\n\n");
  printf("Number of cells with a certain cell classification:\n\n");
  for(int i = 1; i < 6; i++){
    printf("Classification number %u in %f cells.\n", i, no_cells_with_classification_no[i]);
  }

  */
  
  printf("\n\n");
  printf("Number of cells with an ignored coordinate:\n\n");
  for(int i = 1; i < 4; i++){
    printf("Ignored coordinate no. %u in %f cells.\n", i, no_of_coordinates_ignored[i]);
  }

 
 
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
