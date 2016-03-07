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

struct axis {
  double axis_vector[3];
};

static double Euclidean_norm (double a[]) {

  /* Computes the Euclidean norm of a vector given in Cartesian coordinates. */

  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

static struct axis divide_by_scalar (struct axis a, double scalar) {

  /* Component-wise scalar division of a three dimensional axis vector given in Cartesian coordinates. */

  a.axis_vector[0] = a.axis_vector[0]/scalar;
  a.axis_vector[1] = a.axis_vector[1]/scalar;
  a.axis_vector[2] = a.axis_vector[2]/scalar;

  return a;
}

static struct axis normalize_vector(struct axis a){

  /* Normalizes an axis vector a by dividing it though its magnitude. */

  a = divide_by_scalar(a, Euclidean_norm(a.axis_vector));

  return a;
}

static double sign(double x){

  /* Is +1 if x is positive, -1 if x is negative and 0 if x is zero.*/

  return (x > 0) -  (x < 0);
}

static double determinant(double matrix[3][3]){
  return matrix[0][0] * matrix[1][1] * matrix[2][2] + matrix[0][1] * matrix[1][2] * matrix[2][0] + matrix[0][2] * matrix[1][0] * matrix[2][1] - matrix[0][2] * matrix[1][1] * matrix[2][0] - matrix[0][1] * matrix[1][0] * matrix[2][2] - matrix[0][0] * matrix[1][2] * matrix[2][1];
}

static void find_unit_normal(double a[3], double b[3], double c[3], double normal[3]){
  
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



static int vector_in_octant_no(double vector[3]){

  /*
    Analogue to quadrants in two-dimensional space there are eight octants in Euclidean space. This function returns the octant a given vector is in.
    These are the octant numbers and  their characteristic coordinate signs (x, y, z):

    1)    top-front-right      (+, +, +)
    2)    top-back-right       (-, +, +)
    3)    top-back-left        (-, -, +)
    4)    top-front-left       (+, -, +)
    5)    bottom-front-right   (+, +, -)
    6)    bottom-back-right    (-, +, -)
    7)    bottom-back-left     (-, -, -)
    8)    bottom-fron-left     (+, -, -)
   */

  if ((vector[0] >= 0) && (vector[1] >= 0) && (vector[2] >= 0)){
    return 1;
  } 

  if ((vector[0] < 0) && (vector[1] >= 0) && (vector[2] >= 0)){
    return 2;
  } 

  if ((vector[0] < 0) && (vector[1] < 0) && (vector[2] >= 0)){
    return 3;
  } 

  if ((vector[0] >= 0) && (vector[1] < 0) && (vector[2] >= 0)){
    return 4;
  } 

  if ((vector[0] >= 0) && (vector[1] >= 0) && (vector[2] < 0)){
    return 5;
  } 

  if ((vector[0] < 0) && (vector[1] >= 0) && (vector[2] < 0)){
    return 6;
  } 

  if ((vector[0] < 0) && (vector[1] < 0) && (vector[2] < 0)){
    return 7;
  } 

  if ((vector[0] >= 0) && (vector[1] < 0) && (vector[2] < 0)){
    return 8;
  } 
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

static struct axis compute_the_new_z_axis(double cell_corners_in_Euclidean_space[]){

  /* Takes the first three corners/vertices of the cell and computes two edges originating at the first corner/vertex. These two edges are on the same plane as all the other vertices. THERE NEED TO BE AT LEAST THREE CORNERS. */

  double edge_one[3] = {(cell_corners_in_Euclidean_space)[3 + 0] -(cell_corners_in_Euclidean_space)[0],
			(cell_corners_in_Euclidean_space)[3 + 1] -(cell_corners_in_Euclidean_space)[1],
			(cell_corners_in_Euclidean_space)[3 + 2] -(cell_corners_in_Euclidean_space)[2]};
  
  double edge_two[3] = {(cell_corners_in_Euclidean_space)[6 + 0] - (cell_corners_in_Euclidean_space)[0],
			(cell_corners_in_Euclidean_space)[6 + 1] - (cell_corners_in_Euclidean_space)[1],
			(cell_corners_in_Euclidean_space)[6 + 2] - (cell_corners_in_Euclidean_space)[2]};

  printf("Edge one is: (%f, %f, %f)\n", edge_one[0], edge_one[1], edge_one[2]);
  printf("Edge two is: (%f, %f, %f)\n\n", edge_two[0], edge_two[1], edge_two[2]);

  struct axis new_z_axis;

  /* The cross product of the two edges is the surface normal and the new z-axis. crossproduct_d is defined in clipping/geometry.h */
  
  crossproduct_ld(edge_one, edge_two, new_z_axis.axis_vector);

  new_z_axis = normalize_vector(new_z_axis);

  /* Regardless of whether the vertices are arranged in clockwise or counterclockwise order, the z-axis is to point inward. */

  double vector_pointing_from_origin_to_first_cell_vertex[3] = {cell_corners_in_Euclidean_space[0], cell_corners_in_Euclidean_space[1], cell_corners_in_Euclidean_space[2]};
  
  int axis_points_inward = 0;

  if (sign(dotproduct(new_z_axis.axis_vector, vector_pointing_from_origin_to_first_cell_vertex)) < 0){
    axis_points_inward = 1;
  }
  
  if (axis_points_inward == 0){
    new_z_axis.axis_vector[0] = new_z_axis.axis_vector[0] * (-1);
    new_z_axis.axis_vector[1] = new_z_axis.axis_vector[1] * (-1);
    new_z_axis.axis_vector[2] = new_z_axis.axis_vector[2] * (-1);
  }
  
  return new_z_axis;
}



 
static struct axis compute_the_new_y_axis(struct axis new_z_axis){

  /* Then the new y-axis is the result of the cross product of the new z-axis and the old x-axis. crossproduct_d is defined in clipping/geometry.h */
  
  struct axis old_x_axis;
 
  old_x_axis.axis_vector[0] = 1;
  old_x_axis.axis_vector[1] = 0;
  old_x_axis.axis_vector[2] = 0;

  struct axis new_y_axis;
  
  /* The y-axis is always to point 'upward'.*/

  if (new_z_axis.axis_vector[1] >= 0){
    crossproduct_ld(old_x_axis.axis_vector, new_z_axis.axis_vector, new_y_axis.axis_vector);
  } else {
    crossproduct_ld(new_z_axis.axis_vector, old_x_axis.axis_vector, new_y_axis.axis_vector);
  }

  new_y_axis = normalize_vector(new_y_axis);

  return new_y_axis;
}

static struct axis compute_the_new_x_axis(struct axis new_z_axis, struct axis new_y_axis){

  /* The new x-axis is the result of the cross product of the new z-axis and the new y-axis. crossproduct_d is defined in clipping/geometry.h */
  
  struct axis new_x_axis;

  crossproduct_ld(new_z_axis.axis_vector, new_y_axis.axis_vector, new_x_axis.axis_vector);
  
  new_x_axis = normalize_vector(new_x_axis);

  return new_x_axis;
}

static void transform_Euclidean_corner_coordinates_onto_the_cell_plane(struct axis new_x_axis, struct axis new_y_axis, struct axis new_z_axis, double (*p_cell_corners_in_Euclidean_space)[30], double (*p_cell_corners_on_cell_plane)[20], int ncorners){

  /* Since the three axes form an orthonormal base, the inverse and the transpose of their transformation matrix are the same. */

  double transposed_transformation_matrix[3][3];

  transposed_transformation_matrix[0][0] = new_x_axis.axis_vector[0];
  transposed_transformation_matrix[0][1] = new_y_axis.axis_vector[0];
  transposed_transformation_matrix[0][2] = new_z_axis.axis_vector[0];

  transposed_transformation_matrix[1][0] = new_x_axis.axis_vector[1];
  transposed_transformation_matrix[1][1] = new_y_axis.axis_vector[1];
  transposed_transformation_matrix[1][2] = new_z_axis.axis_vector[1];

  transposed_transformation_matrix[2][0] = new_x_axis.axis_vector[2];
  transposed_transformation_matrix[2][1] = new_y_axis.axis_vector[2];
  transposed_transformation_matrix[2][2] = new_z_axis.axis_vector[2];

  /* The points from each corner are now transformed from the world frame to the cell frame.  */

  for (int corner_no = 0; corner_no < ncorners; corner_no++){
    
    /* Multiplying the transformation with the vertex coordinates. */
    
    double first_row_multiplicand_vector[3] = {transposed_transformation_matrix[0][0], transposed_transformation_matrix[1][0], transposed_transformation_matrix[2][0]};    
    double second_row_multiplicand_vector[3] = {transposed_transformation_matrix[0][1], transposed_transformation_matrix[1][1], transposed_transformation_matrix[2][1]};
    double third_row_multiplicand_vector[3] = {transposed_transformation_matrix[0][2], transposed_transformation_matrix[1][2], transposed_transformation_matrix[2][2]};

    double multiplier_vector[3] = {(*p_cell_corners_in_Euclidean_space)[(corner_no * 3) + 0], (*p_cell_corners_in_Euclidean_space)[(corner_no * 3) + 1], (*p_cell_corners_in_Euclidean_space)[(corner_no * 3) + 2]};

    (*p_cell_corners_on_cell_plane)[(corner_no * 2) + 0] = dotproduct(first_row_multiplicand_vector, multiplier_vector);
    (*p_cell_corners_on_cell_plane)[(corner_no * 2) + 1] = dotproduct(second_row_multiplicand_vector, multiplier_vector);
  }
  
  /* The last vertex is set to be the same as the first one in order to close the cell.*/
  
  (*p_cell_corners_on_cell_plane)[(ncorners * 2) + 0] = (*p_cell_corners_on_cell_plane)[0];
  (*p_cell_corners_on_cell_plane)[(ncorners * 2) + 1] = (*p_cell_corners_on_cell_plane)[1];
  
}

static void transform_Euclidean_center_coordinates_onto_the_cell_plane(struct axis new_x_axis, struct axis new_y_axis, struct axis new_z_axis, double (*p_center_point_in_Euclidean_space)[3], double (*p_center_point_on_cell_plane)[2]){
  
  /* Since the three axes form an orthonormal base, the inverse and the transpose of their transformation matrix are the same. */
  
  double transposed_transformation_matrix[3][3];
  
  transposed_transformation_matrix[0][0] = new_x_axis.axis_vector[0];
  transposed_transformation_matrix[0][1] = new_y_axis.axis_vector[0];
  transposed_transformation_matrix[0][2] = new_z_axis.axis_vector[0];

  transposed_transformation_matrix[1][0] = new_x_axis.axis_vector[1];
  transposed_transformation_matrix[1][1] = new_y_axis.axis_vector[1];
  transposed_transformation_matrix[1][2] = new_z_axis.axis_vector[1];

  transposed_transformation_matrix[2][0] = new_x_axis.axis_vector[2];
  transposed_transformation_matrix[2][1] = new_y_axis.axis_vector[2];
  transposed_transformation_matrix[2][2] = new_z_axis.axis_vector[2];

  double first_row_multiplicand_vector[3] = {transposed_transformation_matrix[0][0], transposed_transformation_matrix[1][0], transposed_transformation_matrix[2][0]};    
  double second_row_multiplicand_vector[3] = {transposed_transformation_matrix[0][1], transposed_transformation_matrix[1][1], transposed_transformation_matrix[2][1]};
  double third_row_multiplicand_vector[3] = {transposed_transformation_matrix[0][2], transposed_transformation_matrix[1][2], transposed_transformation_matrix[2][2]};
  
  double multiplier_vector[3] = {(*p_center_point_in_Euclidean_space)[0], (*p_center_point_in_Euclidean_space)[1], (*p_center_point_in_Euclidean_space)[2]};
  
  (*p_center_point_on_cell_plane)[0] = dotproduct(first_row_multiplicand_vector, multiplier_vector);
  (*p_center_point_on_cell_plane)[1] = dotproduct(second_row_multiplicand_vector, multiplier_vector);
  
  /*
  printf("The coordinates of the presumed center point of the cell on the cell plane are: (%f, %f)\n\n", (*p_center_point_on_cell_plane)[0], (*p_center_point_on_cell_plane)[1]);
  */
}


static void project_Euclidean_corner_coordinates_onto_the_cell_plane(struct axis new_x_axis, struct axis new_y_axis, double (*p_cell_corners_in_Euclidean_space)[30], double (*p_cell_corners_on_cell_plane)[20], int ncorners){

  /* All corner points are projected onto the new x- and y-axes. dotproduct is defined in clipping/clipping.c */

  double corner_vector[3];

  for (int corner_no = 0; corner_no < ncorners; corner_no++){
    
    for (int vector_component = 0; vector_component < 3; ++vector_component){
      corner_vector[vector_component] = (*p_cell_corners_in_Euclidean_space)[(corner_no * 3) + vector_component];
    }
    
    (*p_cell_corners_on_cell_plane)[(corner_no * 2) + 0] = dotproduct(new_x_axis.axis_vector, corner_vector);
    (*p_cell_corners_on_cell_plane)[(corner_no * 2) + 1] = dotproduct(new_y_axis.axis_vector, corner_vector);  
    
  }
      
  /* The last vertex is set to be the same as the first one in order to close the cell.*/

  (*p_cell_corners_on_cell_plane)[(ncorners * 2) + 0] = (*p_cell_corners_on_cell_plane)[0];
  (*p_cell_corners_on_cell_plane)[(ncorners * 2) + 1] = (*p_cell_corners_on_cell_plane)[1];
  
  /*
  printf("The coordinates of the cell vertices on the cell plane are: ");

  for (int corner_no =0; corner_no <= ncorners; corner_no++){
    printf("(%f, %f) ", cell_corners_on_cell_plane[(corner_no * 2) + 0], cell_corners_on_cell_plane[(corner_no * 2) + 1]);
  }

  printf("\n\n");
  */
}

static void project_Euclidean_center_coordinates_onto_the_cell_plane(struct axis new_x_axis, struct axis new_y_axis, double (*p_center_point_in_Euclidean_space)[3], double (*p_center_point_on_cell_plane)[2]){
  
  /* The center point is projected onto the new x- and -y-axes. dotproduct is defined in clipping/clipping.c  */

  (*p_center_point_on_cell_plane)[0] = dotproduct((new_x_axis.axis_vector), (*p_center_point_in_Euclidean_space));
  (*p_center_point_on_cell_plane)[1] = dotproduct((new_y_axis.axis_vector), (*p_center_point_in_Euclidean_space));
  
  /*
  printf("The coordinates of the presumed center point of the cell on the cell plane are: (%f, %f)\n\n", (*p_center_point_on_cell_plane)[0], (*p_center_point_on_cell_plane)[1]);
  */
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




static int is_simple_polygon_convex(double cell_corners[], int number_corners){

  /* Uses the perp-dot product to tell if a simple polygon, a cell, is convex. */

  double direction = 0;

  for (int i = 0; i < number_corners - 1; i++){

    /* Tests in which direction edge B winds that is connected to edge A when walking along the polygon edges. Does so for all edges of the polygon. */

    double edge_a[2] = {cell_corners[(i + 1) * 2 + 0] - cell_corners[i * 2 + 0], cell_corners[(i + 1) * 2 + 1] - cell_corners[i * 2 + 1]};
    double edge_b[2] = {cell_corners[(i + 2) * 2 + 0] - cell_corners[(i + 1) * 2 + 0], cell_corners[(i + 2) * 2 + 1] - cell_corners[(i + 1) * 2 + 1]};
    
    double turns_to = sign(perp_dot_product(edge_a, edge_b));
    
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

  /* This algorithm works with non-convex and even self-intersecting polygons. It results in twice the area of the polygon. */
  
  double twice_the_polygon_area = 0;

  for (int i = 0; i < number_corners - 1; i++){
    
    twice_the_polygon_area += (cell_corners[(i + 1) * 2 + 0] - cell_corners[i * 2 + 0]) * (cell_corners[(i + 1) * 2 + 1] - cell_corners[i * 2 + 1]);
  }
  return twice_the_polygon_area / 2;
}

static int are_polygon_vertices_arranged_in_clockwise_order(double cell_area){

  /* A positive area indicates a clockwise arrangement of vertices, a negative area a counterclockwise arrangement. There should be an area to begin with. */
  
  if (cell_area > 0){
    return 1;
  }
  if (cell_area < 0){
    return 0;
  }
}

static double inefficiently_calculate_the_area_of_a_polygon_in_Euclidean_space(int number_of_vertices, double cell_corners_in_Euclidean_space[30]){

  /* number_of_vertices is the number for the open polygon. The array additionally has another vertex which is equal to the first vertex, closing the polygon. */

  double cell_area = 0;
  double surface_normal_of_the_cell[3];
  double result_of_the_cross_product[3];
  double sum_of_cross_product_results[3];
 
  int i, j, k;

  /* If there are less than three corners in the cell, it is unusable and considered degenerate. No area can be computed. */

  if (number_of_vertices < 3){
    return 0;
  }

  /* Takes the first three corners/vertices of the cell and calculates the unit normal via determinants. */

  double vertex_one[3] = {cell_corners_in_Euclidean_space[0],
			  cell_corners_in_Euclidean_space[1],
			  cell_corners_in_Euclidean_space[2]};
  
  double vertex_two[3] = {cell_corners_in_Euclidean_space[3 + 0],
			  cell_corners_in_Euclidean_space[3 + 1],
			  cell_corners_in_Euclidean_space[3 + 2]};

  double vertex_three[3] = {cell_corners_in_Euclidean_space[6 + 0],
			    cell_corners_in_Euclidean_space[6 + 1],
			    cell_corners_in_Euclidean_space[6 + 2]};

			  
  find_unit_normal(vertex_one, vertex_two, vertex_three, surface_normal_of_the_cell);
  
  /*
  printf("The surface normal is (%f, %f, %f).\n\n", surface_normal_of_the_cell[0], surface_normal_of_the_cell[1], surface_normal_of_the_cell[2]);
  */

  /* We are calculating the dot product of the surface normal and the sum of the cross products of the two vertices of each edge. The result is twice the signed area of the polygon. */

  for (i = 0; i < number_of_vertices; i++){
    
    double vertex_one[3] = {cell_corners_in_Euclidean_space[i * 3 + 0], cell_corners_in_Euclidean_space[i * 3 + 1], cell_corners_in_Euclidean_space[i * 3 + 2]};
    double vertex_two[3] = {cell_corners_in_Euclidean_space[(i + 1) * 3 + 0], cell_corners_in_Euclidean_space[(i + 1) * 3 + 1], cell_corners_in_Euclidean_space[(i + 1) * 3 + 2]};
    
    /*
    printf("The vertices used for calculation are (%f, %f, %f) and (%f, %f, %f).\n\n", vertex_one[0], vertex_one[1], vertex_one[2], vertex_two[0], vertex_two[1], vertex_two[2]);
    */

    crossproduct_ld(vertex_one, vertex_two, result_of_the_cross_product);

    /*
    printf("Their cross product is (%f, %f, %f).\n\n", result_of_the_cross_product[0], result_of_the_cross_product[1], result_of_the_cross_product[2]);
    */

    sum_of_cross_product_results[0] += result_of_the_cross_product[0];
    sum_of_cross_product_results[1] += result_of_the_cross_product[1];
    sum_of_cross_product_results[2] += result_of_the_cross_product[2];
  }

  /*
  printf("The sum of the cross product results is (%f, %f, %f).\n\n", sum_of_cross_product_results[0], sum_of_cross_product_results[1], sum_of_cross_product_results[2]);
  */

  cell_area = dotproduct(surface_normal_of_the_cell, sum_of_cross_product_results) / 2;

  /*
  printf("The cell area is %f.\n\n", cell_area);
  */
  return cell_area;
}




static double calculate_the_area_of_a_polygon_in_Euclidean_space(int number_of_vertices, double cell_corners_in_Euclidean_space[30]){

  /* number_of_vertices is the number for the open polygon. The array additionally has another vertex which is equal to the first vertex, closing the polygon. */

  double cell_area = 0;
  double abs_x, abs_y, abs_z;
  double surface_normal_of_the_cell[3];
  double magnitude_of_normal;
  int coordinate_to_ignore;
  int i, j, k;

  /* If there are less than three corners in the cell, it is unusable and considered degenerate. No area can be computed. */

  if (number_of_vertices < 3){
    return 0;
  }

  /* Takes the first three corners/vertices of the cell and calculates the unit normal via determinants. */

  double vertex_one[3] = {cell_corners_in_Euclidean_space[0],
			  cell_corners_in_Euclidean_space[1],
			  cell_corners_in_Euclidean_space[2]};
  
  double vertex_two[3] = {cell_corners_in_Euclidean_space[3 + 0],
			  cell_corners_in_Euclidean_space[3 + 1],
			  cell_corners_in_Euclidean_space[3 + 2]};

  double vertex_three[3] = {cell_corners_in_Euclidean_space[6 + 0],
			    cell_corners_in_Euclidean_space[6 + 1],
			    cell_corners_in_Euclidean_space[6 + 2]};

			  
  find_unit_normal(vertex_one, vertex_two, vertex_three, surface_normal_of_the_cell);
  
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

  /* The area of the projected cell is computed. */

  switch(coordinate_to_ignore){
  case 1:
    for (i = 1, j = 2, k = 0; i < number_of_vertices; i++, j++, k++){
      cell_area += (cell_corners_in_Euclidean_space[i * 3 + 1] * (cell_corners_in_Euclidean_space[j * 3 + 2] - cell_corners_in_Euclidean_space[k * 3 + 2]));
      break;
    }
  case 2:
    for (i = 1, j = 2, k = 0; i < number_of_vertices; i++, j++, k++){
      cell_area += (cell_corners_in_Euclidean_space[i * 3 + 2] * (cell_corners_in_Euclidean_space[j * 3 + 0] - cell_corners_in_Euclidean_space[k * 3 + 0]));
      break;
    }
  case 3:
    for (i = 1, j = 2, k = 0; i < number_of_vertices; i++, j++, k++){
      cell_area += (cell_corners_in_Euclidean_space[i * 3 + 0] * (cell_corners_in_Euclidean_space[j * 3 + 1] - cell_corners_in_Euclidean_space[k * 3 + 1]));
      break;
    }     
  }

  /* Wrap around term. */

  switch(coordinate_to_ignore){
  case 1:
    cell_area += cell_corners_in_Euclidean_space[number_of_vertices * 3 + 1] * (cell_corners_in_Euclidean_space[3 + 2] - cell_corners_in_Euclidean_space[(number_of_vertices - 1) * 3 + 2]);
    break;
  case 2:
    cell_area += cell_corners_in_Euclidean_space[number_of_vertices * 3 + 2] * (cell_corners_in_Euclidean_space[3 + 0] - cell_corners_in_Euclidean_space[(number_of_vertices - 1) * 3 + 0]);
    break;
  case 3:
    cell_area += cell_corners_in_Euclidean_space[number_of_vertices * 3 + 0] * (cell_corners_in_Euclidean_space[3 + 1] - cell_corners_in_Euclidean_space[(number_of_vertices - 1) * 3 + 1]);
    break;
  }

  /* The computed area is rescaled in order to calculate the size before projection. */

  magnitude_of_normal =  sqrt(abs_x * abs_x + abs_y * abs_y + abs_z * abs_z);

  switch(coordinate_to_ignore){
  case 1:
    cell_area *=  (magnitude_of_normal / (2 * surface_normal_of_the_cell[0]));
    break;
  case 2:
    cell_area *=  (magnitude_of_normal / (2 * surface_normal_of_the_cell[1]));
    break;
  case 3:
    cell_area *=  (magnitude_of_normal / (2 * surface_normal_of_the_cell[2]));
  }  
  return cell_area;
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
     This function performs three tests on each cell of a given grid:

     1) it tests whether all cell bounds have the same orientation, i.e. the corners of the cell are in clockwise or counterclockwise order
     2) it tests whether the cell is convex
     3) it tests whether the center point is within the bounds of the cell

     It performs these tests after longitude and latitude on the unit circle have been converted first to Cartesian coordinates in Euclidean space and subsequently to two dimensional coordinates on the plane each cell occupies.  
  */
  
  double center_point_in_Euclidean_space[3];
  double cell_corners_in_Euclidean_space_open_cell[ncorner * 3];
  double cell_corners_in_Euclidean_space[(ncorner + 1) * 3];
  double corners_of_the_final_cell[ncorner * 3];
  double cell_corner_coordinates[3];
  double center_point_on_cell_plane[2];
  double cell_corners_on_cell_plane[(ncorner + 1)  * 2];
  double no_cells_with_a_specific_no_of_corners[ncorner + 1];
  

  int cell_no;
  int corner_no;
  int orthogonal_axes = 0;
  int axes_pointing_inward = 0;
  int no_of_cells_with_duplicates = 0;
  int vector_component;

  int no_of_cells_with_vertices_arranged_in_clockwise_order = 0;
  int no_of_cells_with_vertices_arranged_in_counterclockwise_order = 0;
  int no_of_convex_cells = 0;
  int no_of_degenerate_cells = 0;
  int no_of_cells_with_center_points_within_their_bounds = 0;


  /* 
     Latitude and longitude are spherical coordinates on a unit circle. Each such coordinate tuple is transformed into a triple of Cartesian coordinates in Euclidean space. 
     This is first done for the presumed center point of the cell and then for all the corners of the cell. LLtoXYZ is defined in clipping/geometry.h 
  */

  
  /* The values for the corners of the final cell are retrieved for comparison within the loop. */


  for (corner_no = 0; corner_no < ncorner; corner_no++){
    LLtoXYZ(grid_corner_lon[(gridsize - 1) * ncorner + corner_no], grid_corner_lat[(gridsize - 1) * ncorner + corner_no], cell_corner_coordinates);
    
    corners_of_the_final_cell[corner_no * 3 + 0] = cell_corner_coordinates[0];
    corners_of_the_final_cell[corner_no * 3 + 1] = cell_corner_coordinates[1];
    corners_of_the_final_cell[corner_no * 3 + 2] = cell_corner_coordinates[2];
  }
  

  for (cell_no = 0; cell_no < gridsize; cell_no++)
    {
      printf("\n\n");
      printf("Cell Number %d:\n\n", cell_no + 1);
      printf("Euclidean coordinates are:\n\n");
      
      LLtoXYZ(grid_center_lon[cell_no], grid_center_lat[cell_no], center_point_in_Euclidean_space);
      
      for (corner_no = 0; corner_no < ncorner; corner_no++)
	{
	  LLtoXYZ(grid_corner_lon[cell_no * ncorner + corner_no], grid_corner_lat[cell_no * ncorner + corner_no], cell_corner_coordinates);
	 
	  /* The components of the result vector are appended to the list of cell corner coordinates. */
	  
	  for (vector_component = 0; vector_component < 3; vector_component++){	    
	    cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + vector_component] = cell_corner_coordinates[vector_component];	  
	  }
	  printf("(%f, %f, %f)\n", cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + 0], cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + 1], cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + 2]);
	}
      
      printf("\n");
      
      /* 
	 Not all cells have the same number of corners. The array, however, has ncorner * 3  values for each cell, where ncorner is the maximum number of corners. Unused values have been filled with the values of the final cell.
	 The following identifies the surplus corners and gives the correct length of the cell.
      */
      
      int actual_number_of_corners = ncorner;
      
      
      for (corner_no = ncorner - 1; corner_no >= 0; corner_no--){
	if (cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + 0] == corners_of_the_final_cell[corner_no * 3 + 0]){
	  if (cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + 1] == corners_of_the_final_cell[corner_no * 3 + 1]){
	    if (cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + 2] == corners_of_the_final_cell[corner_no * 3 + 2]){
	      actual_number_of_corners = actual_number_of_corners - 1;
	    }
	  }
	} else {
	  break;
	}	
      }            
      
      no_cells_with_a_specific_no_of_corners[actual_number_of_corners] = no_cells_with_a_specific_no_of_corners[actual_number_of_corners] + 1;
      
      
      /* Checks if there are any duplicate vertices in the list of corners. */
      
      /* Note that the last (additional) corner has not been set yet. */

      if (no_of_duplicates_in_this_list_of_vertices(cell_corners_in_Euclidean_space_open_cell, actual_number_of_corners * 3) > 0){
	no_of_cells_with_duplicates += 1;
      }

      /* We are creating a closed polygon/cell by setting the additional last corner to be the same as the first one. */

      for (corner_no = 0; corner_no < actual_number_of_corners; corner_no++){
	cell_corners_in_Euclidean_space[corner_no * 3 + 0] = cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + 0];
	cell_corners_in_Euclidean_space[corner_no * 3 + 1] = cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + 1];
	cell_corners_in_Euclidean_space[corner_no * 3 + 2] = cell_corners_in_Euclidean_space_open_cell[corner_no * 3 + 2];
      }

      cell_corners_in_Euclidean_space[(ncorner * 3) + 0] = cell_corners_in_Euclidean_space[0];
      cell_corners_in_Euclidean_space[(ncorner * 3) + 1] = cell_corners_in_Euclidean_space[1];
      cell_corners_in_Euclidean_space[(ncorner * 3) + 2] = cell_corners_in_Euclidean_space[2];

      /* The area of this cell is calculated. */
            
      double cell_area = inefficiently_calculate_the_area_of_a_polygon_in_Euclidean_space(ncorner, cell_corners_in_Euclidean_space);

      printf("The cell area is: %f\n\n", cell_area);
      
      /* Checking whether the cell has a size greater zero. */

      if (cell_area == 0){
	printf("This cell is degenerate. Its area equals zero. The cell corners are colinear.\n\n");
	no_of_degenerate_cells += 1;	
	continue;	 
      }
      
      /* Checking whether the cell corners/polygon vertices are arranged in a clockwise or counterclockwise order. This is done by looking at the sign of the cell area. */

      if (are_polygon_vertices_arranged_in_clockwise_order(cell_area) == 1){
	printf("The cell's corners are arranged in clockwise order.\n\n");
	no_of_cells_with_vertices_arranged_in_clockwise_order += 1;
      } else {
	printf("The cell's corners are arranged in counterclockwise order.\n\n");
	no_of_cells_with_vertices_arranged_in_counterclockwise_order += 1;
      }
      
      continue;

      /*
	Each cell corresponds to a two-dimensional polygon now in unknown orientation in three-dimensional space. Each cell and its center point are coplanar. THIS IS A GIVEN !!!
	In order to solve the two-dimensional point-in-polygon problem for each cell, the three-dimensional coordinates of the polygon and its center point are projected onto the two-dimensional plane they form.
        
	This is done in the following steps:

	1) Compute two vectors that lie on the cell plane. The first three corners of the cell are used to do this. THIS MEANS THAT AT LEAST THREE CORNERS MUST BE GIVEN.
	   Then compute the normal of the cell plane the two vectors are on. The normal is the new z-axis.
	2) Compute the new y-axis by computing the cross product of the new z-axis and the old x-axis.
	3) Compute the new x-axis by computing the cross product of the new z-axis and the new y-axis.
	4) Project every corner point onto the new x- and y-axes by using the dot product.

	The result is a xy tuple for each corner and the presumend center point which is a projection onto the plane the xyz corner points form.
      */

      
      struct axis new_z_axis = compute_the_new_z_axis(cell_corners_in_Euclidean_space);

      /* Regardless of whether the vertices are arranged in clockwise or counterclockwise order, the z-axis is to point inward. */
      
      double vector_pointing_from_origin_to_first_cell_vertex[3] = {cell_corners_in_Euclidean_space[0], cell_corners_in_Euclidean_space[1], cell_corners_in_Euclidean_space[2]};
      
      int axis_points_inward = 0;
      
      if (sign(dotproduct(new_z_axis.axis_vector, vector_pointing_from_origin_to_first_cell_vertex)) < 0){
	axis_points_inward = 1;
	axes_pointing_inward += 1;
      }
      
      if (axis_points_inward == 0){
	new_z_axis.axis_vector[0] = new_z_axis.axis_vector[0] * (-1);
	new_z_axis.axis_vector[1] = new_z_axis.axis_vector[1] * (-1);
	new_z_axis.axis_vector[2] = new_z_axis.axis_vector[2] * (-1);
      }

      struct axis new_y_axis = compute_the_new_y_axis(new_z_axis);      
      struct axis new_x_axis = compute_the_new_x_axis(new_z_axis, new_y_axis);
      
      int ortho_zy = dotproduct(new_z_axis.axis_vector, new_y_axis.axis_vector);
      int ortho_yx = dotproduct(new_y_axis.axis_vector, new_x_axis.axis_vector);
      int ortho_zx = dotproduct(new_z_axis.axis_vector, new_x_axis.axis_vector);

      if ((ortho_zy + ortho_yx + ortho_zx) == 0){
	orthogonal_axes = orthogonal_axes + 1;
      } 
      
      transform_Euclidean_corner_coordinates_onto_the_cell_plane(new_x_axis, new_y_axis, new_z_axis, &cell_corners_in_Euclidean_space, &cell_corners_on_cell_plane, ncorner);
      transform_Euclidean_center_coordinates_onto_the_cell_plane(new_x_axis, new_y_axis, new_z_axis, &center_point_in_Euclidean_space, &center_point_on_cell_plane);

      
  
      /* Convexity of the cell is tested. */

      printf("Checking if the cell is convex ... ");

      if (is_simple_polygon_convex(cell_corners_on_cell_plane, ncorner) == 1){
	printf("cell is convex!\n\n");
	no_of_convex_cells += 1;
      } else {
	printf("cell is not convex!\n\n");
	no_of_degenerate_cells += 1;
      }
    
      /* The winding numbers algorithm is used to test whether the presumed center point is within the bounds of the cell. */
        
      int winding_number = winding_numbers_algorithm_without_printfs(cell_corners_on_cell_plane, ncorner, center_point_on_cell_plane);

      if (winding_number == 0){
	printf("The presumed center point lies OUTSIDE the bounds of the cell.\n\n\n\n");
      } else {
	printf("The presumed center point lies INSIDE the bounds of the cell.\n\n\n\n");
	no_of_cells_with_center_points_within_their_bounds += 1;
      }

    }

  

  printf("\n\n");

  printf("The corners of the final cell are:\n\n");
  for (corner_no = 0; corner_no < ncorner; corner_no++){
    printf("(%f, %f, %f)\n", corners_of_the_final_cell[corner_no * 3 + 0], corners_of_the_final_cell[corner_no * 3 + 1], corners_of_the_final_cell[corner_no * 3 + 2]);
  }

  printf("\nNumber of cells with a certain number of corners:\n\n");

  for(int i = 0; i < ncorner + 1; i++){
    printf("%d corners in %f cells.\n", i, no_cells_with_a_specific_no_of_corners[i]);
  }
  

  printf("\n\n");
  
  printf("There are %u cells in all.\n\n", gridsize );
  
  printf("Number of cells with duplicate vertices: %u\n", no_of_cells_with_duplicates);
  printf("Number of cells with clockwise ordering of vertices: %u\n", no_of_cells_with_vertices_arranged_in_clockwise_order);
  printf("Number of cells with counterclockwise ordering of vertices: %u\n", no_of_cells_with_vertices_arranged_in_counterclockwise_order);
  printf("Number of convex cells: %u\n", no_of_convex_cells);
  printf("Number of nonsimple or degenerate cells: %u\n", no_of_degenerate_cells);
  printf("Number of cells with presumed center points within their bounds: %u\n\n", no_of_cells_with_center_points_within_their_bounds);

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
