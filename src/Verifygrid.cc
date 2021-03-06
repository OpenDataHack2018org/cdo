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

/*
   This module contains the following operators:
*/

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "pstream_int.h"
#include "time.h"
extern "C" {
#include "lib/yac/geometry.h"
}

static void
quick_sort(double *array, size_t array_length)
{
  size_t i, j;
  double temp;

  if (array_length < 2) return;
  double p = array[array_length / 2];
  for (i = 0, j = array_length - 1;; i++, j--)
    {
      while (array[i] < p) i++;
      while (p < array[j]) j--;
      if (i >= j) break;
      temp = array[i];
      array[i] = array[j];
      array[j] = temp;
    }
  quick_sort(array, i);
  quick_sort(array + i, array_length - i);
}

/* Quicksort is called with a pointer to the array of center points to be sorted and an integer indicating its length.
 * It sorts the array by its longitude coordinates */
static void
quick_sort_by_lon(double *array, size_t array_length)
{
  if (array_length < 4) return;

  double p = ((array_length / 2) % 2) ? array[(array_length / 2) + 1] : array[array_length / 2];

  double temp_lon, temp_lat;
  size_t i, j;
  for (i = 0, j = array_length - 2;; i += 2, j -= 2)
    {
      while (array[i] < p) i += 2;

      while (p < array[j]) j -= 2;

      if (i >= j) break;

      temp_lon = array[i];
      temp_lat = array[i + 1];
      array[i] = array[j];
      array[i + 1] = array[j + 1];
      array[j] = temp_lon;
      array[j + 1] = temp_lat;
    }

  quick_sort_by_lon(array, i);
  quick_sort_by_lon(array + i, array_length - i);
}

/* This uses quicksort to sort the latitude coordinates in a subarray of all coordinates. */
static void
quick_sort_of_subarray_by_lat(double *array, size_t subarray_start, size_t subarray_end)
{
  size_t subarray_length = (subarray_end - subarray_start) / 2 + 1;
  std::vector<double> subarray(subarray_length);
  size_t subarray_index = 0;

  for (size_t index = subarray_start + 1; index <= subarray_end + 1; index += 2)
    {
      subarray[subarray_index] = array[index];
      subarray_index += 1;
    }

  quick_sort(subarray.data(), subarray_length);

  subarray_index = 0;

  for (size_t index = subarray_start + 1; index <= subarray_end + 1; index += 2)
    {
      array[index] = subarray[subarray_index];
      subarray_index += 1;
    }
}

static double
determinant(double matrix[3][3])
{
  /* Calculates the determinant for a 3 x 3 matrix. */

  return matrix[0][0] * matrix[1][1] * matrix[2][2] + matrix[0][1] * matrix[1][2] * matrix[2][0]
         + matrix[0][2] * matrix[1][0] * matrix[2][1] - matrix[0][2] * matrix[1][1] * matrix[2][0]
         - matrix[0][1] * matrix[1][0] * matrix[2][2] - matrix[0][0] * matrix[1][2] * matrix[2][1];
}

static void
find_unit_normal(double a[3], double b[3], double c[3], double *unit_normal)
{
  /* Calculates the unit normal for a plane defined on three points a, b, c in Euclidean space. */

  double matrix_for_x[3][3] = { { 1, a[1], a[2] }, { 1, b[1], b[2] }, { 1, c[1], c[2] } };

  double x = determinant(matrix_for_x);

  double matrix_for_y[3][3] = { { a[0], 1, a[2] }, { b[0], 1, b[2] }, { c[0], 1, c[2] } };

  double y = determinant(matrix_for_y);

  double matrix_for_z[3][3] = { { a[0], a[1], 1 }, { b[0], b[1], 1 }, { c[0], c[1], 1 } };

  double z = determinant(matrix_for_z);

  double magnitude = sqrt(x * x + y * y + z * z);

  unit_normal[0] = x / magnitude;
  unit_normal[1] = y / magnitude;
  unit_normal[2] = z / magnitude;
}

int
find_coordinate_to_ignore(double *cell_corners_xyz)
{
  /* Takes the first three corners/vertices of the cell and calculates the unit normal via determinants. */

  double *pcorner_coordinates1 = &cell_corners_xyz[0];
  double *pcorner_coordinates2 = &cell_corners_xyz[3];
  double *pcorner_coordinates3 = &cell_corners_xyz[6];

  double surface_normal_of_the_cell[3];
  find_unit_normal(pcorner_coordinates1, pcorner_coordinates2, pcorner_coordinates3, surface_normal_of_the_cell);

  /* The surface normal is used to choose the coordinate to ignore. */

  double abs_x = fabs(surface_normal_of_the_cell[0]);
  double abs_y = fabs(surface_normal_of_the_cell[1]);
  double abs_z = fabs(surface_normal_of_the_cell[2]);

  int coordinate_to_ignore = 3;

  if (abs_x > abs_y)
    {
      if (abs_x > abs_z) coordinate_to_ignore = 1;
    }
  else
    {
      if (abs_y > abs_z) coordinate_to_ignore = 2;
    }

  return coordinate_to_ignore;
}

static double
is_point_left_of_edge(double point_on_line_1[2], double point_on_line_2[2], double point[2])
{
  /*
     Computes whether a point is left of the line through point_on_line_1 and point_on_line_2. This is part of the
     solution to the point in polygon problem. Returns 0 if the point is on the line, > 0 if the point is left of the
     line, and < 0 if the point is right of the line. This algorithm is by Dan Sunday (geomalgorithms.com) and is
     completely free for use and modification.
  */

  double answer = ((point_on_line_2[0] - point_on_line_1[0]) * (point[1] - point_on_line_1[1])
                   - (point[0] - point_on_line_1[0]) * (point_on_line_2[1] - point_on_line_1[1]));

  return answer;
}

int
winding_numbers_algorithm(double cell_corners[], int number_corners, double point[])
{
  /*
     Computes whether a point is inside the bounds of a cell. This is the solution to the point in polygon problem.
     Returns 0 if the point is outside, returns 1 if the point is inside the cell. Based on an algorithm by Dan Sunday
     (geomalgorithms.com). His algorithm is completely free for use and modification.
  */

  int winding_number = 0;

  for (int i = 0; i < number_corners - 1; i++)
    {
      if (cell_corners[i * 2 + 1] <= point[1])
        {
          if (cell_corners[(i + 1) * 2 + 1] > point[1])
            {
              double point_on_edge_1[2] = { cell_corners[i * 2 + 0], cell_corners[i * 2 + 1] };
              double point_on_edge_2[2] = { cell_corners[(i + 1) * 2 + 0], cell_corners[(i + 1) * 2 + 1] };

              if (is_point_left_of_edge(point_on_edge_1, point_on_edge_2, point) > 0) winding_number++;
            }
        }
      else
        {
          if (cell_corners[(i + 1) * 2 + 1] <= point[1])
            {
              double point_on_edge_1[2] = { cell_corners[i * 2 + 0], cell_corners[i * 2 + 1] };
              double point_on_edge_2[2] = { cell_corners[(i + 1) * 2 + 0], cell_corners[(i + 1) * 2 + 1] };

              if (is_point_left_of_edge(point_on_edge_1, point_on_edge_2, point) < 0) winding_number--;
            }
        }
    }

  return winding_number;
}

static double
sign(double x)
{
  /* Is +1 if x is positive, -1 if x is negative and 0 if x is zero.*/

  return (x > 0) - (x < 0);
}

static bool
is_simple_polygon_convex(double cell_corners[], int number_corners)
{
  /* Tests in which direction the polygon winds when walking along its edges. Does so for all edges of the polygon. */

  double direction = 0;

  for (int i = 0; i < number_corners - 2; i++)
    {
      double turns_to = (cell_corners[i * 2 + 0] - cell_corners[(i + 1) * 2 + 0])
                            * (cell_corners[(i + 1) * 2 + 1] - cell_corners[(i + 2) * 2 + 1])
                        - (cell_corners[i * 2 + 1] - cell_corners[(i + 1) * 2 + 1])
                              * (cell_corners[(i + 1) * 2 + 0] - cell_corners[(i + 2) * 2 + 0]);

      /* In the first iteration the direction of winding of the entire polygon is set. Better not be 0.*/

      if (i == 1) direction = turns_to;

      if (IS_NOT_EQUAL(sign(direction), sign(turns_to)))
        {
          if (IS_NOT_EQUAL(direction, 0)) return false;
        }
      else
        {
          direction = turns_to;
        }
    }

  return true;
}

double
calculate_the_polygon_area(double cell_corners[], int number_corners)
{
  /* This algorithm is based on the calculation from Wolfram Mathworld Polygon Area. It results in the area of planar
   * non-self-intersecting polygon. */

  double twice_the_polygon_area = 0;

  for (int i = 0; i < number_corners - 1; i++)
    twice_the_polygon_area
        += (cell_corners[i * 2 + 0] * cell_corners[(i + 1) * 2 + 1]) - (cell_corners[(i + 1) * 2 + 0] * cell_corners[i * 2 + 1]);

  return twice_the_polygon_area / 2;
}

bool
are_polygon_vertices_arranged_in_clockwise_order(double cell_area)
{
  bool status = false;

  /* A negative area indicates a clockwise arrangement of vertices, a positive area a counterclockwise arrangement.
   * There should be an area to begin with.
   */
  if (cell_area < 0) status = true;

  return status;
}

static void
verify_grid(int gridtype, size_t gridsize, int gridno, int ngrids, int ncorner, double *grid_center_lon, double *grid_center_lat,
            double *grid_corner_lon, double *grid_corner_lat)
{
  /*
     First, this function performs the following test:

     1) it tests whether there are duplicate cells in the given grid by
     comparing their center points

     Additionally, on each cell of a given grid:

     2) it tests whether all cells are convex and all cell bounds have the same
     orientation, i.e. the corners of the cell are in clockwise or
     counterclockwise order

     3) it tests whether the center point is within the bounds of the cell

     The results of the tests are printed on stdout.
  */

  constexpr double eps = 0.000000001;
  double center_point_xyz[3];
  std::vector<double> cell_corners_xyz_open_cell(3 * ncorner);

  double corner_coordinates[3];
  double center_point_plane_projection[2];

  size_t no_of_cells_with_duplicates = 0;
  size_t no_usable_cells = 0;
  size_t no_convex_cells = 0;
  size_t no_clockwise_cells = 0;
  size_t no_counterclockwise_cells = 0;
  size_t no_of_cells_with_center_points_out_of_bounds = 0;
  size_t no_unique_center_points = 1;

  std::vector<int> no_cells_with_a_specific_no_of_corners(ncorner);

  for (int i = 0; i < ncorner; i++) no_cells_with_a_specific_no_of_corners[i] = 0;

  if (ngrids == 1)
    cdoPrintBlue("Grid consists of %zu cells (type: %s), of which", gridsize, gridNamePtr(gridtype));
  else
    cdoPrintBlue("Grid no %u (of %u) consists of %zu cells (type: %s), of which", gridno + 1, ngrids, gridsize,
                 gridNamePtr(gridtype));
  // cdoPrint("");

  /* For performing the first test, an array of all center point coordinates is built. */

  double *center_point_array = (double *) Malloc(gridsize * 2 * sizeof(double));

  for (size_t cell_no = 0; cell_no < gridsize; cell_no++)
    {
      center_point_array[cell_no * 2 + 0] = grid_center_lon[cell_no];
      center_point_array[cell_no * 2 + 1] = grid_center_lat[cell_no];
    }

  /* The cell center points are sorted by their first coordinate (lon) with quicksort. */

  quick_sort_by_lon(center_point_array, gridsize * 2);

  /* Now the lat coordinates in subarrays that reflect equal lon coordinates are sorted with quicksort. */

  int subarray_start = 0;
  int subarray_end = 0;

  for (size_t cell_no = 0; cell_no < gridsize - 1; cell_no++)
    {
      if (cell_no == gridsize - 2)
        {
          subarray_end = gridsize * 2 - 2;
          quick_sort_of_subarray_by_lat(center_point_array, subarray_start, subarray_end);
        }

      if (fabs(center_point_array[cell_no * 2 + 0] - center_point_array[(cell_no + 1) * 2 + 0]) > eps)
        {
          subarray_end = cell_no * 2;
          if ((subarray_end - subarray_start) > 1) quick_sort_of_subarray_by_lat(center_point_array, subarray_start, subarray_end);

          subarray_start = subarray_end + 2;
        }
    }

  /* Now checking for the number of unique center point coordinates. */

  for (size_t cell_no = 0; cell_no < gridsize - 1; cell_no++)
    {
      if (fabs(center_point_array[cell_no * 2 + 0] - center_point_array[(cell_no + 1) * 2 + 0]) < eps)
        {
          if (fabs(center_point_array[cell_no * 2 + 1] - center_point_array[(cell_no + 1) * 2 + 1]) < eps)
            continue;
          else
            no_unique_center_points += 1;
        }
      else
        {
          no_unique_center_points += 1;
        }
    }

  Free(center_point_array);

  // used only actual_number_of_corners
  int *marked_duplicate_indices = (int *) Malloc(ncorner * sizeof(int));
  double *cell_corners_xyz_without_duplicates = (double *) Malloc(3 * ncorner * sizeof(double));
  double *cell_corners_xyz = (double *) Malloc(3 * (ncorner + 1) * sizeof(double));
  double *cell_corners_plane_projection = (double *) Malloc(2 * (ncorner + 1) * sizeof(double));

  /*
     Latitude and longitude are spherical coordinates on a unit circle. Each such coordinate tuple is transformed into a
     triple of Cartesian coordinates in Euclidean space. This is first done for the presumed center point of the cell
     and then for all the corners of the cell. LLtoXYZ_deg is defined in clipping/geometry.h
  */

  for (size_t cell_no = 0; cell_no < gridsize; cell_no++)
    {
      /* Conversion of center point spherical coordinates to Cartesian coordinates. */

      LLtoXYZ_deg(grid_center_lon[cell_no], grid_center_lat[cell_no], center_point_xyz);

      for (int corner_no = 0; corner_no < ncorner; corner_no++)
        {
          /* Conversion of corner spherical coordinates to Cartesian coordinates. */

          LLtoXYZ_deg(grid_corner_lon[cell_no * ncorner + corner_no], grid_corner_lat[cell_no * ncorner + corner_no],
                      corner_coordinates);

          /* The components of the result vector are appended to the list of cell corner coordinates. */

          int off = corner_no * 3;
          cell_corners_xyz_open_cell[off + 0] = corner_coordinates[0];
          cell_corners_xyz_open_cell[off + 1] = corner_coordinates[1];
          cell_corners_xyz_open_cell[off + 2] = corner_coordinates[2];
        }

      /*
         Not all cells have the same number of corners. The array, however, has ncorner * 3  values for each cell, where
         ncorner is the maximum number of corners. Unused values have been filled with the values of the final cell. The
         following identifies the surplus corners and gives the correct length of the cell.
      */

      int actual_number_of_corners = ncorner;

      for (int corner_no = ncorner - 1; corner_no > 0; corner_no--)
        {
          int off = corner_no * 3;
          int off2 = (corner_no - 1) * 3;
          if (IS_EQUAL(cell_corners_xyz_open_cell[off + 0], cell_corners_xyz_open_cell[off2 + 0])
              && IS_EQUAL(cell_corners_xyz_open_cell[off + 1], cell_corners_xyz_open_cell[off2 + 1])
              && IS_EQUAL(cell_corners_xyz_open_cell[off + 2], cell_corners_xyz_open_cell[off2 + 2]))
            actual_number_of_corners = actual_number_of_corners - 1;
          else
            break;
        }

      no_cells_with_a_specific_no_of_corners[actual_number_of_corners - 1] += 1;

      /* If there are less than three corners in the cell, it is unusable and considered degenerate. No area can be
       * computed. */

      if (actual_number_of_corners < 3)
        {
          if (cdoVerbose)
            fprintf(stdout, "Less than three vertices found in cell no %zu. This cell is considered degenerate and "
                            "will be omitted from further computation!\n",
                    cell_no + 1);

          continue;
        }

      no_usable_cells++;

      /* Checks if there are any duplicate vertices in the list of corners. Note that the last (additional) corner has
       * not been set yet. */

      for (int i = 0; i < actual_number_of_corners; i++) marked_duplicate_indices[i] = 0;

      int no_duplicates = 0;

      for (int i = 0; i < actual_number_of_corners * 3; i = i + 3)
        for (int j = i + 3; j < actual_number_of_corners * 3; j = j + 3)
          if (fabs(cell_corners_xyz_open_cell[i + 0] - cell_corners_xyz_open_cell[j + 0]) < eps
              && fabs(cell_corners_xyz_open_cell[i + 1] - cell_corners_xyz_open_cell[j + 1]) < eps
              && fabs(cell_corners_xyz_open_cell[i + 2] - cell_corners_xyz_open_cell[j + 2]) < eps)
            {
              if (cdoVerbose)
                fprintf(stdout, "The duplicate vertex %f, %f, %f was found in cell no %zu.\n", cell_corners_xyz_open_cell[j],
                        cell_corners_xyz_open_cell[j + 1], cell_corners_xyz_open_cell[j + 2], cell_no + 1);

              no_duplicates += 1;
              marked_duplicate_indices[j / 3] = 1;
            }

      /* Writes the unique corner vertices in a new array. */

      int unique_corner_number = 0;

      for (int corner_no = 0; corner_no < actual_number_of_corners; corner_no++)
        {
          if (marked_duplicate_indices[corner_no] == 0)
            {
              int off = corner_no * 3;
              int off2 = unique_corner_number * 3;
              cell_corners_xyz_without_duplicates[off2 + 0] = cell_corners_xyz_open_cell[off + 0];
              cell_corners_xyz_without_duplicates[off2 + 1] = cell_corners_xyz_open_cell[off + 1];
              cell_corners_xyz_without_duplicates[off2 + 2] = cell_corners_xyz_open_cell[off + 2];
              unique_corner_number += 1;
            }
        }

      actual_number_of_corners = actual_number_of_corners - no_duplicates;

      if (no_duplicates != 0) no_of_cells_with_duplicates += 1;

      /* If there are less than three corners in the cell left after removing duplicates, it is unusable and considered
       * degenerate. No area can be computed. */

      if (actual_number_of_corners < 3)
        {
          if (cdoVerbose)
            fprintf(stdout, "Less than three vertices found in cell no %zu. This cell is considered degenerate and "
                            "will be omitted from further computation!\n",
                    cell_no + 1);

          continue;
        }

      /* We are creating a closed polygon/cell by setting the additional last corner to be the same as the first one. */

      for (int corner_no = 0; corner_no < actual_number_of_corners; corner_no++)
        {
          int off = corner_no * 3;
          cell_corners_xyz[off + 0] = cell_corners_xyz_without_duplicates[off + 0];
          cell_corners_xyz[off + 1] = cell_corners_xyz_without_duplicates[off + 1];
          cell_corners_xyz[off + 2] = cell_corners_xyz_without_duplicates[off + 2];
        }

      cell_corners_xyz[actual_number_of_corners * 3 + 0] = cell_corners_xyz[0];
      cell_corners_xyz[actual_number_of_corners * 3 + 1] = cell_corners_xyz[1];
      cell_corners_xyz[actual_number_of_corners * 3 + 2] = cell_corners_xyz[2];

      int coordinate_to_ignore = find_coordinate_to_ignore(cell_corners_xyz);

      /* The remaining two-dimensional coordinates are extracted into one array for all the cell's corners and into one
       * array for the center point. */

      /* The following projection on the plane that two coordinate axes lie on changes the arrangement of the polygon
         vertices if the coordinate to be ignored along the third axis is smaller than 0. In this case, the result of
         the computation of the orientation of vertices needs to be  inverted. Clockwise becomes counterclockwise and
         vice versa. */

      bool invert_result = false;
      if (cell_corners_xyz[coordinate_to_ignore - 1] < 0) invert_result = true;

      switch (coordinate_to_ignore)
        {
        case 1:
          for (int corner_no = 0; corner_no <= actual_number_of_corners; corner_no++)
            {
              cell_corners_plane_projection[corner_no * 2 + 0] = cell_corners_xyz[corner_no * 3 + 1];
              cell_corners_plane_projection[corner_no * 2 + 1] = cell_corners_xyz[corner_no * 3 + 2];
            }
          center_point_plane_projection[0] = center_point_xyz[1];
          center_point_plane_projection[1] = center_point_xyz[2];
          break;
        case 2:
          for (int corner_no = 0; corner_no <= actual_number_of_corners; corner_no++)
            {
              cell_corners_plane_projection[corner_no * 2 + 0] = cell_corners_xyz[corner_no * 3 + 2];
              cell_corners_plane_projection[corner_no * 2 + 1] = cell_corners_xyz[corner_no * 3 + 0];
            }
          center_point_plane_projection[0] = center_point_xyz[2];
          center_point_plane_projection[1] = center_point_xyz[0];
          break;
        case 3:
          for (int corner_no = 0; corner_no <= actual_number_of_corners; corner_no++)
            {
              cell_corners_plane_projection[corner_no * 2 + 0] = cell_corners_xyz[corner_no * 3 + 0];
              cell_corners_plane_projection[corner_no * 2 + 1] = cell_corners_xyz[corner_no * 3 + 1];
            }
          center_point_plane_projection[0] = center_point_xyz[0];
          center_point_plane_projection[1] = center_point_xyz[1];
          break;
        }

      /* Checking for convexity of the cell. */

      if (is_simple_polygon_convex(cell_corners_plane_projection, actual_number_of_corners + 1)) no_convex_cells += 1;

      /* Checking the arrangement or direction of cell vertices. */

      long double polygon_area = calculate_the_polygon_area(cell_corners_plane_projection, actual_number_of_corners + 1);
      bool is_clockwise = are_polygon_vertices_arranged_in_clockwise_order(polygon_area);

      /* If the direction of the vertices was flipped during the projection onto the two-dimensional plane, the previous
       * result needs to be inverted now. */

      if (invert_result) is_clockwise = !is_clockwise;

      /* The overall counter of (counter)clockwise cells is increased by one. */

      if (is_clockwise)
        no_clockwise_cells += 1;
      else
        no_counterclockwise_cells += 1;

      /* The winding numbers algorithm is used to test whether the presumed center point is within the bounds of the
       * cell. */

      int winding_number
          = winding_numbers_algorithm(cell_corners_plane_projection, actual_number_of_corners + 1, center_point_plane_projection);

      // if ( winding_number == 0 ) printf("%d,", cell_no+1);
      if (winding_number == 0) no_of_cells_with_center_points_out_of_bounds += 1;

      if (cdoVerbose && winding_number == 0)
        {
          printf("cell_no %zu: ", cell_no + 1);
          printf(" lon=%g lat=%g : ", grid_center_lon[cell_no], grid_center_lat[cell_no]);
          for (int corner_no = 0; corner_no < ncorner; corner_no++)
            printf(" %g/%g ", grid_corner_lon[cell_no * ncorner + corner_no], grid_corner_lat[cell_no * ncorner + corner_no]);
          printf("\n");
        }
    }

  Free(marked_duplicate_indices);
  Free(cell_corners_plane_projection);
  Free(cell_corners_xyz);
  Free(cell_corners_xyz_without_duplicates);

  size_t no_nonunique_cells = gridsize - no_unique_center_points;
  size_t no_nonconvex_cells = gridsize - no_convex_cells;
  size_t no_nonusable_cells = gridsize - no_usable_cells;

  for (int i = 2; i < ncorner; i++)
    if (no_cells_with_a_specific_no_of_corners[i])
      cdoPrintBlue("%9d cells have %d vertices", no_cells_with_a_specific_no_of_corners[i], i + 1);

  if (no_of_cells_with_duplicates) cdoPrintBlue("%9zu cells have duplicate vertices", no_of_cells_with_duplicates);

  if (no_nonusable_cells) cdoPrintRed("%9zu cells have unusable vertices", no_nonusable_cells);

  if (no_nonunique_cells) cdoPrintRed("%9zu cells are not unique", no_nonunique_cells);

  if (no_nonconvex_cells) cdoPrintRed("%9zu cells are non-convex", no_nonconvex_cells);

  if (no_clockwise_cells) cdoPrintRed("%9zu cells have their vertices arranged in a clockwise order", no_clockwise_cells);

  if (no_of_cells_with_center_points_out_of_bounds)
    cdoPrintRed("%9zu cells have their center points located outside their boundaries",
                no_of_cells_with_center_points_out_of_bounds);

  // cdoPrint("");
}

void *
Verifygrid(void *argument)
{
  char units[CDI_MAX_NAME];

  cdoInitialize(argument);

  int VERIFYGRID = cdoOperatorAdd("verifygrid", 0, 0, NULL);

  int operatorID = cdoOperatorID();

  int streamID = cdoStreamOpenRead(cdoStreamName(0));

  int vlistID = pstreamInqVlist(streamID);

  int ngrids = vlistNgrids(vlistID);
  for (int gridno = 0; gridno < ngrids; ++gridno)
    {
      bool lgeo = true;
      bool luse_grid_corner = true;

      int gridID = vlistGrid(vlistID, gridno);
      int gridtype = gridInqType(gridID);

      if (gridtype == GRID_GME) gridID = gridToUnstructured(gridID, 1);

      if (gridtype != GRID_UNSTRUCTURED && gridtype != GRID_CURVILINEAR)
        {
          if (gridtype == GRID_GENERIC || gridtype == GRID_SPECTRAL)
            {
              lgeo = false;
            }
          else
            {
              gridID = gridToCurvilinear(gridID, 1);
            }
        }

      size_t gridsize = gridInqSize(gridID);
      /*
        if ( gridInqMaskGME(gridID, NULL) )
        {
        int *grid_mask = (int*) Malloc(gridsize*sizeof(int));
        gridInqMaskGME(gridID, grid_mask);
        Free(grid_mask);
        }
      */
      if (lgeo)
        {
          std::vector<double> grid_corner_lat, grid_corner_lon;
          int ncorner = 4;
          if (gridInqType(gridID) == GRID_UNSTRUCTURED) ncorner = gridInqNvertex(gridID);

          std::vector<double> grid_center_lat(gridsize);
          std::vector<double> grid_center_lon(gridsize);

          gridInqYvals(gridID, grid_center_lat.data());
          gridInqXvals(gridID, grid_center_lon.data());

          /* Convert lat/lon units if required */
          gridInqXunits(gridID, units);
          grid_to_degree(units, gridsize, grid_center_lon.data(), "grid center lon");
          gridInqYunits(gridID, units);
          grid_to_degree(units, gridsize, grid_center_lat.data(), "grid center lat");

          if (luse_grid_corner)
            {
              if (ncorner == 0) cdoAbort("grid corner missing!");
              size_t nalloc = ncorner * gridsize;
              grid_corner_lat.resize(nalloc);
              grid_corner_lon.resize(nalloc);

              if (gridInqYbounds(gridID, NULL) && gridInqXbounds(gridID, NULL))
                {
                  gridInqYbounds(gridID, grid_corner_lat.data());
                  gridInqXbounds(gridID, grid_corner_lon.data());
                }
              else
                {
                  cdoAbort("Grid corner missing!");
                }

              /* Note: using units from latitude instead from bounds */
              grid_to_degree(units, ncorner * gridsize, grid_corner_lon.data(), "grid corner lon");
              grid_to_degree(units, ncorner * gridsize, grid_corner_lat.data(), "grid corner lat");
            }

          if (operatorID == VERIFYGRID)
            verify_grid(gridtype, gridsize, gridno, ngrids, ncorner, grid_center_lon.data(), grid_center_lat.data(),
                        grid_corner_lon.data(), grid_corner_lat.data());
        }
      else
        {
          if (ngrids == 1)
            cdoPrintBlue("Grid consists of %zu points (type: %s)", gridsize, gridNamePtr(gridtype));
          else
            cdoPrintBlue("Grid no %u (of %u) consists of %zu points (type: %s)", gridno + 1, ngrids, gridsize,
                         gridNamePtr(gridtype));
          // cdoPrint("");
        }
    }

  pstreamClose(streamID);

  cdoFinish();

  return 0;
}
