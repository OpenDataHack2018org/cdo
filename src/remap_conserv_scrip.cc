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

#include "cdo_int.h"
#include "cdo_wtime.h"
#include "grid.h"
#include "remap.h"
#include "remap_store_link_cnsrv.h"
#include "cdoOptions.h"

#define ZERO 0.0
#define ONE 1.0
#define TWO 2.0
#define THREE 3.0
#define HALF 0.5
#define QUART 0.25

#define BABY_STEP 0.001 /* original value */

/* static double north_thresh =  1.45;  */ /* threshold for coord transformation */
/* static double south_thresh = -2.00;  */ /* threshold for coord transformation */
static double north_thresh = 2.00;         /* threshold for coord transformation */
static double south_thresh = -2.00;        /* threshold for coord transformation */

/* threshold for coord transformation */
void
remap_set_threshhold(double threshhold)
{
  north_thresh = threshhold;
  south_thresh = -threshhold;
}

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/*                                                                         */
/*      CONSERVATIVE INTERPOLATION                                         */
/*                                                                         */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

/*
    This routine is identical to the intersection routine except
    that a coordinate transformation (using a Lambert azimuthal
    equivalent projection) is performed to treat polar cells more
    accurately.
*/
static void
pole_intersection(long *location, double *intrsct_lat, double *intrsct_lon, bool *lcoinc, bool *lthresh, double beglat,
                  double beglon, double endlat, double endlon, double *begseg, bool lrevers, long num_srch_cells, long srch_corners,
                  const size_t *restrict srch_add, const double *restrict srch_corner_lat, const double *restrict srch_corner_lon,
                  bool *luse_last, double *intrsct_x, double *intrsct_y, int *avoid_pole_count, double *avoid_pole_offset)
{
  /*
    Intent(in):
    double beglat, beglon,  ! beginning lat/lon endpoints for segment
    double endlat, endlon   ! ending    lat/lon endpoints for segment
    bool   lrevers          ! flag true if segment integrated in reverse

    Intent(inout) :
    double begseg[2] ! begin lat/lon of full segment
    int *location    ! address in destination array containing this
                     ! segment -- also may contain last location on entry
    bool *lthresh    ! flag segment crossing threshold boundary

    intent(out):
    *int lcoinc      ! flag segment coincident with grid line
    double *intrsct_lat, *intrsct_lon ! lat/lon coords of next intersect.
  */
  /* Local variables */
  long n, next_n, cell;
  long ioffset;
  int loutside; /* flags points outside grid */

  double pi4, rns;                           /*  north/south conversion */
  double x1, x2;                             /*  local x variables for segment */
  double y1, y2;                             /*  local y variables for segment */
  double begx, begy;                         /*  beginning x,y variables for segment */
  double endx, endy;                         /*  beginning x,y variables for segment */
  double begsegx, begsegy;                   /*  beginning x,y variables for segment */
  double grdx1, grdx2;                       /*  local x variables for grid cell */
  double grdy1, grdy2;                       /*  local y variables for grid cell */
  double vec1_y, vec1_x;                     /*  vectors and cross products used */
  double vec2_y, vec2_x;                     /*  during grid search */
  double cross_product, eps;                 /*  eps=small offset away from intersect */
  double s1, s2, determ;                     /*  variables used for linear solve to   */
  double mat1, mat2, mat3, mat4, rhs1, rhs2; /* find intersection */

  double *srch_corner_x; /*  x of each corner of srch cells */
  double *srch_corner_y; /*  y of each corner of srch cells */

  /*printf("pole_intersection: %g %g %g %g\n", beglat, beglon, endlat,
   * endlon);*/

  /* Initialize defaults, flags, etc. */

  if (!*lthresh) *location = -1;
  *lcoinc = false;
  *intrsct_lat = endlat;
  *intrsct_lon = endlon;

  loutside = FALSE;
  s1 = ZERO;

  /* Convert coordinates */

  srch_corner_x = (double *) Malloc(srch_corners * num_srch_cells * sizeof(double));
  srch_corner_y = (double *) Malloc(srch_corners * num_srch_cells * sizeof(double));

  if (beglat > ZERO)
    {
      pi4 = QUART * PI;
      rns = ONE;
    }
  else
    {
      pi4 = -QUART * PI;
      rns = -ONE;
    }

  if (*luse_last)
    {
      x1 = *intrsct_x;
      y1 = *intrsct_y;
    }
  else
    {
      x1 = rns * TWO * sin(pi4 - HALF * beglat) * cos(beglon);
      y1 = TWO * sin(pi4 - HALF * beglat) * sin(beglon);
      *luse_last = true;
    }

  x2 = rns * TWO * sin(pi4 - HALF * endlat) * cos(endlon);
  y2 = TWO * sin(pi4 - HALF * endlat) * sin(endlon);

  for (n = 0; n < srch_corners * num_srch_cells; ++n)
    {
      srch_corner_x[n] = rns * TWO * sin(pi4 - HALF * srch_corner_lat[n]) * cos(srch_corner_lon[n]);
      srch_corner_y[n] = TWO * sin(pi4 - HALF * srch_corner_lat[n]) * sin(srch_corner_lon[n]);
    }

  begx = x1;
  begy = y1;
  endx = x2;
  endy = y2;
  begsegx = rns * TWO * sin(pi4 - HALF * begseg[0]) * cos(begseg[1]);
  begsegy = TWO * sin(pi4 - HALF * begseg[0]) * sin(begseg[1]);
  *intrsct_x = endx;
  *intrsct_y = endy;

  /*
     Search for location of this segment in ocean grid using cross
     product method to determine whether a point is enclosed by a cell
  */
  while (TRUE) /* srch_loop */
    {
      /* If last segment crossed threshold, use that location */

      if (*lthresh)
        {
          for (cell = 0; cell < num_srch_cells; ++cell)
            if (srch_add[cell] == (size_t) *location)
              {
                eps = TINY;
                goto after_srch_loop;
              }
        }

      /* Otherwise normal search algorithm */

      for (cell = 0; cell < num_srch_cells; ++cell) /* cell_loop  */
        {
          ioffset = cell * srch_corners;
          for (n = 0; n < srch_corners; ++n) /* corner_loop */
            {
              next_n = (n + 1) % srch_corners;
              /*
                Here we take the cross product of the vector making
                up each cell side with the vector formed by the vertex
                and search point.  If all the cross products are
                positive, the point is contained in the cell.
              */
              vec1_x = srch_corner_x[ioffset + next_n] - srch_corner_x[ioffset + n];
              vec1_y = srch_corner_y[ioffset + next_n] - srch_corner_y[ioffset + n];
              vec2_x = x1 - srch_corner_x[ioffset + n];
              vec2_y = y1 - srch_corner_y[ioffset + n];

              /* If endpoint coincident with vertex, offset the endpoint */

              if (IS_EQUAL(vec2_x, 0) && IS_EQUAL(vec2_y, 0))
                {
                  x1 += 1.e-10 * (x2 - x1);
                  y1 += 1.e-10 * (y2 - y1);
                  vec2_x = x1 - srch_corner_x[ioffset + n];
                  vec2_y = y1 - srch_corner_y[ioffset + n];
                }

              cross_product = vec1_x * vec2_y - vec2_x * vec1_y;

              /*
                If the cross product for a side is ZERO, the point
                  lies exactly on the side or the length of a side
                  is ZERO.  If the length is ZERO set det > 0.
                  otherwise, perform another cross
                  product between the side and the segment itself.
                If this cross product is also ZERO, the line is
                  coincident with the cell boundary - perform the
                  dot product and only choose the cell if the dot
                  product is positive (parallel vs anti-parallel).
              */
              if (IS_EQUAL(cross_product, 0))
                {
                  if (IS_NOT_EQUAL(vec1_x, 0) || IS_NOT_EQUAL(vec1_y, 0))
                    {
                      vec2_x = x2 - x1;
                      vec2_y = y2 - y1;
                      cross_product = vec1_x * vec2_y - vec2_x * vec1_y;
                    }
                  else
                    cross_product = ONE;

                  if (IS_EQUAL(cross_product, 0))
                    {
                      *lcoinc = true;
                      cross_product = vec1_x * vec2_x + vec1_y * vec2_y;
                      if (lrevers) cross_product = -cross_product;
                    }
                }

              /* If cross product is less than ZERO, this cell doesn't work */

              if (cross_product < ZERO) break; /* corner_loop */

            } /* corner_loop */

          /* If cross products all positive, we found the location */

          if (n >= srch_corners)
            {
              *location = srch_add[cell];
              /*
                If the beginning of this segment was outside the
                grid, invert the segment so the intersection found
                will be the first intersection with the grid
              */
              if (loutside)
                {
                  x2 = begx;
                  y2 = begy;
                  *location = -1;
                  eps = -TINY;
                }
              else
                eps = TINY;

              goto after_srch_loop;
            }

          /* Otherwise move on to next cell */

        } /* cell_loop */

      /*
        If no cell found, the point lies outside the grid.
        take some baby steps along the segment to see if any
        part of the segment lies inside the grid.
      */
      loutside = TRUE;
      s1 = s1 + BABY_STEP;
      x1 = begx + s1 * (x2 - begx);
      y1 = begy + s1 * (y2 - begy);

      /* Reached the end of the segment and still outside the grid return no
       * intersection */

      if (s1 >= ONE)
        {
          Free(srch_corner_y);
          Free(srch_corner_x);
          *luse_last = false;
          return;
        }
    } /* srch_loop */

after_srch_loop:

  /*
    Now that a cell is found, search for the next intersection.
    Loop over sides of the cell to find intersection with side
    must check all sides for coincidences or intersections
  */

  ioffset = cell * srch_corners;

  for (n = 0; n < srch_corners; ++n) /* intrsct_loop */
    {
      next_n = (n + 1) % srch_corners;

      grdy1 = srch_corner_y[ioffset + n];
      grdy2 = srch_corner_y[ioffset + next_n];
      grdx1 = srch_corner_x[ioffset + n];
      grdx2 = srch_corner_x[ioffset + next_n];

      /* Set up linear system to solve for intersection */

      mat1 = x2 - x1;
      mat2 = grdx1 - grdx2;
      mat3 = y2 - y1;
      mat4 = grdy1 - grdy2;
      rhs1 = grdx1 - x1;
      rhs2 = grdy1 - y1;

      determ = mat1 * mat4 - mat2 * mat3;

      /*
         If the determinant is ZERO, the segments are either
           parallel or coincident.  Coincidences were detected
           above so do nothing.
         If the determinant is non-ZERO, solve for the linear
           parameters s for the intersection point on each line
           segment.
         If 0<s1,s2<1 then the segment intersects with this side.
           Return the point of intersection (adding a small
           number so the intersection is off the grid line).
      */
      if (fabs(determ) > 1.e-30)
        {
          s1 = (rhs1 * mat4 - mat2 * rhs2) / determ;
          s2 = (mat1 * rhs2 - rhs1 * mat3) / determ;

          /* Uwe Schulzweida: s1 >= ZERO! (bug fix) */
          if (s2 >= ZERO && s2 <= ONE && s1 >= ZERO && s1 <= ONE)
            {
              /*
                Recompute intersection based on full segment
                so intersections are consistent for both sweeps
              */
              if (!loutside)
                {
                  mat1 = x2 - begsegx;
                  mat3 = y2 - begsegy;
                  rhs1 = grdx1 - begsegx;
                  rhs2 = grdy1 - begsegy;
                }
              else
                {
                  mat1 = x2 - endx;
                  mat3 = y2 - endy;
                  rhs1 = grdx1 - endx;
                  rhs2 = grdy1 - endy;
                }

              determ = mat1 * mat4 - mat2 * mat3;

              /*
                Sometimes due to roundoff, the previous
                determinant is non-ZERO, but the lines
                are actually coincident.  If this is the
                case, skip the rest.
              */
              if (IS_NOT_EQUAL(determ, 0))
                {
                  s1 = (rhs1 * mat4 - mat2 * rhs2) / determ;
                  s2 = (mat1 * rhs2 - rhs1 * mat3) / determ;

                  if (!loutside)
                    {
                      *intrsct_x = begsegx + s1 * mat1;
                      *intrsct_y = begsegy + s1 * mat3;
                    }
                  else
                    {
                      *intrsct_x = endx + s1 * mat1;
                      *intrsct_y = endy + s1 * mat3;
                    }

                  /* Convert back to lat/lon coordinates */

                  *intrsct_lon = rns * atan2(*intrsct_y, *intrsct_x);
                  if (*intrsct_lon < ZERO) *intrsct_lon = *intrsct_lon + PI2;

                  if (fabs(*intrsct_x) > 1.e-10)
                    *intrsct_lat = (pi4 - asin(rns * HALF * (*intrsct_x) / cos(*intrsct_lon))) * TWO;
                  else if (fabs(*intrsct_y) > 1.e-10)
                    *intrsct_lat = (pi4 - asin(HALF * (*intrsct_y) / sin(*intrsct_lon))) * TWO;
                  else
                    *intrsct_lat = TWO * pi4;

                  /* Add offset in transformed space for next pass. */

                  if (s1 - eps / determ < ONE)
                    {
                      *intrsct_x = *intrsct_x - mat1 * (eps / determ);
                      *intrsct_y = *intrsct_y - mat3 * (eps / determ);
                    }
                  else
                    {
                      if (!loutside)
                        {
                          *intrsct_x = endx;
                          *intrsct_y = endy;
                          *intrsct_lat = endlat;
                          *intrsct_lon = endlon;
                        }
                      else
                        {
                          *intrsct_x = begsegx;
                          *intrsct_y = begsegy;
                          *intrsct_lat = begseg[0];
                          *intrsct_lon = begseg[1];
                        }
                    }

                  break; /* intrsct_loop */
                }
            }
        }

      /* No intersection this side, move on to next side */

    } /* intrsct_loop */

  Free(srch_corner_y);
  Free(srch_corner_x);

  /*
     If segment manages to cross over pole, shift the beginning
     endpoint in order to avoid hitting pole directly
     (it is ok for endpoint to be pole point)
  */

  if (fabs(*intrsct_x) < 1.e-10 && fabs(*intrsct_y) < 1.e-10 && (IS_NOT_EQUAL(endx, 0) && IS_NOT_EQUAL(endy, 0)))
    {
      if (*avoid_pole_count > 2)
        {
          *avoid_pole_count = 0;
          *avoid_pole_offset = 10. * (*avoid_pole_offset);
        }

      cross_product = begsegx * (endy - begsegy) - begsegy * (endx - begsegx);
      *intrsct_lat = begseg[0];
      if (cross_product * (*intrsct_lat) > ZERO)
        {
          *intrsct_lon = beglon + *avoid_pole_offset;
          begseg[1] = begseg[1] + *avoid_pole_offset;
        }
      else
        {
          *intrsct_lon = beglon - *avoid_pole_offset;
          begseg[1] = begseg[1] - *avoid_pole_offset;
        }

      *avoid_pole_count = *avoid_pole_count + 1;
      *luse_last = false;
    }
  else
    {
      *avoid_pole_count = 0;
      *avoid_pole_offset = TINY;
    }

  /*
     If the segment crosses a pole threshold, reset the intersection
     to be the threshold latitude and do not reuse x,y intersect
     on next entry.  Only check if did not cross threshold last
     time - sometimes the coordinate transformation can place a
     segment on the other side of the threshold again
  */
  if (*lthresh)
    {
      if (*intrsct_lat > north_thresh || *intrsct_lat < south_thresh) *lthresh = false;
    }
  else if (beglat > ZERO && *intrsct_lat < north_thresh)
    {
      mat4 = endlat - begseg[0];
      mat3 = endlon - begseg[1];
      if (mat3 > PI) mat3 = mat3 - PI2;
      if (mat3 < -PI) mat3 = mat3 + PI2;
      *intrsct_lat = north_thresh - TINY;
      s1 = (north_thresh - begseg[0]) / mat4;
      *intrsct_lon = begseg[1] + s1 * mat3;
      *luse_last = false;
      *lthresh = true;
    }
  else if (beglat<ZERO && * intrsct_lat> south_thresh)
    {
      mat4 = endlat - begseg[0];
      mat3 = endlon - begseg[1];
      if (mat3 > PI) mat3 = mat3 - PI2;
      if (mat3 < -PI) mat3 = mat3 + PI2;
      *intrsct_lat = south_thresh + TINY;
      s1 = (south_thresh - begseg[0]) / mat4;
      *intrsct_lon = begseg[1] + s1 * mat3;
      *luse_last = false;
      *lthresh = true;
    }

  /* If reached end of segment, do not use x,y intersect on next entry */

  if (IS_EQUAL(*intrsct_lat, endlat) && IS_EQUAL(*intrsct_lon, endlon)) *luse_last = false;

} /* pole_intersection */

/*
   This routine finds the next intersection of a destination grid line with
   the line segment given by beglon, endlon, etc.
   A coincidence flag is returned if the segment is entirely coincident with
   an ocean grid line.  The cells in which to search for an intersection must
   have already been restricted in the calling routine.
*/
static void
intersection(long *location, double *intrsct_lat, double *intrsct_lon, bool *lcoinc, double beglat, double beglon, double endlat,
             double endlon, double *begseg, bool lbegin, bool lrevers, long num_srch_cells, long srch_corners,
             const size_t *restrict srch_add, const double *restrict srch_corner_lat, const double *restrict srch_corner_lon,
             long *last_loc, bool *lthresh, double *intrsct_lat_off, double *intrsct_lon_off, bool *luse_last, double *intrsct_x,
             double *intrsct_y, int *avoid_pole_count, double *avoid_pole_offset)
{
  /*
    Intent(in):
    bool lbegin,            ! flag for first integration along this segment
    bool lrevers            ! flag whether segment integrated in reverse
    double beglat, beglon,  ! beginning lat/lon endpoints for segment
    double endlat, endlon   ! ending    lat/lon endpoints for segment

    Intent(inout) ::
    double *begseg          ! begin lat/lon of full segment

    intent(out):
    int *location           ! address in destination array containing this
    segment bool *lcoinc            ! flag segments which are entirely
    coincident with a grid line double *intrsct_lat, *intrsct_lon ! lat/lon
    coords of next intersect.
  */
  /* Local variables */
  long n, next_n, cell;
  long ioffset;

  int loutside; /* flags points outside grid */

  double lon1, lon2;         /* local longitude variables for segment */
  double lat1, lat2;         /* local latitude  variables for segment */
  double grdlon1, grdlon2;   /* local longitude variables for grid cell */
  double grdlat1, grdlat2;   /* local latitude  variables for grid cell */
  double vec1_lat, vec1_lon; /* vectors and cross products used */
  double vec2_lat, vec2_lon; /* during grid search */
  double cross_product;
  double eps, offset;                                /* small offset away from intersect */
  double s1, s2, determ;                             /* variables used for linear solve to */
  double mat1 = 0, mat2, mat3 = 0, mat4, rhs1, rhs2; /* find intersection */

  /* Initialize defaults, flags, etc. */

  *location = -1;
  *lcoinc = false;
  *intrsct_lat = endlat;
  *intrsct_lon = endlon;

  if (num_srch_cells == 0) return;

  if (beglat > north_thresh || beglat < south_thresh)
    {
      if (*lthresh) *location = *last_loc;
      pole_intersection(location, intrsct_lat, intrsct_lon, lcoinc, lthresh, beglat, beglon, endlat, endlon, begseg, lrevers,
                        num_srch_cells, srch_corners, srch_add, srch_corner_lat, srch_corner_lon, luse_last, intrsct_x, intrsct_y,
                        avoid_pole_count, avoid_pole_offset);

      if (*lthresh)
        {
          *last_loc = *location;
          *intrsct_lat_off = *intrsct_lat;
          *intrsct_lon_off = *intrsct_lon;
        }
      return;
    }

  loutside = FALSE;
  if (lbegin)
    {
      lat1 = beglat;
      lon1 = beglon;
    }
  else
    {
      lat1 = *intrsct_lat_off;
      lon1 = *intrsct_lon_off;
    }

  lat2 = endlat;
  lon2 = endlon;
  if ((lon2 - lon1) > THREE * PIH)
    lon2 -= PI2;
  else if ((lon2 - lon1) < -THREE * PIH)
    lon2 += PI2;

  s1 = ZERO;

  /*
     Search for location of this segment in ocean grid using cross
     product method to determine whether a point is enclosed by a cell
  */
  while (TRUE) /* srch_loop */
    {
      /* If last segment crossed threshold, use that location */

      if (*lthresh)
        {
          for (cell = 0; cell < num_srch_cells; ++cell)
            if (srch_add[cell] == (size_t) *last_loc)
              {
                *location = *last_loc;
                eps = TINY;
                goto after_srch_loop;
              }
        }

      /* Otherwise normal search algorithm */

      for (cell = 0; cell < num_srch_cells; ++cell) /* cell_loop  */
        {
          ioffset = cell * srch_corners;
          for (n = 0; n < srch_corners; n++) /* corner_loop */
            {
              next_n = (n + 1) % srch_corners;
              /*
                Here we take the cross product of the vector making
                up each cell side with the vector formed by the vertex
                and search point.  If all the cross products are
                positive, the point is contained in the cell.
              */
              vec1_lat = srch_corner_lat[ioffset + next_n] - srch_corner_lat[ioffset + n];
              vec1_lon = srch_corner_lon[ioffset + next_n] - srch_corner_lon[ioffset + n];
              vec2_lat = lat1 - srch_corner_lat[ioffset + n];
              vec2_lon = lon1 - srch_corner_lon[ioffset + n];

              /* If endpoint coincident with vertex, offset the endpoint */

              if (IS_EQUAL(vec2_lat, 0) && IS_EQUAL(vec2_lon, 0))
                {
                  lat1 += 1.e-10 * (lat2 - lat1);
                  lon1 += 1.e-10 * (lon2 - lon1);
                  vec2_lat = lat1 - srch_corner_lat[ioffset + n];
                  vec2_lon = lon1 - srch_corner_lon[ioffset + n];
                }

              /* Check for 0,2pi crossings */

              if (vec1_lon > PI)
                vec1_lon -= PI2;
              else if (vec1_lon < -PI)
                vec1_lon += PI2;

              if (vec2_lon > PI)
                vec2_lon -= PI2;
              else if (vec2_lon < -PI)
                vec2_lon += PI2;

              cross_product = vec1_lon * vec2_lat - vec2_lon * vec1_lat;

              /*
               If the cross product for a side is ZERO, the point
                 lies exactly on the side or the side is degenerate
                 (ZERO length).  If degenerate, set the cross
                 product to a positive number.  Otherwise perform
                 another cross product between the side and the
                 segment itself.
               If this cross product is also ZERO, the line is
                 coincident with the cell boundary - perform the
                 dot product and only choose the cell if the dot
                 product is positive (parallel vs anti-parallel).
              */
              if (IS_EQUAL(cross_product, 0))
                {
                  if (IS_NOT_EQUAL(vec1_lat, 0) || IS_NOT_EQUAL(vec1_lon, 0))
                    {
                      vec2_lat = lat2 - lat1;
                      vec2_lon = lon2 - lon1;

                      if (vec2_lon > PI)
                        vec2_lon -= PI2;
                      else if (vec2_lon < -PI)
                        vec2_lon += PI2;

                      cross_product = vec1_lon * vec2_lat - vec2_lon * vec1_lat;
                    }
                  else
                    cross_product = ONE;

                  if (IS_EQUAL(cross_product, 0))
                    {
                      *lcoinc = true;
                      cross_product = vec1_lon * vec2_lon + vec1_lat * vec2_lat;
                      if (lrevers) cross_product = -cross_product;
                    }
                }

              /* If cross product is less than ZERO, this cell doesn't work */

              if (cross_product < ZERO) break; /* corner_loop */

            } /* corner_loop */

          /* If cross products all positive, we found the location */

          if (n >= srch_corners)
            {
              *location = srch_add[cell];
              /*
                If the beginning of this segment was outside the
                grid, invert the segment so the intersection found
                will be the first intersection with the grid
              */
              if (loutside)
                {
                  lat2 = beglat;
                  lon2 = beglon;
                  *location = -1;
                  eps = -TINY;
                }
              else
                eps = TINY;

              goto after_srch_loop;
            }

          /* Otherwise move on to next cell */

        } /* cell_loop */

      /*
        If still no cell found, the point lies outside the grid.
        Take some baby steps along the segment to see if any
        part of the segment lies inside the grid.
      */
      loutside = TRUE;
      s1 = s1 + BABY_STEP;
      lat1 = beglat + s1 * (endlat - beglat);
      lon1 = beglon + s1 * (lon2 - beglon);

      /* Reached the end of the segment and still outside the grid return no
       * intersection */

      if (s1 >= ONE) return;

    } /* srch_loop */

after_srch_loop:

  /*
    Now that a cell is found, search for the next intersection.
    Loop over sides of the cell to find intersection with side
    must check all sides for coincidences or intersections
  */

  ioffset = cell * srch_corners;

  for (n = 0; n < srch_corners; ++n) /* intrsct_loop */
    {
      next_n = (n + 1) % srch_corners;

      grdlon1 = srch_corner_lon[ioffset + n];
      grdlon2 = srch_corner_lon[ioffset + next_n];
      grdlat1 = srch_corner_lat[ioffset + n];
      grdlat2 = srch_corner_lat[ioffset + next_n];

      /* Set up linear system to solve for intersection */

      mat1 = lat2 - lat1;
      mat2 = grdlat1 - grdlat2;
      mat3 = lon2 - lon1;
      mat4 = grdlon1 - grdlon2;
      rhs1 = grdlat1 - lat1;
      rhs2 = grdlon1 - lon1;

      if (mat3 > PI)
        mat3 -= PI2;
      else if (mat3 < -PI)
        mat3 += PI2;

      if (mat4 > PI)
        mat4 -= PI2;
      else if (mat4 < -PI)
        mat4 += PI2;

      if (rhs2 > PI)
        rhs2 -= PI2;
      else if (rhs2 < -PI)
        rhs2 += PI2;

      determ = mat1 * mat4 - mat2 * mat3;

      /*
         If the determinant is ZERO, the segments are either
           parallel or coincident.  Coincidences were detected
           above so do nothing.
         If the determinant is non-ZERO, solve for the linear
           parameters s for the intersection point on each line
           segment.
         If 0<s1,s2<1 then the segment intersects with this side.
           Return the point of intersection (adding a small
           number so the intersection is off the grid line).
      */
      if (fabs(determ) > 1.e-30)
        {
          s1 = (rhs1 * mat4 - mat2 * rhs2) / determ;
          s2 = (mat1 * rhs2 - rhs1 * mat3) / determ;

          if (s2 >= ZERO && s2 <= ONE && s1 >= ZERO && s1 <= ONE)
            {
              /*
                Recompute intersection based on full segment
                so intersections are consistent for both sweeps
              */
              if (!loutside)
                {
                  mat1 = lat2 - begseg[0];
                  mat3 = lon2 - begseg[1];
                  rhs1 = grdlat1 - begseg[0];
                  rhs2 = grdlon1 - begseg[1];
                }
              else
                {
                  mat1 = begseg[0] - endlat;
                  mat3 = begseg[1] - endlon;
                  rhs1 = grdlat1 - endlat;
                  rhs2 = grdlon1 - endlon;
                }

              if (mat3 > PI)
                mat3 -= PI2;
              else if (mat3 < -PI)
                mat3 += PI2;

              if (rhs2 > PI)
                rhs2 -= PI2;
              else if (rhs2 < -PI)
                rhs2 += PI2;

              determ = mat1 * mat4 - mat2 * mat3;

              /*
                Sometimes due to roundoff, the previous
                determinant is non-ZERO, but the lines
                are actually coincident.  If this is the
                case, skip the rest.
              */
              if (IS_NOT_EQUAL(determ, 0))
                {
                  s1 = (rhs1 * mat4 - mat2 * rhs2) / determ;
                  s2 = (mat1 * rhs2 - rhs1 * mat3) / determ;

                  offset = s1 + eps / determ;
                  if (offset > ONE) offset = ONE;

                  if (!loutside)
                    {
                      *intrsct_lat = begseg[0] + mat1 * s1;
                      *intrsct_lon = begseg[1] + mat3 * s1;
                      *intrsct_lat_off = begseg[0] + mat1 * offset;
                      *intrsct_lon_off = begseg[1] + mat3 * offset;
                    }
                  else
                    {
                      *intrsct_lat = endlat + mat1 * s1;
                      *intrsct_lon = endlon + mat3 * s1;
                      *intrsct_lat_off = endlat + mat1 * offset;
                      *intrsct_lon_off = endlon + mat3 * offset;
                    }
                  break; /* intrsct_loop */
                }
            }
        }

      /* No intersection this side, move on to next side */

    } /* intrsct_loop */

  /*
     If the segment crosses a pole threshold, reset the intersection
     to be the threshold latitude.  Only check if this was not a
     threshold segment since sometimes coordinate transform can end
     up on other side of threshold again.
  */
  if (*lthresh)
    {
      if (*intrsct_lat < north_thresh || *intrsct_lat > south_thresh) *lthresh = false;
    }
  else if (lat1 > ZERO && *intrsct_lat > north_thresh)
    {
      *intrsct_lat = north_thresh + TINY;
      *intrsct_lat_off = north_thresh + eps * mat1;
      s1 = (*intrsct_lat - begseg[0]) / mat1;
      *intrsct_lon = begseg[1] + s1 * mat3;
      *intrsct_lon_off = begseg[1] + (s1 + eps) * mat3;
      *last_loc = *location;
      *lthresh = true;
    }
  else if (lat1 < ZERO && *intrsct_lat < south_thresh)
    {
      *intrsct_lat = south_thresh - TINY;
      *intrsct_lat_off = south_thresh + eps * mat1;
      s1 = (*intrsct_lat - begseg[0]) / mat1;
      *intrsct_lon = begseg[1] + s1 * mat3;
      *intrsct_lon_off = begseg[1] + (s1 + eps) * mat3;
      *last_loc = *location;
      *lthresh = true;
    }

} /* intersection */

/*
   This routine computes the line integral of the flux function
   that results in the interpolation weights.  The line is defined
   by the input lat/lon of the endpoints.
*/
static double
phi_gradient(double in_phi1, double in_phi2, double dphi, double f1, double f2, double grid_lon)
{
  double fint, fac;
  double weight;
  double phi1, phi2;

  phi1 = in_phi1 - grid_lon;
  if (phi1 > PI)
    phi1 -= PI2;
  else if (phi1 < -PI)
    phi1 += PI2;

  phi2 = in_phi2 - grid_lon;
  if (phi2 > PI)
    phi2 -= PI2;
  else if (phi2 < -PI)
    phi2 += PI2;

  if ((phi2 - phi1) < PI && (phi2 - phi1) > -PI)
    weight = dphi * (phi1 * f1 + phi2 * f2);
  else
    {
      if (phi1 > ZERO)
        fac = PI;
      else
        fac = -PI;

      fint = f1 + (f2 - f1) * (fac - phi1) / fabs(dphi);
      weight = HALF * phi1 * (phi1 - fac) * f1 - HALF * phi2 * (phi2 + fac) * f2 + HALF * fac * (phi1 + phi2) * fint;
    }

  return weight;
}

static void
line_integral(double *weights, double in_phi1, double in_phi2, double theta1, double theta2, double grid1_lon, double grid2_lon)
{
  /*
    Intent(in):
    double in_phi1, in_phi2,     ! Longitude endpoints for the segment
    double theta1, theta2,       ! Latitude  endpoints for the segment
    double grid1_lon,            ! Reference coordinates for each
    double grid2_lon             ! Grid (to ensure correct 0,2pi interv.)

    Intent(out):
    double weights[6]            ! Line integral contribution to weights
  */

  /*  Local variables  */
  double dphi, sinth1, sinth2, costh1, costh2;
  double f1, f2;

  /*  Weights for the general case based on a trapezoidal approx to the integrals. */

  sinth1 = sin(theta1);
  sinth2 = sin(theta2);
  costh1 = cos(theta1);
  costh2 = cos(theta2);

  dphi = in_phi1 - in_phi2;
  if (dphi > PI)
    dphi -= PI2;
  else if (dphi < -PI)
    dphi += PI2;

  dphi = HALF * dphi;

  /*
     The first weight is the area overlap integral. The second and
     fourth are second-order latitude gradient weights.
  */
  weights[0] = dphi * (sinth1 + sinth2);
  weights[1] = dphi * (costh1 + costh2 + (theta1 * sinth1 + theta2 * sinth2));
  weights[3] = weights[0];
  weights[4] = weights[1];

  /*
     The third and fifth weights are for the second-order phi gradient
     component.  Must be careful of longitude range.
  */
  f1 = HALF * (costh1 * sinth1 + theta1);
  f2 = HALF * (costh2 * sinth2 + theta2);

  weights[2] = phi_gradient(in_phi1, in_phi2, dphi, f1, f2, grid1_lon);
  weights[5] = phi_gradient(in_phi1, in_phi2, dphi, f1, f2, grid2_lon);

} /* line_integral */

static void
correct_pole(RemapGrid *src_grid, RemapGrid *tgt_grid, RemapVars &rv, double *src_centroid_lat, double *src_centroid_lon,
             double *tgt_centroid_lat, double *tgt_centroid_lon, grid_store_t *grid_store)
{
  /*
     Correct for situations where N/S pole not explicitly included in
     grid (i.e. as a grid corner point). If pole is missing from only
     one grid, need to correct only the area and centroid of that
     grid.  If missing from both, do complete weight calculation.
  */
  long n;
  long num_wts;
  long src_grid_size;
  long tgt_grid_size;
  long src_cell_add; /* current linear address for source grid cell   */
  long tgt_cell_add; /* current linear address for target grid cell   */
  double weights[6]; /* local wgt array */

  num_wts = rv.num_wts;

  src_grid_size = src_grid->size;
  tgt_grid_size = tgt_grid->size;

  /* North Pole */
  weights[0] = PI2;
  weights[1] = PI * PI;
  weights[2] = ZERO;
  weights[3] = PI2;
  weights[4] = PI * PI;
  weights[5] = ZERO;

  src_cell_add = src_grid_size;
  /* pole_loop1 */
  for (n = 0; n < src_grid_size; ++n)
    if (src_grid->cell_area[n] < -THREE * PIH && src_grid->cell_center_lat[n] > ZERO)
      {
        src_cell_add = n;
        break;
      }

  tgt_cell_add = tgt_grid_size;
  /* pole_loop2 */
  for (n = 0; n < tgt_grid_size; ++n)
    if (tgt_grid->cell_area[n] < -THREE * PIH && tgt_grid->cell_center_lat[n] > ZERO)
      {
        tgt_cell_add = n;
        break;
      }

  if (src_cell_add != src_grid_size)
    {
      src_grid->cell_area[src_cell_add] += weights[0];
      src_centroid_lat[src_cell_add] += weights[1];
      src_centroid_lon[src_cell_add] += weights[2];
    }

  if (tgt_cell_add != tgt_grid_size)
    {
      tgt_grid->cell_area[tgt_cell_add] += weights[3];
      tgt_centroid_lat[tgt_cell_add] += weights[4];
      tgt_centroid_lon[tgt_cell_add] += weights[5];
    }

  if (src_cell_add != src_grid_size && tgt_cell_add != tgt_grid_size)
    {
      store_link_cnsrv(rv, src_cell_add, tgt_cell_add, num_wts, weights, grid_store);

      src_grid->cell_frac[src_cell_add] += weights[0];
      tgt_grid->cell_frac[tgt_cell_add] += weights[3];
    }

  /* South Pole */
  weights[0] = PI2;
  weights[1] = -PI * PI;
  weights[2] = ZERO;
  weights[3] = PI2;
  weights[4] = -PI * PI;
  weights[5] = ZERO;

  src_cell_add = src_grid_size;
  /* pole_loop3 */
  for (n = 0; n < src_grid_size; ++n)
    if (src_grid->cell_area[n] < -THREE * PIH && src_grid->cell_center_lat[n] < ZERO)
      {
        src_cell_add = n;
        break;
      }

  tgt_cell_add = tgt_grid_size;
  /* pole_loop4 */
  for (n = 0; n < tgt_grid_size; ++n)
    if (tgt_grid->cell_area[n] < -THREE * PIH && tgt_grid->cell_center_lat[n] < ZERO)
      {
        tgt_cell_add = n;
        break;
      }

  if (src_cell_add != src_grid_size)
    {
      src_grid->cell_area[src_cell_add] += weights[0];
      src_centroid_lat[src_cell_add] += weights[1];
      src_centroid_lon[src_cell_add] += weights[2];
    }

  if (tgt_cell_add != tgt_grid_size)
    {
      tgt_grid->cell_area[tgt_cell_add] += weights[3];
      tgt_centroid_lat[tgt_cell_add] += weights[4];
      tgt_centroid_lon[tgt_cell_add] += weights[5];
    }

  if (src_cell_add != src_grid_size && tgt_cell_add != tgt_grid_size)
    {
      store_link_cnsrv(rv, src_cell_add, tgt_cell_add, num_wts, weights, grid_store);

      src_grid->cell_frac[src_cell_add] += weights[0];
      tgt_grid->cell_frac[tgt_cell_add] += weights[3];
    }
}

static inline void
norm_weight(double norm_factor, double *weights, double src_centroid_lat, double src_centroid_lon)
{
  double weight0 = weights[0];

  weights[0] = weight0 * norm_factor;
  weights[1] = (weights[1] - weight0 * src_centroid_lat) * norm_factor;
  weights[2] = (weights[2] - weight0 * src_centroid_lon) * norm_factor;
}

static void
normalize_weights(RemapGrid *tgt_grid, RemapVars &rv, double *src_centroid_lat, double *src_centroid_lon)
{
  /* Include centroids in weights and normalize using destination area if requested */
  long num_links = rv.num_links;
  long src_cell_add; /* current linear address for source grid cell   */
  long tgt_cell_add; /* current linear address for target grid cell   */
  double *weights = &rv.wts[0];
  double norm_factor = 0; /* factor for normalizing wts */

  if (rv.normOpt == NormOpt::DESTAREA)
    {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(num_links, rv, weights, tgt_grid, src_centroid_lat, \
                                              src_centroid_lon) private(src_cell_add, tgt_cell_add, norm_factor)
#endif
      for (long n = 0; n < num_links; ++n)
        {
          src_cell_add = rv.src_cell_add[n];
          tgt_cell_add = rv.tgt_cell_add[n];
          norm_factor = IS_NOT_EQUAL(tgt_grid->cell_area[tgt_cell_add], 0) ? ONE / tgt_grid->cell_area[tgt_cell_add] : ZERO;
          norm_weight(norm_factor, &weights[n * 3], src_centroid_lat[src_cell_add], src_centroid_lon[src_cell_add]);
        }
    }
  else if (rv.normOpt == NormOpt::FRACAREA)
    {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(num_links, rv, weights, tgt_grid, src_centroid_lat, \
                                              src_centroid_lon) private(src_cell_add, tgt_cell_add, norm_factor)
#endif
      for (long n = 0; n < num_links; ++n)
        {
          src_cell_add = rv.src_cell_add[n];
          tgt_cell_add = rv.tgt_cell_add[n];
          norm_factor = IS_NOT_EQUAL(tgt_grid->cell_frac[tgt_cell_add], 0) ? ONE / tgt_grid->cell_frac[tgt_cell_add] : ZERO;
          norm_weight(norm_factor, &weights[n * 3], src_centroid_lat[src_cell_add], src_centroid_lon[src_cell_add]);
        }
    }
  else if (rv.normOpt == NormOpt::NONE)
    {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(num_links, rv, weights, tgt_grid, src_centroid_lat, \
                                              src_centroid_lon) private(src_cell_add, norm_factor)
#endif
      for (long n = 0; n < num_links; ++n)
        {
          src_cell_add = rv.src_cell_add[n];
          norm_factor = ONE;
          norm_weight(norm_factor, &weights[n * 3], src_centroid_lat[src_cell_add], src_centroid_lon[src_cell_add]);
        }
    }
}

/*
  -----------------------------------------------------------------------

   This routine traces the perimeters of every grid cell on each
   grid checking for intersections with the other grid and computing
   line integrals for each subsegment.

  -----------------------------------------------------------------------
*/
void
remapConservWeightsScrip(RemapSearch &rsearch, RemapVars &rv)
{
  RemapGrid *src_grid = rsearch.srcGrid;
  RemapGrid *tgt_grid = rsearch.tgtGrid;

  bool lcheck = true;

  long max_subseg = 100000; /* max number of subsegments per segment to prevent infinite loop */
  /* 1000 is too small!!! */

  long srch_corners; /* num of corners of srch cells           */
  long i;

  double findex = 0;

  if (cdoVerbose) cdoPrint("Called %s()", __func__);

  progressInit();

  long num_wts = rv.num_wts;

  grid_store_t *grid_store = (grid_store_t *) Malloc(sizeof(grid_store_t));
  grid_store_init(grid_store, tgt_grid->size);

  if (cdoVerbose)
    {
      cdoPrint("North threshhold: %g", north_thresh);
      cdoPrint("South threshhold: %g", south_thresh);
    }

  double start = cdoVerbose ? cdo_get_wtime() : 0;

  long src_grid_size = src_grid->size;
  long tgt_grid_size = tgt_grid->size;

  long src_num_cell_corners = src_grid->num_cell_corners;
  long tgt_num_cell_corners = tgt_grid->num_cell_corners;

  /* Initialize centroid arrays */

  double *src_centroid_lat = (double *) Malloc(src_grid_size * sizeof(double));
  double *src_centroid_lon = (double *) Malloc(src_grid_size * sizeof(double));
  double *tgt_centroid_lat = (double *) Malloc(tgt_grid_size * sizeof(double));
  double *tgt_centroid_lon = (double *) Malloc(tgt_grid_size * sizeof(double));

  for (long n = 0; n < src_grid_size; ++n)
    {
      src_centroid_lat[n] = 0;
      src_centroid_lon[n] = 0;
    }

  for (long n = 0; n < tgt_grid_size; ++n)
    {
      tgt_centroid_lat[n] = 0;
      tgt_centroid_lon[n] = 0;
    }

  std::vector<double *> srch_corner_lat(Threading::ompNumThreads);  // lat of each corner of srch cells
  std::vector<double *> srch_corner_lon(Threading::ompNumThreads);  // lon of each corner of srch cells
  std::vector<long> max_srch_cells(Threading::ompNumThreads);       // num cells in restricted search arrays

  /*  Integrate around each cell on source grid */

  for (i = 0; i < Threading::ompNumThreads; ++i)
    {
      srch_corner_lat[i] = NULL;
      srch_corner_lon[i] = NULL;
      max_srch_cells[i] = 0;
    }

  std::vector<size_t *> srch_add(Threading::ompNumThreads);  // global address of cells in srch arrays
  for (i = 0; i < Threading::ompNumThreads; ++i) srch_add[i] = new size_t[tgt_grid_size];

  srch_corners = tgt_num_cell_corners;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(findex, rsearch, num_wts, src_centroid_lon, src_centroid_lat, grid_store, rv, \
                                              cdoVerbose, max_subseg, srch_corner_lat, srch_corner_lon, max_srch_cells,     \
                                              src_num_cell_corners, srch_corners, src_grid, tgt_grid, tgt_grid_size,        \
                                              src_grid_size, srch_add)
#endif
  for (long src_cell_add = 0; src_cell_add < src_grid_size; ++src_cell_add)
    {
      int ompthID = cdo_omp_get_thread_num();

#ifdef _OPENMP
#pragma omp atomic
#endif
      findex++;
      if (ompthID == 0) progressStatus(0, 0.5, findex / src_grid_size);

      bool lcoinc;  // flag for coincident segments
      bool lthresh = false;
      bool luse_last = false;
      int avoid_pole_count = 0;                         // count attempts to avoid pole
      double avoid_pole_offset = TINY;                  // endpoint offset to avoid pole
      double intrsct_lat, intrsct_lon;                  // lat/lon of next intersect
      double intrsct_lat_off = 0, intrsct_lon_off = 0;  // lat/lon coords offset for next search
      double intrsct_x, intrsct_y;                      // x,y for intersection
      long last_loc = -1;                               // save location when crossing threshold
      long tgt_cell_add;

      // Get search cells
      long num_srch_cells = get_srch_cells(src_cell_add, rsearch.srcBins, rsearch.tgtBins,
                                           &rsearch.srcBins.cell_bound_box[src_cell_add * 4], srch_add[ompthID]);

      if (num_srch_cells == 0) continue;

      // Create search arrays
      if (num_srch_cells > max_srch_cells[ompthID])
        {
          srch_corner_lat[ompthID] = (double *) realloc(srch_corner_lat[ompthID], srch_corners * num_srch_cells * sizeof(double));
          srch_corner_lon[ompthID] = (double *) realloc(srch_corner_lon[ompthID], srch_corners * num_srch_cells * sizeof(double));
          max_srch_cells[ompthID] = num_srch_cells;
        }

      // gather1
      long ioffset;
      for (long n = 0; n < num_srch_cells; ++n)
        {
          tgt_cell_add = srch_add[ompthID][n];
          ioffset = tgt_cell_add * srch_corners;

          long nsrch_corners = n * srch_corners;
          for (long k = 0; k < srch_corners; k++)
            {
              srch_corner_lat[ompthID][nsrch_corners + k] = tgt_grid->cell_corner_lat[ioffset + k];
              srch_corner_lon[ompthID][nsrch_corners + k] = tgt_grid->cell_corner_lon[ioffset + k];
            }
        }

      // Integrate around this cell
      ioffset = src_cell_add * src_num_cell_corners;
      for (long corner = 0; corner < src_num_cell_corners; ++corner)
        {
          long next_corn = (corner + 1) % src_num_cell_corners;

          /* Define endpoints of the current segment */

          double beglat = src_grid->cell_corner_lat[ioffset + corner];
          double beglon = src_grid->cell_corner_lon[ioffset + corner];
          double endlat = src_grid->cell_corner_lat[ioffset + next_corn];
          double endlon = src_grid->cell_corner_lon[ioffset + next_corn];
          bool lrevers = false;

          /*  To ensure exact path taken during both sweeps, always integrate
           * segments in the same direction (SW to NE). */
          if ((endlat < beglat) || (IS_EQUAL(endlat, beglat) && endlon < beglon))
            {
              beglat = src_grid->cell_corner_lat[ioffset + next_corn];
              beglon = src_grid->cell_corner_lon[ioffset + next_corn];
              endlat = src_grid->cell_corner_lat[ioffset + corner];
              endlon = src_grid->cell_corner_lon[ioffset + corner];
              lrevers = true;
            }

          /*
            If this is a constant-longitude segment, skip the rest
            since the line integral contribution will be ZERO.
          */
          if (IS_EQUAL(endlon, beglon)) continue;

          double weights[6];  // local wgt array
          double begseg[2];   // begin lat/lon for full segment
          begseg[0] = beglat;
          begseg[1] = beglon;
          bool lbegin = true;

          long num_subseg = 0;
          /*
            Integrate along this segment, detecting intersections
            and computing the line integral for each sub-segment
          */
          while (IS_NOT_EQUAL(beglat, endlat) || IS_NOT_EQUAL(beglon, endlon))
            {
              /*  Prevent infinite loops if integration gets stuck near cell or threshold boundary */
              num_subseg++;
              if (num_subseg >= max_subseg)
                cdoAbort("Integration stalled: num_subseg exceeded limit (grid1[%zu]: lon1=%g lon2=%g lat1=%g lat2=%g)!",
                         src_cell_add, beglon, endlon, beglat, endlat);

              /* Uwe Schulzweida: skip very small regions */
              if (num_subseg % 1000 == 0)
                {
                  if (fabs(beglat - endlat) < 1.e-10 || fabs(beglon - endlon) < 1.e-10)
                    {
                      if (cdoVerbose)
                        cdoPrint("Skip very small region (grid1[%ld]): lon=%g dlon=%g lat=%g dlat=%g", src_cell_add, beglon,
                                 endlon - beglon, beglat, endlat - beglat);
                      break;
                    }
                }

              /* Find next intersection of this segment with a gridline on grid 2. */

              intersection(&tgt_cell_add, &intrsct_lat, &intrsct_lon, &lcoinc, beglat, beglon, endlat, endlon, begseg, lbegin,
                           lrevers, num_srch_cells, srch_corners, srch_add[ompthID], srch_corner_lat[ompthID],
                           srch_corner_lon[ompthID], &last_loc, &lthresh, &intrsct_lat_off, &intrsct_lon_off, &luse_last,
                           &intrsct_x, &intrsct_y, &avoid_pole_count, &avoid_pole_offset);

              lbegin = false;

              /* Compute line integral for this subsegment. */

              if (tgt_cell_add != -1)
                line_integral(weights, beglon, intrsct_lon, beglat, intrsct_lat, src_grid->cell_center_lon[src_cell_add],
                              tgt_grid->cell_center_lon[tgt_cell_add]);
              else
                line_integral(weights, beglon, intrsct_lon, beglat, intrsct_lat, src_grid->cell_center_lon[src_cell_add],
                              src_grid->cell_center_lon[src_cell_add]);

              /* If integrating in reverse order, change sign of weights */

              if (lrevers)
                for (long k = 0; k < 6; ++k) weights[k] = -weights[k];

              /*
                Store the appropriate addresses and weights.
                Also add contributions to cell areas and centroids.
              */
              if (tgt_cell_add != -1)
                if (src_grid->mask[src_cell_add])
                  {
#ifdef _OPENMP
#pragma omp critical
#endif
                    {
                      store_link_cnsrv(rv, src_cell_add, tgt_cell_add, num_wts, weights, grid_store);

                      tgt_grid->cell_frac[tgt_cell_add] += weights[3];
                    }
                    src_grid->cell_frac[src_cell_add] += weights[0];
                  }

              src_grid->cell_area[src_cell_add] += weights[0];
              src_centroid_lat[src_cell_add] += weights[1];
              src_centroid_lon[src_cell_add] += weights[2];

              /* Reset beglat and beglon for next subsegment. */
              beglat = intrsct_lat;
              beglon = intrsct_lon;
            }
          /* End of segment */
        }
    }

  /* Finished with all cells: deallocate search arrays */

  for (i = 0; i < Threading::ompNumThreads; ++i)
    {
      free(srch_corner_lon[i]);
      free(srch_corner_lat[i]);
    }

  for (i = 0; i < Threading::ompNumThreads; ++i) delete[] srch_add[i];

  /* Integrate around each cell on target grid */

  for (i = 0; i < Threading::ompNumThreads; ++i)
    {
      srch_corner_lat[i] = NULL;
      srch_corner_lon[i] = NULL;
    }

  for (i = 0; i < Threading::ompNumThreads; ++i) max_srch_cells[i] = 0;
  for (i = 0; i < Threading::ompNumThreads; ++i) srch_add[i] = new size_t[src_grid_size];

  srch_corners = src_num_cell_corners;

  findex = 0;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(findex, rsearch, num_wts, tgt_centroid_lon, tgt_centroid_lat, grid_store, rv, \
                                              cdoVerbose, max_subseg, srch_corner_lat, srch_corner_lon, max_srch_cells,     \
                                              tgt_num_cell_corners, srch_corners, src_grid, tgt_grid, tgt_grid_size,        \
                                              src_grid_size, srch_add)
#endif
  for (long tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add)
    {
      int ompthID = cdo_omp_get_thread_num();

#ifdef _OPENMP
#pragma omp atomic
#endif
      findex++;
      if (ompthID == 0) progressStatus(0.5, 0.5, findex / tgt_grid_size);

      bool lcoinc;  // flag for coincident segments
      bool lthresh = false;
      bool luse_last = false;
      int avoid_pole_count = 0;                         // count attempts to avoid pole
      double avoid_pole_offset = TINY;                  // endpoint offset to avoid pole
      double intrsct_lat, intrsct_lon;                  // lat/lon of next intersect
      double intrsct_lat_off = 0, intrsct_lon_off = 0;  // lat/lon coords offset for next search
      double intrsct_x, intrsct_y;                      // x,y for intersection
      long last_loc = -1;                               // save location when crossing threshold
      long src_cell_add;

      // Get search cells
      long num_srch_cells = get_srch_cells(tgt_cell_add, rsearch.tgtBins, rsearch.srcBins,
                                           &rsearch.tgtBins.cell_bound_box[tgt_cell_add * 4], srch_add[ompthID]);

      if (num_srch_cells == 0) continue;

      // Create search arrays
      if (num_srch_cells > max_srch_cells[ompthID])
        {
          srch_corner_lat[ompthID] = (double *) realloc(srch_corner_lat[ompthID], srch_corners * num_srch_cells * sizeof(double));
          srch_corner_lon[ompthID] = (double *) realloc(srch_corner_lon[ompthID], srch_corners * num_srch_cells * sizeof(double));
          max_srch_cells[ompthID] = num_srch_cells;
        }

      // gather2
      long ioffset;
      for (long n = 0; n < num_srch_cells; ++n)
        {
          src_cell_add = srch_add[ompthID][n];
          ioffset = src_cell_add * srch_corners;

          long nsrch_corners = n * srch_corners;
          for (long k = 0; k < srch_corners; ++k)
            {
              srch_corner_lat[ompthID][nsrch_corners + k] = src_grid->cell_corner_lat[ioffset + k];
              srch_corner_lon[ompthID][nsrch_corners + k] = src_grid->cell_corner_lon[ioffset + k];
            }
        }

      // Integrate around this cell
      ioffset = tgt_cell_add * tgt_num_cell_corners;
      for (long corner = 0; corner < tgt_num_cell_corners; ++corner)
        {
          long next_corn = (corner + 1) % tgt_num_cell_corners;

          /* Define endpoints of the current segment */
          double beglat = tgt_grid->cell_corner_lat[ioffset + corner];
          double beglon = tgt_grid->cell_corner_lon[ioffset + corner];
          double endlat = tgt_grid->cell_corner_lat[ioffset + next_corn];
          double endlon = tgt_grid->cell_corner_lon[ioffset + next_corn];
          bool lrevers = false;

          /* To ensure exact path taken during both sweeps, always integrate in the same direction */
          if ((endlat < beglat) || (IS_EQUAL(endlat, beglat) && endlon < beglon))
            {
              beglat = tgt_grid->cell_corner_lat[ioffset + next_corn];
              beglon = tgt_grid->cell_corner_lon[ioffset + next_corn];
              endlat = tgt_grid->cell_corner_lat[ioffset + corner];
              endlon = tgt_grid->cell_corner_lon[ioffset + corner];
              lrevers = true;
            }

          /*
            If this is a constant-longitude segment, skip the rest
            since the line integral contribution will be ZERO.
          */
          if (IS_EQUAL(endlon, beglon)) continue;

          double weights[6];  // local wgt array
          double begseg[2];   // begin lat/lon for full segment
          begseg[0] = beglat;
          begseg[1] = beglon;
          bool lbegin = true;

          long num_subseg = 0;
          /*
            Integrate along this segment, detecting intersections
            and computing the line integral for each sub-segment
          */
          while (IS_NOT_EQUAL(beglat, endlat) || IS_NOT_EQUAL(beglon, endlon))
            {
              /*  Prevent infinite loops if integration gets stuck near cell or threshold boundary */
              num_subseg++;
              if (num_subseg >= max_subseg)
                cdoAbort("Integration stalled: num_subseg exceeded limit (grid2[%zu]: lon1=%g lon2=%g lat1=%g lat2=%g)!",
                         tgt_cell_add, beglon, endlon, beglat, endlat);

              /* Uwe Schulzweida: skip very small regions */
              if (num_subseg % 1000 == 0)
                {
                  if (fabs(beglat - endlat) < 1.e-10 || fabs(beglon - endlon) < 1.e-10)
                    {
                      if (cdoVerbose)
                        cdoPrint("Skip very small region (grid2[%ld]): lon=%g dlon=%g lat=%g dlat=%g", tgt_cell_add, beglon,
                                 endlon - beglon, beglat, endlat - beglat);
                      break;
                    }
                }

              /* Find next intersection of this segment with a gridline on grid 2. */

              intersection(&src_cell_add, &intrsct_lat, &intrsct_lon, &lcoinc, beglat, beglon, endlat, endlon, begseg, lbegin,
                           lrevers, num_srch_cells, srch_corners, srch_add[ompthID], srch_corner_lat[ompthID],
                           srch_corner_lon[ompthID], &last_loc, &lthresh, &intrsct_lat_off, &intrsct_lon_off, &luse_last,
                           &intrsct_x, &intrsct_y, &avoid_pole_count, &avoid_pole_offset);

              lbegin = false;

              /* Compute line integral for this subsegment. */

              if (src_cell_add != -1)
                line_integral(weights, beglon, intrsct_lon, beglat, intrsct_lat, src_grid->cell_center_lon[src_cell_add],
                              tgt_grid->cell_center_lon[tgt_cell_add]);
              else
                line_integral(weights, beglon, intrsct_lon, beglat, intrsct_lat, tgt_grid->cell_center_lon[tgt_cell_add],
                              tgt_grid->cell_center_lon[tgt_cell_add]);

              /* If integrating in reverse order, change sign of weights */

              if (lrevers)
                for (long k = 0; k < 6; ++k) weights[k] = -weights[k];

              /*
                Store the appropriate addresses and weights.
                Also add contributions to cell areas and centroids.
                If there is a coincidence, do not store weights
                because they have been captured in the previous loop.
                The source grid mask is the master mask
              */
              if (!lcoinc && src_cell_add != -1)
                if (src_grid->mask[src_cell_add])
                  {
#ifdef _OPENMP
#pragma omp critical
#endif
                    {
                      store_link_cnsrv(rv, src_cell_add, tgt_cell_add, num_wts, weights, grid_store);

                      src_grid->cell_frac[src_cell_add] += weights[0];
                    }
                    tgt_grid->cell_frac[tgt_cell_add] += weights[3];
                  }

              tgt_grid->cell_area[tgt_cell_add] += weights[3];
              tgt_centroid_lat[tgt_cell_add] += weights[4];
              tgt_centroid_lon[tgt_cell_add] += weights[5];

              /* Reset beglat and beglon for next subsegment. */
              beglat = intrsct_lat;
              beglon = intrsct_lon;
            }
          /* End of segment */
        }
    }

  progressStatus(0, 1, 1);

  /* Finished with all cells: deallocate search arrays */

  for (i = 0; i < Threading::ompNumThreads; ++i)
    {
      free(srch_corner_lon[i]);
      free(srch_corner_lat[i]);
    }

  for (i = 0; i < Threading::ompNumThreads; ++i) delete[] srch_add[i];

  /*
     Correct for situations where N/S pole not explicitly included in
     grid (i.e. as a grid corner point). If pole is missing from only
     one grid, need to correct only the area and centroid of that
     grid.  If missing from both, do complete weight calculation.
  */
  correct_pole(src_grid, tgt_grid, rv, src_centroid_lat, src_centroid_lon, tgt_centroid_lat, tgt_centroid_lon, grid_store);

  grid_store_delete(grid_store);
  Free(grid_store);

  if (rv.num_links != rv.max_links) remapVarsResize(rv, rv.num_links);

  /* Finish centroid computation */

  for (long n = 0; n < src_grid_size; ++n)
    if (IS_NOT_EQUAL(src_grid->cell_area[n], 0))
      {
        src_centroid_lat[n] /= src_grid->cell_area[n];
        src_centroid_lon[n] /= src_grid->cell_area[n];
      }

  for (long n = 0; n < tgt_grid_size; ++n)
    if (IS_NOT_EQUAL(tgt_grid->cell_area[n], 0))
      {
        tgt_centroid_lat[n] /= tgt_grid->cell_area[n];
        tgt_centroid_lon[n] /= tgt_grid->cell_area[n];
      }

  /* 2010-10-08 Uwe Schulzweida: remove all links with weights < 0 */

  /*
  if ( 1 )
    {
      num_links = rv.num_links;

      if ( cdoVerbose )
        for ( long n = 0; n < num_links; n++ )
          printf("wts1: %d %g\n", n, rv.wts[3*n]);

      for ( long n = 0; n < num_links; n++ )
        {
          if ( rv.wts[3*n] < 0 )
            {
              int i, n2, nd;

              for ( n2 = n+1; n2 < num_links; n2++ )
                if ( rv.wts[3*n2] >= 0 ) break;

              nd = n2-n;
              num_links -= nd;
              for ( i = n; i < num_links; i++ )
                {
                  rv.wts[3*i]   = rv.wts[3*(i+nd)];
                  rv.wts[3*i+1] = rv.wts[3*(i+nd)+1];
                  rv.wts[3*i+2] = rv.wts[3*(i+nd)+2];

                  rv.src_cell_add[i] = rv.src_cell_add[i+nd];
                  rv.tgt_cell_add[i] = rv.tgt_cell_add[i+nd];
                }
            }
        }

     if ( cdoVerbose ) cdoPrint("Removed number of links = %zu", rv.num_links -
  num_links);

      rv.num_links = num_links;
    }
  */

  /* Include centroids in weights and normalize using destination area if requested */
  normalize_weights(tgt_grid, rv, src_centroid_lat, src_centroid_lon);

  long num_links = rv.num_links;

  if (cdoVerbose) cdoPrint("Total number of links = %zu", rv.num_links);

  for (long n = 0; n < src_grid_size; ++n)
    if (IS_NOT_EQUAL(src_grid->cell_area[n], 0)) src_grid->cell_frac[n] /= src_grid->cell_area[n];

  for (long n = 0; n < tgt_grid_size; ++n)
    if (IS_NOT_EQUAL(tgt_grid->cell_area[n], 0)) tgt_grid->cell_frac[n] /= tgt_grid->cell_area[n];

  /* Perform some error checking on final weights  */

  if (lcheck)
    {
      remapCheckArea(src_grid_size, &src_grid->cell_area[0], "Source");
      remapCheckArea(tgt_grid_size, &tgt_grid->cell_area[0], "Target");

      for (long n = 0; n < src_grid_size; ++n)
        {
          if (src_centroid_lat[n] < -PIH - .01 || src_centroid_lat[n] > PIH + .01)
            cdoPrint("Source grid centroid lat error: %ld %g", n, src_centroid_lat[n]);

          src_centroid_lat[n] = 0;
          src_centroid_lon[n] = 0;
        }

      for (long n = 0; n < tgt_grid_size; ++n)
        {
          if (tgt_centroid_lat[n] < -PIH - .01 || tgt_centroid_lat[n] > PIH + .01)
            cdoPrint("Target grid centroid lat error: %ld %g", n, tgt_centroid_lat[n]);

          tgt_centroid_lat[n] = 0;
          tgt_centroid_lon[n] = 0;
        }

      remapVarsCheckWeights(rv);

      for (long n = 0; n < num_links; ++n)
        {
          long tgt_cell_add = rv.tgt_cell_add[n];
          tgt_centroid_lat[tgt_cell_add] += rv.wts[3 * n];
        }

      /* 2012-01-24 Uwe Schulzweida: changed [tgt_cell_add] to [n] (bug fix) */
      double norm_factor = 0;  // factor for normalizing wts
      for (long n = 0; n < tgt_grid_size; ++n)
        {
          if (rv.normOpt == NormOpt::DESTAREA)
            norm_factor = tgt_grid->cell_frac[n];
          else if (rv.normOpt == NormOpt::FRACAREA)
            norm_factor = ONE;
          else if (rv.normOpt == NormOpt::NONE)
            norm_factor = tgt_grid->cell_area[n];

          if (tgt_centroid_lat[n] > 0 && fabs(tgt_centroid_lat[n] - norm_factor) > .01)
            cdoPrint("Error: sum of wts for map1 %ld %g %g", n, tgt_centroid_lat[n], norm_factor);
        }
    }  // lcheck

  Free(src_centroid_lat);
  Free(src_centroid_lon);
  Free(tgt_centroid_lat);
  Free(tgt_centroid_lon);

  if (cdoVerbose) cdoPrint("Cells search: %.2f seconds", cdo_get_wtime() - start);
}  // remapConservWeightsScrip
