/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2009 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
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

      Gridboxstat    gridboxmin          Gridbox minimum
      Gridboxstat    gridboxmax          Gridbox maximum
      Gridboxstat    gridboxsum          Gridbox sum
      Gridboxstat    gridboxmean         Gridbox mean
      Gridboxstat    gridboxavg          Gridbox average
      Gridboxstat    gridboxstd          Gridbox standard deviation
      Gridboxstat    gridboxvar          Gridbox variance
*/


#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "functs.h"
#define fast 1

static
int genBoxGrid(int gridID1, int xinc, int yinc)
{
  static char func[] = "genBoxGrid";
  int i, j, i1;
  int gridID2 = -1, gridtype;
  int gridsize1 = 0, nlon1 = 0, nlat1 = 0;
  int gridsize2 = 0, nlon2 = 0, nlat2 = 0;
  int x1, y1, x2, y2;
  int use_x1, use_y1;
  int corner, add, g2_add; 
  int g1_add;
  int circular = gridIsCircular(gridID1) ;
  double *xvals1, *yvals1, *xvals2, *yvals2;
  double *grid1_corner_lon = NULL, *grid1_corner_lat = NULL;
  double *grid2_corner_lon = NULL, *grid2_corner_lat = NULL;
  double on_up, on_lo, ox_up, ox_lo, an_le, an_ri, ax_le, ax_ri;
  double xvals2_0 = 0;
  double area_norm;  

 
  gridtype  = gridInqType(gridID1);
  gridsize1 = gridInqSize(gridID1);
 
  if ( cdoVerbose )
    {
      if ( ! circular )
        fprintf(stderr, "grid is not circular\n");
      else
        fprintf(stderr, "grid is circular\n");
    }

  if ( xinc < 1 || yinc < 1 )
    cdoAbort("xinc and yinc must not be smaller than 1");
  
  if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT || gridtype == GRID_CURVILINEAR )
    {
      nlon1 = gridInqXsize(gridID1);
      nlat1 = gridInqYsize(gridID1);

      nlon2 = nlon1/xinc;
      nlat2 = nlat1/yinc;
      if ( nlon1%xinc ) nlon2++;
      if ( nlat1%yinc ) nlat2++;
      gridsize2 = nlon2*nlat2;

      gridID2 = gridCreate(gridtype, gridsize2);
      gridDefXsize(gridID2, nlon2);
      gridDefYsize(gridID2, nlat2);
    }

  if ( xinc > nlon1 || yinc > nlat1 )
    cdoAbort("xinc and/or yinc exceeds gridsize");


  
  if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT )
    {
      xvals1 = (double *) malloc(nlon1*sizeof(double));
      yvals1 = (double *) malloc(nlat1*sizeof(double));
      xvals2 = (double *) malloc(nlon2*sizeof(double));
      yvals2 = (double *) malloc(nlat2*sizeof(double));
      gridInqXvals(gridID1, xvals1);
      gridInqYvals(gridID1, yvals1);

      if ( gridInqYbounds(gridID1, NULL) && gridInqXbounds(gridID1, NULL) )
        {
          grid1_corner_lon = (double *) malloc(2*nlon1*sizeof(double));
          grid1_corner_lat = (double *) malloc(2*nlat1*sizeof(double));
          grid2_corner_lon = (double *) malloc(2*nlon2*sizeof(double));
          grid2_corner_lat = (double *) malloc(2*nlat2*sizeof(double));
          gridInqXbounds(gridID1, grid1_corner_lon);
          gridInqYbounds(gridID1, grid1_corner_lat);
        }

      j = 0;
      for ( i = 0; i < nlon1; i += xinc )
        {
          i1 = i+(xinc-1);
          if ( i1 >= nlon1-1 ) i1 = nlon1-1; 
          xvals2[j] = xvals1[i] + (xvals1[i1] - xvals1[i])/2.;
          if ( isnan (xvals2[i]) )
            fprintf(stderr, "NaN in xvals2[%i] (i=%i i1=%i) xvals1[i1]=%7.5f xvals1[i]=%7.5f\n", j, i, i1, xvals1[i1], xvals1[i]);
          if ( grid2_corner_lon )
            {
              grid2_corner_lon[2*j] = grid1_corner_lon[2*i];
              grid2_corner_lon[2*j+1] = grid1_corner_lon[2*i1+1];
            }
          j++;        
        }
      j = 0;
      for ( i = 0; i < nlat1; i += yinc )
        {
          i1 = i+(yinc-1);
          if ( i1 >= nlat1-1 ) i1 = nlat1-1; 
          yvals2[j] = yvals1[i] + (yvals1[i1] - yvals1[i])/2;
          if ( isnan (yvals2[i]) )
            fprintf(stderr, "NaN in yvals2[%i] (i=%i i1=%i) yvals1[i1]=%7.5f yvals1[i]=%7.5f\n", j, i, i1, yvals1[i1], yvals1[i]);
          if ( grid2_corner_lat )
            {
              grid2_corner_lat[2*j] = grid1_corner_lat[2*i];
              grid2_corner_lat[2*j+1] = grid1_corner_lat[2*i1+1];
            }
          j++;
        }
      
      gridDefXvals(gridID2, xvals2);
      gridDefYvals(gridID2, yvals2);

      free(xvals2);
      free(yvals2);

      free(xvals1);
      free(yvals1);
      
      if ( grid2_corner_lon && grid2_corner_lat )
        {
          gridDefNvertex(gridID2, 2);
          gridDefXbounds(gridID2, grid2_corner_lon);
          gridDefYbounds(gridID2, grid2_corner_lat);
          free(grid2_corner_lon);
          free(grid2_corner_lat);

          free(grid1_corner_lon);
          free(grid1_corner_lat);
        }
    } /* if ( gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT ) */
  else if ( gridtype == GRID_CURVILINEAR )
    {
      int gridHasBounds = FALSE;
      
      if ( gridInqXbounds(gridID1, NULL) && gridInqYbounds(gridID1, NULL) )
        gridHasBounds = TRUE;
      
      xvals1 = (double *) malloc(nlon1*nlat1*sizeof(double));
      yvals1 = (double *) malloc(nlon1*nlat1*sizeof(double));
      xvals2 = (double *) malloc(nlon2*nlat2*sizeof(double));
      yvals2 = (double *) malloc(nlon2*nlat2*sizeof(double));
      
      gridInqXvals(gridID1, xvals1);
      gridInqYvals(gridID1, yvals1);
      
      if ( gridHasBounds )
        {
          grid1_corner_lon = (double *) malloc(4*nlon1*nlat1*sizeof(double));
          grid1_corner_lat = (double *) malloc(4*nlon1*nlat1*sizeof(double));
          grid2_corner_lon = (double *) malloc(4*nlon2*nlat2*sizeof(double));
          grid2_corner_lat = (double *) malloc(4*nlon2*nlat2*sizeof(double));
          
          gridInqXbounds(gridID1, grid1_corner_lon);
          gridInqYbounds(gridID1, grid1_corner_lat);
        }
      
      /* Process grid2 bounds */
      area_norm = xinc*yinc;
      for ( y2 = 0; y2 < nlat2; y2++ )
        {
          for ( x2 = 0; x2 < nlon2; x2++ )
            {
              g2_add = (y2*nlon2+x2); 
              on_up = on_lo = 360.; ox_up = ox_lo = -360.;
              an_ri = an_le =  90.; ax_ri = ax_le =  -90.; 
              
              for ( y1 = y2*yinc; y1 < yinc*(y2+1); y1++ )
                {
                  if ( y1 >= nlat1 ) use_y1 = nlat1-1;
                  else use_y1 = y1;
                  for ( x1 = x2*xinc; x1 < xinc*(x2+1) ; x1++ )
                    {                    
                      use_x1 = x1;
                      if ( x1 >= nlon1  ) 
                        {
                          if ( circular && use_y1 == y1 ) 
                            use_y1 -= 1;
                          else            
                            use_x1  = nlon1-1;
                        }
                      
                      g1_add= (use_y1*nlon1)+use_x1;                      

                      if ( y1 == y2*yinc && x1 == x2*xinc )
                        {
                          xvals2_0 = xvals1[g1_add];
                          xvals2[g2_add] = xvals1[g1_add]/area_norm;                        
                          if ( isnan(xvals2[g2_add]) ) fprintf(stderr, "isnan 1: xvals2_0 %7.5f xvals2[%i] %7.5f\n", xvals2_0, g2_add, xvals2[g2_add]);
                           yvals2[g2_add] = yvals1[g1_add]/area_norm;  
                        }
                      else if ( fabs(xvals1[g1_add] - xvals2_0) > 270. )
                        {
                          if ( (xvals1[g1_add] - xvals2_0) > 270. )
                            xvals2[g2_add] += (xvals1[g1_add]-360.)/area_norm;
                          else if ( ( xvals1[g1_add] - xvals2_0 ) < -270. )
                            xvals2[g2_add] += (xvals1[g1_add]+360.)/area_norm;
                          if ( isnan(xvals2[g2_add]) ) fprintf(stderr, "isnan 2\n");
                           yvals2[g2_add] += yvals1[g1_add]/area_norm;  
                        }
                      else                      
                        {
                          xvals2[g2_add] += xvals1[g1_add]/area_norm;                                                                    
                          if ( isnan(xvals2[g2_add]) ) fprintf(stderr, "isnan 3\n");
                          yvals2[g2_add] += yvals1[g1_add]/area_norm;  
                        }
                     

                      if ( isnan(xvals2[g2_add]) )
                        fprintf(stderr, "NaN in xvals2[%i]; last g1add %i (%7.5f) area-norm %7.5f\n", g2_add, g1_add, xvals1[g1_add], area_norm);                          
                      if ( isnan(yvals2[g2_add]) )
                        fprintf(stderr, "NaN in yvals2[%i]; last g1add %i (%7.5f) area-norm %7.5f\n", g2_add, g1_add, yvals1[g1_add], area_norm); 
                      //fprintf(stderr, "%9g %9g\n", xvals2[g2_add], yvals2[g2_add], xvals1[g1_add], yvals1[g1_add]);
                      if ( gridHasBounds )
                        {/* 
                         =====================================================================
                         ** O L D
                         */
                          if ( FALSE )
                            {
                              /* upper left cell */
                              if ( y1 == y2*yinc && x1 == x2*xinc)
                                {
                                  for ( corner = 0; corner < 4; corner++ )
                                    {
                                      add = 4*g1_add + corner;
                                      if ( grid1_corner_lon[add] < on_up )
                                        on_up = grid1_corner_lon[add];
                                      if (grid1_corner_lat[add] > ax_le )
                                        ax_le = grid1_corner_lat[add];
                                    }
                                }
                              
                              /* upper right cell */
                              if ( ( y1 == y2*yinc ) && ( x1 == (x2+1)*xinc - 1 ) )
                                {    
                                  for ( corner = 0; corner < 4; corner++ )
                                    {
                                      add = 4*g1_add + corner;
                                      if ( grid1_corner_lon[add] > ox_up )
                                        ox_up = grid1_corner_lon[add];
                                      if (grid1_corner_lat[add] > ax_ri )
                                        ax_ri = grid1_corner_lat[add];
                                    }
                                }
                              
                              /* lower right cell */
                              if ( ( y1 == (y2+1)*yinc -1 ) && (x1 == (x2+1)*xinc -1) )
                                {
                                  for ( corner = 0; corner < 4; corner ++ )
                                    {                             
                                  add = 4*g1_add + corner;
                                      if ( grid1_corner_lon[add] > ox_lo )
                                        ox_lo = grid1_corner_lon[add];
                                      if ( grid1_corner_lat[add] < an_ri )
                                        an_ri = grid1_corner_lat[add];
                                    }
                                }
                              
                              /* lower left cell */
                              if ( ( y1 == (y2+1)*yinc -1 ) && ( x1 == x2*xinc ) )
                                {    
                                  for ( corner = 0; corner < 4; corner ++ )
                                    {
                                      add = 4*g1_add + corner;
                                      if ( grid1_corner_lon[add] < on_lo )
                                    on_lo = grid1_corner_lon[add];
                                      if ( grid1_corner_lat[add] < an_le )
                                        an_le = grid1_corner_lat[add];
                                    }
                                }
                            }  /* O  L  D */
                          //=================================================================
                          else /* N  E  W */
                            {
                              int c_flag[4], corner2, g1_add2, g1_add3;
                              double lon, lat, lon2, lat2, lon3, lat3;
                              /* upper left cell */
                              if ( y1 == y2*yinc && x1 == x2*xinc)
                                {                                 
                                  c_flag[0] = c_flag[1] = c_flag[2] = c_flag[3] = 0;
                                  for ( corner = 0; corner < 4; corner++ )
                                    {                                      
                                      add = 4*g1_add + corner;
                                      lon = grid1_corner_lon[add];
                                      lat = grid1_corner_lat[add]; 
                                      g1_add2 = g1_add+1;
                                      if ( g1_add+nlon1 > gridsize1 ) 
                                        {
                                          cdoWarning("Can't find cell below upper left");
                                          continue; 
                                        }
                                      g1_add3 = g1_add+nlon1;
                                      for ( corner2 = 0; corner2 < 4; corner2++ )
                                        {                                          
                                          lon2 = grid1_corner_lon[4*g1_add2+corner2];
                                          lat2 = grid1_corner_lat[4*g1_add2+corner2];
                                          lon3 = grid1_corner_lon[4*g1_add3+corner2];
                                          lat3 = grid1_corner_lat[4*g1_add3+corner2];
                                          if ((lon2 == lon && lat2 == lat)  ||
                                              (lon3 == lon && lat3 == lat) )
                                            c_flag[corner] = 1;
                                        }                                                                                                                                                                      
                                    }
                                  for ( corner = 0; corner<4; corner++ )
                                    if ( !c_flag[corner] ) break;                                 
                                  on_up = grid1_corner_lon[4*g1_add + corner];
                                  ax_le = grid1_corner_lat[4*g1_add + corner];                                  
                                  if ( c_flag[0] + c_flag[1] + c_flag[2] + c_flag[3] < 3 )
                                    cdoWarning("found two matching corners!");                                  
                                }
                              
                              /* upper right cell */
                              if ( ( y1 == y2*yinc ) && ( x1 == (x2+1)*xinc - 1 ) )
                                {
                                  c_flag[0] = c_flag[1] = c_flag[2] = c_flag[3] = 0;
                                  for ( corner = 0; corner < 4; corner++ )
                                    {                                      
                                      add = 4*g1_add + corner;
                                      lon = grid1_corner_lon[add];
                                      lat = grid1_corner_lat[add]; 
                                      g1_add2 = g1_add-1;                                      
                                      if ( g1_add+nlon1 > gridsize1 ) 
                                        {
                                          cdoWarning("Can't find cell below upper left");
                                          continue; 
                                        }
                                      g1_add3 = g1_add+nlon1;
                                      for ( corner2 = 0; corner2 < 4; corner2++ )
                                        {                                          
                                          lon2 = grid1_corner_lon[4*g1_add2+corner2];
                                          lat2 = grid1_corner_lat[4*g1_add2+corner2];
                                          lon3 = grid1_corner_lon[4*g1_add3+corner2];
                                          lat3 = grid1_corner_lat[4*g1_add3+corner2];
                                          if ((lon2 == lon && lat2 == lat)  ||
                                              (lon3 == lon && lat3 == lat) )
                                            c_flag[corner] = 1;
                                        }                                                                                                                                                                      
                                    }                                 
                                  for ( corner = 0; corner<4; corner++ )
                                    if ( !c_flag[corner] ) break;                                 
                                  ox_up = grid1_corner_lon[4*g1_add + corner];
                                  ax_ri = grid1_corner_lat[4*g1_add + corner];                                  
                                  if ( c_flag[0] + c_flag[1] + c_flag[2] + c_flag[3] < 3 )
                                    cdoWarning("found two matching corners!");                                     
                                }
                            
                            
                              
                              /* lower right cell */
                              if ( ( y1 == (y2+1)*yinc -1 ) && (x1 == (x2+1)*xinc -1) )
                                {                                  
                                  c_flag[0] = c_flag[1] = c_flag[2] = c_flag[3] = 0;
                                  for ( corner = 0; corner < 4; corner++ )
                                    {                                      
                                      add = 4*g1_add + corner;
                                      lon = grid1_corner_lon[add];
                                      lat = grid1_corner_lat[add]; 
                                      g1_add2 = g1_add-1;
                                      if ( g1_add-nlon1 < 0 ) 
                                        {
                                          cdoWarning("Can't find cell above lower right left");
                                          continue; 
                                        }
                                      g1_add3 = g1_add-nlon1;
                                      for ( corner2 = 0; corner2 < 4; corner2++ )
                                        {                                          
                                          lon2 = grid1_corner_lon[4*g1_add2+corner2];
                                          lat2 = grid1_corner_lat[4*g1_add2+corner2];
                                          lon3 = grid1_corner_lon[4*g1_add3+corner2];
                                          lat3 = grid1_corner_lat[4*g1_add3+corner2];
                                          if ((lon2 == lon && lat2 == lat)  ||
                                              (lon3 == lon && lat3 == lat) )
                                            c_flag[corner] = 1;
                                        }                                                                                                                                                                      
                                    }                                  
                                  for ( corner = 0; corner<4; corner++ )
                                    if ( !c_flag[corner] ) break;                                 
                                  ox_lo = grid1_corner_lon[4*g1_add + corner];
                                  an_ri = grid1_corner_lat[4*g1_add + corner];                                  
                                  if ( c_flag[0] + c_flag[1] + c_flag[2] + c_flag[3] < 3 )
                                    cdoWarning("found two matching corners!");        
                                }
                              
                              /* lower left cell */
                              if ( ( y1 == (y2+1)*yinc -1 ) && ( x1 == x2*xinc ) )
                                {    
                                  c_flag[0] = c_flag[1] = c_flag[2] = c_flag[3] = 0;
                                  for ( corner = 0; corner < 4; corner++ )
                                    {                                      
                                      add = 4*g1_add + corner;
                                      lon = grid1_corner_lon[add];
                                      lat = grid1_corner_lat[add]; 
                                      g1_add2 = g1_add+1;
                                      if ( g1_add-nlon1 < 0 ) 
                                        {
                                          cdoWarning("Can't find cell above lower right left");
                                          continue; 
                                        }
                                      g1_add3 = g1_add-nlon1;
                                      for ( corner2 = 0; corner2 < 4; corner2++ )
                                        {                                          
                                          lon2 = grid1_corner_lon[4*g1_add2+corner2];
                                          lat2 = grid1_corner_lat[4*g1_add2+corner2];
                                          lon3 = grid1_corner_lon[4*g1_add3+corner2];
                                          lat3 = grid1_corner_lat[4*g1_add3+corner2];
                                          if ((lon2 == lon && lat2 == lat)  ||
                                              (lon3 == lon && lat3 == lat) )
                                            c_flag[corner] = 1;
                                        }                                                                                                                                                                      
                                    }                                 
                                  for ( corner = 0; corner<4; corner++ )
                                    if ( !c_flag[corner] ) break;                                 
                                  on_lo = grid1_corner_lon[4*g1_add + corner];
                                  an_le = grid1_corner_lat[4*g1_add + corner];                                  
                                  if ( c_flag[0] + c_flag[1] + c_flag[2] + c_flag[3] < 3 )
                                    cdoWarning("found two matching corners!");                                      
                                }
                            }/* else to if(fast) */
                        }/* if ( gridHasBounds) */                  
                    }/* for ( x1 = x2*xinc; x1 < xinc*(x2+1) ; x1++ ) */
                }/* for ( y1 = y2*yinc; y1 < yinc*(y2+1); y1++ ) */
 
              
              if ( gridHasBounds )
                {
                  /* upper left corner */
                  grid2_corner_lon[4*g2_add+3] = on_up;
                  grid2_corner_lat[4*g2_add+3] = ax_le;
                  /* upper right corner */
                  grid2_corner_lon[4*g2_add+2] = ox_up;
                  grid2_corner_lat[4*g2_add+2] = ax_ri;
                  /* lower right corner */
                  grid2_corner_lon[4*g2_add+1] = ox_lo;
                  grid2_corner_lat[4*g2_add+1] = an_ri;
                  /* lower left corner */
                  grid2_corner_lon[4*g2_add+0] = on_lo;
                  grid2_corner_lat[4*g2_add+0] = an_le;
                }
              
              while ( xvals2[g2_add] >  180. ) xvals2[g2_add] -= 360.;
              while ( xvals2[g2_add] < -180. ) xvals2[g2_add] += 360.;
            } /* for ( x2 = 0; x2 < nlon2; x2++ ) */
        } /* for ( y2 = 0; y2 < nlat2; y2++ ) */
      
      gridDefXvals(gridID2, xvals2);
      gridDefYvals(gridID2, yvals2);

      free(xvals2);
      free(yvals2);

      free(xvals1);
      free(yvals1);

      if ( gridHasBounds )
        {
          gridDefNvertex(gridID2, 4);
          gridDefXbounds(gridID2, grid2_corner_lon);
          gridDefYbounds(gridID2, grid2_corner_lat);
          
          free(grid2_corner_lon);
          free(grid2_corner_lat);
          
          free(grid1_corner_lon);
          free(grid1_corner_lat);
        }
    } /* else if ( gridtype == GRID_CURVILINEAR ) */
  else
    {
      cdoAbort("Unsupported grid: %s", gridNamePtr(gridtype));
    }
  
  return gridID2;
}


static
void gridboxstat(FIELD *field1, FIELD *field2, int xinc, int yinc, int statfunc)
{
  static char func[] = "gridboxstat";
  int nlon1, nlat1;
  int nlon2, nlat2;
  int ilat, ilon;
  int gridID1, gridID2;
  int nmiss;
  double **xfield1;
  double **weight = NULL;
  double *array1, *array2;
  double missval;
  int i, j, ii, jj;
  double **xfield2;
  FIELD field;
  int isize;
  int useWeight = FALSE;

  if ( field1->weight ) useWeight = TRUE;

  field.size    = xinc*yinc;
  field.ptr     = (double *) malloc(field.size*sizeof(double));
  field.weight  = NULL;
  if ( useWeight )
    field.weight  = (double *) malloc(field.size*sizeof(double));
  field.missval = field1->missval;
  field.nmiss   = 0;
  
  gridID1 = field1->grid;
  gridID2 = field2->grid;
  array1  = field1->ptr;
  array2  = field2->ptr;
  missval = field1->missval;

  nlon1 = gridInqXsize(gridID1);
  nlat1 = gridInqYsize(gridID1);

  nlon2 = gridInqXsize(gridID2);
  nlat2 = gridInqYsize(gridID2);

  xfield1 = (double **) malloc(nlat1*sizeof(double *));
  if ( useWeight )
    weight  = (double **) malloc(nlat1*sizeof(double *));

  for ( ilat = 0; ilat < nlat1; ilat++ )
    xfield1[ilat] = array1 + ilat*nlon1;

  if ( useWeight )
    for ( ilat = 0; ilat < nlat1; ilat++ )
      weight[ilat] = field1->weight + ilat*nlon1;


  xfield2 = (double **) malloc(nlat2 * sizeof(double *));

  for ( ilat = 0; ilat < nlat2; ilat++ )
    xfield2[ilat] = array2 + ilat*nlon2;

  for ( ilat = 0; ilat < nlat2; ilat++ )
    for ( ilon = 0; ilon < nlon2; ilon++ )
      {
        isize = 0;
        field.nmiss = 0;
        for ( j = 0; j < yinc; ++j )
          {
            jj = ilat*yinc+j;
            if ( jj >= nlat1 ) break;
            for ( i = 0; i < xinc; ++i )
              {
                ii = ilon*xinc+i;
                if ( ii >= nlon1 ) break;
                field.ptr[isize] = xfield1[jj][ii];
                if ( useWeight ) field.weight[isize] = weight[jj][ii];
                if ( DBL_IS_EQUAL(field.ptr[isize], field.missval) ) field.nmiss++;
                isize++;
              }
          }
        
        field.size = isize;
        xfield2[ilat][ilon] = fldfun(field, statfunc);
      }
  
  nmiss = 0;
  for ( i = 0; i < nlat2*nlon2; i++ )
    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;
  
  field2->nmiss = nmiss;
  
  free(xfield2);
  free(xfield1);
  free(field.ptr);
  if ( useWeight )
    free(field.weight);
}


void *Gridboxstat(void *argument)
{
  static char func[] = "Gridboxstat";
  int operatorID;
  int operfunc;
  int streamID1, streamID2;
  int vlistID1, vlistID2;
  int wstatus = FALSE;
  int code = 0, oldcode = 0;
  int index, ngrids;
  int recID, nrecs;
  int tsID, varID, levelID;
  int needWeights = FALSE;
  int gridID1, gridID2;
  int gridsize1, gridsize2;
  FIELD field1, field2;
  int taxisID1, taxisID2;
  int xinc, yinc;

  cdoInitialize(argument);

  operatorInputArg("xinc, yinc");
  operatorCheckArgc(2);
  xinc = atoi(operatorArgv()[0]);
  yinc = atoi(operatorArgv()[1]);

  cdoOperatorAdd("gridboxmin",  func_min,  0, NULL);
  cdoOperatorAdd("gridboxmax",  func_max,  0, NULL);
  cdoOperatorAdd("gridboxsum",  func_sum,  0, NULL);
  cdoOperatorAdd("gridboxmean", func_mean, 0, NULL);
  cdoOperatorAdd("gridboxavg",  func_avg,  0, NULL);
  cdoOperatorAdd("gridboxvar",  func_var,  0, NULL);
  cdoOperatorAdd("gridboxstd",  func_std,  0, NULL);

  operatorID = cdoOperatorID();
  operfunc = cdoOperatorFunc(operatorID);

  if ( operfunc == func_mean || operfunc == func_avg ||
       operfunc == func_var  || operfunc == func_std )
    needWeights = TRUE;

  streamID1 = streamOpenRead(cdoStreamName(0));
  if ( streamID1 < 0 ) cdiError(streamID1, "Open failed on %s", cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  ngrids = vlistNgrids(vlistID1);
  if ( ngrids > 1 )  cdoAbort("Too many different grids!");

  gridID1 = vlistGrid(vlistID1, 0);
  if ( gridInqType(gridID1) == GRID_GAUSSIAN_REDUCED )
    cdoAbort("Gaussian reduced grid found. Use option -R to convert it to a regular grid!");

  gridID2 = genBoxGrid(gridID1, xinc, yinc);
  for ( index = 0; index < ngrids; index++ )
    vlistChangeGridIndex(vlistID2, index, gridID2);

  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID2 < 0 ) cdiError(streamID2, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID2, vlistID2);

  gridsize1 = gridInqSize(gridID1);
  field1.ptr    = (double *) malloc(gridsize1*sizeof(double));
  field1.weight = NULL;
  if ( needWeights )
    field1.weight = (double *) malloc(gridsize1*sizeof(double));

  gridsize2 = gridInqSize(gridID2);
  field2.ptr    = (double *) malloc(gridsize2*sizeof(double));
  field2.weight = NULL;

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID1, &varID, &levelID);
          streamReadRecord(streamID1, field1.ptr, &field1.nmiss);

          field1.grid = vlistInqVarGrid(vlistID1, varID);
          field1.size = gridInqSize(field1.grid);
          field1.missval = vlistInqVarMissval(vlistID1, varID);
          
          field2.grid = gridID2;
          field2.size = gridsize2;
          field2.missval = field1.missval;

          if ( needWeights )
            {
              wstatus = gridWeights(field1.grid, field1.weight);
            }
          code = vlistInqVarCode(vlistID1, varID);
          if ( wstatus != 0 && tsID == 0 && code != oldcode )
            cdoWarning("Using constant grid cell area weights for code %d!", oldcode=code);
          
          gridboxstat(&field1, &field2, xinc, yinc, operfunc);
          
          streamDefRecord(streamID2, varID,  levelID);
          streamWriteRecord(streamID2, field2.ptr, field2.nmiss);
        }
      tsID++;
    }
  
  streamClose(streamID2);
  streamClose(streamID1);
  
  if ( field1.ptr )    free(field1.ptr);
  if ( field1.weight ) free(field1.weight);
  
  if ( field2.ptr )    free(field2.ptr);
  if ( field2.weight ) free(field2.weight);
  
  cdoFinish();
  
  return (0);
}
