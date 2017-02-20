#include <cdi.h>
#include "cdo_int.h"
#include "grid.h"

/*
@Function  define_sample_grid
@Title     Define a sampled grid of another grid

@Prototype int define_sample_grid(int gridSrcID, int sampleFactor)
@Parameter
    @Item  gridSrcID       Source grid
    @Item  sampleFactor    sampleFactor; typically 2,3,4 ...

@Description
The function @func{define_sample_grid} defines a sampled grid of another grid

@EndFunction
*/
int cdo_define_sample_grid(int gridSrcID, int sampleFactor)
{
/* Example of horizontal grids (Harmonie HARM36_L25):
            #
            # gridID 2
            #
            gridtype  = lcc
            gridsize  = 622521
            xsize     = 789
            ysize     = 789
            originLon = -7.89
            originLat = 42.935
            lonParY   = 0
            lat1      = 52.5
            lat2      = 52.5
            xinc      = 2500
            yinc      = 2500
            projection = northpole
=>   RESULT:
            #
            # gridID 2
            #
            gridtype  = lcc
            gridsize  = 156025
            xsize     = 395
            ysize     = 395
            originLon = -7.89
            originLat = 42.935
            lonParY   = 0
            lat1      = 52.5
            lat2      = 52.5
            xinc      = 5000
            yinc      = 5000
            projection = northpole
*/
  if ( cdoDebugExt )
    cdoPrint("cdo_define_sample_grid(gridSrcID=%d, sampleFactor=%d) ...", gridSrcID, sampleFactor);

  int gridtype = gridInqType(gridSrcID);
  if ( ! (gridtype == GRID_GAUSSIAN || gridtype == GRID_LONLAT || gridtype == GRID_PROJECTION ||
          gridtype == GRID_CURVILINEAR || gridtype == GRID_GENERIC || gridtype == GRID_LCC) )
    cdoAbort("Unsupported gridtype: %s", gridNamePtr(gridtype));
    
  int gridXsize = gridInqXsize(gridSrcID);
  int gridYsize = gridInqYsize(gridSrcID);

  if ( (sampleFactor<1) || (gridXsize<1) || (gridYsize<1) || (sampleFactor > (gridXsize/4) ) || (sampleFactor > (gridYsize/4)) )
    cdoAbort("%s(): Unsupported sampleFactor (%d)! Note that: gridXsize = %d, gridYsize = %d",
             __func__, sampleFactor, gridXsize, gridYsize);

  if ( cdoDebugExt>20 ) cdo_print_grid(gridSrcID, 1);

  int xsize = (gridXsize + (sampleFactor-1)) / sampleFactor; // HARM36_L25: (789 + 2-1) / 2 = 395
  int ysize = (gridYsize + (sampleFactor-1)) / sampleFactor;

  int gridID_sampled = gridCreate(gridtype, xsize*ysize);

  // TODO
  /*
    grid_sampled->scanningMode          = grid_src->scanningMode;
    grid_sampled->iScansNegatively      = grid_src->iScansNegatively;
    grid_sampled->jScansPositively      = grid_src->jScansPositively;
    grid_sampled->jPointsAreConsecutive = grid_src->jPointsAreConsecutive;
    grid_sampled->uvRelativeToGrid      = grid_src->uvRelativeToGrid;
  */
  gridDefXsize(gridID_sampled, xsize);
  gridDefYsize(gridID_sampled, ysize);

  gridDefNP(gridID_sampled, gridInqNP(gridSrcID));
  gridDefPrec(gridID_sampled, gridInqPrec(gridSrcID));

  grid_copy_attributes(gridSrcID, gridID_sampled);
  
  if ( gridtype == GRID_PROJECTION ) grid_copy_mapping(gridSrcID, gridID_sampled);

  if ( gridtype == GRID_LCC )
    {
      double originLon, originLat, lonParY, lat1, lat2, xinc, yinc;
      int projflag, scanflag;

      gridInqParamLCC(gridSrcID, &originLon, &originLat, &lonParY, &lat1, &lat2, &xinc, &yinc, &projflag, &scanflag);

      xinc *= sampleFactor;
      yinc *= sampleFactor;
        
      gridDefParamLCC(gridID_sampled, originLon, originLat, lonParY, lat1, lat2, xinc, yinc, projflag, scanflag);
    }
  else if ( gridInqXvals(gridSrcID, NULL) && gridInqYvals(gridSrcID, NULL) )
    {
      if ( gridtype == GRID_CURVILINEAR )
        {
          double *vals = (double *) Malloc(gridXsize*gridYsize*sizeof(double));
          gridInqXvals(gridSrcID, vals);
          double *pvals = vals;
          for ( int j = 0; j < gridYsize; j += sampleFactor )
            for ( int i = 0; i < gridXsize; i += sampleFactor )
              *pvals++ = vals[j*gridXsize+i];
          gridDefXvals(gridID_sampled, vals);
          
          gridInqYvals(gridSrcID, vals);
          pvals = vals;
          for ( int j = 0; j < gridYsize; j += sampleFactor )
            for ( int i = 0; i < gridXsize; i += sampleFactor )
              *pvals++ = vals[j*gridXsize+i];
          gridDefYvals(gridID_sampled, vals);
          Free(vals);
        }
      else
        {
          double *xvals = (double *) Malloc(gridXsize*sizeof(double));
          gridInqXvals(gridSrcID, xvals);
          for ( int i = 0, j = 0; i < gridXsize; i += sampleFactor ) xvals[j++] = xvals[i];
          gridDefXvals(gridID_sampled, xvals);
          Free(xvals);

          double *yvals = (double *) Malloc(gridYsize*sizeof(double));
          gridInqYvals(gridSrcID, yvals);
          for ( int i = 0, j = 0; i < gridYsize; i += sampleFactor ) yvals[j++] = yvals[i];
          gridDefYvals(gridID_sampled, yvals);
          Free(yvals);
        }
    }

  if ( cdoDebugExt>20 )
    {
      cdoPrint("cdo SampleGrid: define_sample_grid(): ");
      cdo_print_grid(gridID_sampled, 1);
    }

  return gridID_sampled;
}


/*
@Function  define_subgrid_grid
@Title     Define a sub-grid of another grid (LCC)

@Prototype int define_subgrid_grid(int gridIDsrc, int subI0, int subI1, int subJ0, int subJ1)
@Parameter
    @Item  gridSrcID                    Source grid
    @Item  subI0,subI1, subJ0, subJ1    Sub-grid indices

@Description
The function @func{define_subgrid_grid} defines a sub-grid of another grid (LCC)

@EndFunction
*/
int cdo_define_subgrid_grid(int gridSrcID, int subI0, int subI1, int subJ0, int subJ1)
{
/* Example of horizontal grids (Harmonie HARM36_L25):
            #
            # gridID 2
            #
            gridtype  = lcc
            gridsize  = 622521
            xsize     = 789
            ysize     = 789
            originLon = -7.89
            originLat = 42.935
            lonParY   = 0
            lat1      = 52.5
            lat2      = 52.5
            xinc      = 2500
            yinc      = 2500
            projection = northpole
=>   RESULT:
            #
            # gridID 2
            #
            gridtype  = lcc
            gridsize  = 156025
            xsize     = 350
            ysize     = 350
            originLon = ...
            originLat = ...
            lonParY   = 0
            lat1      = 52.5
            lat2      = 52.5do SampleGrid: define_subgrid_grid
            xinc      = 2500
            yinc      = 2500
            projection = northpole
*/
  if ( cdoDebugExt )
    cdoPrint("%s(gridSrcID=%d, (subI0,subI1,subJ0,subJ1) = (%d,%d,%d,%d) ...",
             __func__, gridSrcID, subI0,subI1, subJ0, subJ1 );

  int gridXsize = gridInqXsize(gridSrcID);
  int gridYsize = gridInqYsize(gridSrcID);
  int maxIndexI = gridXsize-1;
  int maxIndexJ = gridYsize-1;

  if ( (subI0<0) || (subI0>maxIndexI) ||
       (subI1<=subI0) || (subI1>maxIndexI) ||
       (subJ0<0) || (subJ0>maxIndexJ) ||
       (subJ1<=subJ0) || (subJ1>maxIndexJ) )
    cdoAbort("%s() Incorrect subgrid specified!  (subI0,subI1,subJ0,subJ1) =(%d,%d,%d,%d) Note that: gridXsize = %d, gridYsize = %d", __func__, subI0,subI1, subJ0, subJ1, gridXsize, gridYsize);

  int gridtype = gridInqType(gridSrcID);
  if ( gridtype != GRID_LCC )
    cdoAbort("%s() Error; Only LCC grid is supported; use selindexbox!", __func__);

  double originLon, originLat, lonParY, lat1, lat2, xinc, yinc;
  int projflag, scanflag;

  gridInqParamLCC(gridSrcID, &originLon, &originLat, &lonParY, &lat1, &lat2, &xinc, &yinc, &projflag, &scanflag);

  if ( cdoDebugExt>20 ) cdo_print_grid(gridSrcID, 1);

  if ( cdoDebugExt )
    {
      cdoPrint("%s() Original LCC grid:", __func__);
      cdoPrint("grid Xsize   %d, grid Ysize   %d", gridXsize, gridYsize);
      cdoPrint("originLon %4.3f, originLat %4.3f", originLon, originLat);
      cdoPrint("grid Xinc   %4.3f, grid Yinc   %4.3f", xinc, yinc);
    }
  
  int gridIDcurvl = gridToCurvilinear(gridSrcID, 1);

  originLon = gridInqXval(gridIDcurvl, 0);
  originLat = gridInqYval(gridIDcurvl, 0);

  if ( cdoDebugExt )
    {
      cdoPrint("%s() Original LCC grid as curvilinear (with lats-lons computed):", __func__);
      cdoPrint("grid Xsize   %d, grid Ysize   %d", gridInqXsize(gridIDcurvl), gridInqYsize(gridIDcurvl));
      cdoPrint("grid Xfirst  %4.3f, grid Yfirst  %4.3f", gridInqXval(gridIDcurvl, 0), gridInqYval(gridIDcurvl, 0));
      cdoPrint("grid Xlast   %4.3f, grid Ylast   %4.3f", gridInqXval(gridIDcurvl, gridInqSize(gridIDcurvl) -1), gridInqYval(gridIDcurvl, gridInqSize(gridIDcurvl) -1));
      cdoPrint("originLon %4.3f, originLat %4.3f", originLon, originLat);
    }


  int xsize = subI1 - subI0 + 1;
  int ysize = subJ1 - subJ0 + 1;

  int gridID_sampled = gridCreate(gridtype, xsize*ysize);
  
  // TODO
  /*
    grid_sampled->scanningMode          = grid_src->scanningMode;
    grid_sampled->iScansNegatively      = grid_src->iScansNegatively;
    grid_sampled->jScansPositively      = grid_src->jScansPositively;
    grid_sampled->jPointsAreConsecutive = grid_src->jPointsAreConsecutive;
    grid_sampled->uvRelativeToGrid      = grid_src->uvRelativeToGrid;
  */

  gridDefXsize(gridID_sampled, xsize);
  gridDefYsize(gridID_sampled, ysize);

  gridDefNP(gridID_sampled, gridInqNP(gridSrcID));
  gridDefPrec(gridID_sampled, gridInqPrec(gridSrcID));

  grid_copy_attributes(gridSrcID, gridID_sampled);
  
  if ( gridtype == GRID_PROJECTION ) grid_copy_mapping(gridSrcID, gridID_sampled);

  originLon = gridInqXval(gridIDcurvl, subJ0*gridXsize + subI0 );
  originLat = gridInqYval(gridIDcurvl, subJ0*gridXsize + subI0 );

  if ( cdoDebugExt )
    {
      cdoPrint("%s()  Sub-grid:", __func__);
      cdoPrint("grid Xsize   %d, grid Ysize   %d", gridInqXsize(gridID_sampled), gridInqYsize(gridID_sampled));
      cdoPrint("originLon %4.3f, originLat %4.3f", originLon, originLat);
    }

  gridDefParamLCC(gridID_sampled, originLon, originLat, lonParY, lat1, lat2, xinc, yinc, projflag, scanflag);
    
  gridDestroy(gridIDcurvl);

  if ( cdoDebugExt>20 )
    {
      cdoPrint("%s(): ", __func__);
      cdo_print_grid(gridID_sampled, 1);
    }
    
  return gridID_sampled;
}
