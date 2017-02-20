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
int cdo_define_subgrid_grid(int gridSrcID, int gridIDcurvl, int subI0, int subI1, int subJ0, int subJ1)
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
            lat2      = 52.5
            xinc      = 2500
            yinc      = 2500
            projection = northpole
*/
    if ( cdoDebugExt )
        printf("cdo SampleGrid: define_subgrid_grid(gridSrcID=%d, (subI0,subI1,subJ0,subJ1) =(%d,%d,%d,%d) ...\n",
               gridSrcID, subI0,subI1, subJ0, subJ1 );


    int gridXsize = gridInqXsize(gridSrcID);
    int gridYsize = gridInqYsize(gridSrcID);
    int maxIndexI = gridXsize-1;
    int maxIndexJ = gridYsize-1;

    if ( (subI0<0) || (subI0>maxIndexI) ||
         (subI1<=subI0) || (subI1>maxIndexI) ||
         (subJ0<0) || (subJ0>maxIndexJ) ||
         (subJ1<=subJ0) || (subJ1>maxIndexJ) )
      cdoAbort("cdo_define_subgrid_grid() Incorrect subgrid specified!  (subI0,subI1,subJ0,subJ1) =(%d,%d,%d,%d) Note that: gridXsize = %d, gridYsize = %d", subI0,subI1, subJ0, subJ1, gridXsize, gridYsize);


    // TODO
    cdoAbort("cdo_define_subgrid_grid not implement!");
    /*
    grid_t *grid_src;
    grid_t *grid_sampled;
    int gridID_sampled;

    double originLon, originLat;

    grid_src =  grid_to_pointer(gridSrcID);
    grid_check_ptr(gridSrcID, grid_src);

    if ( grid_src->type != GRID_LCC )
        Error("cdo SampleGrid: define_subgrid_grid() Error; Only LCC grid is supported; use selindexbox");


    double lonParY; double lat1; double lat2; double xinc; double yinc; int projflag; int scanflag;
    gridInqLCC(gridSrcID, &originLon, &originLat, &lonParY, &lat1, &lat2, &xinc, &yinc, &projflag, &scanflag);

    if ( cdoDebugExt>20 ) {
        gridPrint(gridSrcID,1,0);
    }

    if (cdoDebugExt)
    {
        Message("cdo SampleGrid: define_subgrid_grid() Original LCC grid:");
        Message("grid Xsize   %d, grid Ysize   %d", gridInqXsize(gridSrcID), gridInqYsize(gridSrcID));
        Message("originLon %4.3f, originLat %4.3f", originLon, originLat);
        Message("grid Xinc   %4.3f, grid Yinc   %4.3f", xinc, yinc);
    }
    originLon = gridInqXval(gridIDcurvl, 0);
    originLat = gridInqYval(gridIDcurvl, 0);

    if (cdoDebugExt)
    {
        Message("\ncdo SampleGrid: define_subgrid_grid() Original LCC grid as curvilinear (with lats-lons computed):");
        Message("grid Xsize   %d, grid Ysize   %d", gridInqXsize(gridIDcurvl), gridInqYsize(gridIDcurvl));
        Message("grid Xfirst  %4.3f, grid Yfirst  %4.3f", gridInqXval(gridIDcurvl, 0), gridInqYval(gridIDcurvl, 0));
        Message("grid Xlast   %4.3f, grid Ylast   %4.3f", gridInqXval(gridIDcurvl, gridInqSize(gridIDcurvl) -1), gridInqYval(gridIDcurvl, gridInqSize(gridIDcurvl) -1));
        Message("originLon %4.3f, originLat %4.3f", originLon, originLat);
    }


    gridID_sampled = gridDuplicate(gridSrcID);
    grid_sampled = grid_to_pointer(gridID_sampled);
    grid_check_ptr(gridID_sampled, grid_sampled);

    grid_sampled->scanningMode          = grid_src->scanningMode;
    grid_sampled->iScansNegatively      = grid_src->iScansNegatively;
    grid_sampled->jScansPositively      = grid_src->jScansPositively;
    grid_sampled->jPointsAreConsecutive = grid_src->jPointsAreConsecutive;
    grid_sampled->uvRelativeToGrid      = grid_src->uvRelativeToGrid;

    grid_sampled->xsize   = subI1 - subI0 + 1;
    grid_sampled->ysize   = subJ1 - subJ0 + 1;
    grid_sampled->size    = grid_sampled->xsize * grid_sampled->ysize;
    grid_sampled->xinc   = grid_src->xinc;
    grid_sampled->yinc   = grid_src->yinc;


    originLon = gridInqXval(gridIDcurvl, subJ0*gridInqXsize(gridIDcurvl) + subI0 );
    originLat = gridInqYval(gridIDcurvl, subJ0*gridInqXsize(gridIDcurvl) + subI0 );

    if (cdoDebugExt)
    {
        Message("\ncdo SampleGrid: define_subgrid_grid()  Sub-grid:");
        Message("grid Xsize   %d, grid Ysize   %d", gridInqXsize(gridID_sampled), gridInqYsize(gridID_sampled));
        Message("originLon %4.3f, originLat %4.3f", originLon, originLat);
    }

    grid_sampled->lcc_originLon   = originLon;
    grid_sampled->lcc_originLat   = originLat;
    grid_sampled->lcc_lat1   = grid_src->lcc_lat1;
    grid_sampled->lcc_lat2   = grid_src->lcc_lat2;
    grid_sampled->lcc_xinc   = grid_src->lcc_xinc;
    grid_sampled->lcc_yinc   = grid_src->lcc_yinc;
    grid_sampled->lcc_projflag   = grid_src->lcc_projflag;
    grid_sampled->lcc_scanflag   = grid_src->lcc_scanflag;
    grid_sampled->lcc_defined   = grid_src->lcc_defined;

    if ( grid_sampled->type == GRID_LCC )
      gridDefLCC(gridID_sampled, grid_sampled->lcc_originLon, grid_sampled->lcc_originLat, grid_sampled->lcc_lonParY,
             grid_sampled->lcc_lat1, grid_sampled->lcc_lat2, grid_sampled->lcc_xinc, grid_sampled->lcc_yinc,
             grid_sampled->lcc_projflag, grid_sampled->lcc_scanflag);


    if ( cdoDebugExt>20 )
    {
        printf("cdo SampleGrid: define_subgrid_grid(): \n");
        gridPrint(gridID_sampled, 1,0);
    }
    return gridID_sampled;
    */
}
