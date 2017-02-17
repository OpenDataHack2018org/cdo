#include <cdi.h>
#include "cdo_int.h"
#include "grid.h"

extern int cdoDebugExt;
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
       cdoPrint("cdo_define_sample_grid(gridSrcID=%d, sampleFactor=%d) ...\n",gridSrcID, sampleFactor);

    int gridtype = gridInqType(gridSrcID);
    int gridXsize = gridInqXsize(gridSrcID);
    int gridYsize = gridInqYsize(gridSrcID);

    if ( (sampleFactor<1) || (gridXsize<1) || (gridYsize<1) || (sampleFactor > (gridXsize/4) ) || (sampleFactor > (gridYsize/4)) )
        cdoAbort("cdo_define_sample_grid() Unsupported sampleFactor (%d)! Note that: gridXsize = %d, gridYsize = %d",
                 sampleFactor, gridXsize, gridYsize);

    // TODO if ( cdoDebugExt>20 )  gridPrint(gridSrcID,1,0);

    //const double *xvals   = gridInqXvalsPtr(gridSrcID);
    //const double *yvals   = gridInqYvalsPtr(gridSrcID);
    /*
    double xfirst = gridInqXval(gridSrcID,0);   // staggered grid of u-wind
    double yfirst = gridInqYval(gridSrcID,0);
    double xlast  = gridInqXval(gridSrcID, gridXsize-1);
    double ylast  = gridInqYval(gridSrcID, gridYsize-1);
    double xinc     = gridInqXinc(gridSrcID);
    double yinc     = gridInqYinc(gridSrcID);
    */

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

    // for the case of Lambert projection ...
    // TODO
    /*
    grid_sampled->lcc_xinc   = grid_src->lcc_xinc * sampleFactor;
    grid_sampled->lcc_yinc   = grid_src->lcc_yinc * sampleFactor;

    grid_sampled->xinc   = grid_src->xinc * sampleFactor;
    grid_sampled->yinc   = grid_src->yinc * sampleFactor;
    */

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

    // TODO
    /*
    if ( grid_sampled->type == GRID_LCC )
        gridDefLCC( gridID_sampled, grid_sampled->lcc_originLon, grid_sampled->lcc_originLat, grid_sampled->lcc_lonParY,
                    grid_sampled->lcc_lat1, grid_sampled->lcc_lat2, grid_sampled->lcc_xinc, grid_sampled->lcc_yinc,
                    grid_sampled->lcc_projflag, grid_sampled->lcc_scanflag);
    */
    if ( cdoDebugExt>20 )
      {
        printf("cdo SampleGrid: define_sample_grid(): \n");
        // TODO gridPrint(gridID_sampled, 1,0);
      }

    return gridID_sampled;
}
