/*
HIRLAM extensions ..
*/

/*
   This module "SampleGrid" contains the following operators:

    samplegrid      Resample current grid with given factor, typically 2 (which will half the resolution);
                    tested on curvilinear and LCC grids;
    subgrid         Similar to selindexbox but this operator works for LCC grids (tested on HARMONIE NWP model).
*/

#include <ctype.h>
#include "cdo.h"
#include "cdi.h"
#include "cdo_int.h"
#include "grid.h"

#include "griddes.h"
#include "pstream.h"
#include "specspace.h"
#include "list.h"
#include "math.h"

#ifdef HIRLAM_EXTENSIONS

extern int cdoDebugExt; // defined in cdo.c


static
void sampleData(int nwpv, double *array1, int gridID1, double *array2, int gridID2, int resampleFactor )
{
    long nlon1, nlat1;
    long nlon2, nlat2;
    long ilat1, ilon1;
    long ilat2, ilon2;

    nlon1 = gridInqXsize(gridID1);
    nlat1 = gridInqYsize(gridID1);

    nlon2 = gridInqXsize(gridID2);
    nlat2 = gridInqYsize(gridID2);

    if ( cdoDebugExt >= 100) cdoPrint("sampleData():: (nlon1: %d; nlat1: %d) => (nlon2: %d; nlat2: %d); gridID1: %d; gridID2: %d; resampleFactor: %d)",nlon1,nlat1, nlon2,nlat2, gridID1, gridID2, resampleFactor);

    if ( nwpv == 1 )
    {
        for ( ilat1 = 0; ilat1 < nlat1; ilat1+=resampleFactor )
        {
            for ( ilon1 = 0; ilon1 < nlon1; ilon1+=resampleFactor )
                *array2++ = array1[ilat1*nlon1 + ilon1];
        }
    }
    else
    if ( nwpv == 2 ) // complex numbers ... unsupported yet ...
    {
    /*
    for ( ilat = lat1; ilat <= lat2; ilat++ )
    {
      for ( ilon = lon21; ilon <= lon22; ilon++ )
        {
          *array2++ = array1[ilat*nlon1*2 + ilon*2];
          *array2++ = array1[ilat*nlon1*2 + ilon*2+1];
        }
      for ( ilon = lon11; ilon <= lon12; ilon++ )
        {
          *array2++ = array1[ilat*nlon1*2 + ilon*2];
          *array2++ = array1[ilat*nlon1*2 + ilon*2+1];
        }
    } */
    }
}

static
void cropData(int nwpv, double *array1, int gridID1, double *array2, int gridID2, int subI0, int subI1, int  subJ0, int  subJ1 )
{
    long nlon1, nlat1;
    long nlon2, nlat2;
    long ilat1, ilon1;
    long ilat2, ilon2;
    long array2Idx=0;

    nlon1 = gridInqXsize(gridID1);
    nlat1 = gridInqYsize(gridID1);

    nlon2 = gridInqXsize(gridID2);
    nlat2 = gridInqYsize(gridID2);
    long rowLen;
    rowLen =  subI1 - subI0 +1; // must be same as   nlon1

    if (rowLen!= nlon2)
        cdoAbort("cropData() rowLen!= nlon2 [%d != %d]", rowLen, nlon2);


    if ( nwpv == 1 )
    {
        if ( cdoDebugExt>=10 ) cdoPrint("cropData(%d,%d,%d,%d) ...\n",subI0,subI1, subJ0, subJ1 );

        for ( ilat1 = subJ0; ilat1 <= subJ1; ilat1++ ) // copy the last row as well..
        {
            //if ( cdoDebugExt>20 ) cdoPrint("cropData(): ilat1=%d; subJ0=%d; subJ1=%d; rowLen=%d ", ilat1, subJ0, subJ1, rowLen );
            memcpy((void*)&array2[array2Idx], (void*)&array1[ilat1*nlon1 + subI0], rowLen*sizeof(double));
            array2Idx += rowLen;
        }

    }
    else
    if ( nwpv == 2 ) // complex numbers ... unsupported yet ...
    {
    }
}


void *SampleGrid(void *argument)
{
    int SAMPLEGRID;
    int SUBGRID;
    int operatorID;
    int streamID1, streamID2;
    int nrecs, nvars;
    int tsID, recID, varID, levelID;
    int gridsize, gridsize2;
    int vlistID1, vlistID2;
    int gridSrcID = -1, gridIDsampled;
    int resampleFactor;
    int subI0,subI1, subJ0, subJ1;
    int index, ngrids, gridtype = -1;
    int nmiss;
    int *vars = NULL;
    int i;
    int nwpv; // number of words per value; real:1  complex:2
    double missval;
    double *array1 = NULL, *array2 = NULL;
    int taxisID1, taxisID2;
    typedef struct {
        int gridSrcID, gridIDsampled;
        int *cellidx, nvals;
        int subI0,subI1, subJ0, subJ1;
    } sbox_t;
    sbox_t *sbox = NULL;


    cdoInitialize(argument);

    SAMPLEGRID  = cdoOperatorAdd("samplegrid",  0, 0, "resample factor, typically 2 (which will half the resolution)");
    SUBGRID  = cdoOperatorAdd("subgrid",  0, 0, " sub-grid indices: i0,i1,j0,j1");

    operatorID = cdoOperatorID();


    int nch = operatorArgc();

    if (operatorID == SAMPLEGRID)
    {
        if ( cdoDebugExt ) cdoPrint("samplegrid operator requested..");
        if ( nch<1 ) cdoAbort("Number of input arguments < 1; At least 1 argument needed: resample-factor (2,3,4, .. etc)");
        if ( ! isdigit(*operatorArgv()[0]) )
            cdoAbort("The input argument is not a number !");
        resampleFactor = atoi(operatorArgv()[0]);

        if ( cdoDebugExt ) cdoPrint("resampleFactor = %d", resampleFactor);
    }
    else
    if (operatorID == SUBGRID)
    {
        if ( cdoDebugExt ) cdoPrint("subgrid operator requested..");
        if ( nch<4 ) cdoAbort("Number of input arguments < 4; Must specify sub-grid indices: i0,i1,j0,j1; This works only with LCC grid. For other grids use: selindexbox");
        if ( ! isdigit(*operatorArgv()[0]) )
            cdoAbort("The input argument is not a number !");
        subI0 = atoi(operatorArgv()[0]);
        subI1 = atoi(operatorArgv()[1]);
        subJ0 = atoi(operatorArgv()[2]);
        subJ1 = atoi(operatorArgv()[3]);

        resampleFactor = atoi(operatorArgv()[0]);

        if ( cdoDebugExt ) cdoPrint("resampleFactor = %d", resampleFactor);
    }
    else
        cdoAbort("Unknown operator ...");

    streamID1 = streamOpenRead(cdoStreamName(0));

    vlistID1 = streamInqVlist(streamID1);
    vlistID2 = vlistDuplicate(vlistID1);

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    nvars = vlistNvars(vlistID1);
    vars  = (int *) malloc(nvars*sizeof(int));
    for ( varID = 0; varID < nvars; varID++ ) vars[varID] = FALSE;

    ngrids = vlistNgrids(vlistID1);

    if ( cdoDebugExt ) cdoPrint("ngrids = %d", ngrids);

    sbox = (sbox_t *) malloc(ngrids*sizeof(sbox_t));

    for ( index = 0; index < ngrids; index++ )
    {
        gridSrcID  = vlistGrid(vlistID1, index);
        gridtype = gridInqType(gridSrcID);
        if ( ( gridtype != GRID_CURVILINEAR ) && (   gridInqXsize(gridSrcID) > 0 && gridInqYsize(gridSrcID) > 0 ) )
        {
            if (operatorID == SAMPLEGRID)
                gridIDsampled = define_sample_grid(gridSrcID, resampleFactor);
            else
                if (operatorID == SUBGRID)
                {
                    if ( gridtype != GRID_LCC )
                        cdoAbort("Unsupported grid type: %s; This works only with LCC grid. For other grids use: selindexbox", gridNamePtr(gridtype));

                    int gridIDcurvl;

                    gridIDcurvl = gridToCurvilinear(gridSrcID, 1);
                    if ( gridInqType(gridIDcurvl) != GRID_CURVILINEAR )
                    {
                        gridDestroy(gridIDcurvl);
                        cdoAbort("cdo SampleGrid: define_subgrid_grid() Creation of curvilinear grid definition failed: type != GRID_CURVILINEAR");
                    }

                    gridIDsampled = define_subgrid_grid(gridSrcID, gridIDcurvl, subI0,subI1, subJ0, subJ1);

                    gridDestroy(gridIDcurvl);
                }
            sbox[index].gridSrcID = gridSrcID;
            sbox[index].gridIDsampled = gridIDsampled;

            if ( cdoDebugExt>=10 ) gridPrint(gridSrcID, 1,0);
            if ( cdoDebugExt>=10 ) gridPrint(gridIDsampled, 1,0);

            vlistChangeGridIndex(vlistID2, index, gridIDsampled);
            for ( varID = 0; varID < nvars; varID++ )
                if ( gridSrcID == vlistInqVarGrid(vlistID1, varID) )
                    vars[varID] = TRUE;
        }
        else
        {
          cdoAbort("Unsupported grid type: %s", gridNamePtr(gridtype));
        }
    }


    if ( cdoDebugExt )
    {
        if (operatorID == SAMPLEGRID)
            cdoPrint("Resampled grid has been created.");
        if (operatorID == SUBGRID)
            cdoPrint("Sub-grid has been created.");
    }

    streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

    streamDefVlist(streamID2, vlistID2);

    gridsize = vlistGridsizeMax(vlistID1);
    if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
    array1 = (double *) malloc(gridsize*sizeof(double));

    gridsize2 = vlistGridsizeMax(vlistID2);
    if ( vlistNumber(vlistID2) != CDI_REAL ) gridsize2 *= 2;
    array2 = (double *) malloc(gridsize2*sizeof(double));

    if ( cdoDebugExt )
    {
        cdoPrint("gridsize = %ld, gridsize2 = %ld, ", gridsize, gridsize2);
    }

    tsID = 0;
    while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
        taxisCopyTimestep(taxisID2, taxisID1);

        streamDefTimestep(streamID2, tsID);

        for ( recID = 0; recID < nrecs; recID++ )
        {
            streamInqRecord(streamID1, &varID, &levelID);
            streamReadRecord(streamID1, array1, &nmiss);

            streamDefRecord(streamID2, varID, levelID);

            if ( cdoDebugExt>=20 ) cdoPrint("Processing record (%d) of %d.",recID, nrecs);

            if ( vars[varID] )
            {
                if (  vlistInqVarDatatype(vlistID1, varID) == DATATYPE_CPX32 ||
                    vlistInqVarDatatype(vlistID1, varID) == DATATYPE_CPX64 )
                    nwpv = 2;
                else
                    nwpv = 1;

                gridSrcID = vlistInqVarGrid(vlistID1, varID);

                for ( index = 0; index < ngrids; index++ )
                    if ( gridSrcID == sbox[index].gridSrcID ) break;

                if ( index == ngrids ) cdoAbort("Internal problem, grid not found!");

                gridsize2 = gridInqSize(sbox[index].gridIDsampled);
                gridIDsampled = sbox[index].gridIDsampled;

                if (operatorID == SAMPLEGRID) {
                    //if ( cdoDebugExt ) cdoPrint("Calling sampleData gridSrcID: %d; gridIDsampled: %d",gridSrcID, gridIDsampled);
                    sampleData(nwpv, array1, gridSrcID, array2, gridIDsampled, resampleFactor);
                }
                else
                    if (operatorID == SUBGRID) {
                        cropData(nwpv, array1, gridSrcID, array2, gridIDsampled, subI0,subI1, subJ0, subJ1);
                    }


                if ( nmiss )
                {
                    nmiss = 0;
                    missval = vlistInqVarMissval(vlistID2, varID);
                    for ( i = 0; i < gridsize2; i++ )
                        if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;
                }
//    if ( cdoDebugExt ) cdoPrint("AAA1");

                streamWriteRecord(streamID2, array2, nmiss);
                //streamWriteRecord(streamID2, array1, nmiss);
            }
            else
            {
                streamWriteRecord(streamID2, array1, nmiss);
            }
        }
        tsID++;
    }

    streamClose(streamID2);
    streamClose(streamID1);

    vlistDestroy(vlistID2);

    if ( vars   ) free(vars);
    if ( array2 ) free(array2);
    if ( array1 ) free(array1);

    if ( sbox )
    {
        free(sbox);
    }

    cdoFinish();

    return (0);
}

#endif
