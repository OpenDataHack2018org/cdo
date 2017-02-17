/*
   This module "SampleGrid" contains the following operators:

    samplegrid      Resample current grid with given factor, typically 2 (which will half the resolution);
                    tested on curvilinear and LCC grids;
    subgrid         Similar to selindexbox but this operator works for LCC grids (tested on HARMONIE NWP model).
*/

#include "cdi.h"
#include "cdo_int.h"
#include "grid.h"
#include "griddes.h"
#include "pstream.h"


extern int cdoDebugExt; // defined in cdo.c

static
void sampleData(double *array1, int gridID1, double *array2, int gridID2, int resampleFactor)
{
  long nlon1 = gridInqXsize(gridID1);
  long nlat1 = gridInqYsize(gridID1);

  long nlon2 = gridInqXsize(gridID2);
  long nlat2 = gridInqYsize(gridID2);

  if ( cdoDebugExt >= 100 )
    cdoPrint("sampleData():: (nlon1: %d; nlat1: %d) => (nlon2: %d; nlat2: %d); gridID1: %d; gridID2: %d; resampleFactor: %d)",
             nlon1,nlat1, nlon2,nlat2, gridID1, gridID2, resampleFactor);

  for ( long ilat1 = 0; ilat1 < nlat1; ilat1+=resampleFactor )
    for ( long ilon1 = 0; ilon1 < nlon1; ilon1+=resampleFactor )
      *array2++ = array1[ilat1*nlon1 + ilon1];
}

static
void cropData(double *array1, int gridID1, double *array2, int gridID2, int subI0, int subI1, int  subJ0, int  subJ1 )
{
  long nlon1 = gridInqXsize(gridID1);
  long nlon2 = gridInqXsize(gridID2);
  long rowLen = subI1 - subI0 + 1; // must be same as   nlon1

  if ( rowLen!= nlon2 )
    cdoAbort("cropData() rowLen!= nlon2 [%d != %d]", rowLen, nlon2);

  if ( cdoDebugExt>=10 ) cdoPrint("cropData(%d,%d,%d,%d) ...\n",subI0,subI1, subJ0, subJ1 );

  long array2Idx = 0;
  for ( long ilat1 = subJ0; ilat1 <= subJ1; ilat1++ ) // copy the last row as well..
    {
      //if ( cdoDebugExt>20 ) cdoPrint("cropData(): ilat1=%d; subJ0=%d; subJ1=%d; rowLen=%d ", ilat1, subJ0, subJ1, rowLen );
      memcpy((void*)&array2[array2Idx], (void*)&array1[ilat1*nlon1 + subI0], rowLen*sizeof(double));
      array2Idx += rowLen;
    }
}


void *SampleGrid(void *argument)
{
  int nrecs;
  int varID, levelID;
  int resampleFactor;
  int subI0 = 0, subI1 = 0, subJ0 = 0, subJ1 = 0;
  int index;
  int nmiss;
  typedef struct {
    int gridSrcID, gridIDsampled;
    int *cellidx, nvals;
    int subI0,subI1, subJ0, subJ1;
  } sbox_t;

  cdoInitialize(argument);

  int SAMPLEGRID  = cdoOperatorAdd("samplegrid",  0, 0, "resample factor, typically 2 (which will half the resolution)");
  int SUBGRID  = cdoOperatorAdd("subgrid",  0, 0, " sub-grid indices: i0,i1,j0,j1");

  int operatorID = cdoOperatorID();

  int nch = operatorArgc();

  if ( operatorID == SAMPLEGRID )
    {
      if ( cdoDebugExt ) cdoPrint("samplegrid operator requested..");
      if ( nch<1 ) cdoAbort("Number of input arguments < 1; At least 1 argument needed: resample-factor (2,3,4, .. etc)");
      if ( ! isdigit(*operatorArgv()[0]) )
        cdoAbort("The input argument is not a number !");
      resampleFactor = atoi(operatorArgv()[0]);

      if ( cdoDebugExt ) cdoPrint("resampleFactor = %d", resampleFactor);
    }
  else if ( operatorID == SUBGRID )
    {
      if ( cdoDebugExt ) cdoPrint("subgrid operator requested..");
      if ( nch<4 ) cdoAbort("Number of input arguments < 4; Must specify sub-grid indices: i0,i1,j0,j1; This works only with LCC grid. For other grids use: selindexbox");
      if ( ! isdigit(*operatorArgv()[0]) )
        cdoAbort("The input argument is not a number !");
      subI0 = atoi(operatorArgv()[0]);
      subI1 = atoi(operatorArgv()[1]);
      subJ0 = atoi(operatorArgv()[2]);
      subJ1 = atoi(operatorArgv()[3]);
    }
  else
    cdoAbort("Unknown operator ...");

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int nvars = vlistNvars(vlistID1);
  bool *vars  = (bool *) Malloc(nvars*sizeof(bool));
  for ( varID = 0; varID < nvars; varID++ ) vars[varID] = false;

  int ngrids = vlistNgrids(vlistID1);

  if ( cdoDebugExt ) cdoPrint("ngrids = %d", ngrids);

  sbox_t *sbox = (sbox_t *) Malloc(ngrids*sizeof(sbox_t));

  for ( int index = 0; index < ngrids; index++ )
    {
      int gridSrcID = vlistGrid(vlistID1, index);
      int gridIDsampled = -1;
      int gridtype = gridInqType(gridSrcID);
      if ( (gridtype != GRID_CURVILINEAR) && (gridInqXsize(gridSrcID) > 0 && gridInqYsize(gridSrcID) > 0 ) )
        {
          if ( operatorID == SAMPLEGRID )
            gridIDsampled = cdo_define_sample_grid(gridSrcID, resampleFactor);
          else if ( operatorID == SUBGRID )
            {
              if ( gridtype != GRID_LCC )
                cdoAbort("Unsupported grid type: %s; This works only with LCC grid. For other grids use: selindexbox", gridNamePtr(gridtype));

              int gridIDcurvl = gridToCurvilinear(gridSrcID, 1);
              if ( gridInqType(gridIDcurvl) != GRID_CURVILINEAR )
                {
                  gridDestroy(gridIDcurvl);
                  cdoAbort("cdo SampleGrid: define_subgrid_grid() Creation of curvilinear grid definition failed: type != GRID_CURVILINEAR");
                }

              // TODO gridIDsampled = define_subgrid_grid(gridSrcID, gridIDcurvl, subI0,subI1, subJ0, subJ1);
              cdoAbort("Call to define_subgrid_grid() missing!");
                    
              gridDestroy(gridIDcurvl);
            }

          sbox[index].gridSrcID = gridSrcID;
          sbox[index].gridIDsampled = gridIDsampled;

          // TODO if ( cdoDebugExt>=10 ) gridPrint(gridSrcID, 1,0);
          // if ( cdoDebugExt>=10 ) gridPrint(gridIDsampled, 1,0);

          vlistChangeGridIndex(vlistID2, index, gridIDsampled);
          for ( varID = 0; varID < nvars; varID++ )
            if ( gridSrcID == vlistInqVarGrid(vlistID1, varID) )
              vars[varID] = true;
        }
      else
        {
          cdoAbort("Unsupported grid type: %s", gridNamePtr(gridtype));
        }
    }

  if ( cdoDebugExt )
    {
      if ( operatorID == SAMPLEGRID ) cdoPrint("Resampled grid has been created.");
      if ( operatorID == SUBGRID    ) cdoPrint("Sub-grid has been created.");
    }

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
  double *array1 = (double *) Malloc(gridsize*sizeof(double));

  int gridsize2 = vlistGridsizeMax(vlistID2);
  if ( vlistNumber(vlistID2) != CDI_REAL ) gridsize2 *= 2;
  double *array2 = (double *) Malloc(gridsize2*sizeof(double));

  if ( cdoDebugExt ) cdoPrint("gridsize = %ld, gridsize2 = %ld", gridsize, gridsize2);

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
        {
          streamInqRecord(streamID1, &varID, &levelID);
          streamReadRecord(streamID1, array1, &nmiss);

          streamDefRecord(streamID2, varID, levelID);

          if ( cdoDebugExt>=20 ) cdoPrint("Processing record (%d) of %d.",recID, nrecs);

          if ( vars[varID] )
            {
              int gridSrcID = vlistInqVarGrid(vlistID1, varID);

              for ( index = 0; index < ngrids; index++ )
                if ( gridSrcID == sbox[index].gridSrcID ) break;

              if ( index == ngrids ) cdoAbort("Internal problem, grid not found!");

              int gridIDsampled = sbox[index].gridIDsampled;
              gridsize2 = gridInqSize(gridIDsampled);

              if (operatorID == SAMPLEGRID)
                {
                  //if ( cdoDebugExt ) cdoPrint("Calling sampleData gridSrcID: %d; gridIDsampled: %d",gridSrcID, gridIDsampled);
                  sampleData(array1, gridSrcID, array2, gridIDsampled, resampleFactor);
                }
              else if (operatorID == SUBGRID)
                {
                  cropData(array1, gridSrcID, array2, gridIDsampled, subI0,subI1, subJ0, subJ1);
                }

              if ( nmiss )
                {
                  nmiss = 0;
                  double missval = vlistInqVarMissval(vlistID2, varID);
                  for ( int i = 0; i < gridsize2; i++ )
                    if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;
                }

              streamWriteRecord(streamID2, array2, nmiss);
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

  if ( vars ) Free(vars);
  if ( array2 ) Free(array2);
  if ( array1 ) Free(array1);

  if ( sbox ) Free(sbox);

  cdoFinish();

  return 0;
}
