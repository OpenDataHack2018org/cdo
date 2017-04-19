#include <cdi.h>
#include "cdo_int.h"
#include "grid.h"


typedef struct {
  int ncells;
  int *neighbor; // neighbor cell index
  int *parent;   // parent cell index
  int *child;    // child cell index
  const char *filename;
} cellindex_type;


static
void copy_data_to_index(int ncells, const double *restrict data, int *restrict cellindex)
{
  for ( int i = 0; i < ncells; ++i )
    cellindex[i] = (int) lround(data[i]);
}

static
void free_cellindex(cellindex_type *cellindex)
{
  if ( cellindex->neighbor ) Free(cellindex->neighbor);
  if ( cellindex->parent ) Free(cellindex->parent);
  if ( cellindex->child ) Free(cellindex->child);
  Free(cellindex);
}

static
cellindex_type *read_cellindex(const char *filename)
{
  openLock();
  int streamID = streamOpenRead(filename);
  openUnlock();

  if ( streamID < 0 ) cdiOpenError(streamID, "Open failed on >%s<", filename);

  int vlistID = streamInqVlist(streamID);
  int ngrids = vlistNgrids(vlistID);
  int gridID = -1;
  for ( int index = 0; index < ngrids; ++index )
    {
      gridID = vlistGrid(vlistID, index);
      if ( gridInqType(gridID) == GRID_UNSTRUCTURED && gridInqNvertex(gridID) == 3 ) break;
    }

  if ( gridID == -1 ) cdoAbort("No ICON grid found in %s!", filename);

  int nid = CDI_UNDEFID;
  int pid = CDI_UNDEFID;
  // int cid = CDI_UNDEFID;
  int nvars = vlistNvars(vlistID);
  char varname[CDI_MAX_NAME];
  for ( int varID = 0; varID < nvars; ++varID)
    {
      vlistInqVarName(vlistID, varID, varname);
      if      ( strcmp(varname, "neighbor_cell_index") == 0 ) nid = varID;
      else if ( strcmp(varname, "parent_cell_index") == 0 ) pid = varID;
      // else if ( strcmp(varname, "child_cell_index") == 0 ) cid = varID;
    }

  if ( nid == CDI_UNDEFID ) cdoAbort("neighbor_cell_index not found in %s!", filename);
  if ( pid == CDI_UNDEFID ) cdoAbort("parent_cell_index not found in %s!", filename);
  // if ( cid == CDI_UNDEFID ) cdoAbort("child_cell_index not found in %s!", filename);

  int ncells = gridInqSize(gridID);

  cellindex_type *cellindex = (cellindex_type*) Malloc(sizeof(cellindex_type));
  cellindex->ncells = ncells;

  cellindex->neighbor = (int*) Malloc(3*ncells*sizeof(int));
  cellindex->parent   = (int*) Malloc(  ncells*sizeof(int));
  cellindex->child    = NULL;
  // cellindex->child    = (cid != CDI_UNDEFID) ? (int*) Malloc(4*ncells*sizeof(int)) : NULL;
  double *data = (double *) Malloc(ncells*sizeof(double));

  int nrecs = streamInqTimestep(streamID, 0);
  for ( int recID = 0; recID < nrecs; recID++ )
    {
      int varID, levelID, nmiss;
      streamInqRecord(streamID, &varID, &levelID);
      if ( varID == nid || varID == pid /*|| varID == cid*/ )
        {
          streamReadRecord(streamID, data, &nmiss);
          if      ( varID == nid ) copy_data_to_index(ncells, data, cellindex->neighbor+levelID*ncells);
          else if ( varID == pid ) copy_data_to_index(ncells, data, cellindex->parent);
          //  else if ( varID == cid ) copy_data_to_index(ncells, data, cellindex->child+levelID*ncells);
        }
    }

  // Fortran to C index
  for ( int i = 0; i < 3*ncells; ++i ) cellindex->neighbor[i] -= 1;
  for ( int i = 0; i < ncells; ++i ) cellindex->parent[i] -= 1;

  streamClose(streamID);

  Free(data);

  return cellindex;
}

static
int read_grid(const char *filename)
{
  openLock();
  int streamID = streamOpenRead(filename);
  openUnlock();

  if ( streamID < 0 ) cdiOpenError(streamID, "Open failed on >%s<", filename);

  int vlistID = streamInqVlist(streamID);
  int ngrids = vlistNgrids(vlistID);
  int gridID = -1;
  for ( int index = 0; index < ngrids; ++index )
    {
      gridID = vlistGrid(vlistID, index);
      if ( gridInqType(gridID) == GRID_UNSTRUCTURED && gridInqNvertex(gridID) == 3 ) break;
    }

  if ( gridID == -1 ) cdoAbort("No ICON grid found in %s!", filename);

  int gridID2 = gridDuplicate(gridID);

  streamClose(streamID);

  return gridID2;
}

/**
* Return the first index of element x fits.
*
* If no interval can be found return -1.

* @param *array ascending sorted list
* @param n      length of the sorted list
* @param search the element to find a position for 
*/
static
int find_index(int search, int n, const int *restrict array)
{
  int first = 0;
  int last = n - 1;
  int middle = (first+last)/2;
 
  while ( first <= last )
    {
      if ( array[middle] < search )
        first = middle + 1;    
      else if ( array[middle] == search )
        {
          for ( int i = middle; i >= 0; i-- )
            {
              if ( array[i] == search ) middle = i;
              else break;
            }
          return middle;
        }
      else
        last = middle - 1;
 
      middle = (first + last)/2;
    }

  return -1;
}

static
void compute_child(cellindex_type *cellindex1, cellindex_type *cellindex2)
{
  int ncells1 = cellindex1->ncells;
  int *parent1 = cellindex1->parent;
  int ncells2 = cellindex2->ncells;
  int *child2 = (int*) Malloc(4*ncells2*sizeof(int));
  cellindex2->child = child2;
  for ( int i = 0; i< ncells2; ++i )
    {
      for ( int k = 0; k < 4; ++k ) child2[i*4+k] = -1;
      int j = find_index(i, ncells1, parent1);
      if ( j < 0 ) continue;
      for ( int k = 0; k < 4; ++k )
        {
          if ( i != parent1[j+k] ) break;
          child2[i*4+k] = j+k;
        }
     // if ( i%10000 == 0 ) printf("%d %d %d %d %d %d\n", i, j, parent1[j], parent1[j+1], parent1[j+2], parent1[j+3]);
    }
}

static
void compute_sum(int i, int *n, double *sum, double *sumq, int kci, cellindex_type **cellindex, double *array)
{
  // printf("compute: i, kci %d %d\n", i, kci);
  int ncells2 = cellindex[kci]->ncells;
  if ( i < 0 || i > ncells2 ) cdoAbort("Child grid cell index %d out of bounds %d!", i, ncells2);

  for ( int k = 0; k < 4; ++k )
    {
      int index = cellindex[kci]->child[i*4+k];
      if ( index == -1 ) break;
      if ( kci == 1 )
        {
          *sum += array[index];
          *sumq += array[index]*array[index];
          *n += 1;
        }
      else compute_sum(index, n, sum, sumq, kci-1, cellindex, array);
    }
}

static
void samplegrid(double missval, int nci, cellindex_type **cellindex, double *array1, double *array2, double *array3)
{
  int kci = nci-1;
  int ncells2 = cellindex[kci]->ncells;
#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(missval, ncells2, kci, cellindex, array1, array2, array3)
#endif
  for ( int i = 0; i < ncells2; ++i )
    {
      int n = 0;
      double sum = 0, sumq = 0;
      compute_sum(i, &n, &sum, &sumq, kci, cellindex, array1);
      array2[i] = n ? sum/n : missval;  // mean
      double var1 = (n*n > n) ? (sumq*n - sum*sum) / (n*n - n) : missval;
      if ( var1 < 0 && var1 > -1.e-5 ) var1 = 0;
      array3[i] = var_to_std(var1, missval); // std1
    }
}


#include "pstream.h"


void *Samplegridicon(void *argument)
{
  int nrecs;
  int varID, levelID;

  cdoInitialize(argument);

  cdoOperatorAdd("samplegridicon",  0, 0, "sample grids");

  int nsamplegrids = operatorArgc();
  if ( nsamplegrids < 2 ) cdoAbort("Parameter missing!");

  cellindex_type *cellindex[nsamplegrids];
  
  for ( int i = 0; i < nsamplegrids; ++i )
    {
      cellindex[i] = read_cellindex(operatorArgv()[i]);
      cellindex[i]->filename = operatorArgv()[i];
      if ( cdoVerbose ) cdoPrint("Found %d grid cells in %s", cellindex[i]->ncells, cellindex[i]->filename);
    }

  for ( int i = 0; i < nsamplegrids-1; ++i )
    compute_child(cellindex[i], cellindex[i+1]);
  
  int gridID2 = read_grid(operatorArgv()[nsamplegrids-1]);

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);
  int vlistID3 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  int taxisID3 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);
  vlistDefTaxis(vlistID3, taxisID3);

  int ngrids = vlistNgrids(vlistID1);
  for ( int index = 0; index < ngrids; ++index )
    {
      int gridID = vlistGrid(vlistID1, index);
      int gridtype = gridInqType(gridID);
      if ( !(gridtype == GRID_UNSTRUCTURED && gridInqNvertex(gridID) == 3) )
        cdoAbort("Unsupported gridtype: %s with %d corners", gridNamePtr(gridtype), gridInqNvertex(gridID));

      vlistChangeGridIndex(vlistID2, index, gridID2);
      vlistChangeGridIndex(vlistID3, index, gridID2);
    }

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  streamDefVlist(streamID2, vlistID2);

  int streamID3 = streamOpenWrite(cdoStreamName(2), cdoFiletype());
  streamDefVlist(streamID3, vlistID3);

  int gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
  double *array1 = (double *) Malloc(gridsize*sizeof(double));

  int gridsize2 = gridInqSize(gridID2);
  if ( cdoVerbose ) cdoPrint("Target gridsize = %d", gridsize2);
  if ( vlistNumber(vlistID2) != CDI_REAL ) gridsize2 *= 2;
  double *array2 = (double *) Malloc(gridsize2*sizeof(double));
  double *array3 = (double *) Malloc(gridsize2*sizeof(double));

  int tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);
      taxisCopyTimestep(taxisID3, taxisID1);
      streamDefTimestep(streamID2, tsID);
      streamDefTimestep(streamID3, tsID);

      for ( int recID = 0; recID < nrecs; recID++ )
        {
          int nmiss;
          streamInqRecord(streamID1, &varID, &levelID);
          streamReadRecord(streamID1, array1, &nmiss);

          double missval = vlistInqVarMissval(vlistID1, varID);

          samplegrid(missval, nsamplegrids, cellindex, array1, array2, array3);

          nmiss = 0;
          for ( int i = 0; i < gridsize2; ++i )
            if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;

          streamDefRecord(streamID2, varID, levelID);
          streamWriteRecord(streamID2, array2, nmiss);

          nmiss = 0;
          for ( int i = 0; i < gridsize2; ++i )
            if ( DBL_IS_EQUAL(array3[i], missval) ) nmiss++;

          streamDefRecord(streamID3, varID, levelID);
          streamWriteRecord(streamID3, array3, nmiss);
        }

      tsID++;
    }

  streamClose(streamID3);
  streamClose(streamID2);
  streamClose(streamID1);

  vlistDestroy(vlistID2);
  gridDestroy(gridID2);

  if ( array3 ) Free(array3);
  if ( array2 ) Free(array2);
  if ( array1 ) Free(array1);

  for ( int i = 0; i < nsamplegrids; ++i ) free_cellindex(cellindex[i]);

  cdoFinish();

  return 0;
}
