#include <cdi.h>
#include "cdo_int.h"
#include "grid.h"

typedef struct {
  int ncells;
  int *neighbor; // neighbor cell index
  int *parent;   // parent cell index
  int *child;    // child cell index
} cellindex_type;

static
void copy_data_to_index(int ncells, const double *restrict data, int *restrict cellindex)
{
  for ( int i = 0; i < ncells; ++i )
    cellindex[i] = (int) lround(data[i]);
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
* Find the interval i-1 .. i in which an element x fits and return i, the 
* bigger one of the interval borders or x itself if it is an interval border.
*
* If no interval can be found return the length of the array.

* @param *array ascending or descending sorted list
* @param nelem  length of the sorted list
* @param x      the element to find a position for 
*/
static
int find_element(int x, int nelem, const int *restrict array)
{
  int ii;
  int mid = 0;
  int first = 1;
  int last = nelem;

  if ( array[0] < array[nelem-1] ) // ascending order
    {
      /* return the length of the array if x is out of bounds */
      if ( x < array[0] || x > array[nelem-1] ) return nelem;

      /* search for the interval in which x fits */
      // implementation: binary search algorithm
      for ( ii = 1; ii < nelem; ++ii )
	{
	  // binary search: divide search room in the middle
	  // mid = first + ((last - first) >> 1);
	  // faster!
	  mid = (first + last) >> 1;
      
	  /* return the bigger interval border of the interval in which x fits */
	  // if ( x >= array[mid-1] && x <= array[mid] ) break;
	  // faster!
	  if ( !(x < array[mid-1] || x > array[mid]) ) break;

	  // binary search: ignore half of the search room
	  if ( x > array[mid] )
	    first = mid;
	  else
	    last = mid;
	}
    }
  else
    {
      /* return the length of the array if x is out of bounds */
      if ( x < array[nelem-1] || x > array[0] ) return nelem;

      /* search for the interval in which x fits */
      // implementation: binary search algorithm
      for ( ii = 1; ii < nelem; ++ii )
	{
	  // binary search: divide search room in the middle
	  // mid = first + ((last - first) >> 1);
	  // faster!
	  mid = (first + last) >> 1;
      
	  /* return the bigger interval border of the interval in which x fits */
	  // if ( x >= array[mid] && x <= array[mid-1] ) break;
	  // faster!
	  if ( !(x < array[mid] || x > array[mid-1]) ) break;

	  // binary search: ignore half of the search room
	  if ( x < array[mid] )
	    first = mid;
	  else
	    last = mid;
	}
    }

  if ( mid > 1 && IS_EQUAL(x,array[mid-1]) ) mid--;

  return mid;
}

static
void samplegrid(cellindex_type *cellindex1, double *array1, cellindex_type *cellindex2, double *array2, int *samp2)
{
  int ncells1 = cellindex1->ncells;
  int ncells2 = cellindex2->ncells;
  for ( int i = 0; i< ncells2; ++i )
    {
      int j = find_element(i+1, ncells1, cellindex1->parent);
      if ( i%10000 == 0 )
          printf("%d %d %d %d %d %d %d\n", i, j, cellindex1->parent[j-2], cellindex1->parent[j-1], cellindex1->parent[j],
                 cellindex1->parent[j+1], cellindex1->parent[j+2]);
      if ( j == i+1 )
        {
        }
    }
}

static
void compute_child(cellindex_type *cellindex1, cellindex_type *cellindex2)
{
  int ncells2 = cellindex2->ncells;
  cellindex2->child = (int*) Malloc(4*ncells2*sizeof(int));
}

/*
static
void map_index(int num_links, int dst_size, int *src_idx, int *dst_idx, double *src_array, double *dst_array, double missval)
{
  for ( int n = 0; n < dst_size; ++n ) dst_array[n] = missval;
  for ( int n = 0; n < num_links; ++n ) dst_array[dst_idx[n]] = 0.;
  for ( int n = 0; n < num_links; ++n ) dst_array[dst_idx[n]] += src_array[src_idx[n]];
}
*/

#include "pstream.h"


void *Samplegridicon(void *argument)
{
  int nrecs;
  int varID, levelID;
  int nmiss;

  cdoInitialize(argument);

  cdoOperatorAdd("samplegridicon",  0, 0, "sample grids");

  int nsamplegrids = operatorArgc();
  if ( nsamplegrids < 2 ) cdoAbort("Parameter missing!");

  cellindex_type *cellindex1 = read_cellindex(operatorArgv()[0]);
  printf("ncells %d\n", cellindex1->ncells);
  cellindex_type *cellindex2 = read_cellindex(operatorArgv()[1]);
  printf("ncells %d\n", cellindex2->ncells);

  compute_child(cellindex1, cellindex2);
  
  int gridID2 = read_grid(operatorArgv()[1]);

  int streamID1 = streamOpenRead(cdoStreamName(0));

  int vlistID1 = streamInqVlist(streamID1);
  int vlistID2 = vlistDuplicate(vlistID1);

  int taxisID1 = vlistInqTaxis(vlistID1);
  int taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  int ngrids = vlistNgrids(vlistID1);
  for ( int index = 0; index < ngrids; ++index )
    {
      int gridID = vlistGrid(vlistID1, index);
      int gridtype = gridInqType(gridID);
      if ( !(gridtype == GRID_UNSTRUCTURED && gridInqNvertex(gridID) == 3) )
        cdoAbort("Unsupported gridtype: %s with %d corners", gridNamePtr(gridtype), gridInqNvertex(gridID));

      vlistChangeGridIndex(vlistID2, index, gridID2);
    }

  int streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  int gridsize = vlistGridsizeMax(vlistID1);
  if ( vlistNumber(vlistID1) != CDI_REAL ) gridsize *= 2;
  double *array1 = (double *) Malloc(gridsize*sizeof(double));

  int gridsize2 = gridInqSize(gridID2);
  if ( vlistNumber(vlistID2) != CDI_REAL ) gridsize2 *= 2;
  double *array2 = (double *) Malloc(gridsize2*sizeof(double));
  int *samp2 = (int *) Malloc(gridsize2*sizeof(int));
  for ( int i = 0; i < gridsize2; ++i ) samp2[i] = 0;
  for ( int i = 0; i < gridsize2; ++i ) array2[i] = 0;

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

          int gridID = vlistInqVarGrid(vlistID1, varID);

          //for ( int i = 0; i < gridsize2; ++i ) array2[i] = array1[i];
          samplegrid(cellindex1, array1, cellindex2, array2, samp2);

          if ( nmiss )
            {
              nmiss = 0;
              double missval = vlistInqVarMissval(vlistID2, varID);
              for ( int i = 0; i < gridsize2; i++ )
                if ( DBL_IS_EQUAL(array2[i], missval) ) nmiss++;
            }

          streamWriteRecord(streamID2, array2, nmiss);
        }

      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  vlistDestroy(vlistID2);
  gridDestroy(gridID2);

  if ( samp2 ) Free(samp2);
  if ( array2 ) Free(array2);
  if ( array1 ) Free(array1);

  if ( cellindex1->neighbor ) Free(cellindex1->neighbor);
  if ( cellindex1->parent ) Free(cellindex1->parent);
  if ( cellindex1->child ) Free(cellindex1->child);

  if ( cellindex2->neighbor ) Free(cellindex2->neighbor);
  if ( cellindex2->parent ) Free(cellindex2->parent);
  if ( cellindex2->child ) Free(cellindex2->child);

  cdoFinish();

  return 0;
}
