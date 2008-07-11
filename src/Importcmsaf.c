#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif

#if  defined  (HAVE_LIBHDF5)
#  include "hdf5.h"
#endif

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

#define  MAX_DSETS  1024

typedef struct {
  char *name;
  int dtype;
  int nx;
  int ny;
  int gridsize;
  double scale;
  double offset;
  double missval;
  double *array;
}
DSET_OBJ;

typedef struct {
  int nsets;
  DSET_OBJ obj[MAX_DSETS];
}
DSETS;


#if  defined  (HAVE_LIBHDF5)
static
void read_dataset(hid_t loc_id, const char *name, void *opdata)
{
  static char func[] = "read_dataset";
  hid_t dset_id, type_id;
  H5T_class_t t_class;
  H5O_info_t oinfo;           /* Object info */
  hid_t   dataspace;   
  hsize_t dims_out[9];  /* dataset dimensions           */
  herr_t  status;	/* Generic return value		*/
  hid_t   attr, atype, atype_mem;
  H5T_class_t  type_class;
  int     rank;
  int gridsize;
  double *array = NULL;
  int nsets;
  int i;
  char    string_out[80];     /* Buffer to read string attribute back */

  dset_id = H5Dopen2(loc_id, name, H5P_DEFAULT);
  printf("dsetid: %d %s\n", dset_id, name);
  
  type_id = H5Dget_type(dset_id);  /* get datatype*/

  t_class = H5Tget_class(type_id);
  if(t_class < 0) {
    puts(" Invalid datatype.\n");
  }
  else {
    if(t_class == H5T_INTEGER)
      puts(" Datatype is 'H5T_NATIVE_INTEGER'.\n");
    if(t_class == H5T_FLOAT)
      puts(" Datatype is 'H5T_NATIVE_FLOAT'.\n");
    if(t_class == H5T_STRING)
      puts(" Datatype is 'H5T_NATIVE_STRING'.\n");
    if(t_class == H5T_BITFIELD)
      puts(" Datatype is 'H5T_NATIVE_BITFIELD'.\n");
    if(t_class == H5T_OPAQUE)
      puts(" Datatype is 'H5T_NATIVE_OPAQUE'.\n");
    if(t_class == H5T_COMPOUND)
      puts(" Datatype is 'H5T_NATIVE_COMPOUND'.\n");
  }

  dataspace = H5Dget_space(dset_id);    /* dataspace handle */
  rank      = H5Sget_simple_extent_ndims(dataspace);
  status    = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
  if ( rank != 2 )
    {
      cdoWarning("Unexpected rank = %d!", rank);
      goto RETURN;
    }
  printf("\nRank: %d\nDimensions: %lu x %lu \n", rank,
	 (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]));

  gridsize = dims_out[0]*dims_out[1];

  array = (double *) malloc(gridsize*sizeof(double));

  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, array);

  /* free(array); */

  nsets = ((DSETS *) opdata)->nsets;
  if ( nsets < MAX_DSETS )
    {
      ((DSETS *) opdata)->obj[nsets].name     = strdup(name);
      ((DSETS *) opdata)->obj[nsets].nx       = dims_out[0];
      ((DSETS *) opdata)->obj[nsets].ny       = dims_out[1];
      ((DSETS *) opdata)->obj[nsets].gridsize = gridsize;
      ((DSETS *) opdata)->obj[nsets].array    = array;

      H5Oget_info(dset_id, &oinfo);
      for(i = 0; i < (unsigned)oinfo.num_attrs; i++)
	{
	  attr = H5Aopen_by_idx(dset_id, ".", H5_INDEX_CRT_ORDER, H5_ITER_INC, (hsize_t)i, H5P_DEFAULT, H5P_DEFAULT);
	  atype = H5Aget_type(attr);
	  type_class = H5Tget_class(atype);
	  if (type_class == H5T_STRING)
	    {
	      atype_mem = H5Tget_native_type(atype, H5T_DIR_ASCEND);
	      H5Aread(attr, atype_mem, string_out);
	      printf("Found string attribute; its index is %d , value =   %s \n", i, string_out);
	      H5Tclose(atype_mem);
	    } 
	  H5Aclose(attr);
	  H5Tclose(atype);
	}
      
      ((DSETS *) opdata)->nsets++;
    }
  else
    {
      cdoWarning("Too many datasets (MAX = %d)!", MAX_DSETS);
      goto RETURN;
    }

 RETURN:

  H5Sclose(dataspace);
  H5Dclose(dset_id);
  H5Tclose(type_id);
}
#endif

#if  defined  (HAVE_LIBHDF5)
static herr_t
file_info(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata)
{
  H5O_info_t obj_info;

  /* avoid compiler warnings */
  loc_id = loc_id;
  opdata = opdata;
  linfo = linfo;

  /*
   * Display group name. The name is passed to the function by
   * the Library. Some magic :-)
   */
  H5Oget_info_by_name(loc_id, name, &obj_info, H5P_DEFAULT);

  switch (obj_info.type) {
  case H5G_GROUP: 
    printf(" Object with name %s is a group \n", name);
    break;
  case H5G_DATASET: 
    printf(" Object with name %s is a dataset \n", name);
    if ( strstr(name, "PALETTE") )
      {
	printf(" Skip dataset: %s\n", name);
      }
    else if ( strstr(name, "egion") )
      {
	printf(" Skip dataset: %s\n", name);
      }
    else
      {
	printf(" Read dataset: %s\n", name);
	read_dataset(loc_id, name, opdata);
      }
    break;
  case H5G_TYPE: 
    printf(" Object with name %s is a named datatype \n", name);
    break;
  default:
    printf(" Unable to identify an object ");
  }

  return 0;
}
#endif

void dsets_init(DSETS *dsets)
{
  int i;

  dsets->nsets = 0;

  for ( i = 0; i < MAX_DSETS; ++i )
    {
      dsets->obj[i].name    = NULL;
      dsets->obj[i].dtype   = cdoDefaultDataType;
      dsets->obj[i].scale   = 1;
      dsets->obj[i].offset  = 0;
      dsets->obj[i].missval = cdiInqMissval();
      dsets->obj[i].array   = NULL;   
    }
}


void *Importcmsaf(void *argument)
{
#if  defined  (HAVE_LIBHDF5)
  static char func[] = "Importcmsaf";
  int streamID;
  int gridID, zaxisID, taxisID, vlistID;
  int i;
  int nvars;
  int nmiss;
  int ivar;
  int varID, levelID;
  int nx, ny, gridsize;
  double *array;
  hid_t	  file_id;	/* HDF5 File ID	        	*/
  herr_t  status;	/* Generic return value		*/
  DSETS dsets;

  cdoInitialize(argument);

  dsets_init(&dsets);

  /* Open an existing file. */
  file_id = H5Fopen(cdoStreamName(0), H5F_ACC_RDONLY, H5P_DEFAULT);

  H5Literate(file_id, H5_INDEX_NAME, H5_ITER_INC, NULL, file_info, (void *) &dsets);

  if ( dsets.nsets == 0 ) cdoAbort("No dataset found!");

  if ( cdoVerbose )
    for ( ivar = 0; ivar < dsets.nsets; ++ivar )
      {
	cdoPrint("%d %-20s %dx%d", ivar+1, dsets.obj[ivar].name, dsets.obj[ivar].nx, dsets.obj[ivar].ny);
      }

  nx = dsets.obj[0].nx;
  ny = dsets.obj[0].ny;
  gridsize = dsets.obj[0].gridsize;

  for ( ivar = 1; ivar < dsets.nsets; ++ivar )
    {
      if ( nx != dsets.obj[0].nx || ny != dsets.obj[0].ny )
	cdoAbort("Gridsize must not change!");
    }

  gridID = gridCreate(GRID_GENERIC, gridsize);
  gridDefXsize(gridID, nx);
  gridDefYsize(gridID, ny);

  zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);

  vlistID = vlistCreate();

  taxisID = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID, taxisID);

  for ( ivar = 0; ivar < dsets.nsets; ++ivar )
    {
      varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_CONSTANT);
      vlistDefVarName(vlistID, varID,  dsets.obj[ivar].name);
      /*
      vlistDefVarUnits(vlistID, varID, units[i]);
      */
      vlistDefVarDatatype(vlistID, varID, dsets.obj[ivar].dtype);
      vlistDefVarMissval(vlistID, varID, dsets.obj[ivar].missval);
      vlistDefVarScalefactor(vlistID, varID, dsets.obj[ivar].scale);
      vlistDefVarAddoffset(vlistID, varID, dsets.obj[ivar].offset);
    }

  streamID = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID, vlistID);

  streamDefTimestep(streamID, 0);

  for ( ivar = 0; ivar < dsets.nsets; ++ivar )
    {
      varID = ivar;
      levelID = 0;
      nmiss = 0;
      array = dsets.obj[ivar].array;

      streamDefRecord(streamID,  varID,  levelID);
      streamWriteRecord(streamID, array, nmiss);
    }

  /* Close file */
  status = H5Fclose(file_id);

  processDefVarNum(vlistNvars(vlistID), streamID);

  streamClose(streamID);

  vlistDestroy(vlistID);
  gridDestroy(gridID);
  zaxisDestroy(zaxisID);
  taxisDestroy(taxisID);

  for ( ivar = 0; ivar < dsets.nsets; ++ivar )
    free(dsets.obj[ivar].array);

  cdoFinish();
#else
  cdoWarning("HDF5 support not compiled in!");
#endif

  return (0);
}
