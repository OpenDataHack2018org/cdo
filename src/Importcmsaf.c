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
  char *description;
  int dtype;
  int nx;
  int ny;
  int gridsize;
  int lscale;
  int loffset;
  int lmissval;
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
  H5O_info_t oinfo;           /* Object info */
  hid_t   dataspace;   
  hsize_t dims_out[9];  /* dataset dimensions           */
  herr_t  status;	/* Generic return value		*/
  hid_t   attr, atype, atype_mem;
  hid_t   native_type;
  char attname[256];
  H5T_class_t  type_class;
  int     rank;
  int gridsize;
  int attint;
  double attflt;
  double *array = NULL;
  double offset = 0, scale = 1, missval = cdiInqMissval();
  int loffset = 0, lscale = 0, lmissval = 0;
  int nsets;
  int i;
  int laddoffset, lscalefactor;
  int dtype = DATATYPE_FLT32;
  char    attstring[4096];     /* Buffer to read string attribute back */

  dset_id = H5Dopen2(loc_id, name, H5P_DEFAULT);
  
  type_id = H5Dget_type(dset_id);  /* get datatype*/

  type_class = H5Tget_class(type_id);
  if(type_class < 0) {
    puts(" Invalid datatype.\n");
  }
  else {
    if(type_class == H5T_INTEGER)
      puts("   Datatype is 'H5T_NATIVE_INTEGER'.\n");
    if(type_class == H5T_FLOAT)
      puts("   Datatype is 'H5T_NATIVE_FLOAT'.\n");
    if(type_class == H5T_STRING)
      puts("   Datatype is 'H5T_NATIVE_STRING'.\n");
    if(type_class == H5T_BITFIELD)
      puts("   Datatype is 'H5T_NATIVE_BITFIELD'.\n");
    if(type_class == H5T_OPAQUE)
      puts("   Datatype is 'H5T_NATIVE_OPAQUE'.\n");
    if(type_class == H5T_COMPOUND)
      puts("   Datatype is 'H5T_NATIVE_COMPOUND'.\n");
  }

  native_type = H5Tget_native_type(type_id, H5T_DIR_ASCEND);
  if      ( H5Tequal(native_type, H5T_NATIVE_CHAR)   > 0 ) dtype = DATATYPE_INT8;
  else if ( H5Tequal(native_type, H5T_NATIVE_UCHAR)  > 0 ) dtype = DATATYPE_UINT8;
  else if ( H5Tequal(native_type, H5T_NATIVE_SHORT)  > 0 ) dtype = DATATYPE_INT16;
  else if ( H5Tequal(native_type, H5T_NATIVE_USHORT) > 0 ) dtype = DATATYPE_UINT16;
  else if ( H5Tequal(native_type, H5T_NATIVE_INT)    > 0 ) dtype = DATATYPE_INT32;
  else if ( H5Tequal(native_type, H5T_NATIVE_UINT)   > 0 ) dtype = DATATYPE_UINT32;
  else if ( H5Tequal(native_type, H5T_NATIVE_FLOAT)  > 0 ) dtype = DATATYPE_FLT32;
  else if ( H5Tequal(native_type, H5T_NATIVE_DOUBLE) > 0 ) dtype = DATATYPE_FLT64;
  else
    {
      cdoWarning("%s skipped, unsupported native datatype!", name);
      goto RETURN;
    }
  H5Tclose(native_type);

  dataspace = H5Dget_space(dset_id);    /* dataspace handle */
  rank      = H5Sget_simple_extent_ndims(dataspace);
  status    = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
  if ( rank != 2 )
    {
      cdoWarning("%s skipped, unsupported rank (=%d)!", rank);
      goto RETURN;
    }

  gridsize = dims_out[0]*dims_out[1];

  array = (double *) malloc(gridsize*sizeof(double));

  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, array);

  /* free(array); */

  nsets = ((DSETS *) opdata)->nsets;
  if ( nsets < MAX_DSETS )
    {
      H5Oget_info(dset_id, &oinfo);
      for( i = 0; i < (int)oinfo.num_attrs; i++ )
	{
	  attr = H5Aopen_by_idx(dset_id, ".", H5_INDEX_CRT_ORDER, H5_ITER_INC, 
				(hsize_t)i, H5P_DEFAULT, H5P_DEFAULT);
	  atype = H5Aget_type(attr);
	  H5Aget_name(attr, sizeof(attname), attname);

	  if ( strcmp(attname, "CLASS") == 0 ||
	       strcmp(attname, "IMAGE_VERSION") == 0 ||
	       strcmp(attname, "PALETTE") == 0 ) continue;

	  type_class = H5Tget_class(atype);

	  atype_mem = H5Tget_native_type(atype, H5T_DIR_ASCEND);

	  if ( strcmp(attname, "intercept") == 0 )
	    {
	      H5Aread(attr, H5T_NATIVE_DOUBLE, &offset);
	      loffset = 1;
	    }
	  else if ( strcmp(attname, "gain") == 0 )
	    {
	      H5Aread(attr, H5T_NATIVE_DOUBLE, &scale);
	      lscale = 1;
	    }
	  else if ( strcmp(attname, "no_data_value") == 0 )
	    {
	      H5Aread(attr, H5T_NATIVE_DOUBLE, &missval);
	      lmissval = 1;
	    }
	  else if ( strcmp(attname, "description") == 0 )
	    {
	      H5Aread(attr, atype_mem, attstring);
	    }

	  H5Tclose(atype_mem);
	  H5Aclose(attr);
	  H5Tclose(atype);
	}
      
      ((DSETS *) opdata)->obj[nsets].name     = strdup(name);
      ((DSETS *) opdata)->obj[nsets].nx       = dims_out[0];
      ((DSETS *) opdata)->obj[nsets].ny       = dims_out[1];
      ((DSETS *) opdata)->obj[nsets].gridsize = gridsize;
      ((DSETS *) opdata)->obj[nsets].array    = array;

      ((DSETS *) opdata)->obj[nsets].dtype    = dtype;

      ((DSETS *) opdata)->obj[nsets].loffset  = loffset;
      ((DSETS *) opdata)->obj[nsets].lscale   = lscale;
      ((DSETS *) opdata)->obj[nsets].lmissval = lmissval;
      ((DSETS *) opdata)->obj[nsets].offset   = offset;
      ((DSETS *) opdata)->obj[nsets].scale    = scale;
      ((DSETS *) opdata)->obj[nsets].missval  = missval;

      if ( attstring[0] )
	((DSETS *) opdata)->obj[nsets].description = strdup(attstring);

      ((DSETS *) opdata)->nsets++;

      laddoffset   = !DBL_IS_EQUAL(offset, 0);
      lscalefactor = !DBL_IS_EQUAL(scale,  1);

      if ( laddoffset || lscalefactor )
	{	  
	  for ( i = 0; i < gridsize; i++ )
	    if ( !DBL_IS_EQUAL(array[i], missval) )
	      {
		if ( lscalefactor ) array[i] *= scale;
		if ( laddoffset )   array[i] += offset;
	      }
	}
    }
  else
    {
      cdoWarning("Too many datasets (MAX = %d)!", MAX_DSETS);
      goto RETURN;
    }

  H5Sclose(dataspace);

 RETURN:

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
	printf("   Skip dataset: %s\n", name);
      }
    else if ( strstr(name, "egion") )
      {
	printf("   Skip dataset: %s\n", name);
      }
    else
      {
	printf("   Read dataset: %s\n", name);
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


#if  defined  (HAVE_LIBHDF5)
void get_global_att(hid_t file_id, int vlistID)
{
  static char func[] = "get_global_att";
  H5O_info_t oinfo;           /* Object info */
  herr_t  status;	/* Generic return value		*/
  hid_t   attr, atype, atype_mem;
  char attname[256];
  H5T_class_t  type_class;
  int attint;
  double attflt;
  int i;
  char attstring[4096];     /* Buffer to read string attribute back */

  attstring[0] = 0;

  H5Oget_info(file_id, &oinfo);
  for( i = 0; i < (int)oinfo.num_attrs; i++ )
    {
      attr = H5Aopen_by_idx(file_id, ".", H5_INDEX_CRT_ORDER, H5_ITER_INC, 
			    (hsize_t)i, H5P_DEFAULT, H5P_DEFAULT);
      atype = H5Aget_type(attr);
      H5Aget_name(attr, sizeof(attname), attname);

      type_class = H5Tget_class(atype);
      if ( type_class == H5T_STRING )
	{
	  atype_mem = H5Tget_native_type(atype, H5T_DIR_ASCEND);
	  H5Aread(attr, atype_mem, attstring);
	  H5Tclose(atype_mem);
	  vlistDefAttTxt(vlistID, CDI_GLOBAL, attname, (int)strlen(attstring), attstring);
	} 
      else if ( type_class == H5T_INTEGER )
	{
	  atype_mem = H5Tget_native_type(atype, H5T_DIR_ASCEND);
	  H5Aread(attr, H5T_NATIVE_INT, &attint);
	  H5Tclose(atype_mem);
	  vlistDefAttInt(vlistID, CDI_GLOBAL, attname, 1, &attint);
	} 
      else if ( type_class == H5T_FLOAT )
	{
	  atype_mem = H5Tget_native_type(atype, H5T_DIR_ASCEND);
	  H5Aread(attr, H5T_NATIVE_DOUBLE, &attflt);
	  H5Tclose(atype_mem);
	  vlistDefAttFlt(vlistID, CDI_GLOBAL, attname, 1, &attflt);
	} 
      H5Aclose(attr);
      H5Tclose(atype);
    }
}
#endif


void dsets_init(DSETS *dsets)
{
  int i;

  dsets->nsets = 0;

  for ( i = 0; i < MAX_DSETS; ++i )
    {
      dsets->obj[i].name        = NULL;
      dsets->obj[i].description = NULL;
      dsets->obj[i].dtype       = cdoDefaultDataType;
      dsets->obj[i].lscale      = 0;
      dsets->obj[i].loffset     = 0;
      dsets->obj[i].lmissval    = 0;
      dsets->obj[i].missval     = cdiInqMissval();
      dsets->obj[i].array       = NULL;   
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
  double missval, minval, maxval;
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
      if ( dsets.obj[ivar].description )
	vlistDefVarLongname(vlistID, varID,  dsets.obj[ivar].description);
	
      /*
      vlistDefVarUnits(vlistID, varID, units[i]);
      */
      vlistDefVarDatatype(vlistID, varID, dsets.obj[ivar].dtype);
      if ( dsets.obj[ivar].lmissval )
	vlistDefVarMissval(vlistID, varID, dsets.obj[ivar].missval);
      if ( dsets.obj[ivar].lscale )
	vlistDefVarScalefactor(vlistID, varID, dsets.obj[ivar].scale);
      if ( dsets.obj[ivar].loffset )
	vlistDefVarAddoffset(vlistID, varID, dsets.obj[ivar].offset);
    }

  get_global_att(file_id, vlistID);

  streamID = streamOpenWrite(cdoStreamName(1), cdoFiletype());
  if ( streamID < 0 ) cdiError(streamID, "Open failed on %s", cdoStreamName(1));

  streamDefVlist(streamID, vlistID);

  streamDefTimestep(streamID, 0);

  for ( ivar = 0; ivar < dsets.nsets; ++ivar )
    {
      varID   = ivar;
      levelID = 0;

      gridsize = dsets.obj[ivar].gridsize;
      missval  = dsets.obj[ivar].missval;

      nmiss = 0;
      array = dsets.obj[ivar].array;

      minval =  1e35;
      maxval = -1e35;

      for ( i = 0; i < gridsize; i++ )
	{
	  if ( !DBL_IS_EQUAL(array[i], missval) )
	    {
	      if ( array[i] < minval ) minval = array[i];
	      if ( array[i] > maxval ) maxval = array[i];
	    }
	}

      for ( i = 0; i < gridsize; i++ )
	if ( DBL_IS_EQUAL(array[i], missval) ) nmiss++;

      printf(" Write var %d,  nmiss %d, missval %g, minval %g, maxval %g\n",
	     varID, nmiss, missval, minval, maxval);

      if ( ! (missval < minval || missval > maxval) )
	printf(" Warning: missval is inside of valid values\n");

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
