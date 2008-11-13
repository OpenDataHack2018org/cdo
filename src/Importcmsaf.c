#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif

#define H5_USE_16_API

#if  defined  (HAVE_LIBHDF5)
#  include "hdf5.h"
#endif

#include <ctype.h>

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define  MAX_DSETS  1024

typedef struct {
  char *name;
  char *description;
  char *units;
  int dtype;
  int nx;
  int ny;
  int nz;
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
void print_filter(hid_t dset_id, char *varname)
{
  hid_t plist;
  H5Z_filter_t filter;
  unsigned int flags;
  int idx;
  unsigned int cd_values;
  int nfilter;
  size_t cd_nelmts = 1;
  size_t pnamelen = 64;
  char pname[64];

  /* get filter */
  plist = H5Dget_create_plist(dset_id);
  nfilter = H5Pget_nfilters(plist);

  for ( idx = 0; idx < nfilter; idx++ )
    {
      filter = H5Pget_filter(plist, idx, &flags, &cd_nelmts, &cd_values, 
			     pnamelen, pname);
      cdoPrint("Dataset %s: filter %d =  %s", varname, idx+1, pname);
    }

  H5Pclose(plist);
}
#endif


#if  defined  (HAVE_LIBHDF5)
static
void read_dataset(hid_t loc_id, const char *name, void *opdata)
{
  static char func[] = "read_dataset";
  hid_t dset_id, type_id;
  hid_t   dataspace;   
  hsize_t dims_out[9];  /* dataset dimensions           */
  herr_t  status;	/* Generic return value		*/
  hid_t   attr, atype, atype_mem;
  hid_t   native_type;
  int iattr;
  char attname[256];
  H5T_class_t  type_class;
  H5T_class_t atype_class;
  int     rank;
  int nx, ny, nz;
  int gridsize, offset;
  double *array;
  double addoffset = 0, scalefactor = 1, missval = cdiInqMissval();
  int lmissval = 0;
  int nset;
  int i;
  int ftype = 0;
  int len;
  int laddoffset, lscalefactor;
  int dtype = DATATYPE_FLT32;
  char attstring[4096];     /* Buffer to read string attribute back */
  char varname[256];
  short *mask = NULL;
  double minval, maxval;
  int nmiss;
  int num_attrs;

  attstring[0] = 0;
  strcpy(varname, name);

  dset_id = H5Dopen(loc_id, varname);
  
  type_id = H5Dget_type(dset_id);  /* get datatype*/

  type_class = H5Tget_class(type_id);
  if ( type_class < 0 )
    {
      cdoAbort(" Invalid datatype for %s", varname);
    }
  /*
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
  */
  native_type = H5Tget_native_type(type_id, H5T_DIR_ASCEND);
  if      ( H5Tequal(native_type, H5T_NATIVE_CHAR)   > 0 ) {ftype=0; dtype = DATATYPE_INT8;}
  else if ( H5Tequal(native_type, H5T_NATIVE_UCHAR)  > 0 ) {ftype=0; dtype = DATATYPE_UINT8;}
  else if ( H5Tequal(native_type, H5T_NATIVE_SHORT)  > 0 ) {ftype=0; dtype = DATATYPE_INT16;}
  else if ( H5Tequal(native_type, H5T_NATIVE_USHORT) > 0 ) {ftype=0; dtype = DATATYPE_UINT16;}
  else if ( H5Tequal(native_type, H5T_NATIVE_INT)    > 0 ) {ftype=0; dtype = DATATYPE_INT32;}
  else if ( H5Tequal(native_type, H5T_NATIVE_UINT)   > 0 ) {ftype=0; dtype = DATATYPE_UINT32;}
  else if ( H5Tequal(native_type, H5T_NATIVE_FLOAT)  > 0 ) {ftype=1; dtype = DATATYPE_FLT32;}
  else if ( H5Tequal(native_type, H5T_NATIVE_DOUBLE) > 0 ) {ftype=1; dtype = DATATYPE_FLT64;}
  else
    {
      cdoWarning("Dataset %s skipped, unsupported native datatype!", varname);
      goto RETURN;
    }
  H5Tclose(native_type);

  dataspace = H5Dget_space(dset_id);    /* dataspace handle */
  rank      = H5Sget_simple_extent_ndims(dataspace);
  status    = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
  if ( rank != 2 )
    {
      cdoWarning("Dataset %s skipped, unsupported rank (=%d)!", rank);
      goto RETURN;
    }

  nx = dims_out[1];
  ny = dims_out[0];
  nz = 1;

  len = (int) strlen(varname);
  if ( isdigit(varname[len-1]) )
    {
      nz = atoi(&varname[len-1]);
      varname[len-1] = 0;
    }

  gridsize = nx*ny;

  if ( nz == 1 )
    nset = ((DSETS *) opdata)->nsets;
  else
    {
      for ( nset = 0; nset < ((DSETS *) opdata)->nsets; ++nset )
	{
	  if ( strcmp(varname, ((DSETS *) opdata)->obj[nset].name) == 0 ) break;
	}

      if ( nset >= ((DSETS *) opdata)->nsets )
	cdoAbort("3D var %s not found!", varname);
    }

  if ( nset < MAX_DSETS )
    {
      if ( cdoVerbose ) print_filter(dset_id, varname);

      num_attrs = H5Aget_num_attrs(dset_id);
      for( i = 0; i < num_attrs; i++ )
	{
	  attr = H5Aopen_idx(dset_id, i);
	  atype = H5Aget_type(attr);
	  H5Aget_name(attr, sizeof(attname), attname);

	  if ( strcmp(attname, "CLASS") == 0 ||
	       strcmp(attname, "IMAGE_VERSION") == 0 ||
	       strcmp(attname, "PALETTE") == 0 ) continue;

	  atype_mem = H5Tget_native_type(atype, H5T_DIR_ASCEND);
	  atype_class = H5Tget_class(atype);

	  if ( strcmp(attname, "intercept") == 0 )
	    {
	      if ( atype_class == H5T_FLOAT )
		{
		  H5Aread(attr, H5T_NATIVE_DOUBLE, &addoffset);
		  laddoffset = 1;
		}
	      else if ( atype_class == H5T_INTEGER )
		{
		  H5Aread(attr, H5T_NATIVE_INT, &iattr);
		  addoffset = iattr;
		  laddoffset = 1;
		}
	      else
		cdoWarning("Attribute %s has unsupported data type!", attname);
	    }
	  else if ( strcmp(attname, "gain") == 0 )
	    {
	      if ( atype_class == H5T_FLOAT )
		{
		  H5Aread(attr, H5T_NATIVE_DOUBLE, &scalefactor);
		  lscalefactor = 1;
		}
	      else if ( atype_class == H5T_INTEGER )
		{
		  H5Aread(attr, H5T_NATIVE_INT, &iattr);
		  scalefactor = iattr;
		  lscalefactor = 1;
		}
	      else
		cdoWarning("Attribute %s has unsupported data type!", attname);
	    }
	  else if ( strcmp(attname, "no_data_value") == 0 ||
		    strcmp(attname, "nodata value") == 0  ||
		    strcmp(attname, "nodata") == 0 )
	    {
	      if ( atype_class == H5T_FLOAT )
		{
		  H5Aread(attr, H5T_NATIVE_DOUBLE, &missval);
		  lmissval = 1;
		}
	      else if ( atype_class == H5T_INTEGER )
		{
		  H5Aread(attr, H5T_NATIVE_INT, &iattr);
		  missval = iattr;
		  lmissval = 1;
		}
	      else
		cdoWarning("Attribute %s has unsupported data type!", attname);
	    }
	  else if ( strcmp(attname, "description") == 0 )
	    {
	      H5Aread(attr, atype_mem, attstring);
	      if ( ((DSETS *) opdata)->obj[nset].description )
		free(((DSETS *) opdata)->obj[nset].description);
	      ((DSETS *) opdata)->obj[nset].description = strdup(attstring);
	    }
	  /*
	  else if ( strcmp(attname, "title") == 0 )
	    {
	      H5Aread(attr, atype_mem, attstring);
	      if ( ((DSETS *) opdata)->obj[nset].description == NULL )
		((DSETS *) opdata)->obj[nset].description = strdup(attstring);
	    }
	  */
	  else if ( strcmp(attname, "unit") == 0 )
	    {
	      H5Aread(attr, atype_mem, attstring);
	      ((DSETS *) opdata)->obj[nset].units = strdup(attstring);
	    }

	  H5Tclose(atype_mem);
	  H5Aclose(attr);
	  H5Tclose(atype);
	}
      
      offset = gridsize*(nz-1);
      array  = ((DSETS *) opdata)->obj[nset].array;
      array  = (double *) realloc(array, gridsize*nz*sizeof(double));
      ((DSETS *) opdata)->obj[nset].array    = array;
      array  = array+offset;

      mask = (short *) malloc(gridsize*sizeof(short));
      memset(mask, 0, gridsize*sizeof(short));

      if ( ftype )
	{
	  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, array);
	}
      else
	{
	  int *iarray, i;
	  iarray = (int *) malloc(gridsize*sizeof(int));
	  status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, iarray);
	  for ( i = 0; i < gridsize; ++i ) array[i] = iarray[i];
	  free(iarray);
	}

      ((DSETS *) opdata)->obj[nset].name     = strdup(varname);
      ((DSETS *) opdata)->obj[nset].nx       = nx;
      ((DSETS *) opdata)->obj[nset].ny       = ny;
      ((DSETS *) opdata)->obj[nset].nz       = nz;
      ((DSETS *) opdata)->obj[nset].gridsize = gridsize;

      ((DSETS *) opdata)->obj[nset].dtype    = dtype;

      ((DSETS *) opdata)->obj[nset].loffset  = laddoffset;
      ((DSETS *) opdata)->obj[nset].lscale   = lscalefactor;
      ((DSETS *) opdata)->obj[nset].lmissval = lmissval;
      ((DSETS *) opdata)->obj[nset].offset   = addoffset;
      ((DSETS *) opdata)->obj[nset].scale    = scalefactor;
      ((DSETS *) opdata)->obj[nset].missval  = missval;

      if ( nz == 1 ) ((DSETS *) opdata)->nsets++;

      nmiss  = 0;
 
      minval =  1e35;
      maxval = -1e35;

      for ( i = 0; i < gridsize; i++ )
	{
	  if ( array[i] < minval ) minval = array[i];
	  if ( array[i] > maxval ) maxval = array[i];
	}

      if ( cdoVerbose )
	cdoPrint("Dataset %s: dtype = %d  minval = %g  maxval = %g",
		 varname, dtype,  minval, maxval);

      if ( dtype == DATATYPE_UINT8 )
	{
	  if ( minval >= 0 && maxval <= 127 )
	    {
	      dtype = DATATYPE_INT8;
	      ((DSETS *) opdata)->obj[nset].dtype    = dtype;	      
	    }
	}
      else if ( dtype == DATATYPE_UINT16 )
	{
	  if ( minval >= 0 && maxval <= 32767 )
	    {
	      dtype = DATATYPE_INT16;
	      ((DSETS *) opdata)->obj[nset].dtype    = dtype;	      
	    }
	}

      laddoffset   = !DBL_IS_EQUAL(addoffset, 0);
      lscalefactor = !DBL_IS_EQUAL(scalefactor,  1);

      if ( laddoffset || lscalefactor )
	{	  
	  for ( i = 0; i < gridsize; i++ )
	    if ( !DBL_IS_EQUAL(array[i], missval) )
	      {
		mask[i] = 0;

		if ( lscalefactor ) array[i] *= scalefactor;
		if ( laddoffset   ) array[i] += addoffset;
	      }
	    else
	      {
		nmiss++;
		mask[i] = 1;
	      }
	}
	    
      minval =  1e35;
      maxval = -1e35;

      for ( i = 0; i < gridsize; i++ )
	if ( mask[i] == 0 )
	  {
	    if ( array[i] < minval ) minval = array[i];
	    if ( array[i] > maxval ) maxval = array[i];
	  }

      if ( cdoVerbose )
	cdoPrint("Dataset %s: dtype = %d  minval = %g  maxval = %g",
		 varname, dtype,  minval, maxval);

      if ( nmiss > 0 )
	{
	  if ( ! (missval < minval || missval > maxval) )
	    {
	      if ( DBL_IS_EQUAL(missval, 255) && dtype == DATATYPE_UINT8 )
		{
		  missval = -255;
		  dtype   = DATATYPE_INT16;
		  ((DSETS *) opdata)->obj[nset].dtype    = dtype;
		  ((DSETS *) opdata)->obj[nset].missval  = missval;
		  cdoPrint("Dataset %s: change missval to %g and datatype to INT16!",
			   varname, missval);

		  for ( i = 0; i < gridsize; i++ )
		    if ( mask[i] ) array[i] = missval;
		}
	      else
		cdoWarning(" Missval is inside of valid values!\n"
                           " Name: %s  Range: %g - %g  Missval: %g\n",
			   varname, minval, maxval, missval);
	    }
	}

      free(mask); 
      mask = NULL;
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
obj_info(hid_t loc_id, const char *name, void *opdata)
{
  H5G_obj_t obj_type;
  H5G_stat_t statbuf;

  /* avoid compiler warnings */
  loc_id = loc_id;
  opdata = opdata;

  H5Gget_objinfo(loc_id, name, FALSE, &statbuf);

  obj_type = statbuf.type;

  switch (obj_type) {
  case H5G_GROUP:
    if ( cdoVerbose ) cdoPrint(" Object with name %s is a group", name);
    if ( strcmp(name, "Data") == 0 )
      {
	H5Giterate(loc_id, name, NULL, obj_info, opdata);
      }
    break;
  case H5G_DATASET: 
    if ( cdoVerbose ) cdoPrint(" Object with name %s is a dataset", name);
    if ( strstr(name, "PALETTE") )
      {
	if ( cdoVerbose ) cdoPrint("   Skip dataset: %s", name);
      }
    else if ( strstr(name, "egion") )
      {
	if ( cdoVerbose ) cdoPrint("   Skip dataset: %s", name);
      }
    else
      {
	if ( cdoVerbose ) cdoPrint("   Read dataset: %s", name);
	read_dataset(loc_id, name, opdata);
      }
    break;
  case H5G_TYPE: 
    if ( cdoVerbose ) cdoPrint(" Object with name %s is a named datatype", name);
    break;
  default:
    cdoAbort(" Unable to identify an object %s", name);
  }

  return 0;
}
#endif


#if  defined  (HAVE_LIBHDF5)
void get_global_att(hid_t file_id, int vlistID)
{
  hid_t   attr, atype, atype_mem, obj_id, grp_id = -1;
  char attname[256];
  H5T_class_t  type_class;
  int attint;
  double attflt;
  int i, pos;
  int num_attrs;
  char attstring[4096];     /* Buffer to read string attribute back */

  attstring[0] = 0;

  obj_id = H5Gopen(file_id, "/");

  num_attrs = H5Aget_num_attrs(obj_id);
  if ( num_attrs == 0 )
    {
      grp_id = H5Gopen(obj_id, "Metadata");

      if ( grp_id < 0 ) return;

      obj_id = grp_id;

      num_attrs = H5Aget_num_attrs(obj_id);
    }

  for( i = 0; i < num_attrs; i++ )
    {
      attr = H5Aopen_idx(obj_id, i);
      atype = H5Aget_type(attr);
      H5Aget_name(attr, sizeof(attname), attname);

      /* remove illegal characters */
      for ( pos = 0; pos < (int)strlen(attname); ++pos )
	if ( attname[pos] == '&' ) attname[pos] = '_';

      atype_mem = H5Tget_native_type(atype, H5T_DIR_ASCEND);
      type_class = H5Tget_class(atype);
      if ( type_class == H5T_STRING )
	{
	  H5Aread(attr, atype_mem, attstring);
	  vlistDefAttTxt(vlistID, CDI_GLOBAL, attname, (int)strlen(attstring), attstring);
	} 
      else if ( type_class == H5T_INTEGER )
	{
	  H5Aread(attr, H5T_NATIVE_INT, &attint);
	  vlistDefAttInt(vlistID, CDI_GLOBAL, attname, 1, &attint);
	} 
      else if ( type_class == H5T_FLOAT )
	{
	  H5Aread(attr, H5T_NATIVE_DOUBLE, &attflt);
	  vlistDefAttFlt(vlistID, CDI_GLOBAL, attname, 1, &attflt);
	} 
      H5Tclose(atype_mem);
      H5Aclose(attr);
      H5Tclose(atype);
    }

  if ( grp_id >= 0 ) H5Gclose(grp_id);	

}
#endif


void dsets_init(DSETS *dsets)
{
  int i;

  dsets->nsets = 0;

  for ( i = 0; i < MAX_DSETS; ++i )
    {
      dsets->obj[i].nx          = 0;
      dsets->obj[i].ny          = 0;
      dsets->obj[i].nz          = 0;
      dsets->obj[i].name        = NULL;
      dsets->obj[i].description = NULL;
      dsets->obj[i].units       = NULL;
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
  int i, offset;
  int nmiss;
  int ivar;
  int varID, levelID;
  int nx, ny, nz, gridsize;
  double *array;
  double missval, minval, maxval;
  hid_t	  file_id;	/* HDF5 File ID	        	*/
  herr_t  status;	/* Generic return value		*/
  DSETS dsets;

  cdoInitialize(argument);

  dsets_init(&dsets);

  /* Open an existing file. */
  file_id = H5Fopen(cdoStreamName(0), H5F_ACC_RDONLY, H5P_DEFAULT);

  /* cmsaf_type = get_cmsaf_type(file_id); */

  H5Giterate(file_id, "/", NULL, obj_info, (void *) &dsets);

  if ( dsets.nsets == 0 ) cdoAbort("No dataset found!");

  gridsize = dsets.obj[0].gridsize;
  nx = dsets.obj[0].nx;
  ny = dsets.obj[0].ny;
  nz = dsets.obj[0].nz;

  if ( cdoVerbose )
    for ( ivar = 0; ivar < dsets.nsets; ++ivar )
      cdoPrint(" Var %d %-20s %dx%d nlev = %d", ivar, dsets.obj[ivar].name, nx, ny, nz);

  for ( ivar = 1; ivar < dsets.nsets; ++ivar )
    {
      if ( nx != dsets.obj[0].nx || ny != dsets.obj[0].ny )
	cdoAbort("Gridsize must not change!");
      if ( nz != dsets.obj[0].nz )
	cdoAbort("Number of levels must not change!");
    }

  gridID = gridCreate(GRID_GENERIC, gridsize);
  gridDefXsize(gridID, nx);
  gridDefYsize(gridID, ny);

  if ( nz == 1 )
    zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);
  else
    {
      double *levels;
      levels = (double *) malloc(nz*sizeof(double));
      for ( i = 0; i < nz; ++i ) levels[i] = i+1;
      zaxisID = zaxisCreate(ZAXIS_GENERIC, nz);
      zaxisDefLevels(zaxisID, levels);
      free(levels);
    }

  vlistID = vlistCreate();

  taxisID = taxisCreate(TAXIS_ABSOLUTE);
  vlistDefTaxis(vlistID, taxisID);

  for ( ivar = 0; ivar < dsets.nsets; ++ivar )
    {
      varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_CONSTANT);
      vlistDefVarName(vlistID, varID,  dsets.obj[ivar].name);
      if ( dsets.obj[ivar].description )
	vlistDefVarLongname(vlistID, varID,  dsets.obj[ivar].description);
      if ( dsets.obj[ivar].units )
	vlistDefVarUnits(vlistID, varID,  dsets.obj[ivar].units);
	
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

      gridsize = dsets.obj[ivar].gridsize;
      missval  = dsets.obj[ivar].missval;

      for ( levelID = 0; levelID < nz; ++levelID )
	{
	  offset = gridsize*levelID;
	  array  = dsets.obj[ivar].array+offset;

	  nmiss  = 0;
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

	  if ( cdoVerbose )
	    cdoPrint(" Write var %d,  level %d, nmiss %d, missval %g, minval %g, maxval %g",
		 varID, levelID, nmiss, missval, minval, maxval);
	  /*
	  if ( ! (missval < minval || missval > maxval) )
	    cdoWarning(" Missval is inside of valid values! Name: %s  Range: %g - %g  Missval: %g\n",
		       dsets.obj[ivar].name, minval, maxval, missval);
	  */
	  streamDefRecord(streamID,  varID,  levelID);
	  streamWriteRecord(streamID, array, nmiss);
	}
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
