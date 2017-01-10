#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <cdi.h>
#include "cdo_int.h"

static void
printDblsPrefixAutoBrk(FILE *fp, int dig, const char prefix[], size_t nbyte0,
                       size_t n, const double vals[])
{
  fputs(prefix, fp);
  size_t nbyte = nbyte0;
  for ( size_t i = 0; i < n; i++ )
    {
      if ( nbyte > 80 )
        {
          fprintf(fp, "\n%*s", (int)nbyte0, "");
          nbyte = nbyte0;
        }
      nbyte += (size_t)fprintf(fp, "%.*g ", dig, vals[i]);
    }
  fputs("\n", fp);
}

static void
printIntsPrefixAutoBrk(FILE *fp, const char prefix[], size_t nbyte0,
                       size_t n, const int vals[])
{
  fputs(prefix, fp);
  size_t nbyte = nbyte0;
  for ( size_t i = 0; i < n; i++ )
    {
      if ( nbyte > 80 )
        {
          fprintf(fp, "\n%*s", (int)nbyte0, "");
          nbyte = nbyte0;
        }
      nbyte += (size_t)fprintf(fp, "%d ", vals[i]);
    }
  fputs("\n", fp);
}

static void
printBounds(FILE *fp, int dig, const char prefix[], size_t nbyte0,
            size_t n, size_t nvertex, const double bounds[])
{
  fputs(prefix, fp);
  if ( n > 0 )
    {
      for ( size_t iv = 0; iv < nvertex; iv++ )
        fprintf(fp, "%.*g ", dig, bounds[iv]);
      for ( size_t i = 1; i < (size_t)n; i++ )
        {
          fprintf(fp, "\n%*s", (int)nbyte0, "");
          for ( size_t iv = 0; iv < nvertex; iv++ )
            fprintf(fp, "%.*g ", dig, bounds[i*nvertex+iv]);
        }
      fputs("\n", fp);
    }
}

static void
printMask(FILE *fp, const char prefix[], size_t nbyte0,
          size_t n, const int mask[])
{
  fputs(prefix, fp);
  size_t nbyte = nbyte0;
  for ( size_t i = 0; i < n; i++ )
    {
      if ( nbyte > 80 )
        {
          fprintf(fp, "\n%*s", (int)nbyte0, "");
          nbyte = nbyte0;
        }
      nbyte += (size_t)fprintf(fp, "%d ", mask[i]);
    }
  fputs("\n", fp);
}

static inline
void *resizeBuffer(void **buf, size_t *bufSize, size_t reqSize)
{
  if (reqSize > *bufSize)
    {
      *buf = Realloc(*buf, reqSize);
      *bufSize = reqSize;
    }
  return *buf;
}

static
void gridPrintAttributes(FILE *fp, int gridID)
{
  int cdiID = gridID;
  int varID = CDI_GLOBAL;
  int atttype, attlen;
  char attname[CDI_MAX_NAME+1];
  void *attBuf = NULL;
  size_t attBufSize = 0;

  int natts;
  cdiInqNatts(cdiID, varID, &natts);

  for ( int iatt = 0; iatt < natts; ++iatt )
    {
      cdiInqAtt(cdiID, varID, iatt, attname, &atttype, &attlen);

      if ( attlen == 0 ) continue;

      if ( atttype == CDI_DATATYPE_TXT )
        {
          size_t attSize = (size_t)(attlen+1)*sizeof(char);
          char *atttxt = (char *)resizeBuffer(&attBuf, &attBufSize, attSize);
          cdiInqAttTxt(cdiID, varID, attname, attlen, atttxt);
          atttxt[attlen] = 0;
          fprintf(fp, "ATTR_TXT: %s = \"%s\"\n", attname, atttxt);
        }
      else if ( atttype == CDI_DATATYPE_INT8  || atttype == CDI_DATATYPE_UINT8  ||
                atttype == CDI_DATATYPE_INT16 || atttype == CDI_DATATYPE_UINT16 ||
                atttype == CDI_DATATYPE_INT32 || atttype == CDI_DATATYPE_UINT32 )
        {
          size_t attSize = (size_t)attlen*sizeof(int);
          int *attint = (int *)resizeBuffer(&attBuf, &attBufSize, attSize);
          cdiInqAttInt(cdiID, varID, attname, attlen, &attint[0]);
          if ( attlen == 1 )
            fprintf(fp, "ATTR_INT: %s =", attname);
          else
            fprintf(fp, "ATTR_INT_%d: %s =", attlen, attname);
          for ( int i = 0; i < attlen; ++i ) fprintf(fp, " %d", attint[i]);
          fprintf(fp, "\n");
        }
      else if ( atttype == CDI_DATATYPE_FLT32 || atttype == CDI_DATATYPE_FLT64 )
        {
          size_t attSize = (size_t)attlen * sizeof(double);
          double *attflt = (double *)resizeBuffer(&attBuf, &attBufSize, attSize);
          int dig = (atttype == CDI_DATATYPE_FLT64) ? 15 : 7;
          cdiInqAttFlt(cdiID, varID, attname, attlen, attflt);
          if ( attlen == 1 )
            fprintf(fp, "ATTR_FLT: %s =", attname);
          else
            fprintf(fp, "ATTR_FLT_%d: %s =", attlen, attname);
          for ( int i = 0; i < attlen; ++i ) fprintf(fp, " %.*g", dig, attflt[i]);
          fprintf(fp, "\n");
        }
    }

  Free(attBuf);
}

static
void gridPrintKernel(int gridID, int opt, FILE *fp)
{
  int xdim, ydim;
  char attstr[CDI_MAX_NAME];
  char attstr2[CDI_MAX_NAME];
  size_t nxvals = (size_t) gridInqXvals(gridID, NULL);
  size_t nyvals = (size_t) gridInqYvals(gridID, NULL);
  size_t nxbounds = (size_t) gridInqXbounds(gridID, NULL);
  size_t nybounds = (size_t) gridInqYbounds(gridID, NULL);

  int type     = gridInqType(gridID);
  int gridsize = gridInqSize(gridID);
  int xsize    = gridInqXsize(gridID);
  int ysize    = gridInqYsize(gridID);
  int nvertex  = gridInqNvertex(gridID);
  int prec     = gridInqPrec(gridID);

  int dig = (prec == CDI_DATATYPE_FLT64) ? 15 : 7;

  fprintf(fp, "gridtype  = %s\n" "gridsize  = %d\n", gridNamePtr(type), gridsize);

  if ( type != GRID_GME )
    {
      if ( type != GRID_UNSTRUCTURED && type != GRID_SPECTRAL && type != GRID_FOURIER )
        {
          if ( xsize > 0 ) fprintf(fp, "xsize     = %d\n", xsize);
          if ( ysize > 0 ) fprintf(fp, "ysize     = %d\n", ysize);
        }

      if ( nxvals > 0 )
        {
          attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_XNAME, CDI_MAX_NAME, attstr);
          if ( attstr[0] )  fprintf(fp, "xname     = %s\n", attstr);
          attstr2[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_XDIMNAME, CDI_MAX_NAME, attstr2);
          if ( attstr2[0] && strcmp(attstr, attstr2) )  fprintf(fp, "xdimname  = %s\n", attstr2);
          attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_XLONGNAME, CDI_MAX_NAME, attstr);
          if ( attstr[0] )  fprintf(fp, "xlongname = %s\n", attstr);
          attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_XUNITS, CDI_MAX_NAME, attstr);
          if ( attstr[0] )  fprintf(fp, "xunits    = %s\n", attstr);
        }

      if ( nyvals > 0 )
        {
          attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_YNAME, CDI_MAX_NAME, attstr);
          if ( attstr[0] )  fprintf(fp, "yname     = %s\n", attstr);
          attstr2[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_YDIMNAME, CDI_MAX_NAME, attstr2);
          if ( attstr2[0] && strcmp(attstr, attstr2) )  fprintf(fp, "ydimname  = %s\n", attstr2);
          attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_YLONGNAME, CDI_MAX_NAME, attstr);
          if ( attstr[0] )  fprintf(fp, "ylongname = %s\n", attstr);
          attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_YUNITS, CDI_MAX_NAME, attstr);
          if ( attstr[0] )  fprintf(fp, "yunits    = %s\n", attstr);
        }

      if ( type == GRID_UNSTRUCTURED || type == GRID_CURVILINEAR )
        {
          attstr[0] = 0; cdiGridInqKeyStr(gridID, CDI_KEY_VDIMNAME, CDI_MAX_NAME, attstr);
          if ( attstr[0] ) fprintf(fp, "vdimname  = %s\n", attstr);
        }
      if ( type == GRID_UNSTRUCTURED && nvertex > 0 ) fprintf(fp, "nvertex   = %d\n", nvertex);
    }

  switch (type)
    {
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_GAUSSIAN_REDUCED:
    case GRID_GENERIC:
    case GRID_PROJECTION:
    case GRID_CURVILINEAR:
    case GRID_UNSTRUCTURED:
      {
        if ( type == GRID_GAUSSIAN || type == GRID_GAUSSIAN_REDUCED ) fprintf(fp, "np        = %d\n", gridInqNP(gridID));

	if ( type == GRID_CURVILINEAR || type == GRID_UNSTRUCTURED )
	  {
	    xdim = gridsize;
	    ydim = gridsize;
	  }
        else if ( type == GRID_GAUSSIAN_REDUCED )
          {
	    xdim = 2;
	    ydim = ysize;
          }
	else
	  {
	    xdim = xsize;
	    ydim = ysize;
	  }

	if ( type == GRID_UNSTRUCTURED )
          {
            int number = gridInqNumber(gridID);
            int position = gridInqPosition(gridID);
            // const unsigned char *d;
            if ( number > 0 )
              {
                fprintf(fp, "number    = %d\n", number);
                if ( position >= 0 ) fprintf(fp, "position  = %d\n", position);
              }

            if ( gridInqReference(gridID, NULL) )
              {
                char reference_link[8192];
                gridInqReference(gridID, reference_link);
                fprintf(fp, "uri       = %s\n", reference_link);
              }
          }

	if ( nxvals > 0 )
	  {
	    double xfirst = 0.0, xinc = 0.0;

	    if ( type == GRID_LONLAT     || type == GRID_GAUSSIAN ||
		 type == GRID_PROJECTION || type == GRID_GENERIC )
	      {
		xfirst = gridInqXval(gridID, 0);
		xinc   = gridInqXinc(gridID);
	      }

	    if ( IS_NOT_EQUAL(xinc, 0) && opt )
	      {
                fprintf(fp, "xfirst    = %.*g\n"
                        "xinc      = %.*g\n", dig, xfirst, dig, xinc);
	      }
	    else
	      {
                double *xvals = (double*) Malloc(nxvals*sizeof(double));
                gridInqXvals(gridID, xvals);
                static const char prefix[] = "xvals     = ";
                printDblsPrefixAutoBrk(fp, dig, prefix, sizeof(prefix)-1, nxvals, xvals);
                Free(xvals);
	      }
	  }

	if ( nxbounds )
	  {
            double *xbounds = (double*) Malloc(nxbounds*sizeof(double));
            gridInqXbounds(gridID, xbounds);
            static const char prefix[] = "xbounds   = ";
            printBounds(fp, dig, prefix, sizeof(prefix)-1, xdim, nvertex, xbounds);
            Free(xbounds);
	  }

	if ( nyvals > 0 )
	  {
	    double yfirst = 0.0, yinc = 0.0;

	    if ( type == GRID_LONLAT || type == GRID_GENERIC ||
                 type == GRID_PROJECTION || type == GRID_GENERIC )
	      {
		yfirst = gridInqYval(gridID, 0);
		yinc   = gridInqYinc(gridID);
	      }

	    if ( IS_NOT_EQUAL(yinc, 0) && opt )
	      {
	  	fprintf(fp, "yfirst    = %.*g\n"
                        "yinc      = %.*g\n", dig, yfirst, dig, yinc);
	      }
	    else
	      {
                double *yvals = (double*) Malloc(nyvals*sizeof(double));
                gridInqYvals(gridID, yvals);
                static const char prefix[] = "yvals     = ";
                printDblsPrefixAutoBrk(fp, dig, prefix, sizeof(prefix)-1, nyvals, yvals);
                Free(yvals);
	      }
	  }

	if ( nybounds )
	  {
            double *ybounds = (double*) Malloc(nybounds*sizeof(double));
            gridInqYbounds(gridID, ybounds);
            static const char prefix[] = "ybounds   = ";
            printBounds(fp, dig, prefix, sizeof(prefix)-1, ydim, nvertex, ybounds);
            Free(ybounds);
	  }

	if ( gridHasArea(gridID) )
	  {
            double *area = (double*) Malloc(gridsize*sizeof(double));
            gridInqArea(gridID, area);
            static const char prefix[] = "area      = ";
            printDblsPrefixAutoBrk(fp, dig, prefix, sizeof(prefix)-1, gridsize, area);
            Free(area);
	  }

        if ( type == GRID_GAUSSIAN_REDUCED )
          {
            static const char prefix[] = "rowlon    = ";
            int *rowlon = (int *)Malloc((size_t)ysize*sizeof(int));
            gridInqRowlon(gridID, rowlon);
            printIntsPrefixAutoBrk(fp, prefix, sizeof(prefix)-1,
                                   (size_t)(ysize > 0 ? ysize : 0), rowlon);
            Free(rowlon);
          }

        if ( type == GRID_PROJECTION ) gridPrintAttributes(fp, gridID);

	break;
      }
    case GRID_LCC:
      {
	double originLon = 0, originLat = 0, lonParY = 0, lat1 = 0, lat2 = 0, xincm = 0, yincm = 0;
	int projflag = 0, scanflag = 0;
	gridInqParamLCC(gridID, &originLon, &originLat, &lonParY, &lat1, &lat2, &xincm, &yincm,
                        &projflag, &scanflag);

	fprintf(fp,
                "originLon = %.*g\n"
                "originLat = %.*g\n"
                "lonParY   = %.*g\n"
                "lat1      = %.*g\n"
                "lat2      = %.*g\n"
                "xinc      = %.*g\n"
                "yinc      = %.*g\n"
                "projection = %s\n",
                dig, originLon, dig, originLat, dig, lonParY,
                dig, lat1, dig, lat2, dig, xincm, dig, yincm,
                (projflag & 128) == 0 ? "northpole" : "southpole");
	break;
      }
    case GRID_SPECTRAL:
      {
        fprintf(fp, "truncation = %d\n"
                "complexpacking = %d\n", gridInqTrunc(gridID), gridInqComplexPacking(gridID) );
        break;
      }
    case GRID_FOURIER:
      {
	fprintf(fp, "truncation = %d\n", gridInqTrunc(gridID));
	break;
      }
    case GRID_GME:
      {
        int nd, ni, ni2, ni3;
        gridInqParamGME(gridID, &nd, &ni, &ni2, &ni3);
        fprintf(fp, "ni        = %d\n", ni );
        break;
      }
   default:
      {
	fprintf(stderr, "Unsupported grid type: %s\n", gridNamePtr(type));
        break;
      }
    }
  /* TODO !!!
  unsigned char uuidOfHGrid[CDI_UUID_SIZE];
  gridInqUUID(gridID, uuidOfHGrid);
  if ( !cdiUUIDIsNull(uuidOfHGrid) )
    {
      char uuidOfHGridStr[37];
      cdiUUID2Str(uuidOfHGrid, uuidOfHGridStr);
      if ( uuidOfHGridStr[0] != 0 && strlen(uuidOfHGridStr) == 36 )
        fprintf(fp, "uuid      = %s\n", uuidOfHGridStr);
    }
  */
  if ( gridInqMask(gridID, NULL) )
    {
      int *mask = (gridsize>0) ? (int*) Malloc((size_t)gridsize*sizeof(int)) : NULL;
      gridInqMask(gridID, mask);
      static const char prefix[] = "mask      = ";
      printMask(fp, prefix, sizeof(prefix)-1,
                (size_t)(gridsize > 0 ? gridsize : 0), mask);
      if ( mask ) Free(mask);
    }
}


void cdo_print_grid(int gridID, int opt)
{
  gridPrintKernel(gridID, opt, stdout);
}
