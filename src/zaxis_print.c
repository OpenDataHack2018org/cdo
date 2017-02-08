#include <cdi.h>
#include "cdi_uuid.h"
#include "cdo_int.h"


static
void printDblsPrefixAutoBrk(FILE *fp, int dig, const char prefix[], int nbyte0, size_t n, const double vals[], size_t extbreak)
{
  fputs(prefix, fp);
  int nbyte = nbyte0;
  for ( size_t i = 0; i < n; i++ )
    {
      if ( nbyte > 80  || (i && i == extbreak) )
        {
          fprintf(fp, "\n%*s", nbyte0, "");
          nbyte = nbyte0;
        }
      nbyte += fprintf(fp, "%.*g ", dig, vals[i]);
    }
  fputs("\n", fp);
}

static
void zaxis_print_kernel(int zaxisID, FILE *fp)
{
  char attstr[CDI_MAX_NAME];
  int type    = zaxisInqType(zaxisID);
  int nlevels = zaxisInqSize(zaxisID);
  int prec    = zaxisInqPrec(zaxisID);
  size_t nvals = (size_t) zaxisInqLevels(zaxisID, NULL);

  int dig = (prec == CDI_DATATYPE_FLT64) ? 15 : 7;

  fprintf(fp, "zaxistype = %s\n", zaxisNamePtr(type));
  fprintf(fp, "size      = %d\n", nlevels);

  if ( nlevels == 1 && zaxisInqScalar(zaxisID) ) fprintf(fp, "scalar    = true\n");

  attstr[0] = 0; cdiZaxisInqKeyStr(zaxisID, CDI_KEY_NAME, CDI_MAX_NAME, attstr);
  if ( attstr[0] )  fprintf(fp, "name      = %s\n", attstr);
  attstr[0] = 0; cdiZaxisInqKeyStr(zaxisID, CDI_KEY_LONGNAME, CDI_MAX_NAME, attstr);
  if ( attstr[0] )  fprintf(fp, "longname  = \"%s\"\n", attstr);
  attstr[0] = 0; cdiZaxisInqKeyStr(zaxisID, CDI_KEY_UNITS, CDI_MAX_NAME, attstr);
  if ( attstr[0] )  fprintf(fp, "units     = \"%s\"\n", attstr);

  double *vals = (double*) Malloc(nvals*sizeof(double));

  if ( nvals )
    {                
      zaxisInqLevels(zaxisID, vals);
      static const char prefix[] = "levels    = ";
      printDblsPrefixAutoBrk(fp, dig, prefix, (int)sizeof(prefix)-1, nvals, vals, 0);
    }

  if ( zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL) )
    {
      {
        zaxisInqLbounds(zaxisID, vals);
        static const char prefix[] = "lbounds   = ";
        printDblsPrefixAutoBrk(fp, dig, prefix, (int)sizeof(prefix)-1, nvals, vals, 0);
      }

      {
        zaxisInqUbounds(zaxisID, vals);
        static const char prefix[] = "ubounds   = ";
        printDblsPrefixAutoBrk(fp, dig, prefix, (int)sizeof(prefix)-1, nvals, vals, 0);
      }
    }

  Free(vals);

  if ( type == ZAXIS_HYBRID || type == ZAXIS_HYBRID_HALF )
    {
      int vctsize = zaxisInqVctSize(zaxisID);
      fprintf(fp, "vctsize   = %d\n", vctsize);
      if ( vctsize )
        {
          double *vct = (double*) Malloc(vctsize*sizeof(double));
          zaxisInqVct(zaxisID, vct);
          static const char prefix[] = "vct       = ";
          printDblsPrefixAutoBrk(fp, dig, prefix, (int)sizeof(prefix)-1, vctsize, vct, vctsize/2);
          Free(vct);
        }
    }

  if ( type == ZAXIS_REFERENCE )
    {
      unsigned char uuid[CDI_UUID_SIZE];
      zaxisInqUUID(zaxisID, uuid);
      if ( !cdiUUIDIsNull(uuid) )
        {
          char uuidStr[37];
          cdiUUID2Str(uuid, uuidStr);
          if ( uuidStr[0] != 0 && strlen(uuidStr) == 36 )
            fprintf(fp, "uuid      = %s\n", uuidStr);
        }
    }
}


void cdo_print_zaxis(int zaxisID)
{
  zaxis_print_kernel(zaxisID, stdout);
}
