/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2018 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/
#include <cdi.h>
#include "cdi_uuid.h"
#include "cdo_int.h"

void cdo_print_attributes(FILE *fp, int cdiID, int varID, int nblanks);

static void
printDblsPrefixAutoBrk(FILE *fp, int dig, const char prefix[], int nbyte0, size_t n, const double vals[], size_t extbreak)
{
  fputs(prefix, fp);
  int nbyte = nbyte0;
  for (size_t i = 0; i < n; i++)
    {
      if (nbyte > 80 || (i && i == extbreak))
        {
          fprintf(fp, "\n%*s", nbyte0, "");
          nbyte = nbyte0;
        }
      nbyte += fprintf(fp, "%.*g ", dig, vals[i]);
    }
  fputs("\n", fp);
}

static void
zaxis_print_kernel(int zaxisID, FILE *fp)
{
  char attstr[CDI_MAX_NAME];
  int type = zaxisInqType(zaxisID);
  int nlevels = zaxisInqSize(zaxisID);
  int prec = zaxisInqDatatype(zaxisID);
  size_t nvals = (size_t) zaxisInqLevels(zaxisID, NULL);

  int dig = (prec == CDI_DATATYPE_FLT64) ? CDO_dbl_digits : CDO_flt_digits;

  fprintf(fp, "zaxistype = %s\n", zaxisNamePtr(type));
  fprintf(fp, "size      = %d\n", nlevels);

  if (nlevels == 1 && zaxisInqScalar(zaxisID)) fprintf(fp, "scalar    = true\n");

  attstr[0] = 0;
  cdiZaxisInqKeyStr(zaxisID, CDI_KEY_NAME, CDI_MAX_NAME, attstr);
  if (attstr[0]) fprintf(fp, "name      = %s\n", attstr);
  attstr[0] = 0;
  cdiZaxisInqKeyStr(zaxisID, CDI_KEY_LONGNAME, CDI_MAX_NAME, attstr);
  if (attstr[0]) fprintf(fp, "longname  = \"%s\"\n", attstr);
  attstr[0] = 0;
  cdiZaxisInqKeyStr(zaxisID, CDI_KEY_UNITS, CDI_MAX_NAME, attstr);
  if (attstr[0]) fprintf(fp, "units     = \"%s\"\n", attstr);

  char **cvals = NULL;
  std::vector<double> vals;
  if (nvals) vals.resize(nvals);

  if (nvals)
    {
      zaxisInqLevels(zaxisID, vals.data());
      static const char prefix[] = "levels    = ";
      printDblsPrefixAutoBrk(fp, dig, prefix, (int) sizeof(prefix) - 1, nvals, vals.data(), 0);
    }
  else if (type == ZAXIS_CHAR)
    {
      int clen = zaxisInqCLen(zaxisID);
      zaxisInqCVals(zaxisID, &cvals);
      fprintf(fp, "levels    = \n");
      for (int i = 0; i < nlevels; i++)
        {
          fprintf(fp, "     [%2d] = %.*s\n", i, clen, cvals[i]);
          Free(cvals[i]);
        }
    }

  if (zaxisInqLbounds(zaxisID, NULL) && zaxisInqUbounds(zaxisID, NULL))
    {
      {
        zaxisInqLbounds(zaxisID, vals.data());
        static const char prefix[] = "lbounds   = ";
        printDblsPrefixAutoBrk(fp, dig, prefix, (int) sizeof(prefix) - 1, nvals, vals.data(), 0);
      }

      {
        zaxisInqUbounds(zaxisID, vals.data());
        static const char prefix[] = "ubounds   = ";
        printDblsPrefixAutoBrk(fp, dig, prefix, (int) sizeof(prefix) - 1, nvals, vals.data(), 0);
      }
    }

  if (cvals) Free(cvals);

  if (type == ZAXIS_HYBRID || type == ZAXIS_HYBRID_HALF)
    {
      int vctsize = zaxisInqVctSize(zaxisID);
      if (vctsize)
        {
          fprintf(fp, "vctsize   = %d\n", vctsize);
          std::vector<double> vct(vctsize);
          zaxisInqVct(zaxisID, vct.data());
          static const char prefix[] = "vct       = ";
          printDblsPrefixAutoBrk(fp, dig, prefix, (int) sizeof(prefix) - 1, vctsize, vct.data(), vctsize / 2);
        }
    }

  if (type == ZAXIS_REFERENCE)
    {
      unsigned char uuid[CDI_UUID_SIZE];
      zaxisInqUUID(zaxisID, uuid);
      if (!cdiUUIDIsNull(uuid))
        {
          char uuidStr[37];
          cdiUUID2Str(uuid, uuidStr);
          if (uuidStr[0] != 0 && strlen(uuidStr) == 36) fprintf(fp, "uuid      = %s\n", uuidStr);
        }
    }

  cdo_print_attributes(fp, zaxisID, CDI_GLOBAL, 0);
}

void
cdo_print_zaxis(int zaxisID)
{
  zaxis_print_kernel(zaxisID, stdout);
}
