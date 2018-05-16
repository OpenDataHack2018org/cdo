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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cdi.h>

#ifdef HAVE_HDF5_H
#include <hdf5.h>
#endif

#ifdef HAVE_ZLIB_H
#include <zlib.h>
#endif

#ifdef HAVE_LIBXML2
#include <libxml/xmlversion.h>
#endif

#ifdef HAVE_CURL_CURL_H
#include <curl/curl.h>
#endif

#ifdef HAVE_PROJ_API_H
#include <proj_api.h>
#endif

#ifdef HAVE_LIBCMOR
extern "C" {
#include "cmor.h"
}
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "cdo_int.h"  // HAVE_OPENMP4

extern "C" {
size_t getMemorySize(void);
}

void
printFeatures(void)
{
  fprintf(stderr, "Features:");
  size_t memory_size = getMemorySize() / (1024 * 1024 * 1024);
  if (memory_size > 0) fprintf(stderr, " %zuGB", memory_size);
#ifdef __cplusplus
  fprintf(stderr, " C++%d", (int) ((__cplusplus - (__cplusplus / 10000) * 10000) / 100));
#endif
#ifdef HAVE_CF_INTERFACE
  fprintf(stderr, " Fortran");
#endif
#ifdef ENABLE_DATA
  fprintf(stderr, " DATA");
#endif
#ifdef HAVE_LIBPTHREAD
  fprintf(stderr, " PTHREADS");
#endif
#ifdef _OPENMP
  fprintf(stderr, " OpenMP");
#if defined(HAVE_OPENMP45)
  fprintf(stderr, "45");
#elif defined(HAVE_OPENMP4)
  fprintf(stderr, "4");
#elif defined(HAVE_OPENMP3)
  fprintf(stderr, "3");
#endif
#endif
#ifdef HAVE_LIBHDF5
  fprintf(stderr, " HDF5");
#endif
#ifdef HAVE_NETCDF4
  fprintf(stderr, " NC4");
#ifdef HAVE_NC4HDF5
  fprintf(stderr, "/HDF5");
#ifdef HAVE_NC4HDF5_THREADSAFE
  fprintf(stderr, "/threadsafe");
#endif
#endif
#endif
#ifdef HAVE_LIBNC_DAP
  fprintf(stderr, " OPeNDAP");
#endif
#ifdef HAVE_LIBSZ
  fprintf(stderr, " SZ");
#endif
/*
#ifdef HAVE_LIBZ
fprintf(stderr, " Z");
#endif
*/
#ifdef HAVE_LIBUDUNITS2
  fprintf(stderr, " UDUNITS2");
#endif
#ifdef HAVE_LIBPROJ
  fprintf(stderr, " PROJ.4");
#endif
#ifdef HAVE_LIBXML2
  fprintf(stderr, " XML2");
#endif
#ifdef HAVE_LIBMAGICS
  fprintf(stderr, " MAGICS");
#endif
#ifdef HAVE_LIBDRMAA
  fprintf(stderr, " DRMAA");
#endif
#ifdef HAVE_LIBCURL
  fprintf(stderr, " CURL");
#endif
#ifdef HAVE_LIBFFTW3
  fprintf(stderr, " FFTW3");
#endif
#ifdef HAVE_LIBCMOR
  fprintf(stderr, " CMOR");
#endif
#if defined(__AVX2__)
  fprintf(stderr, " AVX2");
#elif defined(__AVX__)
  fprintf(stderr, " AVX");
#elif defined(__SSE4_2__)
  fprintf(stderr, " SSE4_2");
#elif defined(__SSE4_1__)
  fprintf(stderr, " SSE4_1");
#elif defined(__SSE3__)
  fprintf(stderr, " SSE3");
#elif defined(__SSE2__)
  fprintf(stderr, " SSE2");
#endif
  fprintf(stderr, "\n");
}

void
printLibraries(void)
{
  fprintf(stderr, "Libraries:");
#ifdef HAVE_LIBHDF5
  fprintf(stderr, " HDF5");
#ifdef H5_VERS_MAJOR
  unsigned h5h_majnum = H5_VERS_MAJOR, h5h_minnum = H5_VERS_MINOR, h5h_relnum = H5_VERS_RELEASE;
  fprintf(stderr, "/%u.%u.%u", h5h_majnum, h5h_minnum, h5h_relnum);

  unsigned h5l_majnum, h5l_minnum, h5l_relnum;
  H5get_libversion(&h5l_majnum, &h5l_minnum, &h5l_relnum);
  if ((h5h_majnum != h5l_majnum) || (h5h_minnum != h5l_minnum) || (h5h_relnum != h5l_relnum))
    fprintf(stderr, "(%u.%u.%u)", h5l_majnum, h5l_minnum, h5l_relnum);
#endif
#endif
/*
#ifdef HAVE_LIBZ
{
  fprintf(stderr, " zlib/%s", zlibVersion());
#ifdef ZLIB_VERSION
  if ( strcmp(ZLIB_VERSION, zlibVersion()) != 0 )
    fprintf(stderr, "(h%s)", ZLIB_VERSION);
#else
  fprintf(stderr, "(header not found)");
#endif
}
#endif
*/
#ifdef HAVE_LIBPROJ
  fprintf(stderr, " proj");
#ifdef PJ_VERSION
  fprintf(stderr, "/%g", PJ_VERSION * 0.01);
#endif
#endif

#ifdef HAVE_LIBCMOR
  fprintf(stderr, " CMOR");
#ifdef CMOR_VERSION_MAJOR
  fprintf(stderr, "/%u.%u.%u", CMOR_VERSION_MAJOR, CMOR_VERSION_MINOR, CMOR_VERSION_PATCH);
#endif
#endif

#ifdef HAVE_LIBXML2
  fprintf(stderr, " xml2");
#ifdef LIBXML_DOTTED_VERSION
  fprintf(stderr, "/%s", LIBXML_DOTTED_VERSION);
#endif
#endif

#ifdef HAVE_LIBCURL
  {
    curl_version_info_data *version_data = curl_version_info(CURLVERSION_NOW);
    fprintf(stderr, " curl/%s", version_data->version);
#ifdef LIBCURL_VERSION
    if (strcmp(LIBCURL_VERSION, version_data->version) != 0) fprintf(stderr, "(h%s)", LIBCURL_VERSION);
#else
    fprintf(stderr, "(header not found)");
#endif
  }
#endif

  fprintf(stderr, "\n");
}

void
cdoConfig(const char *option)
{
  static const char *YN[] = { "no", "yes" };
  const char *has_srv = YN[cdiHaveFiletype(CDI_FILETYPE_SRV)];
  const char *has_ext = YN[cdiHaveFiletype(CDI_FILETYPE_EXT)];
  const char *has_ieg = YN[cdiHaveFiletype(CDI_FILETYPE_IEG)];
  const char *has_grb = YN[cdiHaveFiletype(CDI_FILETYPE_GRB)];
  const char *has_grb2 = YN[cdiHaveFiletype(CDI_FILETYPE_GRB2)];
  const char *has_nc = YN[cdiHaveFiletype(CDI_FILETYPE_NC)];
  const char *has_nc2 = YN[cdiHaveFiletype(CDI_FILETYPE_NC2)];
  const char *has_nc4 = YN[cdiHaveFiletype(CDI_FILETYPE_NC4)];
  const char *has_nc4c = YN[cdiHaveFiletype(CDI_FILETYPE_NC4C)];
  const char *has_nc5 = YN[cdiHaveFiletype(CDI_FILETYPE_NC5)];
  const char *has_hdf5 = YN[0];
  const char *has_cgribex = YN[0];
  const char *has_cmor = YN[0];
  const char *has_proj = YN[0];
  const char *has_threads = YN[0];
  const char *has_openmp = YN[0];
  const char *has_wordexp = YN[0];

#ifdef HAVE_LIBHDF5
  has_hdf5 = YN[1];
#endif

#ifdef HAVE_LIBCGRIBEX
  has_cgribex = YN[1];
#endif

#ifdef HAVE_LIBCMOR
  has_cmor = YN[1];
#endif

#ifdef HAVE_LIBPROJ
  has_proj = YN[1];
#endif

#ifdef HAVE_LIBPTHREAD
  has_threads = YN[1];
#endif

#ifdef _OPENMP
  has_openmp = YN[1];
#endif

#ifdef HAVE_WORDEXP_H
  has_wordexp = YN[1];
#endif

  if (STR_IS_EQ("all-json", option) || STR_IS_EQ("all", option))
    {
      fprintf(stdout, "{");
      fprintf(stdout, "\"has-srv\":\"%s\"", has_srv);
      fprintf(stdout, ",\"has-ext\":\"%s\"", has_ext);
      fprintf(stdout, ",\"has-ieg\":\"%s\"", has_ieg);
      fprintf(stdout, ",\"has-grb\":\"%s\"", has_grb);
      fprintf(stdout, ",\"has-grb1\":\"%s\"", has_grb);
      fprintf(stdout, ",\"has-grb2\":\"%s\"", has_grb2);
      fprintf(stdout, ",\"has-nc\":\"%s\"", has_nc);
      fprintf(stdout, ",\"has-nc2\":\"%s\"", has_nc2);
      fprintf(stdout, ",\"has-nc4\":\"%s\"", has_nc4);
      fprintf(stdout, ",\"has-nc4c\":\"%s\"", has_nc4c);
      fprintf(stdout, ",\"has-nc5\":\"%s\"", has_nc5);
      fprintf(stdout, ",\"has-hdf5\":\"%s\"", has_hdf5);
      fprintf(stdout, ",\"has-cgribex\":\"%s\"", has_cgribex);
      fprintf(stdout, ",\"has-cmor\":\"%s\"", has_cmor);
      fprintf(stdout, ",\"has-proj\":\"%s\"", has_proj);
      fprintf(stdout, ",\"has-threads\":\"%s\"", has_threads);
      fprintf(stdout, ",\"has-openmp\":\"%s\"", has_openmp);
      fprintf(stdout, ",\"has-wordexp\":\"%s\"", has_wordexp);
      fprintf(stdout, "}\n");
    }
  else
    {
      if (STR_IS_EQ("has-srv", option))
        fprintf(stdout, "%s\n", has_srv);
      else if (STR_IS_EQ("has-ext", option))
        fprintf(stdout, "%s\n", has_ext);
      else if (STR_IS_EQ("has-ieg", option))
        fprintf(stdout, "%s\n", has_ieg);
      else if (STR_IS_EQ("has-grb", option))
        fprintf(stdout, "%s\n", has_grb);
      else if (STR_IS_EQ("has-grb1", option))
        fprintf(stdout, "%s\n", has_grb);
      else if (STR_IS_EQ("has-grb2", option))
        fprintf(stdout, "%s\n", has_grb2);
      else if (STR_IS_EQ("has-nc", option))
        fprintf(stdout, "%s\n", has_nc);
      else if (STR_IS_EQ("has-nc2", option))
        fprintf(stdout, "%s\n", has_nc2);
      else if (STR_IS_EQ("has-nc4", option))
        fprintf(stdout, "%s\n", has_nc4);
      else if (STR_IS_EQ("has-nc4c", option))
        fprintf(stdout, "%s\n", has_nc4c);
      else if (STR_IS_EQ("has-nc5", option))
        fprintf(stdout, "%s\n", has_nc5);
      else if (STR_IS_EQ("has-hdf5", option))
        fprintf(stdout, "%s\n", has_hdf5);
      else if (STR_IS_EQ("has-cgribex", option))
        fprintf(stdout, "%s\n", has_cgribex);
      else if (STR_IS_EQ("has-cmor", option))
        fprintf(stdout, "%s\n", has_cmor);
      else if (STR_IS_EQ("has-proj", option))
        fprintf(stdout, "%s\n", has_proj);
      else if (STR_IS_EQ("has-threads", option))
        fprintf(stdout, "%s\n", has_threads);
      else if (STR_IS_EQ("has-openmp", option))
        fprintf(stdout, "%s\n", has_openmp);
      else if (STR_IS_EQ("has-wordexp", option))
        fprintf(stdout, "%s\n", has_wordexp);
      else
        {
          fprintf(stdout, "unknown config option: %s\n", option);
          fprintf(stdout, "\n");
          fprintf(stdout, "Available config option:\n");
          fprintf(stdout, "\n");
          fprintf(stdout, "  has-srv     whether SERVICE is enabled\n");
          fprintf(stdout, "  has-ext     whether EXTRA is enabled\n");
          fprintf(stdout, "  has-ieg     whether IEG is enabled\n");
          fprintf(stdout, "  has-grb     whether GRIB 1 is enabled\n");
          fprintf(stdout, "  has-grb1    whether GRIB 1 is enabled\n");
          fprintf(stdout, "  has-grb2    whether GRIB 2 is enabled\n");
          fprintf(stdout, "  has-nc      whether NetCDF is enabled\n");
          fprintf(stdout, "  has-nc2     whether NetCDF 2 is enabled\n");
          fprintf(stdout, "  has-nc4     whether NetCDF 4 is enabled\n");
          fprintf(stdout, "  has-nc4c    whether NetCDF 4 classic is enabled\n");
          fprintf(stdout, "  has-nc5     whether NetCDF5 is enabled\n");
          fprintf(stdout, "  has-hdf5    whether HDF5 is enabled\n");
          fprintf(stdout, "  has-cgribex whether CGRIBEX is enabled\n");
          fprintf(stdout, "  has-cmor    whether CMOR is enabled\n");
          fprintf(stdout, "  has-proj    whether PROJ is enabled\n");
          fprintf(stdout, "  has-threads whether PTHREADS is enabled\n");
          fprintf(stdout, "  has-openmp  whether OPENMP is enabled\n");
          fprintf(stdout, "  has-wordexp whether WORDEXP is enabled\n");

          exit(EXIT_FAILURE);
        }
    }

  exit(EXIT_SUCCESS);
}
