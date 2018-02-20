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

#if defined(HAVE_HDF5_H)
#include <hdf5.h>
#endif

#if defined(HAVE_ZLIB_H)
#include <zlib.h>
#endif

#ifdef HAVE_LIBXML2
#include <libxml/xmlversion.h>
#endif

#if defined(HAVE_CURL_CURL_H)
#include <curl/curl.h>
#endif

#if defined(HAVE_PROJ_API_H)
#include <proj_api.h>
#endif

#ifdef HAVE_LIBCMOR
extern "C" {
#include "cmor.h"
}
#endif

#include <stdio.h>
#include <string.h>

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
  fprintf(stderr, " C++%d",
          (int) ((__cplusplus - (__cplusplus / 10000) * 10000) / 100));
#endif
#if defined(HAVE_CF_INTERFACE)
  fprintf(stderr, " Fortran");
#endif
#if defined(ENABLE_DATA)
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
#if defined(HAVE_LIBHDF5)
  fprintf(stderr, " HDF5");
#endif
#if defined(HAVE_NETCDF4)
  fprintf(stderr, " NC4");
#if defined(HAVE_NC4HDF5)
  fprintf(stderr, "/HDF5");
#if defined(HAVE_NC4HDF5_THREADSAFE)
  fprintf(stderr, "/threadsafe");
#endif
#endif
#endif
#if defined(HAVE_LIBNC_DAP)
  fprintf(stderr, " OPeNDAP");
#endif
#if defined(HAVE_LIBSZ)
  fprintf(stderr, " SZ");
#endif
  /*
#if defined(HAVE_LIBZ)
  fprintf(stderr, " Z");
#endif
  */
#if defined(HAVE_LIBUDUNITS2)
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
#if defined(HAVE_LIBDRMAA)
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
#if defined(HAVE_LIBHDF5)
  fprintf(stderr, " HDF5");
#if defined(H5_VERS_MAJOR)
  unsigned h5h_majnum = H5_VERS_MAJOR, h5h_minnum = H5_VERS_MINOR,
           h5h_relnum = H5_VERS_RELEASE;
  fprintf(stderr, "/%u.%u.%u", h5h_majnum, h5h_minnum, h5h_relnum);

  unsigned h5l_majnum, h5l_minnum, h5l_relnum;
  H5get_libversion(&h5l_majnum, &h5l_minnum, &h5l_relnum);
  if ((h5h_majnum != h5l_majnum) || (h5h_minnum != h5l_minnum)
      || (h5h_relnum != h5l_relnum))
    fprintf(stderr, "(%u.%u.%u)", h5l_majnum, h5l_minnum, h5l_relnum);
#endif
#endif
    /*
  #if defined(HAVE_LIBZ)
    {
      fprintf(stderr, " zlib/%s", zlibVersion());
  #if defined(ZLIB_VERSION)
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
#if defined(PJ_VERSION)
  fprintf(stderr, "/%g", PJ_VERSION * 0.01);
#endif
#endif

#ifdef HAVE_LIBCMOR
  fprintf(stderr, " CMOR");
#if defined(CMOR_VERSION_MAJOR)
  fprintf(stderr, "/%u.%u.%u", CMOR_VERSION_MAJOR, CMOR_VERSION_MINOR,
          CMOR_VERSION_PATCH);
#endif
#endif

#ifdef HAVE_LIBXML2
  fprintf(stderr, " xml2");
#if defined(LIBXML_DOTTED_VERSION)
  fprintf(stderr, "/%s", LIBXML_DOTTED_VERSION);
#endif
#endif

#ifdef HAVE_LIBCURL
  {
    curl_version_info_data *version_data = curl_version_info(CURLVERSION_NOW);
    fprintf(stderr, " curl/%s", version_data->version);
#if defined(LIBCURL_VERSION)
    if (strcmp(LIBCURL_VERSION, version_data->version) != 0)
      fprintf(stderr, "(h%s)", LIBCURL_VERSION);
#else
    fprintf(stderr, "(header not found)");
#endif
  }
#endif

  fprintf(stderr, "\n");
}
