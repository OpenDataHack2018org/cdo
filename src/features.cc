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

#include <thread> // std::thread::hardware_concurrency()

extern "C" {
size_t getMemorySize(void);
}

void
printFeatures(void)
{
  fprintf(stderr, "Features:");
  size_t memory_size = getMemorySize() / (1024 * 1024 * 1024);
  if (memory_size > 0) fprintf(stderr, " %zuGB", memory_size);
  fprintf(stderr, " %uthreads", std::thread::hardware_concurrency());
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


#include <iostream>
#include <utility>

void
cdoConfig(const char *option)
{
  std::map<std::string, std::pair<std::string, bool>> configMap;

  // clang-format off
  configMap["has-srv"]      = {"SERVICE",          cdiHaveFiletype(CDI_FILETYPE_SRV)};
  configMap["has-ext"]      = {"EXTRA",            cdiHaveFiletype(CDI_FILETYPE_EXT)};
  configMap["has-ieg"]      = {"IEG",              cdiHaveFiletype(CDI_FILETYPE_IEG)};
  configMap["has-grb"]      = {"GRIB 1",           cdiHaveFiletype(CDI_FILETYPE_GRB)};
  configMap["has-grb1"]     = {"GRIB 1",           cdiHaveFiletype(CDI_FILETYPE_GRB)};
  configMap["has-grb2"]     = {"GRIB 2",           cdiHaveFiletype(CDI_FILETYPE_GRB2)};
  configMap["has-nc"]       = {"NetCDF",           cdiHaveFiletype(CDI_FILETYPE_NC)};
  configMap["has-nc2"]      = {"NetCDF 2",         cdiHaveFiletype(CDI_FILETYPE_NC2)};
  configMap["has-nc4"]      = {"NetCDF 4",         cdiHaveFiletype(CDI_FILETYPE_NC4)};
  configMap["has-nc4c"]     = {"NetCDF 4 classic", cdiHaveFiletype(CDI_FILETYPE_NC4C)};
  configMap["has-nc5"]      = {"NetCDF 5",         cdiHaveFiletype(CDI_FILETYPE_NC5)};
  configMap["has-hdf5"]     = {"HDF5",             false};
  configMap["has-cgribex"]  = {"CGRIBEX",          false};
  configMap["has-cmor"]     = {"CMOR",             false};
  configMap["has-proj"]     = {"PROJ",             false};
  configMap["has-threads"]  = {"PTHREADS",         false};
  configMap["has-openmp"]   = {"OPENMP",           false};
  configMap["has-wordexp"]  = {"WORDEXP",          false};
  // clang-format on

#ifdef HAVE_LIBHDF5
  configMap["has-hdf5"].second = true;
#endif

#ifdef HAVE_LIBCGRIBEX
  configMap["has-cgribex"].second = true;
#endif

#ifdef HAVE_LIBCMOR
  configMap["has-cmor"].second = true;
#endif

#ifdef HAVE_LIBPROJ
  configMap["has-proj"].second = true;
#endif

#ifdef HAVE_LIBPTHREAD
  configMap["has-threads"].second = true;
#endif

#ifdef _OPENMP
  configMap["has-openmp"].second = true;
#endif

#ifdef HAVE_WORDEXP_H
  configMap["has-wordexp"].second = true;
#endif

  if (STR_IS_EQ("all-json", option) || STR_IS_EQ("all", option))
    {
      std::cout << "{";
      int i = 0;
      for (auto entry : configMap)
        {
          if (i++) fprintf(stdout, ",");
          std::cout << "\"" << entry.first << "\":\"" << (entry.second.second ? "yes" : "no") << "\"";
        }
      std::cout << "}\n";
    }
  else
    {
      bool foundOption = false;
      for (auto entry : configMap)
        {
          if (STR_IS_EQ(entry.first.c_str(), option))
            {
              foundOption = true;
              std::cout << (entry.second.second ? "yes" : "no") << "\n";
            }
        }
      
      if ( !foundOption )
        {
          fprintf(stdout, "unknown config option: %s\n", option);
          fprintf(stdout, "\n");
          fprintf(stdout, "Available config option:\n");
          fprintf(stdout, "\n");
          for (auto entry : configMap)
            fprintf(stdout, "  %-12s  whether %s is enabled\n", entry.first.c_str(), entry.second.first.c_str());

          exit(EXIT_FAILURE);
        }
    }

  exit(EXIT_SUCCESS);
}
