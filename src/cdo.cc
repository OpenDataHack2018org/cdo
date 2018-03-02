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

#include "module_definitions.h"
#include <iostream>
#include <vector>
#include "operator_help.h"
#if defined(HAVE_EXECINFO_H)
#include <execinfo.h>
#endif

#if defined(HAVE_WORDEXP_H)
#include <wordexp.h>
#endif

#include <signal.h>
#include <fenv.h>
/*#include <malloc.h>*/ /* mallopt and malloc_stats */
#include <sys/stat.h>
#ifdef HAVE_GETRLIMIT
#ifdef HAVE_SYS_RESOURCE_H
#include <sys/time.h>     /* getrlimit */
#include <sys/resource.h> /* getrlimit */
#endif
#endif
#include <unistd.h> /* sysconf, gethostname */
#include <thread>
#include "timer.h"

#if defined(SX)
#define RLIM_T long long
#else
#define RLIM_T rlim_t
#endif

#include <cdi.h>

#include "cdo_int.h"
#include "cdo_task.h"

#include "cdo_getopt.h"
#include "cdoDebugOutput.h"

#ifdef HAVE_LIBPTHREAD
#include "pthread_debug.h"
#endif

#include "modules.h"
#include "error.h"
#include "grid_proj.h"
#include "percentiles.h"
#include "util_wildcards.h"
#include "util_string.h"
#include "process_int.h"
#include "cdoOptions.h"
#include "timer.h"
#include "commandline.h"
#include "text.h"
#include "datetime.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#if !defined(VERSION)
#define VERSION "0.0.1"
#endif

#define MAX_NUM_VARNAMES 256

#include <string>
#include <cstring>

static int Debug = 0;
static int Version = 0;
static int Help = 0;
static int DebugLevel = 0;
static int numThreads = 0;
static int timer_total;
static int CDO_netcdf_hdr_pad = 0;
static int CDO_Rusage = 0;
const char *CDO_username;

extern "C" {
void streamGrbDefDataScanningMode(int scanmode);
}

void gridsearch_set_method(const char *methodstr);

#define PRINT_RLIMIT(resource)                                                           \
  {                                                                                      \
    int status;                                                                          \
    struct rlimit rlim;                                                                  \
    status = getrlimit(resource, &rlim);                                                 \
    if (status == 0)                                                                     \
      {                                                                                  \
        if (sizeof(RLIM_T) > sizeof(long))                                               \
          {                                                                              \
            fprintf(stderr, "CUR %-15s = %llu\n", #resource, (long long) rlim.rlim_cur); \
            fprintf(stderr, "MAX %-15s = %llu\n", #resource, (long long) rlim.rlim_max); \
          }                                                                              \
        else                                                                             \
          {                                                                              \
            fprintf(stderr, "CUR %-15s = %lu\n", #resource, (long) rlim.rlim_cur);       \
            fprintf(stderr, "MAX %-15s = %lu\n", #resource, (long) rlim.rlim_max);       \
          }                                                                              \
      }                                                                                  \
  }

#define ITSME (STR_IS_EQ(CDO_username, "\x6d\x32\x31\x34\x30\x30\x33"))

static void
cdo_stackframe(void)
{
#if defined HAVE_EXECINFO_H && defined HAVE_BACKTRACE
  void *callstack[32];
  int frames = backtrace(callstack, 32);
  char **messages = backtrace_symbols(callstack, frames);

  fprintf(stderr, "[bt] Execution path:\n");
  if (messages)
    {
      for (int i = 0; i < frames; ++i) fprintf(stderr, "[bt] %s\n", messages[i]);
      free(messages);
    }
#endif
}

static int
cdo_feenableexcept(int excepts)
{
#if defined HAVE_FEENABLEEXCEPT
  int feenableexcept(int);
  int old_excepts = feenableexcept(excepts);
  return old_excepts;
#else
  static fenv_t fenv;
  unsigned new_excepts = ((unsigned) excepts) & FE_ALL_EXCEPT;
  int old_excepts = -1;  // previous masks

  if (fegetenv(&fenv)) return -1;
#if defined(HAVE_FENV_T___CONTROL) && defined(HAVE_FENV_T___MXCSR)
  old_excepts = (int) (fenv.__control & FE_ALL_EXCEPT);

  // unmask
  fenv.__control &= ~new_excepts;
  fenv.__mxcsr &= ~(new_excepts << 7);
#endif

  return (fesetenv(&fenv) ? -1 : (int) old_excepts);
#endif
}

static void
cdo_sig_handler(int signo)
{
  if (signo == SIGFPE)
    {
      cdo_stackframe();
      cdoAbort("floating-point exception!");
    }
}

static void
cdo_set_digits(const char *optarg)
{
  char *ptr1 = 0;
  if (optarg != 0 && (int) strlen(optarg) > 0 && optarg[0] != ',') CDO_flt_digits = (int) strtol(optarg, &ptr1, 10);

  if (CDO_flt_digits < 1 || CDO_flt_digits > 20)
    cdoAbort("Unreasonable value for float significant digits: %d", CDO_flt_digits);

  if (ptr1 && *ptr1 == ',')
    {
      char *ptr2 = 0;
      CDO_dbl_digits = (int) strtol(ptr1 + 1, &ptr2, 10);
      if (ptr2 == ptr1 + 1 || CDO_dbl_digits < 1 || CDO_dbl_digits > 20)
        cdoAbort("Unreasonable value for double significant digits: %d", CDO_dbl_digits);
    }
}

static void
cdo_version(void)
{
  const int filetypes[] = { CDI_FILETYPE_SRV, CDI_FILETYPE_EXT, CDI_FILETYPE_IEG, CDI_FILETYPE_GRB,  CDI_FILETYPE_GRB2,
                            CDI_FILETYPE_NC,  CDI_FILETYPE_NC2, CDI_FILETYPE_NC4, CDI_FILETYPE_NC4C, CDI_FILETYPE_NC5 };
  const char *typenames[] = { "srv", "ext", "ieg", "grb1", "grb2", "nc1", "nc2", "nc4", "nc4c", "nc5" };

  fprintf(stderr, "%s\n", CDO_version);
#ifdef SYSTEM_TYPE
  fprintf(stderr, "System: %s\n", SYSTEM_TYPE);
#endif
#if defined(CXX_COMPILER)
  fprintf(stderr, "CXX Compiler: %s\n", CXX_COMPILER);
#if defined(CXX_VERSION)
  fprintf(stderr, "CXX version : %s\n", CXX_VERSION);
#endif
#endif
#if defined(C_COMPILER)
  fprintf(stderr, "C Compiler: %s\n", C_COMPILER);
#if defined(C_VERSION)
  fprintf(stderr, "C version : %s\n", C_VERSION);
#endif
#endif
#if defined(F77_COMPILER)
  fprintf(stderr, "F77 Compiler: %s\n", F77_COMPILER);
#if defined(F77_VERSION)
  fprintf(stderr, "F77 version : %s\n", F77_VERSION);
#endif
#endif

  printFeatures();
  printLibraries();

  fprintf(stderr, "Filetypes: ");
  set_text_color(stderr, BRIGHT, GREEN);
  for (size_t i = 0; i < sizeof(filetypes) / sizeof(int); ++i)
    if (cdiHaveFiletype(filetypes[i])) fprintf(stderr, "%s ", typenames[i]);
  reset_text_color(stderr);
  fprintf(stderr, "\n");

  cdiPrintVersion();
  fprintf(stderr, "\n");
}

static void
cdo_usage(void)
{
  const char *name;

  /*  fprintf(stderr, "%s\n", CDO_version);*/
  /*  fprintf(stderr, "\n");*/
  fprintf(stderr, "usage : cdo  [Options]  Operator1  [-Operator2  [-OperatorN]]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  Options:\n");
  set_text_color(stderr, RESET, BLUE);
  fprintf(stderr, "    -a             Generate an absolute time axis\n");
  fprintf(stderr, "    -b <nbits>     Set the number of bits for the output precision\n");
  fprintf(stderr, "                   (I8/I16/I32/F32/F64 for "
                  "nc1/nc2/nc4/nc4c/nc5; F32/F64 for grb2/srv/ext/ieg; P1 - "
                  "P24 for grb1/grb2)\n");
  fprintf(stderr, "                   Add L or B to set the byteorder to "
                  "Little or Big endian\n");
  fprintf(stderr, "    --cmor         CMOR conform NetCDF output\n");
  fprintf(stderr, "    -C, --color    Colorized output messages\n");
  fprintf(stderr, "    --enableexcept <except>\n");
  fprintf(stderr, "                   Set individual floating-point traps "
                  "(DIVBYZERO, INEXACT, INVALID, OVERFLOW, UNDERFLOW, "
                  "ALL_EXCEPT)\n");
  fprintf(stderr, "    -f, --format <format>\n");
  fprintf(stderr, "                   Format of the output file. "
                  "(grb1/grb2/nc1/nc2/nc4/nc4c/nc5/srv/ext/ieg)\n");
  fprintf(stderr, "    -g <grid>      Set default grid name or file. Available grids: \n");
  fprintf(stderr, "                   n<N>, t<RES>, tl<RES>, global_<DXY>, "
                  "r<NX>x<NY>, g<NX>x<NY>, gme<NI>, lon=<LON>/lat=<LAT>\n");
  fprintf(stderr, "    -h, --help     Help information for the operators\n");
  fprintf(stderr, "    --history      Do not append to NetCDF \"history\" "
                  "global attribute\n");
  fprintf(stderr, "    --netcdf_hdr_pad, --hdr_pad, --header_pad <nbr>\n");
  fprintf(stderr, "                   Pad NetCDF output header with nbr bytes\n");
  /*
  fprintf(stderr, "    -i <inst>      Institution name/file\n");
  fprintf(stderr, "                   Predefined instituts: ");
  for ( int id = 0; id < institutInqNumber; id++ )
    if ( (name = institutInqNamePtr(id)) )
      fprintf(stderr, " %s", name);
  fprintf(stderr, "\n");
  */
  /* fprintf(stderr, "    -l <level>     Level file\n"); */
  fprintf(stderr, "    -k <chunktype> NetCDF4 chunk type: auto, grid or lines\n");
  fprintf(stderr, "    -L             Lock IO (sequential access)\n");
  fprintf(stderr, "    -M             Switch to indicate that the I/O streams "
                  "have missing values\n");
  fprintf(stderr, "    -m <missval>   Set the missing value of non NetCDF files "
                  "(default: %g)\n",
          cdiInqMissval());
  fprintf(stderr, "    --no_warnings  Inhibit warning messages\n");
  fprintf(stderr, "    -O             Overwrite existing output file, if checked\n");
  fprintf(stderr, "    --operators    List of all operators\n");
#ifdef _OPENMP
  fprintf(stderr, "    -P <nthreads>  Set number of OpenMP threads\n");
#endif
  fprintf(stderr, "    --percentile <method>\n");
  fprintf(stderr, "                   Percentile method: nrank, nist, numpy, "
                  "numpy_lower, numpy_higher, numpy_nearest\n");
  fprintf(stderr, "    --precision <float_digits[,double_digits]>\n");
  fprintf(stderr, "                   Precision to use in displaying "
                  "floating-point data (default: 7,15)\n");
  fprintf(stderr, "    --reduce_dim   Reduce NetCDF dimensions\n");
  if (ITSME) fprintf(stderr, "    --remap_genweights 0/1\n");
  fprintf(stderr, "    -R, --regular  Convert GRIB1 data from reduced to "
                  "regular grid (cgribex only)\n");
  fprintf(stderr, "    -r             Generate a relative time axis\n");
  fprintf(stderr, "    -S             Create an extra output stream for the "
                  "module TIMSTAT. This stream\n");
  fprintf(stderr, "                   contains the number of non missing "
                  "values for each output period.\n");
  fprintf(stderr, "    -s, --silent   Silent mode\n");
  fprintf(stderr, "    --sortname     Alphanumeric sorting of NetCDF parameter names\n");
  fprintf(stderr, "    -t <codetab>   Set GRIB1 default parameter code table "
                  "name or file (cgribex only)\n");
  fprintf(stderr, "                   Predefined tables: ");
  for (int id = 0; id < tableInqNumber(); id++)
    if ((name = tableInqNamePtr(id))) fprintf(stderr, " %s", name);
  fprintf(stderr, "\n");

  fprintf(stderr, "    --timestat_date <srcdate>\n");
  fprintf(stderr, "                   Target timestamp (time statistics): "
                  "first, middle, midhigh or last source timestep.\n");
  fprintf(stderr, "    -V, --version  Print the version number\n");
  fprintf(stderr, "    -v, --verbose  Print extra details for some operators\n");
  fprintf(stderr, "    -W             Print extra warning messages\n");
  fprintf(stderr, "    -z szip        SZIP compression of GRIB1 records\n");
  fprintf(stderr, "       aec         AEC compression of GRIB2 records\n");
  fprintf(stderr, "       jpeg        JPEG compression of GRIB2 records\n");
  fprintf(stderr, "        zip[_1-9]  Deflate compression of NetCDF4 variables\n");
#ifdef HIRLAM_EXTENSIONS
  fprintf(stderr, "    --Dkext <debLev>   Setting debugLevel for extensions\n");
  fprintf(stderr, "    --outputGribDataScanningMode <mode>   Setting grib "
                  "scanning mode for data in output file <0, 64, 96>; Default "
                  "is 64\n");
#endif  // HIRLAM_EXTENSIONS
  reset_text_color(stderr);
  fprintf(stderr, "\n");

  fprintf(stderr, "  Operators:\n");
  fprintf(stderr, "    Use option --operators for a list of all operators.\n");
  /*
  set_text_color(stderr, RESET, GREEN);
  operatorPrintAll();
  reset_text_color(stderr);
  */

  fprintf(stderr, "\n");
  fprintf(stderr, "  CDO version %s, Copyright (C) 2003-2018 Uwe Schulzweida\n", VERSION);
  //  fprintf(stderr, "  Available from <http://mpimet.mpg.de/cdo>\n");
  fprintf(stderr, "  This is free software and comes with ABSOLUTELY NO WARRANTY\n");
  fprintf(stderr, "  Report bugs to <http://mpimet.mpg.de/cdo>\n");
}

static void
cdo_init_is_tty(void)
{
  struct stat statbuf;
  fstat(0, &statbuf);
  if (S_ISCHR(statbuf.st_mode)) stdin_is_tty = 1;
  fstat(1, &statbuf);
  if (S_ISCHR(statbuf.st_mode)) stdout_is_tty = 1;
  fstat(2, &statbuf);
  if (S_ISCHR(statbuf.st_mode)) stderr_is_tty = 1;
}

static void
cdoPrintHelp(std::vector<std::string> help /*, char *xoperator*/)
{
  if (help.empty())
    fprintf(stderr, "No help available for this operator!\n");
  else
    {
      bool lprint;
      for (unsigned long i = 0; i < help.size(); i++)
        {
          lprint = !(help[i][0] == '\0' && help[i + 1][0] == ' ');

          if (lprint)
            {
              if (COLOR_STDOUT)
                {
                  if ((help[i].compare("NAME") == 0) || (help[i].compare("SYNOPSIS") == 0)
                      || (help[i].compare("DESCRIPTION") == 0) || (help[i].compare("OPERATORS") == 0)
                      || (help[i].compare("NAMELIST") == 0) || (help[i].compare("PARAMETER") == 0)
                      || (help[i].compare("ENVIRONMENT") == 0) || (help[i].compare("NOTE") == 0)
                      || (help[i].compare("EXAMPLES") == 0))
                    {
                      set_text_color(stdout, BRIGHT, BLACK);
                      fprintf(stdout, "%s", help[i].c_str());
                      reset_text_color(stdout);
                      fprintf(stdout, "\n");
                    }
                  else
                    fprintf(stdout, "%s\n", help[i].c_str());
                }
              else
                {
                  fprintf(stdout, "%s\n", help[i].c_str());
                }
            }
        }
    }
}

#undef IsBigendian
#define IsBigendian() (u_byteorder.c[sizeof(long) - 1])

static void
setDefaultDataType(const char *datatypestr)
{
  static union
  {
    unsigned long l;
    unsigned char c[sizeof(long)];
  } u_byteorder = { 1 };
  int nbits = -1;
  enum
  {
    D_UINT,
    D_INT,
    D_FLT,
    D_CPX
  };
  int dtype = -1;

  int datatype = tolower(*datatypestr);
  if (datatype == 'i')
    {
      dtype = D_INT;
      datatypestr++;
    }
  else if (datatype == 'u')
    {
      dtype = D_UINT;
      datatypestr++;
    }
  else if (datatype == 'f')
    {
      dtype = D_FLT;
      datatypestr++;
    }
  else if (datatype == 'c')
    {
      dtype = D_CPX;
      datatypestr++;
    }
  else if (datatype == 'p')
    {
      datatypestr++;
    }

  if (isdigit((int) *datatypestr))
    {
      nbits = atoi(datatypestr);
      datatypestr += 1;
      if (nbits >= 10) datatypestr += 1;

      if (dtype == -1)
        {
          if (nbits > 0 && nbits < 32)
            cdoDefaultDataType = nbits;
          else if (nbits == 32)
            {
              if (cdoDefaultFileType == CDI_FILETYPE_GRB)
                cdoDefaultDataType = CDI_DATATYPE_PACK32;
              else
                cdoDefaultDataType = CDI_DATATYPE_FLT32;
            }
          else if (nbits == 64)
            cdoDefaultDataType = CDI_DATATYPE_FLT64;
          else
            {
              fprintf(stderr, "Unsupported number of bits %d!\n", nbits);
              fprintf(stderr, "Use I8/I16/I32/F32/F64 for "
                              "nc1/nc2/nc4/nc4c/nc5; F32/F64 for "
                              "grb2/srv/ext/ieg; P1 - P24 for grb1/grb2.\n");
              exit(EXIT_FAILURE);
            }
        }
      else
        {
          if (dtype == D_INT)
            {
              if (nbits == 8)
                cdoDefaultDataType = CDI_DATATYPE_INT8;
              else if (nbits == 16)
                cdoDefaultDataType = CDI_DATATYPE_INT16;
              else if (nbits == 32)
                cdoDefaultDataType = CDI_DATATYPE_INT32;
              else
                {
                  fprintf(stderr, "Unsupported number of bits = %d for datatype INT!\n", nbits);
                  exit(EXIT_FAILURE);
                }
            }
          else if (dtype == D_UINT)
            {
              if (nbits == 8)
                cdoDefaultDataType = CDI_DATATYPE_UINT8;
              else if (nbits == 16)
                cdoDefaultDataType = CDI_DATATYPE_UINT16;
              else if (nbits == 32)
                cdoDefaultDataType = CDI_DATATYPE_UINT32;
              else
                {
                  fprintf(stderr, "Unsupported number of bits = %d for datatype UINT!\n", nbits);
                  exit(EXIT_FAILURE);
                }
            }
          else if (dtype == D_FLT)
            {
              if (nbits == 32)
                cdoDefaultDataType = CDI_DATATYPE_FLT32;
              else if (nbits == 64)
                cdoDefaultDataType = CDI_DATATYPE_FLT64;
              else
                {
                  fprintf(stderr, "Unsupported number of bits = %d for datatype FLT!\n", nbits);
                  exit(EXIT_FAILURE);
                }
            }
          else if (dtype == D_CPX)
            {
              if (nbits == 32)
                cdoDefaultDataType = CDI_DATATYPE_CPX32;
              else if (nbits == 64)
                cdoDefaultDataType = CDI_DATATYPE_CPX64;
              else
                {
                  fprintf(stderr, "Unsupported number of bits = %d for datatype CPX!\n", nbits);
                  exit(EXIT_FAILURE);
                }
            }
        }
    }

  if (*datatypestr != 0)
    {
      if (*datatypestr == 'l' || *datatypestr == 'L')
        {
          if (IsBigendian()) cdoDefaultByteorder = CDI_LITTLEENDIAN;
          datatypestr++;
        }
      else if (*datatypestr == 'b' || *datatypestr == 'B')
        {
          if (!IsBigendian()) cdoDefaultByteorder = CDI_BIGENDIAN;
          datatypestr++;
        }
      else
        {
          fprintf(stderr, "Unsupported character in number of bytes: >%s< !\n", datatypestr);
          exit(EXIT_FAILURE);
        }
    }
}
/*
static
void setDefaultDataTypeByte(char *datatypestr)
{
  static union {unsigned long l; unsigned char c[sizeof(long)];} u_byteorder =
{1}; int datatype = -1;

  if ( isdigit((int) *datatypestr) )
    {
      datatype = atoi(datatypestr);
      datatypestr++;

      if      ( datatype == 1 ) cdoDefaultDataType = CDI_DATATYPE_PACK8;
      else if ( datatype == 2 ) cdoDefaultDataType = CDI_DATATYPE_PACK16;
      else if ( datatype == 3 ) cdoDefaultDataType = CDI_DATATYPE_PACK24;
      else if ( datatype == 4 ) cdoDefaultDataType = CDI_DATATYPE_FLT32;
      else if ( datatype == 8 ) cdoDefaultDataType = CDI_DATATYPE_FLT64;
      else
        {
          fprintf(stderr, "Unsupported datatype %d!\n", datatype);
          fprintf(stderr, "Use 4/8 for filetype nc/srv/ext/ieg and 1/2/3 for
grb1/grb2.\n"); exit(EXIT_FAILURE);
        }
    }

  if ( *datatypestr != 0 )
    {
      if ( *datatypestr == 'l' || *datatypestr == 'L' )
        {
          if ( IsBigendian() ) cdoDefaultByteorder = CDI_LITTLEENDIAN;
          datatypestr++;setDefaultDataTypeByte
        }
      else if ( *datatypestr == 'b' || *datatypestr == 'B' )
        {
          if ( ! IsBigendian() ) cdoDefaultByteorder = CDI_BIGENDIAN;
          datatypestr++;
        }
      else
        {
          fprintf(stderr, "Unsupported character in number of bytes: %s!\n",
datatypestr); exit(EXIT_FAILURE);
        }
    }
}
*/
static void
setDefaultFileType(const char *filetypestr, int labort)
{
  if (filetypestr)
    {
      const char *ftstr = filetypestr;
      size_t len;

      // clang-format off
      if      ( cmpstrlen(filetypestr, "grb2", len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_GRB2;}
      else if ( cmpstrlen(filetypestr, "grb1", len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_GRB; }
      else if ( cmpstrlen(filetypestr, "grb",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_GRB; }
      else if ( cmpstrlen(filetypestr, "nc2",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_NC2; }
      else if ( cmpstrlen(filetypestr, "nc4c", len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_NC4C;}
      else if ( cmpstrlen(filetypestr, "nc4",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_NC4; }
      else if ( cmpstrlen(filetypestr, "nc5",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_NC5; }
      else if ( cmpstrlen(filetypestr, "nc1",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_NC;  }
      else if ( cmpstrlen(filetypestr, "nc",   len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_NC2; }
      else if ( cmpstrlen(filetypestr, "srv",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_SRV; }
      else if ( cmpstrlen(filetypestr, "ext",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_EXT; }
      else if ( cmpstrlen(filetypestr, "ieg",  len)  == 0 ) { ftstr += len; cdoDefaultFileType = CDI_FILETYPE_IEG; }
      else
        {
          if ( labort )
            {
              fprintf(stderr, "Unsupported filetype %s!\n", filetypestr);
              fprintf(stderr, "Available filetypes: grb1/grb2/nc1/nc2/nc4/nc4c/nc5/srv/ext/ieg\n");
              exit(EXIT_FAILURE);
            }
          else
            {
              return;
            }
        }
      // clang-format on

      if (cdoDefaultFileType != CDI_UNDEFID && *ftstr != 0)
        {
          if (*ftstr == '_')
            {
              setDefaultDataType(++ftstr);
            }
          else
            {
              fprintf(stderr, "Unexpected character >%c< in file type >%s<!\n", *ftstr, filetypestr);
              fprintf(stderr, "Use format[_nbits] with:\n");
              fprintf(stderr, "    format = grb1, grb2, nc1, nc2, nc4, nc4c, "
                              "nc5, srv, ext or ieg\n");
              fprintf(stderr, "    nbits  = 32/64 for "
                              "grb2/nc1/nc2/nc4/nc4c/nc5/srv/ext/ieg; 1 - 24 "
                              "for grb1/grb2\n");
              exit(EXIT_FAILURE);
            }
        }
    }
}

#define NTESTS 11
#include <inttypes.h>
static int
getMemAlignment(void)
{
  int ma = -1;
  double *ptr[NTESTS];
  int64_t iptr;
  size_t tsize[NTESTS] = { 1, 3, 5, 9, 17, 33, 69, 121, 251, 510, 1025 };
  size_t ma_check[4] = { 8, 16, 32, 64 };
  int ma_result[4] = { 1, 1, 1, 1 };

  for (int i = 0; i < NTESTS; ++i)
    {
      ptr[i] = (double *) malloc(tsize[i]);
      iptr = (int64_t) ptr[i];
      for (int k = 0; k < 4; ++k)
        if (iptr % ma_check[k]) ma_result[k] = 0;
    }
  for (int i = 0; i < NTESTS; ++i) free(ptr[i]);

  for (int i = NTESTS - 1; i >= 0; i--)
    {
      ptr[i] = (double *) malloc(tsize[i] + 5);
      iptr = (int64_t) ptr[i];
      for (int k = 0; k < 4; ++k)
        if (iptr % ma_check[k]) ma_result[k] = 0;
    }
  for (int i = 0; i < NTESTS; ++i) free(ptr[i]);

  for (int k = 0; k < 4; ++k)
    if (ma_result[k]) ma = ma_check[k];

  return ma;
}

static void
defineCompress(const char *arg)
{
  size_t len = strlen(arg);

  if (strncmp(arg, "szip", len) == 0)
    {
      Options::cdoCompType = CDI_COMPRESS_SZIP;
      Options::cdoCompLevel = 0;
    }
  else if (strncmp(arg, "aec", len) == 0)
    {
      Options::cdoCompType = CDI_COMPRESS_AEC;
      Options::cdoCompLevel = 0;
    }
  else if (strncmp(arg, "jpeg", len) == 0)
    {
      Options::cdoCompType = CDI_COMPRESS_JPEG;
      Options::cdoCompLevel = 0;
    }
  else if (strncmp(arg, "zip", 3) == 0)
    {
      Options::cdoCompType = CDI_COMPRESS_ZIP;
      if (len == 5 && arg[3] == '_' && isdigit(arg[4]))
        Options::cdoCompLevel = atoi(&arg[4]);
      else
        Options::cdoCompLevel = 1;
    }
  else
    {
      fprintf(stderr, "Compression type '%s' unsupported!\n", arg);
      exit(EXIT_FAILURE);
    }
}

static void
defineChunktype(const char *arg)
{
  if (STR_IS_EQ("auto", arg))
    cdoChunkType = CDI_CHUNK_AUTO;
  else if (STR_IS_EQ("grid", arg))
    cdoChunkType = CDI_CHUNK_GRID;
  else if (STR_IS_EQ("lines", arg))
    cdoChunkType = CDI_CHUNK_LINES;
  else
    {
      fprintf(stderr, "Chunk type '%s' unsupported!\n", arg);
      exit(EXIT_FAILURE);
    }
}

static void
defineVarnames(const char *arg)
{
  size_t len = strlen(arg);
  size_t istart = 0;
  while (istart < len && (arg[istart] == ' ' || arg[istart] == ',')) istart++;

  len -= istart;

  if (len)
    {
      cdoVarnames = (char **) Malloc(MAX_NUM_VARNAMES * sizeof(char *));

      char *pbuf = strdup(arg + istart);
      cdoVarnames[cdoNumVarnames++] = pbuf;

      char *commapos = pbuf;
      while ((commapos = strchr(commapos, ',')) != NULL)
        {
          *commapos++ = '\0';
          if (strlen(commapos))
            {
              if (cdoNumVarnames >= MAX_NUM_VARNAMES) cdoAbort("Too many variable names (limit=%d)!", MAX_NUM_VARNAMES);

              cdoVarnames[cdoNumVarnames++] = commapos;
            }
        }
      /*
      for ( int i = 0; i < cdoNumVarnames; ++i )
        printf("varname %d: %s\n", i+1, cdoVarnames[i]);
      */
    }
}

static void
get_env_vars(void)
{
  CDO_username = getenv("LOGNAME");
  if (CDO_username == NULL)
    {
      CDO_username = getenv("USER");
      if (CDO_username == NULL) CDO_username = "unknown";
    }

  char *envstr = getenv("CDO_GRID_SEARCH_DIR");
  if (envstr)
    {
      size_t len = strlen(envstr);
      if (len > 0)
        {
          len += 2;
          cdoGridSearchDir = (char *) Malloc(len);
          memcpy(cdoGridSearchDir, envstr, len - 1);
          if (cdoGridSearchDir[len - 3] != '/')
            {
              cdoGridSearchDir[len - 2] = '/';
              cdoGridSearchDir[len - 1] = 0;
            }
        }
    }

  envstr = getenv("CDO_DISABLE_HISTORY");
  if (envstr)
    {
      if (atoi(envstr) == 1)
        {
          CDO_Reset_History = TRUE;
          if (cdoVerbose) fprintf(stderr, "CDO_DISABLE_HISTORY = %s\n", envstr);
        }
    }

  envstr = getenv("CDO_RESET_HISTORY");
  if (envstr)
    {
      if (atoi(envstr) == 1)
        {
          CDO_Reset_History = TRUE;
          if (cdoVerbose) fprintf(stderr, "CDO_RESET_HISTORY = %s\n", envstr);
        }
    }

  envstr = getenv("CDO_HISTORY_INFO");
  if (envstr)
    {
      int ival = atoi(envstr);
      if (ival == 0 || ival == 1)
        {
          CDO_Append_History = ival;
          if (cdoVerbose) fprintf(stderr, "CDO_HISTORY_INFO = %s\n", envstr);
        }
    }

  CDO_File_Suffix[0] = 0;

  envstr = getenv("CDO_FILE_SUFFIX");
  if (envstr)
    {
      if (envstr[0])
        {
          strncat(CDO_File_Suffix, envstr, sizeof(CDO_File_Suffix) - 1);
          if (cdoVerbose) fprintf(stderr, "CDO_FILE_SUFFIX = %s\n", envstr);
        }
    }

  envstr = getenv("CDO_DISABLE_FILESUFFIX");
  if (envstr)
    {
      if (atoi(envstr) == 1)
        {
          strcat(CDO_File_Suffix, "NULL");
          if (cdoVerbose) fprintf(stderr, "CDO_DISABLE_FILESUFFIX = %s\n", envstr);
        }
    }

  envstr = getenv("CDO_DIAG");
  if (envstr)
    {
      if (atoi(envstr) == 1)
        {
          cdoDiag = TRUE;
          if (cdoVerbose) fprintf(stderr, "CDO_DIAG = %s\n", envstr);
        }
    }

  envstr = getenv("CDO_USE_FFTW");
  if (envstr)
    {
      int ival = atoi(envstr);
      if (ival == 0 || ival == 1)
        {
          CDO_Use_FFTW = ival;
          if (cdoVerbose) fprintf(stderr, "CDO_Use_FFTW = %s\n", envstr);
        }
    }

  envstr = getenv("CDO_VERSION_INFO");
  if (envstr)
    {
      int ival = atoi(envstr);
      if (ival == 0 || ival == 1)
        {
          CDO_Version_Info = ival;
          if (cdoVerbose) fprintf(stderr, "CDO_Version_Info = %s\n", envstr);
        }
    }

  envstr = getenv("CDO_COLOR");
  if (envstr)
    {
      int ival = atoi(envstr);
      if (ival == 0 || ival == 1)
        {
          CDO_Color = ival;
          if (cdoVerbose) fprintf(stderr, "CDO_COLOR = %s\n", envstr);
        }
    }
  else if (CDO_Color == FALSE && ITSME)
    CDO_Color = TRUE;
}

static void
print_system_info()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "CDO_Color           = %d\n", CDO_Color);
  fprintf(stderr, "CDO_Reset_History   = %d\n", CDO_Reset_History);
  fprintf(stderr, "CDO_File_Suffix     = %s\n", CDO_File_Suffix);
  fprintf(stderr, "cdoDefaultFileType  = %d\n", cdoDefaultFileType);
  fprintf(stderr, "cdoDefaultDataType  = %d\n", cdoDefaultDataType);
  fprintf(stderr, "cdoDefaultByteorder = %d\n", cdoDefaultByteorder);
  fprintf(stderr, "cdoDefaultTableID   = %d\n", cdoDefaultTableID);
  fprintf(stderr, "\n");

  const char *envstr;
  envstr = getenv("HOSTTYPE");
  if (envstr) fprintf(stderr, "HOSTTYPE            = %s\n", envstr);
  envstr = getenv("VENDOR");
  if (envstr) fprintf(stderr, "VENDOR              = %s\n", envstr);
  envstr = getenv("OSTYPE");
  if (envstr) fprintf(stderr, "OSTYPE              = %s\n", envstr);
  envstr = getenv("MACHTYPE");
  if (envstr) fprintf(stderr, "MACHTYPE            = %s\n", envstr);
  fprintf(stderr, "\n");

#if defined(_ARCH_PWR6)
  fprintf(stderr, "Predefined: _ARCH_PWR6\n");
#elif defined(_ARCH_PWR7)
  fprintf(stderr, "Predefined: _ARCH_PWR7\n");
#endif

#if defined(__AVX2__)
  fprintf(stderr, "Predefined: __AVX2__\n");
#elif defined(__AVX__)
  fprintf(stderr, "Predefined: __AVX__\n");
#elif defined(__SSE4_2__)
  fprintf(stderr, "Predefined: __SSE4_2__\n");
#elif defined(__SSE4_1__)
  fprintf(stderr, "Predefined: __SSE4_1__\n");
#elif defined(__SSE3__)
  fprintf(stderr, "Predefined: __SSE3__\n");
#elif defined(__SSE2__)
  fprintf(stderr, "Predefined: __SSE2__\n");
#endif
  fprintf(stderr, "\n");

  fprintf(stderr, "sizeof(size_t)      = %zu\n", sizeof(size_t));
  fprintf(stderr, "mem alignment       = %d\n\n", getMemAlignment());

#if defined(HAVE_MMAP)
  fprintf(stderr, "HAVE_MMAP\n");
#endif
#if defined(HAVE_MEMORY_H)
  fprintf(stderr, "HAVE_MEMORY_H\n");
#endif
  fprintf(stderr, "\n");

#if defined(_OPENACC)
  fprintf(stderr, "OPENACC VERSION     = %d\n", _OPENACC);
#endif
/* OPENMP 3:  201107 */
/* OPENMP 4:  201307 gcc 4.9 */
#ifdef _OPENMP
  fprintf(stderr, "OPENMP VERSION      = %d\n", _OPENMP);
#endif
#if defined(__cplusplus)
  fprintf(stderr, "__cplusplus         = %ld\n", __cplusplus);
#endif
#if defined(__GNUC__)
  fprintf(stderr, "GNUC VERSION        = %d\n", __GNUC__);
#endif
#if defined(__GNUC_MINOR__)
  fprintf(stderr, "GNUC MINOR          = %d\n", __GNUC_MINOR__);
#endif
#if defined(__ICC)
  fprintf(stderr, "ICC VERSION         = %d\n", __ICC);
#endif
#if defined(__STDC__)
  fprintf(stderr, "STD ANSI C          = %d\n", __STDC__);
#endif
#if defined(__STD_VERSION__)
  fprintf(stderr, "STD VERSION         = %ld\n", __STD_VERSION__);
#endif
#if defined(__STDC_VERSION__)
  fprintf(stderr, "STDC VERSION        = %ld\n", __STDC_VERSION__);
#endif
#if defined(__STD_HOSTED__)
  fprintf(stderr, "STD HOSTED          = %d\n", __STD_HOSTED__);
#endif
#if defined(FLT_EVAL_METHOD)
  fprintf(stderr, "FLT_EVAL_METHOD     = %d\n", FLT_EVAL_METHOD);
#endif
#if defined(FP_FAST_FMA)
  fprintf(stderr, "FP_FAST_FMA         = defined\n");
#endif
#if defined(__FAST_MATH__)
  fprintf(stderr, "__FAST_MATH__       = defined\n");
#endif
  fprintf(stderr, "\n");

#if defined(_SC_VERSION)
  fprintf(stderr, "POSIX.1 VERSION     = %ld\n", sysconf(_SC_VERSION));
#endif
#if defined(_SC_ARG_MAX)
  fprintf(stderr, "POSIX.1 ARG_MAX     = %ld\n", sysconf(_SC_ARG_MAX));
#endif
#if defined(_SC_CHILD_MAX)
  fprintf(stderr, "POSIX.1 CHILD_MAX   = %ld\n", sysconf(_SC_CHILD_MAX));
#endif
#if defined(_SC_STREAM_MAX)
  fprintf(stderr, "POSIX.1 STREAM_MAX  = %ld\n", sysconf(_SC_STREAM_MAX));
#endif
#if defined(_SC_OPEN_MAX)
  fprintf(stderr, "POSIX.1 OPEN_MAX    = %ld\n", sysconf(_SC_OPEN_MAX));
#endif
#if defined(_SC_PAGESIZE)
  fprintf(stderr, "POSIX.1 PAGESIZE    = %ld\n", sysconf(_SC_PAGESIZE));
#endif

  fprintf(stderr, "\n");

#if defined(HAVE_GETRLIMIT)
#if defined(RLIMIT_FSIZE)
  PRINT_RLIMIT(RLIMIT_FSIZE);
#endif
#if defined(RLIMIT_NOFILE)
  PRINT_RLIMIT(RLIMIT_NOFILE);
#endif
#if defined(RLIMIT_STACK)
  PRINT_RLIMIT(RLIMIT_STACK);
#endif
#endif
  fprintf(stderr, "\n");
}

static void
check_stacksize()
{
#if defined(HAVE_GETRLIMIT)
#if defined(RLIMIT_STACK)
  {
    struct rlimit rlim;
    int status = getrlimit(RLIMIT_STACK, &rlim);
    if (status == 0)
      {
#define MIN_STACK_SIZE 67108864L /* 64MB */
        RLIM_T min_stack_size = MIN_STACK_SIZE;
        if (min_stack_size > rlim.rlim_max) min_stack_size = rlim.rlim_max;
        if (rlim.rlim_cur < min_stack_size)
          {
            rlim.rlim_cur = min_stack_size;

            status = setrlimit(RLIMIT_STACK, &rlim);
            if (Debug)
              {
                if (status == 0)
                  {
                    fprintf(stderr, "Set stack size to %ld\n", (long) min_stack_size);
                    PRINT_RLIMIT(RLIMIT_STACK);
                  }
                else
                  fprintf(stderr, "Set stack size to %ld failed!\n", (long) min_stack_size);
                fprintf(stderr, "\n");
              }
          }
      }
  }
#endif
#endif
}

static void
cdo_set_options(void)
{
  if (Debug)
    {
      fprintf(stderr, "CDO_CMOR_Mode       = %d\n", CDO_CMOR_Mode);
      fprintf(stderr, "CDO_netcdf_hdr_pad  = %d\n", CDO_netcdf_hdr_pad);
      fprintf(stderr, "\n");
    }

  if (CDO_CMOR_Mode) cdiDefGlobal("CMOR_MODE", CDO_CMOR_Mode);
  if (CDO_Reduce_Dim) cdiDefGlobal("REDUCE_DIM", CDO_Reduce_Dim);
  if (CDO_netcdf_hdr_pad > 0) cdiDefGlobal("NETCDF_HDR_PAD", CDO_netcdf_hdr_pad);
}

static long
str_to_int(const char *intstring)
{
  long intval = -1;

  if (intstring)
    {
      long fact = 1;
      int len = (int) strlen(intstring);
      for (int loop = 0; loop < len; loop++)
        {
          if (!isdigit((int) intstring[loop]))
            {
              switch (tolower((int) intstring[loop]))
                {
                case 'k': fact = 1024; break;
                case 'm': fact = 1048576; break;
                case 'g': fact = 1073741824; break;
                default: fact = 0; break;
                }
              break;
            }
        }

      if (fact) intval = fact * atol(intstring);
    }

  return intval;
}

static int
parse_options_long(int argc, char *argv[])
{
  int c;
  int lnetcdf_hdr_pad;
  int luse_fftw;
  int lgridsearchnn;
  int lgridsearchradius;
  int lremap_genweights;
  int lprecision;
  int lpercentile;
  int lprintoperatorsno = 0;
  int lprintoperators = 0;
  int lenableexcept;
  int ltimestat_date;
  int ltimestat_bounds;
  int lsortname;
  int lsortparam;
  int ldebLevel;
  int lscmode;

  // clang-format off
  struct cdo_option opt_long[] =
    {
      { "precision",         required_argument,        &lprecision,   1  },
      { "percentile",        required_argument,        &lpercentile,  1  },
      { "netcdf_hdr_pad",    required_argument,    &lnetcdf_hdr_pad,  1  },
      { "header_pad",        required_argument,    &lnetcdf_hdr_pad,  1  },
      { "hdr_pad",           required_argument,    &lnetcdf_hdr_pad,  1  },
      { "use_fftw",          required_argument,          &luse_fftw,  1  },
      { "gridsearchnn",      required_argument,      &lgridsearchnn,  1  },
      { "gridsearchradius",  required_argument,  &lgridsearchradius,  1  },
      { "remap_genweights",  required_argument,  &lremap_genweights,  1  },
      { "enableexcept",      required_argument,      &lenableexcept,  1  },
      { "timestat_date",     required_argument,     &ltimestat_date,  1  },
      { "timestat_bounds",         no_argument,   &ltimestat_bounds,  1  },
      { "cmor",                    no_argument,      &CDO_CMOR_Mode,  1  },
      { "reduce_dim",              no_argument,     &CDO_Reduce_Dim,  1  },
      { "float",                   no_argument,        &CDO_Memtype,  MEMTYPE_FLOAT  },
      { "rusage",                  no_argument,         &CDO_Rusage,  1  },
      { "operators_no_output",     no_argument,  &lprintoperatorsno,  1  },
      { "operators",               no_argument,    &lprintoperators,  1  },
      { "no_warnings",             no_argument,           &_Verbose,  0  },
      { "color",                   no_argument,                NULL, 'C' },
      { "format",            required_argument,                NULL, 'f' },
      { "help",                    no_argument,                NULL, 'h' },
      { "history",                 no_argument, &CDO_Append_History,  0  },
      { "no_history",              no_argument, &CDO_Append_History,  0  },
      { "regular",                 no_argument,                NULL, 'R' },
      { "silent",                  no_argument,                NULL, 's' },
      { "sort",                    no_argument,                NULL, 'Q' },
      { "sortname",                no_argument,          &lsortname,  1  },
      { "sortparam",               no_argument,         &lsortparam,  1  },
      { "table",             required_argument,                NULL, 't' },
      { "verbose",                 no_argument,                NULL, 'v' },
      { "version",                 no_argument,                NULL, 'V' },
      { "Dkext",             required_argument,          &ldebLevel,  1  },
      { "outputGribDataScanningMode", required_argument,  &lscmode,   1  },
      { "seperateDebugFromLog", required_argument,             NULL,  2  },
      { NULL,                                0,                NULL,  0  }
    };
  // clang-format on

  CDO_opterr = 1;

  while (1)
    {
      // IMPORTANT: BY EVERY OPTION that takes arguments you MUST set its
      // trigger variable to ZERO; otherwise the parameters of other options get
      // wrongly assigned.
      lprecision = 0;
      lpercentile = 0;
      lnetcdf_hdr_pad = 0;
      luse_fftw = 0;
      lgridsearchnn = 0;
      lgridsearchradius = 0;
      lremap_genweights = 0;
      lenableexcept = 0;
      ltimestat_date = 0;
      ltimestat_bounds = 0;
      lsortname = 0;
      lsortparam = 0;
      ldebLevel = 0;
      lscmode = 0;

      c = cdo_getopt_long(argc, argv, "f:b:e:P:g:i:k:l:m:n:t:D:z:aBCcdhLMOpQRrsSTuVvWXZ", opt_long, NULL);
      if (c == -1) break;

      switch (c)
        {
        case '?':
          // cdo_usage();
          // fprintf(stderr, "Illegal option!\n");
          return -1;
        // break;
        case ':':
          // cdo_usage();
          // fprintf(stderr, "Option requires an argument!\n");
          return -1;
        // break;
        case 0:
          if (lnetcdf_hdr_pad)
            {
              int netcdf_hdr_pad = str_to_int(CDO_optarg);
              if (netcdf_hdr_pad >= 0) CDO_netcdf_hdr_pad = netcdf_hdr_pad;
            }
          else if (lprecision)
            {
              cdo_set_digits(CDO_optarg);
            }
          else if (lpercentile)
            {
              percentile_set_method(CDO_optarg);
            }
          else if (lenableexcept)
            {
              int except = -1;
              if (STR_IS_EQ(CDO_optarg, "DIVBYZERO"))
                except = FE_DIVBYZERO;
              else if (STR_IS_EQ(CDO_optarg, "INEXACT"))
                except = FE_INEXACT;
              else if (STR_IS_EQ(CDO_optarg, "INVALID"))
                except = FE_INVALID;
              else if (STR_IS_EQ(CDO_optarg, "OVERFLOW"))
                except = FE_OVERFLOW;
              else if (STR_IS_EQ(CDO_optarg, "UNDERFLOW"))
                except = FE_UNDERFLOW;
              else if (STR_IS_EQ(CDO_optarg, "ALL_EXCEPT"))
                except = FE_ALL_EXCEPT;
              if (except < 0) cdoAbort("option --%s: unsupported argument: %s", "enableexcept", CDO_optarg);
              cdo_feenableexcept(except);
              if (signal(SIGFPE, cdo_sig_handler) == SIG_ERR) cdoWarning("can't catch SIGFPE!");
            }
          else if (ltimestat_date)
            {
              int timestatdate = -1;
              if (STR_IS_EQ(CDO_optarg, "first"))
                timestatdate = TIMESTAT_FIRST;
              else if (STR_IS_EQ(CDO_optarg, "last"))
                timestatdate = TIMESTAT_LAST;
              else if (STR_IS_EQ(CDO_optarg, "middle"))
                timestatdate = TIMESTAT_MEAN;
              else if (STR_IS_EQ(CDO_optarg, "midhigh"))
                timestatdate = TIMESTAT_MIDHIGH;
              if (timestatdate < 0) cdoAbort("option --%s: unsupported argument: %s", "timestat_date", CDO_optarg);
              extern int CDO_Timestat_Date;
              CDO_Timestat_Date = timestatdate;
            }
          else if (ltimestat_bounds)
            {
              extern bool CDO_Timestat_Bounds;
              CDO_Timestat_Bounds = true;
            }
          else if (luse_fftw)
            {
              int intarg = parameter2int(CDO_optarg);
              if (intarg != 0 && intarg != 1)
                cdoAbort("Unsupported value for option --use_fftw=%d [range: 0-1]", intarg);
              CDO_Use_FFTW = intarg;
            }
          else if (lgridsearchnn)
            {
              gridsearch_set_method(CDO_optarg);
            }
          else if (lgridsearchradius)
            {
              extern double gridsearch_radius;
              double fval = radius_str_to_deg(CDO_optarg);
              if (fval < 0 || fval > 180) cdoAbort("%s=%g out of bounds (0-180 deg)!", "gridsearchradius", fval);
              gridsearch_radius = fval;
            }
          else if (lremap_genweights)
            {
              int intarg = parameter2int(CDO_optarg);
              if (intarg != 0 && intarg != 1)
                cdoAbort("Unsupported value for option --remap_genweights %d [0/1]", intarg);
              REMAP_genweights = intarg;
            }
          else if (lsortname)
            {
              cdiDefGlobal("SORTNAME", TRUE);
            }
          else if (lsortparam)
            {
              cdiDefGlobal("SORTPARAM", TRUE);
            }
#ifdef HIRLAM_EXTENSIONS
          else if (ldebLevel)
            {
              int newDebLevelVal = parameter2int(CDO_optarg);
              if (newDebLevelVal > 0)
                {
                  extern int cdiDebugExt;
                  CdoDebug::cdoDebugExt = newDebLevelVal;
                  cdiDebugExt = newDebLevelVal;
                }
            }
          else if (lscmode)
            {
              int scanningModeValue = atoi(CDO_optarg);
              if (CdoDebug::cdoDebugExt) printf("scanningModeValue=%d\n", scanningModeValue);

              if ((scanningModeValue == 0) || (scanningModeValue == 64) || (scanningModeValue == 96))
                {
                  streamGrbDefDataScanningMode(scanningModeValue);  // -1: not used; allowed modes: <0,
                                                                    // 64, 96>; Default is 64
                }
              else
                {
                  cdoAbort("Warning: %d not in allowed modes: <0, 64, 96>; "
                           "Using default: 64\n",
                           scanningModeValue);
                  streamGrbDefDataScanningMode(64);
                }
            }
#endif
          break;
        case 'a': cdoDefaultTimeType = TAXIS_ABSOLUTE; break;
        case 'b': setDefaultDataType(CDO_optarg); break;
        case 'B': cdoBenchmark = TRUE; break;
        case 'C': CDO_Color = TRUE; break;
        case 'c': cdoCheckDatarange = TRUE; break;
        case 'd': Debug = 1; break;
        case 'D':
          Debug = 1;
          DebugLevel = atoi(CDO_optarg);
          break;
        case 'e': {
#if defined(HAVE_GETHOSTNAME)
            char host[1024];
            gethostname(host, sizeof(host));
            cdoExpName = CDO_optarg;
            /* printf("host: %s %s\n", host, cdoExpName); */
            cdoExpMode = STR_IS_EQ(host, cdoExpName) ? CDO_EXP_REMOTE : CDO_EXP_LOCAL;
#else
            fprintf(stderr, "Function gethostname not available!\n");
            exit(EXIT_FAILURE);
#endif
            break;
          }
        case 'f': setDefaultFileType(CDO_optarg, 1); break;
        case 'g': cdo_set_grids(CDO_optarg); break;
        case 'h': Help = 1; break;
        case 'i': defineInstitution(CDO_optarg); break;
        case 'k': defineChunktype(CDO_optarg); break;
        case 'L': Threading::cdoLockIO = TRUE; break;
        case 'l': defineZaxis(CDO_optarg); break;
        case 'm': cdiDefMissval(atof(CDO_optarg)); break;
        case 'M': cdiDefGlobal("HAVE_MISSVAL", TRUE); break;
        case 'n': defineVarnames(CDO_optarg); break;
        case 'O': cdoOverwriteMode = TRUE; break;
        case 'P':
          if (*CDO_optarg < '1' || *CDO_optarg > '9')
            {
              fprintf(stderr, "Unexpected character in number of OpenMP threads (-P "
                              "<nthreads>): %s!\n",
                      CDO_optarg);
              exit(EXIT_FAILURE);
            }
          numThreads = atoi(CDO_optarg);
          break;
        case 'p':
          CDO_Parallel_Read = TRUE;
          CDO_task = true;
          break;
        case 'Q': cdiDefGlobal("SORTNAME", TRUE); break;
        case 'R':
          cdoRegulargrid = TRUE;
          cdiDefGlobal("REGULARGRID", TRUE);
          break;
        case 'r': cdoDefaultTimeType = TAXIS_RELATIVE; break;
        case 'S': cdoDiag = TRUE; break;
        case 's': Options::silentMode = TRUE; break;
        case 'T': cdoTimer = TRUE; break;
        case 't': cdoDefaultTableID = defineTable(CDO_optarg); break;
        case 'u': Options::cdoInteractive = TRUE; break;
        case 'V': Version = 1; break;
        case 'v':
          cdoVerbose = TRUE;
          _Verbose = 1;
          break;
        case 'W': /* Warning messages */ _Verbose = 1; break;
        case 'X': /* multi threaded I/O */ cdoParIO = TRUE; break;
        case 'Z': Options::cdoCompress = true; break;
        case 'z': defineCompress(CDO_optarg); break;
        case 2:
#ifdef DEBUG
          CdoDebug::outfile = CDO_optarg;
          CdoDebug::print_to_seperate_file = true;
#endif
          break;
        }
    }

  if (lprintoperators || lprintoperatorsno)
    {
      set_text_color(stderr, RESET, GREEN);
      bool print_no_output = lprintoperatorsno > 0;
      operatorPrintList(print_no_output);
      // operatorPrintAll();
      reset_text_color(stderr);
      return 1;
    }

  return 0;
}

static void
cdo_rusage(void)
{
#if defined(HAVE_SYS_RESOURCE_H) && defined(RUSAGE_SELF)
  struct rusage ru;
  int status = getrusage(RUSAGE_SELF, &ru);

  if (status == 0)
    {
      double ut = ru.ru_utime.tv_sec + 0.000001 * ru.ru_utime.tv_usec;
      double st = ru.ru_stime.tv_sec + 0.000001 * ru.ru_stime.tv_usec;

      fprintf(stderr, "  User time:     %.3f seconds\n", ut);
      fprintf(stderr, "  System time:   %.3f seconds\n", st);
      fprintf(stderr, "  Total time:    %.3f seconds\n", ut + st);
      fprintf(stderr, "  Memory usage:  %.2f MBytes\n", ru.ru_maxrss / (1024. * 1024.));
      fprintf(stderr, "  Page reclaims: %5ld page%s\n", ru.ru_minflt, ADD_PLURAL(ru.ru_minflt));
      fprintf(stderr, "  Page faults:   %5ld page%s\n", ru.ru_majflt, ADD_PLURAL(ru.ru_majflt));
      fprintf(stderr, "  Swaps:         %5ld\n", ru.ru_nswap);
      fprintf(stderr, "  Disk read:     %5ld block%s\n", ru.ru_inblock, ADD_PLURAL(ru.ru_inblock));
      fprintf(stderr, "  Disk Write:    %5ld block%s\n", ru.ru_oublock, ADD_PLURAL(ru.ru_oublock));
    }
#endif
}

/* clang-format off */
// stream in  -1 means: unlimited number of input streams
// stream out -1 means: usage of obase
/***
 * Initializes all hardcoded modules.
 */
void init_modules()
{
/*                             function        help function      operator names          mode number     num streams
                                                                                                  type       in out      */
  add_module("Adisit"        , {Adisit        , AdisitHelp        , AdisitOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Afterburner"   , {Afterburner   , AfterburnerHelp   , AfterburnerOperators   , 1 , CDI_REAL , -1 , 1  });
  add_module("Arith"         , {Arith         , ArithHelp         , ArithOperators         , 1 , CDI_REAL , 2  , 1  });
  add_module("Arithc"        , {Arithc        , ArithcHelp        , ArithcOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Arithdays"     , {Arithdays     , ArithdaysHelp     , ArithdaysOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Arithlat"      , {Arithlat      , {}                , ArithlatOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Cat"           , {Cat           , CopyHelp          , CatOperators           , 1 , CDI_REAL , -1 , 1  });
  add_module("CDItest"       , {CDItest       , {}                , CDItestOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("CDIread"       , {CDIread       , {}                , CDIreadOperators       , 1 , CDI_REAL , 1  , 0  });
  add_module("CDIwrite"      , {CDIwrite      , {}                , CDIwriteOperators      , 1 , CDI_REAL , 0  , 1  });
  add_module("Change"        , {Change        , ChangeHelp        , ChangeOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Change_e5slm"  , {Change_e5slm  , {}                , Change_e5slmOperators  , 0 , CDI_REAL , 1  , 1  });
  add_module("Cloudlayer"    , {Cloudlayer    , {}                , CloudlayerOperators    , 1 , CDI_REAL , 1  , 1  });
  add_module("CMOR"          , {CMOR          , CMORHelp          , CMOROperators          , 1 , CDI_REAL , 1  , 0  });
  add_module("CMOR_lite"     , {CMOR_lite     , CMORliteHelp      , CMORliteOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("CMOR_table"    , {CMOR_table    , {}                , CMORtableOperators     , 1 , CDI_REAL , 0  , 0  });
  add_module("Collgrid"      , {Collgrid      , CollgridHelp      , CollgridOperators      , 1 , CDI_REAL , -1 , 1  });
  add_module("Command"       , {Command       , {}                , CommandOperators       , 0 , CDI_REAL , 1  , 0  });
  add_module("Comp"          , {Comp          , CompHelp          , CompOperators          , 1 , CDI_REAL , 2  , 1  });
  add_module("Compc"         , {Compc         , CompcHelp         , CompcOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("Complextorect" , {Complextorect , {}                , ComplextorectOperators , 1 , CDI_COMP , 1  , 2  });
  add_module("Cond"          , {Cond          , CondHelp          , CondOperators          , 1 , CDI_REAL , 2  , 1  });
  add_module("Cond2"         , {Cond2         , Cond2Help         , Cond2Operators         , 1 , CDI_REAL , 3  , 1  });
  add_module("Condc"         , {Condc         , CondcHelp         , CondcOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("Consecstat"    , {Consecstat    , ConsecstatHelp    , ConsecstatOperators    , 1 , CDI_REAL , 1  , 1  });
  add_module("Copy"          , {Copy          , CopyHelp          , CopyOperators          , 1 , CDI_REAL , -1 , 1  });
  add_module("Deltat"        , {Deltat        , {}                , DeltatOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Deltime"       , {Deltime       , {}                , DeltimeOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Derivepar"     , {Derivepar     , DeriveparHelp     , DeriveparOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Detrend"       , {Detrend       , DetrendHelp       , DetrendOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Diff"          , {Diff          , DiffHelp          , DiffOperators          , 1 , CDI_REAL , 2  , 0  });
  add_module("Distgrid"      , {Distgrid      , DistgridHelp      , DistgridOperators      , 1 , CDI_REAL , 1  , -1  });
  add_module("Duplicate"     , {Duplicate     , DuplicateHelp     , DuplicateOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Echam5ini"     , {Echam5ini     , {}                , Echam5iniOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Enlarge"       , {Enlarge       , EnlargeHelp       , EnlargeOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Enlargegrid"   , {Enlargegrid   , {}                , EnlargegridOperators   , 0 , CDI_REAL , 1  , 1  });
  add_module("Ensstat"       , {Ensstat       , EnsstatHelp       , EnsstatOperators       , 1 , CDI_REAL , -1 , 1  });
  add_module("Ensstat3"      , {Ensstat3      , Ensstat2Help      , Ensstat3Operators      , 1 , CDI_REAL , -1 , 1  });
  add_module("Ensval"        , {Ensval        , EnsvalHelp        , EnsvalOperators        , 1 , CDI_REAL , -1 , 1  });
  add_module("Eofcoeff"      , {Eofcoeff      , EofcoeffHelp      , EofcoeffOperators      , 1 , CDI_REAL , 2  , -1 });
  add_module("Eofcoeff3d"    , {Eofcoeff3d    , EofcoeffHelp      , Eofcoeff3dOperators    , 1 , CDI_REAL , 2  , -1 });
  add_module("EOFs"          , {EOFs          , EOFsHelp          , EOFsOperators          , 1 , CDI_REAL , 1  , 2  });
  add_module("EOF3d"         , {EOF3d         , EOFsHelp          , EOF3dOperators         , 1 , CDI_REAL , 1  , 2  });
  add_module("EstFreq"       , {EstFreq       , {}                , EstFreqOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Expr"          , {Expr          , ExprHelp          , ExprOperators          , 1 , CDI_REAL , 1  , 1  });
  add_module("FC"            , {FC            , {}                , FCOperators            , 1 , CDI_REAL , 1  , 1  });
  add_module("Filedes"       , {Filedes       , FiledesHelp       , FiledesOperators       , 1 , CDI_BOTH , 1  , 0  });
  add_module("Fillmiss"      , {Fillmiss      , {}                , FillmissOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Filter"        , {Filter        , FilterHelp        , FilterOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Fldrms"        , {Fldrms        , {}                , FldrmsOperators        , 1 , CDI_REAL , 2  , 1  });
  add_module("Fldstat"       , {Fldstat       , FldstatHelp       , FldstatOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Fldstatcor"    , {Fldstat2      , FldcorHelp        , FldcorOperators        , 1 , CDI_REAL , 2  , 1  });
  add_module("Fldstatvar"    , {Fldstat2      , FldcovarHelp      , FldcovarOperators      , 1 , CDI_REAL , 2  , 1  });
  add_module("Fourier"       , {Fourier       , {}                , FourierOperators       , 1 , CDI_COMP , 1  , 1  });
  add_module("Gengrid"       , {Gengrid       , {}                , GengridOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("Gradsdes"      , {Gradsdes      , GradsdesHelp      , GradsdesOperators      , 1 , CDI_REAL , 1  , 0  });
  add_module("Gridboxstat"   , {Gridboxstat   , GridboxstatHelp   , GridboxstatOperators   , 1 , CDI_REAL , 1  , 1  });
  add_module("Gridcell"      , {Gridcell      , GridcellHelp      , GridcellOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Gridsearch"    , {Gridsearch    , {}                , GridsearchOperators    , 0 , CDI_REAL , 0  , 0  });
  add_module("Harmonic"      , {Harmonic      , {}                , HarmonicOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Histogram"     , {Histogram     , HistogramHelp     , HistogramOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Importamsr"    , {Importamsr    , ImportamsrHelp    , ImportamsrOperators    , 1 , CDI_REAL , 1  , 1  });
  add_module("Importbinary"  , {Importbinary  , ImportbinaryHelp  , ImportbinaryOperators  , 1 , CDI_REAL , 1  , 1  });
  add_module("Importcmsaf"   , {Importcmsaf   , ImportcmsafHelp   , ImportcmsafOperators   , 1 , CDI_REAL , 1  , 1  });
  add_module("Importobs"     , {Importobs     , {}                , ImportobsOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Info"          , {Info          , InfoHelp          , InfoOperators          , 1 , CDI_BOTH , -1 , 0  });
  add_module("Input"         , {Input         , InputHelp         , InputOperators         , 1 , CDI_REAL , 0  , 1  });
  add_module("Intgrid"       , {Intgrid       , {}                , IntgridOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Intgridtraj"   , {Intgridtraj   , {}                , IntgridtrajOperators   , 1 , CDI_REAL , 1  , 1  });
  add_module("Intlevel"      , {Intlevel      , IntlevelHelp      , IntlevelOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Intlevel3d"    , {Intlevel3d    , Intlevel3dHelp    , Intlevel3dOperators    , 1 , CDI_REAL , 2  , 1  });
  add_module("Inttime"       , {Inttime       , InttimeHelp       , InttimeOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Intntime"      , {Intntime      , InttimeHelp       , IntntimeOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Intyear"       , {Intyear       , IntyearHelp       , IntyearOperators       , 1 , CDI_REAL , 2  , -1 });
  add_module("Invert"        , {Invert        , InvertHelp        , InvertOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Invertlev"     , {Invertlev     , InvertlevHelp     , InvertlevOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Isosurface"    , {Isosurface    , {}                , IsosurfaceOperators    , 1 , CDI_REAL , 1  , 1  });
  add_module("MapReduce"     , {MapReduce     , MapReduceHelp     , MapReduceOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Maskbox"       , {Maskbox       , MaskboxHelp       , MaskboxOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Maskregion"    , {Maskbox       , MaskregionHelp    , MaskregionOperators    , 1 , CDI_REAL , 1  , 1  });
  add_module("Mastrfu"       , {Mastrfu       , MastrfuHelp       , MastrfuOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Math"          , {Math          , MathHelp          , MathOperators          , 1 , CDI_BOTH , 1  , 1  });
  add_module("Merge"         , {Merge         , MergeHelp         , MergeOperators         , 1 , CDI_REAL , -1 , 1  });
  add_module("Mergetime"     , {Mergetime     , MergeHelp         , MergetimeOperators     , 1 , CDI_REAL , -1 , 1  });
  add_module("Mergegrid"     , {Mergegrid     , MergegridHelp     , MergegridOperators     , 1 , CDI_REAL , 2  , 1  });
  add_module("Merstat"       , {Merstat       , MerstatHelp       , MerstatOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Monarith"      , {Monarith      , MonarithHelp      , MonarithOperators      , 1 , CDI_REAL , 2  , 1  });
  add_module("Mrotuv"        , {Mrotuv        , {}                , MrotuvOperators        , 1 , CDI_REAL , 1  , 2  });
  add_module("Mrotuvb"       , {Mrotuvb       , {}                , MrotuvbOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("NCL_wind"      , {NCL_wind      , NCL_windHelp      , NCL_windOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Ninfo"         , {Ninfo         , NinfoHelp         , NinfoOperators         , 1 , CDI_BOTH , 1  , 0  });
  add_module("Nmldump"       , {Nmldump       , {}                , NmldumpOperators       , 0 , CDI_REAL , 0  , 0  });
  add_module("Output"        , {Output        , OutputHelp        , OutputOperators        , 1 , CDI_REAL , -1 , 0  });
  add_module("Outputtab"     , {Output        , OutputtabHelp     , OutputtabOperators     , 1 , CDI_REAL , -1 , 0  });
  add_module("Outputgmt"     , {Outputgmt     , OutputgmtHelp     , OutputgmtOperators     , 1 , CDI_REAL , 1  , 0  });
  add_module("Pack"          , {Pack          , {}                , PackOperators          , 1 , CDI_REAL , 1  , 1  });
  add_module("Pardup"        , {Pardup        , {}                , PardupOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Pinfo"         , {Pinfo         , {}                , PinfoOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("Pressure"      , {Pressure      , {}                , PressureOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Regres"        , {Regres        , RegresHelp        , RegresOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Remap"         , {Remap         , RemapHelp         , RemapOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("Remapbil"      , {Remap         , RemapbilHelp      , RemapbilOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Remapbic"      , {Remap         , RemapbicHelp      , RemapbicOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Remapnn"       , {Remap         , RemapnnHelp       , RemapnnOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Remapdis"      , {Remap         , RemapdisHelp      , RemapdisOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Remapycon"     , {Remap         , RemapyconHelp     , RemapyconOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Remapcon"      , {Remap         , RemapconHelp      , RemapconOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Remapcon2"     , {Remap         , Remapcon2Help     , Remapcon2Operators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Remaplaf"      , {Remap         , RemaplafHelp      , RemaplafOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Remapgrid"     , {Remap         , {}                , RemapgridOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Remapeta"      , {Remapeta      , RemapetaHelp      , RemapetaOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Replace"       , {Replace       , ReplaceHelp       , ReplaceOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("Replacevalues" , {Replacevalues , ReplacevaluesHelp , ReplacevaluesOperators , 1 , CDI_REAL , 1  , 1  });
  add_module("Rhopot"        , {Rhopot        , RhopotHelp        , RhopotOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Rotuv"         , {Rotuv         , RotuvbHelp        , RotuvOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("Runpctl"       , {Runpctl       , RunpctlHelp       , RunpctlOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Runstat"       , {Runstat       , RunstatHelp       , RunstatOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Samplegridicon", {Samplegridicon, {}                , SamplegridiconOperators, 1,  CDI_REAL,  1  , 2  });
  add_module("Seascount"     , {Seascount     , {}                , SeascountOperators     , 1 , CDI_BOTH , 1  , 1  });
  add_module("Seaspctl"      , {Seaspctl      , SeaspctlHelp      , SeaspctlOperators      , 1 , CDI_REAL , 3  , 1  });
  add_module("Seasstat"      , {Seasstat      , SeasstatHelp      , SeasstatOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Selbox"        , {Selbox        , SelboxHelp        , SelboxOperators        , 1 , CDI_BOTH , 1  , 1  });
  add_module("Selgridcell"   , {Selgridcell   , {}                , SelgridcellOperators   , 1 , CDI_BOTH , 1  , 1  });
  add_module("Select"        , {Select        , SelectHelp        , SelectOperators        , 1 , CDI_BOTH , -1 , 1  });
  add_module("Selvar"        , {Selvar        , SelvarHelp        , SelvarOperators        , 1 , CDI_BOTH , 1  , 1  });
  add_module("Selrec"        , {Selrec        , SelvarHelp        , SelrecOperators        , 1 , CDI_BOTH , 1  , 1  });
  add_module("Seloperator"   , {Seloperator   , {}                , SeloperatorOperators   , 1 , CDI_REAL , 1  , 1  });
  add_module("Seltime"       , {Seltime       , SeltimeHelp       , SeltimeOperators       , 1 , CDI_BOTH , 1  , 1  });
  add_module("Set"           , {Set           , SetHelp           , SetOperators           , 1 , CDI_BOTH , 1  , 1  });
  add_module("Setattribute"  , {Setattribute  , SetattributeHelp  , SetattributeOperators  , 1 , CDI_REAL , 1  , 1  });
  add_module("Setbox"        , {Setbox        , SetboxHelp        , SetboxOperators        , 1 , CDI_REAL , 1  , 1  });
  //add_module("Setgatt"       , {Setgatt       , SetgattHelp       , SetgattOperators       , 1 , CDI_BOTH , 1  , 1  });  
  add_module("Setgrid"       , {Setgrid       , SetgridHelp       , SetgridOperators       , 1 , CDI_BOTH , 1  , 1  });
  add_module("Sethalo"       , {Sethalo       , SethaloHelp       , SethaloOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Setmiss"       , {Setmiss       , SetmissHelp       , SetmissOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Setmisstonn"   , {Fillmiss      , SetmissHelp       , SetmisstonnOperators   , 1 , CDI_REAL , 1  , 1  });
  add_module("Setcodetab"    , {Setpartab     , SetHelp           , SetcodetabOperators    , 1 , CDI_REAL , 1  , 1  });
  add_module("Setpartab"     , {Setpartab     , SetpartabHelp     , SetpartabOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Setrcaname"    , {Setrcaname    , {}                , SetrcanameOperators    , 1 , CDI_REAL , 1  , 1  });
  add_module("Settime"       , {Settime       , SettimeHelp       , SettimeOperators       , 1 , CDI_BOTH , 1  , 1  });
  add_module("Setzaxis"      , {Setzaxis      , SetzaxisHelp      , SetzaxisOperators      , 1 , CDI_BOTH , 1  , 1  });
  add_module("Shiftxy"       , {Shiftxy       , {}                , ShiftxyOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Showinfo"      , {Showinfo      , ShowinfoHelp      , ShowinfoOperators      , 1 , CDI_BOTH , 1  , 0  });
  add_module("Showattribute" , {Showattribute , ShowattributeHelp , ShowattributeOperators , 1 , CDI_REAL , 1  , 0  });
  add_module("Sinfo"         , {Sinfo         , SinfoHelp         , SinfoOperators         , 1 , CDI_BOTH , -1 , 0  });
  add_module("Smooth"        , {Smooth        , SmoothHelp        , SmoothOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Sort"          , {Sort          , {}                , SortOperators          , 1 , CDI_REAL , 1  , 1  });
  add_module("Sorttimestamp" , {Sorttimestamp , {}                , SorttimestampOperators , 1 , CDI_REAL , -1 , 1  });
  add_module("Specinfo"      , {Specinfo      , {}                , SpecinfoOperators      , 1 , CDI_REAL , 0  , 0  });
  add_module("Spectral"      , {Spectral      , SpectralHelp      , SpectralOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Spectrum"      , {Spectrum      , {}                , SpectrumOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Split"         , {Split         , SplitHelp         , SplitOperators         , 1 , CDI_BOTH , 1  , -1 });
  add_module("Splitrec"      , {Splitrec      , SplitHelp         , SplitrecOperators      , 1 , CDI_BOTH , 1  , -1 });
  add_module("Splitsel"      , {Splitsel      , SplitselHelp      , SplitselOperators      , 1 , CDI_BOTH , 1  , -1 });
  add_module("Splittime"     , {Splittime     , SplittimeHelp     , SplittimeOperators     , 1 , CDI_BOTH , 1  , -1 });
  add_module("Splityear"     , {Splityear     , SplittimeHelp     , SplityearOperators     , 1 , CDI_BOTH , 1  , -1 });
  add_module("Subtrend"      , {Subtrend      , SubtrendHelp      , SubtrendOperators      , 1 , CDI_REAL , 3  , 1  });
  add_module("Tee"           , {Tee           , TeeHelp           , TeeOperators           , 1 , CDI_REAL , 2  , 1  });
  add_module("Template1"     , {Template1     , {}                , Template1Operators     , 0 , CDI_REAL , 1  , 1  });
  add_module("Template2"     , {Template2     , {}                , Template2Operators     , 0 , CDI_REAL , 1  , 1  });
  add_module("Test"          , {Test          , {}                , TestOperators          , 0 , CDI_REAL , 1  , 1  });
  add_module("Test2"         , {Test2         , {}                , Test2Operators         , 0 , CDI_REAL , 2  , 1  });
  add_module("Testdata"      , {Testdata      , {}                , TestdataOperators      , 0 , CDI_REAL , 1  , 1  });
  add_module("Tests"         , {Tests         , {}                , TestsOperators         , 0 , CDI_REAL , 1  , 1  });
  add_module("Timcount"      , {Timcount      , {}                , TimcountOperators      , 1 , CDI_BOTH , 1  , 1  });
  add_module("Yearcount"     , {Timcount      , {}                , YearcountOperators     , 1 , CDI_BOTH , 1  , 1  });
  add_module("Moncount"      , {Timcount      , {}                , MoncountOperators      , 1 , CDI_BOTH , 1  , 1  });
  add_module("Daycount"      , {Timcount      , {}                , DaycountOperators      , 1 , CDI_BOTH , 1  , 1  });
  add_module("Hourcount"     , {Timcount      , {}                , HourcountOperators     , 1 , CDI_BOTH , 1  , 1  });
  add_module("Timcumsum"     , {Timcumsum     , TimcumsumHelp     , TimcumsumOperators     , 1 , CDI_BOTH , 1  , 1  });
  add_module("Timpctl"       , {Timpctl       , TimpctlHelp       , TimpctlOperators       , 1 , CDI_REAL , 3  , 1  });
  add_module("Yearpctl"      , {Timpctl       , YearpctlHelp      , YearpctlOperators      , 1 , CDI_REAL , 3  , 1  });
  add_module("Monpctl"       , {Timpctl       , MonpctlHelp       , MonpctlOperators       , 1 , CDI_REAL , 3  , 1  });
  add_module("Daypctl"       , {Timpctl       , DaypctlHelp       , DaypctlOperators       , 1 , CDI_REAL , 3  , 1  });
  add_module("Hourpctl"      , {Timpctl       , HourpctlHelp      , HourpctlOperators      , 1 , CDI_REAL , 3  , 1  });
  add_module("Timselpctl"    , {Timselpctl    , TimselpctlHelp    , TimselpctlOperators    , 1 , CDI_REAL , 3  , 1  });
  add_module("Timsort"       , {Timsort       , TimsortHelp       , TimsortOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Timselstat"    , {Timselstat    , TimselstatHelp    , TimselstatOperators    , 1 , CDI_REAL , 1  , 1  });
  add_module("XTimstat"      , {XTimstat      , {}                , XTimstatOperators      , 0 , CDI_BOTH , 1  , 1  });
  add_module("Timstat"       , {Timstat       , TimstatHelp       , TimstatOperators       , 1 , CDI_BOTH , 1  , 1  });
  add_module("Yearstat"      , {Timstat       , YearstatHelp      , YearstatOperators      , 1 , CDI_BOTH , 1  , 1  });
  add_module("Monstat"       , {Timstat       , MonstatHelp       , MonstatOperators       , 1 , CDI_BOTH , 1  , 1  });
  add_module("Daystat"       , {Timstat       , DaystatHelp       , DaystatOperators       , 1 , CDI_BOTH , 1  , 1  });
  add_module("Hourstat"      , {Timstat       , HourstatHelp      , HourstatOperators      , 1 , CDI_BOTH , 1  , 1  });
  add_module("Timcor"        , {Timstat2      , TimcorHelp        , TimcorOperators        , 1 , CDI_REAL , 2  , 1  });
  add_module("Timscorvar"    , {Timstat2      , TimcovarHelp      , TimcovarOperators      , 1 , CDI_REAL , 2  , 1  });
  add_module("Timstat3"      , {Timstat3      , {}                , Timstat3Operators      , 1 , CDI_REAL , 2  , 1  });
  add_module("Tinfo"         , {Tinfo         , {}                , TinfoOperators         , 1 , CDI_BOTH , 1  , 0  });
  add_module("Tocomplex"     , {Tocomplex     , {}                , TocomplexOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Transpose"     , {Transpose     , {}                , TransposeOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Trend"         , {Trend         , TrendHelp         , TrendOperators         , 1 , CDI_REAL , 1  , 2  });
  add_module("Trms"          , {Trms          , {}                , TrmsOperators          , 0 , CDI_REAL , 2  , 1  });
  add_module("Tstepcount"    , {Tstepcount    , {}                , TstepcountOperators    , 1 , CDI_REAL , 1  , 1  });
  add_module("Vargen"        , {Vargen        , VargenHelp        , VargenOperators        , 1 , CDI_REAL , 0  , 1  });
  add_module("Varrms"        , {Varrms        , {}                , VarrmsOperators        , 0 , CDI_REAL , 2  , 1  });
  add_module("Vertintml"     , {Vertintml     , VertintmlHelp     , VertintmlOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Vertintap"     , {Vertintap     , VertintapHelp     , VertintapOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Vertstat"      , {Vertstat      , VertstatHelp      , VertstatOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Vertcum"       , {Vertcum       , {}                , VertcumOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Vertwind"      , {Vertwind      , {}                , VertwindOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Verifygrid"    , {Verifygrid    , {}                , VerifygridOperators    , 1 , CDI_REAL , 1  , 0  });
  add_module("Wind"          , {Wind          , WindHelp          , WindOperators          , 1 , CDI_REAL , 1  , 1  });
  add_module("Writegrid"     , {Writegrid     , {}                , WritegridOperators     , 1 , CDI_REAL , 1  , 1  }); // no cdi output
  add_module("Writerandom"   , {Writerandom   , {}                , WriterandomOperators   , 1 , CDI_REAL , 1  , 1  });
  add_module("YAR"           , {YAR           , {}                , YAROperators           , 0 , CDI_REAL , 1  , 1  });
  add_module("Yearmonstat"   , {Yearmonstat   , YearmonstatHelp   , YearmonstatOperators   , 1 , CDI_REAL , 1  , 1  });
  add_module("Ydayarith"     , {Ydayarith     , YdayarithHelp     , YdayarithOperators     , 1 , CDI_REAL , 2  , 1  });
  add_module("Ydaypctl"      , {Ydaypctl      , YdaypctlHelp      , YdaypctlOperators      , 1 , CDI_REAL , 3  , 1  });
  add_module("Ydaystat"      , {Ydaystat      , YdaystatHelp      , YdaystatOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Ydrunpctl"     , {Ydrunpctl     , YdrunpctlHelp     , YdrunpctlOperators     , 1 , CDI_REAL , 3  , 1  });
  add_module("Ydrunstat"     , {Ydrunstat     , YdrunstatHelp     , YdrunstatOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Yhourarith"    , {Yhourarith    , YhourarithHelp    , YhourarithOperators    , 1 , CDI_REAL , 2  , 1  });
  add_module("Yhourstat"     , {Yhourstat     , YhourstatHelp     , YhourstatOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Ymonarith"     , {Ymonarith     , YmonarithHelp     , YmonarithOperators     , 1 , CDI_REAL , 2  , 1  });
  add_module("Yseasarith"    , {Ymonarith     , YseasarithHelp    , YseasarithOperators    , 1 , CDI_REAL , 2  , 1  });
  add_module("Ymonpctl"      , {Ymonpctl      , YmonpctlHelp      , YmonpctlOperators      , 1 , CDI_REAL , 3  , 1  });
  add_module("Ymonstat"      , {Ymonstat      , YmonstatHelp      , YmonstatOperators      , 1 , CDI_REAL , 1  , 1  });
  add_module("Yseaspctl"     , {Yseaspctl     , YseaspctlHelp     , YseaspctlOperators     , 1 , CDI_REAL , 3  , 1  });
  add_module("Yseasstat"     , {Yseasstat     , YseasstatHelp     , YseasstatOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Zonstat"       , {Zonstat       , ZonstatHelp       , ZonstatOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaCfd"        , {EcaCfd        , EcaCfdHelp        , EcaCfdOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaCsu"        , {EcaCsu        , EcaCsuHelp        , EcaCsuOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaCwdi"       , {EcaCwdi       , EcaCwdiHelp       , EcaCwdiOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaCwfi"       , {EcaCwfi       , EcaCwfiHelp       , EcaCwfiOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaEtr"        , {EcaEtr        , EcaEtrHelp        , EcaEtrOperators        , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaFd"         , {EcaFd         , EcaFdHelp         , EcaFdOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaGsl"        , {EcaGsl        , EcaGslHelp        , EcaGslOperators        , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaHd"         , {EcaHd         , EcaHdHelp         , EcaHdOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaHwdi"       , {EcaHwdi       , EcaHwdiHelp       , EcaHwdiOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaHwfi"       , {EcaHwfi       , EcaHwfiHelp       , EcaHwfiOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaId"         , {EcaId         , EcaIdHelp         , EcaIdOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaSu"         , {EcaSu         , EcaSuHelp         , EcaSuOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaTr"         , {EcaTr         , EcaTrHelp         , EcaTrOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaTg10p"      , {EcaTg10p      , EcaTg10pHelp      , EcaTg10pOperators      , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaTg90p"      , {EcaTg90p      , EcaTg90pHelp      , EcaTg90pOperators      , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaTn10p"      , {EcaTn10p      , EcaTn10pHelp      , EcaTn10pOperators      , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaTn90p"      , {EcaTn90p      , EcaTn90pHelp      , EcaTn90pOperators      , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaTx10p"      , {EcaTx10p      , EcaTx10pHelp      , EcaTx10pOperators      , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaTx90p"      , {EcaTx90p      , EcaTx90pHelp      , EcaTx90pOperators      , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaCdd"        , {EcaCdd        , EcaCddHelp        , EcaCddOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaCwd"        , {EcaCwd        , EcaCwdHelp        , EcaCwdOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaRr1"        , {EcaRr1        , EcaRr1Help        , EcaRr1Operators        , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaPd"         , {EcaPd         , EcaPdHelp         , EcaPdOperators         , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaR75p"       , {EcaR75p       , EcaR75pHelp       , EcaR75pOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaR75ptot"    , {EcaR75ptot    , EcaR75ptotHelp    , EcaR75ptotOperators    , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaR90p"       , {EcaR90p       , EcaR90pHelp       , EcaR90pOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaR90ptot"    , {EcaR90ptot    , EcaR90ptotHelp    , EcaR90ptotOperators    , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaR95p"       , {EcaR95p       , EcaR95pHelp       , EcaR95pOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaR95ptot"    , {EcaR95ptot    , EcaR95ptotHelp    , EcaR95ptotOperators    , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaR99p"       , {EcaR99p       , EcaR99pHelp       , EcaR99pOperators       , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaR99ptot"    , {EcaR99ptot    , EcaR99ptotHelp    , EcaR99ptotOperators    , 1 , CDI_REAL , 2  , 1  });
  add_module("EcaRx1day"     , {EcaRx1day     , EcaRx1dayHelp     , EcaRx1dayOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaRx5day"     , {EcaRx5day     , EcaRx5dayHelp     , EcaRx5dayOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("EcaSdii"       , {EcaSdii       , EcaSdiiHelp       , EcaSdiiOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Fdns"          , {Fdns          , FdnsHelp          , FdnsOperators          , 1 , CDI_REAL , 2  , 1  });
  add_module("Strwin"        , {Strwin        , StrwinHelp        , StrwinOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Strbre"        , {Strbre        , StrbreHelp        , StrbreOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Strgal"        , {Strgal        , StrgalHelp        , StrgalOperators        , 1 , CDI_REAL , 1  , 1  });
  add_module("Hurr"          , {Hurr          , HurrHelp          , HurrOperators          , 1 , CDI_REAL , 1  , 1  });
  // add_module("Hi"         , { Hi           , {}                , HiOperators            , 1 , CDI_REAL , 3  , 1   });
  add_module("Wct"           , {Wct           , WctHelp           , WctOperators           , 1 , CDI_REAL , 2  , 1  });
  add_module("Magplot"       , {Magplot       , MagplotHelp       , MagplotOperators       , 1 , CDI_REAL , 1  , 1  });
  add_module("Magvector"     , {Magvector     , MagvectorHelp     , MagvectorOperators     , 1 , CDI_REAL , 1  , 1  });
  add_module("Maggraph"      , {Maggraph      , MaggraphHelp      , MaggraphOperators      , 1 , CDI_REAL , -1 , 1  });
  // HIRLAM_EXTENSIONS
  add_module( "Samplegrid"   , { Samplegrid   , SamplegridHelp    , SamplegridOperators    , 1 , CDI_REAL , 1  , 1 });
  add_module( "Selmulti "    , { Selmulti     , SelmultiHelp      , SelmultiOperators      , 1 , CDI_REAL , 1  , 1 });
  add_module( "WindTrans"    , { WindTrans    , WindTransHelp     , WindTransOperators     , 1 , CDI_REAL , 1  , 1 });

  init_aliases();
}

/**
 * Initializes all hardcoded aliases
 */
void init_aliases()
{
  add_alias("afterburner"     , "after");
  add_alias("anomaly"         , "ymonsub");
  add_alias("deltap_fl"       , "deltap");
  add_alias("diffv"           , "diffn");
  add_alias("covar0"          , "timcovar");
  add_alias("covar0r"         , "fldcovar");
  add_alias("gather"          , "collgrid");
  add_alias("geopotheight"    , "gheight");
  add_alias("globavg"         , "fldavg");
  add_alias("import_grads"    , "import_binary");
  add_alias("infos"           , "sinfo");
  add_alias("infov"           , "infon");
  add_alias("intgrid"         , "intgridbil");
  add_alias("log"             , "ln");
  add_alias("lmean"           , "ymonmean");
  add_alias("lmmean"          , "ymonmean");
  add_alias("lmavg"           , "ymonavg");
  add_alias("lmstd"           , "ymonstd");
  add_alias("lsmean"          , "yseasmean");
  add_alias("chvar"           , "chname");
  add_alias("nvar"            , "npar");
  add_alias("outputkey"       , "outputtab");
  add_alias("vardes"          , "codetab");
  add_alias("pardes"          , "codetab");
  add_alias("selvar"          , "selname");
  add_alias("delvar"          , "delname");
  add_alias("selindex"        , "selgridcell");
  add_alias("remapcon1"       , "remaplaf");
  add_alias("remapdis1"       , "remapnn");
  add_alias("scatter"         , "distgrid");
  add_alias("showvar"         , "showname");
  add_alias("selgridname"     , "selgrid");
  add_alias("setvar"          , "setname");
  add_alias("setpartabv"      , "setpartabn");
  add_alias("setpartab"       , "setcodetab");
  add_alias("sinfov"          , "sinfon");
  add_alias("sortvar"         , "sortname");
  add_alias("splitvar"        , "splitname");
  add_alias("sort"            , "timsort");
  add_alias("eca_r1mm"        , "eca_rr1");
  add_alias("fpressure"       , "pressure_fl");
  add_alias("hpressure"       , "pressure_hl");
  add_alias("ensrkhist_space" , "ensrkhistspace");
  add_alias("ensrkhist_time"  , "ensrkhisttime");
  add_alias("gridverify"      , "verifygrid");
  add_alias("outputcenter"    , "gmtxyz");
  add_alias("outputbounds"    , "gmtcells");
  add_alias("selseas"         , "selseason");
  add_alias("selmon"          , "selmonth");
}


int main(int argc, char *argv[])
{
  int lstop = FALSE;
  int noff = 0;
  int status = 0;
  const char *operatorArg = NULL;

  cdo_init_is_tty();

  memExitOnError();

  _Verbose = 1;
  CDO_Reduce_Dim = 0;

  /* mallopt(M_MMAP_MAX, 0); */

  setCommandLine(argc, argv);

  Cdo::progname = getProgname(argv[0]);
  if ( strncmp(Cdo::progname, "cdo", 3) == 0 && strlen(Cdo::progname) > 3 ) noff = 3;
  if ( noff ) setDefaultFileType(Cdo::progname+noff, 0);

  get_env_vars();
  init_modules();
  status = parse_options_long(argc, argv);

  if ( status != 0 ) return -1;

  cdo_set_options();
#ifdef DEBUG
    CdoDebug::CdoStartMessage();
    MESSAGE(CdoDebug::argvToString(argc,(const char**) argv));
#endif
  if ( Debug || Version ) cdo_version();

  if ( Debug )
    {
      fprintf(stderr, "stdin_is_tty:   %d\n", stdin_is_tty);
      fprintf(stderr, "stdout_is_tty:  %d\n", stdout_is_tty);
      fprintf(stderr, "stderr_is_tty:  %d\n", stderr_is_tty);
      print_system_info();
    }

  check_stacksize();

  if ( Debug ) 
  { 
      print_pthread_info();
      //      fprintf(stderr, "C++ max thread      = %u\n", std::thread::hardware_concurrency());
  }

#ifdef  _OPENMP
  if ( numThreads <= 0 ) numThreads = 1;
  omp_set_num_threads(numThreads);

  if ( Debug )
    {
      fprintf(stderr, "OMP num procs       = %d\n", omp_get_num_procs());
      fprintf(stderr, "OMP max threads     = %d\n", omp_get_max_threads());
      fprintf(stderr, "OMP num threads     = %d\n", omp_get_num_threads());
#if defined(HAVE_OPENMP3)
      fprintf(stderr, "OMP thread limit    = %d\n", omp_get_thread_limit());
      omp_sched_t kind;
      int modifer;
      omp_get_schedule(&kind, &modifer);
      fprintf(stderr, "OMP schedule        = %d (1:static; 2:dynamic; 3:guided; 4:auto)\n", (int) kind);
#endif
#ifdef  HAVE_OPENMP4
      fprintf(stderr, "OMP proc bind       = %d (0:false; 1:true; 2:master; 3:close; 4:spread)\n", (int) omp_get_proc_bind());
#if !defined(__ICC)
      fprintf(stderr, "OMP num devices     = %d\n", omp_get_num_devices());
#endif
#endif
    }

  Threading::ompNumThreads = omp_get_max_threads();
  if ( omp_get_max_threads() > omp_get_num_procs() )
    fprintf(stderr, "Warning: Number of OMP threads is greater than number of Cores=%d!\n", omp_get_num_procs());
  if ( Threading::ompNumThreads < numThreads )
    fprintf(stderr, "Warning: omp_get_max_threads() returns %d!\n", Threading::ompNumThreads);
  if ( cdoVerbose )
    {
      fprintf(stderr, " OpenMP:  num_procs=%d  max_threads=%d", omp_get_num_procs(), omp_get_max_threads());
#ifdef  HAVE_OPENMP4
#if !defined(__ICC)
      fprintf(stderr, "  num_devices=%d", omp_get_num_devices());
#endif
#endif
      fprintf(stderr, "\n");
    }
#else
  if ( numThreads > 0 )
    {
      fprintf(stderr, "Option -P failed, OpenMP support not compiled in!\n");
      return -1;
    }
#endif
      std::vector<std::string> new_argv(&argv[CDO_optind], argv + argc);
      new_argv = expandWildCards(new_argv);

      ///*TEMP*/ // should not be needed when std::string is standart string 
      std::vector<char*> new_cargv(new_argv.size());
      for(unsigned long i = 0; i < new_argv.size(); i++)
      {
          new_cargv[i] = strdup(new_argv[i].c_str());
      }
      //temprorary end

  if ( CDO_optind >= argc )
      {
      if ( ! Version && ! Help )
        {
          fprintf(stderr, "\nNo operator given!\n\n");
          cdo_usage();
          status = 1;
        }

      if ( Help ) cdo_usage();
      lstop = TRUE;
    }

  if ( lstop ) return status;

  if ( cdoDefaultTableID != CDI_UNDEFID ) cdiDefTableID(cdoDefaultTableID);

  extern int (*proj_lonlat_to_lcc_func)();
  proj_lonlat_to_lcc_func = (int (*)()) proj_lonlat_to_lcc;
  extern int (*proj_lcc_to_lonlat_func)();
  proj_lcc_to_lonlat_func = (int (*)()) proj_lcc_to_lonlat;
  
  const char *operatorName = get_original(getOperatorName(argv[CDO_optind]));

  if ( Help )
    {
      cdoPrintHelp(operatorHelp(operatorName));
    }
  else if ( cdoExpMode == CDO_EXP_LOCAL )
    {
      exp_run(argc, argv, cdoExpName);
    }
  else
    {
      if ( Debug )
        {
          if ( DebugLevel == 0 ) DebugLevel = 1;
          cdiDebug(DebugLevel);
          CdoDebug::SetDebug(DebugLevel);
        }

      timer_total  = timer_new("total");
      timer_read   = timer_new("read");
      timer_write  = timer_new("write");

      timer_start(timer_total);


#ifdef CUSTOM_MODULES
      load_custom_modules("custom_modules");
      getProcess(0)->m_module.func(getProcess(0));
      close_library_handles();
#else
      createProcesses(new_argv.size(),(const char**) &new_cargv[0] );
      cdoValidateProcesses();
      runProcesses();
#endif
      clearProcesses();

      timer_stop(timer_total);

      if ( cdoTimer ) timer_report();
    }

  if ( cdoVarnames )
    {
      if ( cdoNumVarnames ) Free(cdoVarnames[0]);
      Free(cdoVarnames);
    }

  /* problems with alias!!! if ( operatorName ) Free(operatorName); */ 

  /* malloc_stats(); */

  if ( cdoGridSearchDir ) Free(cdoGridSearchDir);

  if ( CDO_Rusage ) cdo_rusage();
#ifdef DEBUG
    CdoDebug::CdoEndMessage();
#endif
  return status;
}
