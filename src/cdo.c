/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2008 Uwe Schulzweida, Uwe.Schulzweida@zmaw.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#if  defined  (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
/*#include <malloc.h>*/ /* mallopt and malloc_stats */
#if ! defined (SX)
#include <sys/time.h>       /* getrlimit */
#include <sys/resource.h>   /* getrlimit */
#endif
#include <unistd.h>         /* getrlimit sysconf */

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"


#if  defined  (HAVE_LIBPTHREAD)
#include "pstream_int.h"
#include "pthread_debug.h"
#endif

#include "modules.h"
#include "util.h"

#if defined (_OPENMP)
#  include <omp.h>
#endif

#if ! defined (VERSION)
#  define  VERSION  "0.0.1"
#endif

char CDO_Version[] = "Climate Data Operators version "VERSION" (http://www.mpimet.mpg.de/cdo)";


char *Progname;

int ompNumThreads = 1;

int cdoDefaultFileType  = CDI_UNDEFID;
int cdoDefaultDataType  = CDI_UNDEFID;
int cdoDefaultTimeType  = CDI_UNDEFID;
int cdoDefaultByteorder = CDI_UNDEFID;
int cdoDefaultTableID   = CDI_UNDEFID;

int cdoDisableFilesuffix = FALSE;
int cdoDisableHistory = FALSE;
int cdoZtype       = COMPRESS_NONE;
int cdoZlevel      = 0;
int cdoLogOff      = FALSE;
int cdoSilentMode  = FALSE;
int cdoRegulargrid = FALSE;
int cdoBenchmark   = FALSE;
int cdoTimer       = FALSE;
int cdoVerbose     = FALSE;
int cdoDebug       = 0;

int cdoExpMode     = -1;
char *cdoExpName   = NULL;
void exp_run(int argc, char *argv[], char *cdoExpName);


int timer_total, timer_read, timer_write;
int timer_remap, timer_remap_con, timer_remap_con2, timer_remap_con3;


#define PRINT_RLIMIT(resource) \
      { \
	int status; \
	struct rlimit rlim; \
	status = getrlimit(resource, &rlim); \
	if ( status == 0 ) \
	  { \
	    if ( sizeof(rlim_t) > sizeof(long) ) \
	      { \
		fprintf(stderr, "CUR %-15s = %llu\n", #resource, (long long) rlim.rlim_cur); \
		fprintf(stderr, "MAX %-15s = %llu\n", #resource, (long long) rlim.rlim_max); \
	      } \
	    else \
	      { \
		fprintf(stderr, "CUR %-15s = %lu\n", #resource, (long) rlim.rlim_cur); \
		fprintf(stderr, "MAX %-15s = %lu\n", #resource, (long) rlim.rlim_max); \
	      } \
	  } \
      }


static void version(void)
{
  fprintf(stderr, "%s\n", CDO_Version);
#if defined (COMPILER)
  fprintf(stderr, "Compiler: %s\n", COMPILER);
#endif
#if defined (COMP_VERSION)
  fprintf(stderr, " version: %s\n", COMP_VERSION);
#endif
#if defined (HAVE_LIBSZ) || defined (HAVE_LIBDRMAA) || defined (HAVE_LIBCURL)
  fprintf(stderr, "    with:");
#if defined (HAVE_LIBSZ)
  fprintf(stderr, " libsz");
#endif
#if defined (HAVE_LIBDRMAA)
  fprintf(stderr, " libdrmaa");
#endif
#if defined (HAVE_LIBCURL)
  fprintf(stderr, " libcurl");
#endif
  fprintf(stderr, "\n");
#endif
#if defined (USER_NAME) && defined(HOST_NAME)
  fprintf(stderr, "Compiled: by %s on %s %s %s\n", USER_NAME, HOST_NAME, __DATE__, __TIME__);
#endif
  cdiPrintVersion();
  fprintf(stderr, "\n");
}


static void usage(void)
{
  int id = 0;
  char *name;

  /*  fprintf(stderr, "%s\n", CDO_Version);*/
  /*  fprintf(stderr, "\n");*/
  fprintf(stderr, "usage : cdo  [Options]  [Operators]  [Operands]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  Options:\n");
  fprintf(stderr, "    -a             Convert from a relative to an absolute time axis\n");
  fprintf(stderr, "    -b <nbits>     Set the number of bits for the output precision\n");
  fprintf(stderr, "                   (32/64 for nc, nc2, srv, ext, ieg; 1 - 32 for grb)\n");
  fprintf(stderr, "    -f <format>    Format of the output file. (grb, nc, nc2, srv, ext or ieg)\n");
  fprintf(stderr, "    -g <grid>      Grid name or file. Available grids: \n");
  fprintf(stderr, "                   t<RES>grid, t<RES>spec, r<NX>x<NY>, g<NX>x<NY>, gme<NI>\n");
  fprintf(stderr, "    -h             Help information for the operators\n");
  /*
  fprintf(stderr, "    -i <inst>      Institution name/file\n");
  fprintf(stderr, "                   Predefined instituts: ");
  for ( id = 0; id < institutInqNumber; id++ )
    if ( (name = institutInqNamePtr(id)) )
      fprintf(stderr, " %s", name);
  fprintf(stderr, "\n");
  */
  /* fprintf(stderr, "    -l <level>     Level file\n"); */
  fprintf(stderr, "    -m <missval>   Set the default missing value (default: %g)\n", cdiInqMissval());
  /*
  fprintf(stderr, "    -p <prec>      Set the precision of the output data in bytes\n");
  fprintf(stderr, "                   (4/8 for nc, nc2, srv, ext, ieg; 1/2/3 for grb)\n");
  */
  fprintf(stderr, "    -R             Convert GRIB data from reduced to regular grid\n");
  fprintf(stderr, "    -r             Convert from an absolute to a relative time axis\n");
  fprintf(stderr, "    -s             Silent mode\n");
  fprintf(stderr, "    -t <partab>    Parameter table name or file\n");
  fprintf(stderr, "                   Predefined tables: ");
  for ( id = 0; id < tableInqNumber(); id++ )
    if ( (name = tableInqNamePtr(id)) )
      fprintf(stderr, " %s", name);
  fprintf(stderr, "\n");

  fprintf(stderr, "    -V             Print the version number\n");
  fprintf(stderr, "    -v             Print extra details for some operators\n");
  fprintf(stderr, "    -z szip        Compress GRIB records with szip\n");
  fprintf(stderr, "\n");

  fprintf(stderr, "  Operators:\n");
  operatorPrintAll();

  fprintf(stderr, "\n");
  fprintf(stderr, "  Operands:\n");
  fprintf(stderr, "    The path name of the input/output files for the operator.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  CDO version %s, Copyright (C) 2003-2008 Uwe Schulzweida\n", VERSION);
  fprintf(stderr, "  Available from http://www.mpimet.mpg.de/cdo\n");
  fprintf(stderr, "  This is free software and comes with ABSOLUTELY NO WARRANTY\n");
  fprintf(stderr, "  Report bugs to Uwe.Schulzweida@zmaw.de\n");
}


static void cdoPrintHelp(char *phelp[]/*, char *xoperator*/)
{
  if ( phelp == NULL )
    printf("No help available for this operator!\n");
  else
    {
      int lprint;
      while ( *phelp )
	{
	  lprint = TRUE;
	  if ( *phelp[0] == '\0' )
	    if ( *(phelp+1) )
	      if ( *(phelp+1)[0] == ' ' ) lprint = FALSE;
	  
	  if ( lprint ) printf("%s\n", *phelp);

	  phelp++;
	}
    }
}


static void cdoSetDebug(int level)
{
  /*
    level   0: off
    level   1: on
    level   2: cdi
    level   4: memory
    level   8: file
    level  16: format
    level  32: cdo
    level  64: stream
    level 128: pipe
    level 256: pthread
   */
  cdiDebug(level);

  if ( level == 1 || level &  32 ) cdoDebug = 1;
  if ( level == 1 || level &  64 ) pstreamDebug(1);
#if  defined  (HAVE_LIBPTHREAD)
  if ( level == 1 || level & 128 ) pipeDebug(1);
  if ( level == 1 || level & 256 ) Pthread_debug(1);
#endif
}


static int cdoOptind = 1;
static char *cdoOptarg;

static int cdoGetopt(int argc, char * const argv[], const char *optstring)
{
  static int optpos = 0;
  int optval = -1, value;
  int opthasarg = 0;
  int optstrlen = strlen(optstring);
  int iargc;

  cdoOptarg = NULL;

  while ( optpos < optstrlen && cdoOptind < argc )
    {
      value = optstring[optpos];
      optpos++;
      if ( optstring[optpos] == ':' )
	{
	  opthasarg = 1;
	  optpos++;
	}
      else
	opthasarg = 0;

      for ( iargc = 1; iargc < argc; iargc++ )
	{
	  if ( *argv[iargc] == '-' && strlen(argv[iargc]) == 2 )
	    {
	      if ( (argv[iargc][1]) == value )
		{
		  optval = value;
		  cdoOptind++;
		  if ( opthasarg )
		    {
		      cdoOptarg = argv[iargc+1];
		      cdoOptind++;
		    }
		  break;
		}
	    }
	}
      if ( iargc < argc ) break;
    }

  if ( opthasarg && cdoOptarg == NULL ) optval = ':';

  return (optval);
}


#undef  IsBigendian
#define IsBigendian()  ( u_byteorder.c[sizeof(long) - 1] )


static void setDefaultDataType(char *datatypestr)
{
  static union {unsigned long l; unsigned char c[sizeof(long)];} u_byteorder = {1};
  int datatype = -1;

  if ( isdigit((int) *datatypestr) )
    {
      datatype = atoi(datatypestr);
      if ( datatype < 10 )
	datatypestr += 1;
      else
	datatypestr += 2;

      if      ( datatype > 0 && datatype < 32 ) cdoDefaultDataType = datatype;
      else if ( datatype == 32 )
	{
	  if ( cdoDefaultFileType == FILETYPE_GRB )
	    cdoDefaultDataType = DATATYPE_PACK32;
	  else
	    cdoDefaultDataType = DATATYPE_FLT32;
	}
      else if ( datatype == 64 ) cdoDefaultDataType = DATATYPE_FLT64;
      else
	{
	  fprintf(stderr, "Unsupported datatype %d!\n", datatype);
	  fprintf(stderr, "Use 32/64 for filetype nc, srv, ext, ieg and 1-32 for grb.\n");
	  exit(EXIT_FAILURE);
	}
    }

  if ( *datatypestr == 'l' || *datatypestr == 'L' )
    {
      if ( IsBigendian() ) cdoDefaultByteorder = CDI_LITTLEENDIAN;
      datatypestr++;
    }

  if ( *datatypestr == 'b' || *datatypestr == 'B' )
    {
      if ( ! IsBigendian() ) cdoDefaultByteorder = CDI_BIGENDIAN;
      datatypestr++;
    }
}


static void setDefaultDataTypeByte(char *datatypestr)
{
  static union {unsigned long l; unsigned char c[sizeof(long)];} u_byteorder = {1};
  int datatype = -1;

  if ( isdigit((int) *datatypestr) )
    {
      datatype = atoi(datatypestr);
      datatypestr++;

      if      ( datatype == 1 ) cdoDefaultDataType = DATATYPE_PACK8;
      else if ( datatype == 2 ) cdoDefaultDataType = DATATYPE_PACK16;
      else if ( datatype == 3 ) cdoDefaultDataType = DATATYPE_PACK24;
      else if ( datatype == 4 ) cdoDefaultDataType = DATATYPE_FLT32;
      else if ( datatype == 8 ) cdoDefaultDataType = DATATYPE_FLT64;
      else
	{
	  fprintf(stderr, "Unsupported datatype %d!\n", datatype);
	  fprintf(stderr, "Use 4/8 for filetype nc, srv, ext, ieg and 1/2/3 for grb.\n");
	  exit(EXIT_FAILURE);
	}
    }

  if ( *datatypestr == 'l' || *datatypestr == 'L' )
    {
      if ( IsBigendian() ) cdoDefaultByteorder = CDI_LITTLEENDIAN;
      datatypestr++;
    }

  if ( *datatypestr == 'b' || *datatypestr == 'B' )
    {
      if ( ! IsBigendian() ) cdoDefaultByteorder = CDI_BIGENDIAN;
      datatypestr++;
    }
}


void setDefaultFileType(char *filetypestr, int labort)
{
  if ( filetypestr )
    {
      char *ftstr = filetypestr;

      if      ( strncmp(filetypestr, "grb", 3)  == 0 ) { ftstr += 3; cdoDefaultFileType = FILETYPE_GRB; }
      else if ( strncmp(filetypestr, "nc2", 3)  == 0 ) { ftstr += 3; cdoDefaultFileType = FILETYPE_NC2; }
      else if ( strncmp(filetypestr, "nc",  2)  == 0 ) { ftstr += 2; cdoDefaultFileType = FILETYPE_NC;  }
      else if ( strncmp(filetypestr, "srv", 3)  == 0 ) { ftstr += 3; cdoDefaultFileType = FILETYPE_SRV; }
      else if ( strncmp(filetypestr, "ext", 3)  == 0 ) { ftstr += 3; cdoDefaultFileType = FILETYPE_EXT; }
      else if ( strncmp(filetypestr, "ieg", 3)  == 0 ) { ftstr += 3; cdoDefaultFileType = FILETYPE_IEG; }
      else
	{
	  if ( labort )
	    {
	      fprintf(stderr, "Unsupported filetype %s!\n", filetypestr);
	      fprintf(stderr, "Available filetypes: grb, nc, nc2, srv, ext and ieg\n");
	      exit(EXIT_FAILURE);
	    }
	  else
	    {
	      return;
	    }
	}

      if ( cdoDefaultFileType != CDI_UNDEFID && *ftstr != 0 )
	{
	  if ( *ftstr == '_' )
	    {
	      ftstr++;

	      setDefaultDataType(ftstr);
	    }
	  else
	    {
	      fprintf(stderr, "Unexpected character >%c< in file type >%s<!\n", *ftstr, filetypestr);
	      fprintf(stderr, "Use format[_nbits] with:\n");
	      fprintf(stderr, "    format = grb, nc, nc2, srv, ext or ieg\n");
	      fprintf(stderr, "    nbits  = 32/64 for nc, nc2, srv, ext, ieg; 1 - 32 for grb\n");
	      exit(EXIT_FAILURE);
	    }
	}
    }
}


int cdoFiletype(void)
{
  if ( cdoDefaultFileType == CDI_UNDEFID )
    {
      cdoDefaultFileType = FILETYPE_GRB;
      if ( ! cdoSilentMode )
	cdoPrint("Set default filetype to GRIB");
    }

  return (cdoDefaultFileType);
}


void defineCompress(const char *arg)
{
  size_t len = strlen(arg);

  if      ( strncmp(arg, "szip", len) == 0 )
    {
      cdoZtype  = COMPRESS_SZIP;
      cdoZlevel = 0;
    }
  else if ( strncmp(arg, "gzip", len) == 0 )
    {
      cdoZtype  = COMPRESS_GZIP;
      cdoZlevel = 6;
    }
  else
    fprintf(stderr, "Compression %s unsupported!\n", arg);
}


int main(int argc, char *argv[])
{
  static char func[] = "main";
  int c;
  int Debug = 0;
  int Version = 0;
  int Help = 0;
  int DebugLevel = 0;
  int lstop = FALSE;
  int noff = 0;
  int status = 0;
  char *operatorName = NULL;
  char *operatorArg = NULL;
  char *argument = NULL;
  extern int dmemory_ExitOnError;

  dmemory_ExitOnError = 1;

  /* mallopt(M_MMAP_MAX, 0); */
 
  setCommandLine(argc, argv);

  Progname = getProgname(argv[0]);

  if ( strncmp(Progname, "cdo", 3) == 0 && strlen(Progname) > 3 ) noff = 3;

  /* old versions !!!! */
  if ( strncmp(Progname, "gdo", 3) == 0 && strlen(Progname) > 3 ) noff = 3;
  if ( strncmp(Progname, "gm",  2) == 0 && strlen(Progname) > 2 ) noff = 2;

  if ( noff ) setDefaultFileType(Progname+noff, 0);

  while ( (c = cdoGetopt(argc, argv, "f:b:e:p:g:i:l:m:t:D:z:aBdhRrsTVvZ")) != -1 )
    {
      switch (c)
	{
	case 'a':
	  cdoDefaultTimeType = TAXIS_ABSOLUTE;
	  break;
	case 'b':
	  setDefaultDataType(cdoOptarg);
	  break;
	case 'B':
	  cdoBenchmark = TRUE;
	  break;
	case 'd':
	  Debug = 1;
	  break;
	case 'D':
	  Debug = 1;
	  DebugLevel = atoi(cdoOptarg);
	  break;
	case 'e':
	  {
	  char host[1024];
	  gethostname(host, sizeof(host));
	  cdoExpName = cdoOptarg;
	  /* printf("host: %s %s\n", host, cdoExpName); */
	  if ( strcmp(host, cdoExpName) == 0 )
	    cdoExpMode = CDO_EXP_REMOTE;
	  else
            cdoExpMode = CDO_EXP_LOCAL;
          break;
	  }
	case 'f':
	  setDefaultFileType(cdoOptarg, 1);
	  break;
	case 'g':
	  defineGrid(cdoOptarg);
	  break;
	case 'h':	
	  Help = 1;
	  break;
	case 'i':
	  defineInstitution(cdoOptarg);
	  break;
	case 'l':
	  defineZaxis(cdoOptarg);
	  break;
	case 'm':
	  cdiDefMissval(atof(cdoOptarg));
	  break;
	case 'p':
	  setDefaultDataTypeByte(cdoOptarg);
	  break;
	case 'R':
	  cdoRegulargrid = TRUE;
	  cdiDefGlobal("CD_REGULARGRID", TRUE);
	  break;
	case 'r':
	  cdoDefaultTimeType = TAXIS_RELATIVE;
	  break;
	case 's':
	  cdoSilentMode = TRUE;
	  break;
	case 'T':
	  cdoTimer = TRUE;
	  break;
	case 't':
	  cdoDefaultTableID = defineTable(cdoOptarg);
	  break;
	case 'V':
	  Version = 1;
	  break;
	case 'v':
	  cdoVerbose = TRUE;
	  break;
	case 'z':
	  defineCompress(cdoOptarg);
          break;
	case ':':
	  fprintf(stderr, "\nmissing parameter for one of the options\n\n");	  
	  Help = 1;
	  break;
	}
    }

  {
      char *envstr;

      envstr = getenv("CDO_LOG_OFF");
      if ( envstr )
	{
	  if ( atoi(envstr) == 1 )
	    {
	      cdoLogOff = TRUE;
	      if ( cdoVerbose )
		fprintf(stderr, "CDO_LOG_OFF         = %s\n", envstr);
	    }
	}

      envstr = getenv("CDO_DISABLE_HISTORY");
      if ( envstr )
	{
	  if ( atoi(envstr) == 1 )
	    {
	      cdoDisableHistory = TRUE;
	      if ( cdoVerbose )
		fprintf(stderr, "CDO_DISABLE_HISTORY = %s\n", envstr);
	    }
	}

      envstr = getenv("CDO_DISABLE_FILESUFFIX");
      if ( envstr )
	{
	  if ( atoi(envstr) == 1 )
	    {
	      cdoDisableFilesuffix = TRUE;
	      if ( cdoVerbose )
		fprintf(stderr, "CDO_DISABLE_FILESUFFIX = %s\n", envstr);
	    }
	}
  }

  if ( Debug || Version ) version();

  if ( Debug )
    {
      char *envstr;

      if ( DebugLevel == 0 ) DebugLevel = 1;
      cdoSetDebug(DebugLevel);
      fprintf(stderr, "\n");
      fprintf(stderr, "cdoDefaultFileType  = %d\n", cdoDefaultFileType);
      fprintf(stderr, "cdoDefaultDataType  = %d\n", cdoDefaultDataType);
      fprintf(stderr, "cdoDefaultByteorder = %d\n", cdoDefaultByteorder);
      fprintf(stderr, "cdoDefaultTableID   = %d\n", cdoDefaultTableID);
      fprintf(stderr, "\n");

      envstr = getenv("HOSTTYPE");
      if ( envstr ) fprintf(stderr, "HOSTTYPE            = %s\n", envstr);
      envstr = getenv("VENDOR");
      if ( envstr ) fprintf(stderr, "VENDOR              = %s\n", envstr);
      envstr = getenv("OSTYPE");
      if ( envstr ) fprintf(stderr, "OSTYPE              = %s\n", envstr);
      envstr = getenv("MACHTYPE");
      if ( envstr ) fprintf(stderr, "MACHTYPE            = %s\n", envstr);
      fprintf(stderr, "\n");

#if defined (HAVE_MMAP)
      fprintf(stderr, "HAVE_MMAP\n");
#endif
#if defined (HAVE_MEMORY_H)
      fprintf(stderr, "HAVE_MEMORY_H\n");
#endif
      fprintf(stderr, "\n");

#if defined (__STDC__)
      fprintf(stderr, "STD ANSI C          = %d\n", __STDC__);
#endif
#if defined (__STD_VERSION__)
      fprintf(stderr, "STD VERSION         = %ld\n", __STD_VERSION__);
#endif
#if defined (__STDC_VERSION__)
      fprintf(stderr, "STDC VERSION        = %ld\n", __STDC_VERSION__);
#endif
#if defined (__STD_HOSTED__)
      fprintf(stderr, "STD HOSTED          = %d\n", __STD_HOSTED__);
#endif
#if defined (FLT_EVAL_METHOD)
      fprintf(stderr, "FLT_EVAL_METHOD     = %d\n", FLT_EVAL_METHOD);
#endif
#if defined (FP_FAST_FMA)
      fprintf(stderr, "FP_FAST_FMA         = defined\n");
#endif
      fprintf(stderr, "\n");

#if defined (_SC_VERSION)
      fprintf(stderr, "POSIX.1 VERSION     = %ld\n", sysconf(_SC_VERSION));
#endif
#if defined (_SC_ARG_MAX)
      fprintf(stderr, "POSIX.1 ARG_MAX     = %ld\n", sysconf(_SC_ARG_MAX));
#endif
#if defined (_SC_CHILD_MAX)
      fprintf(stderr, "POSIX.1 CHILD_MAX   = %ld\n", sysconf(_SC_CHILD_MAX));
#endif
#if defined (_SC_STREAM_MAX)
      fprintf(stderr, "POSIX.1 STREAM_MAX  = %ld\n", sysconf(_SC_STREAM_MAX));
#endif
#if defined (_SC_OPEN_MAX)
      fprintf(stderr, "POSIX.1 OPEN_MAX    = %ld\n", sysconf(_SC_OPEN_MAX));
#endif
#if defined (_SC_PAGESIZE)
      fprintf(stderr, "POSIX.1 PAGESIZE    = %ld\n", sysconf(_SC_PAGESIZE));
#endif

      fprintf(stderr, "\n");

#if defined (RLIMIT_FSIZE)
      PRINT_RLIMIT(RLIMIT_FSIZE);
#endif
#if defined (RLIMIT_NOFILE)
      PRINT_RLIMIT(RLIMIT_NOFILE);
#endif
#if defined (RLIMIT_STACK)
      PRINT_RLIMIT(RLIMIT_STACK);
#endif
      fprintf(stderr, "\n");
    }

#if defined (RLIMIT_STACK)
  {
#define  MIN_STACK_SIZE  67108864L  /* 64MB */
    int status;
    struct rlimit rlim;
    rlim_t min_stack_size = MIN_STACK_SIZE;

    status = getrlimit(RLIMIT_STACK, &rlim);

    if ( status == 0 )
      {
	if ( min_stack_size > rlim.rlim_max ) min_stack_size = rlim.rlim_max;
	if ( rlim.rlim_cur < min_stack_size )
	  {
	    rlim.rlim_cur = min_stack_size;

	    status = setrlimit(RLIMIT_STACK, &rlim);
	    if ( Debug )
	      {
		if ( status == 0 )
		  {
		    fprintf(stderr, "Set stack size to %ld\n", MIN_STACK_SIZE);
		    PRINT_RLIMIT(RLIMIT_STACK);
		  }
		else
		  fprintf(stderr, "Set stack size to %ld failed!\n", MIN_STACK_SIZE);
	      }
	  }
      }
  }
#endif

  if ( Debug )
    {
      print_pthread_info();
    }

#if defined (_OPENMP)
  ompNumThreads = omp_get_max_threads();
  fprintf(stderr, " OpenMP:  num_procs = %d  max_threads = %d\n",
	  omp_get_num_procs(), omp_get_max_threads());
#endif


  if ( cdoOptind < argc )
    {
      operatorArg = argv[cdoOptind];
      argument = makeArgument(argc-cdoOptind, &argv[cdoOptind]);
    }
  else
    {
      if ( ! Version && ! Help )
	{
	  fprintf(stderr, "\nno operator given\n\n");
	  usage();
	  status = 1;
	}

      if ( Help ) usage();
      lstop = TRUE;
    }

  if ( lstop ) return (status);

  if ( cdoDefaultTableID != CDI_UNDEFID ) cdiDefTableID(cdoDefaultTableID);

  operatorName = getOperatorName(operatorArg);

  if ( Help )
    {
      cdoPrintHelp(operatorHelp(operatorName)/*, operatorName*/);
    }
  else if ( cdoExpMode == CDO_EXP_LOCAL )
    {
      exp_run(argc, argv, cdoExpName);
    }
  else
    {
      if ( cdoTimer )
	{
	  timer_total      = timer_new("total");
	  timer_read       = timer_new("read");
	  timer_write      = timer_new("write");
	  timer_remap      = timer_new("remap");
	  timer_remap_con  = timer_new("remap con");
	  timer_remap_con2 = timer_new("remap con2");
	  timer_remap_con3 = timer_new("remap con3");
	}

      if ( cdoTimer ) timer_start(timer_total);

      operatorModule(operatorName)(argument);

      if ( cdoTimer ) timer_stop(timer_total);

      if ( cdoTimer ) timer_report();
    }

  if ( argument ) free(argument);
  /* problems with alias!!! if ( operatorName ) free(operatorName); */ 

  /* malloc_stats(); */


  return (status);
}
