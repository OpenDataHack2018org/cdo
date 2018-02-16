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
#ifdef  HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdarg.h>
#include <errno.h>
#include <cdi.h>

#include "cdo_int.h"
#include "process_int.h"
#include "error.h"
#include "text.h"
#include "cdoOptions.h"

void pstreamCloseAll(void);

void cdiOpenError(int cdiErrno, const char *fmt, const char *path)
{	
  printf("\n");
  set_text_color(stderr, RESET, RED);
   fprintf(stderr, "%s: ", processInqPrompt());
  reset_text_color(stderr);
  set_text_color(stderr, RESET, BLACK);
  fprintf(stderr, fmt, path);
  reset_text_color(stderr);
   fprintf(stderr, "\n");

  fprintf(stderr, "%s\n", cdiStringError(cdiErrno));

  if ( cdiErrno == CDI_ELIBNAVAIL )
    {
      int byteorder;
      int filetype = cdiGetFiletype(path, &byteorder);

      switch (filetype)
	{
	case CDI_FILETYPE_GRB:
          break;
	case CDI_FILETYPE_GRB2:
          fprintf(stderr, "To create a CDO application with GRIB2 support use: ./configure --with-eccodes=<ECCODES root directory> ...\n");
          break;
	case CDI_FILETYPE_SRV:
          break;
	case CDI_FILETYPE_EXT:
          break;
	case CDI_FILETYPE_IEG:
          break;
	case CDI_FILETYPE_NC:
	case CDI_FILETYPE_NC2:
	case CDI_FILETYPE_NC4:
	case CDI_FILETYPE_NC4C:
	case CDI_FILETYPE_NC5:
          {
            const char *ncv = (filetype == CDI_FILETYPE_NC4 || filetype == CDI_FILETYPE_NC4C) ? "4" :
              ((filetype == CDI_FILETYPE_NC2) ? "2" :
               ((filetype == CDI_FILETYPE_NC5) ? "5" : ""));
#if defined HAVE_LIBNETCDF
            fprintf(stderr, "CDO was build with a NetCDF version which doesn't support NetCDF%s data!\n", ncv);
#else
            fprintf(stderr, "To create a CDO application with NetCDF%s support use: ./configure --with-netcdf=<NetCDF%s root directory> ...\n", ncv, ncv);
#endif
            break;
          }
	default:
          break;
	}
    }
  
  if ( _ExitOnError )
    {
      pstreamCloseAll();
      exit(EXIT_FAILURE);
    }
}

void cdoAbort(const char *fmt, ...)
{
  printf("\n");
  set_text_color(stderr, RESET, RED);
   fprintf(stderr, "%s (Abort): ", processInqPrompt());
  reset_text_color(stderr);
  
  set_text_color(stderr, RESET, BLACK);
  va_list args;	
  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  reset_text_color(stderr);
  fprintf(stderr, "\n");

  if ( _ExitOnError )
    {
      pstreamCloseAll();
      exit(EXIT_FAILURE);
    }
}


void cdoWarning(const char *fmt, ...)
{
  if ( _Verbose )
    {
      set_text_color(stderr, BRIGHT, YELLOW);
      fprintf(stderr, "%s (Warning): ", processInqPrompt());
      reset_text_color(stderr);

      set_text_color(stderr, RESET, BLACK);
      va_list args;
      va_start(args, fmt);
      vfprintf(stderr, fmt, args);
      va_end(args);
      reset_text_color(stderr);
      fprintf(stderr, "\n");
    }
}


void cdoPrint(const char *fmt, ...)
{
  if ( ! Options::silentMode )
    {
      set_text_color(stderr, RESET, GREEN);
      fprintf(stderr, "%s: ", processInqPrompt());
      reset_text_color(stderr);

      set_text_color(stderr, RESET, BLACK);
      va_list args;
      va_start(args, fmt);
      vfprintf(stderr, fmt, args);
      va_end(args);
      reset_text_color(stderr);
      fprintf(stderr, "\n");
    }
}


void cdoPrintBlue(const char *fmt, ...)
{
  if ( ! Options::silentMode )
    {
      set_text_color(stderr, RESET, GREEN);
      fprintf(stderr, "%s: ", processInqPrompt());
      reset_text_color(stderr);

      set_text_color(stderr, RESET, BLUE);
      va_list args;
      va_start(args, fmt);
      vfprintf(stderr, fmt, args);
      va_end(args);
      reset_text_color(stderr);
      fprintf(stderr, "\n");
    }
}


void cdoPrintRed(const char *fmt, ...)
{
  if ( ! Options::silentMode )
    {
      set_text_color(stderr, RESET, GREEN);
      fprintf(stderr, "%s: ", processInqPrompt());
      reset_text_color(stderr);

      set_text_color(stderr, RESET, RED);
      va_list args;
      va_start(args, fmt);
      vfprintf(stderr, fmt, args);
      va_end(args);
      reset_text_color(stderr);
      fprintf(stderr, "\n");
    }
}
