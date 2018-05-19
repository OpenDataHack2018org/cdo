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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>

#include "namelist.h"

int
main(int argc, char *argv[])
{
  if (argc != 2)
    {
      fprintf(stderr, "Usage: %s namelist\n", argv[0]);
      return -1;
    }

  const char *filename = argv[1];
  printf("Parse namelist %s:\n", filename);

  struct stat sbuf;
  size_t filesize = (stat(filename, &sbuf) == 0) ? sbuf.st_size : 0;

  if (filesize == 0)
    {
      fprintf(stderr, "Empty table file: %s\n", filename);
      return -1;
    }

  FILE *fp = fopen(filename, "r");
  if (fp == NULL)
    {
      fprintf(stderr, "Open failed on %s: %s\n", filename, strerror(errno));
      return -1;
    }

  char *buffer = (char *) malloc(filesize);
  size_t nitems = fread(buffer, 1, filesize, fp);

  fclose(fp);

  if (nitems != filesize)
    {
      fprintf(stderr, "Read failed on %s!\n", filename);
      return -1;
    }

  NamelistParser p;

  NamelistError status = p.parse(buffer, filesize);
  printf("Processed number of lines: %d\n", p.lineno - 1);
  if (status != 0)
    {
      switch (status)
        {
        case NamelistError::INVAL:
          fprintf(stderr, "Namelist error: Invalid character in %s (line=%d character='%c')!\n", filename, p.lineno,
                  buffer[p.pos]);
          break;
        case NamelistError::PART:
          fprintf(stderr, "Namelist error: End of string not found in %s (line=%d)!\n", filename, p.lineno);
          break;
        case NamelistError::INKEY:
          fprintf(stderr, "Namelist error: Invalid key word in %s (line=%d)!\n", filename, p.lineno);
          break;
        case NamelistError::INTYP:
          fprintf(stderr, "Namelist error: Invalid key word type in %s (line=%d)!\n", filename, p.lineno);
          break;
        case NamelistError::INOBJ: fprintf(stderr, "Namelist error: Invalid object in %s (line=%d)!\n", filename, p.lineno); break;
        case NamelistError::EMKEY: fprintf(stderr, "Namelsit error: Emtry key name in %s (line=%d)!\n", filename, p.lineno); break;
        default: fprintf(stderr, "Namelsit error in %s (line=%d)!\n", filename, p.lineno); break;
        }
    }

  p.dump(buffer);

  free(buffer);

  return 0;
}
