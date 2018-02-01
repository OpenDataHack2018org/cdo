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
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "dmemory.h"
#include "listbuf.h"

listbuf_t *listbuf_new(void)
{
  listbuf_t *listbuf = (listbuf_t *) malloc(sizeof(listbuf_t));
  if ( listbuf )
    {
      listbuf->size = 0;
      listbuf->buffer = NULL;
      listbuf->name = NULL;
    }

  return listbuf;
}


int listbuf_read(listbuf_t *listbuf, FILE *fp, const char *name)
{
  int filedes = fileno(fp);
  struct stat buf;
  size_t filesize = 0;
  if ( fstat(filedes, &buf) == 0 ) filesize = (size_t) buf.st_size;

  if ( filesize == 0 )
    {
      fprintf(stderr, "%s: empty stream: %s\n", __func__, name);
      return -1;
    }

  char *buffer = (char*) Malloc(filesize);
  size_t nitems = fread(buffer, 1, filesize, fp);

  if ( nitems != filesize )
    {
      Free(buffer);
      fprintf(stderr, "%s: read failed on %s!\n", __func__, name);
      return -1;
    }

  listbuf->size = filesize;
  listbuf->buffer = buffer;

  if ( name ) listbuf->name = strdup(name);

  return 0;
}


void listbuf_destroy(listbuf_t *listbuf)
{
  if ( listbuf )
    {
      if ( listbuf->buffer ) free(listbuf->buffer);
      if ( listbuf->name ) free(listbuf->name);
      listbuf->size = 0;
      listbuf->buffer = NULL;
      listbuf->name = NULL;
      free(listbuf);
    }
}


