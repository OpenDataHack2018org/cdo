#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>

#include "namelist.h"


int main(int argc, char *argv[])
{
  if ( argc != 2 )
    {
      fprintf(stderr, "Usage: %s namelist\n", argv[0]);
      return -1;
    }

  const char *filename = argv[1];
  printf("Parse namelist %s:\n", filename);

  struct stat sbuf;
  size_t filesize = (stat(filename, &sbuf) == 0) ? sbuf.st_size : 0;

  if ( filesize == 0 )
    {
      fprintf(stderr, "Empty table file: %s\n", filename);
      return -1;
    }

  FILE *fp = fopen(filename, "r");
  if ( fp == NULL )
    {
      fprintf(stderr, "Open failed on %s: %s\n", filename, strerror(errno));
      return -1;
    }

  char *buffer = (char*) malloc(filesize);
  size_t nitems = fread(buffer, 1, filesize, fp);

  fclose(fp);

  if ( nitems != filesize )
    {
      fprintf(stderr, "Read failed on %s!\n", filename);
      return -1;
    }

  namelist_parser *p = namelist_new();

  int status = namelist_parse(p, buffer, filesize);
  printf("Processed number of lines: %d\n", p->lineno-1);
  if ( status != 0 )
    {
      switch (status)
        {
        case NAMELIST_ERROR_INVAL: fprintf(stderr, "Namelist error: Invalid character in %s (line=%d character='%c')!\n", filename, p->lineno, buffer[p->pos]); break;
        case NAMELIST_ERROR_PART:  fprintf(stderr, "Namelist error: End of string not found in %s (line=%d)!\n", filename, p->lineno); break;
        case NAMELIST_ERROR_INKEY: fprintf(stderr, "Namelist error: Invalid key word in %s (line=%d)!\n", filename, p->lineno); break;
        case NAMELIST_ERROR_INTYP: fprintf(stderr, "Namelist error: Invalid key word type in %s (line=%d)!\n", filename, p->lineno); break;
        case NAMELIST_ERROR_INOBJ: fprintf(stderr, "Namelist error: Invalid object in %s (line=%d)!\n", filename, p->lineno); break;
        case NAMELIST_ERROR_EMKEY: fprintf(stderr, "Namelsit error: Emtry key name in %s (line=%d)!\n", filename, p->lineno); break;
        default:                   fprintf(stderr, "Namelsit error in %s (line=%d)!\n", filename, p->lineno); break;
        }
    }

  namelist_dump(p, buffer);

  namelist_destroy(p);

  free(buffer);

  return 0;
}
