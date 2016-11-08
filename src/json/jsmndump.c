#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include "jsmn.h"


/*
 * An example of reading JSON from stdin and printing its content to stdout.
 * The output looks like YAML, but I'm not sure if it's really compatible.
 */
static
int dump(const char *js, jsmntok_t *t, size_t count, int indent)
{
  int i, j, k;
  if (count == 0)  return 0;

  if (t->type == JSMN_PRIMITIVE)
    {
      printf("%.*s", t->end - t->start, js+t->start);
      return 1;
    }
  else if (t->type == JSMN_STRING)
    {
      printf("'%.*s'", t->end - t->start, js+t->start);
      return 1;
    }
  else if (t->type == JSMN_OBJECT)
    {
      printf("\n");
      printf("Object: size %d\n", t->size);
      j = 0;
      for (i = 0; i < t->size; i++)
        {
          for (k = 0; k < indent; k++) printf("  ");
          j += dump(js, t+1+j, count-j, indent+1);
          printf(": ");
          j += dump(js, t+1+j, count-j, indent+1);
          printf("\n");
        }
      return j+1;
    }
  else if (t->type == JSMN_ARRAY)
    {
      j = 0;
      printf("\n");
      for (i = 0; i < t->size; i++)
        {
          for (k = 0; k < indent-1; k++) printf("  ");
          printf("   - ");
          j += dump(js, t+1+j, count-j, indent+1);
          printf("\n");
        }
      return j+1;
    }
  return 0;
}


int main()
{
  FILE *fp = stdin;
  int filedes = fileno(fp);
  struct stat buf;
  size_t filesize = 0;
  if ( fstat(filedes, &buf) == 0 ) filesize = (size_t) buf.st_size;

  if ( filesize == 0 )
    {
      fprintf(stderr, "Empty stream!\n");
      return -1;
    }

  char *buffer = (char*) malloc(filesize);
  size_t nitems = fread(buffer, 1, filesize, fp);

  if ( nitems != filesize )
    {
      free(buffer);
      fprintf(stderr, "Read failed on stdin!\n");
      return -1;
    }

  /* Prepare parser */
  jsmn_parser *p = jsmn_new();

  int status = jsmn_parse(p, buffer, filesize);
  if ( status != 0 )
    {
      const char *filename = "stdin";
      switch (status)
        {
        case JSMN_ERROR_INVAL: fprintf(stderr, "JSON error: Invalid character in %s (line=%d character='%c')!\n", filename, p->lineno, buffer[p->pos]); break;
        case JSMN_ERROR_PART:  fprintf(stderr, "JSON error: End of string not found in %s (line=%d)!\n", filename, p->lineno); break;
        default:               fprintf(stderr, "JSON error in %s (line=%d)\n", filename, p->lineno); break;
        }
    }

  dump(buffer, p->tokens, p->toknext, 0);
  free(buffer);

  jsmn_destroy(p);

  return EXIT_SUCCESS;
}
