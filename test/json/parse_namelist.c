#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>

typedef enum {
  NAMELIST_UNDEFINED = 0,
  NAMELIST_OBJECT    = 1,
  NAMELIST_KEY       = 2,
  NAMELIST_STRING    = 3,
  NAMELIST_WORD      = 4
} namelisttype_t;

/**
 * NAMELIST token description.
 * type         type (object, array, string etc.)
 * start        start position in NAMELIST buffer
 * end          end position in NAMELIST buffer
 */
typedef struct {
  namelisttype_t type;
  int start;
  int end;
  int size;
} namelisttok_t;


typedef struct {
  unsigned int pos;
  unsigned int toknext;
  int toksuper;
} namelist_parser;


static
void namelist_init(namelist_parser *parser)
{
  parser->pos = 0;
  parser->toknext = 0;
  parser->toksuper = -1;
}


namelist_parser *namelist_new(void)
{
  namelist_parser *parser = (namelist_parser *) malloc(sizeof(namelist_parser));

  namelist_init(parser);

  return parser;
}


void namelist_destroy(namelist_parser *parser)
{
  if ( parser )
    {
      free(parser);
    }
}


int namelist_parse(namelist_parser *parser, const char *buf, size_t len)
{
  int status = 0;
  int lineno = 0;
  int count = parser->toknext;
  namelisttok_t *token;

  for ( ; parser->pos < len && buf[parser->pos] != '\0'; parser->pos++ )
    {
      namelisttype_t type;

      char c = buf[parser->pos];
      switch (c)
        {
        case '&':
          break;
        case '/':
          break;
        case '\t': case ' ':
          break;
        case '\r':
          lineno++;
          if ( parser->pos+1 < len && buf[parser->pos+1] == '\n' ) parser->pos++;
          break;
        case '\n':
          lineno++;
          break;
        case '#': case '!': // Skip to end of line
          for (; parser->pos < len && buf[parser->pos] != '\0'; parser->pos++)
            if ( buf[parser->pos] == '\r' || buf[parser->pos] == '\n' ) break;
        case '\"':
          break;
        }
    }

  printf("Processed number of lines: %d\n", lineno);
  return status;
}


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
  namelist_parse(p, buffer, filesize);
  /*
  namelist_dump(p, buffer, filesize);
  */
  namelist_destroy(p);

  free(buffer);

  return 0;
}
