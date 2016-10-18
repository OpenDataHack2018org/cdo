#include <errno.h>
#include "cdo_int.h"
#include "pmlist.h"
#include "json/jsmn.h"


static
char *readLineFromBuffer(char *buffer, size_t *buffersize, char *line, size_t len)
{
  int ichar;
  size_t ipos = 0;

  while ( *buffersize )
    {
      ichar = *buffer;
      (*buffersize)--;
      buffer++;
      if ( ichar == '\r' )
        {
          if ( *buffersize )
            {
              ichar = *buffer;
              if ( ichar == '\n' )
                {
                  (*buffersize)--;
                  buffer++;
                }
            }
          break;
        }
      if ( ichar == '\n' ) break;
      line[ipos++] = ichar;
      if ( ipos >= len )
        {
          fprintf(stderr, "readLineFromBuffer: end of line not found (maxlen = %ld)!\n", len);
          break;
        }
    }
  line[ipos] = 0;

  if ( *buffersize == 0 && ipos == 0 ) buffer = NULL;

  return buffer;
}

static
char *skipSeparator(char *pline)
{
  while ( isspace((int) *pline) ) pline++;
  if ( *pline == '=' || *pline == ':' ) pline++;
  while ( isspace((int) *pline) ) pline++;

  return pline;
}

static
char *getElementName(char *pline, char *name)
{
  while ( isspace((int) *pline) ) pline++;
  size_t len = strlen(pline);
  size_t pos = 0;
  while ( pos < len && !isspace((int) *(pline+pos)) && *(pline+pos) != '=' && *(pline+pos) != ':' ) pos++;

  strncpy(name, pline, pos);
  name[pos] = 0;

  pline += pos;
  return pline;
}

static
char *getElementValue(char *pline)
{
  while ( isspace((int) *pline) ) pline++;
  size_t len = strlen(pline);
  while ( isspace((int) *(pline+len-1)) && len ) { *(pline+len-1) = 0; len--;}

  return pline;
}


void parse_buffer_to_pml(list_t *pml, size_t buffersize, char *buffer)
{
  char line[4096];
  char name[256];
  char *pline;
  char listkey1[] = "axis_entry:";
  char listkey2[] = "variable_entry:";
  int linenumber = 0;
  int listtype = 0;
  list_t *kvl = NULL;

  while ( (buffer = readLineFromBuffer(buffer, &buffersize, line, sizeof(line))) )
    {
      linenumber++;
      pline = line;
      while ( isspace((int) *pline) ) pline++;
      if ( *pline == '#' || *pline == '!' || *pline == '\0' ) continue;
      //  len = (int) strlen(pline);
      if ( listtype == 0 && *pline == '&' ) listtype = 1;
      
      if ( strncmp(pline, listkey1, strlen(listkey1)) == 0 )
	{
	  pline += strlen(listkey1);

	  listtype = 2;

          kvl = list_new(sizeof(keyValues_t *), free_keyval, "axis");
          list_append(pml, &kvl);

	  pline = skipSeparator(pline);
	  pline = getElementValue(pline);

	  if ( *pline ) kvlist_append(kvl, "name", (const char **)&pline, 1);
	}
      else if ( strncmp(pline, listkey2, strlen(listkey2)) == 0 )
	{
	  pline += strlen(listkey2);

	  listtype = 2;

          kvl = list_new(sizeof(keyValues_t *), free_keyval, "variable");
          list_append(pml, &kvl);

	  pline = skipSeparator(pline);
	  pline = getElementValue(pline);

	  if ( *pline ) kvlist_append(kvl, "name", (const char **)&pline, 1);
	}
      else
	{
	  pline = getElementName(pline, name);
	  pline = skipSeparator(pline);
	  pline = getElementValue(pline);

	  if ( kvl == NULL )
            {
              kvl = list_new(sizeof(keyValues_t *), free_keyval, "global");
              list_append(pml, &kvl);
            }

	  if ( *pline ) kvlist_append(kvl, name, (const char **)&pline, 1);
	}
    }
}

static
int dump_json(const char *js, jsmntok_t *t, size_t count, int level)
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
      //  printf("Object: size %d\n", t->size);
      printf("Object: size %d count %d level %d\n", t->size, (int)count, level);
      j = 0;
      for (i = 0; i < t->size; i++)
        {
          for (k = 0; k < level; k++) printf("  ");
          j += dump_json(js, t+1+j, count-j, level+1);
          printf(": ");
          j += dump_json(js, t+1+j, count-j, level+1);
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
          for (k = 0; k < level-1; k++) printf("  ");
          printf("   - ");
          j += dump_json(js, t+1+j, count-j, level+1);
          printf("\n");
        }
      return j+1;
    }
  return 0;
}

static
int json_to_pml(list_t *pml, const char *js, jsmntok_t *t, size_t count, int level)
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
      //  printf("Object: size %d\n", t->size);
      printf("Object: size %d count %d level %d\n", t->size, (int)count, level);
      j = 0;
      for (i = 0; i < t->size; i++)
        {
          for (k = 0; k < level; k++) printf("  ");
          j += json_to_pml(pml, js, t+1+j, count-j, level+1);
          printf(": ");
          j += json_to_pml(pml, js, t+1+j, count-j, level+1);
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
          for (k = 0; k < level-1; k++) printf("  ");
          printf("   - ");
          j += json_to_pml(pml, js, t+1+j, count-j, level+1);
          printf("\n");
        }
      return j+1;
    }
  return 0;
}


void parse_json_buffer_to_pml(list_t *pml, size_t buffersize, char *buffer)
{
  int r;
  char *js = buffer;
  size_t jslen = buffersize;
        
  jsmn_parser p;
  jsmntok_t *tok;
  size_t tokcount = 2;

  /* Prepare parser */
  jsmn_init(&p);

  /* Allocate some tokens as a start */
  tok = Malloc(sizeof(*tok) * tokcount);

 again:
  r = jsmn_parse(&p, js, jslen, tok, tokcount);
  if ( r < 0 )
    {
      if ( r == JSMN_ERROR_NOMEM )
        {
          tokcount = tokcount * 2;
          tok = Realloc(tok, sizeof(*tok) * tokcount);
          goto again;
        }
    }

  dump_json(js, tok, p.toknext, 0);
}


list_t *cdo_parse_cmor_file(const char *filename)
{
  assert(filename != NULL);

  size_t filesize = fileSize(filename);

  FILE *fp = fopen(filename, "r");
  if ( fp == NULL )
    {
      fprintf(stderr, "Open failed on %s: %s\n", filename, strerror(errno));
      return NULL;
    }

  char *buffer = (char*) Malloc(filesize);
  size_t nitems = fread(buffer, 1, filesize, fp);

  fclose(fp);

  if ( nitems != filesize )
    {
      fprintf(stderr, "Read failed on %s!\n", filename);
      return NULL;
    }
 
  list_t *pml = list_new(sizeof(list_t *), free_kvlist, filename);

  if ( buffer[0] == '{' )
    parse_json_buffer_to_pml(pml, filesize, buffer);
  else
    parse_buffer_to_pml(pml, filesize, buffer);
  
  return pml;
}
