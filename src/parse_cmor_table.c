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

          kvl = list_new(sizeof(keyValues_t *), free_keyval, listkey1);
          list_append(pml, &kvl);

	  pline = skipSeparator(pline);
	  pline = getElementValue(pline);

	  if ( *pline ) kvlist_append(kvl, "name", (const char **)&pline, 1);
	}
      else if ( strncmp(pline, listkey2, strlen(listkey2)) == 0 )
	{
	  pline += strlen(listkey2);

	  listtype = 2;

          kvl = list_new(sizeof(keyValues_t *), free_keyval, listkey2);
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
              kvl = list_new(sizeof(keyValues_t *), free_keyval, "Header");
              list_append(pml, &kvl);
            }

	  if ( *pline ) kvlist_append(kvl, name, (const char **)&pline, 1);
	}
    }
}

// not used
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
void kvlist_append_json(list_t *kvl, const char *key, const char *js, jsmntok_t *t, int nvalues)
{
  char name[1024];
  keyValues_t *keyval = (keyValues_t *) malloc(sizeof(keyValues_t));
  keyval->key = strdup(key);
  keyval->nvalues = nvalues;
  keyval->values = (char **) malloc(nvalues*sizeof(char*));
  for ( int i = 0; i < nvalues; ++i )
    {
      sprintf(name, "%.*s", t[i].end - t[i].start, js+t[i].start);
      // printf("set %s: '%s'\n", key, name);
      keyval->values[i] = strdup(name);
    }
  list_append(kvl, &keyval);
}

static
int json_to_pml(list_t *pml, const char *js, jsmntok_t *t, int count)
{
  bool debug = false;
  char name[1024];
  list_t *kvl = NULL;
  int pmlname = -1;
  int i = 0;
  if ( t[0].type == JSMN_OBJECT )
    for ( int ib = 0; ib < t[0].size; ++ib )
      {
        ++i;
        pmlname = i;
        if ( debug ) printf("  object: %.*s\n", t[i].end - t[i].start, js+t[i].start);
        ++i;
        if ( t[i].type == JSMN_OBJECT )
          {
            int ic = 0;
          NEXT:
            sprintf(name, "%.*s", t[pmlname].end - t[pmlname].start, js+t[pmlname].start);
            // printf("new object: %s\n", name);
            kvl = list_new(sizeof(keyValues_t *), free_keyval, name);
            list_append(pml, &kvl);
                
            if ( t[i+2].type == JSMN_OBJECT )
              {
                if ( ic == 0 ) ic = t[i].size;
                else           ic--;
                
                ++i;
                kvlist_append_json(kvl, "name", js, &t[i], 1);
                if ( debug ) printf("    name: '%.*s'\n", t[i].end - t[i].start, js+t[i].start);
                ++i;
              }
            int n = t[i].size;
            for ( int jb = 0; jb < n; ++jb )
              {
                ++i;
                sprintf(name, "%.*s", t[i].end - t[i].start, js+t[i].start);
                if ( debug ) printf("    %.*s:", t[i].end - t[i].start, js+t[i].start);
                ++i;
                if ( t[i].type == JSMN_ARRAY )
                  {
                    int nk = t[i].size;
                    kvlist_append_json(kvl, name, js, &t[i+1], nk);
                    for ( int k = 0; k < nk; ++k )
                      {
                        ++i;
                        if ( debug ) printf(" '%.*s'", t[i].end - t[i].start, js+t[i].start);
                      }
                  }
                else
                  {
                    kvlist_append_json(kvl, name, js, &t[i], 1);
                    if ( debug ) printf(" '%.*s'", t[i].end - t[i].start, js+t[i].start);
                  }
                if ( debug ) printf("\n");
              }
            if ( ic > 1 ) goto NEXT;
          }
      }

  if ( debug ) printf("Processed %d of %d tokens!\n", i, count-1);

  return 0;
}


void parse_json_buffer_to_pml(list_t *pml, size_t buffersize, char *buffer)
{
  char *js = buffer;
  size_t jslen = buffersize;

  /* Prepare parser */
  jsmn_parser p;
  jsmn_init(&p);

  /* Allocate some tokens as a start */
  size_t tokcount = 2;
  jsmntok_t *tok = Malloc(sizeof(*tok) * tokcount);

  int r;
 AGAIN:
  r = jsmn_parse(&p, js, jslen, tok, tokcount);
  if ( r < 0 )
    {
      if ( r == JSMN_ERROR_NOMEM )
        {
          tokcount = tokcount * 2;
          tok = Realloc(tok, sizeof(*tok) * tokcount);
          goto AGAIN;
        }
    }

  json_to_pml(pml, js, tok, (int)p.toknext);

  Free(tok);
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

  Free(buffer);
  
  return pml;
}
