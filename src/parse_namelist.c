#include <errno.h>
#include "cdo_int.h"
#include "pmlist.h"
#include "namelist.h"


static
void kvlist_append_namelist(list_t *kvl, const char *key, const char *buffer, namelisttok_t *t, int nvalues)
{
  keyValues_t *keyval = (keyValues_t *) malloc(sizeof(keyValues_t));
  keyval->key = strdup(key);
  keyval->nvalues = nvalues;
  keyval->values = NULL;
  if ( nvalues > 0 )
    {
      keyval->values = (char **) malloc(nvalues*sizeof(char*));
      for ( int i = 0; i < nvalues; ++i )
        {
          size_t len = t[i].end - t[i].start;
          char *value = (char*) malloc((len+1)*sizeof(char));
          //printf(" value >%.*s<\n", len, buffer+t[i].start);
          snprintf(value, len+1, "%.*s", (int)len, buffer+t[i].start);
          value[len] = 0;
          keyval->values[i] = value;
        }
      list_append(kvl, &keyval);
    }
}

static
int get_number_of_values(int ntok, namelisttok_t *tokens)
{
  int it;
  
  for ( it = 0; it < ntok; ++it )
    {
      namelisttok_t *t = &tokens[it];
      if ( t->type != NAMELIST_WORD && t->type != NAMELIST_STRING ) break;
    }

  if ( it == ntok ) it = 0;
  
  return it;
}

static
int namelist_to_pml(list_t *pml, namelist_parser *parser, char *buf)
{
  char name[4096];
  list_t *kvl = NULL;
  namelisttok_t *t;
  namelisttok_t *tokens = parser->tokens;
  unsigned int ntok = parser->toknext;
  // printf("Number of tokens %d\n", ntok);

  for ( unsigned int it = 0; it < ntok; ++it )
    {
      t = &tokens[it];
      // printf("Token %u", it+1);
      if ( t->type == NAMELIST_OBJECT )
        {
          name[0] = 0;
          if ( it+1 < ntok && tokens[it+1].type == NAMELIST_WORD )
            {
              it++;
              t = &tokens[it];
              snprintf(name, sizeof(name), "%.*s", t->end - t->start, buf+t->start);
              name[sizeof(name)-1] = 0;
            }
          kvl = kvlist_new(name);
          list_append(pml, &kvl);
        }
      else if ( t->type == NAMELIST_KEY )
        {
          // printf(" key >%.*s<\n", t->end - t->start, buf+t->start);
          snprintf(name, sizeof(name), "%.*s", t->end - t->start, buf+t->start);
          name[sizeof(name)-1] = 0;
          if ( kvl == NULL ) printf("kvl not defined\n");
          int nvalues = get_number_of_values(ntok-it, &tokens[it+1]);
          kvlist_append_namelist(kvl, name, buf, &tokens[it+1], nvalues);
          it += nvalues;
        }
      else
        {
          // printf(" token >%.*s<\n", t->end - t->start, buf+t->start);
          break;
        }
    }

  return 0;
}


void parse_namelist_buffer_to_pml(list_t *pml, const char *filename, size_t buffersize, char *buffer)
{
  namelist_parser *p = namelist_new();

  int status = namelist_parse(p, buffer, buffersize);
  if ( status )
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

  // namelist_dump(p, buffer);

  namelist_to_pml(pml, p, buffer);

  namelist_destroy(p);
}



list_t *cdo_parse_namelist(const char *filename)
{
  assert(filename != NULL);

  size_t filesize = fileSize(filename);
  if ( filesize == 0 )
    {
      fprintf(stderr, "Empty table file: %s\n", filename);
      return NULL;
    }

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

  parse_namelist_buffer_to_pml(pml, filename, filesize, buffer);

  Free(buffer);
  
  return pml;
}
