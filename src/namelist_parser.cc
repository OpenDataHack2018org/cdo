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
#include "cdo_int.h"
#include "namelist.h"

static void
kvlist_append_namelist(list_t *kvlist, const char *key, const char *buffer, NamelistToken *t, int nvalues)
{
  char vbuf[4096];
  keyValues_t *keyval = (keyValues_t *) Malloc(sizeof(keyValues_t));
  keyval->key = strdup(key);
  keyval->nvalues = nvalues;
  keyval->values = NULL;

  if (nvalues > 0) keyval->values = (char **) Malloc(nvalues * sizeof(char *));
  for (int i = 0; i < nvalues; ++i)
    {
      size_t len = t[i].end - t[i].start;
      char *value = (char *) Malloc((len + 1) * sizeof(char));
      // printf(" value[%d] >%.*s<\n", i, (int)len, buffer+t[i].start);
      const char *pval = buffer + t[i].start;
      if (len < sizeof(vbuf))  // snprintf seems to call strlen(pval)
        {
          memcpy(vbuf, buffer + t[i].start, len);
          vbuf[len] = 0;
          pval = vbuf;
        }
      snprintf(value, len + 1, "%.*s", (int) len, pval);
      value[len] = 0;
      keyval->values[i] = value;
    }

  list_append(kvlist, &keyval);
}

static int
get_number_of_values(int ntok, NamelistToken *tokens)
{
  int it;

  for (it = 0; it < ntok; ++it)
    {
      NamelistToken *t = &tokens[it];
      if (t->type != NamelistType::WORD && t->type != NamelistType::STRING) break;
    }

  return it;
}

static int
namelist_to_pml(list_t *pmlist, NamelistParser &parser, char *buf)
{
  char name[4096];
  list_t *kvlist = NULL;
  NamelistToken *t;
  NamelistToken *tokens = parser.tokens.data();
  unsigned int ntok = parser.toknext;
  // printf("Number of tokens %d\n", ntok);

  for (unsigned int it = 0; it < ntok; ++it)
    {
      t = &tokens[it];
      // printf("Token %u", it+1);
      if (t->type == NamelistType::OBJECT)
        {
          name[0] = 0;
          if (it + 1 < ntok && tokens[it + 1].type == NamelistType::WORD)
            {
              it++;
              t = &tokens[it];
              snprintf(name, sizeof(name), "%.*s", t->end - t->start, buf + t->start);
              name[sizeof(name) - 1] = 0;
            }
          kvlist = kvlist_new(name);
          list_append(pmlist, &kvlist);
        }
      else if (t->type == NamelistType::KEY)
        {
          if (kvlist == NULL)
            {
              kvlist = kvlist_new(NULL);
              list_append(pmlist, &kvlist);
            }
          // printf(" key >%.*s<\n", t->end - t->start, buf+t->start);
          snprintf(name, sizeof(name), "%.*s", t->end - t->start, buf + t->start);
          name[sizeof(name) - 1] = 0;
          int nvalues = get_number_of_values(ntok - it - 1, &tokens[it + 1]);
          // printf("nvalues %d\n", nvalues);
          kvlist_append_namelist(kvlist, name, buf, &tokens[it + 1], nvalues);
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

list_t *
namelistbuf_to_pmlist(listbuf_t *listbuf)
{
  const char *name = listbuf->name;
  NamelistParser p;

  NamelistError status = p.parse(listbuf->buffer, listbuf->size);
  if (status!=NamelistError::UNDEFINED)
    {
      switch (status)
        {
        case NamelistError::INVAL:
          fprintf(stderr, "Namelist error: Invalid character in %s (line=%d character='%c')!\n", name, p.lineno,
                  listbuf->buffer[p.pos]);
          break;
        case NamelistError::PART:
          fprintf(stderr, "Namelist error: End of string not found in %s (line=%d)!\n", name, p.lineno);
          break;
        case NamelistError::INKEY: fprintf(stderr, "Namelist error: Invalid key word in %s (line=%d)!\n", name, p.lineno); break;
        case NamelistError::INTYP:
          fprintf(stderr, "Namelist error: Invalid key word type in %s (line=%d)!\n", name, p.lineno);
          break;
        case NamelistError::INOBJ: fprintf(stderr, "Namelist error: Invalid object in %s (line=%d)!\n", name, p.lineno); break;
        case NamelistError::EMKEY: fprintf(stderr, "Namelist error: Emtry key name in %s (line=%d)!\n", name, p.lineno); break;
        default: fprintf(stderr, "Namelist error in %s (line=%d)!\n", name, p.lineno); break;
        }
      cdoAbort("Namelist error!");
    }

  // p.dump(listbuf->buffer);
  if (p.verify())
    {
      fprintf(stderr, "Namelist error: Invalid contents in %s!\n", name);
      cdoAbort("Namelist error!");
    }

  list_t *pmlist = list_new(sizeof(list_t *), free_kvlist, listbuf->name);

  namelist_to_pml(pmlist, p, listbuf->buffer);

  return pmlist;
}

list_t *
namelist_to_pmlist(FILE *fp, const char *name)
{
  listbuf_t *listbuf = listbuf_new();
  if (listbuf_read(listbuf, fp, name)) cdoAbort("Read error on namelist %s!", name);

  list_t *pmlist = namelistbuf_to_pmlist(listbuf);
  if (pmlist == NULL) cdoAbort("Namelist not found!");

  listbuf_destroy(listbuf);

  return pmlist;
}
