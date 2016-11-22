#include <cdo_int.h>
#include "sellist.h"


sellist_t *sellist_create(list_t *kvlist)
{
  printf("kvlist size %d\n", list_size(kvlist));
  sellist_t *sellist = (sellist_t *) malloc(sizeof(sellist_t));
  sellist->size = list_size(kvlist);
  sellist->entry = (selentry_t *) malloc(sellist->size*sizeof(selentry_t));

  int i = 0;
  for ( listNode_t *kvnode = kvlist->head; kvnode; kvnode = kvnode->next )
    {
      keyValues_t *kv = *(keyValues_t **)kvnode->data;
      sellist->entry[i].key = kv->key;
      sellist->entry[i].values = kv->values;
      sellist->entry[i].nvalues = kv->nvalues;
      const char *key = kv->key;
      char **values = kv->values;
      int nvalues = kv->nvalues;
      printf("%s =", key);
      for ( int i = 0; i < nvalues; ++i ) printf(" '%s'", values[i]);
      printf("\n");
      ++i;
    }

  for ( int i = 0; i < sellist->size; ++i )
    {
      const char *key = sellist->entry[i].key;
      char **values = sellist->entry[i].values;
      int nvalues = sellist->entry[i].nvalues;
      sellist->entry[i].flag = NULL;
      sellist->entry[i].cvalues = NULL;
      if ( nvalues )
        {
          sellist->entry[i].flag = (bool*) calloc(nvalues, sizeof(bool));
          sellist->entry[i].cvalues = (cvalues_t*) calloc(nvalues, sizeof(cvalues_t));
        }
      printf("%s =", key);
      for ( int i = 0; i < nvalues; ++i ) printf(" '%s'", values[i]);
      printf("\n");
    }

  return sellist;
}


void sellist_destroy(sellist_t *sellist)
{
  if ( sellist )
    {
      for ( int i = 0; i < sellist->size; ++i )
        {
          const char *key = sellist->entry[i].key;
          char **values = sellist->entry[i].values;
          int nvalues = sellist->entry[i].nvalues;
          sellist->entry[i].flag = NULL;
          if ( sellist->entry[i].flag ) free(sellist->entry[i].flag);
          if ( sellist->entry[i].txt ) free(sellist->entry[i].txt);
          if ( sellist->entry[i].cvalues ) free(sellist->entry[i].cvalues);
        }

      free(sellist);
    }
}


int sellist_add(sellist_t *sellist, const char *txt, const char *name, int type)
{
  int idx = -1;

  if ( sellist )
    {
      for ( int i = 0; i < sellist->size; ++i )
        {
          const char *key = sellist->entry[i].key;
          if ( strcmp(key, name) == 0 )
            {
              idx = i;
              break;
            }
        }

      if ( idx >= 0 && idx < sellist->size )
        {
          selentry_t *e = &(sellist->entry[idx]);
          e->type = type;
          e->txt = strdup(txt);
          int nvalues = e->nvalues;
          for ( int i = 0; i < nvalues; ++i )
            switch (type)
              {
              case SELLIST_INT:  e->cvalues[i].ival = parameter2int(e->values[i]); break;
              case SELLIST_FLT:  e->cvalues[i].dval = parameter2double(e->values[i]); break;
              case SELLIST_WORD: e->cvalues[i].cval = parameter2word(e->values[i]); break;
              }
          printf("add %s idx=%d\n", e->key, idx);
          printf("%s =", e->key);
          for ( int i = 0; i < nvalues; ++i )
            switch (type)
              {
              case SELLIST_INT:  printf(" %d", e->cvalues[i].ival); break;
              case SELLIST_FLT:  printf(" %g", e->cvalues[i].dval); break;
              case SELLIST_WORD: printf(" %s", e->cvalues[i].cval); break;
              }
          printf("\n");
        }
    }

  if ( idx >= 0 )
    {
      selentry_t *e = &(sellist->entry[idx]);
      printf("add %s idx=%d nvalues=%d\n", e->key, idx, e->nvalues);
    }
  return idx;
}


int sellist_nvalues(sellist_t *sellist, int idx)
{
  int nvalues = 0;

  if ( sellist && idx >= 0 && idx < sellist->size ) nvalues = sellist->entry[idx].nvalues;

  return nvalues;
}


bool sellist_check(sellist_t *sellist, int idx, void *par)
{
  bool found = false;

  if ( idx < 0 || idx >= sellist->size ) return found;

  int nvalues = sellist_nvalues(sellist, idx);

  if ( nvalues )
    {
      selentry_t *e = &(sellist->entry[idx]);
      printf("check for %s idx=%d nvalues=%d\n", e->key, idx, nvalues);
      int type = e->type;
      for ( int i = 0; i < nvalues; ++i )
        {
          switch (type)
            {
            case SELLIST_INT:  if ( *(int*)par == e->cvalues[i].ival )                  found = true; break;
            case SELLIST_FLT:  if ( fabs(*(double*)par - e->cvalues[i].dval) < 1.e-4 )  found = true; break;
            case SELLIST_WORD: if ( wildcardmatch(e->cvalues[i].cval, *(char**)par) )   found = true; break;
            }

          if ( found ) e->flag[i] = true;
        }
    }

  return found;
}


void sellist_check_flag(sellist_t *sellist, int idx)
{
  if ( idx < 0 || idx >= sellist->size ) return;

  int nvalues = sellist_nvalues(sellist, idx);

  if ( nvalues )
    {
      selentry_t *e = &(sellist->entry[idx]);
      printf("check flag for %s idx=%d\n", e->key, idx);
      for ( int i = 0; i < nvalues; ++i )
        if ( e->flag[i] == false ) cdoWarning("%s >%s< not found!", e->txt, e->values[i]);
    }
}

