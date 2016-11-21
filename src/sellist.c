#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
      if ( nvalues ) sellist->entry[i].flag = (bool*) calloc(nvalues, sizeof(bool));
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
        }

      free(sellist);
    }
}


int sellist_add(sellist_t *sellist, const char *txt, const char *name, int type)
{
  int entry = -1;

  if ( sellist )
    {
      for ( int i = 0; i < sellist->size; ++i )
        {
          const char *key = sellist->entry[i].key;
          if ( strcmp(key, name) == 0 )
            {
              sellist->entry[i].type = type;
              sellist->entry[i].txt = strdup(txt);
              entry = i;
              break;
            }
        }
    }

  return entry;
}


int sellist_num_par(sellist_t *sellist, int entry)
{
  int num_par = 0;

  if ( sellist && entry >= 0 && entry < sellist->size ) num_par = sellist->entry[entry].nvalues;

  return num_par;
}
