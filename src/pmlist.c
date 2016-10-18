#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "pmlist.h"

keyValues_t *kvlist_search(list_t *kvl, const char *key)
{
  if ( key )
    {
      listNode_t *node = kvl->head;
      while ( node )
        {
          keyValues_t *kv = *(keyValues_t **)node->data;
          if ( kv->key && *(kv->key) == *key && strcmp(kv->key, key) == 0 ) return kv;
          node = node->next;
        }
    }

  return NULL;
}


list_t *pml_search_kvl(list_t *pml, const char *key, const char *value)
{
  if ( pml && key && value )
    {
      listNode_t *node = pml->head;
      while ( node )
        {
          if ( node->data )
            {
              list_t *kvl = *(list_t **)node->data;
              keyValues_t *kv = kvlist_search(kvl, key);
              if ( kv && kv->nvalues > 0 && *(kv->values[0]) == *value && strcmp(kv->values[0], value) == 0 ) return kvl;
            }
          node = node->next;
        }
    }

  return NULL;
}


bool kvl_print_iter(void *data)
{
  keyValues_t *keyval = *(keyValues_t **)data;
  //  printf("Found %s list with %d keys: \n", list_name(keyval), list_size(keyval));
  printf("  key=%s  val0=%s\n", keyval->key, keyval->values[0]);
  return true;
}


bool pml_print_iter(void *data)
{
  list_t *kvl = *(list_t **)data;
  printf("Found %s list with %d keys: \n", list_name(kvl), list_size(kvl));
  list_for_each(kvl, kvl_print_iter);
  return true;
}


void free_keyval(void *data)
{
  keyValues_t *keyval = *(keyValues_t **)data;
  if ( keyval->key ) free(keyval->key);
  int nvalues = keyval->nvalues;
  for ( int i = 0; i < nvalues; ++i )
    if ( keyval->values[i] ) free(keyval->values[i]);
  free(keyval->values);
  free(keyval);
}


void free_kvlist(void *data)
{
  list_t *kvl = *(list_t **)data;
  //int n = list_size(kvl);
  list_destroy(kvl);
  //printf("Successfully freed %d keyvalues...\n", n);
}


void kvlist_append(list_t *kvl, const char *key, const char **values, int nvalues)
{
  keyValues_t *keyval = (keyValues_t *) malloc(sizeof(keyValues_t));
  keyval->key = strdup(key);
  keyval->nvalues = nvalues;
  keyval->values = (char **) malloc(nvalues*sizeof(char*));
  for ( int i = 0; i < nvalues; ++i ) keyval->values[i] = strdup(values[i]);
  list_append(kvl, &keyval);
}

/*
int main(void)
{
  int numLists = 0;
  printf("Generating list with lists of keyValues...\n");

  list_t *pml = list_new(sizeof(list_t *), free_kvlist, "parameter");

  {
    const char *k1 = "longname", *k1vals[] = {"surface temperature"};
    const char *k2 = "name",     *k2vals[] = {"temperature"};
    const char *k3 = "values",   *k3vals[] = {"273.15", "292.5", "301.4"};

    list_t *kvl = list_new(sizeof(keyValues_t *), free_keyval, "p1");

    kvlist_append(kvl, k1, k1vals, sizeof(k1vals)/sizeof(k1vals[0]));
    kvlist_append(kvl, k2, k2vals, sizeof(k2vals)/sizeof(k2vals[0]));
    kvlist_append(kvl, k3, k3vals, sizeof(k3vals)/sizeof(k3vals[0]));

    list_append(pml, &kvl);
    numLists++;
  }
  {
    const char *k1 = "longname", *k1vals[] = {"surface pressure"};
    const char *k2 = "name",     *k2vals[] = {"pressure"};
    const char *k3 = "values",   *k3vals[] = {"1000", "850", "500"};
    const char *k4 = "units",    *k4vals[] = {"hPa"};

    list_t *kvl = list_new(sizeof(keyValues_t *), free_keyval, "p1");

    kvlist_append(kvl, k1, k1vals, sizeof(k1vals)/sizeof(k1vals[0]));
    kvlist_append(kvl, k2, k2vals, sizeof(k2vals)/sizeof(k2vals[0]));
    kvlist_append(kvl, k3, k3vals, sizeof(k3vals)/sizeof(k3vals[0]));
    kvlist_append(kvl, k4, k4vals, sizeof(k4vals)/sizeof(k4vals[0]));

    list_append(pml, &kvl);
    numLists++;
  }
  {
    const char *k1 = "longname", *k1vals[] = {"Air Temperature"};
    const char *k2 = "name",     *k2vals[] = {"ta"};
    const char *k3 = "valid_max",*k3vals[] = {"336"};
    const char *k4 = "units",    *k4vals[] = {"K"};

    list_t *kvl = list_new(sizeof(keyValues_t *), free_keyval, "p1");

    kvlist_append(kvl, k1, k1vals, sizeof(k1vals)/sizeof(k1vals[0]));
    kvlist_append(kvl, k2, k2vals, sizeof(k2vals)/sizeof(k2vals[0]));
    kvlist_append(kvl, k3, k3vals, sizeof(k3vals)/sizeof(k3vals[0]));
    kvlist_append(kvl, k4, k4vals, sizeof(k4vals)/sizeof(k4vals[0]));

    list_append(pml, &kvl);
    numLists++;
  }
  
 
  printf("pmlist=%s n=%d:\n", list_name(pml), list_size(pml));
  list_for_each(pml, pml_print_iter);

  list_t *kvl = pml_search_kvl(pml, "name", "pressure");
  printf("kvlist of name=pressure:\n");
  if ( kvl ) list_for_each(kvl, kvl_print_iter);

  printf("\n");
  int n = list_size(kvl);
  for ( int i = 0; i < n; ++i )
    {
      keyValues_t *kv = (keyValues_t *) list_entry(kvl, i);
      if ( kv ) printf("  %d: key=%s  val0=%s\n", i+1, kv->key, kv->values[0]);
    }

  printf("\n");
  int i = 0;
  for ( listNode_t *node = kvl->head; node; node = node->next )
    {
      keyValues_t *kv = *(keyValues_t **)node->data;
      if ( kv ) printf("  %d: key=%s  val0=%s\n", i+1, kv->key, kv->values[0]);
      ++i;
    }

  list_destroy(pml);
  printf("Successfully freed %d kvlists...\n", numLists);

  return 0;
}
*/
