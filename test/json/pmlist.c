#include <stdio.h>
#include <string.h>
#include <stdlib.h>


struct keyvalues
{
  int nvalues;
  char **values;
  char *key;
};


struct kvlist_node {
  struct keyvalues kv;
  struct kvlist_node *next;
};


struct kvlist {
  struct kvlist_node *head;
  char *name;
  int length;
};


struct kvlist_node *kvlist_search(struct kvlist *kvl, const char *key)
{
  if ( key )
    {
      struct kvlist_node *node = kvl->head;
      while ( node )
        {
          if ( node->kv.key && strcmp(node->kv.key, key) == 0 ) return node;
          node = node->next;
        }
    }

  return NULL;
}


struct keyvalues *kvlist_entry(struct kvlist *kvl, int index)
{
  if ( kvl )
    {
      int i = 0;
      struct kvlist_node *node = kvl->head;
      while ( node )
        {
          if ( i == index ) return &(node->kv);
          node = node->next;
          ++i;
        }
    }

  return NULL;
}


void kvlist_node_print(struct kvlist_node *node)
{
  if ( node )
    {
      printf("\t%s\t= ", node->kv.key);
      int nvalues = node->kv.nvalues;
      for ( int i = 0; i < nvalues; ++i ) printf(" \"%s\"", node->kv.values[i]);
      printf("\n");
    }
}

void kvlist_node_fill(struct kvlist_node *node, const char *key, char *values[], int nvalues)
{
  node->next = NULL;
  node->kv.key = strdup(key);
  node->kv.nvalues = nvalues;
  node->kv.values = (char **) malloc(nvalues*sizeof(char*));
  for ( int i = 0; i < nvalues; ++i ) node->kv.values[i] = strdup(values[i]);
}


struct kvlist_node *kvlist_node_new()
{
  return (struct kvlist_node *) malloc(sizeof(struct kvlist_node));
}


void kvlist_append(struct kvlist *kvl, struct kvlist_node *node)
{
  if ( kvl )
    {
      struct kvlist_node **list = &(kvl->head);
      while ( *list ) list = &(*list)->next;
      *list = node;
      kvl->length++;
    }
}


void kvlist_delete(struct kvlist **pkvl)
{
  struct kvlist *kvl = *pkvl;
  if ( kvl )
    {
      struct kvlist_node *node = kvl->head;
      struct kvlist_node *next;
      struct keyvalues *kv;

      while ( node )
        {
          kv = &(node->kv);
          int nvalues = kv->nvalues;
          if ( kv->values )
            {
              for ( int i = 0; i < nvalues; ++i ) if ( kv->values[i] ) free(kv->values[i]);
              free(kv->values);
            }
          if ( kv->key ) free(kv->key);
          next = node->next;
          free(node);
          node = next;
        }

      if ( kvl->name ) free(kvl->name);
      free(kvl);
      *pkvl = NULL;
    }
}


void kvlist_print(struct kvlist *kvl)
{
  if ( kvl )
    {
      printf("  kvlist=%s n=%d:\n", kvl->name, kvl->length);
      struct kvlist_node *node = kvl->head;
      while ( node )
        {
          kvlist_node_print(node);
          node = node->next;
        }
    }
}


struct kvlist *kvlist_new(const char *name)
{
  struct kvlist *kvl = (struct kvlist *) malloc(sizeof(struct kvlist));
  kvl->head = NULL;
  kvl->name = name ? strdup(name) : NULL;
  kvl->length = 0;
  return kvl;
};


struct param_node {
  struct param_node *next;
  struct kvlist *kvl;
};


struct pmlist{
  struct param_node *head;
  char *name;
  int length;
};


struct pmlist *pmlist_new(const char *name)
{
  struct pmlist *pml = (struct pmlist *) malloc(sizeof(struct pmlist));
  pml->head = NULL;
  pml->name = name ? strdup(name) : NULL;
  pml->length = 0;
  return pml;
};


void pmlist_delete(struct pmlist **ppml)
{
  struct pmlist *pml = *ppml;
  if ( pml )
    {
      struct param_node *node = pml->head;
      struct param_node *next;

      while ( node )
        {
          if ( node->kvl ) kvlist_delete(&(node->kvl));
          next = node->next;
          free(node);
          node = next;
        }

      if ( pml->name ) free(pml->name);
      free(pml);
      *ppml = NULL;
    }
}


void pmlist_print(struct pmlist *pml)
{
  if ( pml )
    {
      printf("pmlist=%s n=%d:\n", pml->name, pml->length);
      struct param_node *node = pml->head;
      while ( node )
        {
          kvlist_print(node->kvl);
          node = node->next;
        }
    }
}


void pmlist_append(struct pmlist *pml, struct kvlist *kvl)
{
  if ( pml )
    {
      struct param_node **list = &(pml->head);
      while ( *list ) list = &(*list)->next;
      struct param_node *node = (struct param_node *) malloc(sizeof(struct param_node));
      node->next = NULL;
      node->kvl = kvl;
      *list = node;
      pml->length++;
    }
}


struct kvlist *pmlist_search(struct pmlist *pml, const char *key, const char *value)
{
  if ( pml && key && value )
    {
      struct param_node *node = pml->head;
      while ( node )
        {
          if ( node->kvl )
            {
              struct kvlist_node *kvnode = kvlist_search(node->kvl, key);
              if ( kvnode && kvnode->kv.nvalues > 0 && strcmp(kvnode->kv.values[0], value) == 0 ) return node->kvl;
            }
          node = node->next;
        }
    }

  return NULL;
}


int main(void)
{
  struct pmlist *pml = pmlist_new("list1");

  {
    struct kvlist *kvl = kvlist_new("parameter");
    struct kvlist_node *node;

    char *k1 = "longname", *k1vals[] = {"surface temperature"};
    char *k2 = "name",     *k2vals[] = {"temperature"};
    char *k3 = "values",   *k3vals[] = {"273.15", "292.5", "301.4"};

    node = kvlist_node_new();
    kvlist_node_fill(node, k1, k1vals, sizeof(k1vals)/sizeof(k1vals[0]));
    kvlist_append(kvl, node);

    node = kvlist_node_new();
    kvlist_node_fill(node, k2, k2vals, sizeof(k2vals)/sizeof(k2vals[0]));
    kvlist_append(kvl, node);

    node = kvlist_node_new();
    kvlist_node_fill(node, k3, k3vals, sizeof(k3vals)/sizeof(k3vals[0]));
    kvlist_append(kvl, node);

    pmlist_append(pml, kvl);
  }
  {
    struct kvlist *kvl = kvlist_new(NULL);
    struct kvlist_node *node;

    char *k1 = "longname", *k1vals[] = {"surface pressure"};
    char *k2 = "name",     *k2vals[] = {"pressure"};
    char *k3 = "values",   *k3vals[] = {"1000", "850", "500"};
    char *k4 = "units",    *k4vals[] = {"hPa"};

    node = kvlist_node_new();
    kvlist_node_fill(node, k1, k1vals, sizeof(k1vals)/sizeof(k1vals[0]));
    kvlist_append(kvl, node);

    node = kvlist_node_new();
    kvlist_node_fill(node, k2, k2vals, sizeof(k2vals)/sizeof(k2vals[0]));
    kvlist_append(kvl, node);

    node = kvlist_node_new();
    kvlist_node_fill(node, k3, k3vals, sizeof(k3vals)/sizeof(k3vals[0]));
    kvlist_append(kvl, node);

    node = kvlist_node_new();
    kvlist_node_fill(node, k4, k4vals, sizeof(k4vals)/sizeof(k4vals[0]));
    kvlist_append(kvl, node);

    pmlist_append(pml, kvl);
  }
  {
    struct kvlist *kvl = kvlist_new(NULL);
    struct kvlist_node *node;

    char *k1 = "longname", *k1vals[] = {"Air Temperature"};
    char *k2 = "name",     *k2vals[] = {"ta"};
    char *k3 = "valid_max",*k3vals[] = {"336"};
    char *k4 = "units",    *k4vals[] = {"K"};

    node = kvlist_node_new();
    kvlist_node_fill(node, k1, k1vals, sizeof(k1vals)/sizeof(k1vals[0]));
    kvlist_append(kvl, node);

    node = kvlist_node_new();
    kvlist_node_fill(node, k2, k2vals, sizeof(k2vals)/sizeof(k2vals[0]));
    kvlist_append(kvl, node);

    node = kvlist_node_new();
    kvlist_node_fill(node, k3, k3vals, sizeof(k3vals)/sizeof(k3vals[0]));
    kvlist_append(kvl, node);

    node = kvlist_node_new();
    kvlist_node_fill(node, k4, k4vals, sizeof(k4vals)/sizeof(k4vals[0]));
    kvlist_append(kvl, node);

    pmlist_append(pml, kvl);
  }
  

  pmlist_print(pml);

  struct kvlist *kvl = pmlist_search(pml, "name", "pressure");
  printf("kvlist of name=pressure:\n");
  kvlist_print(kvl);

  printf("\n");
  int n = kvl->length;
  for ( int i = 0; i < n; ++i )
    {
      struct keyvalues *kv = kvlist_entry(kvl, i);
      if ( kv ) printf("  %d: key=%s  val0=%s\n", i+1, kv->key, kv->values[0]);
    }

  printf("\n");
  int i = 0;
  for ( struct kvlist_node *node = kvl->head; node; node=node->next )
    {
      struct keyvalues *kv = &(node->kv);
      if ( kv ) printf("  %d: key=%s  val0=%s\n", i+1, kv->key, kv->values[0]);
      ++i;
    }
  
  pmlist_delete(&pml);

  return 0;
}
