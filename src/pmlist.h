#ifndef _PMLIST_H
#define _PMLIST_H

#include "list.h"

typedef struct {
  int nvalues;
  char *key;
  char **values;
} keyValues_t;

keyValues_t *kvlist_search(list_t *kvl, const char *key);
list_t *pml_search_kvl(list_t *pml, const char *key, const char *value);

bool kvl_print_iter(void *data);
bool pml_print_iter(void *data);

void kvlist_print(list_t *kvl);

void free_keyval(void *data);
void free_kvlist(void *data);

void kvlist_append(list_t *kvl, const char *key, const char **values, int nvalues);

int kvlist_parse_cmdline(list_t *kvl, int nparams, char **params);

#endif
