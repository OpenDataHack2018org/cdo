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

list_t *kvlist_new(const char *name);
void kvlist_destroy(list_t *list);
void kvlist_append(list_t *kvl, const char *key, const char **values, int nvalues);
int kvlist_parse_cmdline(list_t *kvl, int nparams, char **params);

list_t *pml_search_kvl_ventry(list_t *pml, const char *key, const char *value, int nentry, const char **entry);
list_t *pml_get_kvl_ventry(list_t *pml, int nentry, const char **entry);

#endif
