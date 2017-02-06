#include <stdlib.h>
#include "uthash.h"

struct kv
{
  char *key;
  char *value;
  UT_hash_handle hh;
};

void hinsert(struct kv **ht, const char *key, const char *value)
{
  /* Insert new keys. Do not overwrite values of existing keys. */
  struct kv *e, *s;
  HASH_FIND_STR(*ht, key, s);
  if ( s == NULL)
    {
      e = malloc(sizeof(struct kv));
      e->key = malloc(strlen(key) + 1);
      e->value = malloc(strlen(value) + 1);
      strcpy(e->key, key);
      strcpy(e->value, value);
      HASH_ADD_KEYPTR(hh, *ht, e->key, strlen(e->key), e);
    }
}
