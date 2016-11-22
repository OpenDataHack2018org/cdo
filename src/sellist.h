#ifndef _SELLIST_H
#define _SELLIST_H

#include "pmlist.h"

typedef union {
  int ival;
  double dval;
  const char *cval;
} cvalues_t;

typedef struct {
  int nvalues;
  char *key;
  char **values;
  bool *flag;
  int type;
  char *txt;
  cvalues_t *cvalues;
} selentry_t;


typedef struct {
  int size;
  selentry_t *entry;
} sellist_t;


#define  SELLIST_INT         1
#define  SELLIST_FLT         2
#define  SELLIST_WORD        3

#define  SELLIST_DEF_INT(name)                  int x_##name = 0
#define  SELLIST_DEF_FLT(name)                  double x_##name = 0
#define  SELLIST_DEF_WORD(name)                 const char *x_##name = 0
#define  SELLIST_ADD_INT(sellist, name, txt)    SELLIST_DEF_INT(name);  int xpid_##name = sellist_add(sellist, txt, #name, SELLIST_INT)
#define  SELLIST_ADD_FLT(sellist, name, txt)    SELLIST_DEF_FLT(name);  int xpid_##name = sellist_add(sellist, txt, #name, SELLIST_FLT)
#define  SELLIST_ADD_WORD(sellist, name, txt)   SELLIST_DEF_WORD(name); int xpid_##name = sellist_add(sellist, txt, #name, SELLIST_WORD)
#define  SELLIST_NOCC(name)                     sellist_nvalues(sellist, xpid_##name)
#define  SELLIST_CHECK(name)                    sellist_check(sellist, xpid_##name, &x_##name)
#define  SELLIST_CHECK_FLAG(name)               sellist_check_flag(sellist, xpid_##name)


sellist_t *sellist_create(list_t *kvlist);
void sellist_destroy(sellist_t *sellist);
int sellist_add(sellist_t *sellist, const char *txt, const char *name, int type);
int sellist_nvalues(sellist_t *sellist, int idx);
bool sellist_check(sellist_t *sellist, int idx, void *par);
void sellist_check_flag(sellist_t *sellist, int idx);

#endif
