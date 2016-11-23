#include <cdo_int.h>
#include "sellist.h"

//#define SELDEBUG 1

sellist_t *sellist_create(list_t *kvlist)
{
  sellist_t *sellist = (sellist_t *) malloc(sizeof(sellist_t));
  sellist->size = list_size(kvlist);
  sellist->entry = (selentry_t *) malloc(sellist->size*sizeof(selentry_t));

  int i = 0;
  for ( listNode_t *kvnode = kvlist->head; kvnode; kvnode = kvnode->next )
    {
      keyValues_t *kv = *(keyValues_t **)kvnode->data;
      selentry_t *e = &(sellist->entry[i]);
      e->key = kv->key;
      e->values = kv->values;
      e->nvalues = kv->nvalues;
#ifdef SELDEBUG
      printf("%s =", e->key);
      for ( int ii = 0; ii < e->nvalues; ++ii ) printf(" '%s'", e->values[ii]);
      printf("\n");
#endif
      ++i;
    }

  for ( int i = 0; i < sellist->size; ++i )
    {
      selentry_t *e = &(sellist->entry[i]);
      e->flag = NULL;
      e->cvalues = NULL;
      if ( e->nvalues )
        {
          e->flag = (bool*) calloc(e->nvalues, sizeof(bool));
          e->cvalues = (cvalues_t*) calloc(e->nvalues, sizeof(cvalues_t));
        }
#ifdef SELDEBUG
      printf("%s =", e->key);
      for ( int ii = 0; ii < e->nvalues; ++ii ) printf(" '%s'", e->values[ii]);
      printf("\n");
#endif
    }

  return sellist;
}


void sellist_destroy(sellist_t *sellist)
{
  if ( sellist )
    {
      for ( int i = 0; i < sellist->size; ++i )
        {
          selentry_t *e = &(sellist->entry[i]);
          if ( e->flag ) free(e->flag);
          if ( e->txt ) free(e->txt);
          if ( e->cvalues ) free(e->cvalues);
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
#ifdef SELDEBUG          
          printf("%s =", e->key);
          for ( int i = 0; i < nvalues; ++i )
            switch (type)
              {
              case SELLIST_INT:  printf(" %d", e->cvalues[i].ival); break;
              case SELLIST_FLT:  printf(" %g", e->cvalues[i].dval); break;
              case SELLIST_WORD: printf(" %s", e->cvalues[i].cval); break;
              }
          printf("\n");
#endif
        }
    }

  return idx;
}


int sellist_nvalues(sellist_t *sellist, int idx)
{
  int nvalues = 0;

  if ( sellist && idx >= 0 && idx < sellist->size ) nvalues = sellist->entry[idx].nvalues;

  return nvalues;
}


void sellist_check_flag(sellist_t *sellist, int idx)
{
  if ( idx < 0 || idx >= sellist->size ) return;

  int nvalues = sellist_nvalues(sellist, idx);

  if ( nvalues )
    {
      selentry_t *e = &(sellist->entry[idx]);
      for ( int i = 0; i < nvalues; ++i )
        if ( e->flag[i] == false ) cdoWarning("%s >%s< not found!", e->txt, e->values[i]);
    }
}


bool sellist_check(sellist_t *sellist, int idx, void *par)
{
  bool found = false;

  if ( idx < 0 || idx >= sellist->size ) return found;

  int nvalues = sellist_nvalues(sellist, idx);

  if ( nvalues )
    {
      selentry_t *e = &(sellist->entry[idx]);
      int type = e->type;
      for ( int i = 0; i < nvalues; ++i )
        {
          switch (type)
            {
            case SELLIST_INT:  if ( *(int*)par == e->cvalues[i].ival )                     { found = true; e->flag[i] = true; } break;
            case SELLIST_FLT:  if ( fabs(*(double*)par - e->cvalues[i].dval) < 1.e-4 )     { found = true; e->flag[i] = true; } break;
            case SELLIST_WORD: if ( wildcardmatch(e->cvalues[i].cval, *(char**)par) == 0 ) { found = true; e->flag[i] = true; } break;
            }
        }
    }

  return found;
}


bool sellist_check_date(sellist_t *sellist, int idx, const char *par)
{
  bool found = false;

  if ( idx < 0 || idx >= sellist->size ) return found;

  int nvalues = sellist_nvalues(sellist, idx);

  if ( nvalues )
    {
      char wcdate[512];
      selentry_t *e = &(sellist->entry[idx]);

      if ( *par == ' ' ) ++par;

      for ( int i = 0; i < nvalues; ++i )
        {
          strcpy(wcdate, e->values[i]);
          strcat(wcdate, "*");
          if ( wildcardmatch(wcdate, par) == 0 ) { found = true; e->flag[i] = true; }
        }
    }

  return found;
}

void season_to_months(const char *season, int *imonths);

bool sellist_check_season(sellist_t *sellist, int idx, int month)
{
  assert(month>=1&&month<=12);
  bool found = false;

  if ( idx < 0 || idx >= sellist->size ) return found;

  int nvalues = sellist_nvalues(sellist, idx);

  if ( nvalues )
    {
      int imon[13]; /* 1-12 ! */
      selentry_t *e = &(sellist->entry[idx]);

      for ( int i = 0; i < nvalues; ++i )
        {
          for ( int m = 0; m < 13; ++m ) imon[m] = 0;
          season_to_months(e->values[i], imon);
          if ( imon[month] ) { found = true; e->flag[i] = true; }
        }
    }

  return found;
}


void sellist_def_flag(sellist_t *sellist, int idx, int vindex, bool flag)
{
  if ( idx < 0 || idx >= sellist->size ) return;

  int nvalues = sellist_nvalues(sellist, idx);

  if ( nvalues )
    {
      selentry_t *e = &(sellist->entry[idx]);
      if ( vindex >= 0 && vindex < nvalues ) e->flag[vindex] = flag;
    }  
}


void sellist_get_par(sellist_t *sellist, int idx, int vindex, void *par)
{
  if ( idx < 0 || idx >= sellist->size ) return;

  int nvalues = sellist_nvalues(sellist, idx);

  if ( nvalues )
    {
      selentry_t *e = &(sellist->entry[idx]);
      int type = e->type;
      if ( vindex >= 0 && vindex < nvalues )
        {
          switch (type)
            {
            case SELLIST_INT:  *(int*)par = e->cvalues[vindex].ival; break;
            case SELLIST_FLT:  *(double*)par = e->cvalues[vindex].dval; break;
            case SELLIST_WORD: *(const char**)par = e->cvalues[vindex].cval; break;
            }
        }
    }
}


void sellist_def_par(sellist_t *sellist, int idx, int vindex, void *par)
{
  if ( idx < 0 || idx >= sellist->size ) return;

  int nvalues = sellist_nvalues(sellist, idx);

  if ( nvalues )
    {
      selentry_t *e = &(sellist->entry[idx]);
      int type = e->type;
      if ( vindex >= 0 && vindex < nvalues )
        {
          switch (type)
            {
            case SELLIST_INT:  e->cvalues[vindex].ival = *(int*)par; break;
            case SELLIST_FLT:  e->cvalues[vindex].dval = *(double*)par; break;
            case SELLIST_WORD: e->cvalues[vindex].cval = *(const char**)par; break;
            }
        }
    }
}
