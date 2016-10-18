#include "cdo_int.h"
#include "pml.h"

static
void pml_init(pml_t *pml, const char *name)
{
  pml->size = 0;
  pml->name = strdup(name);
}


pml_t *pml_create(const char *name)
{
  pml_t *pml = (pml_t*) Malloc(sizeof(pml_t));

  pml_init(pml, name);

  return pml;
}


void pml_destroy(pml_t *pml)
{
  if ( pml == NULL ) return;

  for ( int i = 0; i < pml->size; ++i )
    {
      if ( pml->entry[i]->txt ) free(pml->entry[i]->txt);
      if ( pml->entry[i]->name ) free(pml->entry[i]->name);
      if ( pml->entry[i] ) Free(pml->entry[i]);
    }

  if ( pml->name ) free(pml->name);

  Free(pml);
}


void pml_print(pml_t *pml)
{
  if ( pml == NULL ) return;

  fprintf(stdout, "Parameter list: %s\n", pml->name);
  fprintf(stdout, " Num  Name             Type  Size   Occ  Entries\n");

  for ( int i = 0; i < pml->size; ++i )
    {
      pml_entry_t *entry = pml->entry[i];
      fprintf(stdout, "%4d  %-16s %4d  %4d  %4d ",
	      i+1, pml->entry[i]->name, pml->entry[i]->type, (int)pml->entry[i]->size, pml->entry[i]->occ);
      int nout = pml->entry[i]->occ;
      if ( nout > 8 ) nout = 8;

      for ( int j = 0; j < nout; j++ )
        {
          switch (entry->type)
            {
            case PML_WORD: fprintf(stdout, " %s", ((char **)entry->ptr)[j]); break;
            case PML_INT:  fprintf(stdout, " %d", ((int *)entry->ptr)[j]); break;
            case PML_FLT:  fprintf(stdout, " %g", ((double *)entry->ptr)[j]); break;
            }
        }
      fprintf(stdout, "\n");
    }
}


int pml_add(pml_t *pml, const char *txt, const char *name, int type, void *ptr, size_t size)
{
  if ( pml->size >= MAX_PML_ENTRY )
    {
      fprintf(stderr, "Too many entries in parameter list %s! (Max = %d)\n", pml->name, MAX_PML_ENTRY);
      return -1;
    }

  pml_entry_t *pml_entry = (pml_entry_t*) Malloc(sizeof(pml_entry_t));

  pml_entry->txt  = strdup(txt);
  pml_entry->name = strdup(name);
  pml_entry->len  = strlen(name);
  pml_entry->type = type;
  pml_entry->ptr  = ptr;
  pml_entry->size = size;
  pml_entry->occ  = 0;

  int entry = pml->size;
  pml->entry[pml->size++] = pml_entry;

  return entry;
}


int pml_num_par(pml_t *pml, int pid)
{
  int nocc = 0;

  if ( pml && pid >= 0 && pid < pml->size ) nocc = pml->entry[pid]->occ;

  return nocc;
}


int pml_num(pml_t *pml, const char *name)
{
  pml_entry_t *entry;
  int i, nocc = 0;

  if ( pml == NULL ) return nocc;

  for ( i = 0; i < pml->size; i++ )
    {
      entry = pml->entry[i];
      if ( strcmp(name, entry->name) == 0 )
	{
	  nocc = entry->occ;
	  break;
	}
    }

  if ( i == pml->size )
    fprintf(stderr, "Parameter list entry %s not found in %s\n", name, pml->name);

  return nocc;
}

void split_intstring(const char *intstr, int *first, int *last, int *inc);

int pml_add_entry(pml_entry_t *entry, char *arg)
{
  int status = 0;

  if ( entry->type == PML_INT )
    {
      int ival, first, last, inc;

      split_intstring(arg, &first, &last, &inc);

      if ( inc >= 0 )
	{
	  for ( ival = first; ival <= last; ival += inc )
	    if ( entry->occ < (int) entry->size )
	      ((int *) entry->ptr)[entry->occ++] = ival;
	}
      else
	{
	  for ( ival = first; ival >= last; ival += inc )
	    if ( entry->occ < (int) entry->size )
	      ((int *) entry->ptr)[entry->occ++] = ival;
	}
    }
  else if ( entry->type == PML_FLT )
    {
      if ( entry->occ < (int) entry->size )
	((double *) entry->ptr)[entry->occ++] = atof(arg);
    }
  else if ( entry->type == PML_WORD )
    {
      if ( entry->occ < (int) entry->size )
	((char **) entry->ptr)[entry->occ++] = strdupx(arg);
    }
  else
    {
      fprintf(stderr, "Unsupported type!\n");
    }

  return status;
}

static
void pml_process(pml_entry_t *entry, int argc, char **argv)
{
  for ( int i = 0; i < argc; ++i )
    {
      char *parg = argv[i];
      if ( i == 0 )
	{
	  char *epos = strchr(parg, '=');
	  if ( epos == NULL ) fprintf(stderr, "Internal problem, keyword not found!\n");
	  parg += epos-parg+1;
	}

      pml_add_entry(entry, parg);
    }
}


int pml_read(pml_t *pml, int argc, char **argv)
{
  pml_entry_t *entry = NULL;
  pml_entry_t *pentry[MAX_PML_ENTRY];
  int params[MAX_PML_ENTRY];
  int num_par[MAX_PML_ENTRY];
  int nparams = 0;
  int i;
  char *epos;
  size_t len;
  int bufsize = 0;
  int status = 0;

  if ( argc == 0 ) return 0;

  for ( i = 0; i < argc; ++i )
    {
      len = strlen(argv[i]);
      bufsize += len+1;
    }

  char *parbuf = (char*) Malloc(bufsize*sizeof(char));
  memset(parbuf, 0, bufsize*sizeof(char));

  int istart = 0;
  while ( istart < argc )
    {
      epos = strchr(argv[istart], '=');
      if ( epos == NULL )
	{
	  fprintf(stderr, "Parameter >%s< has no keyword!\n", argv[istart]);
	  status = 1;
	  goto END_LABEL;
	}

      len = epos - argv[istart];
      for ( i = 0; i < pml->size; ++i )
	{
	  entry = pml->entry[i];
	  if ( entry->len == len )
	    if ( memcmp(entry->name, argv[istart], len) == 0 ) break;
	}

      if ( i == pml->size )
	{
	  fprintf(stderr, "Parameter >%s< has an invalid keyword!\n", argv[istart]);
	  status = 2;
	  goto END_LABEL;
	}

      num_par[nparams] = 0;
      pentry[nparams]  = entry;
      params[nparams]  = istart;
      num_par[nparams] = 1;
      
      istart++;
      for ( i = istart; i < argc; ++i )
	{
	  if ( *argv[i] == 0 ) { i++; break;}
	  epos = strchr(argv[i], '=');
	  if ( epos != NULL ) break;

	  num_par[nparams]++;
	}

      istart = i;

      nparams++;
    }

  for ( i = 0; i < nparams; ++i )
    {
      pml_process(pentry[i], num_par[i], &argv[params[i]]);
    }


 END_LABEL:

  Free(parbuf);

  return status;
}


bool pml_check_int(pml_t *pml, int pid, int *parlist, bool *flaglist, int par)
{
  bool found = false;
  int npar = pml_num_par(pml, pid);
  for ( int i = 0; i < npar; i++ )
    if ( par == parlist[i] ) { found = true; flaglist[i] = true;/* break;*/}

  return found;
}


bool pml_check_flt(pml_t *pml, int pid, double *parlist, bool *flaglist, double par)
{
  bool found = false;
  int npar = pml_num_par(pml, pid);
  for ( int i = 0; i < npar; i++ )
    if ( fabs(par - parlist[i]) < 1.e-4 ) { found = true; flaglist[i] = true;/* break;*/}

  return found;
}


bool pml_check_word(pml_t *pml, int pid, char **parlist, bool *flaglist, const char *par)
{
  bool found = false;
  int npar = pml_num_par(pml, pid);
  for ( int i = 0; i < npar; i++ )
    if ( wildcardmatch(parlist[i], par) == 0 ) { found = true; flaglist[i] = true;/* break;*/}

  return found;
}


bool pml_check_date(pml_t *pml, int pid, char **parlist, bool *flaglist, const char *par)
{
  bool found = false;
  int npar = pml_num_par(pml, pid);
  char wcdate[512];

  if ( *par == ' ' ) ++par;

  for ( int i = 0; i < npar; i++ )
    {
      strcpy(wcdate, parlist[i]);
      strcat(wcdate, "*");
      if ( wildcardmatch(wcdate, par) == 0 ) { found = true; flaglist[i] = true;/* break;*/}
    }

  return found;
}

void season_to_months(const char *season, int *imonths);

bool pml_check_season(pml_t *pml, int pid, char **parlist, bool *flaglist, int month)
{
  assert(month>=1&&month<=12);
  bool found = false;
  int npar = pml_num_par(pml, pid);
  int imon[13]; /* 1-12 ! */

  for ( int i = 0; i < npar; i++ )
    {
      for ( int m = 0; m < 13; ++m ) imon[m] = 0;
      season_to_months(parlist[i], imon);
      if ( imon[month] ) { found = true; flaglist[i] = true;/* break;*/}
    }

  return found;
}


void pml_check_int_flag(pml_t *pml, int pid, int *parlist, bool *flaglist)
{
  int npar = pml_num_par(pml, pid);
  for ( int i = 0; i < npar; ++i )
    if ( flaglist[i] == false )
      cdoWarning("%s >%d< not found!", pml->entry[pid]->txt, parlist[i]);
}


void pml_check_flt_flag(pml_t *pml, int pid, double *parlist, bool *flaglist)
{
  int npar = pml_num_par(pml, pid);
  for ( int i = 0; i < npar; ++i )
    if ( flaglist[i] == false )
      cdoWarning("%s >%g< not found!", pml->entry[pid]->txt, parlist[i]);
}


void pml_check_word_flag(pml_t *pml, int pid, char **parlist, bool *flaglist)
{
  int npar = pml_num_par(pml, pid);
  for ( int i = 0; i < npar; ++i )
    if ( flaglist[i] == false )
      cdoWarning("%s >%s< not found!", pml->entry[pid]->txt, parlist[i]);
}
