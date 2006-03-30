/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2006 Uwe Schulzweida, schulzweida@dkrz.de
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>

#include "namelist.h"

#if ! defined (strdup)
char *strdup(const char *s);
#endif
int readline(FILE *fp, char *line, int len);


#define  func_1         -1 /* nptype */
#define  func_2         -2 /* nptype */
#define  func_3         -3 /* nptype */
#define  NML_NIX         0 /* nptype */
#define  NML_TEXTU       5
#define  NML_TEXTL       6
#define  NML_NPR         7
#define  NML_NUMBER      8 /* nptype */
#define  NML_KEYWORD     9 /* nptype */

#define  PRINT_NOT       1
#define  PRINT_ALL       2

#undef   TRUE
#define  TRUE   1
#undef   FALSE
#define  FALSE  0

#undef   MIN
#define  MIN(a,b)  ((a) < (b) ? (a) : (b))
#undef   MAX
#define  MAX(a,b)  ((a) > (b) ? (a) : (b))
#undef   NINT
#define  NINT(x)   ((x) < 0 ? (int)((x)-0.5) : (int)((x)+0.5))

struct PGMSTAT
{
  int pdis;
  int intract;
};

struct PGMSTAT pgmstat;

#define MAXLINL 1024

struct NAMELINE
{
  int nptype, namitf, namitl;
  char lineac[MAXLINL], lineuc[MAXLINL], linelc[MAXLINL];
};

struct NAMELINE nameline;


static void namelist_init(NAMELIST *namelist, const char *name)
{
  namelist->size = 0;
  namelist->name = strdup(name);
}


NAMELIST *namelistNew(const char *name)
{
  NAMELIST *namelist;

  namelist = (NAMELIST *) malloc(sizeof(NAMELIST));

  namelist_init(namelist, name);

  return (namelist);
}


void namelistDelete(NAMELIST *nml)
{
  int i;

  if ( nml )
    {
      for ( i = 0; i < nml->size; i++ )
	{
	  if ( nml->entry[i]->name ) free(nml->entry[i]->name);
	  free(nml->entry[i]);
	}
      
      if ( nml->name ) free(nml->name);
      free(nml);
    }
}


void namelistPrint(NAMELIST *nml)
{
  NML_ENTRY *entry;
  int i, j, nout;

  if ( nml == NULL ) return;

  fprintf(stdout, "Namelist: %s\n", nml->name);
  fprintf(stdout, " Num  Name       Type  Size   Dis   Occ  Entries\n");

  for ( i = 0; i < nml->size; i++ )
    {
      entry = nml->entry[i];
      fprintf(stdout, "%4d  %-10s %4d  %4d  %4d  %4d ",
	      i+1, nml->entry[i]->name, nml->entry[i]->type, (int)nml->entry[i]->size,
	      nml->entry[i]->dis, nml->entry[i]->occ);
      nout = nml->entry[i]->occ;
      if ( nout > 8 ) nout = 8;
      if      ( entry->type >= NML_TEXT )
	fprintf(stdout, " '%s'", ((char *)entry->ptr));
      else if ( entry->type == NML_WORD )
	for ( j = 0; j < nout; j++ )
	  fprintf(stdout, " %s", ((char **)entry->ptr)[j]);
      else if ( entry->type == NML_INT )
	for ( j = 0; j < nout; j++ )
	  fprintf(stdout, " %d", ((int *)entry->ptr)[j]);
      else if ( entry->type == NML_DOUBLE )
	for ( j = 0; j < nout; j++ )
	  fprintf(stdout, " %g", ((double *)entry->ptr)[j]);
      
      fprintf(stdout, "\n");
    }
}


void namelistAdd(NAMELIST *nml, const char *name, int type, int dis, void *ptr, size_t size)
{
  NML_ENTRY *nml_entry;

  if ( nml->size >= MAX_NML_ENTRY )
    {
      fprintf(stderr, "Too many namelist entries in %s! (Max = %d)\n", nml->name, MAX_NML_ENTRY);
      return;
    }

  nml_entry = (NML_ENTRY *) malloc(sizeof(NML_ENTRY));

  nml_entry->name = strdup(name);
  nml_entry->type = type;
  nml_entry->ptr  = ptr;
  nml_entry->size = size;
  nml_entry->dis  = dis;
  nml_entry->occ  = 0;

  nml->entry[nml->size++] = nml_entry;
}


int namelistNum(NAMELIST *nml, const char *name)
{
  NML_ENTRY *entry;
  int i, nocc = 0;

  if ( nml == NULL ) return (nocc);

  for ( i = 0; i < nml->size; i++ )
    {
      entry = nml->entry[i];
      if ( strcmp(name, entry->name) == 0 )
	{
	  nocc = entry->occ;
	  break;
	}
    }

  if ( i == nml->size )
    fprintf(stderr, "Namelist entry %s not found in %s\n", name, nml->name);

  return (nocc);
}


static void getnite(void)
{
  int nst, i, j;

  nst = nameline.namitl + 1;

  while ( TRUE )
    {
      for ( i = nst; i < MAXLINL; i++ )
	{
	  if ( nameline.linelc[i] == 0 ) break;

          if      (   nameline.linelc[i] == ' ' ) continue;
	  else if (   nameline.linelc[i] == '=' ) continue;
          else if (   nameline.linelc[i] == ',' ) continue;
          else if ( ((nameline.linelc[i] >= 'a')  &&
		     (nameline.linelc[i] <= 'z')) ||
	             (nameline.linelc[i] == '$')  ||
		     (nameline.linelc[i] == '&') )
	    {
	      nameline.nptype = NML_KEYWORD;
	      nameline.namitf = i;
	      for ( j = nameline.namitf+1; j < MAXLINL; j++ )
		{
		  if ( !(islower((int) nameline.linelc[j]) ||
			 isdigit((int) nameline.linelc[j])) )
		    {
		      nameline.namitl = j - 1;
		      return;
		    }
		}
	      nameline.namitl = MAXLINL;
	      return;
	    }
          else if ( nameline.linelc[i] == '\'' ||
		    nameline.linelc[i] == '\"' ||
		    nameline.linelc[i] == '`'  ||
		    nameline.linelc[i] == '/' )
	    {
	      nameline.nptype = NML_TEXT;
	      nameline.namitf = i;
	      for ( j = nameline.namitf+1; j < MAXLINL; j++ )
		if (nameline.linelc[j] == nameline.linelc[nameline.namitf])
		  {
 		    nameline.namitl = j;
		    return;
		  }
	      nameline.namitl = MAXLINL + 1;
	      return;
	    }
          else if ( (nameline.linelc[i] >= '0'  &&
		     nameline.linelc[i] <= '9') ||
		     nameline.linelc[i] == '+'  ||
		     nameline.linelc[i] == '-'  ||
		     nameline.linelc[i] == '.' )
	    {
	      nameline.nptype = NML_NUMBER;
	      nameline.namitf = i;
	      for ( j = i+1; j < MAXLINL; j++)
		{
		  if ( nameline.linelc[j] >= '0' && nameline.linelc[j] <= '9' ) continue;
	          else if ( nameline.linelc[j] == '+' ) continue;
		  else if ( nameline.linelc[j] == '-' ) continue;
         	  else if ( nameline.linelc[j] == '.' ) continue;
		  else if ( nameline.linelc[j] == 'e' ) continue;
                  else
		    {
		      nameline.namitl = j - 1;
		      return;
		    }
		}
	      nameline.namitl = MAXLINL;
	      return;
	    }
        }

      if ( ! readline(stdin, nameline.lineac, MAXLINL) ) break;

      for ( i = 0; i < MAXLINL; i++ )
	{
	  nameline.linelc[i] = tolower(nameline.lineac[i]);
	  nameline.lineuc[i] = toupper(nameline.lineac[i]);
        }
      nst = 0;
    }

  nameline.nptype = NML_NIX;
}


static void rdnlsgl(void *var, int ntyp, int nlen, int *nocc)
{
  if ( nameline.nptype == NML_NUMBER )
    {
      if ( *nocc >= nlen )
	{
	  nameline.nptype = func_1;
          return;
	}
      else if ( ntyp == NML_INT )
	{
	  ((int *)var)[*nocc] = atoi(&nameline.lineac[nameline.namitf]);
          *nocc += 1;
	}
      else if ( ntyp == NML_DOUBLE )
	{
	  ((double *)var)[*nocc] = atof(&nameline.lineac[nameline.namitf]);
          *nocc += 1;
	}
      else
	{
          nameline.nptype = func_2;
          return;
        }
    }
  else if ( nameline.nptype == NML_TEXT )
    {
      int i, j=0, newnocc;

      newnocc = MIN(nlen, *nocc+nameline.namitl-nameline.namitf-1);

      if      ( ntyp == NML_TEXT )
	for (i=*nocc; i<newnocc; i++)
	  ((char *)var)[i] = nameline.lineac[nameline.namitf+1+j++];
      else if ( ntyp == NML_TEXTU )
	for (i=*nocc; i<newnocc; i++)
	  ((char *)var)[i] = nameline.lineuc[nameline.namitf+1+j++];
      else if ( ntyp == NML_TEXTL )
	for (i=*nocc; i<newnocc; i++)
	  ((char *)var)[i] = nameline.linelc[nameline.namitf+1+j++];
      else
	{
	  nameline.nptype = func_3;
	  return;
	}

      *nocc = newnocc;
    }
  else if ( nameline.nptype == NML_WORD )
    {
      int i, len;

      if ( *nocc < nlen )
	{
	  len = nameline.namitl - nameline.namitf + 1;
	  ((char **)var)[*nocc] = (char*) calloc((size_t)len+1, sizeof(char));
	  for ( i = 0; i < len; i++ )
	    ((char **)var)[*nocc][i] = nameline.lineac[nameline.namitf+i];
	  *nocc += 1;
	}
    }
  else
    {
      fprintf(stderr, "Namelist parameter type %d unknown!\n", nameline.nptype);
      return;
    }

  nameline.nptype = ntyp;
}


static void nml_print_entry(NML_ENTRY *entry, int ife)
{
  int nout, j;

  if ( entry->size == 0 ) return;

  if ( entry->type == NML_NPR ) return;

  if ( ife == PRINT_ALL )
    nout = entry->occ;
  else
    {
      nout = MAX(entry->occ, entry->dis);
      if ( nout == 0 ) return;
    }

  printf(" %-24s", entry->name);

  if      ( entry->type >= NML_TEXT )
    printf("'%s'", ((char *)entry->ptr));
  else if ( entry->type == NML_WORD )
    for ( j = 0; j < nout; j++ )
      printf(" %s", ((char **)entry->ptr)[j]);
  else if ( entry->type == NML_INT )
    for ( j = 0; j < nout; j++ )
      printf(" %d", ((int *)entry->ptr)[j]);
  else if ( entry->type == NML_DOUBLE )
    for ( j = 0; j < nout; j++ )
      printf(" %g", ((double *)entry->ptr)[j]);

  printf("\n");
}


static void nml_print(NAMELIST *nml, int ife)
{
  int i;

  for ( i = 0; i < nml->size; i++ )
    nml_print_entry(nml->entry[i], ife);
}


void namelistRead(NAMELIST *nml)
{
  /*
    cn  name
    nt  type
    nl  size length
    nc  occ count
    no  dis list
  */
  int clear = FALSE;
  int j, jj, match = -1, wordmatch = -1;
  size_t len;
  char namecx[16], *pnamecx = NULL;
  int nparam;

  nparam = nml->size;

  nameline.namitl = MAXLINL;

 L2000:

  getnite();

  if ( nameline.nptype == NML_NIX )
    {
      goto L3000;
    }
  else if ( nameline.nptype == NML_KEYWORD )
    {
      memset(namecx, '\0', 16);
      len = (size_t) (nameline.namitl - nameline.namitf + 1);

      if ( nameline.lineac[nameline.namitf] == '$' || 
	   nameline.lineac[nameline.namitf] == '&' )
	{
          if ( nameline.namitl-nameline.namitf > 16 ) goto L3000;

          if ( nameline.namitf < nameline.namitl)
	    pnamecx = &nameline.linelc[nameline.namitf+1];
          if ( strncmp(pnamecx, "select", 6) == 0 || 
	       strncmp(pnamecx, "params", 6) == 0 )
	    goto L2000;

	  goto L3000;
        }
      if ( nameline.namitl-nameline.namitf >= 16 ) goto L3000;

      pnamecx = &nameline.linelc[nameline.namitf];
      strncpy(namecx, pnamecx, len);

      if ( len == 5 )
	if ( strncmp(pnamecx, "clear", len) == 0 )
	  {
	    clear = TRUE;
	    goto L2000;
	  }

      match = -1;
      for ( j = 0; j < nparam; j++ )
	{
	  if ( strlen(nml->entry[j]->name) == len )
	    if ( strncmp(pnamecx, nml->entry[j]->name, len) == 0 )
	      {
		jj = j;
		while ( nml->entry[jj]->type == NML_NPR ) jj--;
		if ( match == -1 )
		  match = jj;
		else if ( match != jj )
		  match = -2;
		break;
	      }
	}

      if ( match < 0 )
	{
	  if ( wordmatch >= 0 )
	    {
	      match = wordmatch;
	      nameline.nptype = NML_WORD;
	      goto L777;
	    }

          printf(" * unidentified or ambiguous parameter <%s>\n", namecx);
	  printf(" * valid parameters and values specified so far are\n");

	  nml_print(nml, PRINT_ALL);

          if ( ! pgmstat.intract )
	    {
	      fprintf(stderr, "Namelist error!\n");
	      return;
	    }
        }
      else
	{
          if ( clear ) nml->entry[match]->occ = 0;
	  if ( nml->entry[match]->type == NML_WORD )
	    wordmatch = match;
	  else
	    wordmatch = -1;
        }
      clear = FALSE;
    }
  else
    {
    L777:
      j = match;

      rdnlsgl(nml->entry[j]->ptr, nml->entry[j]->type, (int)nml->entry[j]->size, &nml->entry[j]->occ);

      if ( nameline.nptype != nml->entry[j]->type )
	{
          printf(" * value ignored for parameter <%s> %5d\n", namecx, nameline.nptype);
	  printf(" * valid parameters and values specified so far are\n");

	  nml_print(nml, PRINT_ALL);

          if ( ! pgmstat.intract )
	    {
	      fprintf(stderr, "Namelist error!\n");
	    }
        }
    }

  goto L2000;

 L3000:

  if ( pgmstat.pdis == PRINT_NOT ) return;

  nml_print(nml, pgmstat.pdis);
}

