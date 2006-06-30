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

#define  PRINT_ALL       2
#define  PRINT_MIN       3

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
  int intract;
};

struct PGMSTAT pgmstat;


static void namelist_init(NAMELIST *namelist, const char *name)
{
  namelist->size = 0;
  namelist->dis  = 1;
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
  int i, iocc;

  if ( nml )
    {
      for ( i = 0; i < nml->size; i++ )
	{
	  if ( nml->entry[i]->name ) free(nml->entry[i]->name);
	  if ( nml->entry[i]->type == NML_WORD )
	    for ( iocc = 0; iocc < nml->entry[i]->occ; iocc++ )
	      {
		if ( ((char **)nml->entry[i]->ptr)[iocc] )
		  free(((char **)nml->entry[i]->ptr)[iocc]);
	      }

	  free(nml->entry[i]);
	}
      
      if ( nml->name ) free(nml->name);
      free(nml);
    }
}


void namelistClear(NAMELIST *nml)
{
  int i, iocc;

  if ( nml )
    {
      for ( i = 0; i < nml->size; i++ )
	{
	  if ( nml->entry[i]->type == NML_WORD )
	    for ( iocc = 0; iocc < nml->entry[i]->occ; iocc++ )
	      {
		if ( ((char **)nml->entry[i]->ptr)[iocc] )
		  free(((char **)nml->entry[i]->ptr)[iocc]);
	      }

	  nml->entry[i]->occ = 0;
	}
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
      if ( entry->type >= NML_TEXT )
	{
	  if ( nout > 32 ) nout = 32;
	}
      else
	{
	  if ( nout > 8 ) nout = 8;
	}

      if      ( entry->type >= NML_TEXT )
	{
	  fprintf(stdout, " '");
	  for ( j = 0; j < nout; j++ )
	    fprintf(stdout, "%c", ((char *)entry->ptr)[j]);
	  fprintf(stdout, "'");
	}
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


static void getnite(FILE *nmlfp, NAMELIST *nml)
{
  int nst, i, j;

  nst = nml->line.namitl + 1;

  while ( TRUE )
    {
      for ( i = nst; i < MAX_LINE_LEN; i++ )
	{
	  if ( nml->line.linelc[i] == 0 ) break;

          if      (   nml->line.linelc[i] == ' ' ) continue;
	  else if (   nml->line.linelc[i] == '=' ) continue;
          else if (   nml->line.linelc[i] == ',' ) continue;
          else if ( ((nml->line.linelc[i] >= 'a')  &&
		     (nml->line.linelc[i] <= 'z')) ||
	             (nml->line.linelc[i] == '_')  ||
	             (nml->line.linelc[i] == '$')  ||
		     (nml->line.linelc[i] == '&') )
	    {
	      nml->line.nptype = NML_KEYWORD;
	      nml->line.namitf = i;
	      for ( j = nml->line.namitf+1; j < MAX_LINE_LEN; j++ )
		{
		  if ( !(islower((int) nml->line.linelc[j]) ||
			 (((int) nml->line.linelc[j]) == '_') ||
			 isdigit((int) nml->line.linelc[j])) )
		    {
		      nml->line.namitl = j - 1;
		      return;
		    }
		}
	      nml->line.namitl = MAX_LINE_LEN;
	      return;
	    }
          else if ( nml->line.linelc[i] == '\'' ||
		    nml->line.linelc[i] == '\"' ||
		    nml->line.linelc[i] == '`'  ||
		    nml->line.linelc[i] == '/' )
	    {
	      nml->line.nptype = NML_TEXT;
	      nml->line.namitf = i;
	      for ( j = nml->line.namitf+1; j < MAX_LINE_LEN; j++ )
		if (nml->line.linelc[j] == nml->line.linelc[nml->line.namitf])
		  {
 		    nml->line.namitl = j;
		    return;
		  }
	      nml->line.namitl = MAX_LINE_LEN + 1;
	      return;
	    }
          else if ( (nml->line.linelc[i] >= '0'  &&
		     nml->line.linelc[i] <= '9') ||
		     nml->line.linelc[i] == '+'  ||
		     nml->line.linelc[i] == '-'  ||
		     nml->line.linelc[i] == '.' )
	    {
	      nml->line.nptype = NML_NUMBER;
	      nml->line.namitf = i;
	      for ( j = i+1; j < MAX_LINE_LEN; j++)
		{
		  if ( nml->line.linelc[j] >= '0' && nml->line.linelc[j] <= '9' ) continue;
	          else if ( nml->line.linelc[j] == '+' ) continue;
		  else if ( nml->line.linelc[j] == '-' ) continue;
         	  else if ( nml->line.linelc[j] == '.' ) continue;
		  else if ( nml->line.linelc[j] == 'e' ) continue;
                  else
		    {
		      nml->line.namitl = j - 1;
		      return;
		    }
		}
	      nml->line.namitl = MAX_LINE_LEN;
	      return;
	    }
        }

      if ( ! readline(nmlfp, nml->line.lineac, MAX_LINE_LEN) ) break;

      for ( i = 0; i < MAX_LINE_LEN; i++ )
	{
	  nml->line.linelc[i] = tolower(nml->line.lineac[i]);
	  nml->line.lineuc[i] = toupper(nml->line.lineac[i]);
        }
      nst = 0;
    }

  nml->line.nptype = NML_NIX;
}


static void rdnlsgl(NAMELIST *nml, void *var, int ntyp, int nlen, int *nocc)
{
  if ( nml->line.nptype == NML_NUMBER )
    {
      if ( *nocc >= nlen )
	{
	  nml->line.nptype = func_1;
          return;
	}
      else if ( ntyp == NML_INT )
	{
	  ((int *)var)[*nocc] = atoi(&nml->line.lineac[nml->line.namitf]);
          *nocc += 1;
	}
      else if ( ntyp == NML_DOUBLE )
	{
	  ((double *)var)[*nocc] = atof(&nml->line.lineac[nml->line.namitf]);
          *nocc += 1;
	}
      else
	{
          nml->line.nptype = func_2;
          return;
        }
    }
  else if ( nml->line.nptype == NML_TEXT )
    {
      int i, j=0, newnocc;

      newnocc = MIN(nlen, *nocc+nml->line.namitl-nml->line.namitf-1);

      if      ( ntyp == NML_TEXT )
	for (i=*nocc; i<newnocc; i++)
	  ((char *)var)[i] = nml->line.lineac[nml->line.namitf+1+j++];
      else if ( ntyp == NML_TEXTU )
	for (i=*nocc; i<newnocc; i++)
	  ((char *)var)[i] = nml->line.lineuc[nml->line.namitf+1+j++];
      else if ( ntyp == NML_TEXTL )
	for (i=*nocc; i<newnocc; i++)
	  ((char *)var)[i] = nml->line.linelc[nml->line.namitf+1+j++];
      else
	{
	  nml->line.nptype = func_3;
	  return;
	}

      *nocc = newnocc;
    }
  else if ( nml->line.nptype == NML_WORD )
    {
      int i, len;

      if ( *nocc < nlen )
	{
	  len = nml->line.namitl - nml->line.namitf + 1;
	  ((char **)var)[*nocc] = (char*) calloc((size_t)len+1, sizeof(char));
	  for ( i = 0; i < len; i++ )
	    ((char **)var)[*nocc][i] = nml->line.lineac[nml->line.namitf+i];
	  *nocc += 1;
	}
    }
  else
    {
      fprintf(stderr, "Namelist parameter type %d unknown!\n", nml->line.nptype);
      return;
    }

  nml->line.nptype = ntyp;
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
    {
      printf("'");
      for ( j = 0; j < nout; j++ )
	printf("%c", ((char *)entry->ptr)[j]);
      printf("'");
    }
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

#define  MAX_WORD_LEN  256

void namelistRead(FILE *nmlfp, NAMELIST *nml)
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
  char namecx[MAX_WORD_LEN], *pnamecx = NULL;
  int nparam;

  nparam = nml->size;

  nml->line.namitl = MAX_LINE_LEN;

 L2000:

  getnite(nmlfp, nml);

  if ( nml->line.nptype == NML_NIX )
    {
      goto L3000;
    }
  else if ( nml->line.nptype == NML_KEYWORD )
    {
      memset(namecx, '\0', MAX_WORD_LEN);
      len = (size_t) (nml->line.namitl - nml->line.namitf + 1);

      if ( nml->line.lineac[nml->line.namitf] == '$' || 
	   nml->line.lineac[nml->line.namitf] == '&' )
	{
          if ( nml->line.namitl-nml->line.namitf > MAX_WORD_LEN ) goto L3000;

          if ( nml->line.namitf < nml->line.namitl)
	    pnamecx = &nml->line.linelc[nml->line.namitf+1];

          if ( strncmp(pnamecx, "select", 6) == 0 || 
	       strncmp(pnamecx, "params", 6) == 0 || 
	       strncmp(pnamecx, nml->name, strlen(nml->name)) == 0 ) goto L2000;

	  goto L3000;
        }

      if ( nml->line.namitl-nml->line.namitf >= MAX_WORD_LEN ) goto L3000;

      pnamecx = &nml->line.linelc[nml->line.namitf];
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
	      nml->line.nptype = NML_WORD;
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

      rdnlsgl(nml, nml->entry[j]->ptr, nml->entry[j]->type, (int)nml->entry[j]->size, &nml->entry[j]->occ);

      if ( nml->line.nptype != nml->entry[j]->type )
	{
          printf(" * value ignored for parameter <%s> %5d\n", namecx, nml->line.nptype);
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

  if ( nml->dis == 0 ) return;

  nml_print(nml, PRINT_MIN);
}

