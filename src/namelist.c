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

int readline(FILE *fp, char *line, int len);

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
  int nfpm;
  char *prog;
};

struct PGMSTAT pgmstat;

#define MAXLINL 267

struct NAMELINE
{
  int nptype, namitf, namitl;
  char lineac[MAXLINL], lineuc[MAXLINL], linelc[MAXLINL];
};

struct NAMELINE nameline;


static void namelist_init(NAMELIST *namelist)
{
  static char func[] = "namelist_init";

  namelist->size = 0;
}


NAMELIST *namelistNew()
{
  static char func[] = "namelistNew";  
  NAMELIST *namelist;

  namelist = (NAMELIST *) malloc(sizeof(NAMELIST));

  namelist_init(namelist);

  return (namelist);
}


void namelistDelete(NAMELIST *namelist)
{
  static char func[] = "namelistDelete";  

  if ( namelist )
    {
      /*     if ( namelist->name )     free(namelist->name);*/
      free(namelist);
    }
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

          if      (   nameline.linelc[i] == ' ' ) ;
	  else if (   nameline.linelc[i] == '=' ) ;
          else if (   nameline.linelc[i] == ',' ) ;
          else if ( ((nameline.linelc[i] >= 'a')  &&
		     (nameline.linelc[i] <= 'z')) ||
	             (nameline.linelc[i] == '$')  ||
		     (nameline.linelc[i] == '&') )
	    {
	      nameline.nptype = func_nkw;
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
	      nameline.nptype = func_nir;
	      nameline.namitf = i;
	      for ( j = i+1; j < MAXLINL; j++)
		{
		  if ( nameline.linelc[j] >= '0' && nameline.linelc[j] <= '9' ) ;
	          else if ( nameline.linelc[j] == '+' ) ;
		  else if ( nameline.linelc[j] == '-' ) ;
         	  else if ( nameline.linelc[j] == '.' ) ;
		  else if ( nameline.linelc[j] == 'e' ) ;
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

  nameline.nptype = func_nix;
}


static void rdnlsgl(void *var, int ntyp, int nlen, int *nocc)
{
  if ( nameline.nptype == func_nir )
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
      else if ( ntyp == func_ntu )
	for (i=*nocc; i<newnocc; i++)
	  ((char *)var)[i] = nameline.lineuc[nameline.namitf+1+j++];
      else if ( ntyp == func_ntl )
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
	  ((char **)var)[*nocc] = (char*) calloc(len+1, sizeof(char));
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


static void rdnlscdo(void *var, int ntyp, int nlen, int nocc, int nlst, char *cpnam, int ife)
{
  char namout[24];
  int nout, j;

  if ( nlst == -1 ) return;

  strcpy(namout, cpnam);

  if ( ntyp == func_npr ) return;

  if ( ife == func_all )
    nout = nocc;
  else
    {
      nout = MAX(nocc, nlst);
      if ( nout == 0 ) return;
    }

  printf(" %-24s", namout);

  if      ( ntyp >= NML_TEXT )
    printf("'%s'", ((char *)var));
  else if ( ntyp == NML_WORD )
    for ( j = 0; j < nout; j++ )
      printf(" %s", ((char **)var)[j]);
  else if ( ntyp == NML_INT )
    for ( j = 0; j < nout; j++ )
      printf(" %d", ((int *)var)[j]);
  else if ( ntyp == NML_DOUBLE )
    for ( j = 0; j < nout; j++ )
      printf(" %g", ((double *)var)[j]);

  printf("\n");
}


static void rdnlout(int ipl, char *cn[], int nt[], int nl[], int nc[], int no[], int ife, void *vparam[])
{
  int i;

  for ( i = 0; i < ipl; i++ )
    rdnlscdo(vparam[i], nt[i], nl[i], nc[i], no[i], cn[i], ife);
}


void namelist(int nparam, char *cn[], int nt[], int nl[], int nc[], int no[], ...)
{
  /*
    cn  name
    nt  type
    nl  length
    nc  count
    no  list
  */
  const int numpar = 99;
  int clear = FALSE;
  int j, jj, match = -1, wordmatch = -1;
  size_t len;
  int index;
  char namecx[16], *pnamecx = NULL;
  void *vparam[99];
  va_list ap;

  if ( nparam > numpar )
    {
      fprintf(stderr, "Too much parameter in namelist!\n");
      return;
    }

  va_start(ap, no);

  for ( index = 0; index < nparam; index++ )
    vparam[index] = va_arg(ap, void *);

  va_end(ap);

  nameline.namitl = MAXLINL;

  if ( pgmstat.nfpm < 0 ) goto L3000;

 L2000:

  getnite();

  if ( nameline.nptype == func_nix )
    {
      goto L3000;
    }
  else if ( nameline.nptype == func_nkw )
    {
      memset(namecx, '\0', 16);
      len = nameline.namitl - nameline.namitf + 1;

      if ( nameline.lineac[nameline.namitf] == '$' || 
	   nameline.lineac[nameline.namitf] == '&' )
	{
          if ( nameline.namitl-nameline.namitf > 16 ) goto L3000;

          if ( nameline.namitf < nameline.namitl)
	    pnamecx = &nameline.linelc[nameline.namitf+1];
          if ( strncmp(pnamecx, "select", 6) == 0 || 
	       strncmp(pnamecx, "params", 6) == 0 || 
	       strncmp(pnamecx, pgmstat.prog, strlen(pgmstat.prog)) == 0 )
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
	  if ( strlen(cn[j]) == len )
	    if ( strncmp(pnamecx, cn[j], len) == 0 )
	      {
		jj = j;
		while ( nt[jj] == func_npr ) jj--;
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

	  rdnlout(nparam, cn, nt, nl, nc, no, func_all, vparam);

          if ( ! pgmstat.intract )
	    {
	      fprintf(stderr, "Namelist error!\n");
	      return;
	    }
        }
      else
	{
          if ( clear ) nc[match] = 0;
	  if ( nt[match] == NML_WORD )
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

      rdnlsgl(vparam[j], nt[j], nl[j], &nc[j]);

      if ( nameline.nptype != nt[j] )
	{
          printf(" * value ignored for parameter <%s> %5d\n", namecx, nameline.nptype);
	  printf(" * valid parameters and values specified so far are\n");

	  rdnlout(nparam, cn, nt, nl, nc, no, func_all, vparam);

          if ( ! pgmstat.intract )
	    {
	      fprintf(stderr, "Namelist error!\n");
	    }
        }
    }

  goto L2000;

 L3000:

  if ( pgmstat.pdis == func_not ) return;

  rdnlout(nparam, cn, nt, nl, nc, no, pgmstat.pdis, vparam);
}
