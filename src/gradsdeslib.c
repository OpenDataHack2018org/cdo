#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "cdi.h"
/* #include "cdo.h" */
/* #include "cdo_int.h" */
#include "gradsdeslib.h"

extern int cdoDefaultDataType;

void dsets_init(dsets_t *dsets)
{
  int i;

  dsets->name[0]    = 0;
  dsets->dnam[0]    = 0;
  dsets->title[0]   = 0;
  dsets->bswap      = 0;
  dsets->fhdr       = 0;
  dsets->xyhdr      = 0;
  dsets->seqflg     = 0;
  dsets->yrflg      = 0;
  dsets->zrflg      = 0;
  dsets->tmplat     = 0;
  dsets->pa2mb      = 0;
  dsets->calendar   = 0;

  for ( i = 0; i < 5; ++i ) dsets->dnum[i]    = 0;
  

  dsets->nsets      = 0;
  dsets->mergelevel = 0;
  dsets->lgeoloc    = 0;
  dsets->lregion    = 0;
  dsets->lprojtype  = 0;
  dsets->lmetadata  = 0;

  for ( i = 0; i < MAX_DSETS; ++i )
    {
      dsets->obj[i].nx          = 0;
      dsets->obj[i].ny          = 0;
      dsets->obj[i].nz          = 0;
      dsets->obj[i].name        = NULL;
      dsets->obj[i].description = NULL;
      dsets->obj[i].units       = NULL;
      dsets->obj[i].title       = NULL;
      dsets->obj[i].time        = NULL;
      dsets->obj[i].dtype       = cdoDefaultDataType;
      dsets->obj[i].lscale      = 0;
      dsets->obj[i].loffset     = 0;
      dsets->obj[i].lmissval    = 0;
      dsets->obj[i].missval     = cdiInqMissval();
      dsets->obj[i].array       = NULL;   
    }
}

/*mf version
  convert all upper case alphabetic characters to lower case.
  The GrADS system is case insensitive, and assumes lower case
  internally in most cases. Does not turn to lower case if in "'s
*/
void lowcas (char *ch) {
int i;
int qflag=0;

  while (*ch!='\0' && *ch!='\n') {
    i = *ch;
    if(*ch == '\"' && qflag == 0 ) {
      qflag=1;
      } else if(*ch == '\"' && qflag == 1 ) {
	qflag=0;
      }
    if (i>64 && i<91 && qflag==0) {
      i+=32;
      *ch = i;
    } else if(i == 95) {
      *ch=i;
    }
    ch++;
  }
}

/* Compares two strings.  A match occurs if the leading
   blank-delimited words in the two strings match.  CR and NULL also
   serve as delimiters.                                               */

int cmpwrd (char *ch1, char *ch2) {

  while (*ch1==' '||*ch1=='\t') ch1++;  /* Advance past leading blanks.     */
  while (*ch2==' '||*ch2=='\t') ch2++;

  while (*ch1 == *ch2) {
    if (*ch1==' '||*ch1=='\t'||*ch1=='\0'||*ch1=='\n'||*ch1=='\r' ) return (1);
    ch1++; ch2++;
  }

  if ( (*ch1==' '||*ch1=='\t'||*ch1=='\0'||*ch1=='\n'||*ch1=='\r') &&
       (*ch2==' '||*ch2=='\t'||*ch2=='\0'||*ch2=='\n'||*ch2=='\r') ) return (1);
  return (0);
}


/* Parses a number in a character string.
   This routine will detect numbers of the form:
       nnnn
       -nnnn

   Args:    ch     - pointer to the number, in character form.
            val    - integer value returned
            return value  - address of 1st character past the
                            number parsed.  NULL if no number found
                            at pointer ch or if the number is an
                            invalid format.
             */

char *intprs (char *ch, int *val) {

int nflag,flag;

  nflag = 0;
  if (*ch=='-') { nflag = 1; ch++; }
  else if (*ch=='+') ch++;

  *val = 0;
  flag = 1;

  while (*ch>='0' && *ch<='9') {
    *val = *val*10 + (int)(*ch-'0');
    flag = 0;
    ch++;
  }

  if (flag) return (NULL);

  if (nflag) *val = -1 * *val;
  return (ch);
}

char *longprs (char *ch, long *val) {

int nflag,flag;

  nflag = 0;
  if (*ch=='-') { nflag = 1; ch++; }
  else if (*ch=='+') ch++;

  *val = 0;
  flag = 1;

  while (*ch>='0' && *ch<='9') {
    *val = *val*10 + (int)(*ch-'0');
    flag = 0;
    ch++;
  }

  if (flag) return (NULL);

  if (nflag) *val = -1 * *val;
  return (ch);
}

/* Moves a pointer to the start of the next blank-delimited word
   in a string.  If not found, NULL is returned.                     */

char * nxtwrd (char *ch) {

  while (*ch!=' '&&*ch!='\t') {                     /* Skip 1st word  */
    if (*ch == '\0' || *ch == '\n' || *ch == '\r') return (NULL);
    ch++;
  }
  while (*ch==' '||*ch=='\t') ch++;                 /* Find next word */
  if (*ch == '\0' || *ch == '\n' || *ch == '\r') return (NULL);
  return (ch);
}

/* Copies a string of a specified length, or when \0 or \n is hit.
   Trailing blanks are removed, and the output string is terminated
   with '\0'.                                                         */

void getstr (char *ch1, char *ch2, int len) {
char *ch;

  ch = ch1;
  while (len>0 && *ch2!='\n' && *ch2!='\0') {
    *ch1 = *ch2;
    len--;
    ch1++;  ch2++;
  }
  ch1--;
  while (ch1>=ch && *ch1==' ') ch1--;
  ch1++;
  *ch1 = '\0';
}

/* Copies a word of a specified length, or when \0 or \n or \r or ' ' is
   encountered.  The word is terminated with '\0'. ch2 is src, ch1 is dest */

void getwrd (char *ch1, char *ch2, int len) {
char *ch;

  ch = ch1;
  while (len>0 && *ch2!='\n' && *ch2!='\0' && *ch2!='\r' && *ch2!=' ' ) {
    *ch1 = *ch2;
    len--;
    ch1++;  ch2++;
  }
  *ch1 = '\0';
}

/* Converts strings to double */
char * getdbl(char *ch, double *val) {
  char * pos;
  double res;

  res = strtod(ch, &pos);
  if (pos==ch) {
    return NULL;
  } else {
    *val = res;
    return pos;
  }
}

/* Converts strings to double */
char * getflt(char *ch, float *val) {
char * pos;
  *val = (float)strtod(ch, &pos);
  if (pos==ch) {
    return NULL;
  } else {
    return pos;
  }
}

/* Expand file names prefixed with '^' from data descriptor
   files */

void fnmexp (char *out, char *in1, char *in2) {
char *pos, *ch, envv[20], *envr, CR=13;
int i,j;

  if (*in1=='$') {
    in1++;
    i = 0;
    while (*in1!='/' && *in1!='\0' && i<16) {
      envv[i] = *in1;
      i++; in1++;
    }
    envv[i] = '\0';
    envr = getenv(envv);
    if (envr) {
      i = 0; j = 0;
      while (*(envr+j)) {
        *(out+i) = *(envr+j);
        i++; j++;
      }
      /* handle CR for descriptor files created under MS Windows */
      while (*in1!='\0' && *in1!=' ' && *in1!='\n' && *in1!=CR) {
        *(out+i) = *in1;
        i++; in1++;
      }
      *(out+i) = '\0';
    }
    return;
  }
  ch = in2;
  pos=NULL;
  while (*ch!='\0' && *ch!=' ' && *ch!='\n') {
    if (*ch=='/') pos=ch;
    ch++;
  }
  if (pos) pos++;
  while (pos!=NULL && in2<pos) {
    *out = *in2;
    out++; in2++;
  }
  in1++;
  while (*in1!='\0' && *in1!=' ' && *in1!='\n' && *in1!=CR) {
    *out = *in1;
    out++; in1++;
  }
  *out = '\0';
}


/* Linear conversion routine for dimension conversions.               */

gadouble liconv (gadouble *vals, gadouble v) {
  return ( (*vals * v) + *(vals+1) );
}

/* Non-linear scaling routine for discrete levels.  Linear interp
   between levels.  Scaling beyond upper and lower bounds is
   linear based on the last and first grid spacing, respectively.
   In each case a pointer to a list of values is provided.  The
   list contains in its first element the number of values
   in the list.    */

/* Convert a grid value to a world coordinate value.
   This operation needs to be efficient, since it gets done
   very often.  */

gadouble gr2lev (gadouble *vals, gadouble gr) {
gaint i;
  if (gr<1.0) return ( *(vals+1) + (1.0-gr)*(*(vals+1)-*(vals+2)) );
  if (gr>*vals) {
    i = (gaint)(*vals+0.1);
    return ( *(vals+i) + (gr-*vals)*(*(vals+i)-*(vals+i-1)) );
  }
  i = (gaint)gr;
  return (*(vals+i)+((gr-(gadouble)i)*(*(vals+i+1)-*(vals+i))));
}

/* Convert from world coordinate value to grid value.  This operation
   is not set up to be efficient, under the assumption that it won't
   get done all that often.  */

gadouble lev2gr (gadouble *vals, gadouble lev) {
gaint i,num;
gadouble gr;
  num = (gaint)(*vals+0.1);
  for (i=1; i<num; i++) {
    if ( (lev >= *(vals+i) && lev <= *(vals+i+1)) ||
         (lev <= *(vals+i) && lev >= *(vals+i+1)) ) {
      gr = (gadouble)i + (lev - *(vals+i))/(*(vals+i+1) - *(vals+i));
      return (gr);
    }
  }
  if (*(vals+1)<*(vals+num)) {
    if (lev<*(vals+1)) {
      gr = 1.0 + ((lev-*(vals+1))/(*(vals+2)-*(vals+1)));
      return (gr);
    }
    gr = (gadouble)i + ((lev-*(vals+i))/(*(vals+i)-*(vals+i-1)));
    return (gr);
  } else {
    if (lev>*(vals+1)) {
      gr = 1.0 + ((lev-*(vals+1))/(*(vals+2)-*(vals+1)));
      return (gr);
    }
    gr = (gadouble)i + ((lev-*(vals+i))/(*(vals+i)-*(vals+i-1)));
    return (gr);
  }
}

/* Process linear scaling args */
/*
gaint deflin (char *ch, dsets_t *pfi, gaint dim, gaint flag) {
gadouble *vals,v1,v2;

  vals = (gadouble *)galloc(sizeof(gadouble)*6,"vals1");
  if (vals==NULL) return (-1);

  if ((ch = nxtwrd(ch))==NULL) goto err1;
  if (getdbl(ch,&v1)==NULL) goto err1;
  if (flag) v2 = 1.0;
  else {
    if ((ch = nxtwrd(ch))==NULL) goto err2;
    if (getdbl(ch,&v2)==NULL) goto err2;
  }
  if (dim!=3 && v2<=0.0) goto err2;
  *(vals)   = v2;
  *(vals+1) = v1 - v2;
  *(vals+2) = -999.9;
  pfi->grvals[dim] = vals;
  *(vals+4) = -1.0 * ( (v1-v2)/v2 );
  *(vals+3) = 1.0/v2;
  *(vals+5) = -999.9;
  pfi->abvals[dim] = vals+3;
  pfi->ab2gr[dim] = liconv;
  pfi->gr2ab[dim] = liconv;
  pfi->linear[dim] = 1;
  return (0);

err1:
  gaprnt (0,"Open Error:  Missing or invalid dimension");
  gaprnt (0," starting value\n");
  gree(vals,"f178");
  return (1);

err2:
  gaprnt (0,"Open Error:  Missing or invalid dimension");
  gaprnt (0," increment value\n");
  gree(vals,"179");
  return (1);
}
*/
/* Process levels values in def record */
/* Return codes:  -1 is memory allocation error, 1 is other error */
/*
gaint deflev (char *ch, char *rec, dsets_t *pfi, gaint dim) {
gadouble *vvs,*vals,v1;
gaint i;

  if (pfi->dnum[dim]==1) {
    i = deflin (ch, pfi, dim, 1);
    return (i);
  }

  vals = (gadouble *)galloc((pfi->dnum[dim]+5)*sizeof(gadouble),"vals2");
  if (vals==NULL) return (-1);

  vvs = vals;
  *vvs = (gadouble)pfi->dnum[dim];
  vvs++;
  for (i=0; i<pfi->dnum[dim]; i++) {
    if ( (ch = nxtwrd(ch))==NULL) {
      if (fgets(rec,256,descr)==NULL) goto err2;
      ch = rec;
      while (*ch==' ' || *ch=='\t') ch++;
      if (*ch=='\0' || *ch=='\n') goto err3;
    }
    if (getdbl(ch,&v1)==NULL) goto err1;
    *vvs = v1;
    vvs++;
  }
  *vvs = -999.9;
  pfi->abvals[dim] = vals;
  pfi->grvals[dim] = vals;
  pfi->ab2gr[dim] = lev2gr;
  pfi->gr2ab[dim] = gr2lev;
  pfi->linear[dim] = 0;
  return (0);

err1:
  gaprnt (0,"Open Error:  Invalid value in LEVELS data\n");
  gree(vals,"f180");
  return (1);

err2:
  gaprnt (0,"Open Error:  Unexpected EOF reading descriptor file\n");
  gaprnt (0,"   EOF occurred reading LEVELS values\n");
  gree(vals,"f181");
  return (1);

err3:
  gaprnt (0,"Open Error:  Blank Record found in LEVELS data\n");
  gree(vals,"f182");
  return (1);
}
*/

void gaprnt (int i, char *ch)
{
  printf ("%s",ch);
}



#undef  IsBigendian
#define IsBigendian()  ( u_byteorder.c[sizeof(long) - 1] )

int read_gradsdes(char *filename, dsets_t *pfi)
{
  /* IsBigendian returns 1 for big endian byte order */
  static union {unsigned long l; unsigned char c[sizeof(long)];} u_byteorder = {1};
  int status = 0;
  int reclen;
  int ichar;
  char rec[MAX_RECLEN], mrec[MAX_RECLEN];
  char *ch, *pos;
  FILE *descr;
  int i, ii, jj;
  int hdrb, trlb;
  int flgs[8];
  int BYTEORDER = IsBigendian();
  gadouble v1,v2,ev1,ev2,temp;
 
  hdrb = 0;
  trlb = 0;

  /* Try to open descriptor file */
  descr = fopen (filename, "r");
  if ( descr == NULL )
    {
      fprintf(stderr, "fopen failed on %s", filename);
      return (-1);
    }

  /* Copy descriptor file name into gafile structure */
  getwrd (pfi->dnam,filename,MAX_NAMELEN);

  /* initialize error flags */
  for (i=0;i<8;i++) flgs[i] = 1;

  /* Parse the data descriptor file */
  while ( fgets(rec, MAX_RECLEN, descr) != NULL )
    {

      /* Remove any leading blanks from rec */
      reclen = (int)strlen(rec);
      jj = 0;
      while ( jj<reclen && rec[0]==' ' ) 
	{
	  for (ii=0; ii<reclen; ii++) rec[ii] = rec[ii+1];
	  jj++;
	}
      /* replace newline with null at end of record */
      for (ichar = (int)strlen(rec) - 1 ;  ichar >= 0 ;  --ichar)
	{
	  if (rec[ichar] == '\n') {
	    rec[ichar] = '\0' ;
	    break ; 
	  }
	}
      /* Keep mixed case and lower case versions of rec handy */
      strcpy (mrec,rec);   
      lowcas(rec);
      
      if (!isalnum(mrec[0]))
	{
	  /* check if comment contains attribute metadata */
	  /*
	  if ((strncmp("*:attr",mrec,6)==0) || (strncmp("@",mrec,1)==0))
	    {
	      if ((ddfattr(mrec,pfi)) == -1) goto retrn;
	    }
	  */
	}
      else if (cmpwrd("byteswapped",rec))
	{
	  pfi->bswap = 1;
	}
      else if (cmpwrd("fileheader",rec))
	{
	  if ( (ch=nxtwrd(rec))==NULL )
	    {
	      gaprnt (1,"Description file warning: Missing fileheader length\n");
	    }
	  else
	    {
	      ch = longprs(ch,&(pfi->fhdr));
	      if (ch==NULL)
		{
		  gaprnt (1,"Fileheader record invalid\n");
		  pfi->fhdr = 0;
		}
	    }
	} 
      else if (cmpwrd("xyheader",rec)) 
	{
	  if ( (ch=nxtwrd(rec))==NULL )
	    {
	      gaprnt (1,"Description file warning: Missing xy grid header length\n");
	    } 
	  else 
	    {
	      ch = longprs(ch,&(pfi->xyhdr));
	      if (ch==NULL)
		{
		  gaprnt (1,"xy grid header length invalid\n");
		  pfi->xyhdr = 0;
		}
	      else
		{
		  pfi->xyhdr = pfi->xyhdr/4;
		}
	    }
	} 
      else if (cmpwrd("format",rec) || cmpwrd("options",rec))
	{
	  if ( (ch=nxtwrd(rec))==NULL )
	    {
	      gaprnt (1,"Description file warning: Missing options keyword\n");
	    }
	  else
	    {
	      while (ch!=NULL)
		{
		  if (cmpwrd("sequential",ch)) pfi->seqflg = 1;
		  else if (cmpwrd("yrev",ch)) pfi->yrflg = 1;
		  else if (cmpwrd("zrev",ch)) pfi->zrflg = 1;
		  else if (cmpwrd("template",ch)) pfi->tmplat = 1;
		  else if (cmpwrd("byteswapped",ch)) pfi->bswap = 1;
#if GRIB2
		  else if (cmpwrd("pascals",ch)) pfi->pa2mb = 1;
#endif
		  else if (cmpwrd("365_day_calendar",ch))
		    {
		      pfi->calendar=1;
		    }
		  else if (cmpwrd("big_endian",ch))
		    {
		      if (!BYTEORDER) pfi->bswap = 1;
		    }
		  else if (cmpwrd("little_endian",ch))
		    {
		      if (BYTEORDER) pfi->bswap = 1;
		    }
		  else {
		    gaprnt (0,"Open Error:  Data file type invalid\n");
		    goto err9;
		  }
		  ch = nxtwrd(ch);
		}
	    }
	}
      else if (cmpwrd("trailerbytes",rec))
	{
	  if ( (ch=nxtwrd(rec))==NULL )
	    {
	      gaprnt (1,"Trailerbytes record invalid\n");
	    }
	  else 
	    {
	      ch = intprs(ch,&trlb);
	      if (ch==NULL)
		{
		  gaprnt (1,"Trailerbytes record invalid\n");
		  trlb = 0;
		}
	      else 
		{
		  trlb = trlb/4;
		}
	    }
	}
      else if (cmpwrd("headerbytes",rec)|| cmpwrd("theader",rec))
	{
	  if ( (ch=nxtwrd(rec))==NULL )
	    {
	      gaprnt (1,"headerbytes/theader record invalid\n");
	    }
	  else 
	    {
	      ch = intprs(ch,&hdrb);
	      if (ch==NULL)
		{
		  gaprnt (1,"headerbytes/theader record invalid\n");
		  hdrb = 0;
		}
	      else
		{
		  hdrb = hdrb/4;
		}
	    }
	}
      else if (cmpwrd("title",rec))
	{
	  if ( (ch=nxtwrd(mrec))==NULL )
	    {
	      gaprnt (1,"Description file warning: Missing title string\n");
	    } 
	  else
	    {
	      getstr (pfi->title,ch,MAX_NAMELEN);
	      flgs[7] = 0;
	    }
	} 
      else if (cmpwrd("dset",rec))
	{
	  ch = nxtwrd(mrec);
	  if (ch==NULL) 
	    {
	      gaprnt (0,"Descriptor File Error:  Data file name is missing\n");
	      goto err9;
	    }
	  if (*ch=='^' || *ch=='$')
	    {
	      fnmexp (pfi->name,ch,filename);
	    } 
	  else 
	    {
	      getwrd (pfi->name,ch,MAX_NAMELEN);
	    }
	  flgs[5] = 0;
	}   
      else if (cmpwrd("undef",rec))
	{
	  ch = nxtwrd(mrec);
	  if (ch==NULL)
	    {
	      gaprnt (0,"Open Error:  Missing undef value\n");
	      goto err9;
	    }
      
	  pos = getdbl(ch,&(pfi->undef));
	  if (pos==NULL)
	    {
	      gaprnt (0,"Open Error:  Invalid undef value\n");
	      goto err9;
	    } 

	  pfi->ulow = fabs(pfi->undef/EPSILON);
	  pfi->uhi  = pfi->undef + pfi->ulow;
	  pfi->ulow = pfi->undef - pfi->ulow;
	  flgs[4] = 0;
	}
      /*
      else if (cmpwrd("xdef",rec))
	{
	  if (pfi->type == 2) continue;
	  if ( (ch = nxtwrd(rec)) == NULL) goto err1;
	  if ( (pos = intprs(ch,&(pfi->dnum[0])))==NULL) goto err1;
	  if (pfi->dnum[0]<1) 
	    {
	      sprintf(pout,"Warning: Invalid XDEF syntax in %s -- Changing size of X axis from %d to 1 \n",
		      pfi->dnam,pfi->dnum[0]);
	      gaprnt (1,pout);
	      pfi->dnum[0] = 1;
	    }
	  if (*pos!=' ') goto err1;
	  if ( (ch = nxtwrd(ch))==NULL) goto err2;
	  if (cmpwrd("linear",ch)) 
	    {
	      rc = deflin(ch, pfi, 0, 0);
	      if (rc==-1) goto err8; 
	      if (rc) goto err9;
	      v2 = *(pfi->grvals[0]);
	      v1 = *(pfi->grvals[0]+1) + v2;
	      temp = v1+((gadouble)(pfi->dnum[0]))*v2;
	      temp=temp-360.0;
	      if (fabs(temp-v1)<0.01) pfi->wrap = 1;
	    }
	  else if (cmpwrd("levels",ch))
	    {
	      rc = deflev (ch, rec, pfi, 0);
	      if (rc==-1)  goto err8; 
	      if (rc) goto err9;
	    } else goto err2;
	  flgs[0] = 0;
	} 
      else if (cmpwrd("ydef",rec))
	{
	  if (pfi->type == 2) continue;
	  if ( (ch = nxtwrd(rec)) == NULL) goto err1;
	  if ( (pos = intprs(ch,&(pfi->dnum[1])))==NULL) goto err1;
	  if (pfi->dnum[1]<1)
	    {
	      sprintf(pout,"Warning: Invalid YDEF syntax in %s -- Changing size of Y axis from %d to 1 \n",
		      pfi->dnam,pfi->dnum[1]);
	      gaprnt (1,pout);
	      pfi->dnum[1] = 1;
	    }
	  if (*pos!=' ') goto err1;
	  if ( (ch = nxtwrd(ch))==NULL) goto err2;
	  if (cmpwrd("linear",ch))
	    {
	      rc = deflin(ch, pfi, 1, 0);
	      if (rc==-1) goto err8; 
	      if (rc) goto err9;
	    }
	  else if (cmpwrd("levels",ch))
	    {
	      rc = deflev (ch, rec, pfi, 1);
	      if (rc==-1) goto err8;
	      if (rc) goto err9;
	    }
	  flgs[1] = 0;
	}
      else if (cmpwrd("zdef",rec))
	{
	  if (pfi->type == 2) continue;
	  if ( (ch = nxtwrd(rec)) == NULL) goto err1;
	  if ( (pos = intprs(ch,&(pfi->dnum[2])))==NULL) goto err1;
	  if (pfi->dnum[2]<1)
	    {
	      sprintf(pout,"Warning: Invalid ZDEF syntax in %s -- Changing size of Z axis from %d to 1 \n",
		      pfi->dnam,pfi->dnum[2]);
	      gaprnt (1,pout);
	      pfi->dnum[2] = 1;
	    }
	  if (*pos!=' ') goto err1;
	  if ( (ch = nxtwrd(ch))==NULL) goto err2;
	  if (cmpwrd("linear",ch))
	    {
	      rc = deflin(ch, pfi, 2, 0);
	      if (rc==-1) goto err8; 
	      if (rc) goto err9;
	    }
	  else if (cmpwrd("levels",ch))
	    {
	      rc = deflev (ch, rec, pfi, 2);
	      if (rc==-1) goto err8; 
	      if (rc) goto err9;
	    } else goto err2;
	  flgs[2] = 0;
	}
      else if (cmpwrd("tdef",rec))
	{
	  if ( (ch = nxtwrd(rec)) == NULL) goto err1;
	  if ( (pos = intprs(ch,&(pfi->dnum[3])))==NULL) goto err1;
	  if (pfi->dnum[3]<1)
	    {
	      sprintf(pout,"Warning: Invalid TDEF syntax in %s -- Changing size of T axis from %d to 1 \n",
		      pfi->dnam,pfi->dnum[3]);
	      gaprnt (1,pout);
	      pfi->dnum[3] = 1;
	    }
	  if (*pos!=' ') goto err1;
	  if ( (ch = nxtwrd(ch))==NULL) goto err2;
	  if (cmpwrd("linear",ch))
	    {
	      if ( (ch = nxtwrd(ch))==NULL) goto err3a_tdef;
	      tdef.yr = -1000;
	      tdef.mo = -1000;
	      tdef.dy = -1000;
	      if ( (pos = adtprs(ch,&tdef,&dt1))==NULL) goto err3b_tdef;
	      if (*pos!=' ' || dt1.yr == -1000 || dt1.mo == -1000.0 ||
		  dt1.dy == -1000) goto err3c_tdef;
	      if ( (ch = nxtwrd(ch))==NULL) goto err4a_tdef;
	      if ( (pos = rdtprs(ch,&dt2))==NULL) goto err4b_tdef;
	      v1 = (dt2.yr * 12) + dt2.mo;
	      v2 = (dt2.dy * 1440) + (dt2.hr * 60) + dt2.mn;
      *//* check if 0 dt *//*
	      if ( (v1 == 0) && (v2 == 0) ) goto err4c_tdef;  
	      if ((vals = (gadouble *)galloc(sizeof(gadouble)*8,"tvals5")) == NULL) goto err8; 
	      *(vals) = dt1.yr;
	      *(vals+1) = dt1.mo;
	      *(vals+2) = dt1.dy;
	      *(vals+3) = dt1.hr;
	      *(vals+4) = dt1.mn;
	      *(vals+5) = v1;
	      *(vals+6) = v2;
	      *(vals+7) = -999.9;
	      pfi->grvals[3] = vals;
	      pfi->abvals[3] = vals;
	      pfi->linear[3] = 1;
	    } else goto err2;
	  flgs[3] = 0;
	}
*/
    }
  
    /* Handle the chsub records.  time1, time2, then a string,  multiple times */
      /*
    } else if (cmpwrd("chsub",rec)) {
    }
      */

  fclose(descr);

  return (status);

 err1:
  gaprnt (0,"Open Error:  Missing or invalid dimension size.\n");
  goto err9;

 err2:
  gaprnt (0,"Open Error:  Missing or invalid dimension");
  gaprnt (0," scaling type\n");
  goto err9;

 err3a_tdef:
  gaprnt (0,"Open Error:  Start Time missing in tdef card");
  gaprnt (0," starting value\n");
  goto err9;

 err3b_tdef:
  gaprnt (0,"Open Error:  Invalid start time in tdef card");
  gaprnt (0," starting value\n");
  goto err9;

 err3c_tdef:
  gaprnt (0,"Open Error:  Missing or invalid dimension");
  gaprnt (0," starting value\n");
  goto err9;

 err3:
  gaprnt (0,"Open Error:  Missing or invalid dimension");
  gaprnt (0," starting value\n");
  goto err9;

 err4a_tdef:
  gaprnt (0,"Open Error:  Time increment missing in tdef\n");
  gaprnt (0," use 1 for single time data\n");
  goto err9;

 err4b_tdef:
  gaprnt (0,"Open Error:  Invalid time increment in tdef\n");
  gaprnt (0," use 1 for single time data\n");
  goto err9;

 err4c_tdef:
  gaprnt (0,"Open Error:  0 time increment in tdef\n");
  gaprnt (0," use 1 for single time data\n");
  goto err9;

 err5:
  gaprnt (0,"Open Error:  Missing or invalid variable");
  gaprnt (0," count\n");
  goto err9;

 err6:
  gaprnt (0,"Open Error:  Invalid variable record\n");
  goto err9;

 err6a:
  gaprnt (0,"Open Error:  Invalid x,y pair\n");
  goto err9;

 err7a: 
  gaprnt (0,"Open Error:  EOF occurred reading ensemble names\n");
  goto err9;

 err7b:
  gaprnt (0,"Open Error:  Blank record found in EDEF data\n");
  goto err9;

 err7c:
  gaprnt (0,"Open Error:  Invalid ensemble grib codes\n");
  goto err9;

 err7d:
  gaprnt (0,"Open Error:  Invalid ensemble name\n");
  goto err9;

 err7e:
  gaprnt (0,"Open Error:  Invalid ensemble record\n");
  goto err9;

 err8:
  gaprnt (0,"Open Error:  Memory allocation Error in gaddes.c\n");
  goto retrn;

 err9:
  gaprnt (0,"  --> The invalid description file record is: \n");
  gaprnt (0,"  --> ");
  gaprnt (0,mrec);
  gaprnt (0,"\n");

 retrn:
  gaprnt (0,"  The data file was not opened. \n");
  fclose (descr);
  return(1);
}
