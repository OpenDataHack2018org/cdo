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

static char pout[512];
FILE *descr;             /* File descriptor pointer */
int cal365 = 0;
int fullyear = -999;

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
  dsets->type       = 1;      /* Assume grid unless told otherwise */
  dsets->ncflg      = 0;      /* Assume not netcdf */

  dsets->pchsub1    = NULL;

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

/* Date/Time manipulation routines.  Note that these routines
   are not particularly efficient, thus Date/Time conversions
   should be kept to a minimum.                                      */

static gaint mosiz[13] = {0,31,28,31,30,31,30,31,31,30,31,30,31};

/* Test for leap year.  Rules are:

      Divisible by 4, it is a leap year, unless....
      Divisible by 100, it is not a leap year, unless...
      Divisible by 400, it is a leap year.                           */

gaint qleap (gaint year)  {
gaint i,y;

/*mf - disable if 365 day calendar mf*/

 if(/*mfcmn.*/cal365 == 1) return(0);

  y = year;

  i = y / 4;
  i = (i*4) - y;
  if (i!=0) return (0);

  i = y / 100;
  i = (i*100) - y;
  if (i!=0) return (1);

  i = y / 400;
  i = (i*400) - y;
  if (i!=0) return (0);

  return (1);


}

static char *mons[12] = {"jan","feb","mar","apr","may","jun",
			 "jul","aug","sep","oct","nov","dec"};

/* Parse an absolute date/time value.  Format is:

   12:00z 1jan 1989 (jan,feb,mar,apr,may,jun,jul,aug,sep,oct,nov,dec)

   Must have Z or Month abbrev, or value is invalid.  'def' contains
   higher order missing values (usually from tmin in pst).  Lower order
   values are defaulted to be: dy = 1, hr = 0, mn = 0.              */

char *adtprs (char *ch, struct dt *def, struct dt *dtim) {
gaint val,flag,i;
char *pos;
char monam[5];

  pos = ch;

  dtim->mn = 0;
  dtim->hr = 0;
  dtim->dy = 1;

  if (*ch>='0' && *ch<='9') {
  flag = 0;
    ch = intprs (ch,&val);
    if (*ch == ':' || tolower(*ch) == 'z') {
      if (val>23) {
        gaprnt (0,"Syntax Error:  Invalid Date/Time value.\n");
        sprintf (pout,"  Hour = %i -- greater than 23\n",val);
        gaprnt (0,pout);
        return (NULL);
      }
      dtim->hr = val;
      if (*ch == ':') {
        ch++;
        if (*ch>='0' && *ch<='9') {
          ch = intprs (ch,&val);
          if (val>59) {
            gaprnt (0,"Syntax Error:  Invalid Date/Time value.\n");
            sprintf (pout,"  Minute = %i -- greater than 59\n",val);
            gaprnt (0,pout);
            return (NULL);
          }
          if (tolower(*ch)!='z') {
            gaprnt (0,"Syntax Error:  Invalid Date/Time value.\n");
            gaprnt (0,"  'z' delimiter is missing \n");
            return (NULL);
          }
          dtim->mn = val;
          ch++;
          if (*ch>='0' && *ch<='9') ch = intprs (ch,&val);
          else val = def->dy;
        } else {
          gaprnt (0,"Syntax Error:  Invalid Date/Time value.\n");
          gaprnt (0,"  Missing minute value \n");
          return (NULL);
        }
      } else {
        ch++;
        if (*ch>='0' && *ch<='9') ch = intprs (ch,&val);
        else val = def->dy;
      }
    } else flag = 2;
    dtim->dy = val;
  } else flag = 1;

  monam[0] = tolower(*ch);
  monam[1] = tolower(*(ch+1));
  monam[2] = tolower(*(ch+2));
  monam[3] = '\0';

  i = 0;
  while (i<12 && !cmpwrd(monam,mons[i]) ) i++;
  i++;

  if (i==13) {
    if (flag==1) {
      gaprnt (0,"Syntax Error:  Invalid Date/Time value.\n");
      gaprnt (0,"  Expected month abbreviation, none found\n");
      return (NULL);
    }
    if (flag==2) {
      gaprnt (0,"Syntax Error:  Invalid Date/Time value.\n");
      gaprnt (0,"  Missing month abbreviation or 'z' delimiter\n");
      return (NULL);
    }
    dtim->mo = def->mo;
    dtim->yr = def->yr;
  } else {
    dtim->mo = i;
    ch+=3;
    /* parse year */
    if (*ch>='0' && *ch<='9') {
      /* use fullyear only if year 1 = 0001*/
      if(*(ch+2)>='0' && *(ch+2)<='9') {
	/*mfcmn.*/fullyear=1;
      } else {
	/*mfcmn.*/fullyear=0;
      }
      ch = intprs (ch,&val);
    } else {
      val = def->yr;
    }

    /* turn off setting of < 100 years to 1900 or 2000 */
    if(/*mfcmn.*/fullyear == 0) {
      if (val<50) val+=2000;
      else if (val<100) val+=1900;
    }
    dtim->yr = val;
  }

  i = mosiz[dtim->mo];
  if (dtim->mo==2 && qleap(dtim->yr)) i = 29;
  if (dtim->dy > i) {
    gaprnt (0,"Syntax Error:  Invalid Date/Time value.\n");
    sprintf (pout,"  Day = %i -- greater than %i \n",dtim->dy,i);
    gaprnt (0,pout);
    return (NULL);
  }
  return (ch);
}

/* Parse a relative date/time (offset).  Format is:

   nn (yr/mo/dy/hr/mn)

   Examples:  5mo
              1dy12hr
              etc.

   Missing values are filled in with 0s.                             */

char *rdtprs (char *ch, struct dt *dtim) {
gaint flag,val;
char *pos;
char id[3];

  pos = ch;

  dtim->yr = 0;
  dtim->mo = 0;
  dtim->dy = 0;
  dtim->hr = 0;
  dtim->mn = 0;

  flag = 1;

  while (*ch>='0' && *ch<='9') {
    flag = 0;
    ch = intprs(ch,&val);
    id[0] = *ch; id[1] = *(ch+1); id[2] = '\0';
    if (cmpwrd("yr",id)) dtim->yr = val;
    else if (cmpwrd("mo",id)) dtim->mo = val;
    else if (cmpwrd("dy",id)) dtim->dy = val;
    else if (cmpwrd("hr",id)) dtim->hr = val;
    else if (cmpwrd("mn",id)) dtim->mn = val;
    else {
      gaprnt (0,"Syntax Error:  Invalid Date/Time offset.\n");
      sprintf (pout,"  Expecting yr/mo/dy/hr/mn, found %s\n",id);
      gaprnt (0,pout);
      return (NULL);
    }
    ch+=2;
  }
  if (flag) {
    gaprnt (0,"Syntax Error:  Invalid Date/Time offset.\n");
    gaprnt (0,"  No offset value given\n");
    return (NULL);
  }
  return (ch);
}

/* Compares two strings.  A match occurs if the leading
   blank-delimited words in the two strings match.  CR and NULL also
   serve as delimiters.                                               */

gaint cmpwrd (char *ch1, char *ch2) {

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


/* Determines word length up to next delimiter */

gaint wrdlen (char *ch2) {
gaint len;
  len = 0;
  while (*ch2!='\n' && *ch2!='\0' && *ch2!=' ' && *ch2!='\t') {
    len++;
    ch2++;
  }
  return(len);
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

/* Process levels values in def record */
/* Return codes:  -1 is memory allocation error, 1 is other error */

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


/*  handle var name of the form longnm=>abbrv
    or just the abbrv with no long name */

gaint getvnm (struct gavar *pvar, char *mrec) {
gaint ib,i,j,k,len,flag;

  ib = 0;
  while (*(mrec+ib)==' ') ib++;

  if (*(mrec+ib)=='\0' || *(mrec+ib)=='\n') return(1);

  /* Scan for the '=>' string */
  len = 0;
  i = ib;
  flag = 0;

  while (1) {
    if (*(mrec+i)==' ' || *(mrec+i)=='\0' || *(mrec+i)=='\n') break;
    if (*(mrec+i)=='=' && *(mrec+i+1)=='>') {
      flag = 1;
      break;
    }
    len++ ; i++; 
  }

  if (flag) {
    for (j=ib; j<i; j++) {
      k = j-ib;
      pvar->longnm[k] = *(mrec+j); 
      /* substitute ~ for spaces in longname */
      if (pvar->longnm[k]=='~') pvar->longnm[k]=' '; 
    }
    pvar->longnm[len] = '\0';
    i+=2;
  } else {
    i = 0;
    pvar->longnm[0] = '\0';
  } 

  if (*(mrec+i)=='\n' || *(mrec+i)=='\0') return (1);

  getwrd (pvar->abbrv, mrec+i, 15);
  lowcas(pvar->abbrv);

  /* Check if 1st character is lower-case alphabetic */
  if (islower(*(pvar->abbrv))) return(0);
  else return (1);
}



#undef  IsBigendian
#define IsBigendian()  ( u_byteorder.c[sizeof(long) - 1] )

int read_gradsdes(char *filename, dsets_t *pfi)
{
  /* IsBigendian returns 1 for big endian byte order */
  static union {unsigned long l; unsigned char c[sizeof(long)];} u_byteorder = {1};
  struct gavar *pvar;
  struct dt tdef,dt1,dt2;
  struct gachsub *pchsub;
  int status = 0;
  int reclen;
  int ichar;
  char rec[MAX_RECLEN], mrec[MAX_RECLEN];
  char *ch, *pos;
  int i, j, ii, jj;
  gaint hdrb, trlb;
  gaint size=0,rc,len,flag,tim1,tim2;
  gaint flgs[8];
  int BYTEORDER = IsBigendian();
  gadouble *vals;
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
      /* Handle the chsub records.  time1, time2, then a string,  multiple times */
      else if (cmpwrd("chsub",rec))
	{
	  /* point to first block in chain */
	  pchsub = pfi->pchsub1;    
	  if (pchsub!=NULL)
	    {
	      while (pchsub->forw!=NULL) {
		pchsub = pchsub->forw;       /* advance to end of chain */
	      }
	    }
	  flag = 0;
	  ch = mrec;
	  while (1) 
	    {
	      if ( (ch=nxtwrd(ch)) == NULL ) break;
	      flag = 1;
	      if ( (ch = intprs(ch,&tim1)) == NULL) break;
	      if ( (ch=nxtwrd(ch)) == NULL ) break;
	      if (*ch=='*' && (*(ch+1)==' '||*(ch+1)=='\t')) tim2 = -99;
	      else if ( (ch = intprs(ch,&tim2)) == NULL) break;
	      if ( (ch=nxtwrd(ch)) == NULL ) break;
	      flag = 0;
	      if (pchsub) 
		{   /* chain exists */
		  pchsub->forw = (struct gachsub *)galloc(sizeof(struct gachsub),"chsubnew");
		  if (pchsub->forw==NULL) {
		    gaprnt(0,"Open Error: memory allocation failed for pchsub\n");
		    goto err8; 
		  }
		  pchsub = pchsub->forw;
		  pchsub->forw = NULL;
		} 
	      else 
		{        /* start a new chain */
		  pfi->pchsub1 = (struct gachsub *)galloc(sizeof(struct gachsub),"chsub1");
		  if (pfi->pchsub1==NULL)  {
		    gaprnt(0,"Open Error: memory allocation failed for pchsub1\n");
		    goto err8; 
		  }
		  pchsub = pfi->pchsub1;
		  pchsub->forw = NULL;
		}
	      len = wrdlen(ch);
	      if ((pchsub->ch = (char *)galloc(len+1,"chsubstr")) == NULL) goto err8;
	      getwrd(pchsub->ch,ch,len);
	      pchsub->t1 = tim1;
	      pchsub->t2 = tim2;
	    }
	  if (flag) 
	    {
	      gaprnt (1,"Description file warning: Invalid chsub record; Ignored\n");
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
	      /* check if 0 dt */
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
      else if (cmpwrd("vars",rec))
	{
	  if ( (ch = nxtwrd(rec)) == NULL) goto err5;
	  if ( (pos = intprs(ch,&(pfi->vnum)))==NULL) goto err5;
	  size = pfi->vnum * (sizeof(struct gavar) + 7 );
	  if ((pvar = (struct gavar *)galloc(size,"pvar2")) == NULL) goto err8;
	  pfi->pvar1 = pvar;
	  i = 0;
	  while (i<pfi->vnum)
	    {
	      /* initialize variables in the pvar structure */
	      pvar->offset = 0; 
	      pvar->recoff = 0;
	      pvar->ncvid = -999;
	      pvar->sdvid = -999;
	      pvar->levels = 0;
	      pvar->dfrm = 0;
	      pvar->var_t = 0;
	      pvar->scale = 1;
	      pvar->add = 0;  
	      pvar->undef= -9.99E33; 
	      pvar->vecpair = -999;
	      pvar->isu = 0;
	      pvar->isdvar = 0;
	      pvar->nvardims = 0; 

	      /* get the complete variable declaration */
	      if (fgets(rec,512,descr)==NULL) 
		{
		  gaprnt (0,"Open Error:  Unexpected EOF reading variables\n");
		  sprintf (pout, "Was expecting %i records.  Found %i.\n", pfi->vnum, i);
		  gaprnt (2,pout);
		  goto retrn;
		}
	      /* remove any leading blanks from rec */
	      reclen = strlen(rec);
	      jj = 0;
	      while (jj<reclen && rec[0]==' ')
		{
		  for (ii=0; ii<reclen; ii++) rec[ii] = rec[ii+1];
		  jj++;
		}
	      /* replace newline with null at end of record */
	      for (ichar = strlen(rec) - 1 ;  ichar >= 0 ;  --ichar)
		{
		  if (rec[ichar] == '\n')
		    {
		      rec[ichar] = '\0' ;
		      break ; 
		    }
		}
	      /* Keep mixed case and lower case versions of rec handy */
	      strcpy (mrec,rec);
	      lowcas(rec);
	      /* Allow comments between VARS and ENDVARS */
	      if (!isalnum(*(mrec)))
		{
		  /* Parse comment if it contains attribute metadata  */
		  /*
		  if ((strncmp("*:attr",mrec,6)==0) || (strncmp("@",mrec,1)==0)) {
		    if ((ddfattr(mrec,pfi)) == -1) goto retrn;
		    else continue;
		  }
		  else */continue; 
		}
	      if (cmpwrd("endvars",rec))
		{
		  gaprnt (0,"Open Error:  Unexpected ENDVARS record\n");
		  sprintf (pout, "Was expecting %i records.  Found %i.\n", pfi->vnum, i);
		  gaprnt (2,pout);
		  goto err9;
		}
	
	      /* get abbrv and full variable name if there */
	      if ((getvnm(pvar, mrec))!=0) goto err6;

	      /* parse the levels fields */
	      if ( (ch=nxtwrd(rec))==NULL) goto err6;
	      /* begin with 8th element of units aray for levels values */
	      for (j=0;j<16;j++) pvar->units[j] = -999;
	      j = 8;          
	      while (1) 
		{
		  if (j==8) {
		    /* first element is num levels */
		    if ((ch=intprs(ch,&(pvar->levels)))==NULL) goto err6;      
		  }
		  else {
		    /* remaining elements are grib2 level codes */
		    if ((ch=getdbl(ch,&(pvar->units[j-1])))==NULL) goto err6;  
		  }
		  /* advance through comma-delimited list of levels args */
		  while (*ch==' ') ch++;
		  if (*ch=='\0' || *ch=='\n') goto err6;
		  if (*ch!=',') break;
		  ch++;
		  while (*ch==',') { ch++; j++;}  /* advance past back to back commas */
		  while (*ch==' ') ch++;
		  if (*ch=='\0' || *ch=='\n') goto err6;
		  j++;
		  if (j>15) goto err6;
		}

	      /* parse the units fields; begin with 0th element for variable units */
	      j = 0;
	      pvar->nvardims=0;
	      while (1)
		{
		  if (*ch=='x'||*ch=='y'||*ch=='z'||*ch=='t'||*ch=='e')
		    { 
		      if (*(ch+1)!=',' && *(ch+1)!=' ') goto err6;
		      if (*ch=='x') { pvar->units[j] = -100; pvar->nvardims++; }
		      if (*ch=='y') { pvar->units[j] = -101; pvar->nvardims++; }
		      if (*ch=='z') { pvar->units[j] = -102; pvar->nvardims++; }
		      if (*ch=='t') { pvar->units[j] = -103; pvar->nvardims++; }
		      if (*ch=='e') { pvar->units[j] = -104; pvar->nvardims++; }
		      ch++;
		    } 
		  else 
		    {
		      if ( (ch=getdbl(ch,&(pvar->units[j])))==NULL ) goto err6;
		      /* no negative array indices for ncflag files */
		      if ((pfi->ncflg) && (pvar->units[j] < 0))  goto err6;   
		    }
		  while (*ch==' ') ch++;
		  if (*ch=='\0' || *ch=='\n') goto err6;
		  if (*ch!=',') break;
		  ch++;
		  while (*ch==' ') ch++;
		  if (*ch=='\0' || *ch=='\n') goto err6;
		  j++;
		  if (j>8) goto err6;
		}

	      /* parse the variable description */
	      getstr (pvar->varnm,mrec+(ch-rec),127);


	      /* var_t is for data files with dimension sequence: X, Y, Z, T, V */
	      if ((pvar->units[0]==-1) && 
		  (pvar->units[1]==20)) 
		pvar->var_t = 1;

	      /* non-float data types */
	      if ((pvar->units[0]==-1) && 
		  (pvar->units[1]==40))
		{

		  if (pvar->units[2]== 1) pvar->dfrm = 1;
		  if (pvar->units[2]== 2)
		    {
		      pvar->dfrm = 2;
		      if (pvar->units[3]==-1) pvar->dfrm = -2;
		    }
		  if (pvar->units[2]== 4) pvar->dfrm = 4;
		}

	      i++; pvar++;
	    }

	  /* Get ENDVARS statement and any additional comments */
	  if (fgets(rec,512,descr)==NULL) {
	    gaprnt (0,"Open Error:  Missing ENDVARS statement.\n");
	    goto retrn;
	  }
	  /* Remove any leading blanks from rec */
	  reclen = strlen(rec);
	  jj = 0;
	  while (jj<reclen && rec[0]==' ') {
	    for (ii=0; ii<reclen; ii++) rec[ii] = rec[ii+1];
	    jj++;
	  }
	  /* replace newline with null at end of record */
	  for (ichar = strlen(rec) - 1 ;  ichar >= 0 ;  --ichar) {
	    if (rec[ichar] == '\n') {
	      rec[ichar] = '\0' ;
	      break ; 
	    }
	  }
	  /* Keep mixed case and lower case versions handy */
	  strcpy (mrec,rec);
	  lowcas(rec);
	  while (!cmpwrd("endvars",rec)) 
	    {
	      /* see if it's an attribute comment */
	      if (!isalnum(*(mrec))) {
		/*
		if ((strncmp("*:attr",mrec,6)==0) || (strncmp("@",mrec,1)==0)) {
		  if ((ddfattr(mrec,pfi)) == -1) goto retrn;
		}
		*/
	      }
	      else {
		sprintf(pout,"Open Error:  Looking for \"endvars\", found \"%s\" instead.\n",rec);
		gaprnt (0,pout);
		goto err9;
	      }
	      /* get a new record */
	      if (fgets(rec,512,descr)==NULL) {
		gaprnt (0,"Open Error:  Missing ENDVARS statement.\n");
		goto retrn;
	      }
	      /* Remove any leading blanks from new record */
	      reclen = strlen(rec);
	      jj = 0;
	      while (jj<reclen && rec[0]==' ') {
		for (ii=0; ii<reclen; ii++) rec[ii] = rec[ii+1];
		jj++;
	      }
	      /* replace newline with null at end of record */
	      for (ichar = strlen(rec) - 1 ;  ichar >= 0 ;  --ichar) {
		if (rec[ichar] == '\n') {
		  rec[ichar] = '\0' ;
		  break ; 
		}
	      }
	      /* Keep mixed case and lower case versions handy */
	      strcpy (mrec,rec);
	      lowcas(rec);
	    }
	  /* vars block parsed without error */
	  flgs[6] = 0;

	} 
      else
	{
	  /* parse error of .ctl file */
	  gaprnt (0,"Open Error:  Unknown keyword in description file\n");
	  goto err9;
	}
    }

  /* Done scanning!
     Check if scanned stuff makes sense, and then set things up correctly */


  /* Make sure there are no conflicting options and data types */
  pvar=pfi->pvar1;
  for (j=1; j<=pfi->vnum; j++) {
    if (pvar->units[0]==-1 && pvar->units[1]==20) {
      if (pfi->tmplat) {
	gaprnt(0,"Open Error: Variables with transposed VAR-T dimensions cannot be templated together\n");
	err=1;
      }
      if (hdrb>0) {
	gaprnt(0,"Open Error: Variables with transposed VAR-T dimensions are incompatible with time headers\n");
	err=1;
      }
      if (trlb>0) {
	gaprnt(0,"Open Error: Variables with transposed VAR-T dimensions are incompatible with TRAILERBYTES\n");
	err=1;
      }
    }
    pvar++;
  }
  if (err) goto retrn;


  /* Figure out locations of variables within a time group */
  pvar = pfi->pvar1;

  /* Grid data */
  if (pfi->type==1) {
    pfi->gsiz = pfi->dnum[0] * pfi->dnum[1];
    if (pfi->ppflag) pfi->gsiz = pfi->ppisiz * pfi->ppjsiz;
    /* add the XY header to gsiz */
    if (pfi->xyhdr) {
      if (pvar->dfrm == 1) {
	pfi->xyhdr = pfi->xyhdr*4/1;          
      } 
      else if (pvar->dfrm ==  2 || pvar->dfrm == -2 ) {
	pfi->xyhdr = pfi->xyhdr*4/2;
      } 
      pfi->gsiz = pfi->gsiz + pfi->xyhdr;
    }

    /* adjust the size of hdrb and trlb for non-float data */
    if (pvar->dfrm == 1) {
      hdrb = hdrb*4/1;
      trlb = trlb*4/1;
    } 
    else if (pvar->dfrm == 2 || pvar->dfrm == -2 ) {
      hdrb = hdrb*4/2;
      trlb = trlb*4/2;
    } 
    
    if (pfi->seqflg) {
      /* pad the grid size with 2 4-byte chunks */
      if (pvar->dfrm == 1) {
	pfi->gsiz += 8;
      } 
      else if (pvar->dfrm == 2 || pvar->dfrm == -2 ) {
	pfi->gsiz += 4;
      } 
      else {
	pfi->gsiz += 2;             
      }
      /* pad the header with 2 4-byte chunks*/
      if (hdrb>0) {
	if (pvar->dfrm == 1) {
	  hdrb = hdrb + 8;
	} 
	else if (pvar->dfrm == 2 || pvar->dfrm == -2 ) {
	  hdrb = hdrb + 4;
	} 
	else {
	  hdrb += 2; 
	}
      }
      /* how far we have to go into the file before getting to 1st var */
      if (pvar->dfrm == 1) {
	pvar->offset = 4+hdrb;
	acum = 4+hdrb;
      } 
      else if (pvar->dfrm == 2 || pvar->dfrm == -2 ) {
	pvar->offset = 2+hdrb;
	acum = 2+hdrb;
      } 
      else {
	pvar->offset = 1+hdrb;
	acum = 1+hdrb;
      } 
    }
    else {
      /* how far we have to go into the file before getting to 1st var */
      pvar->offset = hdrb;
      acum = hdrb;
    }

    levs = pvar->levels;
    if (levs==0) levs=1;
    pvar->recoff = 0;
    recacm = 0;
    pvar++;
    acumvz=acum;

    for (i=1; i<pfi->vnum; i++) {
      if (pvar->var_t) {   
	acum = acum + levs*(pfi->gsiz)*(pfi->dnum[3]); 
      } else {                              
	acum = acum + (levs*pfi->gsiz);
	acumstride = acum ;
      }
      recacm += levs;
      pvar->offset = acum;
      pvar->recoff = recacm;
      levs = pvar->levels;
      if (levs==0) levs=1;
      pvar++;
    }

    recacm += levs;

    /* last variable */
    acum = acum + (levs*pfi->gsiz);

    pfi->tsiz = acum;
    pfi->trecs = recacm;
    if (pfi->seqflg) pfi->tsiz-=1;
    pfi->tsiz += trlb;
    
  } 
  else {
    fprintf(stderr, "Grid data type unsupported!");
    return (-1);
  }

/* set the global calendar and check if we are trying to change with a new file...
   we do this here to set the calandar for templating */

  if (mfcmn.cal365<0) {
    mfcmn.cal365=pfi->calendar;
  } else {
    if (pfi->calendar != mfcmn.cal365) {
      gaprnt(0,"Attempt to change the global calendar...\n");
      if (mfcmn.cal365) {
	gaprnt(0,"The calendar is NOW 365 DAYS and you attempted to open a standard calendar file\n");
      } else {
	gaprnt(0,"The calendar is NOW STANDARD and you attempted to open a 365-day calendar file\n");
      }
      goto retrn;
    }
  }



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
