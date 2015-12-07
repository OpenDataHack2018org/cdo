/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2015 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#ifndef _FIELD_H
#define _FIELD_H

double var_to_std(double rvar, double missval);

enum field_flag {
  FIELD_NONE  =  1,
  FIELD_PTR   =  2,
  FIELD_WGT   =  4,
  FIELD_PTR2  =  8,
  FIELD_FLT   = 16,
  FIELD_ALL   = FIELD_PTR | FIELD_WGT
};


#define  MADDMN(x,y)  (DBL_IS_EQUAL((x),missval1) || DBL_IS_EQUAL((y),missval2) ? missval1 : (x)+(y))
#define  MSUBMN(x,y)  (DBL_IS_EQUAL((x),missval1) || DBL_IS_EQUAL((y),missval2) ? missval1 : (x)-(y))
#define  MMULMN(x,y)  (DBL_IS_EQUAL((x),missval1) || DBL_IS_EQUAL((y),missval2) ? missval1 : (x)*(y))
#define  MDIVMN(x,y)  (DBL_IS_EQUAL((x),missval1) || DBL_IS_EQUAL((y),missval2) || DBL_IS_EQUAL((y),0.) ? missval1 : (x)/(y))
#define  MPOWMN(x,y)  (DBL_IS_EQUAL((x),missval1) || DBL_IS_EQUAL((y),missval2) ? missval1 : pow((x),(y)))
#define  MSQRTMN(x)   (DBL_IS_EQUAL((x),missval1) || (x)<0 ? missval1 : sqrt(x))


#define  ADD(x,y)  ((x)+(y))
#define  SUB(x,y)  ((x)-(y))
#define  MUL(x,y)  ((x)*(y))
#define  DIV(x,y)  (IS_EQUAL((y),0.) ? missval1 : (x)/(y))
#define  POW(x,y)  pow((x),(y))
#define  SQRT(x)   sqrt(x)


#define  ADDM(x,y)  (IS_EQUAL((x),missval1) || IS_EQUAL((y),missval2) ? missval1 : (x)+(y))
#define  SUBM(x,y)  (IS_EQUAL((x),missval1) || IS_EQUAL((y),missval2) ? missval1 : (x)-(y))
#define  MULM(x,y)  (IS_EQUAL((x),missval1) || IS_EQUAL((y),missval2) ? missval1 : (x)*(y))
#define  DIVM(x,y)  (IS_EQUAL((x),missval1) || IS_EQUAL((y),missval2) || IS_EQUAL((y),0.) ? missval1 : (x)/(y))
#define  POWM(x,y)  (IS_EQUAL((x),missval1) || IS_EQUAL((y),missval2) ? missval1 : pow((x),(y)))
#define  SQRTM(x)   (IS_EQUAL((x),missval1) || (x)<0 ? missval1 : sqrt(x))


#define  ADDMN(x,y)  FADDMN(x, y, missval1, missval2)
#define  SUBMN(x,y)  FSUBMN(x, y, missval1, missval2)
#define  MULMN(x,y)  FMULMN(x, y, missval1, missval2)
#define  DIVMN(x,y)  FDIVMN(x, y, missval1, missval2)
#define  POWMN(x,y)  FPOWMN(x, y, missval1, missval2)
#define  SQRTMN(x)   FSQRTMN(x, missval1)


static inline
double FADDMN(double x, double y, double missval1, double missval2) { return MADDMN(x,y);}
static inline
double FSUBMN(double x, double y, double missval1, double missval2) { return MSUBMN(x, y);}
static inline
double FMULMN(double x, double y, double missval1, double missval2) { return MMULMN(x, y);}
static inline
double FDIVMN(double x, double y, double missval1, double missval2) { return MDIVMN(x, y);}
static inline
double FPOWMN(double x, double y, double missval1, double missval2) { return MPOWMN(x, y);}
static inline
double FSQRTMN(double x, double missval1) { return MSQRTMN(x);}

typedef struct {
  int      nwpv; // number of words per value; real:1  complex:2
  int      grid;
  int      zaxis;
  size_t   size;
  int      nsamp;
  int      nmiss;
  int      nmiss2;
  int      memtype;
  double   missval;
  double  *weight;
  double  *ptr;
  float   *ptrf;
  void    *ptr2;
}
field_t;


/* fieldmem.c */

void      field_init(field_t *field);
field_t **field_malloc(const int vlistID, const int ptype);
field_t **field_calloc(const int vlistID, const int ptype);
void      field_free(field_t **field, const int vlistID);

/* field.c */

double fldfun(field_t field, int function);
double fldmin(field_t field);
double fldmax(field_t field);
double fldsum(field_t field);
double fldavg(field_t field);
double fldmean(field_t field);
double fldstd(field_t field);
double fldstd1(field_t field);
double fldvar(field_t field);
double fldvar1(field_t field);
double fldpctl(field_t field, const double pn);
void   fldunm(field_t *field);
int    fldhvs(field_t *field, const size_t nlevels);

/* ENS VALIDATION */
double fldcrps(field_t field);
double fldbrs(field_t field);
double fldrank(field_t field);
double fldroc(field_t field);

/* fieldzon.c */

void zonfun(field_t field1, field_t *field2, const int function);
void zonmin(field_t field1, field_t *field2);
void zonmax(field_t field1, field_t *field2);
void zonrange(field_t field1, field_t *field2);
void zonsum(field_t field1, field_t *field2);
void zonavg(field_t field1, field_t *field2);
void zonmean(field_t field1, field_t *field2);
void zonstd(field_t field1, field_t *field2);
void zonstd1(field_t field1, field_t *field2);
void zonvar(field_t field1, field_t *field2);
void zonvar1(field_t field1, field_t *field2);
void zonpctl(field_t field1, field_t *field2, const int k);

/* fieldmer.c */

void merfun(field_t field1, field_t *field2, const int function);
void mermin(field_t field1, field_t *field2);
void mermax(field_t field1, field_t *field2);
void mersum(field_t field1, field_t *field2);
void meravg(field_t field1, field_t *field2);
void mermean(field_t field1, field_t *field2);
void merstd(field_t field1, field_t *field2);
void merstd1(field_t field1, field_t *field2);
void mervar(field_t field1, field_t *field2);
void mervar1(field_t field1, field_t *field2);
void merpctl(field_t field1, field_t *field2, const int k);

void fldrms(field_t field1, field_t field2, field_t *field3);

void varrms(field_t field1, field_t field2, field_t *field3);

/* fieldc.c */

void farcfun(field_t *field, const double rconst, const int function);

void farcmul(field_t *field, const double rconst);
void farcdiv(field_t *field, const double rconst);
void farcadd(field_t *field, const double rconst);
void farcsub(field_t *field, const double rconst);

void farmod(field_t *field, const double divisor);

void farinv(field_t *field);
void farround(field_t *field);

/* field2.c */

void farfun(field_t *field1, field_t field2, int function);

void farcpy(field_t *field1, field_t field2);
void faradd(field_t *field1, field_t field2);
void farsum(field_t *field1, field_t field2);
void farsumw(field_t *field1, field_t field2, double w);
void farsumq(field_t *field1, field_t field2);
void farsumqw(field_t *field1, field_t field2, double w);
void farsumtr(field_t *field1, field_t field2, const double refval);
void farsub(field_t *field1, field_t field2);
void farmul(field_t *field1, field_t field2);
void fardiv(field_t *field1, field_t field2);
void farmin(field_t *field1, field_t field2);
void farmax(field_t *field1, field_t field2);
void farvar(field_t *field1, field_t field2, field_t field3, int divisor);
void farstd(field_t *field1, field_t field2, field_t field3, int divisor);
void farcvar(field_t *field1, field_t field2, int nsets, int divisor);
void farcstd(field_t *field1, field_t field2, int nsets, int divisor);
void farmoq(field_t *field1, field_t field2);
void farmoqw(field_t *field1, field_t field2, double w);
void faratan2(field_t *field1, field_t field2);

void farcount(field_t *field1, field_t field2);

#endif  /* _FIELD_H */
