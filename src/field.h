#/*
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

#ifndef _FIELD_H
#define _FIELD_H


#define  ADD(x,y)  ( DBL_IS_EQUAL((x),missval1) || DBL_IS_EQUAL((y),missval2) ? missval1 : (x)+(y) )
#define  SUB(x,y)  ( DBL_IS_EQUAL((x),missval1) || DBL_IS_EQUAL((y),missval2) ? missval1 : (x)-(y) )
#define  MUL(x,y)  ( DBL_IS_EQUAL((x),0)||DBL_IS_EQUAL((y),0) ? 0 : DBL_IS_EQUAL((x),missval1) || DBL_IS_EQUAL((y),missval2) ? missval1 : (x)*(y) )
#define  DIV(x,y)  ( DBL_IS_EQUAL((x),missval1) || DBL_IS_EQUAL((y),missval2) || DBL_IS_EQUAL((y),0) ? missval1 : (x)/(y) )
#define  ROOT(x)   ( DBL_IS_EQUAL((x),missval1) || (x)<0 ? missval1 : sqrt(x) )


typedef struct {
  int      grid;
  int      zaxis;
  int      size;
  int      nsamp;
  int      nmiss;
  double   missval;
  double  *weight;
  double  *ptr;
}
FIELD;

/* field.c */

double fldfun(FIELD field, int function);
double fldmin(FIELD field);
double fldmax(FIELD field);
double fldsum(FIELD field);
double fldavg(FIELD field);
double fldmean(FIELD field);
double fldstd(FIELD field);
double fldvar(FIELD field);
/* RQ */
double fldpctl(FIELD field, int k);
/* QR */

/* fieldzon.c */

void zonfun(FIELD field1, FIELD *field2, int function);
void zonmin(FIELD field1, FIELD *field2);
void zonmax(FIELD field1, FIELD *field2);
void zonsum(FIELD field1, FIELD *field2);
void zonavg(FIELD field1, FIELD *field2);
void zonmean(FIELD field1, FIELD *field2);
void zonstd(FIELD field1, FIELD *field2);
void zonvar(FIELD field1, FIELD *field2);
/* RQ */
void zonpctl(FIELD field1, FIELD *field2, int k);
/* QR */

/* fieldmer.c */

void merfun(FIELD field1, FIELD *field2, int function);
void mermin(FIELD field1, FIELD *field2);
void mermax(FIELD field1, FIELD *field2);
void mersum(FIELD field1, FIELD *field2);
void meravg(FIELD field1, FIELD *field2);
void mermean(FIELD field1, FIELD *field2);
void merstd(FIELD field1, FIELD *field2);
void mervar(FIELD field1, FIELD *field2);
/* RQ */
void merpctl(FIELD field1, FIELD *field2, int k);
/* QR */

void fldrms(FIELD field1, FIELD field2, FIELD *field3);

void varrms(FIELD field1, FIELD field2, FIELD *field3);

/* fieldc.c */

void farcfun(FIELD *field, double rconst, int function);

void farcmul(FIELD *field, double rconst);
void farcdiv(FIELD *field, double rconst);
void farcadd(FIELD *field, double rconst);
void farcsub(FIELD *field, double rconst);

void farinv(FIELD *field);

/* field2.c */

void farfun(FIELD *field1, FIELD field2, int function);

void faradd(FIELD *field1, FIELD field2);
void farsum(FIELD *field1, FIELD field2);
void farsumq(FIELD *field1, FIELD field2);
void farsub(FIELD *field1, FIELD field2);
void farmul(FIELD *field1, FIELD field2);
void fardiv(FIELD *field1, FIELD field2);
void farmin(FIELD *field1, FIELD field2);
void farmax(FIELD *field1, FIELD field2);
void farvar(FIELD *field1, FIELD field2, FIELD field3);
void farstd(FIELD *field1, FIELD field2, FIELD field3);
void farcvar(FIELD *field1, FIELD field2, double rconst1);
void farcstd(FIELD *field1, FIELD field2, double rconst1);
void farmoq(FIELD *field1, FIELD field2);
void faratan2(FIELD *field1, FIELD field2);

#endif  /* _FIELD_H */
