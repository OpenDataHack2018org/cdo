/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2018 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#ifndef FIELD_H
#define FIELD_H

#include "compare.h"
#include "array.h"

double varToStd(double rvar, double missval);

enum field_flag
{
  FIELD_NONE = 1,
  FIELD_PTR = 2,
  FIELD_WGT = 4,
  FIELD_PTR2 = 8,
  FIELD_FLT = 16,
  FIELD_ALL = FIELD_PTR | FIELD_WGT
};

typedef struct
{
  int fpeRaised;
  int nwpv;  // number of words per value; real:1  complex:2
  int memtype;
  int grid;
  int zaxis;
  size_t size;
  size_t nsamp;
  size_t nmiss;
  size_t nmiss2;
  double missval;
  double *weight;
  double *ptr;
  float *ptrf;
  void *ptr2;
} field_type;

typedef struct
{
  short varID;
  short levelID;
  bool lconst;
} recinfo_type;

/* fieldmem.cc */

void field_init(field_type *field);
field_type **field_malloc(const int vlistID, const int ptype);
field_type **field_calloc(const int vlistID, const int ptype);
void field_free(field_type **field, const int vlistID);

/* field.cc */

double fldfun(field_type field, int function);
double fldrange(field_type field);
double fldmin(field_type field);
double fldmax(field_type field);
double fldsum(field_type field);
double fldavg(field_type field);
double fldmean(field_type field);
double fldstd(field_type field);
double fldstd1(field_type field);
double fldvar(field_type field);
double fldvar1(field_type field);
double fldavgw(field_type field);
double fldmeanw(field_type field);
double fldstdw(field_type field);
double fldstd1w(field_type field);
double fldvarw(field_type field);
double fldvar1w(field_type field);
double fldpctl(field_type field, const double pn);
void fldunm(field_type *field);
int fldhvs(field_type *field, const size_t nlevels);
double fldskew(field_type field);
double fldkurt(field_type field);

/* ENS VALIDATION */
double fldbrs(field_type field);
double fldrank(field_type field);
double fldroc(field_type field);

/* fieldzon.cc */

void zonfun(field_type field1, field_type *field2, const int function);
void zonmin(field_type field1, field_type *field2);
void zonmax(field_type field1, field_type *field2);
void zonrange(field_type field1, field_type *field2);
void zonsum(field_type field1, field_type *field2);
void zonavg(field_type field1, field_type *field2);
void zonmean(field_type field1, field_type *field2);
void zonstd(field_type field1, field_type *field2);
void zonstd1(field_type field1, field_type *field2);
void zonvar(field_type field1, field_type *field2);
void zonvar1(field_type field1, field_type *field2);
void zonpctl(field_type field1, field_type *field2, const int k);

/* fieldmer.cc */

void merfun(field_type field1, field_type *field2, const int function);
void mermin(field_type field1, field_type *field2);
void mermax(field_type field1, field_type *field2);
void merrange(field_type field1, field_type *field2);
void mersum(field_type field1, field_type *field2);
void meravgw(field_type field1, field_type *field2);
void mermeanw(field_type field1, field_type *field2);
void merstdw(field_type field1, field_type *field2);
void merstd1w(field_type field1, field_type *field2);
void mervarw(field_type field1, field_type *field2);
void mervar1w(field_type field1, field_type *field2);
void merpctl(field_type field1, field_type *field2, const int k);

void fldrms(field_type field1, field_type field2, field_type *field3);

void varrms(field_type field1, field_type field2, field_type *field3);

/* fieldc.cc */

void farcfun(field_type *field, const double rconst, const int function);

void farcmul(field_type *field, const double rconst);
void farcdiv(field_type *field, const double rconst);
void farcadd(field_type *field, const double rconst);
void farcsub(field_type *field, const double rconst);

void farmod(field_type *field, const double divisor);

void farinv(field_type *field);
void farround(field_type *field);

/* field2.cc */

void farfun(field_type *field1, field_type field2, int function);

void farcpy(field_type *field1, field_type field2);
void faradd(field_type *field1, field_type field2);
void farsum(field_type *field1, field_type field2);
void farsumw(field_type *field1, field_type field2, double w);
void farsumq(field_type *field1, field_type field2);
void farsumqw(field_type *field1, field_type field2, double w);
void farsumtr(field_type *field1, field_type field2, const double refval);
void farsub(field_type *field1, field_type field2);
void farmul(field_type *field1, field_type field2);
void fardiv(field_type *field1, field_type field2);
void farmin(field_type *field1, field_type field2);
void farmax(field_type *field1, field_type field2);
void farvar(field_type *field1, field_type field2, field_type field3, int divisor);
void farstd(field_type *field1, field_type field2, field_type field3, int divisor);
void farcvar(field_type *field1, field_type field2, int nsets, int divisor);
void farcstd(field_type *field1, field_type field2, int nsets, int divisor);
void farmoq(field_type *field1, field_type field2);
void farmoqw(field_type *field1, field_type field2, double w);
void faratan2(field_type *field1, field_type field2);
void farsetmiss(field_type *field1, field_type field2);

void farcount(field_type *field1, field_type field2);

#endif /* FIELD_H */
