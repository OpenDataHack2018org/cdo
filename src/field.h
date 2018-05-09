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

struct Field
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
};

struct RecordInfo
{
  short varID;
  short levelID;
  bool lconst;
};

// fieldmem.cc
void field_init(Field *field);
Field **field_malloc(const int vlistID, const int ptype);
Field **field_calloc(const int vlistID, const int ptype);
void field_free(Field **field, const int vlistID);

// field.cc
double fldfun(const Field &field, int function);
double fldrange(const Field &field);
double fldmin(const Field &field);
double fldmax(const Field &field);
double fldsum(const Field &field);
double fldmean(const Field &field);
double fldmeanw(const Field &field);
double fldavg(const Field &field);
double fldavgw(const Field &field);
double fldstd(const Field &field);
double fldstd1(const Field &field);
double fldvar(const Field &field);
double fldvar1(const Field &field);
double fldstdw(const Field &field);
double fldstd1w(const Field &field);
double fldvarw(const Field &field);
double fldvar1w(const Field &field);
double fldskew(const Field &field);
double fldkurt(const Field &field);

// ENS VALIDATION
double fldbrs(const Field &field);
double fldrank(const Field &field);
double fldroc(const Field &field);

double fldpctl(Field field, const double pn);
void fldunm(Field *field);
int fldhvs(Field *field, const size_t nlevels);

// fieldzon.cc
void zonfun(Field field1, Field *field2, const int function);
void zonmin(Field field1, Field *field2);
void zonmax(Field field1, Field *field2);
void zonrange(Field field1, Field *field2);
void zonsum(Field field1, Field *field2);
void zonavg(Field field1, Field *field2);
void zonmean(Field field1, Field *field2);
void zonstd(Field field1, Field *field2);
void zonstd1(Field field1, Field *field2);
void zonvar(Field field1, Field *field2);
void zonvar1(Field field1, Field *field2);
void zonpctl(Field field1, Field *field2, const int k);

/* fieldmer.cc */

void merfun(Field field1, Field *field2, const int function);
void mermin(Field field1, Field *field2);
void mermax(Field field1, Field *field2);
void merrange(Field field1, Field *field2);
void mersum(Field field1, Field *field2);
void meravgw(Field field1, Field *field2);
void mermeanw(Field field1, Field *field2);
void merstdw(Field field1, Field *field2);
void merstd1w(Field field1, Field *field2);
void mervarw(Field field1, Field *field2);
void mervar1w(Field field1, Field *field2);
void merpctl(Field field1, Field *field2, const int k);

void fldrms(Field field1, Field field2, Field *field3);

void varrms(Field field1, Field field2, Field *field3);

// fieldc.cc
void farcfun(Field *field, const double rconst, const int function);

void farcmul(Field *field, const double rconst);
void farcdiv(Field *field, const double rconst);
void farcadd(Field *field, const double rconst);
void farcsub(Field *field, const double rconst);

void farmod(Field *field, const double divisor);

void farinv(Field *field);
void farround(Field *field);

// field2.cc
void farfun(Field *field1, Field field2, int function);

void farcpy(Field *field1, Field field2);
void faradd(Field *field1, Field field2);
void farsum(Field *field1, Field field2);
void farsumw(Field *field1, Field field2, double w);
void farsumq(Field *field1, Field field2);
void farsumqw(Field *field1, Field field2, double w);
void farsumtr(Field *field1, Field field2, const double refval);
void farsub(Field *field1, Field field2);
void farmul(Field *field1, Field field2);
void fardiv(Field *field1, Field field2);
void farmin(Field *field1, Field field2);
void farmax(Field *field1, Field field2);
void farminidx(Field *field1, Field *field2, Field field3, int idx);
void farmaxidx(Field *field1, Field *field2, Field field3, int idx);
void farvar(Field *field1, Field field2, Field field3, int divisor);
void farstd(Field *field1, Field field2, Field field3, int divisor);
void farcvar(Field *field1, Field field2, int nsets, int divisor);
void farcstd(Field *field1, Field field2, int nsets, int divisor);
void farmoq(Field *field1, Field field2);
void farmoqw(Field *field1, Field field2, double w);
void faratan2(Field *field1, Field field2);
void farsetmiss(Field *field1, Field field2);

void farcount(Field *field1, Field field2);

// field2cplx.cc
void farfuncplx(Field *field1, Field field2, int function);

#endif /* FIELD_H */
