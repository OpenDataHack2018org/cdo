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
#ifndef  ARRAY_H
#define  ARRAY_H

#define  MADDMN(x,y)  (DBL_IS_EQUAL((x),missval1) || DBL_IS_EQUAL((y),missval2) ? missval1 : (x)+(y))
#define  MSUBMN(x,y)  (DBL_IS_EQUAL((x),missval1) || DBL_IS_EQUAL((y),missval2) ? missval1 : (x)-(y))
#define  MMULMN(x,y)  (DBL_IS_EQUAL((x),0.)||DBL_IS_EQUAL((y),0.) ? 0 : DBL_IS_EQUAL((x),missval1) || DBL_IS_EQUAL((y),missval2) ? missval1 : (x)*(y))
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
#define  MULM(x,y)  (IS_EQUAL((x),0.)||IS_EQUAL((y),0.) ? 0 : IS_EQUAL((x),missval1) || IS_EQUAL((y),missval2) ? missval1 : (x)*(y))
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


const char *fpe_errstr(int fpeRaised);

int array_minmaxsum_val(size_t len, const double *array, double *rmin, double *rmax, double *rsum);
int array_minmaxmean_val(size_t len, const double *array, double *rmin, double *rmax, double *rmean);

int array_add_array(size_t len, double *restrict array1, const double *restrict array2);

double arrayMin(size_t len, const double *restrict array);
double arrayMax(size_t len, const double *restrict array);
double arrayRange(size_t len, const double *restrict array);
double arrayMinMV(size_t len, const double *restrict array, double missval);
double arrayMaxMV(size_t len, const double *restrict array, double missval);
double arrayRangeMV(size_t len, const double *restrict array, double missval);

double arraySum(size_t len, const double *restrict array);
double arraySumMV(size_t len, const double *restrict array, double missval);

double arrayMean(size_t len, const double *restrict array);
double arrayMeanMV(size_t len, const double *restrict array, double missval);
double arrayWeightedMean(size_t len, const double *restrict array, const double *restrict w, double missval);
double arrayWeightedMeanMV(size_t len, const double *restrict array, const double *restrict w, double missval);

double arrayAvgMV(size_t len, const double *restrict array, double missval);
double arrayWeightedAvgMV(size_t len, const double *restrict array, const double *restrict w, double missval);

#endif //  ARRAY_H

