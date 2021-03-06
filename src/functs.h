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
#ifndef _FUNCTS_H
#define _FUNCTS_H

#define func_fld 7
#define func_all 8
#define func_hrd 9

#define func_min 10
#define func_max 11
#define func_range 12
#define func_sum 13
#define func_avg 14
#define func_mean 15
#define func_std 16
#define func_std1 17
#define func_var 18
#define func_var1 19
#define func_pctl 20
#define func_cor 21
#define func_covar 22
#define func_avgw 23
#define func_meanw 24
#define func_stdw 25
#define func_std1w 26
#define func_varw 27
#define func_var1w 28
#define func_skew 29
#define func_kurt 30

#define func_brs 31
#define func_rank 32
#define func_roc 33
#define func_minidx 34
#define func_maxidx 35


#define func_add 41
#define func_sub 42
#define func_mul 43
#define func_div 44
#define func_mod 45

#define func_atan2 50
#define func_setmiss 51

#define func_read 60
#define func_write 61

#define func_month 84
#define func_year 85
#define func_time 86
#define func_date 87
#define func_step 88
#define func_datetime 89

#define func_lon 98
#define func_lat 99

enum cmp_flag
{
  CMP_NAME = 1,
  CMP_GRID = 2,
  CMP_NLEVEL = 4,
  CMP_GRIDSIZE = 8,
  CMP_HRD = CMP_NAME | CMP_GRIDSIZE,
  CMP_DIM = CMP_GRIDSIZE | CMP_NLEVEL | CMP_GRID,
  CMP_ALL = CMP_NAME | CMP_GRIDSIZE | CMP_NLEVEL | CMP_GRID
};

#endif /* _FUNCTS_H */
