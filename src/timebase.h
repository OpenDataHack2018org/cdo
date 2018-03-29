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
#ifndef _TIMEBASE_H
#define _TIMEBASE_H

#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif

/* date format:  YYYYMMDD */
/* time format:  hhmmss   */

void decode_julday(int calendar, int julday, int *year, int *mon, int *day);
int encode_julday(int calendar, int year, int month, int day);

int date_to_julday(int calendar, int date);
int julday_to_date(int calendar, int julday);

int time_to_sec(int time);
int sec_to_time(int secofday);

void julday_add_seconds(int64_t seconds, int *julday, int *secofday);
void julday_add(int days, int secs, int *julday, int *secofday);
double julday_sub(int julday1, int secofday1, int julday2, int secofday2, int *days, int *secs);

void encode_juldaysec(int calendar, int year, int month, int day, int hour, int minute, int second, int *julday, int *secofday);
void decode_juldaysec(int calendar, int julday, int secofday, int *year, int *month, int *day, int *hour, int *minute,
                      int *second);

#if defined(__cplusplus)
}
#endif

#endif /* _TIMEBASE_H */

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
