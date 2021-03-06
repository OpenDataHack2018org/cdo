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
#include <cdi.h>
#include "cdo_int.h"
#include "datetime.h"
#include "calendar.h"

juldate_t
juldate_encode(int calendar, int64_t date, int time)
{
  int year, month, day, hour, minute, second;
  juldate_t juldate;

  cdiDecodeDate(date, &year, &month, &day);
  cdiDecodeTime(time, &hour, &minute, &second);

  encode_caldaysec(calendar, year, month, day, hour, minute, second, &juldate.julday, &juldate.secofday);

  return juldate;
}

void
juldate_decode(int calendar, juldate_t juldate, int64_t *date, int *time)
{
  int year, month, day, hour, minute, second;

  decode_caldaysec(calendar, juldate.julday, juldate.secofday, &year, &month, &day, &hour, &minute, &second);

  *date = cdiEncodeDate(year, month, day);
  *time = cdiEncodeTime(hour, minute, second);
}

juldate_t
juldate_sub(juldate_t juldate2, juldate_t juldate1)
{
  juldate_t juldate;

  (void) julday_sub(juldate1.julday, juldate1.secofday, juldate2.julday, juldate2.secofday, &juldate.julday, &juldate.secofday);

  return juldate;
}

juldate_t
juldate_add_seconds(int64_t seconds, juldate_t juldate)
{
  juldate_t juldate_new = juldate;

  julday_add_seconds(seconds, &juldate_new.julday, &juldate_new.secofday);

  return juldate_new;
}

double
juldate_to_seconds(juldate_t juldate)
{
  return juldate.julday * 86400. + juldate.secofday;
  ;
}
