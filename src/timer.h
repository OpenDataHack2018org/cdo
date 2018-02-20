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
#ifndef TIMER_H
#define TIMER_H

extern int timer_read, timer_write;  // refactor: both pstream.cc and CDIread.cc
                                     // CDIwrite.cc defined in cdo.cc

void cdoProcessTime(double *utime, double *stime);
int timer_new(const char *text);
double timer_val(int it);
void timer_report(void);
void timer_start(int it);
void timer_stop(int it);

#endif
