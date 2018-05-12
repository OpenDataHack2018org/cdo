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
#include "util_string.h"

TimeStat CDO_Timestat_Date = TimeStat::UNDEF;
bool CDO_Timestat_Bounds = false;

void
setTimestatDate(const char *optarg)
{
  TimeStat timestatdate = TimeStat::UNDEF;

  // clang-format off
  if      (STR_IS_EQ(optarg, "first"))   timestatdate = TimeStat::FIRST;
  else if (STR_IS_EQ(optarg, "last"))    timestatdate = TimeStat::LAST;
  else if (STR_IS_EQ(optarg, "middle"))  timestatdate = TimeStat::MEAN;
  else if (STR_IS_EQ(optarg, "midhigh")) timestatdate = TimeStat::MIDHIGH;
  // clang-format on

  if (timestatdate == TimeStat::UNDEF) cdoAbort("option --%s: unsupported argument: %s", "timestat_date", optarg);

  CDO_Timestat_Date = timestatdate;
}

static void
getTimestatDate(TimeStat *tstat_date)
{
  char *envstr = getenv("CDO_TIMESTAT_DATE");
  if (envstr == NULL) envstr = getenv("RUNSTAT_DATE");
  if (envstr)
    {
      TimeStat env_date = TimeStat::UNDEF;
      char envstrl[8];

      memcpy(envstrl, envstr, 8);
      envstrl[7] = 0;
      strtolower(envstrl);

      // clang-format off
      if      (memcmp(envstrl, "first", 5) == 0)   env_date = TimeStat::FIRST;
      else if (memcmp(envstrl, "last", 4) == 0)    env_date = TimeStat::LAST;
      else if (memcmp(envstrl, "middle", 6) == 0)  env_date = TimeStat::MEAN;
      else if (memcmp(envstrl, "midhigh", 7) == 0) env_date = TimeStat::MIDHIGH;
      // clang-format on

      if (env_date != TimeStat::UNDEF)
        {
          *tstat_date = env_date;
          if (cdoVerbose) cdoPrint("Set CDO_TIMESTAT_DATE to %s", envstr);
        }
    }
}

void
dtlist_init(dtlist_type *dtlist)
{
  dtlist->nalloc = 0;
  dtlist->size = 0;
  dtlist->calendar = -1;
  dtlist->has_bounds = -1;
  dtlist->stat = TimeStat::LAST;
  dtlist->dtinfo = NULL;

  static bool init = false;
  if (!init) getTimestatDate(&CDO_Timestat_Date);
  init = true;
}

dtlist_type *
dtlist_new(void)
{
  dtlist_type *dtlist = (dtlist_type *) Malloc(sizeof(dtlist_type));
  dtlist_init(dtlist);
  return dtlist;
}

void
dtlist_delete(dtlist_type *dtlist)
{
  if (dtlist->nalloc > 0 && dtlist->dtinfo) Free(dtlist->dtinfo);
  Free(dtlist);
}

void
dtlist_taxisInqTimestep(dtlist_type *dtlist, int taxisID, int tsID)
{
  size_t NALLOC = 128;

  if ((size_t) tsID >= dtlist->nalloc)
    {
      dtlist->nalloc += NALLOC;
      dtlist->dtinfo = (dtinfo_type *) Realloc(dtlist->dtinfo, dtlist->nalloc * sizeof(dtinfo_type));
    }

  if ((size_t) tsID >= dtlist->size) dtlist->size = (size_t) tsID + 1;

  dtlist->dtinfo[tsID].v.date = taxisInqVdate(taxisID);
  dtlist->dtinfo[tsID].v.time = taxisInqVtime(taxisID);

  dtlist->dtinfo[tsID].c.date = dtlist->dtinfo[tsID].v.date;
  dtlist->dtinfo[tsID].c.time = dtlist->dtinfo[tsID].v.time;

  if (tsID == 0)
    {
      if (dtlist->has_bounds == -1)
        {
          dtlist->has_bounds = 0;
          if (taxisHasBounds(taxisID)) dtlist->has_bounds = 1;
        }

      if (dtlist->calendar == -1)
        {
          dtlist->calendar = taxisInqCalendar(taxisID);
        }
    }

  if (dtlist->has_bounds)
    {
      taxisInqVdateBounds(taxisID, &(dtlist->dtinfo[tsID].b[0].date), &(dtlist->dtinfo[tsID].b[1].date));
      taxisInqVtimeBounds(taxisID, &(dtlist->dtinfo[tsID].b[0].time), &(dtlist->dtinfo[tsID].b[1].time));

      if (CDO_Timestat_Bounds && dtlist->dtinfo[tsID].v.date == dtlist->dtinfo[tsID].b[1].date
          && dtlist->dtinfo[tsID].v.time == dtlist->dtinfo[tsID].b[1].time)
        {
          int calendar = dtlist->calendar;

          int64_t vdate = dtlist->dtinfo[tsID].b[0].date;
          int vtime = dtlist->dtinfo[tsID].b[0].time;
          juldate_t juldate1 = juldate_encode(calendar, vdate, vtime);

          vdate = dtlist->dtinfo[tsID].b[1].date;
          vtime = dtlist->dtinfo[tsID].b[1].time;
          juldate_t juldate2 = juldate_encode(calendar, vdate, vtime);

          // int hour, minute, second;
          // cdiDecodeTime(vtime, &hour, &minute, &second);

          if (vtime == 0 && juldate_to_seconds(juldate1) < juldate_to_seconds(juldate2))
            {
              juldate_t juldate = juldate_add_seconds(-1, juldate2);
              juldate_decode(calendar, juldate, &vdate, &vtime);

              dtlist->dtinfo[tsID].c.date = vdate;
              dtlist->dtinfo[tsID].c.time = vtime;
            }
        }
    }
  else
    {
      dtlist->dtinfo[tsID].b[0].date = 0;
      dtlist->dtinfo[tsID].b[1].date = 0;
      dtlist->dtinfo[tsID].b[0].time = 0;
      dtlist->dtinfo[tsID].b[1].time = 0;
    }
}

void
dtlist_taxisDefTimestep(dtlist_type *dtlist, int taxisID, int tsID)
{
  if (tsID < 0 || (size_t) tsID >= dtlist->size) cdoAbort("Internal error; tsID out of bounds!");

  taxisDefVdate(taxisID, dtlist->dtinfo[tsID].v.date);
  taxisDefVtime(taxisID, dtlist->dtinfo[tsID].v.time);
  if (dtlist->has_bounds)
    {
      taxisDefVdateBounds(taxisID, dtlist->dtinfo[tsID].b[0].date, dtlist->dtinfo[tsID].b[1].date);
      taxisDefVtimeBounds(taxisID, dtlist->dtinfo[tsID].b[0].time, dtlist->dtinfo[tsID].b[1].time);
    }
}

void
dtlist_mean(dtlist_type *dtlist, int nsteps)
{
  int64_t vdate;
  int vtime;

  if (nsteps % 2 == 0)
    {
      int calendar = dtlist->calendar;

//#define TEST_DTLIST_MEAN 1
#ifdef TEST_DTLIST_MEAN
      vdate = dtlist->dtinfo[0].v.date;
      vtime = dtlist->dtinfo[0].v.time;
      juldate_t juldate0 = juldate_encode(calendar, vdate, vtime);

      juldate_t juldate;
      double seconds = 0;
      for (int i = 1; i < nsteps; ++i)
        {
          vdate = dtlist->dtinfo[i].v.date;
          vtime = dtlist->dtinfo[i].v.time;
          juldate = juldate_encode(calendar, vdate, vtime);

          seconds += juldate_to_seconds(juldate_sub(juldate, juldate0));
        }

      juldate = juldate_add_seconds(lround(seconds / nsteps), juldate0);
      juldate_decode(calendar, juldate, &vdate, &vtime);
#else
      vdate = dtlist->dtinfo[nsteps / 2 - 1].v.date;
      vtime = dtlist->dtinfo[nsteps / 2 - 1].v.time;
      juldate_t juldate1 = juldate_encode(calendar, vdate, vtime);

      vdate = dtlist->dtinfo[nsteps / 2].v.date;
      vtime = dtlist->dtinfo[nsteps / 2].v.time;
      juldate_t juldate2 = juldate_encode(calendar, vdate, vtime);

      double seconds = juldate_to_seconds(juldate_sub(juldate2, juldate1)) / 2;
      juldate_t juldatem = juldate_add_seconds(lround(seconds), juldate1);
      juldate_decode(calendar, juldatem, &vdate, &vtime);
#endif
    }
  else
    {
      vdate = dtlist->dtinfo[nsteps / 2].v.date;
      vtime = dtlist->dtinfo[nsteps / 2].v.time;
    }

  dtlist->timestat.v.date = vdate;
  dtlist->timestat.v.time = vtime;
}

void
dtlist_midhigh(dtlist_type *dtlist, int nsteps)
{
  dtlist->timestat.v.date = dtlist->dtinfo[nsteps / 2].v.date;
  dtlist->timestat.v.time = dtlist->dtinfo[nsteps / 2].v.time;
}

void
dtlist_stat_taxisDefTimestep(dtlist_type *dtlist, int taxisID, int nsteps)
{
  if ((size_t) nsteps > dtlist->size) cdoAbort("Internal error; unexpected nsteps=%d (limit=%ld)!", nsteps, dtlist->size);

  TimeStat stat = dtlist->stat;
  if (CDO_Timestat_Date != TimeStat::UNDEF) stat = CDO_Timestat_Date;

  // clang-format off
  if      (stat == TimeStat::MEAN)    dtlist_mean(dtlist, nsteps);
  else if (stat == TimeStat::MIDHIGH) dtlist_midhigh(dtlist, nsteps);
  else if (stat == TimeStat::FIRST)   dtlist->timestat.v = dtlist->dtinfo[0].v;
  else if (stat == TimeStat::LAST)    dtlist->timestat.v = dtlist->dtinfo[nsteps - 1].v;
  else cdoAbort("Internal error; implementation missing for timestat=%d", (int)stat);
  // clang-format on

  if (dtlist->has_bounds)
    {
      dtlist->timestat.b[0] = dtlist->dtinfo[0].b[0];
      dtlist->timestat.b[1] = dtlist->dtinfo[nsteps - 1].b[1];
    }
  else
    {
      dtlist->timestat.b[0] = dtlist->dtinfo[0].v;
      dtlist->timestat.b[1] = dtlist->dtinfo[nsteps - 1].v;
    }

  taxisDefVdate(taxisID, dtlist->timestat.v.date);
  taxisDefVtime(taxisID, dtlist->timestat.v.time);
  // if ( dtlist->has_bounds )
  {
    taxisDefVdateBounds(taxisID, dtlist->timestat.b[0].date, dtlist->timestat.b[1].date);
    taxisDefVtimeBounds(taxisID, dtlist->timestat.b[0].time, dtlist->timestat.b[1].time);
  }
}

void
dtlist_shift(dtlist_type *dtlist)
{
  for (size_t inp = 0; inp < dtlist->size - 1; inp++)
    {
      dtlist->dtinfo[inp] = dtlist->dtinfo[inp + 1];
    }
}

void
dtlist_set_stat(dtlist_type *dtlist, TimeStat stat)
{
  dtlist->stat = stat;
}

void
dtlist_set_calendar(dtlist_type *dtlist, int calendar)
{
  dtlist->calendar = calendar;
}

int64_t
dtlist_get_vdate(dtlist_type *dtlist, int tsID)
{
  if (tsID < 0 || (size_t) tsID >= dtlist->size) cdoAbort("Internal error; tsID out of bounds!");

  return dtlist->dtinfo[tsID].c.date;
}

int
dtlist_get_vtime(dtlist_type *dtlist, int tsID)
{
  if (tsID < 0 || (size_t) tsID >= dtlist->size) cdoAbort("Internal error; tsID out of bounds!");

  return dtlist->dtinfo[tsID].c.time;
}

void
datetime_avg(int calendar, int ndates, CdoDateTime *datetime)
{
  int64_t vdate;
  int vtime;

  if (ndates % 2 == 0)
    {
      vdate = datetime[ndates / 2 - 1].date;
      vtime = datetime[ndates / 2 - 1].time;
      juldate_t juldate1 = juldate_encode(calendar, vdate, vtime);

      vdate = datetime[ndates / 2].date;
      vtime = datetime[ndates / 2].time;
      juldate_t juldate2 = juldate_encode(calendar, vdate, vtime);

      double seconds = juldate_to_seconds(juldate_sub(juldate2, juldate1)) / 2;
      juldate_t juldatem = juldate_add_seconds(lround(seconds), juldate1);
      juldate_decode(calendar, juldatem, &vdate, &vtime);
    }
  else
    {
      vdate = datetime[ndates / 2].date;
      vtime = datetime[ndates / 2].time;
    }

  datetime[ndates].date = vdate;
  datetime[ndates].time = vtime;
}
