#include <string.h>
#include <math.h>

#include <cdi.h>
#include "cdo.h"
#include "datetime.h"


void taxisInqDTinfo(int taxisID, dtinfo_t *dtinfo)
{
  dtinfo->v.date = taxisInqVdate(taxisID);
  dtinfo->v.time = taxisInqVtime(taxisID);
  if ( taxisHasBounds(taxisID) )
    {
      taxisInqVdateBounds(taxisID, &(dtinfo->b[0].date), &(dtinfo->b[1].date));
      taxisInqVtimeBounds(taxisID, &(dtinfo->b[0].time), &(dtinfo->b[1].time));
    }
}


void taxisDefDTinfo(int taxisID, dtinfo_t dtinfo)
{
  taxisDefVdate(taxisID, dtinfo.v.date);
  taxisDefVtime(taxisID, dtinfo.v.time);
  if ( taxisHasBounds(taxisID) )
    {
      taxisDefVdateBounds(taxisID, dtinfo.b[0].date, dtinfo.b[1].date);
      taxisDefVtimeBounds(taxisID, dtinfo.b[0].time, dtinfo.b[1].time);
    }
}


void dtlist_init(dtlist_type *dtlist)
{
  dtlist->nalloc = 0;
  dtlist->size   = 0;
  dtlist->timestat_date  = TIMESTAT_LAST;
  dtlist->dtinfo = NULL;
}


dtlist_type *dtlist_new(void)
{
  dtlist_type *dtlist = (dtlist_type *) malloc(sizeof(dtlist_type));

  dtlist_init(dtlist);

  return dtlist;
}


void dtlist_delete(dtlist_type *dtlist)
{
  if ( dtlist->nalloc > 0 && dtlist->dtinfo )
    free(dtlist->dtinfo);

  free(dtlist);
}


void dtlist_taxisInqTimestep(int taxisID, int tsID, dtlist_type *dtlist)
{
  size_t NALLOC = 512;

  if ( (size_t)tsID >= dtlist->nalloc )
    {
      dtlist->nalloc += NALLOC;
      dtlist->dtinfo = (dtinfo_type *) realloc(dtlist->dtinfo, dtlist->nalloc*sizeof(dtinfo_type));
    }

  if ( (size_t)tsID >= dtlist->size ) dtlist->size = (size_t)tsID + 1;

  dtlist->dtinfo[tsID].v.date = taxisInqVdate(taxisID);
  dtlist->dtinfo[tsID].v.time = taxisInqVtime(taxisID);

  if ( taxisHasBounds(taxisID) )
    {
      taxisInqVdateBounds(taxisID, &(dtlist->dtinfo[tsID].b[0].date), &(dtlist->dtinfo[tsID].b[1].date));
      taxisInqVtimeBounds(taxisID, &(dtlist->dtinfo[tsID].b[0].time), &(dtlist->dtinfo[tsID].b[1].time));
    }
  else
    {
      dtlist->dtinfo[tsID].b[0].date = 0;
      dtlist->dtinfo[tsID].b[1].date = 0;
      dtlist->dtinfo[tsID].b[0].time = 0;
      dtlist->dtinfo[tsID].b[1].time = 0;
    }
}


void dtlist_stat_taxisDefTimestep(int taxisID, int nsteps, const dtlist_type *dtlist)
{
  if ( (size_t)nsteps != dtlist->size )
    cdoAbort("Internal error; unexpected nsteps!");

  taxisDefVdate(taxisID, dtlist->timestat.v.date);
  taxisDefVtime(taxisID, dtlist->timestat.v.time);
  if ( taxisHasBounds(taxisID) )
    {
      taxisDefVdateBounds(taxisID, dtlist->timestat.b[0].date, dtlist->timestat.b[1].date);
      taxisDefVtimeBounds(taxisID, dtlist->timestat.b[0].time, dtlist->timestat.b[1].time);
    }
}


void dtlist_taxisDefTimestep(int taxisID, int tsID, const dtlist_type *dtlist)
{
  if ( tsID < 0 || (size_t)tsID >= dtlist->size )
    cdoAbort("Internal error; tsID out of bounds!");

  taxisDefVdate(taxisID, dtlist->dtinfo[tsID].v.date);
  taxisDefVtime(taxisID, dtlist->dtinfo[tsID].v.time);
  if ( taxisHasBounds(taxisID) )
    {
      taxisDefVdateBounds(taxisID, dtlist->dtinfo[tsID].b[0].date, dtlist->dtinfo[tsID].b[1].date);
      taxisDefVtimeBounds(taxisID, dtlist->dtinfo[tsID].b[0].time, dtlist->dtinfo[tsID].b[1].time);
    }
}


void datetime_avg_dtinfo(int calendar, int ndates, dtinfo_t *dtinfo)
{
  int vdate, vtime;
  juldate_t juldate1, juldate2, juldatem;
  double seconds;
  /*
  for ( i = 0; i < ndates; i++ )
    fprintf(stdout, "%4d %d %d\n", i+1, dtinfo[i].v.date, dtinfo[i].v.time);
  */
  if ( ndates%2 == 0 )
    {
      /*
      vdate = dtinfo[ndates-1].v.date;
      vtime = dtinfo[ndates-1].v.time;
      */
      vdate = dtinfo[ndates/2-1].v.date;
      vtime = dtinfo[ndates/2-1].v.time;
      juldate1 = juldate_encode(calendar, vdate, vtime);

      vdate = dtinfo[ndates/2].v.date;
      vtime = dtinfo[ndates/2].v.time;
      juldate2 = juldate_encode(calendar, vdate, vtime);

      seconds = juldate_to_seconds(juldate_sub(juldate2, juldate1)) / 2;
      juldatem = juldate_add_seconds((int)lround(seconds), juldate1);
      juldate_decode(calendar, juldatem, &vdate, &vtime);
    }
  else
    {
      vdate = dtinfo[ndates/2].v.date;
      vtime = dtinfo[ndates/2].v.time;
    }

  dtinfo[ndates].v.date = vdate;
  dtinfo[ndates].v.time = vtime;
  /*
  fprintf(stdout, "res: %d %d\n\n", dtinfo[ndates].v.date, dtinfo[ndates].v.time);
  */
}


void datetime_avg(int calendar, int ndates, datetime_t *datetime)
{
  int vdate, vtime;
  juldate_t juldate1, juldate2, juldatem;
  double seconds;
  /*
  for ( i = 0; i < ndates; i++ )
    fprintf(stdout, "%4d %d %d\n", i+1, datetime[i].date, datetime[i].time);
  */
  if ( ndates%2 == 0 )
    {
      /*
      vdate = datetime[ndates-1].date;
      vtime = datetime[ndates-1].time;
      */
      vdate = datetime[ndates/2-1].date;
      vtime = datetime[ndates/2-1].time;
      juldate1 = juldate_encode(calendar, vdate, vtime);

      vdate = datetime[ndates/2].date;
      vtime = datetime[ndates/2].time;
      juldate2 = juldate_encode(calendar, vdate, vtime);

      seconds = juldate_to_seconds(juldate_sub(juldate2, juldate1)) / 2;
      juldatem = juldate_add_seconds((int)lround(seconds), juldate1);
      juldate_decode(calendar, juldatem, &vdate, &vtime);
    }
  else
    {
      vdate = datetime[ndates/2].date;
      vtime = datetime[ndates/2].time;
    }

  datetime[ndates].date = vdate;
  datetime[ndates].time = vtime;
  /*
  fprintf(stdout, "res: %d %d\n\n", datetime[ndates].date, datetime[ndates].time);
  */
}


void get_timestat_date(int *tstat_date)
{
  char *envstr;

  envstr = getenv("TIMESTAT_DATE");
  if ( envstr == NULL ) envstr = getenv("RUNSTAT_DATE");
  if ( envstr )
    {
      int env_date = -1;
      char envstrl[8];

      memcpy(envstrl, envstr, 8);
      envstrl[7] = 0;
      strtolower(envstrl);

      if      ( memcmp(envstrl, "first", 5)  == 0 )  env_date = TIMESTAT_FIRST;
      else if ( memcmp(envstrl, "last", 4)   == 0 )  env_date = TIMESTAT_LAST;
      else if ( memcmp(envstrl, "middle", 6) == 0 )  env_date = TIMESTAT_MEAN;

      if ( env_date >= 0 )
	{
	  *tstat_date = env_date;

	  if ( cdoVerbose ) cdoPrint("Set TIMESTAT_DATE to %s", envstr);
	}
    }
}
