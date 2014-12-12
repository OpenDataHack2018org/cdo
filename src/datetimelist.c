#include "cdi.h"
#include "dmemory.h"
#include "util.h"
#include "datetimelist.h"


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
  printf("init\n");
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

  dtlist->size += 1;

  dtlist->dtinfo[tsID].v.date = taxisInqVdate(taxisID);
  dtlist->dtinfo[tsID].v.time = taxisInqVtime(taxisID);
  if ( taxisHasBounds(taxisID) )
    {
      taxisInqVdateBounds(taxisID, &(dtlist->dtinfo[tsID].b[0].date), &(dtlist->dtinfo[tsID].b[1].date));
      taxisInqVtimeBounds(taxisID, &(dtlist->dtinfo[tsID].b[0].time), &(dtlist->dtinfo[tsID].b[1].time));
    }
}


void dtlist_taxisDefTimestep(int taxisID, int tsID, const dtlist_type *dtlist)
{
  if ( (size_t)tsID >= dtlist->size )
    cdoAbort("Internal error; tsID out of bounds!");

  taxisDefVdate(taxisID, dtlist->dtinfo[tsID].v.date);
  taxisDefVtime(taxisID, dtlist->dtinfo[tsID].v.time);
  if ( taxisHasBounds(taxisID) )
    {
      taxisDefVdateBounds(taxisID, dtlist->dtinfo[tsID].b[0].date, dtlist->dtinfo[tsID].b[1].date);
      taxisDefVtimeBounds(taxisID, dtlist->dtinfo[tsID].b[0].time, dtlist->dtinfo[tsID].b[1].time);
    }
}

