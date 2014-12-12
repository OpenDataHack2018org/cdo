#include <stdio.h>

#define  TIMESTAT_FIRST  1
#define  TIMESTAT_LAST   2
#define  TIMESTAT_MEAN   3


typedef struct {
  int   date;
  int   time;
} datetime_t;

typedef struct
{
  datetime_t v;
  datetime_t b[2];
} dtinfo_t;


typedef struct {
  int   date;
  int   time;
} datetime_type;

typedef struct
{
  datetime_type v;
  datetime_type b[2];
} dtinfo_type;



typedef struct
{
  size_t       nalloc;
  size_t       size;
  int          timestat_date;
  dtinfo_type  timestat;
  dtinfo_type *dtinfo;
} dtlist_type;



void    get_timestat_date(int *tstat_date);
void    datetime_avg(int dpy, int ndates, datetime_t *datetime);
void    datetime_avg_dtinfo(int dpy, int ndates, dtinfo_t *dtinfo);
void    taxisInqDTinfo(int taxisID, dtinfo_t *dtinfo);
void    taxisDefDTinfo(int taxisID, dtinfo_t dtinfo);

dtlist_type *dtlist_new(void);
void dtlist_delete(dtlist_type *dtlist);
void dtlist_taxisInqTimestep(int taxisID, int tsID, dtlist_type *dtlist);
void dtlist_taxisDefTimestep(int taxisID, int tsID, const dtlist_type *dtlist);
