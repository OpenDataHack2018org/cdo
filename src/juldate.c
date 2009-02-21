#include "cdo_int.h"


void encode_caldaysec(int calendar, int year, int month, int day, int hour, int minute, int *julday, int *secofday);
void decode_caldaysec(int calendar, int julday, int secofday, int *year, int *month, int *day, int *hour, int *minute);


juldate_t juldate_encode(int calendar, int date, int time)
{
  int year, month, day, hour, minute;
  juldate_t juldate;

  decode_date(date, &year, &month, &day);
  decode_time(time, &hour, &minute);

  encode_caldaysec(calendar, year, month, day, hour, minute, 
		   &juldate.julday, &juldate.secofday);

  return (juldate);
}


void juldate_decode(int calendar, juldate_t juldate, int *date, int *time)
{
  int year, month, day, hour, minute;
  
  decode_caldaysec(calendar, juldate.julday, juldate.secofday, 
		   &year, &month, &day, &hour, &minute);

  *date = encode_date(year, month, day);
  *time = encode_time(hour, minute);
}


juldate_t juldate_sub(juldate_t juldate2, juldate_t juldate1)
{
  juldate_t juldate;

  (void) julday_sub(juldate1.julday, juldate1.secofday, juldate2.julday, juldate2.secofday, 
		    &juldate.julday, &juldate.secofday);

  return (juldate);
}


juldate_t juldate_add_seconds(int seconds, juldate_t juldate)
{
  juldate_t juldate_new;

  juldate_new = juldate;

  julday_add_seconds(seconds, &juldate_new.julday, &juldate_new.secofday);

  return (juldate_new);
}


double juldate_to_seconds(juldate_t juldate)
{
  double seconds;

  seconds = juldate.julday*86400. + juldate.secofday;

  return (seconds);
}
