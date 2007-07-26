#ifndef  _TIMEBASE_H
#define  _TIMEBASE_H

/* date format:  YYYYMMDD */
/* time format:  hhmm     */

void decode_date(int date, int *year, int *month, int *day);
int encode_date(int year, int month, int day);

void decode_time(int time, int *hour, int *minute);
int encode_time(int hour, int minute);

void decode_julday(int julday, int *year, int *mon, int	*day);
int encode_julday(int year, int month, int day);

int date_to_julday(int date);
int julday_to_date(int julday);

int time_to_sec(int time);
int sec_to_time(int secofday);

void julday_add_seconds(int seconds, int *julday, int *secofday);
void julday_add(int days, int secs, int *julday, int *secofday);
int julday_sub(int julday1, int secofday1, int julday2, int secofday2, int *days, int *secs);

void encode_juldaysec(int year, int month, int day, int hour, int minute, int *julday, int *secofday);
void decode_juldaysec(int julday, int secofday, int *year, int *month, int *day, int *hour, int *minute);

#endif  /* _TIMEBASE_H */
