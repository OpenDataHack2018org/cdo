#ifndef TIMEBASE_H
#define TIMEBASE_H

#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif

/* date format:  YYYYMMDD */
/* time format:  hhmmss   */

void decode_julday(int calendar, int64_t julday, int *year, int *mon, int *day);
int64_t encode_julday(int calendar, int year, int month, int day);

int64_t date_to_julday(int calendar, int64_t date);
int64_t julday_to_date(int calendar, int64_t julday);

int time_to_sec(int time);
int sec_to_time(int secofday);

void   julday_add_seconds(int64_t seconds, int64_t *julday, int *secofday);
void   julday_add(int days, int secs, int64_t *julday, int *secofday);
double julday_sub(int64_t julday1, int secofday1, int64_t julday2, int secofday2, int64_t *days, int *secs);

void encode_juldaysec(int calendar, int year, int month, int day, int hour, int minute, int second, int64_t *julday, int *secofday);
void decode_juldaysec(int calendar, int julday64_t, int secofday, int *year, int *month, int *day, int *hour, int *minute, int *second);

#ifdef __cplusplus
}
#endif

#endif /* TIMEBASE_H */

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
