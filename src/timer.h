#ifndef TIMER_H
#define TIMER_H

void cdoProcessTime(double *utime, double *stime);
int timer_new(const char *text);
double timer_val(int it);
void timer_report(void);
void timer_start(int it);
void timer_stop(int it);

#endif
