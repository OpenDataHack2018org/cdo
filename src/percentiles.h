#ifndef _PERCENTILES_H
#define _PERCENTILES_H

double percentile(double *array, int len, double pn);
void percentile_set_method(const char *methodstr);
void percentile_check_number(double pn);

#endif /* _PERCENTILES_H */
