#ifndef _ARRAY_H
#define _ARRAY_H

const char *fpe_errstr(int fpeRaised);

int array_minmaxmean_val(const double *array, size_t len, double *rmin, double *rmax, double *rmean);

#endif // _ARRAY_H

