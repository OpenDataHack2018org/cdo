#ifndef _GRID_PROJ_H
#define _GRID_PROJ_H


int cdo_lonlat_to_lcc(int gridID, size_t nvals, double *xvals, double *yvals);
int cdo_lcc_to_lonlat(int gridID, size_t nvals, double *xvals, double *yvals);

void cdo_sinu_to_lonlat(size_t nvals, double *xvals, double *yvals);
void cdo_laea_to_lonlat(int gridID, size_t nvals, double *xvals, double *yvals);

void cdo_proj_to_lonlat(char *proj4param, size_t nvals, double *xvals, double *yvals);

#endif  /* _GRID_PROJ_H */
