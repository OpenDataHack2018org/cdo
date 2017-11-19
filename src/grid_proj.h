#ifndef  GRID_PROJ_H
#define  GRID_PROJ_H


int cdo_lonlat_to_lcc(int gridID, size_t nvals, double *xvals, double *yvals);
int cdo_lcc_to_lonlat(int gridID, size_t nvals, double *xvals, double *yvals);

void cdo_sinu_to_lonlat(size_t nvals, double *xvals, double *yvals);
void cdo_laea_to_lonlat(int gridID, size_t nvals, double *xvals, double *yvals);

void cdo_proj_to_lonlat(char *proj4param, size_t nvals, double *xvals, double *yvals);

int proj_lonlat_to_lcc(double missval, double lon_0, double lat_0, double lat_1, double lat_2,
                       double a, double rf, size_t nvals, double *xvals, double *yvals);
int proj_lcc_to_lonlat(double missval, double lon_0, double lat_0, double lat_1, double lat_2,
                       double a, double rf, double x_0, double y_0, size_t nvals, double *xvals, double *yvals);

#endif  /* GRID_PROJ_H */
