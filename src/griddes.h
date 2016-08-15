#ifndef _GRIDDES_H
#define _GRIDDES_H

#include <stdbool.h>

typedef struct {
  int    *mask;
  double *xvals;
  double *yvals;
  double *xbounds;
  double *ybounds;
  double *area;
  double  xfirst, yfirst;
  double  xlast, ylast;
  double  xinc, yinc;
  double  xpole, ypole, angle;    /* rotated north pole             */
  double  originLon;              /* lambert                        */
  double  originLat;
  double  lonParY;
  double  lat1;
  double  lat2;
  int     projflag;
  int     scanflag;
  bool    def_originLon;
  bool    def_originLat;
  bool    def_lonParY;
  bool    def_lat1;
  bool    def_lat2;
  double  a;
  double  lon_0;
  double  lat_0;
  double  lat_1;
  double  lat_2;
  bool    def_lon_0;
  bool    def_lat_0;
  bool    def_lat_1;
  bool    def_lat_2;
  int     prec;
  int     isRotated;              /* TRUE for rotated grids         */
  int     type;
  int     ntr;
  int    *rowlon;
  bool    genBounds;
  int     nvertex;
  long    size;
  int     xsize;
  int     ysize;
  int     np;
  int     lcomplex;
  bool    def_xfirst;
  bool    def_yfirst;
  bool    def_xlast;
  bool    def_ylast;
  bool    def_xinc;
  bool    def_yinc;
  int     nd, ni, ni2, ni3;
  int     number, position;
  unsigned char uuid[CDI_UUID_SIZE];
  char    path[16384];
  char    xname[CDI_MAX_NAME];
  char    xlongname[CDI_MAX_NAME];
  char    xunits[CDI_MAX_NAME];
  char    xdimname[CDI_MAX_NAME];
  char    yname[CDI_MAX_NAME];
  char    ylongname[CDI_MAX_NAME];
  char    yunits[CDI_MAX_NAME];
  char    ydimname[CDI_MAX_NAME];
  char    vdimname[CDI_MAX_NAME];
}
griddes_t;

void gridInit(griddes_t *grid);
int gridDefine(griddes_t grid);

int gridFromNCfile(const char *gridfile);
int gridFromH5file(const char *gridfile);

#endif  /* _GRIDDES_H */
