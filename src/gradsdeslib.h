#ifndef  _GRADSDESLIB_H
#define  _GRADSDESLIB_H

#define  gaint     int
#define  gadouble  double
#define  galloc(x,y)  malloc(x)
#define  gree(x,y)    free(x)

/* Handling of missing data values. After the data I/O is done, 
   grid values are tested to see if they are within a small range 
   (+-value/EPSILON) of the missing value. If true, then the undef 
   mask is set to 0. If false, then the grid data values are good, 
   and the undef mask is set to 1. Everywhere else in the code, 
   undef tests are done on the mask values, not the data. */

#define EPSILON 1e5

#define  MAX_DSETS  1024

typedef struct {
  char *name;
  char *description;
  char *units;
  char *title;
  char *time;
  int dtype;
  int nx;
  int ny;
  int nz;
  int nt;
  int gridsize;
  int lscale;
  int loffset;
  int lmissval;
  double scale;
  double offset;
  double missval;
  double *array;
}
dset_obj_t;


#define  MAX_RECLEN   512
#define  MAX_NAMELEN  512
typedef struct {
  char name[MAX_NAMELEN];
  char dnam[MAX_NAMELEN];
  char title[MAX_NAMELEN];
  int bswap;
  long fhdr;
  long xyhdr;
  int seqflg;
  int yrflg;
  int zrflg;
  int tmplat;
  int pa2mb;
  int calendar;
  /* init !? */
  double undef;
  double ulow;
  double uhi;
  int dnum[5];               /* Dimension sizes for this file.        */
  double (*gr2ab[5]) (double *, double);
                               /* Addresses of routines to do conversion
                                  from grid coordinates to absolute
                                  coordinates for X, Y, Z.  All Date/time
                                  conversions handled by gr2t.          */
  double (*ab2gr[5]) (double *, double);
                               /* Addresses of routines to do conversion
                                  from absolute coordinates to grid
                                  coordinates for X,Y,Z.  All date/time
                                  conversions handled by t2gr.          */
  double *grvals[5];         /* Pointers to conversion information for
                                  grid-to-absolute conversion routines. */
  double *abvals[5];         /* Pointers to conversion information for
                                  absolute-to-grid conversion routines. */

  int nsets;
  int mergelevel;
  int lgeoloc;
  int lregion;
  int lprojtype;
  int lmetadata;
  dset_obj_t obj[MAX_DSETS];
}
dsets_t;

void dsets_init(dsets_t *dsets);

int read_gradsdes(char *filename, dsets_t *pfi);

gadouble liconv (gadouble *, gadouble);
gadouble gr2lev (gadouble *, gadouble);
gadouble lev2gr (gadouble *, gadouble);

#endif  /* _GRADSDESLIB_H */
