#ifndef _GRID_DATA_H
#define _GRID_DATA_H

#ifndef  M_PI
#define  M_PI        3.14159265358979323846264338327950288  /* pi */
#endif

#ifndef  M_PI_2
#define  M_PI_2      1.57079632679489661923132169163975144  /* pi/2 */
#endif

#define EARTH_RADIUS 6371.2290

static double const EarthRadius  = EARTH_RADIUS;
static double const rad          = M_PI / 180.0;
static double const deg          = 180.0 / M_PI;
static double const EarthRadius2 = EARTH_RADIUS * EARTH_RADIUS / 2.0 ;


// forward declaration required by grid_vtable
struct grid;

typedef struct {
   void (*delete)(struct grid *);
} grid_vtable_t;


typedef struct grid {
  grid_vtable_t *vtable;
} grid_t;

void grid_delete(grid_t * grid);

grid_t * reg2d_grid_new(double * coordinates_x, double * coordinates_y, int nx, int ny, int iscyclic);

#endif  /* _GRID_DATA_H */ 
