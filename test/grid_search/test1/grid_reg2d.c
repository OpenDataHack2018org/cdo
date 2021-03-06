#include "grid_search.h"

typedef struct reg2d_grid {

  grid_vtable_t * vtable;
  /*
   struct dep_list cell_to_neigh; //!< dependency list containing cell neighbourhood
                                  //!< information (automatically generated once it is
                                  //!< needed)
				  */
  double * cell_corners_x,       //!< latitude data
         * cell_corners_y,       //!< longitude data
         * cos_cell_corners_x,
         * sin_cell_corners_x,
         * cos_cell_corners_y,
         * sin_cell_corners_y;

   unsigned cell_corners_x_by_user; //!< indicates whether cell_corners_x was provided by
                                    //!< the user or automatically generated by a grid
                                    //!< routine
   unsigned cell_corners_y_by_user; //!< indicates whether cell_corners_y was provided by
                                    //!< the user or automatically generated by a grid
                                    //!< routine

  int nx, ny;
  int iscyclic;

  struct grid_search * grid_search;

} reg2d_grid_t;


static unsigned const * get_cell_x_coord_indices_reg2d(struct grid * grid,
                                                       unsigned cell_index) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   static unsigned cell_corners[4];
   unsigned x_index, y_index;

   y_index = cell_index / reg2d_grid->num_cells[0];
   x_index = cell_index - y_index * reg2d_grid->num_cells[0];

   cell_corners[0] = x_index;
   
   if (reg2d_grid->cyclic[0] &&
       x_index+1 == reg2d_grid->num_cells[0]) {
       
      cell_corners[1] = 0;
      cell_corners[2] = 0;
   } else {
      cell_corners[1] = x_index + 1;
      cell_corners[2] = x_index + 1;
   }
   cell_corners[3] = x_index;

   return cell_corners;
}

static unsigned const * get_cell_y_coord_indices_reg2d(struct grid * grid,
                                                       unsigned cell_index) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   static unsigned cell_corners[4];
   unsigned y_index;

   y_index = cell_index / reg2d_grid->num_cells[0];

   cell_corners[0] = y_index;
   cell_corners[1] = y_index;
   cell_corners[2] = y_index + 1;
   cell_corners[3] = y_index + 1;

   return cell_corners;
}

static void get_grid_cell_reg2d(struct grid_t grid, unsigned cell_index,
                                struct grid_cell * cell) {

   struct reg2d_grid * reg2d_grid;

   reg2d_grid = (struct reg2d_grid *)grid;

   unsigned num_corners;

   //num_corners = get_num_cell_corners_reg2d(grid, cell_index);
   num_corners = 4;

   // if the memory for the coordinates needs to be reallocated
   if (num_corners > cell->array_size) {
      cell->coordinates_x = realloc (cell->coordinates_x, num_corners * sizeof(cell->coordinates_x[0]));
      cell->coordinates_y = realloc (cell->coordinates_y, num_corners * sizeof(cell->coordinates_y[0]));
      cell->coordinates_xyz = realloc (cell->coordinates_xyz, 3 * num_corners * sizeof(cell->coordinates_xyz[0]));
      cell->edge_type = realloc (cell->edge_type, num_corners * sizeof(cell->edge_type[0]));
      cell->array_size = num_corners;
   }

   cell->num_corners = num_corners;

   // set the corners for of the cell
   unsigned const * x_indices, * y_indices;
   unsigned i;

   x_indices = get_cell_x_coord_indices_reg2d(grid, cell_index);
   y_indices = get_cell_y_coord_indices_reg2d(grid, cell_index);

   for (i = 0; i < 4; ++i) {
      cell->coordinates_x[i] = reg2d_grid->cell_corners_x[x_indices[i]];
      cell->coordinates_y[i] = reg2d_grid->cell_corners_y[y_indices[i]];
      /*
      cell->coordinates_xyz[0+i*3] = reg2d_grid->cos_cell_corners_x[x_indices[i]] *
                                     reg2d_grid->cos_cell_corners_y[y_indices[i]];
      cell->coordinates_xyz[1+i*3] = reg2d_grid->sin_cell_corners_x[x_indices[i]] *
                                     reg2d_grid->cos_cell_corners_y[y_indices[i]];
      cell->coordinates_xyz[2+i*3] = reg2d_grid->sin_cell_corners_y[y_indices[i]];
      */
   }


   // set the edge type for the cell
   cell->edge_type[0] = LAT_CIRCLE;
   cell->edge_type[1] = LON_CIRCLE;
   cell->edge_type[2] = LAT_CIRCLE;
   cell->edge_type[3] = LON_CIRCLE;
}

static
void delete_grid_reg2d(grid_t * grid)
{
  reg2d_grid_t * reg2d_grid;

  reg2d_grid = (reg2d_grid_t *)grid;

  // free_dep_list(&(reg2d_grid->cell_to_neigh));
  if ( !reg2d_grid->cell_corners_x_by_user )
    free(reg2d_grid->cell_corners_x);
  if ( !reg2d_grid->cell_corners_y_by_user )
    free(reg2d_grid->cell_corners_y);
  /*
  free(reg2d_grid->sin_cell_corners_x);
  free(reg2d_grid->cos_cell_corners_x);
  free(reg2d_grid->sin_cell_corners_y);
  free(reg2d_grid->cos_cell_corners_y);
  */
  if ( reg2d_grid->grid_search != NULL )
    grid_search_delete(reg2d_grid->grid_search);

  free(reg2d_grid);
}


static grid_vtable_t reg2d_grid_vtable = {
  .get_grid_cell               = get_grid_cell_reg2d,
  .get_cell_x_coord_indices    = get_cell_x_coord_indices_reg2d,
  .get_cell_y_coord_indices    = get_cell_y_coord_indices_reg2d,
  .delete                      = delete_grid_reg2d
};

grid_t * reg2d_grid_new(double * coordinates_x, double * coordinates_y, int nx, int ny, int iscyclic)
{
  reg2d_grid_t * reg2d_grid;

  reg2d_grid = malloc(sizeof(reg2d_grid_t));

  reg2d_grid->vtable = &reg2d_grid_vtable;
   // init_dep_list(&reg2d_grid->cell_to_neigh);

  reg2d_grid->cell_corners_x = coordinates_x;
  reg2d_grid->cell_corners_y = coordinates_y;
  reg2d_grid->nx = nx;
  reg2d_grid->ny = ny;
  reg2d_grid->iscyclic = iscyclic;
  reg2d_grid->cell_corners_x_by_user = 1;
  reg2d_grid->cell_corners_y_by_user = 1;
  reg2d_grid->grid_search = NULL;
  /*
  if (coordinates_x != NULL) {
    unsigned num_x_coords = get_size_x_coords_reg2d((struct grid *) reg2d_grid);
    reg2d_grid->sin_cell_corners_x =
      malloc(num_x_coords * sizeof(*(reg2d_grid->sin_cell_corners_x)));
    reg2d_grid->cos_cell_corners_x = 
      malloc(num_x_coords * sizeof(*(reg2d_grid->cos_cell_corners_x)));

    for (unsigned i = 0; i < num_x_coords; ++i) {
      
      reg2d_grid->sin_cell_corners_x[i] = sin(reg2d_grid->cell_corners_x[i]);
      reg2d_grid->cos_cell_corners_x[i] = cos(reg2d_grid->cell_corners_x[i]);
    }
  } else {
    reg2d_grid->sin_cell_corners_x = NULL;
    reg2d_grid->cos_cell_corners_x = NULL;
  }

  if (coordinates_y != NULL) {
    unsigned num_y_coords = get_size_y_coords_reg2d((struct grid *) reg2d_grid);
    reg2d_grid->sin_cell_corners_y =
      malloc(num_y_coords * sizeof(*(reg2d_grid->sin_cell_corners_y)));
    reg2d_grid->cos_cell_corners_y = 
      malloc(num_y_coords * sizeof(*(reg2d_grid->cos_cell_corners_y)));
    
    for (unsigned i = 0; i < num_y_coords; ++i) {

      reg2d_grid->sin_cell_corners_y[i] = sin(reg2d_grid->cell_corners_y[i]);
      reg2d_grid->cos_cell_corners_y[i] = cos(reg2d_grid->cell_corners_y[i]);
    }
  } else {
    reg2d_grid->sin_cell_corners_y = NULL;
    reg2d_grid->cos_cell_corners_y = NULL;
  }
  */
  return ((grid_t *) reg2d_grid);
}

#if defined(GRID_SEARCH)
int main(void)
{
  double coordinates_x[] = {1, 2, 3, 4, 5};
  double coordinates_y[] = {1, 2, 3, 4, 5};
  int nx = 4, ny = 4;
  int iscyclic = 0;

  grid_t *rgrid = reg2d_grid_new(coordinates_x, coordinates_y, nx, ny, iscyclic);

  grid_delete(rgrid);

  return 0;
}
#endif
