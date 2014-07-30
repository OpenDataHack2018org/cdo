#include "grid_search.h"


typedef struct {
   grid_search_vtable_t * vtable;
  // grid_t * bucket_grid;
   double * bucket_grid_x_coords;
   double * bucket_grid_y_coords;
  // grid_t * grid_data;
   unsigned num_buckets[2];
  // struct dep_list bucket_to_cell;
} bucket_search_t;


static
void bucket_search_delete(grid_search_t * search)
{
  bucket_search_t * bucket_search = (bucket_search_t *) search;

  //free(bucket_search->bucket_grid_x_coords);
  //free(bucket_search->bucket_grid_y_coords);
  //delete_grid(bucket_search->bucket_grid);
  //free_dep_list(&(bucket_search->bucket_to_cell));
  free(search);
}


static grid_search_vtable_t bucket_search_vtable = {
  .grid_search_delete                    = bucket_search_delete
};

grid_search_t * bucket_search_new(grid_t grid_data)
{
   bucket_search_t * search;

   search = malloc(sizeof(bucket_search_t));

   search->vtable = &bucket_search_vtable;

   return ((grid_search_t *) search);
}
