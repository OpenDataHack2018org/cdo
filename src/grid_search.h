#ifndef _GRID_SEARCH_H
#define _GRID_SEARCH_H


struct grid_search {
   void (*grid_search_delete) (struct grid_search * search);
};

typedef struct grid_search grid_search_t;

/*
void do_cell_search(grid_search_t * search, grid_cell_t grid_cell, unsigned * cells_size, unsigned ** cells); 

void do_point_search(grid_search_t * search, grid_point_t grid_point, unsigned * points_size, unsigned ** points);
*/

void grid_search_delete(grid_search_t * search);

grid_search_t *bucket_search_new(grid_t grid_data);

#endif  /* _GRID_SEARCH_H */
