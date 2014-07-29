#include <stdio.h>
#include <assert.h>
#include "grid_search.h"

void grid_search_delete(grid_search_t * search)
{
  assert( search->grid_search_delete != NULL );

  search->grid_search_delete(search);
}
