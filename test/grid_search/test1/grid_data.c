#include <stdio.h>
#include <assert.h>
#include "grid_data.h"

void grid_delete(grid_t * grid)
{
  if ( grid == NULL ) return;

  assert( grid->vtable->delete != NULL );

  grid->vtable->delete(grid);
}
