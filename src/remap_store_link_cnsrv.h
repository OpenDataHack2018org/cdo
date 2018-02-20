/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2018 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/
#ifndef REMAP_STORE_LINK_CNSRV_H
#define REMAP_STORE_LINK_CNSRV_H

extern int remap_store_link_fast;

// used for store_link_fast

#define BLK_SIZE 4096
#define BLK_NUM(x) (x / grid_store->blk_size)
#define BLK_IDX(x) (x % grid_store->blk_size)
struct grid_layer
{
  long *grid2_link;
  struct grid_layer *next;
};

typedef struct grid_layer grid_layer_t;

typedef struct
{
  long blk_size;
  long max_size;
  long nblocks;
  long *blksize;
  long *nlayers;
  grid_layer_t **layers;
} grid_store_t;

void grid_store_init(grid_store_t *grid_store, long gridsize);
void grid_store_delete(grid_store_t *grid_store);

void store_link_cnsrv_fast(remapvars_t *rv, long add1, long add2, long num_wts,
                           double *weights, grid_store_t *grid_store);
void store_link_cnsrv(remapvars_t *rv, long add1, long add2,
                      double *restrict weights, long *link_add1[2],
                      long *link_add2[2]);

#endif /* REMAP_STORE_LINK_CNSRV */
