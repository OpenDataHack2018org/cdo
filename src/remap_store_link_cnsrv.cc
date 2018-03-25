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

#include "cdo_int.h"
#include "remap.h"
#include "remap_store_link_cnsrv.h"

void
grid_store_init(grid_store_t *grid_store, long gridsize)
{
  long iblk;
  long blksize[] = { 128, 256, 512, 1024, 2048, 4096, 8192 };
  long nblks = sizeof(blksize) / sizeof(long);
  long nblocks;

  for (iblk = nblks - 1; iblk >= 0; --iblk)
    if (gridsize / blksize[iblk] > 99) break;

  if (iblk < 0) iblk = 0;

  /* grid_store->blk_size = BLK_SIZE; */
  grid_store->blk_size = blksize[iblk];
  grid_store->max_size = gridsize;

  grid_store->nblocks = grid_store->max_size / grid_store->blk_size;
  if (grid_store->max_size % grid_store->blk_size > 0) grid_store->nblocks++;

  if (cdoVerbose)
    fprintf(stdout, "blksize = %ld  lastblksize = %ld  max_size = %ld  nblocks = %ld\n", grid_store->blk_size,
            grid_store->max_size % grid_store->blk_size, grid_store->max_size, grid_store->nblocks);

  grid_store->blksize = (long *) Malloc(grid_store->nblocks * sizeof(long));
  grid_store->nlayers = (long *) Malloc(grid_store->nblocks * sizeof(long));
  grid_store->layers = (grid_layer_t **) Malloc(grid_store->nblocks * sizeof(grid_layer_t *));

  nblocks = grid_store->nblocks;
  for (iblk = 0; iblk < nblocks; ++iblk)
    {
      grid_store->blksize[iblk] = grid_store->blk_size;
      grid_store->nlayers[iblk] = 0;
      grid_store->layers[iblk] = NULL;
    }
  if (grid_store->max_size % grid_store->blk_size > 0)
    grid_store->blksize[grid_store->nblocks - 1] = grid_store->max_size % grid_store->blk_size;
}

void
grid_store_delete(grid_store_t *grid_store)
{
  long nblocks = grid_store->nblocks;

  for (long iblk = 0; iblk < nblocks; ++iblk)
    {
      long j = 0;
      grid_layer_t *grid_layer = grid_store->layers[iblk];
      long nlayers = grid_store->nlayers[iblk];
      long blksize = grid_store->blksize[iblk];
      for (long ilayer = 0; ilayer < nlayers; ++ilayer)
        {
          if (cdoVerbose)
            {
              for (long i = 0; i < blksize; ++i)
                if (grid_layer->grid2_link[i] != -1) j++;
            }

          grid_layer_t *grid_layer_f = grid_layer;
          Free(grid_layer->grid2_link);
          grid_layer = grid_layer->next;
          Free(grid_layer_f);
        }
      /*
      if ( cdoVerbose )
        {
          fprintf(stderr, "block = %ld nlayers = %ld  allocated = %ld  used =
      %ld\n", iblk+1, nlayers, nlayers*blksize, j);
        }
      */
    }

  Free(grid_store->blksize);
  Free(grid_store->layers);
  Free(grid_store->nlayers);
}

/*
    This routine stores the address and weight for this link in the appropriate
    address and weight arrays and resizes those arrays if necessary.
*/
void
store_link_cnsrv(RemapVars &rv, long add1, long add2, long num_wts, double *weights, grid_store_t *grid_store)
{
  /*
    Input variables:
    int  add1         ! address on source grid
    int  add2         ! address on target grid
    double weights[]  ! array of remapping weights for this link
  */
  /* Local variables */
  long nlink; /* link index */
  grid_layer_t *grid_layer = NULL;

  /*  If all weights are ZERO, do not bother storing the link */
  if (num_wts == 3)
    {
      if (IS_EQUAL(weights[0], 0) && IS_EQUAL(weights[1], 0) && IS_EQUAL(weights[2], 0)) return;
    }
  else
    {
      if (IS_EQUAL(weights[0], 0)) return;
    }

  /* If the link already exists, add the weight to the current weight arrays */

  long iblk = BLK_NUM(add2);
  long iadd2 = BLK_IDX(add2);

  bool lstore_link = false;
  grid_layer_t **grid_layer2 = &grid_store->layers[iblk];
  long nlayer = grid_store->nlayers[iblk];
  long ilayer;
  for (ilayer = 0; ilayer < nlayer; ++ilayer)
    {
      grid_layer = *grid_layer2;
      nlink = grid_layer->grid2_link[iadd2];
      if (nlink == -1)
        {
          break;
        }
      else if ((size_t) add1 == rv.src_cell_add[nlink])
        {
          lstore_link = true;
          break;
        }
      grid_layer2 = &(*grid_layer2)->next;
    }

  if (lstore_link)
    {
      for (long i = 0; i < num_wts; ++i) rv.wts[num_wts * nlink + i] += weights[i];
      return;
    }

  /*
     If the link does not yet exist, increment number of links and
     check to see if remap arrays need to be increased to accomodate
     the new link. Then store the link.
  */
  nlink = rv.num_links;

  if (ilayer < grid_store->nlayers[iblk])
    {
      grid_layer->grid2_link[iadd2] = nlink;
    }
  else
    {
      grid_layer = (grid_layer_t *) Malloc(sizeof(grid_layer_t));
      grid_layer->next = NULL;
      grid_layer->grid2_link = (long *) Malloc(grid_store->blksize[iblk] * sizeof(long));

      long blksize = grid_store->blksize[iblk];
      for (long i = 0; i < blksize; ++i) grid_layer->grid2_link[i] = -1;

      grid_layer->grid2_link[iadd2] = nlink;
      *grid_layer2 = grid_layer;
      grid_store->nlayers[iblk]++;
    }

  rv.num_links++;
  remapVarsEnsureSize(rv, rv.num_links);

  rv.src_cell_add[nlink] = add1;
  rv.tgt_cell_add[nlink] = add2;

  for (long i = 0; i < num_wts; ++i) rv.wts[num_wts * nlink + i] = weights[i];

} /* store_link_cnsrv */
