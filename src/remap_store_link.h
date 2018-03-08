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
#ifndef REMAP_STORE_LINK_H
#define REMAP_STORE_LINK_H

typedef struct
{
  size_t add;
  double weight;
} addweight_t;

typedef struct
{
  size_t add;
  double weight[4];
} addweight4_t;

typedef struct
{
  size_t nlinks;
  size_t offset;
  addweight_t *addweights;
} weightlinks_t;

typedef struct
{
  size_t nlinks;
  size_t offset;
  addweight4_t *addweights;
} weightlinks4_t;

void store_weightlinks(int lalloc, size_t num_weights, size_t *srch_add, double *weights, size_t cell_add,
                       weightlinks_t *weightlinks);
void store_weightlinks4(size_t num_weights, size_t *srch_add, double weights[4][4], size_t cell_add,
                        weightlinks4_t *weightlinks);
void weightlinks2remaplinks(int lalloc, size_t tgt_grid_size, weightlinks_t *weightlinks, RemapVarsType &rv);
void weightlinks2remaplinks4(size_t tgt_grid_size, weightlinks4_t *weightlinks, RemapVarsType &rv);
void sort_add_and_wgts(size_t num_weights, size_t *src_add, double *wgts);
void sort_add_and_wgts4(size_t num_weights, size_t *src_add, double wgts[4][4]);

#endif /* REMAP_STORE_LINK */
