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

#include <algorithm>
#include "cdo_int.h"
#include "remap.h"
#include "remap_store_link.h"

static bool
compareAdds(const Addweight &a, const Addweight &b)
{
  return a.add < b.add;
}

static bool
compareAdds4(const Addweight4 &a, const Addweight4 &b)
{
  return a.add < b.add;
}

static int
qcompareAdds(const void *a, const void *b)
{
  const size_t x = ((const Addweight *) a)->add;
  const size_t y = ((const Addweight *) b)->add;
  return ((x > y) - (x < y)) * 2 + (x > y) - (x < y);
}

static int
qcompareAdds4(const void *a, const void *b)
{
  const size_t x = ((const Addweight4 *) a)->add;
  const size_t y = ((const Addweight4 *) b)->add;
  return ((x > y) - (x < y)) * 2 + (x > y) - (x < y);
}

static void
sortAddweights(size_t numWeights, Addweight *addweights)
{
  size_t n;
  for (n = 1; n < numWeights; ++n)
    if (addweights[n].add < addweights[n - 1].add) break;
  if (n == numWeights) return;

  qsort(addweights, numWeights, sizeof(Addweight), qcompareAdds);
}

static void
sortAddweights4(Addweight4 *addweights)
{
  unsigned n;
  for (n = 1; n < 4; ++n)
    if (addweights[n].add < addweights[n - 1].add) break;
  if (n == 4) return;

  qsort(addweights, 4, sizeof(Addweight4), qcompareAdds4);
}

void
sort_add_and_wgts(size_t numWeights, size_t *src_add, double *wgts)
{
  size_t n;
  for (n = 1; n < numWeights; ++n)
    if (src_add[n] < src_add[n - 1]) break;
  if (n == numWeights) return;

  if (numWeights > 1)
    {
      std::vector<Addweight> addweights(numWeights);

      for (n = 0; n < numWeights; ++n)
        {
          addweights[n].add = src_add[n];
          addweights[n].weight = wgts[n];
        }

      std::sort(addweights.begin(), addweights.end(), compareAdds);

      for (n = 0; n < numWeights; ++n)
        {
          src_add[n] = addweights[n].add;
          wgts[n] = addweights[n].weight;
        }
    }
}

void
sort_add_and_wgts4(size_t numWeights, size_t *src_add, double wgts[4][4])
{
  size_t n;
  for (n = 1; n < numWeights; ++n)
    if (src_add[n] < src_add[n - 1]) break;
  if (n == numWeights) return;

  if (numWeights > 1)
    {
      std::vector<Addweight4> addweights(numWeights);

      for (n = 0; n < numWeights; ++n)
        {
          addweights[n].add = src_add[n];
          for (unsigned k = 0; k < 4; ++k) addweights[n].weight[k] = wgts[n][k];
        }

      std::sort(addweights.begin(), addweights.end(), compareAdds4);

      for (n = 0; n < numWeights; ++n)
        {
          src_add[n] = addweights[n].add;
          for (unsigned k = 0; k < 4; ++k) wgts[n][k] = addweights[n].weight[k];
        }
    }
}

void
storeWeightLinks(int lalloc, size_t numWeights, size_t *srch_add, double *weights, size_t cell_add,
                 std::vector<WeightLinks> &weightLinks)
{
  weightLinks[cell_add].nlinks = 0;
  weightLinks[cell_add].offset = 0;

  if (numWeights)
    {
      Addweight *addweights = NULL;
      if (lalloc)
        addweights = (Addweight *) Malloc(numWeights * sizeof(Addweight));
      else
        addweights = weightLinks[cell_add].addweights;

      for (size_t n = 0; n < numWeights; ++n)
        {
          addweights[n].add = srch_add[n];
          addweights[n].weight = weights[n];
        }

      if (numWeights > 1) sortAddweights(numWeights, addweights);

      weightLinks[cell_add].nlinks = numWeights;

      if (lalloc) weightLinks[cell_add].addweights = addweights;
    }
}

void
storeWeightLinks4(size_t *srch_add, double weights[4][4], size_t cell_add, std::vector<WeightLinks4> &weightLinks)
{
  weightLinks[cell_add].nlinks = 0;
  weightLinks[cell_add].offset = 0;

  Addweight4 *addweights = weightLinks[cell_add].addweights;

  for (unsigned n = 0; n < 4; ++n)
    {
      addweights[n].add = srch_add[n];
      for (unsigned k = 0; k < 4; ++k) addweights[n].weight[k] = weights[n][k];
    }

  sortAddweights4(addweights);

  weightLinks[cell_add].nlinks = 4;
}

void
weightLinksToRemapLinks(int lalloc, size_t gridSize, std::vector<WeightLinks> &weightLinks, RemapVars &rv)
{
  size_t nlinks = 0;

  for (size_t i = 0; i < gridSize; ++i)
    {
      if (weightLinks[i].nlinks)
        {
          weightLinks[i].offset = nlinks;
          nlinks += weightLinks[i].nlinks;
        }
    }

  rv.max_links = nlinks;
  rv.num_links = nlinks;
  if (nlinks)
    {
      rv.src_cell_add.resize(nlinks);
      rv.tgt_cell_add.resize(nlinks);
      rv.wts.resize(nlinks);
      size_t *restrict src_cell_adds = &rv.src_cell_add[0];
      size_t *restrict tgt_cell_adds = &rv.tgt_cell_add[0];
      double *restrict wts = &rv.wts[0];

#ifdef _OPENMP
#pragma omp parallel for schedule(static) default(none) shared(src_cell_adds, tgt_cell_adds, wts, weightLinks, gridSize)
#endif
      for (size_t i = 0; i < gridSize; ++i)
        {
          size_t num_links = weightLinks[i].nlinks;
          if (num_links)
            {
              size_t offset = weightLinks[i].offset;
              Addweight *addweights = weightLinks[i].addweights;
              for (size_t ilink = 0; ilink < num_links; ++ilink)
                {
                  src_cell_adds[offset + ilink] = addweights[ilink].add;
                  tgt_cell_adds[offset + ilink] = i;
                  wts[offset + ilink] = addweights[ilink].weight;
                }
            }
        }

      if (lalloc)
        {
          for (size_t i = 0; i < gridSize; ++i)
            {
              size_t num_links = weightLinks[i].nlinks;
              if (num_links) Free(weightLinks[i].addweights);
            }
        }
      else
        {
          Free(weightLinks[0].addweights);
        }
    }
}

void
weightLinks4ToRemapLinks(size_t gridSize, std::vector<WeightLinks4> &weightLinks, RemapVars &rv)
{
  size_t nlinks = 0;

  for (size_t i = 0; i < gridSize; ++i)
    {
      if (weightLinks[i].nlinks)
        {
          weightLinks[i].offset = nlinks;
          nlinks += weightLinks[i].nlinks;
        }
    }

  rv.max_links = nlinks;
  rv.num_links = nlinks;
  if (nlinks)
    {
      rv.src_cell_add.resize(nlinks);
      rv.tgt_cell_add.resize(nlinks);
      rv.wts.resize(4 * nlinks);
      size_t *restrict src_cell_adds = &rv.src_cell_add[0];
      size_t *restrict tgt_cell_adds = &rv.tgt_cell_add[0];
      double *restrict wts = &rv.wts[0];

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(src_cell_adds, tgt_cell_adds, wts, weightLinks, gridSize)
#endif
      for (size_t i = 0; i < gridSize; ++i)
        {
          size_t num_links = weightLinks[i].nlinks;
          if (num_links)
            {
              size_t offset = weightLinks[i].offset;
              Addweight4 *addweights = weightLinks[i].addweights;
              for (size_t ilink = 0; ilink < num_links; ++ilink)
                {
                  src_cell_adds[offset + ilink] = addweights[ilink].add;
                  tgt_cell_adds[offset + ilink] = i;
                  for (size_t k = 0; k < 4; ++k) wts[(offset + ilink) * 4 + k] = addweights[ilink].weight[k];
                }
            }
        }

      Free(weightLinks[0].addweights);
    }
}

void
weightLinksAlloc(size_t numNeighbors, size_t gridSize, std::vector<WeightLinks> &weightLinks)
{
  weightLinks[0].addweights = (Addweight *) Malloc(numNeighbors * gridSize * sizeof(Addweight));
  for (size_t i = 1; i < gridSize; ++i) weightLinks[i].addweights = weightLinks[0].addweights + numNeighbors * i;
}

void
weightLinks4Alloc(size_t gridSize, std::vector<WeightLinks4> &weightLinks)
{
  weightLinks[0].addweights = (Addweight4 *) Malloc(4 * gridSize * sizeof(Addweight4));
  for (size_t i = 1; i < gridSize; ++i) weightLinks[i].addweights = weightLinks[0].addweights + 4 * i;
}
