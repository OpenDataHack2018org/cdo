#include "cdo_int.h"
#include "cdoOptions.h"
#include "remap_vars.h"
#include "timer.h"

/*
  -----------------------------------------------------------------------

  Performs the remapping based on weights computed elsewhere

  -----------------------------------------------------------------------
*/
void
remap(double *restrict dst_array, double missval, size_t dst_size, const remapVarsType &rv, const double *restrict src_array,
      const double *restrict src_grad1, const double *restrict src_grad2, const double *restrict src_grad3)
{
  /*
    Input arrays:

    int *dst_add         ! destination address for each link
    int *src_add         ! source      address for each link

    int num_wts          ! num of weights used in remapping

    double *map_wts      ! remapping weights for each link

    double *src_array    ! array with source field to be remapped

    optional:

    double *src_grad1    ! gradient arrays on source grid necessary for
    double *src_grad2    ! higher-order remappings
    double *src_grad3

    output variables:

    double *dst_array    ! array for remapped field on destination grid
  */
  size_t num_links = rv.num_links;
  size_t num_wts = rv.num_wts;
  const double *restrict map_wts = &rv.wts[0];
  const size_t *restrict dst_add = &rv.tgt_cell_add[0];
  const size_t *restrict src_add = &rv.src_cell_add[0];
  remaplink_t links = rv.links;
  long links_per_value = rv.links_per_value;
  
  extern int timer_remap;

  // Check the order of the interpolation

  int iorder = (src_grad1 == NULL) ? 1 : 2;

#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(none) shared(dst_size, dst_array, missval)
#endif
  for (size_t n = 0; n < dst_size; ++n) dst_array[n] = missval;

  if (cdoTimer) timer_start(timer_remap);

  if (iorder == 1)  // First order remapping
    {
      if (links.option)
        {
#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(none) shared(num_links, dst_array, dst_add)
#endif
          for (size_t n = 0; n < num_links; ++n) dst_array[dst_add[n]] = 0.;

          for (size_t j = 0; j < links.num_blks; ++j)
            {
              const size_t *restrict dst_addx = links.dst_add[j];
              const size_t *restrict src_addx = links.src_add[j];
              const size_t *restrict windex = links.w_index[j];
              size_t nlinks = links.num_links[j];

#ifdef HAVE_OPENMP4
#pragma omp parallel for simd default(none) shared(nlinks, dst_array, src_array, dst_addx, src_addx, map_wts, num_wts, windex)
#endif
              for (size_t n = 0; n < nlinks; ++n)
                {
                  dst_array[dst_addx[n]] += src_array[src_addx[n]] * map_wts[num_wts * windex[n]];
                }
            }
        }
      else
        {
          long lpv = links_per_value;
          if (lpv > 0)
            {
              size_t nlinks = num_links / lpv;

              if (lpv == 4)
                {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(dst_array, src_array, dst_add, src_add, map_wts, num_wts, nlinks, lpv)
#endif
                  for (size_t n = 0; n < nlinks; ++n)
                    {
                      size_t noff = n * lpv;
                      dst_array[dst_add[noff]] = src_array[src_add[noff]] * map_wts[num_wts * (noff)]
                                                 + src_array[src_add[noff + 1]] * map_wts[num_wts * (noff + 1)]
                                                 + src_array[src_add[noff + 2]] * map_wts[num_wts * (noff + 2)]
                                                 + src_array[src_add[noff + 3]] * map_wts[num_wts * (noff + 3)];
                    }
                }
              else
                {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(dst_array, src_array, dst_add, src_add, map_wts, num_wts, nlinks, lpv)
#endif
                  for (size_t n = 0; n < nlinks; ++n)
                    {
                      size_t noff = n * lpv;
                      dst_array[dst_add[noff]] = src_array[src_add[noff]] * map_wts[num_wts * noff];
                      for (size_t k = 1; k < (size_t) lpv; ++k)
                        dst_array[dst_add[noff]] += src_array[src_add[noff + k]] * map_wts[num_wts * (noff + k)];
                    }
                }
            }
          else
            {
#ifdef SX
#pragma cdir nodep
#endif
              for (size_t n = 0; n < num_links; ++n) dst_array[dst_add[n]] = 0.;

              for (size_t n = 0; n < num_links; ++n)
                {
                  // printf("%5d %5d %5ld %g # dst_add src_add n\n", dst_add[n],
                  // src_add[n], n, map_wts[num_wts*n]);
                  dst_array[dst_add[n]] += src_array[src_add[n]] * map_wts[num_wts * n];
                }
            }
        }
    }
  else  // Second order remapping
    {
#ifdef SX
#pragma cdir nodep
#endif
      for (size_t n = 0; n < num_links; ++n) dst_array[dst_add[n]] = 0.;

      if (num_wts == 3)
        {
          for (size_t n = 0; n < num_links; ++n)
            {
              dst_array[dst_add[n]] += src_array[src_add[n]] * map_wts[3 * n]
                                       + src_grad1[src_add[n]] * map_wts[3 * n + 1]
                                       + src_grad2[src_add[n]] * map_wts[3 * n + 2];
            }
        }
      else if (num_wts == 4)
        {
          for (size_t n = 0; n < num_links; ++n)
            {
              dst_array[dst_add[n]]
                  += src_array[src_add[n]] * map_wts[4 * n] + src_grad1[src_add[n]] * map_wts[4 * n + 1]
                     + src_grad2[src_add[n]] * map_wts[4 * n + 2] + src_grad3[src_add[n]] * map_wts[4 * n + 3];
            }
        }
    }

  if (cdoTimer) timer_stop(timer_remap);
}

static size_t
get_max_add(size_t num_links, size_t size, const size_t *restrict add)
{
  size_t *isum = (size_t *) Calloc(size, sizeof(size_t));

  for (size_t n = 0; n < num_links; ++n) isum[add[n]]++;

  size_t max_add = 0;
  for (size_t i = 0; i < size; ++i)
    if (isum[i] > max_add) max_add = isum[i];
  Free(isum);

  return max_add;
}

static size_t
binary_search_int(const size_t *array, size_t len, size_t value)
{
  int64_t low = 0, high = len - 1, midpoint = 0;

  while (low <= high)
    {
      midpoint = low + (high - low) / 2;

      // check to see if value is equal to item in array
      if (value == array[midpoint]) return midpoint;

      if (value < array[midpoint])
        high = midpoint - 1;
      else
        low = midpoint + 1;
    }

  // item was not found
  return len;
}

/*
  -----------------------------------------------------------------------

  Performs the remapping based on weights computed elsewhere

  -----------------------------------------------------------------------
*/
void
remap_laf(double *restrict dst_array, double missval, size_t dst_size, const remapVarsType &rv, const double *restrict src_array)
{
  /*
    Input arrays:

    int *dst_add         ! destination address for each link
    int *src_add         ! source      address for each link

    int num_wts          ! num of weights used in remapping

    double *map_wts      ! remapping weights for each link

    double *src_array    ! array with source field to be remapped

    output variables:

    double *dst_array    ! array for remapped field on destination grid
  */
  size_t num_links = rv.num_links;
  size_t num_wts = rv.num_wts;
  const double *restrict map_wts = &rv.wts[0];
  const size_t *restrict dst_add = &rv.tgt_cell_add[0];
  const size_t *restrict src_add = &rv.src_cell_add[0];

  arrayFill(dst_size, dst_array, missval);

  if (num_links == 0) return;

  size_t max_cls = get_max_add(num_links, dst_size, dst_add);

#ifdef _OPENMP
  VECTOR_2D(double, src_cls2, Threading::ompNumThreads, max_cls);
  VECTOR_2D(double, src_wts2, Threading::ompNumThreads, max_cls);
#else
  std::vector<double> src_cls(max_cls);
  std::vector<double> src_wts(max_cls);
#endif

  for (size_t n = 0; n < num_links; ++n)
    if (DBL_IS_EQUAL(dst_array[dst_add[n]], missval)) dst_array[dst_add[n]] = 0.0;

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(dst_size, src_cls2, src_wts2, num_links, dst_add, src_add, src_array, \
                                              map_wts, num_wts, dst_array, max_cls) schedule(dynamic, 1)
#endif
  for (size_t i = 0; i < dst_size; ++i)
    {
      size_t k;
      size_t ncls;
#ifdef _OPENMP
      int ompthID = cdo_omp_get_thread_num();
      double *src_cls = &src_cls2[ompthID][0];
      double *src_wts = &src_wts2[ompthID][0];
#endif
      arrayFill(max_cls, &src_cls[0], 0.0);
      arrayFill(max_cls, &src_wts[0], 0.0);
      /*
      ncls = 0;
      for ( n = 0; n < num_links; n++ )
        {
          if ( i == dst_add[n] )
            {
              for ( k = 0; k < ncls; k++ )
                if ( IS_EQUAL(src_array[src_add[n]], src_cls[k]) ) break;

              if ( k == ncls )
                {
                  src_cls[k] = src_array[src_add[n]];
                  ncls++;
                }

              src_wts[k] += map_wts[num_wts*n];
            }
        }
      */
      /* only for sorted dst_add! */
      {
        size_t min_add = 1, max_add = 0;

        size_t n = binary_search_int(dst_add, num_links, i);

        if (n < num_links)
          {
            min_add = n;

            for (n = min_add + 1; n < num_links; ++n)
              if (i != dst_add[n]) break;

            max_add = n;

            for (n = min_add; n > 0; --n)
              if (i != dst_add[n - 1]) break;

            min_add = n;
          }

        ncls = 0;
        for (n = min_add; n < max_add; ++n)
          {
            for (k = 0; k < ncls; ++k)
              if (IS_EQUAL(src_array[src_add[n]], src_cls[k])) break;

            if (k == ncls)
              {
                src_cls[k] = src_array[src_add[n]];
                ncls++;
              }

            src_wts[k] += map_wts[num_wts * n];
          }
      }

      if (ncls)
        {
          size_t imax = 0;
          double wts = src_wts[0];
          for (k = 1; k < ncls; ++k)
            {
              if (src_wts[k] > wts)
                {
                  wts = src_wts[k];
                  imax = k;
                }
            }

          dst_array[i] = src_cls[imax];
        }
    }
}

/*
  -----------------------------------------------------------------------

  Performs the remapping based on weights computed elsewhere

  -----------------------------------------------------------------------
*/
void
remap_sum(double *restrict dst_array, double missval, size_t dst_size, const remapVarsType &rv, const double *restrict src_array)
{
  /*
    Input arrays:

    int *dst_add         ! destination address for each link
    int *src_add         ! source      address for each link

    int num_wts          ! num of weights used in remapping

    double *map_wts      ! remapping weights for each link
    double *src_array    ! array with source field to be remapped

    output variables:

    double *dst_array    ! array for remapped field on destination grid
  */
  size_t num_links = rv.num_links;
  size_t num_wts = rv.num_wts;
  const double *restrict map_wts = &rv.wts[0];
  const size_t *restrict dst_add = &rv.tgt_cell_add[0];
  const size_t *restrict src_add = &rv.src_cell_add[0];

  for (size_t n = 0; n < dst_size; ++n) dst_array[n] = missval;

#ifdef SX
#pragma cdir nodep
#endif
  for (size_t n = 0; n < num_links; ++n)
    if (DBL_IS_EQUAL(dst_array[dst_add[n]], missval)) dst_array[dst_add[n]] = 0.0;

  for (size_t n = 0; n < num_links; ++n)
    {
      /*
        printf("%5d %5d %5d %g # dst_add src_add n\n", dst_add[n], src_add[n],
        n, map_wts[num_wts*n]);
      */
      // dst_array[dst_add[n]] += src_array[src_add[n]]*map_wts[num_wts*n];
      dst_array[dst_add[n]] += src_array[src_add[n]] * map_wts[num_wts * n];
      printf("%zu %zu %zu %g %g %g\n", n, dst_add[n], src_add[n], src_array[src_add[n]], map_wts[num_wts * n],
             dst_array[dst_add[n]]);
    }
}

/*
    This routine initializes some variables and provides an initial
    allocation of arrays (fairly large so frequent resizing unnecessary).
*/
void
remap_vars_init(RemapType mapType, size_t src_grid_size, size_t tgt_grid_size, remapVarsType &rv)
{
  /* Initialize all pointer */
  if (rv.pinit == false) rv.pinit = true;

  /* Determine the number of weights */

  if (mapType == RemapType::CONSERV)
    rv.num_wts = 3;
  else if (mapType == RemapType::CONSERV_YAC)
    rv.num_wts = 1;
  else if (mapType == RemapType::BILINEAR)
    rv.num_wts = 1;
  else if (mapType == RemapType::BICUBIC)
    rv.num_wts = 4;
  else if (mapType == RemapType::DISTWGT)
    rv.num_wts = 1;
  else
    cdoAbort("Unknown mapping method!");

  rv.sort_add = (mapType == RemapType::CONSERV);

  rv.links_per_value = -1;

  /*
   Initialize num_links and set max_links to four times the largest
   of the destination grid sizes initially (can be changed later).
   Set a default resize increment to increase the size of link
   arrays if the number of links exceeds the initial size
 */
  rv.num_links = 0;
  rv.max_links = 4 * tgt_grid_size;

  rv.resize_increment = (size_t)(0.1 * MAX(src_grid_size, tgt_grid_size));

  /*  Allocate address and weight arrays for mapping 1 */
  if (mapType == RemapType::CONSERV)
    {
      rv.src_cell_add.resize(rv.max_links);
      rv.tgt_cell_add.resize(rv.max_links);
      rv.wts.resize(rv.num_wts * rv.max_links);
    }

  rv.links.option = false;
  rv.links.max_links = 0;
  rv.links.num_blks = 0;
  rv.links.num_links = NULL;
  rv.links.src_add = NULL;
  rv.links.dst_add = NULL;
  rv.links.w_index = NULL;

} /* remapVarsInit */

/*
   This routine resizes remapping arrays by increasing(decreasing) the max_links by increment
*/
void
resize_remap_vars(remapVarsType &rv, int64_t increment)
{
  /*
    Input variables:
    int  increment  ! the number of links to add(subtract) to arrays
  */

  /*  Reallocate arrays at new size */

  rv.max_links += increment;

  if (rv.max_links)
    {
      rv.src_cell_add.resize(rv.max_links);
      rv.tgt_cell_add.resize(rv.max_links);
      rv.wts.resize(rv.num_wts * rv.max_links);
    }

} /* resize_remap_vars */

void
remapVarsReorder(remapVarsType &rv)
{
  size_t j, nval = 0, num_blks = 0;
  size_t n;

  size_t num_links = rv.num_links;

  printf("remapVarsReorder\n");
  printf("  num_links %zu\n", num_links);
  rv.links.option = true;

  size_t lastval = -1;
  size_t max_links = 0;
  for (n = 0; n < num_links; n++)
    {
      if (rv.tgt_cell_add[n] == lastval)
        nval++;
      else
        {
          if (nval > num_blks) num_blks = nval;
          nval = 1;
          max_links++;
          lastval = rv.tgt_cell_add[n];
        }
    }

  if (num_blks)
    {
      rv.links.max_links = max_links;
      rv.links.num_blks = num_blks;

      printf("num_links %zu  max_links %zu  num_blks %zu\n", rv.num_links, max_links, num_blks);

      rv.links.num_links = (size_t *) Malloc(num_blks * sizeof(size_t));
      rv.links.dst_add = (size_t **) Malloc(num_blks * sizeof(size_t *));
      rv.links.src_add = (size_t **) Malloc(num_blks * sizeof(size_t *));
      rv.links.w_index = (size_t **) Malloc(num_blks * sizeof(size_t *));
    }

  for (j = 0; j < num_blks; j++)
    {
      rv.links.dst_add[j] = (size_t *) Malloc(max_links * sizeof(size_t));
      rv.links.src_add[j] = (size_t *) Malloc(max_links * sizeof(size_t));
      rv.links.w_index[j] = (size_t *) Malloc(max_links * sizeof(size_t));
    }

  for (j = 0; j < num_blks; j++)
    {
      nval = 0;
      lastval = -1;
      size_t nlinks = 0;

      for (n = 0; n < num_links; n++)
        {
          if (rv.tgt_cell_add[n] == lastval)
            nval++;
          else
            {
              nval = 1;
              lastval = rv.tgt_cell_add[n];
            }

          if (nval == j + 1)
            {
              rv.links.dst_add[j][nlinks] = rv.tgt_cell_add[n];
              rv.links.src_add[j][nlinks] = rv.src_cell_add[n];
              rv.links.w_index[j][nlinks] = n;
              nlinks++;
            }
        }

      rv.links.num_links[j] = nlinks;
      printf("loop %zu  nlinks %zu\n", j + 1, nlinks);
    }
}

void
remapVarsFree(remapVarsType &rv)
{
  if (rv.pinit)
    {
      rv.pinit = false;
      rv.sort_add = false;

      rv.src_cell_add.resize(0);
      rv.tgt_cell_add.resize(0);
      rv.wts.resize(0);

      if (rv.links.option)
        {
          rv.links.option = false;

          if (rv.links.num_blks)
            {
              Free(rv.links.num_links);
              size_t num_blks = rv.links.num_blks;
              for (size_t i = 0; i < num_blks; ++i)
                {
                  Free(rv.links.src_add[i]);
                  Free(rv.links.dst_add[i]);
                  Free(rv.links.w_index[i]);
                }
              Free(rv.links.src_add);
              Free(rv.links.dst_add);
              Free(rv.links.w_index);
            }
        }
    }
  else
    fprintf(stderr, "%s Warning: vars not initialized!\n", __func__);

} /* remapVarsFree */

void
remapVarsCheckWeights(const remapVarsType &rv)
{
  auto num_links = rv.num_links;
  auto num_wts = rv.num_wts;
  auto normOpt = rv.normOpt;
  auto src_cell_add = &rv.src_cell_add[0];
  auto tgt_cell_add = &rv.tgt_cell_add[0];
  auto wts = &rv.wts[0];

  for (size_t n = 0; n < num_links; ++n)
    {
      if (wts[n * num_wts] < -0.01)
        cdoPrint("Map weight < 0! grid1idx=%zu grid2idx=%zu nlink=%zu wts=%g", src_cell_add[n], tgt_cell_add[n], n,
                 wts[n * num_wts]);

      if (normOpt != NormOpt::NONE && wts[n * num_wts] > 1.01)
        cdoPrint("Map weight > 1! grid1idx=%zu grid2idx=%zu nlink=%zu wts=%g", src_cell_add[n], tgt_cell_add[n], n,
                 wts[n * num_wts]);
    }
}
