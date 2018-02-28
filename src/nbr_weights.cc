#include "nbr_weights.h"

void
nbr_store_distance(size_t nadd, double distance, size_t numNeighbors, size_t *restrict nbr_add,
                   double *restrict nbr_dist)
{
  if (numNeighbors == 1)
    {
      if (distance < nbr_dist[0] || (distance <= nbr_dist[0] && nadd < nbr_add[0]))
        {
          nbr_add[0] = nadd;
          nbr_dist[0] = distance;
        }
    }
  else
    {
      for (size_t i = 0; i < numNeighbors; ++i)
        {
          if (distance < nbr_dist[i] || (distance <= nbr_dist[i] && nadd < nbr_add[i]))
            {
              for (size_t n = numNeighbors - 1; n > i; --n)
                {
                  nbr_add[n] = nbr_add[n - 1];
                  nbr_dist[n] = nbr_dist[n - 1];
                }
              nbr_add[i] = nadd;
              nbr_dist[i] = distance;
              break;
            }
        }
    }
}

void
nbr_check_distance(size_t numNeighbors, const size_t *restrict nbr_add, double *restrict nbr_dist)
{
  constexpr double eps = 1.e-14;
  // If distance is zero, set to small number
  for (size_t i = 0; i < numNeighbors; ++i)
    if (nbr_add[i] < SIZE_MAX && nbr_dist[i] <= 0.) nbr_dist[i] = eps;
}
