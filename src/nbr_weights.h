#ifndef NBR_WEIGHTS_H
#define NBR_WEIGHTS_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

class nbrWeightsType
{
 private:
  size_t m_numNeighbors;

 public:
  std::vector<uint8_t> m_mask; // mask at nearest neighbors
  std::vector<size_t> m_add;   // source address at nearest neighbors
  std::vector<double> m_dist;  // angular distance four nearest neighbors

  nbrWeightsType(size_t numNeighbors)
    {
      m_numNeighbors = numNeighbors;
      m_mask.resize(m_numNeighbors);
      m_add.resize(m_numNeighbors);
      m_dist.resize(m_numNeighbors);
    }

  inline void store_distance(size_t nadd, double distance)
  {
    if (m_numNeighbors == 1)
      {
        if (distance < m_dist[0] || (distance <= m_dist[0] && nadd < m_add[0]))
          {
            m_add[0] = nadd;
            m_dist[0] = distance;
          }
      }
    else
      {
        for (size_t i = 0; i < m_numNeighbors; ++i)
          {
            if (distance < m_dist[i] || (distance <= m_dist[i] && nadd < m_add[i]))
              {
                for (size_t n = m_numNeighbors - 1; n > i; --n)
                  {
                    m_add[n] = m_add[n - 1];
                    m_dist[n] = m_dist[n - 1];
                  }
                m_add[i] = nadd;
                m_dist[i] = distance;
                break;
              }
          }
      }
  }

  inline void check_distance()
  {
    constexpr double eps = 1.e-14;
    // If distance is zero, set to small number
    for (size_t i = 0; i < m_numNeighbors; ++i)
      if (m_add[i] < SIZE_MAX && m_dist[i] <= 0.) m_dist[i] = eps;
  }

  inline double compute_weights(const int *restrict src_grid_mask)
  {
    // Compute weights based on inverse distance if mask is false, eliminate those points

    double dist_tot = 0.;  // sum of neighbor distances (for normalizing)

    if (src_grid_mask)
      {
        for (size_t n = 0; n < m_numNeighbors; ++n)
          {
            m_mask[n] = false;
            if (m_add[n] < SIZE_MAX)
              if (src_grid_mask[m_add[n]])
                {
                  m_dist[n] = 1. / m_dist[n];
                  dist_tot += m_dist[n];
                  m_mask[n] = true;
                }
          }
      }
    else
      {
        for (size_t n = 0; n < m_numNeighbors; ++n)
          {
            m_mask[n] = false;
            if (m_add[n] < SIZE_MAX)
              {
                m_dist[n] = 1. / m_dist[n];
                dist_tot += m_dist[n];
                m_mask[n] = true;
              }
          }
      }

    return dist_tot;
  }

  inline size_t normalize_weights(double dist_tot)
  {
    // Normalize weights and store the link
    size_t nadds = 0;
    
    for (size_t n = 0; n < m_numNeighbors; ++n)
      {
        if (m_mask[n])
          {
            m_dist[nadds] = m_dist[n] / dist_tot;
            m_add[nadds] = m_add[n];
            nadds++;
          }
      }

    return nadds;
  }
};

void nbr_store_distance(size_t nadd, double distance, size_t numNeighbors, size_t *restrict nbr_add, double *restrict nbr_dist);
void nbr_check_distance(size_t numNeighbors, const size_t *restrict nbr_add, double *restrict nbr_dist);

#endif
