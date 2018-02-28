#ifndef NBR_WEIGHTS_H
#define NBR_WEIGHTS_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include <iostream>
#include <vector>

class nbrWeightsType
{
 private:
  size_t m_numNeighbors;

 public:
  std::vector<uint8_t> m_mask; // mask at nearest neighbors
  std::vector<size_t> m_addr;   // source address at nearest neighbors
  std::vector<double> m_dist;  // angular distance four nearest neighbors

  inline void init(size_t numNeighbors)
    {
      m_numNeighbors = numNeighbors;
      m_mask.resize(m_numNeighbors);
      m_addr.resize(m_numNeighbors);
      m_dist.resize(m_numNeighbors);
    }

  nbrWeightsType(size_t numNeighbors)
    {
      init(numNeighbors);
    }

  inline size_t numNeighbors()
    {
      return m_numNeighbors;
    }

  inline void init_addr()
    {
      for (size_t i = 0; i < m_numNeighbors; ++i) m_addr[i] = SIZE_MAX;
    }

  inline void init_dist()
    {
      for (size_t i = 0; i < m_numNeighbors; ++i) m_dist[i] = DBL_MAX;
    }

 inline void store_distance(size_t addr, double distance, size_t numNeighbors)
  {
    assert(numNeighbors <= m_numNeighbors);
    if (numNeighbors == 1)
      {
        if (distance < m_dist[0] || (distance <= m_dist[0] && addr < m_addr[0]))
          {
            m_addr[0] = addr;
            m_dist[0] = distance;
          }
      }
    else
      {
        for (size_t i = 0; i < numNeighbors; ++i)
          {
            if (distance < m_dist[i] || (distance <= m_dist[i] && addr < m_addr[i]))
              {
                for (size_t n = numNeighbors - 1; n > i; --n)
                  {
                    m_addr[n] = m_addr[n - 1];
                    m_dist[n] = m_dist[n - 1];
                  }
                m_addr[i] = addr;
                m_dist[i] = distance;
                break;
              }
          }
      }
  }

  inline void check_distance(size_t numNeighbors)
  {
    assert(numNeighbors <= m_numNeighbors);
    constexpr double eps = 1.e-14;
    // If distance is zero, set to small number
    for (size_t i = 0; i < numNeighbors; ++i)
      if (m_addr[i] < SIZE_MAX && m_dist[i] <= 0.) m_dist[i] = eps;
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
            m_addr[nadds] = m_addr[n];
            nadds++;
          }
      }

    return nadds;
  }

  size_t compute_weights()
  {
    // Compute weights based on inverse distance if mask is false, eliminate those points

    double dist_tot = 0.;  // sum of neighbor distances (for normalizing)

    for (size_t n = 0; n < m_numNeighbors; ++n)
      {
        m_mask[n] = false;
        if (m_addr[n] < SIZE_MAX)
          {
            m_dist[n] = 1. / m_dist[n];
            dist_tot += m_dist[n];
            m_mask[n] = true;
          }
      }

    return normalize_weights(dist_tot);
  }

  size_t compute_weights(const int *restrict src_grid_mask)
  {
    // Compute weights based on inverse distance if mask is false, eliminate those points

    double dist_tot = 0.;  // sum of neighbor distances (for normalizing)

    for (size_t n = 0; n < m_numNeighbors; ++n)
      {
        m_mask[n] = false;
        if (m_addr[n] < SIZE_MAX)
          if (src_grid_mask[m_addr[n]])
            {
              m_dist[n] = 1. / m_dist[n];
              dist_tot += m_dist[n];
              m_mask[n] = true;
            }
      }

    return normalize_weights(dist_tot);
  }
};

void nbr_store_distance(size_t nadd, double distance, size_t numNeighbors, size_t *restrict nbr_add, double *restrict nbr_dist);
void nbr_check_distance(size_t numNeighbors, const size_t *restrict nbr_add, double *restrict nbr_dist);

#endif
