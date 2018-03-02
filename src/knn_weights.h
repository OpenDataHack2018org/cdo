#ifndef NBR_WEIGHTS_H
#define NBR_WEIGHTS_H

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include <vector>

double intlin(double x, double y1, double x1, double y2, double x2);

class knnWeightsType
{
private:
  size_t m_numNeighbors;
  size_t m_maxNeighbors;

public:
  std::vector<uint8_t> m_mask;  // mask at nearest neighbors
  std::vector<size_t> m_addr;   // source address at nearest neighbors
  std::vector<double> m_dist;   // angular distance four nearest neighbors
  std::vector<size_t> m_tmpaddr;
  std::vector<double> m_tmpdist;

  inline void
  init(size_t maxNeighbors)
  {
    m_numNeighbors = 0;
    m_maxNeighbors = maxNeighbors;
    m_mask.resize(m_maxNeighbors);
    m_addr.resize(m_maxNeighbors);
    m_dist.resize(m_maxNeighbors);
  }

  knnWeightsType(size_t maxNeighbors) { init(maxNeighbors); }

  inline size_t
  maxNeighbors()
  {
    return m_maxNeighbors;
  }

  inline size_t
  numNeighbors()
  {
    return m_numNeighbors;
  }

  inline void
  init_addr()
  {
    for (size_t i = 0; i < m_maxNeighbors; ++i) m_addr[i] = SIZE_MAX;
  }

  inline void
  init_dist()
  {
    for (size_t i = 0; i < m_maxNeighbors; ++i) m_dist[i] = DBL_MAX;
  }

  inline void
  store_distance(size_t addr, double distance, size_t numNeighbors)
  {
    assert(numNeighbors <= m_maxNeighbors);
    m_numNeighbors = numNeighbors;

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

  inline void
  check_distance()
  {
    constexpr double eps = 1.e-14;
    // If distance is zero, set to small number
    for (size_t i = 0; i < m_numNeighbors; ++i)
      if (m_addr[i] < SIZE_MAX && m_dist[i] <= 0.) m_dist[i] = eps;
  }

  size_t
  normalize_weights(double dist_tot)
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

    m_numNeighbors = nadds;
    return nadds;
  }

  size_t
  compute_weights()
  {
    // Compute weights based on inverse distance if mask is false, eliminate those points

    double dist_tot = 0.;  // sum of neighbor distances (for normalizing)

    for (size_t n = 0; n < m_maxNeighbors; ++n)
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

  size_t
  compute_weights(const int *src_grid_mask)
  {
    // Compute weights based on inverse distance if mask is false, eliminate those points

    double dist_tot = 0.;  // sum of neighbor distances (for normalizing)

    for (size_t n = 0; n < m_maxNeighbors; ++n)
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

  size_t
  compute_weights(const uint8_t *src_grid_mask, double search_radius, double weight0, double weightR)
  {
    // Compute weights based on inverse distance if mask is false, eliminate those points

    double dist_tot = 0.;  // sum of neighbor distances (for normalizing)

    for (size_t n = 0; n < m_maxNeighbors; ++n)
      {
        m_mask[n] = false;
        if (m_addr[n] < SIZE_MAX)
          if (src_grid_mask[m_addr[n]])
            {
              m_dist[n] = intlin(m_dist[n], weight0, 0, weightR, search_radius);
              dist_tot += m_dist[n];
              m_mask[n] = true;
            }
      }

    return normalize_weights(dist_tot);
  }

  double
  array_weights_sum(const double *array)
  {
    double result = 0;
    for (size_t n = 0; n < m_numNeighbors; ++n) result += array[m_addr[n]] * m_dist[n];
    return result;
  }
};

#endif
