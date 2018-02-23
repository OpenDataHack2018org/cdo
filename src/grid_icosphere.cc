// modified code from:
// https://schneide.wordpress.com/2016/07/15/generating-an-icosphere-in-c

#include <stdio.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <array>
#include <map>

static inline double
rad2deg(double v)
{
  return v * 180. / M_PI;
}

// Converts cartesian coordinates to geographical
static inline void
cc2gc(double xyz[], double *lon, double *lat)
{
  *lon = atan2(xyz[1], xyz[0]);
  *lat = M_PI_2 - acos(xyz[2]);
}

using Index = size_t;
using vec3 = std::array<double, 3>;

struct Triangle
{
  Index vertex[3];
};

using TriangleList = std::vector<Triangle>;
using VertexList = std::vector<vec3>;

// icosahedron from ICON
namespace icosahedron
{
// Northern hemisphere are the first 6 elements of vertices[0:5]
// Southern hemisphere are the other 6 elements of vertices[6:11]
// 12 vertices
static VertexList vertices(12);
void
init(void)
{
  constexpr double pi_5 = M_PI * 0.2;
  // first define the vertices of the icosahedron
  const double z_w = 2.0 * acos(1.0 / (2.0 * sin(pi_5)));

  // set poles first - it is simple
  vertices[0] = vec3{ 0.0, 0.0, 1.0 };
  vertices[11] = vec3{ 0.0, 0.0, -1.0 };
  // now set the vertices on the two latitude rings
  int i_mdist[10];
  for (int j = 1; j < 11; ++j)
    {
      if (j % 2 == 0)
        i_mdist[j / 2 + 4] = -1 + (j - 1) - 10 * ((j - 1) / 7);
      else
        i_mdist[(j + 1) / 2 - 1] = -1 + (j - 1) - 10 * ((j - 1) / 7);
    }

  for (int j = 1; j < 11; ++j)
    {
      // toggle the hemisphere
      double i_msgn = (j >= 6) ? -1 : 1;
      // compute the meridian angle for the base vertex.
      double z_rlon = (1.0 + i_mdist[j - 1]) * pi_5;
      // now initialize the coordinates
      vertices[j] = vec3{ sin(z_w) * cos(z_rlon), sin(z_w) * sin(z_rlon), cos(z_w) * i_msgn };
    }
}

// 20 triangles
static const TriangleList triangles
    = { { 0, 1, 2 },  { 0, 2, 3 },  { 0, 3, 4 },  { 0, 4, 5 },  { 0, 5, 1 },   { 6, 2, 1 },  { 7, 3, 2 },
        { 8, 4, 3 },  { 9, 5, 4 },  { 10, 1, 5 }, { 2, 6, 7 },  { 3, 7, 8 },   { 4, 8, 9 },  { 5, 9, 10 },
        { 1, 10, 6 }, { 11, 7, 6 }, { 11, 8, 7 }, { 11, 9, 8 }, { 11, 10, 9 }, { 11, 6, 10 } };
}  // namespace icosahedron

static inline vec3
addVector(const vec3 &a, const vec3 &b)
{
  vec3 c;
  for (unsigned i = 0; i < 3; ++i)
    c[i] = a[i] + b[i];
  return c;
}

static inline vec3
subVector(const vec3 &a, const vec3 &b)
{
  vec3 c;
  for (unsigned i = 0; i < 3; ++i)
    c[i] = a[i] - b[i];
  return c;
}

static inline vec3
normalizeVector(const vec3 &a)
{
  vec3 c;
  double magnitude = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
  for (unsigned i = 0; i < 3; ++i)
    c[i] = a[i] / magnitude;
  return c;
}

using Lookup = std::map<std::pair<Index, Index>, Index>;

Index
vertexForEdge(Lookup &lookup, VertexList &vertices, Index first, Index second)
{
  Lookup::key_type key(first, second);
  if (key.first > key.second) std::swap(key.first, key.second);

  auto inserted = lookup.insert({ key, vertices.size() });
  if (inserted.second)
    {
      auto &edge0 = vertices[first];
      auto &edge1 = vertices[second];
      auto point = normalizeVector(addVector(edge0, edge1));
      vertices.push_back(point);
    }

  return inserted.first->second;
}

TriangleList
subdivide(VertexList &vertices, TriangleList triangles)
{
  Lookup lookup;
  TriangleList result;
  Index mid[3];

  for (auto &&each : triangles)
    {
      for (int edge = 0; edge < 3; ++edge)
        mid[edge] = vertexForEdge(lookup, vertices, each.vertex[edge], each.vertex[(edge + 1) % 3]);

      result.push_back({ each.vertex[0], mid[0], mid[2] });
      result.push_back({ each.vertex[1], mid[1], mid[0] });
      result.push_back({ each.vertex[2], mid[2], mid[1] });
      result.push_back({ mid[0], mid[1], mid[2] });
    }

  return result;
}

using IndexedMesh = std::pair<VertexList, TriangleList>;

IndexedMesh
makeIcosphere(int subdivisions)
{
  icosahedron::init();
  VertexList vertices = icosahedron::vertices;
  TriangleList triangles = icosahedron::triangles;

  for (int i = 0; i < subdivisions; ++i)
    triangles = subdivide(vertices, triangles);

  return { vertices, triangles };
}

static inline vec3
vectorProduct(const vec3 &a, const vec3 &b)
{
  vec3 c;
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
  return c;
}

static inline double
dotProduct(const vec3 &a, const vec3 &b)
{
  double sum = 0;
  for (unsigned i = 0; i < 3; ++i)
    sum += a[i] * b[i];
  return sum;
}

static inline void
d_normalize(vec3 &v)
{
  double dnorm = sqrt(dotProduct(v, v));
  for (unsigned i = 0; i < 3; ++i)
    v[i] /= dnorm;
}

vec3
circumCenterMean(const vec3 &v0, const vec3 &v1, const vec3 &v2)
{
  /*
    v0, v1, v2: the coordinates of the three triangle vertices (unit vectors) in
    counter clockwise order center: the coordinates of the circumcenter unless
    co-linear
  */
  // cu0, cu1, cu2: vector product of center:  e1 x e2
  vec3 e1, e2;  // edges of the underlying planar triangle: v1-v0 ands v2-v0,
                // respectively

  e1 = subVector(v1, v0);
  e2 = subVector(v2, v0);
  vec3 cu0 = vectorProduct(e1, e2);
  if (dotProduct(cu0, v0) < 0.0)
    for (unsigned i = 0; i < 3; ++i)
      cu0[i] = -cu0[i];
  d_normalize(cu0);

  e1 = subVector(v2, v1);
  e2 = subVector(v0, v1);
  vec3 cu1 = vectorProduct(e1, e2);
  if (dotProduct(cu1, v1) < 0.0)
    for (unsigned i = 0; i < 3; ++i)
      cu1[i] = -cu1[i];
  d_normalize(cu1);

  e1 = subVector(v0, v2);
  e2 = subVector(v1, v2);
  vec3 cu2 = vectorProduct(e1, e2);
  if (dotProduct(cu2, v2) < 0.0)
    for (unsigned i = 0; i < 3; ++i)
      cu2[i] = -cu2[i];
  d_normalize(cu2);

  vec3 center = addVector(addVector(cu0, cu1), cu2);
  d_normalize(center);

  return center;
}

size_t
genIcosphereCoords(int subdivisions, bool lbounds, double **xvals, double **yvals, double **xbounds, double **ybounds)
{
  IndexedMesh mesh = makeIcosphere(subdivisions);
  VertexList &vertices = mesh.first;
  TriangleList &triangles = mesh.second;

  size_t ncells = triangles.size();
  *xvals = (double *) malloc(ncells * sizeof(double));
  *yvals = (double *) malloc(ncells * sizeof(double));
  if (lbounds)
    {
      *xbounds = (double *) malloc(3 * ncells * sizeof(double));
      *ybounds = (double *) malloc(3 * ncells * sizeof(double));
    }

  size_t i = 0;
  for (Triangle &t : triangles)
    {
      vec3 center = circumCenterMean(vertices[t.vertex[0]], vertices[t.vertex[1]], vertices[t.vertex[2]]);
      cc2gc(&center[0], &(*xvals)[i], &(*yvals)[i]);
      if (lbounds)
        for (size_t k = 0; k < 3; ++k)
          cc2gc(&vertices[t.vertex[k]][0], &(*xbounds)[i * 3 + k], &(*ybounds)[i * 3 + k]);
      i++;
    }

  return ncells;
}

#ifdef TEST_ICO
int
main(void)
{
  IndexedMesh mesh = makeIcosphere(0);
  VertexList &vertices = mesh.first;
  TriangleList &triangles = mesh.second;
  if (1)
    {
      for (vec3 &v : vertices)
        {
          double lon, lat;
          cc2gc(&v[0], &lon, &lat);
          fprintf(stderr, "xyz:%g %g %g   lon:%g lat:%g\n", v[0], v[1], v[2], rad2deg(lon), rad2deg(lat));
        }
      for (Triangle &t : triangles)
        {
          fprintf(stderr, "index: %d %d %d\n", t.vertex[0], t.vertex[1], t.vertex[2]);
        }
      for (Triangle &t : triangles)
        {
          double lon, lat;
          vec3 center = circumCenterMean(vertices[t.vertex[0]], vertices[t.vertex[1]], vertices[t.vertex[2]]);
          cc2gc(&center[0], &lon, &lat);
          fprintf(stderr, "center: %g %g\n", rad2deg(lon), rad2deg(lat));
        }
      for (Triangle &t : triangles)
        {
          double lon, lat;
          printf(">\n");
          for (int i = 0; i < 3; ++i)
            {
              cc2gc(&vertices[t.vertex[i]][0], &lon, &lat);
              printf("   %g  %g\n", rad2deg(lon), rad2deg(lat));
            }
          cc2gc(&vertices[t.vertex[0]][0], &lon, &lat);
          printf("   %g  %g\n", rad2deg(lon), rad2deg(lat));
        }
    }
  fprintf(stderr, "vertices %zu\n", vertices.size());
  fprintf(stderr, "triangles %zu\n", triangles.size());
}
#endif
