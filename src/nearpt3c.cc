#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#if defined(ENABLE_NEARPT3)
#include "nearpt3x.h"
#endif

#include "nearpt3c.h"

void *nearpt3_preprocess(const int nfixpts, Coord_T **pts)
{
#if defined(ENABLE_NEARPT3)
  return (void *) nearpt3::Preprocess(nfixpts, pts);
#endif
}


int nearpt3_query(void *g, const Coord_T *q)
{
#if defined(ENABLE_NEARPT3)
  return nearpt3::Query((nearpt3::Grid_T<Coord_T> *) g, q);
#endif
}


void nearpt3_destroy(void *g)
{
#if defined(ENABLE_NEARPT3)
  nearpt3::Destroy((nearpt3::Grid_T<Coord_T> *) g);
#endif
}
