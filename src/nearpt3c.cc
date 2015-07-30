#include "nearpt3x.h"
#include "nearpt3c.h"

void *nearpt3_preprocess(const int nfixpts, Coord_T **pts)
{
  return (void *) nearpt3::Preprocess(nfixpts, pts);
}


int nearpt3_query(void *g, const Coord_T *q)
{
  return nearpt3::Query((nearpt3::Grid_T<Coord_T> *) g, q);
}


void nearpt3_destroy(void *g)
{
  nearpt3::Destroy((nearpt3::Grid_T<Coord_T> *) g);
}
