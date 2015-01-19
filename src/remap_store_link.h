#ifndef _REMAP_STORE_LINK_H
#define _REMAP_STORE_LINK_H


typedef struct
{
  int    add;
  double weight;
} addweight_t;

typedef struct {
  int nlinks;
  int offset;
  addweight_t *addweights;
} weightlinks_t;

void store_weightlinks(long num_weights, int *srch_add, double *weights, long cell_add, weightlinks_t *weightlinks);
void weightlinks2remaplinks(long tgt_grid_size, weightlinks_t *weightlinks,  remapvars_t *rv);

#endif  /* _REMAP_STORE_LINK */
