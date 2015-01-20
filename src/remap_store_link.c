#include "cdo.h"
#include "cdo_int.h"
#include "remap.h"
#include "remap_store_link.h"


static
int cmp_adds(const void *s1, const void *s2)
{
  int cmp = 0;
  const addweight_t* c1 = (const addweight_t*) s1;
  const addweight_t* c2 = (const addweight_t*) s2;

  if      ( c1->add < c2->add ) cmp = -1;
  else if ( c1->add > c2->add ) cmp =  1;

  return (cmp);
}

static
void sort_addweights(size_t num_weights, addweight_t *addweights)
{
  size_t n;

  for ( n = 1; n < num_weights; ++n )
    if ( addweights[n].add < addweights[n-1].add ) break;
  if ( n == num_weights ) return;

  qsort(addweights, num_weights, sizeof(addweight_t), cmp_adds);
}


void sort_add_and_wgts(size_t num_weights, int *src_add, double *wgts)
{
  size_t n;

  for ( n = 1; n < num_weights; ++n )
    if ( src_add[n] < src_add[n-1] ) break;
  if ( n == num_weights ) return;

  addweight_t addweights[num_weights];

  for ( n = 0; n < num_weights; ++n )
    {
      addweights[n].add    = src_add[n];
      addweights[n].weight = wgts[n];
    }

  qsort(addweights, num_weights, sizeof(addweight_t), cmp_adds);

  for ( n = 0; n < num_weights; ++n )
    {
      src_add[n] = addweights[n].add;
      wgts[n]    = addweights[n].weight;
    }  
}


void store_weightlinks(long num_weights, int *srch_add, double *weights, long cell_add, weightlinks_t *weightlinks)
{
  weightlinks[cell_add].nlinks = 0;
  weightlinks[cell_add].offset = 0;

  if ( num_weights )
    {
      addweight_t *addweights = (addweight_t *) malloc(num_weights*sizeof(addweight_t));
      for ( long n = 0; n < num_weights; ++n )
	{
	  addweights[n].add    = srch_add[n];
	  addweights[n].weight = weights[n];
	}

      sort_addweights(num_weights, addweights);

      weightlinks[cell_add].addweights = addweights;
      weightlinks[cell_add].nlinks     = num_weights;
    }
}


void weightlinks2remaplinks(long tgt_grid_size, weightlinks_t *weightlinks,  remapvars_t *rv)
{
  long tgt_cell_add;
  long nlinks = 0;

  for ( tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
    {
      if ( weightlinks[tgt_cell_add].nlinks )
	{
	  weightlinks[tgt_cell_add].offset = nlinks;
	  nlinks += weightlinks[tgt_cell_add].nlinks;
	}
    }

  rv->max_links = nlinks;
  rv->num_links = nlinks;
  if ( nlinks )
    {
      rv->src_cell_add = (int*) realloc(rv->src_cell_add, nlinks*sizeof(int));
      rv->tgt_cell_add = (int*) realloc(rv->tgt_cell_add, nlinks*sizeof(int));
      rv->wts          = (double*) realloc(rv->wts, nlinks*sizeof(double));

#if defined(_OPENMP)
#pragma omp parallel for default(shared) \
  shared(rv, weightlinks)		       	 \
  private(tgt_cell_add)
#endif
      for ( tgt_cell_add = 0; tgt_cell_add < tgt_grid_size; ++tgt_cell_add )
	{
	  long num_links = weightlinks[tgt_cell_add].nlinks;
	  long offset    = weightlinks[tgt_cell_add].offset;
	  if ( num_links )
	    {
	      for ( long ilink = 0; ilink < num_links; ++ilink )
		{
		  rv->src_cell_add[offset+ilink] = weightlinks[tgt_cell_add].addweights[ilink].add;
		  rv->tgt_cell_add[offset+ilink] = tgt_cell_add;
		  rv->wts[offset+ilink] = weightlinks[tgt_cell_add].addweights[ilink].weight;
		}
	      free(weightlinks[tgt_cell_add].addweights);
	    }
	}
    }
}
