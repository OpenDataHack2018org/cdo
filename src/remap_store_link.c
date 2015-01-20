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
int cmp_adds4(const void *s1, const void *s2)
{
  int cmp = 0;
  const addweight4_t* c1 = (const addweight4_t*) s1;
  const addweight4_t* c2 = (const addweight4_t*) s2;

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

static
void sort_addweights4(size_t num_weights, addweight4_t *addweights)
{
  size_t n;

  for ( n = 1; n < num_weights; ++n )
    if ( addweights[n].add < addweights[n-1].add ) break;
  if ( n == num_weights ) return;

  qsort(addweights, num_weights, sizeof(addweight4_t), cmp_adds);
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


void sort_add_and_wgts4(size_t num_weights, int *src_add, double *wgts[4])
{
  size_t n;

  for ( n = 1; n < num_weights; ++n )
    if ( src_add[n] < src_add[n-1] ) break;
  if ( n == num_weights ) return;

  addweight4_t addweights[num_weights];

  for ( n = 0; n < num_weights; ++n )
    {
      addweights[n].add       = src_add[n];
      addweights[n].weight[0] = wgts[0][n];
      addweights[n].weight[1] = wgts[1][n];
      addweights[n].weight[2] = wgts[2][n];
      addweights[n].weight[3] = wgts[3][n];
    }

  qsort(addweights, num_weights, sizeof(addweight4_t), cmp_adds4);

  for ( n = 0; n < num_weights; ++n )
    {
      src_add[n] = addweights[n].add;
      wgts[0][n] = addweights[n].weight[0];
      wgts[1][n] = addweights[n].weight[1];
      wgts[2][n] = addweights[n].weight[2];
      wgts[3][n] = addweights[n].weight[3];
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


void store_weightlinks4(long num_weights, int *srch_add, double *weights[4], long cell_add, weightlinks4_t *weightlinks)
{
  weightlinks[cell_add].nlinks = 0;
  weightlinks[cell_add].offset = 0;

  if ( num_weights )
    {
      addweight4_t *addweights = (addweight4_t *) malloc(num_weights*sizeof(addweight4_t));
      for ( long n = 0; n < num_weights; ++n )
	{
	  addweights[n].add       = srch_add[n];
	  addweights[n].weight[0] = weights[0][n];
	  addweights[n].weight[1] = weights[1][n];
	  addweights[n].weight[2] = weights[2][n];
	  addweights[n].weight[3] = weights[3][n];
	}

      sort_addweights4(num_weights, addweights);

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


void weightlinks2remaplinks4(long tgt_grid_size, weightlinks4_t *weightlinks,  remapvars_t *rv)
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
      rv->wts          = (double*) realloc(rv->wts, 4*nlinks*sizeof(double));

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
		  rv->wts[(offset+ilink)*4+0] = weightlinks[tgt_cell_add].addweights[ilink].weight[0];
		  rv->wts[(offset+ilink)*4+1] = weightlinks[tgt_cell_add].addweights[ilink].weight[1];
		  rv->wts[(offset+ilink)*4+2] = weightlinks[tgt_cell_add].addweights[ilink].weight[2];
		  rv->wts[(offset+ilink)*4+3] = weightlinks[tgt_cell_add].addweights[ilink].weight[3];
		}
	      free(weightlinks[tgt_cell_add].addweights);
	    }
	}
    }
}
