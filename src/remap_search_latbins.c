#include "cdo.h"
#include "remap.h"


void calc_bin_addr(long gridsize, long nbins, const restr_t* restrict bin_lats, const restr_t* restrict cell_bound_box, int* restrict bin_addr)
{
  long n, n2, nele, nele4;
  restr_t cell_bound_box_lat1, cell_bound_box_lat2;

  for ( n = 0; n < nbins; ++n )
    {
      n2 = n<<1;
      bin_addr[n2  ] = gridsize;
      bin_addr[n2+1] = 0;
    }

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  private(n, n2, nele4, cell_bound_box_lat1, cell_bound_box_lat2)  \
  shared(gridsize, nbins, bin_lats, cell_bound_box, bin_addr)
#endif
  for ( nele = 0; nele < gridsize; ++nele )
    {
      nele4 = nele<<2;
      cell_bound_box_lat1 = cell_bound_box[nele4  ];
      cell_bound_box_lat2 = cell_bound_box[nele4+1];
      for ( n = 0; n < nbins; ++n )
	{
	  n2 = n<<1;
	  if ( cell_bound_box_lat1 <= bin_lats[n2+1] &&
	       cell_bound_box_lat2 >= bin_lats[n2  ] )
	    {
	      /*
#if defined(_OPENMP)
	      if ( nele < bin_addr[n2  ] || nele > bin_addr[n2+1] )
#pragma omp critical
#endif
	      */
		{
		  bin_addr[n2  ] = MIN(nele, bin_addr[n2  ]);
		  bin_addr[n2+1] = MAX(nele, bin_addr[n2+1]);
		}
	    }
	}
    }
}
/*
static
void calc_bin_addr(long gridsize, long nbins, const restr_t* restrict bin_lats, const restr_t* restrict cell_bound_box, int* restrict bin_addr)
{
  long n, n2, nele, nele4;

  for ( n = 0; n < nbins; ++n )
    {
      n2 = n<<1;
      bin_addr[n2  ] = gridsize;
      bin_addr[n2+1] = 0;

#if defined(_OPENMP)
#pragma omp parallel for default(none) \
  private(nele4)	\
  shared(n2, gridsize, bin_lats, cell_bound_box, bin_addr)
#endif
      for ( nele = 0; nele < gridsize; ++nele )
	{
	  nele4 = nele<<2;

	  if ( cell_bound_box[nele4  ] <= bin_lats[n2+1] &&
	       cell_bound_box[nele4+1] >= bin_lats[n2  ] )
	    {
	      bin_addr[n2  ] = MIN(nele, bin_addr[n2  ]);
	      bin_addr[n2+1] = MAX(nele, bin_addr[n2+1]);
	    }
	}
    }
}
*/

void calc_lat_bins(remapgrid_t* src_grid, remapgrid_t* tgt_grid, int map_type)
{
  long nbins;
  long n;      /* Loop counter                  */
  long n2;
  double dlat;                /* lat/lon intervals for search bins  */
  restr_t *bin_lats = NULL;

  nbins = src_grid->num_srch_bins;
  dlat = PI/nbins;

  if ( cdoVerbose ) cdoPrint("Using %d latitude bins to restrict search.", nbins);

  if ( nbins > 0 )
    {
      bin_lats = src_grid->bin_lats = realloc(src_grid->bin_lats, 2*nbins*sizeof(restr_t));

      for ( n = 0; n < nbins; ++n )
	{
	  n2 = n<<1;
	  bin_lats[n2  ] = RESTR_SCALE((n  )*dlat - PIH);
	  bin_lats[n2+1] = RESTR_SCALE((n+1)*dlat - PIH);
	}

      src_grid->bin_addr = realloc(src_grid->bin_addr, 2*nbins*sizeof(int));

      calc_bin_addr(src_grid->size, nbins, bin_lats, src_grid->cell_bound_box, src_grid->bin_addr);

      if ( map_type == MAP_TYPE_CONSERV || map_type == MAP_TYPE_CONSPHERE )
	{
	  tgt_grid->bin_addr = realloc(tgt_grid->bin_addr, 2*nbins*sizeof(int));

	  calc_bin_addr(tgt_grid->size, nbins, bin_lats, tgt_grid->cell_bound_box, tgt_grid->bin_addr);

	  free(src_grid->bin_lats); src_grid->bin_lats = NULL;
	}
   }

  if ( map_type == MAP_TYPE_CONSPHERE )
    {
      free(tgt_grid->cell_bound_box); tgt_grid->cell_bound_box = NULL;
    }
 
  if ( map_type == MAP_TYPE_DISTWGT )
    {
      free(src_grid->cell_bound_box); src_grid->cell_bound_box = NULL;
    }
}


long get_srch_cells(long tgt_grid_add, long nbins, int *bin_addr1, int *bin_addr2,
		    restr_t *tgt_cell_bound_box, restr_t *src_cell_bound_box, long src_grid_size, int *srch_add)
{
  long num_srch_cells;  /* num cells in restricted search arrays   */
  long min_add;         /* addresses for restricting search of     */
  long max_add;         /* destination grid                        */
  long n, n2;           /* generic counters                        */
  long src_grid_add;    /* current linear address for src cell     */
  long src_grid_addm4;
  restr_t bound_box_lat1, bound_box_lat2, bound_box_lon1, bound_box_lon2;

  /* Restrict searches first using search bins */

  min_add = src_grid_size - 1;
  max_add = 0;

  for ( n = 0; n < nbins; ++n )
    {
      n2 = n<<1;
      if ( tgt_grid_add >= bin_addr1[n2] && tgt_grid_add <= bin_addr1[n2+1] )
	{
	  if ( bin_addr2[n2  ] < min_add ) min_add = bin_addr2[n2  ];
	  if ( bin_addr2[n2+1] > max_add ) max_add = bin_addr2[n2+1];
	}
    }

  /* Further restrict searches using bounding boxes */

  bound_box_lat1 = tgt_cell_bound_box[0];
  bound_box_lat2 = tgt_cell_bound_box[1];
  bound_box_lon1 = tgt_cell_bound_box[2];
  bound_box_lon2 = tgt_cell_bound_box[3];

  num_srch_cells = 0;
  for ( src_grid_add = min_add; src_grid_add <= max_add; ++src_grid_add )
    {
      src_grid_addm4 = src_grid_add<<2;
      if ( (src_cell_bound_box[src_grid_addm4+2] <= bound_box_lon2)  &&
	   (src_cell_bound_box[src_grid_addm4+3] >= bound_box_lon1) )
	{
	  if ( (src_cell_bound_box[src_grid_addm4  ] <= bound_box_lat2)  &&
	       (src_cell_bound_box[src_grid_addm4+1] >= bound_box_lat1) )
	    {
	      srch_add[num_srch_cells] = src_grid_add;
	      num_srch_cells++;
	    }
	}
    }

  if ( bound_box_lon1 < RESTR_SCALE(0.) || bound_box_lon2 > RESTR_SCALE(PI2) )
    {
      if ( bound_box_lon1 < RESTR_SCALE(0.) )
	{
	  bound_box_lon1 += RESTR_SCALE(PI2);
	  bound_box_lon2 += RESTR_SCALE(PI2);
	}
      else
	{
	  bound_box_lon1 -= RESTR_SCALE(PI2);
	  bound_box_lon2 -= RESTR_SCALE(PI2);
	}

      for ( src_grid_add = min_add; src_grid_add <= max_add; ++src_grid_add )
	{
	  src_grid_addm4 = src_grid_add<<2;
	  if ( (src_cell_bound_box[src_grid_addm4+2] <= bound_box_lon2)  &&
	       (src_cell_bound_box[src_grid_addm4+3] >= bound_box_lon1) )
	    {
	      if ( (src_cell_bound_box[src_grid_addm4  ] <= bound_box_lat2)  &&
		   (src_cell_bound_box[src_grid_addm4+1] >= bound_box_lat1) )
		{
		  long ii;
		  for ( ii = 0; ii < num_srch_cells; ++ii )
		    if ( srch_add[ii] == src_grid_add ) break;
		  
		  if ( ii == num_srch_cells )
		    {
		      srch_add[num_srch_cells] = src_grid_add;
		      num_srch_cells++;
		    }
		}
	    }
	}
    }

  return (num_srch_cells);
}
