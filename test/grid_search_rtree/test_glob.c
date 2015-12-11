
#include "index.h"
#include <stdio.h>
#include <time.h>
#include <math.h>
static double const rad          = M_PI / 180.0;

struct Rect rects[] = {
	{0, 0, 2, 2}, // xmin, ymin, xmax, ymax (for 2 dimensional RTree)
	{5, 5, 7, 7},
	{8, 5, 9, 6},
	{7, 1, 9, 2},
};
int nrects = sizeof(rects) / sizeof(rects[0]);
struct Rect search_rect = {
	{6, 4, 10, 6}, // search will find above rects that this one overlaps
};

int MySearchCallback(int id, void* arg) 
{
	// Note: -1 to make up for the +1 when data was inserted
	printf("Hit data rect %d\n", id-1);
	return 1; // keep going
}

int MySearchCallback2(int id, void* arg) 
{
	// Note: -1 to make up for the +1 when data was inserted
  //	printf("Hit data rect %d\n", id-1);
	return 1; // keep going
}

void main()
{
	struct Node* root = RTreeNewIndex();
	int i, nhits;
	printf("nrects = %d\n", nrects);
	/*
	 * Insert all the data rects.
	 * Notes about the arguments:
	 * parameter 1 is the rect being inserted,
	 * parameter 2 is its ID. NOTE: *** ID MUST NEVER BE ZERO ***, hence the +1,
	 * parameter 3 is the root of the tree. Note: its address is passed
	 * because it can change as a result of this call, therefore no other parts
	 * of this code should stash its address since it could change undernieth.
	 * parameter 4 is always zero which means to add from the root.
	 */
	for(i=0; i<nrects; i++)
		RTreeInsertRect(&rects[i], i+1, &root, 0); // i+1 is rect ID. Note: root can change
	nhits = RTreeSearch(root, &search_rect, MySearchCallback, 0);
	printf("Search resulted in %d hits\n", nhits);

      printf("test cell search for cyclic grids\n");
      //---------------
      // setup
      //---------------
      const int nxa = 720*4;
      const int nya = 360*4+1;
      const int nxb = 360*4;
      const int nyb = 180*4+1;
      /*
      const int nxa = 72;
      const int nya = 37;
      const int nxb = 36;
      const int nyb = 19;
      */
      double coords_x_a[nxa];
      double coords_y_a[nya];

      for ( unsigned i = 0; i < nxa; ++i ) coords_x_a[i] = -180 + i*360./nxa;
      for ( unsigned i = 0; i < nya; ++i ) coords_y_a[i] =  -90 + i*360./nxa;

      printf("nxa: %d xa: %g %g ... %g %g\n", nxa, coords_x_a[0], coords_x_a[1], coords_x_a[nxa-2], coords_x_a[nxa-1]);
      printf("nya: %d ya: %g %g ... %g %g\n", nya, coords_y_a[0], coords_y_a[1], coords_y_a[nya-2], coords_y_a[nya-1]);

      for ( unsigned i = 0; i < nxa; ++i ) coords_x_a[i] *= rad;
      for ( unsigned i = 0; i < nya; ++i ) coords_y_a[i] *= rad;   

      double coords_x_b[nxb];
      double coords_y_b[nyb];

      for ( unsigned i = 0; i < nxb; ++i ) coords_x_b[i] = -180 + i*360./nxb;
      for ( unsigned i = 0; i < nyb; ++i ) coords_y_b[i] =  -90 + i*360./nxb;

      printf("nxb: %d xb: %g %g ... %g %g\n", nxb, coords_x_b[0], coords_x_b[1], coords_x_b[nxb-2], coords_x_b[nxb-1]);
      printf("nyb: %d yb: %g %g ... %g %g\n", nyb, coords_y_b[0], coords_y_b[1], coords_y_b[nyb-2], coords_y_b[nyb-1]);

      for ( unsigned i = 0; i < nxb; ++i ) coords_x_b[i] *= rad;
      for ( unsigned i = 0; i < nyb; ++i ) coords_y_b[i] *= rad;
   
      //---------------
      // testing
      //---------------

      clock_t start, finish;

      struct Node* xroot = RTreeNewIndex();
         
      start = clock();

      struct Rect rect;
      for ( unsigned j = 0; j < nya-1; ++j )
        for ( unsigned i = 0; i < nxa-1; ++i )
          {
            rect.boundary[0] = coords_x_a[i];
            rect.boundary[1] = coords_y_a[j];
            rect.boundary[2] = coords_x_a[i+1];
            rect.boundary[3] = coords_y_a[j+1];
            //printf("%u %g %g %g %g\n", j*nxa+i+1, rect.boundary[0], rect.boundary[1], rect.boundary[2], rect.boundary[3]);
            RTreeInsertRect(&rect, j*nxa+i+1, &xroot, 0);
          }

      finish = clock();
      printf("create: %.2f seconds\n", ((double)(finish-start))/CLOCKS_PER_SEC);
      start = clock();

      struct Rect search_rect;
      for ( unsigned j = 0; j < nyb-1; ++j )
        for ( unsigned i = 0; i < nxb-1; ++i )
          {
            search_rect.boundary[0] = coords_x_b[i];
            search_rect.boundary[1] = coords_y_b[j];
            search_rect.boundary[2] = coords_x_b[i+1];
            search_rect.boundary[3] = coords_y_b[j+1];
            // printf("index %u\n", j*nxb+i+1);
            nhits = RTreeSearch(xroot, &search_rect, MySearchCallback2, 0);
          }

      finish = clock();
      printf("search: %.2f seconds\n", ((double)(finish-start))/CLOCKS_PER_SEC);
}
