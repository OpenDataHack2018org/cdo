#include <stdio.h>
#include <stdlib.h>

int get_alignment(double *ptr)
{
  int64_t iptr = (int64_t) ptr;
  int mk[4] = {64, 32, 16, 8};
  int malign = 0;

  for ( int i = 0; i < 4; ++i )
    if ( iptr%mk[i] == 0 )
      {
	malign = mk[4];
	break;
      }
  
  return malign;
}

int main(void)
{
  int nelem = 97;

  double *array = (double *) malloc(nelem*sizeof(double));

  int64_t iptr = (int64_t) array;
  int malign = get_alignment(array);
  printf("mem alignment: %d\n", malign);

  free(array);
  
  return 0;
}
