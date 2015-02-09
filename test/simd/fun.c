void fun1(int nelem, double *restrict array1, double *restrict array2)
{
  for ( int i = 0; i < nelem; ++i )
    array1[i] += array2[i];
}
