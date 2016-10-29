#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>

enum literal_type {E_LTYPE_NONE = 1, E_LTYPE_BYTE, E_LTYPE_SHORT, E_LTYPE_INT, E_LTYPE_FLOAT, E_LTYPE_DOUBLE};

// CHECK RANGE!!!
int literal_type(const char *literal)
{
  if ( literal && *literal )
    {
      char *endptr;
      errno = 0;
      long lval = strtol(literal, &endptr, 10);
      if ( errno == 0 && *endptr == 0 ) return E_LTYPE_INT;
      else if ( errno == 0 && *(endptr+1) == 0 )
        {
          if      ( *endptr == 's' && ( lval >= SHRT_MIN && lval <= SHRT_MAX) )
            return E_LTYPE_SHORT;
          else if ( *endptr == 'b' && ( lval >= SCHAR_MIN && lval <= SCHAR_MAX) )
            return E_LTYPE_BYTE;
        }
      else
        {
          errno = 0;
          float fval = strtof(literal, &endptr);
          (void) fval;
          if ( errno == 0 && (*endptr == 0 || (*(endptr+1) == 0 && *endptr == 'f')) )
            return E_LTYPE_FLOAT;
          else
            {
              double dval = strtod(literal, &endptr);
              (void)dval;
              if ( *endptr == 0 ) return E_LTYPE_DOUBLE;
            }
        }
    }

  return E_LTYPE_NONE;
}

int main(void)
{
  const char *literals[] = {"127b", "-32768s", "-2147483647", "-1.e+36f", "1.e+308"};
  int nliterals = sizeof(literals) / sizeof(literals[0]);
  
  for ( int i = 0; i < nliterals; ++i )
    printf("%d %s type = %d\n", i+1, literals[i], literal_type(literals[i]));

  return 0;
}
