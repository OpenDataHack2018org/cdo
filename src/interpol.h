#ifndef INTERPOL_H
#define INTERPOL_H

#include <stdio.h>

void intgridbil(field_type *field1, field_type *field2);
void intgridcon(field_type *field1, field_type *field2);
void interpolate(field_type *field1, field_type *field2);

void intgriddis(field_type *field1, field_type *field2, size_t nnn);

#endif
