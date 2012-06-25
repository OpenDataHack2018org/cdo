#ifndef STR_UTILITIES_H
#define STR_UTILITIES_H

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include <ctype.h>

int StringSplitWithSeperator( const char *source_string, const char seperator, char*** ptr_split_string );

int IsNumeric (const char *s);

void StrToUpperCase ( char *sPtr );

void StrToLowerCase ( char *sPtr );

#endif
