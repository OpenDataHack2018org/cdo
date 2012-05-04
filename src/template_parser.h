#ifndef TEMPLATE_PARSER_HH
#define TEMPLATE_PARSER_HH

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<locale.h>

#if  defined  (HAVE_LIBXML)
#include<libxml/parser.h>
#include<libxml/tree.h>
#endif

int template_parser( char *Filename, const char *varname );

#endif
