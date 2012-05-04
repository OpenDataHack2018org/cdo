#ifndef RESULTS_TEMPLATE_PARSER_HH
#define RESULTS_TEMPLATE_PARSER_HH

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#if  defined  (HAVE_LIBXML)
#include<libxml/parser.h>
#include<libxml/tree.h>
#endif

#if  defined  (HAVE_LIBXML)
int results_template_parser( xmlNode * a_node, const char *varname ); 
#endif
/* int GetMagicsParameterInfo( const char *user_name, char *magics_name, char *magics_type ); */


#endif
