#ifndef MAGICS_TEMPLATE_PARSER_HH
#define MAGICS_TEMPLATE_PARSER_HH
#endif

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<locale.h>

#include<libxml/parser.h>
#include<libxml/tree.h>

#include "cdo_int.h"

int magics_template_parser( xmlNode * a_node );

int SetMagicsParameterValue( const char *param_name, const char *param_type, const char *param_value );
