/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2018 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cdo_int.h"
#include "template_parser.h"
#include "magics_template_parser.h"
#include "results_template_parser.h"

#ifdef HAVE_LIBXML2
#include <libxml/parser.h>
#include <libxml/tree.h>
#endif

#define DBG_MSG 0

/* extern int GetMagicsParameterInfo( const char *user_name, char **magics_name,
 * char **magics_type ); */

#ifdef HAVE_LIBXML2
extern int GetMagicsParameterInfo(const char *user_name, char *param_value);
#endif

/* Recursive function that sets the results parameters from the XML structure */

int
results_template_parser(void *node, const char *varname)
{
#ifdef HAVE_LIBXML2
  xmlNode *a_node = (xmlNode *) node;
  xmlNode *cur_node = NULL;
  xmlAttrPtr attr = NULL;
  xmlChar *param_value;

  if (a_node == NULL) return 1;

  if (!strcmp((const char *) a_node->name, "results"))
    {
      const char *value = (const char *) xmlGetProp(a_node, (const xmlChar *) "version");

      if (value)
        {
          if (DBG_MSG) printf("Version %s \n", value);

          if (atof(value) > 3.0f)
            {
              return 1;
            }
        }
    }

  for (cur_node = a_node->children; cur_node; cur_node = cur_node->next)
    {
#if 0
	xmlChar *param_name = NULL;
	char *param_type = NULL;
#endif
      param_value = NULL;

      if (cur_node->type == XML_ELEMENT_NODE)
        {

          if (DBG_MSG) printf("Node Name: %s \n", cur_node->name);

          if (cur_node->properties == NULL)
            {
              if (cur_node->children == NULL)
                {
                  printf("NO ATTRIBUTES!!!\n");
                }
            }
          else
            {

/* 	Loop Over the attributes and get the corresponding
        Magics Parameter name and type, set the value
*/

#if 0
	      printf( "Finding varname = %s  result_name = %s\n", varname, xmlGetProp( cur_node,"name") );
#endif

              if (strcmp(varname, (const char *) xmlGetProp(cur_node, (xmlChar *) "name")) == 0)
                {
#if 0
	          printf( "Found varname = %s  result_name = %s\n", varname, xmlGetProp( cur_node,"name") );
#endif

                  for (attr = cur_node->properties; attr; attr = attr->next)
                    {
                      if (attr != NULL)
                        {

                          param_value = xmlNodeGetContent(attr->children);

                          /* if( !GetMagicsParameterInfo( attr->name,
                           * &magics_param_name, &param_type ) ) */
                          if (!GetMagicsParameterInfo((const char *) attr->name, (char *) param_value))
                            {

#if 0
			      printf("Done corresponding Magics Parameter found!\n");
			      printf("Setting corresponding Magics Parameter %s and type %s!\n",magics_param_name, param_type );
	  		      fprintf(stderr, "param_value: %s\n", param_value);
			      SetMagicsParameterValue( magics_param_name, param_type, param_value );
#endif
                            }
                          else
                            {
#if 0
			      printf("No corresponding Magics Parameter found!\n");
#endif
                            }
                        }
                    }

                  break;
                }
              else
                {
                  fprintf(stderr, "Var Name not matching resetting Magics Params!\n");
                  /* Call the Reset functions of all the features to Reset the
                   * magics params to default */
                }
            }
        }
    }
#else

  cdoAbort("XML2 support not compiled in!");

#endif

  return 0;
}
