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
#include "config.h" /* HAVE_LIBMAGICS */
#endif

#include "cdo_int.h"
#include "magics_template_parser.h"
#include "StringUtilities.h"

#ifdef HAVE_LIBXML2
#include <libxml/parser.h>
#include <libxml/tree.h>
#endif

#ifdef HAVE_LIBMAGICS
#include "magics_api.h"
#endif

#define DBG 0

#ifdef HAVE_LIBXML2
extern void *magics_node;
#endif

/* Recursive function that sets the Magics parameters from the XML structure */

int
magics_template_parser(void *node)
{
#ifdef HAVE_LIBXML2
  xmlNode *a_node = (xmlNode *) node;
  int param_set_flag;
  xmlNode *cur_node = NULL;
  const char *param_name, *param_type, *param_value;

  if (a_node == NULL) return 0;

#if 0
    fprintf( stdout,"Parsing the magics Node \n");
#endif

  if (!strcmp((const char *) a_node->name, "magics"))
    {
      const char *value = (const char *) xmlGetProp(a_node, (const xmlChar *) "version");

      if (value)
        {
          if (DBG) printf("Version %s \n", value);

          if (atof(value) > 3.0f)
            {
              return 1;
            }
        }
    }

  for (cur_node = a_node->children; cur_node; cur_node = cur_node->next)
    {
      param_name = NULL;
      param_type = NULL;
      param_value = NULL;

      if (cur_node->type == XML_ELEMENT_NODE)
        {

          if (DBG) printf("Node Name: %s \n", cur_node->name);

#if 0
            fprintf( stdout,"Node Name: %s \n", cur_node->name );
#endif

          if (cur_node->properties == NULL)
            {
              if (cur_node->children == NULL)
                {
                  printf("NO ATTRIBUTES!!!\n");
                }
            }
          else
            {

              param_name = (const char *) xmlGetProp(cur_node, (const xmlChar *) "parameter");
              param_type = (const char *) xmlGetProp(cur_node, (const xmlChar *) "type");
              param_value = (const char *) xmlGetProp(cur_node, (const xmlChar *) "value");
#if 0
    		printf( "\t\tAttr name: %s Type: %s Value: %s \n", param_name,param_type,param_value);
#endif

              param_set_flag = SetMagicsParameterValue(param_name, param_type, param_value);

              if (param_set_flag) printf(" Error in Setting the Parameter %s\n", param_name);
            }
        }
    }
#else

  cdoAbort("XML2 support not compiled in!");

#endif

  return 0;
}

#ifdef HAVE_LIBMAGICS
int
SetMagicsParameterValue(const char *param_name, const char *param_type, const char *param_value)
#else
int
SetMagicsParameterValue(const char *, const char *, const char *)
#endif
{
  int ret_flag = 0;
#ifdef HAVE_LIBMAGICS
  int i;
  int split_str_count = 0;
  char **split_str = NULL;
  const char *sep_char = ",";
  const char *search_char = ";";
  double *float_param_list = NULL;

  if (param_name == NULL)
    {
      ret_flag = 1;
      return ret_flag;
    }

  if (param_value == NULL) ret_flag = 2;

  /*   MAGICS++ ENV RELATED PARAMETERS   */
  if (!strcmp(param_type, "environvar"))
    {
      if (!strcmp(param_name, "quiet_option"))
        {
          if (!strcmp(param_value, "off") || !strcmp(param_value, "OFF"))
            {
#if 0
				printf( "Quiet Option %s \n", param_value );
#endif
              if (!unsetenv("MAGPLUS_QUIET"))
                {
                  if (DBG) fprintf(stderr, "Quiet Option %s is un-set successfully!!! \n", param_value);
                }
              else
                fprintf(stderr, "Quiet Option %s COULDN'T be UNSET!!!\n", param_value);
            }

          if (!strcmp(param_value, "on") || !strcmp(param_value, "ON"))
            {
#if 0
				printf( "Quiet Option %s \n", param_value );
#endif
              if (!setenv("MAGPLUS_QUIET", "1", 1))
                {
                  if (DBG) fprintf(stderr, "Quiet Option %s is set successfully!!! \n", param_value);
                }
              else
                fprintf(stderr, "Quiet Option %s COULDN'T be SET!!!\n", param_value);
            }
        }
    }

  /*   MAGICS++ FLOAT TYPE PARAMETERS   */
  else if (!strcmp(param_type, "float"))
    {
      mag_setr(param_name, atof(param_value));
    }

  /*   MAGICS++ FLOAT ARRAY  TYPE    PARAMETERS   */
  else if (!strcmp(param_type, "floatarray"))
    {

#if 0
	        fprintf(stderr, "param_name : %s\tparam_value: %s\n", param_name, param_value);
#endif
      if (strchr(param_value, ';')) sep_char = ";";
      split_str_count = StringSplitWithSeperator(param_value, sep_char, &split_str);
      if (split_str_count)
        {
          float_param_list = (double *) Malloc(sizeof(double) * split_str_count);
          for (i = 0; i < split_str_count; i++)
            {
#if 0
			        fprintf(stderr, "%d %d %s\n", i, split_str_count, split_str[i]);
#endif
              float_param_list[i] = atof(split_str[i]);
            }
          mag_set1r(param_name, float_param_list, split_str_count);
          Free(float_param_list);
          Free(split_str);
        }
    }

  /*   MAGICS++ INT TYPE    PARAMETERS   */
  else if (!strcmp(param_type, "int"))
    {
      mag_seti(param_name, atoi(param_value));
    }

  /*   MAGICS++ INT ARRAY  TYPE    PARAMETERS   */
  else if (!strcmp(param_type, "intarray"))
    {
      if (strchr(param_value, ';')) sep_char = ";";
      split_str_count = StringSplitWithSeperator(param_value, sep_char, &split_str);
      if (split_str_count)
        {
          int *int_param_list = (int *) Malloc(sizeof(int) * split_str_count);
          for (i = 0; i < split_str_count; i++)
            {
              int_param_list[i] = atoi(split_str[i]);
            }
          mag_set1i(param_name, int_param_list, split_str_count);
          Free(int_param_list);
          Free(split_str);
        }
    }

  /*   MAGICS++ STRING TYPE    PARAMETERS   */
  else if (!strcmp(param_type, "string"))
    {
      mag_setc(param_name, param_value);
    }

  /*   MAGICS++ STRINGARRAY  TYPE    PARAMETERS   */
  else if (!strcmp(param_type, "stringarray"))
    {
      if (DBG) fprintf(stderr, "Input strarr is %s  Sep char is %s Search char is %s\n", param_value, sep_char, search_char);
      if (strstr(param_value, ";"))
        {
          sep_char = ";";
        }

      if (DBG) fprintf(stderr, "Input strarr is %s  Sep char is %s\n", param_value, sep_char);
      split_str_count = StringSplitWithSeperator(param_value, sep_char, &split_str);

      if (DBG) fprintf(stderr, "Input strarr is %s split str count is %d Sep char is %s\n", param_value, split_str_count, sep_char);

      mag_set1c(param_name, (const char **) split_str, split_str_count);
      Free(split_str);
    }
  else
    {
      ret_flag = 3;
      fprintf(stderr, "Unknown Parameter Type\n");
    }
#endif

  return ret_flag;
}
