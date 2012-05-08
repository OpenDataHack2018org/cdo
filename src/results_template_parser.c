#include "template_parser.h"
#include "magics_template_parser.h"
#include "results_template_parser.h"

#define DBG_MSG 1 


extern int GetMagicsParameterInfo( const char *user_name, char **magics_name, char **magics_type );

#if  defined  (HAVE_LIBXML)

extern xmlNode *results_node;


/* Recursive function that sets the results parameters from the XML structure */

int results_template_parser( xmlNode * a_node, const char *varname ) 

{
    xmlNode *cur_node = NULL;
    xmlAttrPtr attr = NULL;
    xmlChar    *param_name,*param_value,*value;
    char       *magics_param_name,*param_type;

    for ( cur_node = a_node; cur_node; cur_node = cur_node->next )
    {
	param_name = NULL;
	param_type = NULL;
	param_value = NULL;

        if ( cur_node->type == XML_ELEMENT_NODE )
        {
		
	    if( DBG_MSG )
            	printf( "Node Name: %s \n", cur_node->name );

	    if( !strcmp( cur_node->name, "results" ) )
	    {
	 	value = xmlGetProp( cur_node, "version" );

		if( value )
		{
	    		if( DBG_MSG )
				printf( "Version %s \n", value ); 

			if( atof( value ) > 3.0f ) 
			{
				return 1;
			}
		}
        	results_template_parser( cur_node->children, varname );
		continue;
	    }
	
	    if( cur_node->properties == NULL )
	    {
		if( cur_node->children == NULL )
		{
			printf( "NO ATTRIBUTES!!!\n" );
		}
		else 
		{
		  results_template_parser( cur_node->children, varname );
		}
	    }
	    else
	    {
		
		/* 	Loop Over the attributes and get the corresponding
                   	Magics Parameter name and type, set the value 
		*/

	      printf("varname = %s  result_name = %s\n", varname, xmlGetProp( cur_node,"name"));

	      if ( strcmp(varname, xmlGetProp( cur_node,"name")) == 0 )
	      {
		  for( attr = cur_node->properties; attr; attr = attr->next )
		  {	
		      if( attr != NULL )
		      {
			  printf( "Attr name: %s Conte: %s \n", attr->name, xmlNodeGetContent( attr->children ) );

			  if( !GetMagicsParameterInfo( attr->name, &magics_param_name, &param_type ) )
			  {
			      printf("Setting corresponding Magics Parameter %s and type %s!\n",magics_param_name, param_type );
			      param_value = xmlNodeGetContent( attr->children );
	  fprintf(stderr, "param_value: %s\n", param_value);
			      SetMagicsParameterValue( magics_param_name, param_type, param_value );
#if 0
#endif
			  }
			  else
			  {
			      printf("No corresponding Magics Parameter found!\n");
			  }
		      }	
		  }	
	      }
	      else
	      {
	   	  printf("Var Name not matching!\n");
	      }
	    }
        }
    }
    return 0;
}
#endif
