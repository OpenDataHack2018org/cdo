#include "template_parser.h"
#include "magics_template_parser.h"
#include "results_template_parser.h"

#define DBG_MSG 0 


extern int magics_template_parser();
extern int results_template_parser();



int template_parser(  char *Filename, const char *varname )

{
#if  defined  (HAVE_LIBXML)
        xmlDoc         *doc = NULL;
        xmlNode        *root_element = NULL;

	printf("XML File Name:\t%s\n",Filename);
        doc = xmlReadFile( Filename, NULL, 0 );
        if ( doc == NULL )
        {
                  printf( "Error: Could not parse the file \"%s\"\n", Filename );
        	  return (1);
        }
        else
        {
                  /* 
		     Get the name of the root element node 
		     If "magics" , call "magics" parser
		     If "results", call "results" parser
		  */                      

                  root_element = xmlDocGetRootElement( doc );
		  
		  if( !strcmp( root_element->name, "magics" ) )
		  {
			printf( "Found Magics Template! \n" );
                  	if ( magics_template_parser( root_element ) == 1 )
		  	{
				printf( "Un-Supported version of Magics++! \n" );
	        		return (2);
			}
		  }
		  else if( !strcmp( root_element->name, "results" ) )
		  {
			printf( "Found Results Template! \n" );
                  	results_template_parser( root_element, varname );
			 /* Needs some error handling */
		  }

 
                  /*** free the document ***/
                  xmlFreeDoc( doc );
        }

        /*** Free the global variables that may
         *   have been allocated by the parser. 
        ***/

        xmlCleanupParser();

        return (0);
#else
	fprintf(stderr, "XML support not compiled in!");
	return (1);
#endif
}
