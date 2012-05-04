#include "CdoMagicsMapper.h"

#define PARAM_COUNT  sizeof( mapper ) / sizeof ( CdoMagicsMapper )


/* Define an array of Mapper structures to sort. */

typedef struct 

{
	char *cdo_name;
	char *magics_name;
	char *magics_type;

} CdoMagicsMapper;


CdoMagicsMapper mapper[] =

{
	{ "clevs","mag_cont_levels","floatarray"},
	{ "ccols","mag_colours","intarray"},
	{ "color_table","mag_colour_tabe","intarray"}

};

int Compare( CdoMagicsMapper  *Parameter_one , CdoMagicsMapper *Parameter_two );

void PrintResult ( const CdoMagicsMapper *c );

/* This is the comparison function used for sorting and searching. */
     
int Compare( CdoMagicsMapper *p1, CdoMagicsMapper *p2 )

{
	return strcmp ( p1->cdo_name, p2->cdo_name );
}
     
     
/* Print information about a critter. */
     
void PrintResult ( const CdoMagicsMapper *c )

{
	printf ( "CDO Name:%s\t MAGICS Name:%s\t MAGICS Type:%s\n", c->cdo_name, c->magics_name, c->magics_type );
}
     
     
/* Do the lookup into the sorted array. */

int GetMagicsParameterInfo( const char *user_name, char **magics_name, char **magics_type )

{
       int ret_flag = 0;
       CdoMagicsMapper target, *result;
       target.cdo_name = user_name; 
       printf ("Finding  %s.\n", user_name);

       printf ("Finding  %s.\n", user_name);
       result = bsearch ( &target, mapper, PARAM_COUNT, sizeof ( CdoMagicsMapper ),
                          ( void * )Compare );
       if ( result )
       {
         	printf ("Found %s.\n", user_name);
         	printf ("Name %s\n", result->magics_name);
		*magics_name = result->magics_name;
		*magics_type = result->magics_type;
       }
       else
       {
         ret_flag = 1;
       }
       return ret_flag;
}
     

     /* Main program. */
  
/*   
int main (void)

{
       int i;
       char *user_name,*param_type,*param_name;
       char *user_name1,*param_type1,*param_name1;
	
     
       for ( i = 0; i < PARAM_COUNT; i++ )
		PrintResult ( &mapper[i] );
       printf ( "\n" );
     
       qsort ( mapper, PARAM_COUNT, sizeof ( CdoMagicsMapper ), ( void * )Compare );
     
       for ( i = 0; i < PARAM_COUNT; i++ )
		PrintResult ( &mapper[i] );
       printf ("\n");

       user_name  = "clevs";

       printf ( "Find UserName %s\n", user_name );
     
       if( !GetMagicsParameterInfo ( user_name, &param_name, &param_type  ) )
       {
	       printf ( "RESULT\t MAGICS Name:%s\t MAGICS Type:%s\n", param_name, param_type );
       }
       else
         printf ("Couldn't find %s.\n", user_name);


       user_name = "ccols";
       printf ( "Find UserName %s\n", user_name );

       if( !GetMagicsParameterInfo ( user_name, &param_name, &param_type  ) )
       {
	       printf ( "RESULT\t MAGICS Name:%s\t MAGICS Type:%s\n", param_name, param_type );
       }
       else
         printf ("Couldn't find %s.\n", user_name);

       for ( i = 0; i < PARAM_COUNT; i++ )
		PrintResult ( &mapper[i] );
       printf ("\n");

       return 0;
}
*/
