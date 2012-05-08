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
	{ "clevs","contour_level_list","floatarray"},
	{ "ccols","mag_2","intarray"},
	{ "color_table","mag_3","intarray"}

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
       static int once = 1;
       int ret_flag = 0;
       CdoMagicsMapper target, *result;
       target.cdo_name = user_name; 

       printf ("Finding  %s.\n", user_name);
       if( once )
       {
       		qsort ( mapper, PARAM_COUNT, sizeof ( CdoMagicsMapper ), ( void * )Compare );
		once = 0;
       }

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
