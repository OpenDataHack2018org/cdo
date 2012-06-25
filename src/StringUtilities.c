#include "StringUtilities.h"

int StringSplitWithSeperator( const char *source_string, const char seperator, char*** ptr_split_string )

{
	char *duplicate_src, **temp_list, *temp_str = NULL;
	int n = 0, nn, i;
	int str_len = 0, sep_count = 0;	

	str_len = strlen( source_string );

	if( !str_len )
	  return 0;

	duplicate_src = strdup( source_string );

	for( i = 0; i < str_len; i++ )
	{
		if( duplicate_src[i] == seperator ) 
			sep_count++;
	}	
 
	temp_list  = ( char** )malloc ( sizeof( char* ) * sep_count+1 );

	temp_str = strtok( duplicate_src, &seperator );
	if( temp_str )
		temp_list[0] = temp_str;

	while( temp_str && n < sep_count ) 
	{
		temp_str = NULL;
		temp_str = strtok( NULL, &seperator );
		if( temp_str )
		  temp_list[++n] = temp_str;
 	}

	*ptr_split_string = temp_list;

	return n+1;
}



int IsNumeric (const char *s)
{
    char *ptr;
    if (s == NULL || *s == '\0' || isspace(*s))
      return 0;
    
    strtod (s, &ptr);
    return *ptr == '\0';
}


void StrToUpperCase ( char *sPtr )
{
    while ( *sPtr != '\0' )
    {
      *sPtr = toupper ( ( unsigned char ) *sPtr );
      ++sPtr;
    }
}


void StrToLowerCase ( char *sPtr )
{
    while ( *sPtr != '\0' )
    {
      *sPtr = tolower ( ( unsigned char ) *sPtr );
      ++sPtr;
    }
}
