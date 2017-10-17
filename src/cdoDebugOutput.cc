#include "cdoDebugOutput.h"



namespace CdoDebug{

#if defined(DEBUG-PSTREAM) || defined(DEBUG)
   int PSTREAM = 1;
#else
   int PSTREAM = 0;
#endif

#if defined(DEBUG-PROCESS) || defined(DEBUG)
   bool PROCESS = true;
#else 
   bool PROCESS = false;
#endif

   std::string outfile = "";
   bool print_to_seperate_file = false;
}

