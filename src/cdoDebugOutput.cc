#include "cdoDebugOutput.h"



namespace CdoDebug{

#ifdef DEBUG
   int PSTREAM = 1;
#else
   int PSTREAM = 0;
#endif
   bool PROCESS = false;
   std::string outfile = "";
   bool print_to_seperate_file = false;
}

