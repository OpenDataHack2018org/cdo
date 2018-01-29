#include <cdi.h>

#include "util.h"

void nospec(int vlistID)
{
  int nvars = vlistNvars(vlistID);
  for ( int varID = 0; varID < nvars; varID++ )
    {
      int gridID = vlistInqVarGrid(vlistID, varID);
      int gridtype = gridInqType(gridID);
      if ( gridtype == GRID_SPECTRAL )
	cdoAbort("Operator not defined for spectral fields");
    }
}
