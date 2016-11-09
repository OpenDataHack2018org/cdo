/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2016 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

/*
   This module contains the following operators:

*/

#include "cdo.h"
#include "cdo_int.h"


void *Nmldump(void *argument)
{
  cdoInitialize(argument);

  list_t *pmlist = namelist_to_pmlist(stdin, "STDIN");

  list_for_each(pmlist, pmlist_print_iter);

  list_destroy(pmlist);

  cdoFinish();

  return 0;
}
