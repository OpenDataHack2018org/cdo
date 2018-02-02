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

#include <cdi.h>
#include "cdoOptions.h"

namespace Cdo
{
    const char *progname;
}
namespace Options
{
    bool benchmark = false;
    bool silentMode = false;

    bool cdoCompress = false;
    int cdoCompType = CDI_COMPRESS_NONE;
    int cdoCompLevel= 0;
}

namespace Threading
{
    int ompNumThreads = 1;
    bool cdoLockIO = false;
}
