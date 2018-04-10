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
#ifndef CDOOPTIONS_H
#define CDOOPTIONS_H

namespace Cdo
{
extern const char *progname;
}

namespace Options
{
extern bool benchmark;
extern bool silentMode;

extern bool cdoCompress;
extern int cdoCompType;
extern int cdoCompLevel;
extern bool cdoInteractive;
}  // namespace Options

namespace Threading
{
extern int ompNumThreads;
extern bool cdoLockIO;
}  // namespace Threading

#endif
