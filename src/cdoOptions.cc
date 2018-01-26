
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
