#ifndef CDOOPTIONS_H
#define CDOOPTIONS_H

/*TEMP*/ // replace with constexpr
#define  CDO_EXP_LOCAL   1
#define  CDO_EXP_REMOTE  2

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
}

namespace Threading
{
    extern int ompNumThreads;
    extern bool cdoLockIO;
}



#endif
