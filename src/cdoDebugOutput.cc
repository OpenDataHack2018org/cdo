#include "cdoDebugOutput.h"

namespace CdoLog
{
    void StdOut(std::stringstream & p_message)
    {
        std::cout << p_message.str();
    }
}

namespace CdoDebug
{
    int PTHREAD;
    int PSTREAM;
    bool PROCESS;
    bool PIPE;
    int ARGUMENT;

    std::string outfile = "";
    bool print_to_seperate_file = false;
    std::fstream outfile_stream;

    std::string
    get_padding(const char *p_func)
    {
      size_t len = strlen(p_func);

      return std::string(30 - len, ' ');
    }

    void
    CdoStartMessage()
    {
      std::stringstream message;
      outfile_stream = std::fstream(outfile, std::fstream::in | std::fstream::app);

      message << std::string(30, ' ') << "  == CDO Start ==" << std::endl;
      printMessage(message);
    }
    void
    CdoEndMessage()
    {
      std::stringstream message;
      message << std::string(30, ' ') << "  == CDO End ==" << std::endl;
      printMessage(message);
      outfile_stream.close();
    }
}
