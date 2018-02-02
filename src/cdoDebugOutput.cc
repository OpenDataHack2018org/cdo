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
#include "cdoDebugOutput.h"

namespace CdoDebug
{
    void SetDebug(int p_debug_level)
    {
        /*
        p_debug_level   0: off
        p_debug_level   1: on
        p_debug_level   2: cdi
        p_debug_level   4: memory
        p_debug_level   8: file
        p_debug_level  16: format
        p_debug_level  32: cdo
        p_debug_level  64: stream
        p_debug_level 128: pipe
        p_debug_level 256: pthread
        p_debug_level 512: process
        */

        if ( p_debug_level == 1 || (p_debug_level &  32) ) cdoDebug = 1;
        if ( p_debug_level == 1 || (p_debug_level &  64) ) PSTREAM = 1;
        if ( p_debug_level == 1 || (p_debug_level &  512) ) PROCESS = 1;
        #ifdef  HAVE_LIBPTHREAD
        if ( p_debug_level == 1 || (p_debug_level & 128) ) PIPE = 1;
        if ( p_debug_level == 1 || (p_debug_level & 256) ) PTHREAD = 1;
        #endif
    }

    //Debug Switches
     int cdoDebug;
     int cdoDebugExt = 0;     //  Debug level for the KNMI extensions
    //Subsystem Debug Switches
     int  PSTREAM = 0;
     bool PROCESS = 0;
     bool PIPE = false;
     int ARGUMENT = 0;
     int PTHREAD = 0;

    //File switches and streams
     std::string outfile;
     bool print_to_seperate_file;
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
      outfile_stream.open(outfile, std::fstream::in | std::fstream::app);

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

    std::string argvToString(int argc, const char ** argv)
    {
        std::string input_string = "";
        for (int i = 0; i < argc; i++)
        {
            input_string += argv[i];
            input_string += " ";
        }
        return input_string;

    }

        void printMessage(std::stringstream &p_message, bool both )
        {
            if(print_to_seperate_file || (print_to_seperate_file && both))
            {
                outfile_stream <<  p_message.str();
            }

            if(!print_to_seperate_file || both)
            {
                std::cout << p_message.str();
            }
        }
}
namespace CdoLog
{
    void StdOut(std::stringstream & p_message)
    {
        std::cout << p_message.str();
    }
}
