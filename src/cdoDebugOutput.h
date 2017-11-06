#ifndef DEBUG_SWITCHES_H
#define DEBUG_SWITCHES_H
#include <iostream>
#include <sstream>

#include <fstream>
#include <string.h>

namespace CdoLog
{
    void StdOut(std::stringstream &message);

    template <typename ...T>
    void  expand(std::stringstream &p_message, T&& ...args)
    {
           //for showing that the dummy array is never used
           using expander = int[];
           //creating dummy array with expanding the parameter pack in its
           //initializer list while writing each element of the pack into message
           expander{0,(void(p_message << std::forward<T>(args)),0)...};
           p_message << std::endl;
    }

    template <typename ...T>
    void StdOut(T&& ...args)
    {
        std::stringstream message;
        expand(message, args...);
        std::cout << message.str();
    }
 
    template <typename ...T>
    void StdErr(T&& ...args)
    {
        std::stringstream message;
        expand(message, args...);
        std::cout << message.str();
    }
}

namespace CdoDebug
{
    //Debug Switches
    extern int cdoDebug;
    extern int cdoDebugExt; //  Debug level for the KNMI extensions
    //Subsystem Debug Switches
    extern int  PSTREAM;
    extern bool PROCESS;
    extern bool PIPE;
    extern int ARGUMENT;
    extern int PTHREAD;

    //File switches and streams
    extern std::string outfile;
    extern bool print_to_seperate_file;
    extern std::fstream outfile_stream;


    std::string get_padding(const char *p_func);

    void CdoStartMessage();
    void CdoEndMessage();
    void SetDebug(int p_debug_level);

   namespace{
        void printMessage(std::stringstream &p_message)
        {
            if(!print_to_seperate_file)
            {
                std::cout << p_message.str();
            }
            else 
            {
                outfile_stream <<  p_message.str();
            }
        }
    }


   
    template <typename ...T>
    void Message_ (const char * p_func, T&& ...args)
    {
        std::stringstream message;
        message << p_func <<": " << get_padding(p_func);
        CdoLog::expand(message, args...);
        printMessage(message);
    }

    
    template <typename ...T>
    void Warning_(T&& ...args)
    {
        std::stringstream message;
        message << "Warning: ";
        CdoLog::expand(message, args...);
        std::cout << message.str();
    }


}

namespace CdoError{
    static int _ExitOnError = 1;

    template <typename ...T>
    void Error_(const char* p_file, const int p_line, const char* caller, T&& ...args)
    {
          std::stringstream message;
          message << "Error in: " << p_file << ":" << p_line << " ";
          CdoLog::expand(message, args...);
          CdoLog::StdOut(message);
          if ( CdoError::_ExitOnError )
          {
              exit(EXIT_FAILURE);
          }
    }
    
    template <typename ...T>
    void SysError_(const char* p_file, const int p_line,  const char* p_func, T&& ...args)
    {
        int saved_errno = errno;
        std::stringstream message;
        message << "SysError in:" << p_file << std::endl;
        message << "    " << "in function: p_func ,line: " << p_line << std::endl;
        CdoLog::StdOut(message, args...);
        CdoLog::StdOut(message);
        if(saved_errno)
        {
            errno = saved_errno;
            perror("Sytem error message");
        }
        exit(EXIT_FAILURE);
    }
}
#if  defined  WITH_CALLER_NAME
#  define  SYS_ERROR(...)    CdoError::SysError_( __FILE__ , __LINE__ , __func__ , __VA_ARGS__)
#  define    ERROR_C(...)    CdoError::Error_( __FILE__ , __LINE__ ,  caller, __VA_ARGS__)
#  define    ERROR(...)    CdoError::Error_( __FILE__ , __LINE__ , __func__ , __VA_ARGS__)
#  define   WARNING(...)    CdoError::Warning_( __func__ , __VA_ARGS__)
#  define  MESSAGE_C(...)    CdoDebug::Message_(   caller , __VA_ARGS__)
#  define   MESSAGE(...)    CdoDebug::Message_( __func__ , __VA_ARGS__)
#else
#  define  SYS_ERROR(...)    CdoError::SysError_(__FILE__, __LINE__,"", __VA_ARGS__)
#  define    ERROR_C(...)    CdoError::Error_(__FILE__, __LINE__,"", __VA_ARGS__)
#  define    ERROR(...)      CdoError::Error_(__FILE__, __LINE__,"", __VA_ARGS__)
#  define   WARNING(...)     CdoError::Warning_(__VA_ARGS__)
#  define  MESSAGE_C(...)    CdoDebug::Message_(__VA_ARGS__)
#  define   MESSAGE(...)     CdoDebug::Message_(__func__,__VA_ARGS__)
#endif

#endif
