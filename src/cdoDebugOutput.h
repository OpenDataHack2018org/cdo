#ifndef DEBUG_SWITCHES_H
#define DEBUG_SWITCHES_H
#include <iostream>
#include <sstream>

#include <fstream>

namespace CdoDebug
{
    extern int  PSTREAM;
    extern bool PROCESS;
    extern std::string outfile;
    extern bool print_to_seperate_file;

    namespace{
        void printMessageToFile(std::stringstream &p_message)
        {
            if(!outfile.empty()){
                std::fstream outfile_stream(outfile,std::fstream::in | std::fstream::app );

                outfile_stream <<  p_message.str();
            }
        }
    }


    template <typename ...T>
    void Message_ (std::stringstream &p_message, T&& ...args)
    {
           //for showing that the dummy array is never used
           using expander = int[];
           //creating dummy array with expanding the parameter pack in its
           //initializer list while writing each element of the pack into message
           expander{0,(void(p_message << std::forward<T>(args)),0)...};
           p_message << std::endl;
    }

    template <typename ...T>
    void Message_ (const char * p_func, T&& ...args)
    {
        std::stringstream message;
        message << p_func <<": ";
        Message_(message, args...);
        if(print_to_seperate_file)
        {
            printMessageToFile(message);
        }
        else
        {
            std::cout << message.str();
        }
    }
}

namespace CdoError{
    static int _ExitOnError = 1;

    template <typename ...T>
    void Error_(const char* p_file, const int p_line, const char* caller, T&& ...args)
    {
          std::stringstream message;
          message << "Error in: " << p_file << ":" << p_line << " ";
          CdoDebug::Message_(message, args...);
          std::cout << message.str();
          if ( CdoError::_ExitOnError )
          {
              exit(EXIT_FAILURE);
          }
    }
    template <typename ...T>
    void Warning_(T&& ...args)
    {
        std::stringstream message;
        message << "Warning: ";
        CdoDebug::Message_(message, args...);
        std::cout << message.str();
    }

    template <typename ...T>
    void SysError_(const char* p_file, const int p_line,  const char* p_func, T&& ...args)
    {
        int saved_errno = errno;
        std::stringstream message;
        message << "SysError in:" << p_file << std::endl;
        message << "    " << "in function: p_func ,line: " << p_line << std::endl;
        CdoDebug::Message_(message, args...);

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
