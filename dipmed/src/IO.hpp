#ifndef INC_IO_hpp
#define INC_IO_hpp

#include <string>
#include <cstdio>

struct IO
{    
     public:
          const static int nparams = 5;
          int nsteps;
          int step;
          int ncell;
          int m2start;
          int m2stop;
          
          std::string input_file;
          std::string output_dir;
          std::string output_file;
          
          IO(int argc, char **argv);
          void copy_input_file();
          void read_input_file();
};

#endif // INC_IO_hpp