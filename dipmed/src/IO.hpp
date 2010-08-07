#ifndef INC_IO_hpp
#define INC_IO_hpp

#include <string>
#include <cstdio>

struct IO
{    
     public:
          /******** Members ***************************************************/
          int nsteps;
          int step;
          int ncell;
          int m2start;
          int m2stop;
          
          double dx;
          double time_scale;
          double dt;
          
          std::string input_file;
          std::string output_dir;
          
          /******** Member functions ******************************************/
          IO(int argc, char **argv);
          void copy_input_file();
          void read_input_file();
          void write_field_to_file(const int n, const double *ex, 
                                    const double *hy);
};

#endif // INC_IO_hpp