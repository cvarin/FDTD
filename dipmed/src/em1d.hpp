#ifndef INC_em1d_hpp
#define INC_em1d_hpp

#include "IO.hpp"

class em1d : public IO
{    
     public:
          /******** Members ***************************************************/
          double *ex;    // e-field container
          double *hy;    // h-field container
          double *cb;    // Stores medium profile

          int bytes_allocated;
          
          double ex_low ; // Temp variables for absorbing boundaries
          double ex_high;

          /******** Member functions ******************************************/
          em1d(const int _argc, const char **_argv);
          ~em1d();
          void advance_a_step(const int _n);
          void print_allocated_memory_in_Kbytes();
          void print_allocated_memory_in_Mbytes();
};

#endif // INC_em1d_hpp