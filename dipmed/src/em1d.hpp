#ifndef INC_em1d_hpp
#define INC_em1d_hpp

#include "IO.hpp"

class em1d : public IO
{    
     public:
          /******** Members ***************************************************/
          static const int number_of_grids = 6;
          double *ex;          // e-field
          double *ex_previous; // e-field (at time step n - 1)
          double *hy;          // h-field
          double *cb;          // medium profile
          double *P;           // polarization
          double *P_previous;  // polarization (at time step n - 1)

          int bytes_allocated;
          
          double ex_low ; // Temp variables for absorbing boundaries
          double ex_high;
          
          static const int source_plane = 5;

          /******** Member functions ******************************************/
          em1d(const int _argc, const char **_argv);
          ~em1d();
          void advance_a_step(const int _n);
          void apply_boundary_E();
          void apply_boundary_H();
          void print_allocated_memory_in_Kbytes();
          void print_allocated_memory_in_Mbytes();
          void update_E(const double t_scale);
          void update_H(const double t_scale);
          void update_medium();
          void update_source_E(const int _n);
          void update_source_H(const int _n);
};

#endif // INC_em1d_hpp