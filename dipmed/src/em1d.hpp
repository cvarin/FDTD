#ifndef INC_em1d_hpp
#define INC_em1d_hpp

#include "IO.hpp"
#include "material.hpp"

class em1d : public material
{    
     public:
          /******** Members ***************************************************/
          double *Ex;          // E-field
          double *Hy;          // H-field
          
          double ex_low ; // Temp variables for absorbing boundaries
          double ex_high;
          
          static const int source_plane = 5;

          double dt_dxeps0; // Parameter for the E-field update
          double dt_dxmu0;  // Parameter for the H-field update
          double dx_co;     // Parameter for the source updates

          /******** Member functions ******************************************/
          em1d(const int _argc, const char **_argv);
          ~em1d();
          void advance_a_step(const int _n);
          void apply_boundary_E();
          void apply_boundary_H();
          double gaussian_pulse(const int _n, const int offset);
          void update_E();
          void update_H();
          void update_source_for_E(const int _n);
          void update_source_for_H(const int _n);
};

#endif // INC_em1d_hpp