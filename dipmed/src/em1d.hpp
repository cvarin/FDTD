#ifndef INC_em1d_hpp
#define INC_em1d_hpp

#include "IO.hpp"
#include "material.hpp"

class em1d : public material
{    
     public:
          /******** Members ***************************************************/
          double *ex;          // E-field
          double *Dx;          // D-field
          double *hy;          // H-field
          double *epsi_rel;    // Relative Permittivity
          double *density_profile; // Goes from 0.0 (no material) to 1.0 (material)
          double *px;          // Polarization
          double *px_previous; // Polarization (at time step n - 1)

          int bytes_allocated;
          
          double ex_low ; // Temp variables for absorbing boundaries
          double ex_high;
          
          static const int source_plane = 5;

          double dt_dxeps0; // Parameter for the E-field update
          double dt_dxmu0; // Parameter for the H-field update

          /******** Member functions ******************************************/
          em1d(const int _argc, const char **_argv);
          ~em1d();
          void advance_a_step(const int _n);
          void apply_boundary_E();
          void apply_boundary_H();
          void print_allocated_memory_in_Kbytes();
          void print_allocated_memory_in_Mbytes();
          void update_E();
          void update_E_with_D();
          void update_E_with_P_and_epsi_rel();
          void update_H();
          void update_source_E(const int _n);
          void update_source_H(const int _n);
};

#endif // INC_em1d_hpp