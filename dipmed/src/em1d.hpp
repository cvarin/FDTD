#ifndef INC_em1d_hpp
#define INC_em1d_hpp

#include "IO.hpp"

class em1d : public IO
{    
     public:
          /******** Members ***************************************************/
          double *ex;          // E-field
          double *Dx;          // D-field
          double *ex_previous; // E-field (at time step n - 1)
          double *hy;          // H-field
          double *epsi_rel;    // Relative Permittivity
          double *N;           // Average number of atom/molecule per m^3
          double *px;          // Polarization
          double *px_previous; // Polarization (at time step n - 1)

          int bytes_allocated;
          
          double ex_low ; // Temp variables for absorbing boundaries
          double ex_high;
          
          static const int source_plane = 5;

          double gam;   // parameter for the polarization differential equation
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
          double static_response(const int k);
          void update_E(const double t_scale);
          void update_E_with_P(const double t_scale);
          void update_H(const double t_scale);
          void update_polarization();
          void update_source_E(const int _n);
          void update_source_H(const int _n);
};

#endif // INC_em1d_hpp