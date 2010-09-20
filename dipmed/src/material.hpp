#ifndef INC_material_hpp
#define INC_material_hpp

#include "IO.hpp"

class material : public IO
{    
     public:
          /******** Members ***************************************************/
          double *epsi_rel;    // Relative Permittivity
          double *density_profile; // Goes from 0.0 (no material) to 1.0 (material)
          double *Px;          // Polarization
          double *Px_previous; // Polarization (at time step n - 1)
          
          double gam;   // Parameter for the Debye differential equation
          double a1;    // Parameters for the Lorentz differential equation
          double a2;
          double a3; 
          
          /******** Member functions ******************************************/
          material(const int _argc, const char **_argv);
          ~material();
          double static_absorption();
          double static_electronic_response();
          void update_polarization(const double *Ex, const int ncell);
};

#endif // INC_material_hpp