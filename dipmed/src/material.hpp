#ifndef INC_material_hpp
#define INC_material_hpp

#include "IO.hpp"

class material : public IO
{    
     public:
          /******** Members ***************************************************/
          double *epsi_rel;    // Relative Permittivity
          double *density_profile; // Goes from 0.0 (no material) to 1.0 (material)
          double *px;          // Polarization
          double *px_previous; // Polarization (at time step n - 1)
          
          double gam;   // parameter for the polarization differential equation
          
          /******** Member functions ******************************************/
          material(const int _argc, const char **_argv);
          ~material();
          double static_response();
          void update_polarization_debye_medium(double *px, double *px_previous,
                                                 const double *ex, 
                                                  const double *density_profile,
                                                   const int ncell);
          void update_polarization_lorentz_medium(double *px, double *px_previous,
                                                   const double *ex, 
                                                    const double *density_profile,
                                                     const int ncell);
};

#endif // INC_material_hpp