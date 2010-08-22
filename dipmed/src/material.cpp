
#include "material.hpp"

#include <omp.h>

#include "constants.hpp"
#include "system.hpp"

/****************** Constructor/Destructor ************************************/
material::material(const int _argc, const char **_argv):IO(_argc,_argv)
{
     bytes_allocated += allocate_1D_array_of_doubles(&epsi_rel,ncell,"epsi_rel");
     bytes_allocated += allocate_1D_array_of_doubles(&density_profile,ncell,"density_profile");
     bytes_allocated += allocate_1D_array_of_doubles(&px_previous,ncell,"P_previous");
     bytes_allocated += allocate_1D_array_of_doubles(&px,ncell,"P");
     gam = 2.0*relaxation_time/dt;
     
     /**************************************************************************/
     // Set free space everywhere
     for(int k=0; k <= ncell-1; k++)
     {
          epsi_rel[k] = 1.0;
          density_profile[k]  = 0.0;
     }
          
     /*************************************************************************/
     // Initialize the medium 2        
     for(int k=m2start; k < m2stop; k++)
     {
          epsi_rel[k] = epsilon;
          density_profile[k] = 1.0;
     }
}

/******************************************************************************/
material::~material()
{
     bytes_allocated -= free_array_of_doubles(px_previous,ncell);
     bytes_allocated -= free_array_of_doubles(px,ncell);
     bytes_allocated -= free_array_of_doubles(density_profile,ncell);
     bytes_allocated -= free_array_of_doubles(epsi_rel,ncell);
}
          
/****************** Member functions ******************************************/
double material::static_response()
{
    return conductivity;
}

/******************************************************************************/
void material::update_polarization_debye_medium(double *px, double *px_previous,
                                                 const double *ex, 
                                                  const double *density_profile,
                                                   const int ncell)
{    
    double pstat;
    #pragma omp parallel for shared(px,px_previous,ex,density_profile) private(pstat) 
    for(int k=0; k < ncell-1; k++)
    {
        px_previous[k] = px[k];
        pstat = density_profile[k]*static_response()*ex[k]*epsi_0;
        px[k] = pstat/(1+gam) - (1-gam)/(1+gam)*px[k];
    }
}

/******************************************************************************/
void material::update_polarization_lorentz_medium(double *px, double *px_previous,
                                                   const double *ex, 
                                                    const double *density_profile,
                                                     const int ncell)
{    
    double pstat;
    #pragma omp parallel for shared(px,px_previous,ex,density_profile) private(pstat) 
    for(int k=0; k < ncell-1; k++)
    {
        px_previous[k] = px[k];
        pstat = density_profile[k]*static_response()*ex[k]*epsi_0;
        px[k] = pstat/(1+gam) - (1-gam)/(1+gam)*px[k];
    }
}

/******************************************************************************/

/****************** End of file ***********************************************/