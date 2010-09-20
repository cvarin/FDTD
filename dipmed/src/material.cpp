
#include "material.hpp"

#include <omp.h>
#include <string.h>

#include "constants.hpp"
#include "system.hpp"

/****************** Constructor/Destructor ************************************/
material::material(const int _argc, const char **_argv):IO(_argc,_argv)
{
     bytes_allocated += allocate_1D_array_of_doubles(&epsi_rel,ncell,"epsi_rel");
     bytes_allocated += allocate_1D_array_of_doubles(&density_profile,ncell,"density_profile");
     bytes_allocated += allocate_1D_array_of_doubles(&Px,ncell,"P");
     bytes_allocated += allocate_1D_array_of_doubles(&Px_previous,ncell,"P_previous");
     
     /**************************************************************************/
     // The the parameters for the differential equations (obtained by 
     // putting finite differences in the differential equations).
     gam = 2.0*relaxation_time/dt;
     double dT = omega0*dt;
     double common = 1.0/(1 + 0.5*lifetime*dt);
     a1 = (2.0 - dT*dT)*common;
     a2 = -(1.0 - 0.5*lifetime*dt)*common;
     a3 = dT*dT*common;
     
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
     bytes_allocated -= free_array_of_doubles(Px_previous,ncell);
     bytes_allocated -= free_array_of_doubles(Px,ncell);
     bytes_allocated -= free_array_of_doubles(density_profile,ncell);
     bytes_allocated -= free_array_of_doubles(epsi_rel,ncell);
}
          
/****************** Member functions ******************************************/
double material::static_absorption()
{
    return conductivity;
}

/******************************************************************************/
double material::static_electronic_response()
{
    return epsilon - 1.0;
}

/******************************************************************************/
void material::update_polarization_debye_medium(const double *ex, const int ncell)
{    
    double Pstat;
    memcpy(Px_previous,Px,ncell*sizeof(double));
    #pragma omp parallel for default(none) shared(ex) private(Pstat) 
    for(int k=0; k < ncell-1; k++)
    {
        Pstat = density_profile[k]*static_absorption()*ex[k]*epsi_0;
        Px[k] = Pstat/(1+gam) - (1-gam)/(1+gam)*Px[k];
    }
}

/******************************************************************************/
void material::update_polarization_lorentz_medium(const double *ex, const int ncell)
{    
    double pstat;
    double p_minus_1;
    double p_minus_2;
    #pragma omp parallel for default(none) shared(ex) private(pstat,p_minus_1,p_minus_2) 
    for(int k=0; k < ncell-1; k++)
    {
        pstat = density_profile[k]*static_electronic_response()*ex[k]*epsi_0;
        p_minus_2 = Px_previous[k];
        Px_previous[k] = Px[k];
        p_minus_1 = Px_previous[k];
        Px[k] = a1*p_minus_1 + a2*p_minus_2 + a3*pstat;
    }
}

/****************** End of file ***********************************************/