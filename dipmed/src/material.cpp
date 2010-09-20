
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
void material::update_polarization(const double *Ex, const int ncell)
{    
    double Pstat;
    memcpy(Px_previous,Px,ncell*sizeof(double));
    #pragma omp parallel for default(none) shared(Ex) private(Pstat) 
    for(int k=0; k < ncell-1; k++)
    {
        Pstat = density_profile[k]*static_absorption()*Ex[k]*epsi_0;
        Px[k] = Pstat/(1+gam) - (1-gam)/(1+gam)*Px[k];
    }
}

/****************** End of file ***********************************************/