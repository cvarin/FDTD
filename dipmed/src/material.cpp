
#include "material.hpp"

#include <omp.h>

#include "constants.hpp"

/****************** Constructor/Destructor ************************************/
material::material(const int _argc, const char **_argv):IO(_argc,_argv)
{
     gam = 2.0*relaxation_time/dt;
}

/******************************************************************************/
material::~material()
{

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