
#include "em1d.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <omp.h>

#include "constants.hpp"
#include "system.hpp"

/****************** Constructor/Destructor ************************************/
em1d::em1d(const int _argc, const char **_argv):material(_argc,_argv)
{ 
    /**************************************************************************/
    // Initialize the E field and all cells to free space
    bytes_allocated += allocate_1D_array_of_doubles(&Ex,ncell,"Ex");
    bytes_allocated += allocate_1D_array_of_doubles(&Dx,ncell,"Dx");
    bytes_allocated += allocate_1D_array_of_doubles(&Hy,ncell,"Hy");
  
    /**************************************************************************/
    print_allocated_memory_in_Kbytes();
    
    /*************************************************************************/
    // Set parameter for the polarization differential equation
    dt_dxeps0 = dt/(dx*epsi_0);
    dt_dxmu0 = dt/(dx*mu_0);
    dx_co = dx/co;
    gam = 2.0*relaxation_time/dt;
    
    /*************************************************************************/
    // Initialize the boundaries
    ex_low  = 0.0;
    ex_high = 0.0;
    
    /*************************************************************************/
    // Save the initial field distribution
    write_field_to_file(0,Ex,Hy);
}

/******************************************************************************/
em1d::~em1d()
{
    bytes_allocated -= free_array_of_doubles(Hy,ncell);
    bytes_allocated -= free_array_of_doubles(Dx,ncell);
    bytes_allocated -= free_array_of_doubles(Ex,ncell);
}
          
/****************** Member functions ******************************************/
void em1d::advance_a_step(const int _n)
{ 
    /**************************************************************************/
    // Update E-field
//     update_E();

//     update_E_with_D();
//     apply_boundary_E();
//     update_source_for_D(_n);
    
    update_E_with_P_and_epsi_rel();
    apply_boundary_E();
    update_source_for_E(_n);
  
    /**************************************************************************/
    // update the material response
//     update_polarization_debye_medium(ex,ncell);
    update_polarization_lorentz_medium(Ex,ncell);
    
    /**************************************************************************/
    // Update H-field
    update_H();
    apply_boundary_H();
    update_source_for_H(_n);
    
    /**************************************************************************/
    // Outputs the progress of the calculation on screen
    if(_n%screenout==0) printf("Progress = %d/%d\n",_n,nsteps);
    
    /**************************************************************************/
    // Writing fields to file                
    if(_n%fileout==0) write_field_to_file(_n,Ex,Hy);
    
    /**************************************************************************/
    flush_output_buffers();
}

/******************************************************************************/
void em1d::apply_boundary_E()
{
    Ex[0]  = ex_low;
    ex_low = Ex[1];
    
    Ex[ncell-1] = ex_high;
    ex_high = Ex[ncell-2];
}

/******************************************************************************/
void em1d::apply_boundary_H()
{

}

/******************************************************************************/
double em1d::gaussian_pulse(const int _n, const int offset)
{
    const double t = (_n - offset*0.5*(1.0-1.0/time_scale))*dt;
    const double T = (t - t0*dx_co)/(spread*dx_co);
    const double carrier = cos(omega_laser*t - ceo_phase);
    const double enveloppe = exp(-T*T);
    return E0*carrier*enveloppe;
}

/******************************************************************************/
void em1d::update_E()
{
    #pragma omp parallel for default(none)
    for(int k=1; k < ncell; k++) Ex[k] += dt_dxeps0/epsi_rel[k]*(Hy[k-1] - Hy[k]);
}

/******************************************************************************/
void em1d::update_E_with_D()
{
    #pragma omp parallel for default(none)
    for(int k=1; k < ncell; k++)
    {
        Dx[k] += dt/dx*(Hy[k-1] - Hy[k]);
        Ex[k] = (Dx[k] - px[k])/(epsi_rel[k]*epsi_0);
    }
}

/******************************************************************************/
void em1d::update_E_with_P_and_epsi_rel()
{
      #pragma omp parallel for default(none)
      for(int k=1; k < ncell; k++)
        Ex[k] += (dt_dxeps0*(Hy[k-1] - Hy[k]) 
                   - 1.0/epsi_0*(px[k] - px_previous[k]))/epsi_rel[k];
}

/******************************************************************************/
void em1d::update_H()
{
    #pragma omp parallel for default(none)
    for(int k=0; k < ncell-1; k++) Hy[k] += dt_dxmu0*(Ex[k] - Ex[k+1]); 
}

/******************************************************************************/
void em1d::update_source_for_D(const int _n)
{   
     Dx[source_plane] += epsi_0*dt_dxeps0*gaussian_pulse(_n,1)/eta_0;
}

/******************************************************************************/
void em1d::update_source_for_E(const int _n)
{   
     Ex[source_plane] += dt_dxeps0*gaussian_pulse(_n,1)/eta_0;
}

/******************************************************************************/
void em1d::update_source_for_H(const int _n)
{
     Hy[source_plane - 1] += dt_dxmu0*gaussian_pulse(_n,0);
}

/****************** End of file ***********************************************/