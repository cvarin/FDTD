
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
    bytes_allocated += allocate_1D_array_of_doubles(&ex,ncell,"ex");
    bytes_allocated += allocate_1D_array_of_doubles(&Dx,ncell,"Dx");
    bytes_allocated += allocate_1D_array_of_doubles(&hy,ncell,"hy");
  
    /**************************************************************************/
    print_allocated_memory_in_Kbytes();
    
    /*************************************************************************/
    // Set parameter for the polarization differential equation
    dt_dxeps0 = dt/(dx*epsi_0);
    dt_dxmu0 = dt/(dx*mu_0);
    gam = 2.0*relaxation_time/dt;
    
    /*************************************************************************/
    // Initialize the boundaries
    ex_low  = 0.0;
    ex_high = 0.0;
    
    /*************************************************************************/
    // Save the initial field distribution
    write_field_to_file(0,ex,hy);
}

/******************************************************************************/
em1d::~em1d()
{
    bytes_allocated -= free_array_of_doubles(hy,ncell);
    bytes_allocated -= free_array_of_doubles(Dx,ncell);
    bytes_allocated -= free_array_of_doubles(ex,ncell);
}
          
/****************** Member functions ******************************************/
void em1d::advance_a_step(const int _n)
{ 
    /**************************************************************************/
    // Update E-field
//     update_E();
//     update_E_with_D();
    update_E_with_P_and_epsi_rel();
    apply_boundary_E();
    update_source_for_E(_n);
  
    /**************************************************************************/
    // update the material response
//     update_polarization_debye_medium(ex,ncell);
    update_polarization_lorentz_medium(ex,ncell);
    
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
    if(_n%fileout==0) write_field_to_file(_n,ex,hy);
    
    /**************************************************************************/
    flush_output_buffers();
}

/******************************************************************************/
void em1d::apply_boundary_E()
{
    ex[0]  = ex_low;
    ex_low = ex[1];
    
    ex[ncell-1] = ex_high;
    ex_high = ex[ncell-2];
}

/******************************************************************************/
void em1d::apply_boundary_H()
{

}

/******************************************************************************/
double em1d::gaussian_pulse(const int _n, const int offset)
{
     double carrier = cos(2.0*Pi*freq_in*dt*(_n + offset*0.5) - ceo_phase);
     double enveloppe = exp(-0.5*pow((t0/time_scale- _n - offset*0.5)*time_scale/spread,2.0));
     return E0*carrier*enveloppe;
}

/******************************************************************************/
void em1d::update_E()
{
    #pragma omp parallel for
    for(int k=1; k < ncell; k++) ex[k] += dt_dxeps0/epsi_rel[k]*(hy[k-1] - hy[k]);
}

/******************************************************************************/
void em1d::update_E_with_D()
{
    #pragma omp parallel for
    for(int k=1; k < ncell; k++)
    {
        Dx[k] += dt/dx*(hy[k-1] - hy[k]);
        ex[k] = (Dx[k] - px[k])/(epsi_rel[k]*epsi_0);
    }
}

/******************************************************************************/
void em1d::update_E_with_P_and_epsi_rel()
{
      #pragma omp parallel for
      for(int k=1; k < ncell; k++)
        ex[k] += (dt_dxeps0*(hy[k-1] - hy[k]) 
                   - 1.0/epsi_0*(px[k] - px_previous[k]))/epsi_rel[k];
}

/******************************************************************************/
void em1d::update_H()
{
    #pragma omp parallel for
    for(int k=0; k < ncell-1; k++) hy[k] += dt_dxmu0*(ex[k] - ex[k+1]); 
}

/******************************************************************************/
void em1d::update_source_for_E(const int _n)
{   
     ex[source_plane + 95] += dt_dxeps0*gaussian_pulse(_n,1)/eta_0;
}

/******************************************************************************/
void em1d::update_source_for_H(const int _n)
{
     hy[source_plane + 95 - 1] += dt_dxmu0*gaussian_pulse(_n,0);
}

/****************** End of file ***********************************************/