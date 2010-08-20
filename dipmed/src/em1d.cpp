
#include "em1d.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <omp.h>

#include "constants.hpp"
#include "system.hpp"

/****************** Constructor/Destructor ************************************/
em1d::em1d(const int _argc, const char **_argv):IO(_argc,_argv)
{
    bytes_allocated = 0;
  
    /**************************************************************************/
    // Initialize the E field and all cells to free space
    bytes_allocated += allocate_1D_array_of_doubles(&ex,ncell,"ex");
    bytes_allocated += allocate_1D_array_of_doubles(&Dx,ncell,"Dx");
    bytes_allocated += allocate_1D_array_of_doubles(&hy,ncell,"hy");
    bytes_allocated += allocate_1D_array_of_doubles(&epsi_rel,ncell,"epsi_rel");
    bytes_allocated += allocate_1D_array_of_doubles(&N,ncell,"N");
    bytes_allocated += allocate_1D_array_of_doubles(&px_previous,ncell,"P_previous");
    bytes_allocated += allocate_1D_array_of_doubles(&px,ncell,"P");
  
    /**************************************************************************/
    print_allocated_memory_in_Kbytes();
    
    /**************************************************************************/
    // Set free space everywhere
    for(int k=0; k <= ncell-1; k++)
    {
        epsi_rel[k] = 1.0;
        N[k]  = 0.0;
    }
        
    /*************************************************************************/
    // Initialize the medium 2        
    for(int k=m2start; k < m2stop; k++)
    {
        epsi_rel[k] = epsilon;
        N[k]  = number_density;
    }
    
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
    bytes_allocated -= free_array_of_doubles(px_previous,ncell);
    bytes_allocated -= free_array_of_doubles(px,ncell);
    bytes_allocated -= free_array_of_doubles(N,ncell);
    bytes_allocated -= free_array_of_doubles(epsi_rel,ncell);
    bytes_allocated -= free_array_of_doubles(hy,ncell);
    bytes_allocated -= free_array_of_doubles(Dx,ncell);
    bytes_allocated -= free_array_of_doubles(ex,ncell);
    print_allocated_memory_in_Kbytes();
}
          
/****************** Member functions ******************************************/
void em1d::advance_a_step(const int _n)
{ 
    /**************************************************************************/
    // Update E-field
//     update_E(time_scale);
//     update_E_with_D(time_scale);
    update_E_with_P(time_scale);
    apply_boundary_E();
    update_source_E(_n);
  
    /**************************************************************************/
    // update the material response
    update_polarization();
    
    /**************************************************************************/
    // Update H-field
    update_H(time_scale);
    apply_boundary_H();
    update_source_H(_n);
    
    /**************************************************************************/
    // Outputs the progress of the calculation on screen
    if(_n%screenout==0) printf("Progress = %d/%d\n",_n,nsteps);
    
    /********************************************************************/
    // Writing fields to file                
    if(_n%fileout==0) write_field_to_file(_n,ex,hy);
    
    /********************************************************************/
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
void em1d::print_allocated_memory_in_Kbytes()
{
    printf("Memory allocated: %f KB\n",(double)bytes_allocated/1024);
}

/******************************************************************************/
void em1d::print_allocated_memory_in_Mbytes()
{
    printf("Memory allocated: %f MB\n",(double)bytes_allocated/1024/1024);
}

/******************************************************************************/
double em1d::static_response(const int k)
{
    return epsi_rel[k] - 1.0;
}

/******************************************************************************/
void em1d::update_E(const double t_scale)
{
    #pragma omp parallel for
    for(int k=1; k < ncell; k++) ex[k] += dt_dxeps0/epsi_rel[k]*(hy[k-1] - hy[k]);
}

/******************************************************************************/
void em1d::update_E_with_D(const double t_scale)
{
    #pragma omp parallel for
    for(int k=1; k < ncell; k++)
    {
        Dx[k] += dt/dx*(hy[k-1] - hy[k]);
        ex[k] = (Dx[k] - px[k])/(epsi_rel[k]*epsi_0);
    }
}

/******************************************************************************/
void em1d::update_E_with_P(const double t_scale)
{
      #pragma omp parallel for
      for(int k=1; k < ncell; k++)
        ex[k] += (dt_dxeps0*(hy[k-1] - hy[k]) 
                   - 1.0/epsi_0*(px[k] - px_previous[k]))/epsi_rel[k];
}

/******************************************************************************/
void em1d::update_H(const double t_scale)
{
    #pragma omp parallel for
    for(int k=0; k < ncell-1; k++) hy[k] += dt_dxmu0*(ex[k] - ex[k+1]); 
}

/******************************************************************************/
void em1d::update_polarization()
{
    double pstat;    
    #pragma omp parallel for
    for(int k=0; k < ncell-1; k++) 
    {
        px_previous[k] = px[k];
        pstat = static_response(k)*ex[k]*epsi_0;
        px[k] = pstat/(1+gam) - (1-gam)/(1+gam)*px[k];
    }
}

/******************************************************************************/
void em1d::update_source_E(const int _n)
{
    double carrier = cos(2.0*Pi*freq_in*dt*_n - ceo_phase);
    double enveloppe = exp(-0.5*pow((t0-_n)/spread,2.0));
    ex[source_plane] += 2.0*E0*carrier*enveloppe;
}

/******************************************************************************/
void em1d::update_source_H(const int _n)
{
  
}

/****************** End of file ***********************************************/