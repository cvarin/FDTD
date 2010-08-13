
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
    ex = (double *)calloc(ncell,sizeof(double)); assert(ex);
    hy = (double *)calloc(ncell,sizeof(double)); assert(hy);
    cb = (double *)calloc(ncell,sizeof(double)); assert(cb);
    
    /**************************************************************************/
    if(ex==NULL) printf("Can't allocate ex\n"),abort();
    if(hy==NULL) printf("Can't allocate hy\n"),abort();
    if(cb==NULL) printf("Can't allocate cb\n"),abort();
  
    /**************************************************************************/
    bytes_allocated += 3*ncell*sizeof(double);
    print_allocated_memory_in_Kbytes();
    
    /**************************************************************************/
    // Set free space everywhere
    for(int k=0; k <= ncell-1; k++) cb[k] = 1.0; 
        
    /*************************************************************************/
    // Initialize the medium 2        
    for(int k=m2start; k < m2stop; k++) cb[k] = 1.0/epsilon;
    
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
    free(cb);
    free(hy);
    free(ex);
    bytes_allocated -= 3*ncell*sizeof(double);
    print_allocated_memory_in_Kbytes();
}
          
/****************** Member functions ******************************************/
void em1d::advance_a_step(const int _n)
{
    /**************************************************************************/
    // Update E-field
    update_E();
    apply_boundary_E();
    update_source_E(_n);
    
    /**************************************************************************/
    // Update H-field
    update_H();
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
void em1d::update_E()
{
    #pragma omp parallel for
    for(int k=1; k < ncell; k++) ex[k] += cb[k]*time_scale*(hy[k-1] - hy[k]); 
}

/******************************************************************************/
void em1d::update_H()
{
    #pragma omp parallel for
    for(int k=0; k < ncell-1; k++) hy[k] += time_scale*(ex[k] - ex[k+1]); 
}

/******************************************************************************/
void em1d::update_source_E(const int _n)
{
    double carrier = sin(2.0*Pi*freq_in*dt*_n);
    double enveloppe = exp(-0.5*pow((t0-_n)/spread,2.0));
    ex[source_plane] += 2.0*carrier*enveloppe;
}

/******************************************************************************/
void em1d::update_source_H(const int _n)
{
  
}

/****************** End of file ***********************************************/