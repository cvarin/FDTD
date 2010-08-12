#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "constants.hpp"
#include "em1d.hpp"

/****************** Main ******************************************************/
int main(int argc, char **argv)
{   
     em1d f(argc,argv);

     /*************************************************************************/     
     f.write_field_to_file(0,f.ex,f.hy);
     for(int n=1; n <= f.nsteps; n++)
     {                   
          // Calculate the Ex field
          for (int k=1; k < f.ncell; k++) f.ex[k] += f.cb[k]*f.time_scale*(f.hy[k-1] - f.hy[k]); 
          
          // Put a Gaussian pulse in the middle
          double carrier = sin(2.0*Pi*f.freq_in*f.dt*n);
          double enveloppe = exp(-0.5*pow((f.t0-n)/f.spread,2.0));
          f.ex[5] += carrier*enveloppe;
          
          // Absorbing boundary conditions for Ex
          f.ex[0]  = f.ex_low;
          f.ex_low = f.ex[1];
          
          f.ex[f.ncell-1] = f.ex_high;
          f.ex_high = f.ex[f.ncell-2];
          
          // Calculate the Hy field
          for(int k=0; k < f.ncell-1; k++) f.hy[k] += f.time_scale*(f.ex[k] - f.ex[k+1]); 
          
          // Outputs the progress of the calculation on screen
          if(n%f.screenout==0) printf("Progress = %d/%d\n",n,f.nsteps);
          
          /********************************************************************/
          // Writing fields to file                
          if(n%f.fileout==0) f.write_field_to_file(n,f.ex,f.hy);
          
          /********************************************************************/
          // Flush output buffers
          fflush(stdout);
          fflush(stderr);
     }
     return 0;
}

/****************** End of file ***********************************************/
