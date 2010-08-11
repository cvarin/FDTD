#include <cmath>

#include "constants.hpp"
#include "IO.hpp"

/****************** Main ******************************************************/
int main(int argc, char **argv)
{   
     IO io(argc,argv);

     /************* Boundary **************************************************/
     double    ex_low_1  = 0.0,// Temp variables for
               ex_low_2  = 0.0,// absorbing boundaries
               ex_high_1 = 0.0,
               ex_high_2 = 0.0;

     /*************************************************************************/
     // Initialize the E field and all cells to free space
     double       ex[io.ncell], hy[io.ncell], // Field containers
                  cb[io.ncell];         // Stores medium profile
     for (int k=0; k <= io.ncell-1; k++)
     { 
          ex[k]   = 0.0;
          hy[k]   = 0.0;
          cb[k]   = 1.0; 
     }
          
     /*************************************************************************/
     // Initialize the medium 2        
     for(int k=io.m2start; k < io.m2stop; k++) cb[k] = 1.0/io.epsilon;

     /*************************************************************************/     
     io.write_field_to_file(0,ex,hy);
     for(int n=1; n <= io.nsteps; n++)
     {                   
          // Calculate the Ex field
          for (int k=1; k < io.ncell; k++) ex[k] += cb[k]*0.5*(hy[k-1] - hy[k]); 
          
          // Put a Gaussian pulse in the middle
          double carrier = sin(2.0*Pi*io.freq_in*io.dt*n);
          double enveloppe = exp(-0.5*pow((io.t0-n)/io.spread,2.0));
          ex[5] += carrier*enveloppe;
          
          // Absorbing boundary conditions for Ex
          ex[0]     = ex_low_2;
          ex_low_2  = ex_low_1;
          ex_low_1  = ex[1];
          
          ex[io.ncell-1] = ex_high_2;
          ex_high_2 = ex_high_1;
          ex_high_1 = ex[io.ncell-2];
          
          // Calculate the Hy field
          for(int k=0; k < io.ncell-1; k++) hy[k] += 0.5*(ex[k] - ex[k+1]); 
          
          // Outputs the progress of the calculation on screen
          if(n%io.screenout==0) printf("Progress = %d/%d\n",n,io.nsteps);
          
          /********************************************************************/
          // Writing fields to file                
          if(n%io.fileout==0) io.write_field_to_file(n,ex,hy);
          
          /********************************************************************/
          // Flush output buffers
          fflush(stdout);
     }
     return 0;
}

/****************** End of file ***********************************************/
