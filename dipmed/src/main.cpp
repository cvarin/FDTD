#include <cmath>

#include "constants.hpp"
#include "IO.hpp"

/****************** Main ******************************************************/
int main(int argc, char **argv)
{   
     IO io(argc,argv);
     
     /************* Declarations **********************************************/
     int          screenout = io.nsteps/10;
     int          fileout  = 10;  

     /************* Medium ****************************************************/
     const int    thickness = 50;
     const int    center  = (int)io.ncell/2;  
     const int    kstart  = center;  
     const int    kstop   = io.ncell; 
//      const int    kstart  = center - (int)center/2;  
//      const int    kstop   = center + (int)center/2; 
     double       epsilon = 4.0;  /* Relative dielectric constant of medium 2 */

     /************* Source ****************************************************/
     const double t0 = 80.0,      // Center of the incident pulse
                  spread = 40.0,  // Width of the incident pulse
                  freq_in = 2.0e9;// Signal Frequency [Hz]
     double       carrier = 0.0,  // Signal carrier
                  enveloppe = 0.0;// Signal enveloppe

     /************* Boundary **************************************************/
     double    ex_low_1  = 0.0,// Temp variables for
               ex_low_2  = 0.0,// absorbing boundaries
               ex_high_1 = 0.0,
               ex_high_2 = 0.0,
               ex_high_3 = 0.0, 
               ex_high_4 = 0.0;

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
     for(int k=kstart; k < kstop; k++) cb[k] = 1.0/epsilon;

     /*************************************************************************/
     // FDTD loop
     double T = 0.0;        // Time
     io.write_field_to_file(0,ex,hy);
     if(io.nsteps <= 1) io.nsteps = 2;
     for(int n=1; n <= io.nsteps; n++)
     {
          T += 1.0;   // T keeps track of the number of times FDTD loop
                    // is executed.
                         
          // Calculate the Ex field
          for (int k=1; k < io.ncell; k++) ex[k] += cb[k]*0.5*(hy[k-1] - hy[k]); 
          
          // Put a Gaussian pulse in the middle
          carrier = sin(2.0*Pi*freq_in*io.dt*T);
          enveloppe = exp(-0.5*pow((t0-T)/spread,2.0));
          ex[5] += carrier*enveloppe;
          
          // Absorbing boundary conditions for Ex
          ex[0]     = ex_low_2;
          ex_low_2  = ex_low_1;
          ex_low_1  = ex[1];
          
          ex[io.ncell-1]  = ex_high_4;
          ex_high_4 = ex_high_3;
          ex_high_3 = ex_high_2;
          ex_high_2 = ex_high_1;
          ex_high_1 = ex[io.ncell-2];
          
          // Calculate the Hy field
          for(int k=0; k < io.ncell-1; k++) hy[k] += 0.5*(ex[k] - ex[k+1]); 
          
          // Outputs the progress of the calculation on screen
          if(n%screenout==0) printf("Progress = %d/%d\n",n,io.nsteps);
          
          /********************************************************************/
          // Writing fields to file                
          if(n%fileout==0) io.write_field_to_file(n,ex,hy);
          
          /********************************************************************/
          // Flush output buffers
          fflush(stdout);
     }
     return 0;
}

/****************** End of file ***********************************************/
