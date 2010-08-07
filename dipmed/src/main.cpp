#include <cmath>
#include <fstream>
#include <iostream>

#include "constants.hpp"
#include "IO.hpp"

void write_field_to_file(const int n, const int ncell, const double *ex, 
                          const double *hy);

/****************** Main ******************************************************/
int main(int argc, char **argv)
{   
     const char *outputdir = "./output/";
     IO(argc,argv);
     
     /************* Declarations **********************************************/
     const int    ncell = 400;     // Number of cells
     const double dx = 0.01,       // Cell size [m]
                  dt = dx/(2*co);  // Time step [s]
     int          NSTEPS = 1000;   // Number of time steps
     int          screenout = NSTEPS/10;  
     int          fileout  = 10;  

     const int    thickness = 50;
     const int    center  = (int)ncell/2;  
     const int    kstart  = center;  
     const int    kstop   = ncell; 
//      const int    kstart  = center - (int)center/2;  
//      const int    kstop   = center + (int)center/2; 
     double       epsilon = 4.0;  /* Relative dielectric constant of medium 2 */

     const double t0 = 80.0,      // Center of the incident pulse
                  spread = 40.0,  // Width of the incident pulse
                  freq_in = 2.0e9;// Signal Frequency [Hz]
     double       carrier = 0.0,  // Signal carrier
                  enveloppe = 0.0;// Signal enveloppe

     double    T = 0.0;        // Time

     double    ex_low_1  = 0.0,// Temp variables for
               ex_low_2  = 0.0,// absorbing boundaries
               ex_high_1 = 0.0,
               ex_high_2 = 0.0,
               ex_high_3 = 0.0, 
               ex_high_4 = 0.0;

     /*************************************************************************/
     // Initialize the E field and all cells to free space
     double       ex[ncell], hy[ncell], // Field containers
                  cb[ncell];         // Stores medium profile
     for (int k=0; k <= ncell-1; k++)
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
     write_field_to_file(0,ncell,ex,hy);
     if(NSTEPS <= 1) NSTEPS = 2;
     for(int n=1; n <= NSTEPS; n++)
     {
          T += 1.0;   // T keeps track of the number of times FDTD loop
                    // is executed.
                         
          // Calculate the Ex field
          for (int k=1; k < ncell; k++) ex[k] += cb[k]*0.5*( hy[k-1] - hy[k] ); 
          
          // Put a Gaussian pulse in the middle
          carrier = sin(2.0*Pi*freq_in*dt*T);
          enveloppe = exp(-0.5*pow((t0-T)/spread,2.0));
          ex[5] += carrier*enveloppe;
          
          // Absorbing boundary conditions for Ex
          ex[0]     = ex_low_2;
          ex_low_2  = ex_low_1;
          ex_low_1  = ex[1];
          
          ex[ncell-1]  = ex_high_4;
          ex_high_4 = ex_high_3;
          ex_high_3 = ex_high_2;
          ex_high_2 = ex_high_1;
          ex_high_1 = ex[ncell-2];
          
          // Calculate the Hy field
          for(int k=0; k < ncell-1; k++) hy[k] += 0.5*( ex[k] - ex[k+1] ); 
          
          // Outputs the progress of the calculation on screen
          if(n%screenout==0) printf("Progress = %d/%d\n",n,NSTEPS);
          
          /********************************************************************/
          // Writing fields to file                
          if(n%fileout==0) write_field_to_file(n,ncell,ex,hy);
          
          /********************************************************************/
          // Flush output buffers
          fflush(stdout);
     }
          
     return 0;
}

/******************************************************************************/
/****************** Local functions *******************************************/
/******************************************************************************/
void write_field_to_file(const int n, const int ncell, const double *ex, 
                          const double *hy)
{
     char filename[200];
     sprintf(filename,"./output/out_%.06d.dat",n);
     std::ofstream fp;
     fp.open(filename,std::ios::out);
     fp.precision(6);
     // Loop goes from 0 to KE-1 (excluded) because hy[KE-1] requires
     // ex[KE] that is undefined.
     for(int k=0; k < ncell-1; k++)
     {   // Typecast to (float) to avoid writing numbers with large exponents in the file
          fp << (float)k     << "\t" 
          << (float)ex[k] << "\t"
          << (float)hy[k] << std::endl;
     }      
     fp.close();
}

/****************** End of file ***********************************************/
