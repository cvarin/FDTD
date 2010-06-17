#include <cmath>
#include <fstream>
#include <iostream>

#include "PhysConsts.hpp"

int main()
{   
    const int    KE = 200,       // Number of cells
                 kc = KE/2;      // Middle of the grid
    const double dx = 0.01,      // Cell size [m]
                 dt = dx/(2*co); // Time step [s]
    double       T = 0.0;        // Time
    int          NSTEPS = 300;   // Number of time steps
    
    const int    kstart  = 100;  // Boundary between mediums 1 and 2
    double       epsilon = 4.0;  /* Relative dielectric constant of 
                                    medium 2 */
    
    const double t0 = 80.0,      // Center of the incident pulse
                 spread = 40.0,  // Width of the incident pulse
                 freq_in = 2.0e9;// Signal Frequency [Hz]
    double       carrier = 0.0,  // Signal carrier
                 enveloppe = 0.0;// Signal enveloppe
    
    double       ex[KE], hy[KE], // Field containers
                 cb[KE];         // Stores medium profile
    int          n = 0,          // Loop iterators
                 k = 0;                    
    
    double       ex_low_1  = 0.0,// Temp variables for
                 ex_low_2  = 0.0,// absorbing boundaries
                 ex_high_1 = 0.0,
                 ex_high_2 = 0.0,
                 ex_high_3 = 0.0, 
                 ex_high_4 = 0.0;
            
    std::ofstream fp;
    
    // Initialize the E field and all cells to free space
    for (k=0; k <= KE-1; k++)
        { ex[k]   = 0.0;
          hy[k]   = 0.0;
          cb[k]   = 1.0; }
          
    // Initialize the medium 2        
    for (k=kstart; k < KE; k++){ cb[k] = 1.0/epsilon; }
    
    // Check if there is at least 2 steps
    if (NSTEPS <= 1){ 
        std::cout << "Error: Number of steps must be larger than 1." 
                  << std::endl;
        getchar();
       }
       
    else{   
        // FDTD loop
        for (n=1; n <= NSTEPS; n++)
          {
            T += 1.0;   // T keeps track of the number of times FDTD loop
                        // is executed.
                           
            // Calculate the Ex field
            for (k=1; k < KE; k++){ 
                  ex[k] += cb[k]*0.5*( hy[k-1] - hy[k] ); 
                }
            
            // Put a Gaussian pulse in the middle
            carrier = sin(2.0*Pi*freq_in*dt*T);
            enveloppe = exp( -0.5*pow((t0-T)/spread,2.0) );
            ex[5] += carrier*enveloppe;
            
            // Absorbing boundary conditions for Ex
            ex[0]     = ex_low_2;
            ex_low_2  = ex_low_1;
            ex_low_1  = ex[1];
            
            ex[KE-1]  = ex_high_4;
            ex_high_4 = ex_high_3;
            ex_high_3 = ex_high_2;
            ex_high_2 = ex_high_1;
            ex_high_1 = ex[KE-2];
            
            // Calculate the Hy field
            for (k=0; k < KE-1; k++){ 
                  hy[k] += 0.5*( ex[k] - ex[k+1] ); 
                }
            
            // Outputs the progress of the calculation on screen
            // Typecast to (int) for an output without digits
            std::cerr << "Progress = " << (int)(100*T/NSTEPS) << "% \r";
                
          } // end of FDTD for loop
        } // end of else
            
    // Writing fields to file                
    fp.open("Out.dat",std::ios::out);
    fp.precision(6);
    
    // Loop goes from 0 to KE-1 (excluded) because hy[KE-1] requires
    // ex[KE] that is undefined.
    for (k=0; k < KE-1; k++){      // Typecast to (float) to avoid
        fp << (float)k     << "\t" // writing numbers with large
           << (float)ex[k] << "\t" // exponents in the file   
           << (float)hy[k] << std::endl;
          }
          
    fp.close();
       
    return 0;
} // end of main
