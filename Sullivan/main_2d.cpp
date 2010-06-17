#include <cmath>
#include <fstream>
#include <iostream>

#include "PhysConsts.hpp"

int main()
{   
    const int    IE = 60,        // Number of cells in x
                 JE = 60,        // Number of cells in y
                 ic = IE/2,      // Middle of the grid in x
                 jc = JE/2;      // Middle of the grid in y
//    const double dx = 0.01,      // Cell size [m]
//                 dt = dx/(2*co); // Time step [s]
    double       T = 0.0;        // Time
    const int    NSTEPS = 50;   // Number of time steps
    
    const double t0 = 20.0,      // Center of the incident pulse
                 spread = 3.0;   // Width of the incident pulse
                 
    double       dz[IE][JE],     // Field containers
                 ez[IE][JE],
                 hx[IE][JE],
                 hy[IE][JE], 
                 ga[IE][JE];     // Stores medium profile
                 
    int  n = 0, i = 0, j = 0;    // Loops iterators
            
    std::ofstream fp;
    
    // Initialize the E field and all cells to free space
    for (j=0; j < JE; j++){
        for (i=0; i < IE; i++){ 
            dz[i][j]   = 0.0;
            ez[i][j]   = 0.0;
            hx[i][j]   = 0.0;
            hy[i][j]   = 0.0;
            ga[i][j]   = 1.0; 
            } // end of i for
        } // end of j for
    
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
                           
            // Calculate the Dz field
            for (j=1; j < JE; j++){
                for (i=1; i < IE; i++){
                    dz[i][j] += 0.5*( hy[i][j] - hy[i-1][j]
                                      - hx[i][j] + hx[i][j-1] ); 
                    } // end of i for
                } // end of j for
            
            // Put a Gaussian pulse in the middle
            dz[ic][jc] = exp( -0.5*pow((t0-T)/spread,2.0) );

            // Calculate the Ez field
            for (j=1; j < JE; j++){
                for (i=1; i < IE; i++){
                    ez[i][j] = ga[i][j]*dz[i][j]; 
                    } // end of i for
                } // end of j for

            // Calculate the Hx field
            for (j=0; j < JE-1; j++){
                for (i=0; i < IE-1; i++){
                    hx[i][j] += 0.5*( ez[i][j] - ez[i][j+1] );
                    } // end of i for
                } // end of j for
            
            // Calculate the Hy field
            for (j=0; j < JE-1; j++){
                for (i=0; i < IE-1; i++){
                    hy[i][j] += 0.5*( ez[i+1][j] - ez[i][j] );
                    } // end of i for
                } // end of j for
            
            
            // Outputs the progress of the calculation on screen
            // Typecast to (int) for an output without digits
            std::cerr << "Progress = " << (int)(100*T/NSTEPS) << "% \r";
                
          } // end of FDTD for loop
        } // end of else
     
            
    // Writing fields to file                
    fp.open("Out.dat",std::ios::out);
    fp.precision(6);
    
    for (i=1; i < IE; i++){
        for (j=1; j < JE; j++){       // Typecast to (float) to avoid
        fp << (float)i        << "\t" // writing numbers with large
           << (float)j        << "\t" // exponents in the file   
           << (float)ez[i][j] << std::endl;
          } // end of i for
        } // end of j for
          
    fp.close();

       
    return 0;
} // end of main
