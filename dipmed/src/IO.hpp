#ifndef INC_IO_hpp
#define INC_IO_hpp

#include <string>

class IO
{    
     public:
          /******** Members ***************************************************/
          int threads;   // Numbers of threads to use for the calculation
          int nsteps;    // Number of time steps
          int step;      // step size
          int ncell;     // Number of cell
          
          int m2start;    // Cell where the medium starts
          int m2stop;     // Cell where the medium stops
          double epsilon; // Medium relative index
          double number_density; // Number of atoms/molecules per m^3
          double relaxation_time; // In seconds.
          
          double dx;         // Cell size
          static const double time_scale = 1.0; // time scale = 1.0 is Courant condition
          double dt;         // Time increment
          
          double t0;         // Center of the gaussian pulse in time steps
          double spread;     // 1/e width of the gaussian pulse
          double freq_in;    // Carrier frequency
          double I;          // Laser intensity in W/cm^2
          double E0;         // Laser E-field amplitude in V/m
          
          std::string input_file;
          std::string output_dir;
          
          int screenout;     // Period for screen output
          int fileout;       // Period for file output
          
          /******** Member functions ******************************************/
          IO(const int _argc, const char **_argv);
          void copy_input_file();
          void read_input_file();
          void write_field_to_file(const int n, const double *ex, 
                                    const double *hy);
};

#endif // INC_IO_hpp