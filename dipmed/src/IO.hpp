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
          double conductivity; // Medium conductivity
          double number_density; // Number of atoms/molecules per m^3
          double relaxation_time; // In seconds.
          double omega0;          // Electronic resonance angular frequency.
          double lifetime;        // Electronic relaxation time.
          
          double dx;         // Cell size
          static const double time_scale = 0.5; // time scale = 1.0 is Courant condition
          // time_scale is fixed to 1 because of the rudimentary boundaries used.
          double dt;         // Time increment
          
          double t0;         // Center of the gaussian pulse in time steps
          double spread;     // 1/e width of the gaussian pulse
          double freq_in;    // Carrier frequency
          double I;          // Laser intensity in W/cm^2
          double E0;         // Laser E-field amplitude in V/m
          double ceo_phase;  // Carrier-envelope offset phase (http://www.rp-photonics.com/carrier_envelope_offset.html)
          
          std::string input_file;
          std::string output_dir;
          
          int screenout;     // Period for screen output
          int fileout;       // Period for file output
          
          int bytes_allocated;
          
          /******** Member functions ******************************************/
          IO(const int _argc, const char **_argv);
          ~IO();
          void copy_input_file();
          void print_allocated_memory_in_Kbytes();
          void print_allocated_memory_in_Mbytes();
          void read_input_file();
          void write_field_to_file(const int n, const double *ex, 
                                    const double *hy);
};

#endif // INC_IO_hpp