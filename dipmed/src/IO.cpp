
#include "IO.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>

#include "constants.hpp"

/****************** Constructor/Destructor ************************************/
IO::IO(int argc, char **argv)
{
     if(argc < 3)
     {
          printf("\n");
          printf("Please provide an input file name and output directory.\n");
          printf("Example: ./program <input file with path> <output directory>\n");
          printf("\n");
          abort();
     }
     
     /************* Read the input file ***************************************/
     this->input_file = argv[1];
     this->output_dir = argv[2];
     printf("\n");
     this->read_input_file();
     this->copy_input_file();
}
          
/****************** Member functions ******************************************/
void IO::copy_input_file()
/*
     Reference : http://hubpages.com/hub/File-Copy-Program-in-C-Language
*/
{
     char ch; 
     std::string output_file = this->output_dir + "input.txt";
     FILE *inF = fopen(this->input_file.c_str(),"r");
     FILE *ouF = fopen(output_file.c_str(),"w");
     if(inF == NULL) perror(this->input_file.c_str()), exit(-1);
     if(ouF == NULL) perror(output_file.c_str()), exit(-1);
     
     while(1)  
     {  
          ch = getc(inF);  
          if(ch==EOF) break;  
          else putc(ch,ouF);  
     }
     
     fclose(inF);
     fclose(ouF);
     printf("Input file copied to %s\n",this->output_dir.c_str());
     printf("\n");
}

/******************************************************************************/
void IO::read_input_file()
{
     FILE *f = fopen(this->input_file.c_str(),"r");
     if(f==NULL)
     {
          printf("Can't open %s\n",this->input_file.c_str());
          abort();
     }
     
     /************* Read the file *********************************************/
     printf("Reading %s\n",this->input_file.c_str());
     int nparams = 6;
     int results[nparams];
     results[0] = fscanf(f,"nsteps = %i\n", &this->nsteps);
     results[1] = fscanf(f,"step = %i\n", &this->step);
     results[2] = fscanf(f,"ncell = %i\n", &this->ncell);
     results[3] = fscanf(f,"m2start = %i\n", &this->m2start);
     results[4] = fscanf(f,"m2stop = %i\n", &this->m2stop);
     results[5] = fscanf(f,"dx = %lf\n", &this->dx);
     results[6] = fscanf(f,"time_scale = %lf\n", &this->time_scale);
     this->dt = this->time_scale*dx/co;
     
     for(int i=nparams;i--;) 
     {
          if(results[i]!=1)
          {
               printf("Input file is incomplete.\n");
               printf("\n");
               printf("Example:\n");
               printf("nsteps = 1000\n");
               printf("step = 10\n");
               printf("ncell = 400\n");
               printf("m2start = 200\n");
               printf("m2stop = 400\n");
               printf("dx = 0.01");
               printf("time_scale = 0.5");
               printf("\n");
               abort();
          }
     }
     
     fclose(f);
     
     /************* Print what was read ***************************************/
     printf("\n");
     printf("Input parameters\n");
     printf("nsteps = %d\n",this->nsteps);
     printf("step = %d\n",this->step);
     printf("ncell = %d\n",this->ncell);
     printf("m2start = %d\n",this->m2start);
     printf("m2stop = %d\n",this->m2stop);
     printf("dx = %f\n",this->dx);
     printf("dt = %e (time scale %f)\n",this->dt,this->time_scale);
     printf("\n");
}

/******************************************************************************/
void IO::write_field_to_file(const int n, const double *ex, const double *hy)
{
     char filename[200];
     sprintf(filename,"./output/out_%.06d.dat",n);
     std::ofstream fp;
     fp.open(filename,std::ios::out);
     fp.precision(6);
     // Loop goes from 0 to KE-1 (excluded) because hy[KE-1] requires
     // ex[KE] that is undefined.
     for(int k=0; k < this->ncell-1; k++)
     {   // Typecast to (float) to avoid writing numbers with large exponents in the file
          fp << (float)k     << "\t" 
             << (float)ex[k] << "\t"
             << (float)hy[k] << std::endl;
     }      
     fp.close();
}
/****************** End of file ***********************************************/