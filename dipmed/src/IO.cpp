#include "IO.hpp"

#include <cstdio>
#include <cstdlib>

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
     this->output_file = this->output_dir + "input.txt";
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
     FILE *inF = fopen(this->input_file.c_str(),"r");
     FILE *ouF = fopen(this->output_file.c_str(),"w");
     if(inF == NULL) perror(this->input_file.c_str()), exit(-1);
     if(ouF == NULL) perror(this->output_file.c_str()), exit(-1);
     
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
     int results[this->nparams];
     results[0] = fscanf(f,"nsteps = %i\n", &this->nsteps);
     results[1] = fscanf(f,"step = %i\n", &this->step);
     results[2] = fscanf(f,"ncell = %i\n", &this->ncell);
     results[3] = fscanf(f,"m2start = %i\n", &this->m2start);
     results[4] = fscanf(f,"m2stop = %i\n", &this->m2stop);
     for(int i=this->nparams;i--;) 
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
     printf("\n");
}

/****************** End of file ***********************************************/