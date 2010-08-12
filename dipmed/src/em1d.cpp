
#include "em1d.hpp"

#include <cassert>
#include <cstdio>
#include <cstdlib>
// #include <fstream>
// 
// #include "constants.hpp"

/****************** Constructor/Destructor ************************************/
em1d::em1d(int _argc, char **_argv):IO(_argc,_argv)
{
    bytes_allocated = 0;
  
    /**************************************************************************/
    // Initialize the E field and all cells to free space
    ex = (double *)calloc(ncell,sizeof(double)); assert(ex);
    hy = (double *)calloc(ncell,sizeof(double)); assert(hy);
    cb = (double *)calloc(ncell,sizeof(double)); assert(cb);
    
    /**************************************************************************/
    if(ex==NULL) printf("Can't allocate ex\n"),abort();
    if(hy==NULL) perror("Can't allocate hy\n"),abort();
    if(cb==NULL) perror("Can't allocate cb\n"),abort();
  
    /**************************************************************************/
    bytes_allocated += 3*ncell*sizeof(double);
    print_allocated_memory_in_Kbytes();
    
    /**************************************************************************/
    // Set free space everywhere
    for(int k=0; k <= ncell-1; k++) cb[k] = 1.0; 
        
    /*************************************************************************/
    // Initialize the medium 2        
    for(int k=m2start; k < m2stop; k++) cb[k] = 1.0/epsilon;
}

/******************************************************************************/
em1d::~em1d()
{
    free(cb);
    free(hy);
    free(ex);
    bytes_allocated -= 3*ncell*sizeof(double);
    print_allocated_memory_in_Kbytes();
}
          
/****************** Member functions ******************************************/
void em1d::print_allocated_memory_in_Kbytes()
{
    printf("Memory allocated: %f KB\n",(double)bytes_allocated/1024);
}

/******************************************************************************/
void em1d::print_allocated_memory_in_Mbytes()
{
    printf("Memory allocated: %f MB\n",(double)bytes_allocated/1024/1024);
}

/******************************************************************************/

/****************** End of file ***********************************************/