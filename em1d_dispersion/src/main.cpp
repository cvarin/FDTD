
#include "system.hpp"
#include "em1d.hpp"

/****************** Main ******************************************************/
int main(const int argc, const char **argv)
{   
     print_character('*',60);
     const time_t start = time(NULL);
     em1d f(argc,argv);

     /*************************************************************************/
     for(int n=1; n <= f.nsteps; n++) f.advance_a_step(n);
     
     /*************************************************************************/
     ShowRunTime(start,time(NULL));
     return EXIT_SUCCESS;
}

/****************** End of file ***********************************************/
