
#include "material.hpp"

// #include <cassert>
// #include <cmath>
// #include <cstdio>
// #include <cstdlib>

#include <omp.h>

// #include "constants.hpp"
// #include "system.hpp"

/****************** Constructor/Destructor ************************************/
material::material(const int _argc, const char **_argv):IO(_argc,_argv)
{

}

/******************************************************************************/
material::~material()
{

}
          
/****************** Member functions ******************************************/
double material::static_response()
{
    return epsilon;
}

/******************************************************************************/

/****************** End of file ***********************************************/