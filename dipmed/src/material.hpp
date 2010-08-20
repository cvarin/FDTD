#ifndef INC_material_hpp
#define INC_material_hpp

#include "IO.hpp"

class material : public IO
{    
     public:
          /******** Members ***************************************************/

          /******** Member functions ******************************************/
          material(const int _argc, const char **_argv);
          ~material();
          double static_response();
};

#endif // INC_material_hpp