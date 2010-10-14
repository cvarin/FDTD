#ifndef INC_em2d_hpp
#define INC_em2d_hpp

#include "system.hpp"

/******************************************************************************/
class em2d
{
     public:
     int IE;
     int JE;
     double **Dz;
     
     /************* Constructor, destructor ***********************************/
     em2d(const int _IE, const int _JE)
     {
          IE = _IE;
          JE = _JE;
          allocate_2D(&Dz,IE,JE);
          
          for(int j=0; j < JE; j++ ) 
               for(int i=0; i < IE; i++ ) 
                    Dz[i][j] = 0.0, printf("Dz[%d][%d] = %f\n",i,j,Dz[i][j]);
     }
     
     /*************************************************************************/
     ~em2d()
     {
          delete_2D(Dz,IE);
     }
};

#endif // #ifndef INC_em2d_hpp