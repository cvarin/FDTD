#include <cmath>
#include <cstdlib>
#include <cstring>
#include <omp.h>                        
#include <stdio.h>                      
#include <time.h>                       

/****************** Global declarations ***********************************************************/
inline void allocate_3D(double ****m, int d1, int d2, int d3) 
{                                                             
     *m = new double** [d1];                                  
     for (int i=0; i<d1; i++)                                 
     {                                                        
          (*m)[i] = new double* [d2];                         
          for (int j=0; j<d2; j++)                            
          {                                                   
              (*m)[i][j] = new double [d3];                   
              for (int k=0; k<d3; k++) (*m)[i][j][k] = 0.0;   
          }                                                   
     }                                                        
}

/**************************************************************************************************/
struct Grid{double Hx,Hy,Hz;double Ex,Ey,Ez;};
// struct Grid{double Ex,Ey,Ez;};
void Allocate_Grid(Grid **g, const unsigned int &size)
{
    *g = (Grid*)calloc(size,sizeof(Grid));
    if (*g == NULL) printf("Can't allocate grid!\n"),abort();

    // If not present, valgrind bails when saving snapshot for unknown reason.
    memset(*g,0,size*sizeof(Grid));
}

/**************************************************************************************************/
void perform_with_3D_vectors(const int &nsteps, const int N[3])
{
     /*********************************************************************************************/
     /************* Première approche... vecteurs 3D indépendants *********************************/
     /*********************************************************************************************/
     double ***Ex;
     double ***Ey;
     double ***Ez;
     double ***Hx;
     double ***Hy;
     double ***Hz;
     allocate_3D(&Ex,N[0],N[1],N[2]);
     allocate_3D(&Ey,N[0],N[1],N[2]);
     allocate_3D(&Ez,N[0],N[1],N[2]);
     allocate_3D(&Hx,N[0],N[1],N[2]);
     allocate_3D(&Hy,N[0],N[1],N[2]);
     allocate_3D(&Hz,N[0],N[1],N[2]);
          
     /************* Do the loop *******************************************************************/
     for(int n=0;n<nsteps;n++)
     {
          #pragma omp parallel for
          for (int i = 1; i < N[0]; i++)
          {
               for (int j = 1; j < N[1]; j++)
               {
                    for (int k = 1; k < N[2]; k++)
                    {
                         Ex[i][j][k] += (Hz[i][j][k] - Hz[i][j-1][k]) - (Hy[i][j][k] - Hy[i][j][k-1]);
                         Ey[i][j][k] += (Hx[i][j][k] - Hx[i][j][k-1]) - (Hz[i][j][k] - Hz[i-1][j][k]);
                         Ez[i][j][k] += (Hy[i][j][k] - Hy[i-1][j][k]) - (Hx[i][j][k] - Hx[i][j-1][k]);
//                          Ex[i][j][k] = exp(i*j*k);
//                          Ey[i][j][k] = exp(i*j*k);
//                          Ez[i][j][k] = exp(i*j*k);
                    }
               }
          }
     }
}

/**************************************************************************************************/
void perform_with_1D_structure(const int &nsteps, const int N[3])
{
     Grid *g;
     Allocate_Grid(&g,N[0]*N[1]*N[2]);
     const int N01 = N[0]*N[1];
     int I,J,K;
     
     /************* Do the loop *******************************************************************/
     for(int n=0;n<nsteps;n++)
     {
          #pragma omp parallel for private(I,J,K) shared(g)
          for (int k = 1; k < N[2]; k++)
          {
               K = k*N[0]*N[1];
               for (int j = 1; j < N[1]; j++)
               {   
                    J = j*N[0] + K;
                    I = 1 + J;
                    for (int i = 1; i < N[0]; i++, I++)
                    {
                         g[I].Ex += (g[I].Hz - g[I-N[0]].Hz) - (g[I].Hy - g[I-N01].Hy);
                         g[I].Ey += (g[I].Hx - g[I-N01].Hx) - (g[I].Hz - g[I-1].Hz);
                         g[I].Ez += (g[I].Hy - g[I-1].Hy) - (g[I].Hx - g[I-N[0]].Hx);
                                        
//                          g[I].Ex = exp(i*j*k);
//                          g[I].Ey = exp(i*j*k);
//                          g[I].Ez = exp(i*j*k);
                    }
               }
          }
     }
}

/****************** Main **************************************************************************/
int main(int argc, char *argv[])
{
     const int N[3] = {200,200,200};
     int nthreads = 1;
     int nsteps = 100;
     if(argc > 1) nsteps = atoi(argv[1]);
     
     /************* Get number of treads and print message ****************************************/
     #if defined (_OPENMP)
     #pragma omp parallel default(none) shared(nthreads)
     {
          nthreads = omp_get_num_threads();
     }
     #endif
     printf("Will perform %d steps ",nsteps);
     #if defined (_OPENMP)
     printf("with %d thread(s)",nthreads);
     #endif
     printf("\n");
     fflush(stdout);
     
     /************* Do it *************************************************************************/
     perform_with_3D_vectors(nsteps,N);
//      perform_with_1D_structure(nsteps,N);
          
     /************* Print the time it took ********************************************************/
     int time = (int)clock();
     printf("Execution time : %d clock cycles",time);
     printf("(%f s) ",(double)time/CLOCKS_PER_SEC);
     printf("by %d threads.\n", nthreads);
     printf("\n\t=> %f s",(double)time/CLOCKS_PER_SEC/nthreads);
     printf("\n");
}

/*************************************************************************************************/