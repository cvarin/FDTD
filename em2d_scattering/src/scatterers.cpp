#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/******************************************************************************/
int main(int argc, char *argv[])
{          
     /************* Simulation parameters *************************************/
     const int IE = 140;
     const int JE = 140;
     const int npml = 8;
     const int outper = 10;
          
     const int ic = (int) (IE/2 - 10);
     const int jc = (int) (JE/2 - 10);
//      const double ddx = .01; /* Cell size */
//      const double dt = ddx/6e8; /* Time steps */

     const double t0 = 40.0;
     const double spread = 12.0;
     
     /*************************************************************************/
     /************* Field arrays **********************************************/
     /*************************************************************************/
     double Dz[IE][JE];
     double Ez[IE][JE];
     double Hx[IE][JE];
     double Hy[IE][JE];
     double ihx[IE][JE];
     double ihy[IE][JE];
     double ga[IE][JE];

     /************* Initialize the arrays *************************************/
     for(int j=0; j < JE; j++ ) 
          for(int i=0; i < IE; i++ ) 
          {
               Ez[i][j] = 0.0;
               Dz[i][j] = 0.0;
               Hx[i][j] = 0.0;
               Hy[i][j] = 0.0;
               ihx[i][j] = 0.0;
               ihy[i][j] = 0.0;
               ga[i][j] = 1.0;
          }

     /*************************************************************************/
     /************* Calculate the PML parameters ******************************/     
     /*************************************************************************/
     double gi2[IE];
     double gi3[IE];
     double fi1[IE];
     double fi2[IE];
     double fi3[IE];
     for(int i=0;i< IE; i++) 
     {
          gi2[i] = 1.0;
          gi3[i] = 1.0;
          fi1[i] = 0.0;
          fi2[i] = 1.0;
          fi3[i] = 1.0;
     }
     
     /*************************************************************************/
     for(int i=0;i<= npml; i++) 
     {
          double xnum = npml - i;
          double xd = npml;
          double xxn = xnum/xd;
          double xn = 0.25*pow(xxn,3.0);
//           printf(" %d %7.4f %7.4f \n",i,xxn,xn);
          gi2[i] = 1.0/(1.0+xn);
          gi2[IE-1-i] = 1.0/(1.0+xn);
          gi3[i] = (1.0 - xn)/(1.0 + xn);
          gi3[IE-i-1] = (1.0 - xn)/(1.0 + xn);
          xxn = (xnum-.5)/xd;
          xn = 0.25*pow(xxn,3.0);
          fi1[i] = xn;
          fi1[IE-2-i] = xn;
          fi2[i] = 1.0/(1.0+xn);
          fi2[IE-2-i] = 1.0/(1.0+xn);
          fi3[i] = (1.0 - xn)/(1.0 + xn);
          fi3[IE-2-i] = (1.0 - xn)/(1.0 + xn);
     }
     
     /*************************************************************************/
     double gj2[JE];
     double gj3[JE];
     double fj1[JE];
     double fj2[JE];
     double fj3[JE];
     for(int j=0;j< JE; j++) 
     {
          gj2[j] = 1.0;
          gj3[j] = 1.0;
          fj1[j] = 0.0;
          fj2[j] = 1.0;
          fj3[j] = 1.0;
     }

     /*************************************************************************/
     for(int j=0;j<= npml; j++) 
     {
          double xnum = npml - j;
          double xd = npml;
          double xxn = xnum/xd;
          double xn = 0.25*pow(xxn,3.0);
//           printf(" %d %7.4f %7.4f \n",i,xxn,xn);
          gj2[j] = 1.0/(1.0+xn);
          gj2[JE-1-j] = 1.0/(1.0+xn);
          gj3[j] = (1.0 - xn)/(1.0 + xn);
          gj3[JE-j-1] = (1.0 - xn)/(1.0 + xn);
          xxn = (xnum-.5)/xd;
          xn = 0.25*pow(xxn,3.0);
          fj1[j] = xn;
          fj1[JE-2-j] = xn;
          fj2[j] = 1.0/(1.0+xn);
          fj2[JE-2-j] = 1.0/(1.0+xn);
          fj3[j] = (1.0 - xn)/(1.0 + xn);
          fj3[JE-2-j] = (1.0 - xn)/(1.0 + xn);
     }

//      printf("gi + fi \n");
//      for (i=0; i< IE; i++) 
//      {
//           printf( "%2d %5.2f %5.2f \n",
//           i,gi2[i],gi3[i]),
//           printf( " %5.2f %5.2f %5.2f \n ",
//           fi1[i],fi2[i],fi3[i]);
//      }
// 
//      printf("gj + fj \n");
//      for (j=0; j< JE; j++) 
//      {
//           printf( "%2d %5.2f %5.2f \n",
//           j,gj2[j],gj3[j]),
//           printf( " %5.2f %5.2f %5.2f \n ",
//           fj1[j],fj2[j],fj3[j]);
//      }

     /*************************************************************************/
     /************* Arrays for the incident wave source ***********************/     
     /*************************************************************************/
     double Ez_inc[JE];
     double Hx_inc[JE];
     for(int j=0; j < JE; j++)
     {
          Ez_inc[j] = 0.0;
          Hx_inc[j] = 0.0;
     }
     
     double Ez_inc_low_m1 = 0.0;
     double Ez_inc_low_m2 = 0.0;
     double Ez_inc_high_m1 = 0.0;
     double Ez_inc_high_m2 = 0.0;
     
     const int ntfsf = npml + 1;
     const int ia = ntfsf;
     const int ib = IE - ia - 1;
     const int ja = ntfsf;
     const int jb = JE - ja - 1;
     
     /*************************************************************************/
     int T = 0;
     int nsteps = 1;
     
     /************* Ask for the number of steps to calculate ******************/
     int result = 0;
     while(1)
     {
          printf( "nsteps --> ");
          result = scanf("%d", &nsteps);
          if(result == 1) break;
          else
          {
               printf("nsteps should be an integer.\n");
               abort();
          }
     }

     /************* Loop over all time steps **********************************/
     for(int n=1; n <=nsteps ; n++) 
     {
          T += 1;

          /******** Update Ez_inc *********************************************/
          for(int j=1; j < JE; j++)
               Ez_inc[j] += 0.5*(Hx_inc[j-1] - Hx_inc[j]);
          
          Ez_inc[0]     = Ez_inc_low_m2;
          Ez_inc_low_m2 = Ez_inc_low_m1;
          Ez_inc_low_m1 = Ez_inc[1];
          
          Ez_inc[JE - 1] = Ez_inc_high_m2;
          Ez_inc_high_m2 = Ez_inc_high_m1;
          Ez_inc_high_m1 = Ez_inc[JE - 2];
          
          /******** Calculate the Dz field ************************************/
          for(int j=1; j < IE; j++ ) 
          {
               for(int i=1; i < IE; i++ ) 
               {
                    Dz[i][j] = gi3[i]*gj3[j]*Dz[i][j]
                    + gi2[i]*gj2[j]*.5*(Hy[i][j] - Hy[i-1][j]
                    - Hx[i][j] + Hx[i][j-1]) ;
               }
          }

          /******** Source ****************************************************/
          double pulse = sin(2.0*3.14159*1.5e9*T);
//           double pulse = exp(-0.5*pow( (T-t0)/spread,2.0));
//           Dz[ic][jc] = 5.0*pulse;
//           Dz[ia+20][jc] = 5.0*pulse;
//           Dz[ic + 50][jc + 50] = pulse;

          /******** Ez TFSF ***************************************************/
          Ez_inc[3] = pulse;
          for(int i = ia; i<= ib; i++)
          {
               Dz[i][ja] = Dz[i][ja] + 0.5*Hx_inc[ja - 1];
               Dz[i][jb] = Dz[i][jb] - 0.5*Hx_inc[jb];
          }
          
          /******** Calculate the Ez field ************************************/
          for(int j=1; j < JE-1; j++ )
          {
               for(int i=1; i < IE-1; i++ )
                    Ez[i][j] = ga[i][j]*Dz[i][j];
          }
     
          /******** Set the PEC at the simulation boundaries ******************/
          for(int j=0; j < JE-1; j++)
          {
               Ez[0][j] = 0.0;
               Ez[IE-1][j] = 0.0;
          }
          for(int i=0; i < IE-1; i++)
          {
               Ez[i][0] = 0.0;
               Ez[i][JE-1] = 0.0;
          }

          /******** Update Hx source ******************************************/
          for(int j = 0; j < JE-1; j++)
               Hx_inc[j] += 0.5*(Ez_inc[j] - Ez_inc[j + 1]);

          /******** Calculate the Hx field ************************************/
          for(int j=0; j < JE-1; j++) 
          {
               for(int i=0; i < IE; i++) 
               {
                    double curl_e = Ez[i][j] - Ez[i][j+1] ;
                    ihx[i][j] = ihx[i][j] + fi1[i]*curl_e ;
                    Hx[i][j] = fj3[j]*Hx[i][j]
                    + fj2[j]*.5*(curl_e + ihx[i][j]);
               }
          }

          /******** Incident Hx values ****************************************/
          for(int i = ia; i <= ib; i++)
          {
               Hx[i][ja-1] += 0.5*Ez_inc[ja];
               Hx[i][jb]   -= 0.5*Ez_inc[jb];
          }

          /******** Calculate the Hy field ************************************/
          for(int j=0; j <= JE-1; j++ ) 
          {
               for(int i=0; i < IE-1; i++ ) 
               {
                    double curl_e = Ez[i+1][j] - Ez[i][j];
                    ihy[i][j] = ihy[i][j] + fj1[j]*curl_e ;
                    Hy[i][j] = fi3[i]*Hy[i][j]
                    + fi2[i]*.5*(curl_e + ihy[i][j]);
               }
          }
          
          /******** Incident Hy values ****************************************/
          for(int j = ja; j <= jb; j++)
          {
               Hy[ia - 1][j] -= 0.5*Ez_inc[j];
               Hy[ib][j]     += 0.5*Ez_inc[j];
          }
          
          /******** Write Ez to file ******************************************/
          if(n%outper==0)
          {
               char filename[200];
               sprintf(filename,"output/Ez_%06d.dat",T);
               FILE *fp = fopen(filename,"w");
               for(int j=0; j < JE; j++ ) 
               {
                    for(int i=0; i < IE; i++ ) 
                    {
                         fprintf( fp,"%6.3f ",Ez[i][j]);
                    }
                    fprintf(fp," \n");
               }
               fclose(fp);
          }
     }

     /************* Write the E field out to a file ***************************/
//      char filename[200];
//      sprintf(filename,"output/Ez_%06d.dat",T);
//      FILE *fp = fopen(filename,"w");
//      for(int j=0; j < JE; j++ ) 
//      {
//           for(int i=0; i < IE; i++ ) 
//           {
//                fprintf( fp,"%6.3f ",Ez[i][j]);
//           }
//           fprintf(fp," \n");
//      }
//      fclose(fp);

     /************* Print message *********************************************/
     printf("T = %d\n ",T);
}

/****************** End of file ***********************************************/
