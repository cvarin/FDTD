# include <math.h>
# include <stdlib.h>
# include <stdio.h>

/******************************************************************************/
int main(int argc, char *argv[])
{
     /************* Simulation parameters *************************************/
     const int IE = 140;
     const int JE = 140;
     const int npml = 8;
     
     const int ic = IE/2-20;
     const int jc = JE/2-20;
//      const double ddx = .01; /* Cell size */
//      const double dt = ddx/6e8; /* Time steps */

     const double t0 = 40.0;
     const double spread = 12.0;
     
     /************* Variables *************************************************/
     double Dz[IE][JE];
     double Ez[IE][JE];
     double Hx[IE][JE];
     double Hy[IE][JE];
     double ga[IE][JE];
     double gi2[IE],gi3[IE];
     double gj2[JE],gj3[IE];
     double fi1[IE],fi2[IE],fi3[JE];
     double fj1[JE],fj2[JE],fj3[JE];
     double ihx[IE][JE],ihy[IE][JE];

     /************* Initialize the arrays *************************************/
     for(int j=0; j < JE; j++ ) 
     {
//           printf( "%2d ",j);
          for(int i=0; i < IE; i++ ) 
          {
               Ez[i][j] = 0.0;
               Dz[i][j] = 0.0;
               Hx[i][j] = 0.0;
               Hy[i][j] = 0.0;
               ihx[i][j] = 0.0;
               ihy[i][j] = 0.0;
               ga[i][j] = 1.0;
//                printf( "%5.2f ",ga[i][j]);
          }
//           printf( " \n");
     }

     /*************************************************************************/
     /************* Calculate the PML parameters ******************************/
     for(int i=0;i< IE; i++) 
     {
          gi2[i] = 1.0;
          gi3[i] = 1.0;
          fi1[i] = 0.0;
          fi2[i] = 1.0;
          fi3[i] = 1.0;
     }
     
     /*************************************************************************/
     for(int j=0;j< JE; j++) 
     {
          gj2[j] = 1.0;
          gj3[j] = 1.0;
          fj1[j] = 0.0;
          fj2[j] = 1.0;
          fj3[j] = 1.0;
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
          /* pulse = sin(2*pi*1500*1e6*dt*T); */
          double pulse = exp(-.5*pow( (T-t0)/spread,2.));
          Dz[ic][jc] = pulse;
//                Dz[ic + 50][jc + 50] = pulse;

          /******** Calculate the Ez field ************************************/
          for(int j=1; j < JE-1; j++ )
          {
               for(int i=1; i < IE-1; i++ )
                    Ez[i][j] = ga[i][j]*Dz[i][j];
          }
     //      printf("%3f %6.2f \n ",T,ez[ic][jc]);
     
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

          /******** Calculate the Hx field ************************************/
          for(int j=0; j < JE-1; j++ ) 
          {
               for(int i=0; i < IE; i++ ) 
               {
                    double curl_e = Ez[i][j] - Ez[i][j+1] ;
                    ihx[i][j] = ihx[i][j] + fi1[i]*curl_e ;
                    Hx[i][j] = fj3[j]*Hx[i][j]
                    + fj2[j]*.5*(curl_e + ihx[i][j]);
               }
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
     }

     /************* Write the E field out to a file ***************************/
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

     /************* Print message *********************************************/
     printf("T = %d\n ",T);
}

/****************** End of file ***********************************************/
