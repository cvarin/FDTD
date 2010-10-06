# include <math.h>
# include <stdlib.h>
# include <stdio.h>

main ()
{
     const int IE = 140;
     const int JE = 140;
     
     const int npml = 8;
     
     float ga[IE][JE],dz[IE][JE],ez[IE][JE],hx[IE][JE],hy[IE][JE];
     int l,n,i,j,ic,jc,nsteps;
     float ddx,dt,T,epsz,pi,epsilon,sigma,eaf;
     float xn,xxn,xnum,xd,curl_e;
     float t0,spread,pulse;
     float gi2[IE],gi3[IE];
     float gj2[JE],gj3[IE];
     float fi1[IE],fi2[IE],fi3[JE];
     float fj1[JE],fj2[JE],fj3[JE];
     float ihx[IE][JE],ihy[IE][JE];
     FILE *fp;

     ic = IE/2-20;
     jc = JE/2-20;
     ddx = .01; /* Cell size */
     dt =ddx/6e8; /* Time steps */
     epsz = 8.8e-12;
     pi=3.14159;

     /* Initialize the arrays */
     for ( j=0; j < JE; j++ ) 
     {
          printf( "%2d ",j);
          for ( i=0; i < IE; i++ ) 
          {
               dz[i][j]= 0.0 ;
               hx[i][j]= 0.0 ;
               hy[i][j]= 0.0 ;
               ihx[i][j]= 0.0 ;
               ihy[i][j]= 0.0 ;
               ga[i][j]= 1.0 ;
               printf( "%5.2f ",ga[i][j]);
          }
          printf( " \n");
     }

     /* Calculate the PML parameters */
     for (i=0;i< IE; i++) 
     {
          gi2[i] = 1.0;
          gi3[i] = 1.0;
          fi1[i] = 0.0;
          fi2[i] = 1.0;
          fi3[i] = 1.0;
     }
     for (j=0;j< IE; j++) 
     {
          gj2[j] = 1.0;
          gj3[j] = 1.0;
          fj1[j] = 0.0;
          fj2[j] = 1.0;
          fj3[j] = 1.0;
     }

//      printf( "Number of PML cells --> ");
//      scanf("%d", &npml);

     for (i=0;i<= npml; i++) 
     {
          xnum = npml - i;
          xd = npml;
          xxn = xnum/xd;
          xn = 0.25*pow(xxn,3.0);
          printf(" %d %7.4f %7.4f \n",i,xxn,xn);
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

     for (j=0;j<= npml; j++) 
     {
          xnum = npml - j;
          xd = npml;
          xxn = xnum/xd;
          xn = 0.25*pow(xxn,3.0);
          printf(" %d %7.4f %7.4f \n",i,xxn,xn);
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

     printf("gi + fi \n");
     for (i=0; i< IE; i++) 
     {
          printf( "%2d %5.2f %5.2f \n",
          i,gi2[i],gi3[i]),
          printf( " %5.2f %5.2f %5.2f \n ",
          fi1[i],fi2[i],fi3[i]);
     }

     printf("gj + fj \n");
     for (j=0; j< JE; j++) 
     {
          printf( "%2d %5.2f %5.2f \n",
          j,gj2[j],gj3[j]),
          printf( " %5.2f %5.2f %5.2f \n ",
          fj1[j],fj2[j],fj3[j]);
     }

     t0 = 40.0;
     spread = 12.0;
     T = 0;
     nsteps = 1;

     while ( nsteps > 0 ) 
     {
          printf( "nsteps --> ");
          scanf("%d", &nsteps);
          printf("%d \n", nsteps);

          for ( n=1; n <=nsteps ; n++) 
          {
               T = T + 1;

               /* ---- Start of the Main FDTD loop ---- */

               /* Calculate the Dz field */
               for ( j=1; j < IE; j++ ) {
               for ( i=1; i < IE; i++ ) {
               dz[i][j] = gi3[i]*gj3[j]*dz[i][j]
               + gi2[i]*gj2[j]*.5*( hy[i][j] - hy[i-1][j]
               - hx[i][j] + hx[i][j-1]) ;
          }
     }

     /* Sinusoidal Source */

     /* pulse = sin(2*pi*1500*1e6*dt*T);; */
     pulse = exp(-.5*pow( (T-t0)/spread,2.));
     dz[ic][jc] = pulse;

     /* Calculate the Ez field */
     /* Leave the Ez edges to 0, as part of the PML */
     for ( j=1; j < JE-1; j++ ) 
     {
          for ( i=1; i < IE-1; i++ ) 
          {
               ez[i][j] = ga[i][j]*dz[i][j] ;
          }
     }

     printf("%3f %6.2f \n ",T,ez[ic][jc]);

     /* Calculate the Hx field */
     for ( j=0; j < JE-1; j++ ) 
     {
          for ( i=0; i < IE; i++ ) 
          {
               curl_e = ez[i][j] - ez[i][j+1] ;
               ihx[i][j] = ihx[i][j] + fi1[i]*curl_e ;
               hx[i][j] = fj3[j]*hx[i][j]
               + fj2[j]*.5*(curl_e + ihx[i][j]);
          }
     }

     /* Calculate the Hy field */
     for ( j=0; j <= JE-1; j++ ) 
     {
          for ( i=0; i < IE-1; i++ ) 
          {
               curl_e = ez[i+1][j] - ez[i][j];
               ihy[i][j] = ihy[i][j] + fj1[j]*curl_e ;
               hy[i][j] = fi3[i]*hy[i][j]
               + fi2[i]*.5*(curl_e + ihy[i][j]);
          }
     }

     }
     /* ---- End of the main FDTD loop ---- */

     for ( j=1; j < JE; j++ ) 
     {
          printf( "%2d ",j);
          for ( i=1; i <= IE; i++ ) 
          {
               printf( "%4.1f",ez[i][j]);
          }
          printf( " \n");
     }

     /* Write the E field out to a file "Ez" */
     fp = fopen( "Ez","w");
     for ( j=0; j < JE; j++ ) 
     {
          for ( i=0; i < IE; i++ ) 
          {
               fprintf( fp,"%6.3f ",ez[i][j]);
          }
          fprintf( fp," \n");
     }

     fclose(fp);

     printf("T = %6.0f \n ",T);

     }
}

/****************** End of file ***********************************************/
