/*  Fd2d_3.4.c. 2D TM simulation of a plane
 wave source impinging on a dielectric cylinder.
 Analysis using Fourier Transforms */
# include <math.h>
# include <stdlib.h>
# include <stdio.h>

#define IE 50
#define JE 50
#define NFREQS 3

main ()
{
    float dz[IE][JE],ez[IE][JE],hx[IE][JE],hy[IE][JE];
    float ga[IE][JE],gb[IE][JE],iz[IE][JE];
    int l,m,n,i,j,ic,jc,nsteps,npml;
    float ddx,dt,T,epsz,pi,epsilon,sigma;
    float xn,xxn,xnum,xd,curl_e;
    float t0,spread,pulse;
    float gi1[IE],gi2[IE],gi3[IE];
    float gj1[JE],gj2[JE],gj3[IE];
    float fi1[IE],fi2[IE],fi3[JE];
    float fj1[JE],fj2[JE],fj3[JE];
    float ihx[IE][JE],ihy[IE][JE];
    FILE *fp, *fopen();
    float ez_inc[JE],hx_inc[JE];
    float ez_inc_low_m1,ez_inc_low_m2;
    float ez_inc_high_m1,ez_inc_high_m2;
    int ia,ib,ja,jb;
    float radius,xdist,ydist,dist;
    float freq[NFREQS],arg[NFREQS];
    float
real_pt[NFREQS][IE][JE],imag_pt[NFREQS][IE][JE],amp[IE][JE],phase[IE][JE];
    float real_in[5],imag_in[5],amp_in[5],phase_in[5];

    ic = (IE/2)-1;
    jc = (JE/2)-1;
    ia = 10;                   /* Total/scattered field boundaries */
    ib = IE-ia-1;
    ja = 10;
    jb = JE-ja-1;
    ddx = .01;                 /* Cell size */
    dt =ddx/6e8;               /* Time step */
    epsz = 8.8e-12;
    pi=3.14159;

    /* Initialize the arrays */
    for ( j=0; j < JE; j++ ) {
       printf( "%2d  ",j);
	   ez_inc[j] = 0.;
	   hx_inc[j] = 0.;
       for ( i=0; i < IE; i++ ) {
           dz[i][j]= 0.0  ;
           ez[i][j]= 0.0  ;
           hx[i][j]= 0.0  ;
           hy[i][j]= 0.0  ;
           iz[i][j]= 0.0  ;
           ihx[i][j]= 0.0  ;
           ihy[i][j]= 0.0  ;
           ga[i][j]= 1.0  ;
           gb[i][j]= 0.0  ;
           printf( "%5.2f ",ga[i][j]);
       }
       printf( " \n");
    }

    /* Parameters for the Fourier Transforms */

    freq[0] =  50.e6;
    freq[1] = 300.e6;
    freq[2] = 700.e6;

   for ( n=0; n < NFREQS; n++)
   { arg[n] = 2*pi*freq[n]*dt;
   printf( "%d %6.2f %7.5f \n",n,freq[n]*1e-6,arg[n]);
   }

    /* Specify the dielectric cylinder */

       printf( "Cylinder radius (cells), epsilon, sigma --> ");
       scanf("%f %f %f", &radius, &epsilon, &sigma);
       printf( "Radius = %5.2f  Eps = %6.2f Sigma = %6.2f \n ",
	 radius,epsilon,sigma);

       for ( j = ja; j < jb; j++ ) {
          for ( i=ia; i < ib; i++ ) {
	      xdist = (ic-i);
	      ydist = (jc-j);
	      dist = sqrt(pow(xdist,2.) + pow(ydist,2.));
		if( dist <= radius) {
		   ga[i][j] = 1./(epsilon + (sigma*dt/epsz));
		   gb[i][j] = sigma*dt/epsz;
                }
           }
        }

      /* Write the ga field out to a file "Ga" */
       fp = fopen( "Ga","w");
       printf( " Ga  \n");
       for ( j=ja; j <= jb; j++ ) {
          for ( i=ia; i <= ib; i++ ) {
             printf( "%5.2f",ga[i][j]);
             fprintf( fp,"%6.3f ",ga[i][j]);
	  }
             printf( " \n");
             fprintf( fp," \n");
       }

       printf( " Gb  \n");
       for ( j=ja; j < jb; j++ ) {
          for ( i=ia; i <= ib; i++ ) {
             printf( "%5.2f",gb[i][j]);
	  }
             printf( " \n");
       }
        fclose(fp);

    /* Calculate the PML parameters */

       for (i=0;i< IE; i++) {
	  gi1[i] = 0.0;
	  gi2[i] = 1.0;
	  gi3[i] = 1.0;
	  fi1[i] = 0.0;
	  fi2[i] = 1.0;
	  fi3[i] = 1.0;
       }
       for (j=0;j< IE; j++) {
	  gj1[j] = 0.0;
	  gj2[j] = 1.0;
	  gj3[j] = 1.0;
	  fj1[j] = 0.0;
	  fj2[j] = 1.0;
	  fj3[j] = 1.0;
       }

       printf( "Number of PML cells --> ");
       scanf("%d", &npml);

       for (i=0;i<= npml; i++) {
       xnum  = npml - i;
       xd  = npml;
       xxn = xnum/xd;
       xn  = 0.25*pow(xxn,3.0);
       printf(" %d %7.4f  %7.4f \n",i,xxn,xn);
          gi1[i] = xn;
          gi1[IE-1-i] = xn;
          gi2[i] = 1.0/(1.0+xn);
          gi2[IE-1-i] = 1.0/(1.0+xn);
          gi3[i] = (1.0 - xn)/(1.0 + xn);
          gi3[IE-i-1] = (1.0 - xn)/(1.0 + xn);
       xxn = (xnum-.5)/xd;
       xn  = 0.25*pow(xxn,3.0);
          fi1[i] = xn;
          fi1[IE-2-i] = xn;
          fi2[i] = 1.0/(1.0+xn);
          fi2[IE-2-i] = 1.0/(1.0+xn);
          fi3[i] = (1.0 - xn)/(1.0 + xn);
          fi3[IE-2-i] = (1.0 - xn)/(1.0 + xn);
       }

       for (j=0;j<= npml; j++) {
       xnum  = npml - j;
       xd  = npml;
       xxn = xnum/xd;
       xn  = 0.25*pow(xxn,3.0);
       printf(" %d %7.4f  %7.4f \n",i,xxn,xn);
          gj1[j] = xn;
          gj1[JE-1-j] = xn;
          gj2[j] = 1.0/(1.0+xn);
          gj2[JE-1-j] = 1.0/(1.0+xn);
          gj3[j] = (1.0 - xn)/(1.0 + xn);
          gj3[JE-j-1] = (1.0 - xn)/(1.0 + xn);
       xxn = (xnum-.5)/xd;
       xn  = 0.25*pow(xxn,3.0);
          fj1[j] = xn;
          fj1[JE-2-j] = xn;
          fj2[j] = 1.0/(1.0+xn);
          fj2[JE-2-j] = 1.0/(1.0+xn);
          fj3[j] = (1.0 - xn)/(1.0 + xn);
          fj3[JE-2-j] = (1.0 - xn)/(1.0 + xn);
       }

       printf("gi + fi \n");
       for (i=0; i< IE; i++) {
           printf( "%2d  %5.2f  %5.2f  %5.2f \n",
	            i,gi1[i],gi2[i],gi3[i]),
           printf( "     %5.2f  %5.2f  %5.2f \n ",
		       fi1[i],fi2[i],fi3[i]);
       }

       printf("gj + fj \n");
       for (j=0; j< JE; j++) {
           printf( "%2d  %5.2f  %5.2f  %5.2f \n",
                    j,gj1[j],gj2[j],gj3[j]),
           printf( "     %5.2f  %5.2f  %5.2f \n ",
		       fj1[j],fj2[j],fj3[j]);
       }

    t0 = 25.0;
    spread = 8.0;
    T = 0;
    nsteps = 1;

while ( nsteps > 0 ) {
       printf( "nsteps --> ");
       scanf("%d", &nsteps);
       printf("%d \n", nsteps);

    for ( n=1; n <=nsteps ; n++)  {
       T = T + 1;

/*  ----   Start of the Main FDTD loop ----  */

      /* Calculate the incidnet Ez  */
      for (j=1; j< JE; j++) {
         ez_inc[j] = ez_inc[j] + .5*(hx_inc[j-1]-hx_inc[j]);
      }

    /* Fourier Tramsform of the incident field */
        for ( m=0; m < NFREQS ; m++ )
        {  real_in[m] = real_in[m] + cos(arg[m]*T)*ez_inc[ja-1] ;
	   imag_in[m] = imag_in[m] - sin(arg[m]*T)*ez_inc[ja-1] ;
	}

      /* ABC for the incident buffer */
      ez_inc[0]      = ez_inc_low_m2;
      ez_inc_low_m2  = ez_inc_low_m1;
      ez_inc_low_m1  = ez_inc[1];

      ez_inc[JE-1]    = ez_inc_high_m2;
      ez_inc_high_m2  = ez_inc_high_m1;
      ez_inc_high_m1  = ez_inc[JE-2];

       /* Calculate the Dz field */
       for ( j=1; j < IE; j++ ) {
          for ( i=1; i < IE; i++ ) {
              dz[i][j] = gi3[i]*gj3[j]*dz[i][j]
	     + gi2[i]*gj2[j]*.5*( hy[i][j] - hy[i-1][j]
			        - hx[i][j] + hx[i][j-1]) ;
          }
       }

       /*  Source */

       /* pulse =  sin(2*pi*400*1e6*dt*T) ; */
       pulse =  exp(-.5*(pow((t0-T)/spread,2.0) ));
         ez_inc[3] = pulse;
  /*     dz[ic-5][jc-5] = pulse; */

       printf("%3.0f  %6.2f \n ",T,ez_inc[3]);
       /* Incident Dz values */

       for (i=ia; i<= ib; i++ )  {
          dz[i][ja] = dz[i][ja] + 0.5*hx_inc[ja-1];
          dz[i][jb] = dz[i][jb] - 0.5*hx_inc[jb];
       }

       /* Calculate the Ez field */
       for ( j=1; j < JE-1; j++ ) {
          for ( i=1; i < IE-1; i++ ) {
              ez[i][j] = ga[i][j]*( dz[i][j] - iz[i][j] ) ;
	      iz[i][j] = iz[i][j] + gb[i][j]*ez[i][j] ;
          }
       }

	      /* Calculate the Fourier transform of Ex. */
     for ( j=0; j < JE; j++ )
     {   for ( i=0; i < JE; i++ )
        {   for ( m=0; m < NFREQS; m++ )
           { real_pt[m][i][j] = real_pt[m][i][j] + cos(arg[m]*T)*ez[i][j] ;
             imag_pt[m][i][j] = imag_pt[m][i][j] + sin(arg[m]*T)*ez[i][j] ;
           }
        }
     }

      /* Calculate the incident Hx */
      for (j=0; j< JE; j++) {
         hx_inc[j] = hx_inc[j] + .5*(ez_inc[j]-ez_inc[j+1]);
      }

       /* Calculate the Hx field */
       for ( j=0; j < JE-1; j++ ) {
          for ( i=0; i < IE; i++ ) {
	     curl_e =   ez[i][j] - ez[i][j+1]  ;
	     ihx[i][j] = ihx[i][j]  + fi1[i]*curl_e ;
	     hx[i][j] = fj3[j]*hx[i][j]
		      + fj2[j]*.5*(curl_e + ihx[i][j]);
          }
       }

       /* Incident Hx values */

       for (i=ia; i<= ib; i++ )  {
          hx[i][ja-1] = hx[i][ja-1] + .5*ez_inc[ja];
          hx[i][jb]   = hx[i][jb] - .5*ez_inc[jb];
       }

       /* Calculate the Hy field */
       for ( j=0; j <= JE-1; j++ ) {
          for ( i=0; i < IE-1; i++ ) {
	     curl_e  = ez[i+1][j] - ez[i][j];
	     ihy[i][j] = ihy[i][j]  + fj1[j]*curl_e ;
	     hy[i][j] = fi3[i]*hy[i][j]
		      + fi2[i]*.5*(curl_e + ihy[i][j]);
          }
       }

       /* Incident Hy values */
       for (j=ja; j<= jb; j++ )  {
          hy[ia-1][j] = hy[ia-1][j] - .5*ez_inc[j];
          hy[ib][j]   = hy[ib][j]   + .5*ez_inc[j];
       }

    }
/*  ----  End of the main FDTD loop ---- */

  /* Calculate the Fouier amplitude and phase of the incident pulse */
       for ( m=0; m < NFREQS; m++ )
       {    amp_in[m]   = sqrt( pow(real_in[m],2.) + pow(imag_in[m],2.));
            phase_in[m] = atan2(imag_in[m],real_in[m]) ;
	     printf( "%d  Input Pulse :  %8.4f %8.4f %8.4f  %7.2f\n",
	     m,real_in[m],imag_in[m],amp_in[m],(180.0/pi)*phase_in[m]);
       }

   /* Calculate the Fouier amplitude and phase of the total field field */
       for ( m=0; m < NFREQS; m++ )
       {
       if( m == 0)       fp = fopen( "amp1","w");
       else if( m == 1)  fp = fopen( "amp2","w");
       else if( m == 2)  fp = fopen( "amp3","w");
       {   printf( "%2d  %7.2f  MHz\n",m,freq[m]*1.e-6);
           for ( j=ja; j <= jb; j++ )
           {  if( ga[ic][j] < 1.00 )
              {  amp[ic][j]  = (1./amp_in[m])
                 *sqrt( pow(real_pt[m][ic][j],2.) + pow(imag_pt[m][ic][j],2.));
                 printf( "%2d %9.4f \n",jc-j,amp[ic][j]);
                 fprintf( fp," %9.4f \n",amp[ic][j]);
               }
           }
        }
        fclose(fp);
	}

     /*  for ( j=1; j < JE; j++ ) {
          printf( "%2d  %5.2f",j,ez_inc[j]);
          for ( i=1; i <= IE; i++ ) {
             printf( "%5.2f",ez[i][j]);
	  }
             printf( " \n");
       } */

        /* Write the E field out to a file "Ez" */
       fp = fopen( "Ez","w");
       for ( j=0; j < jc; j++ ) {
	   for ( i=0; i < ic; i++ ) {
             fprintf( fp,"%6.3f ",ez[2*i][2*j]);
	  }
             fprintf( fp," \n");
        }
       printf("T = %3.0f  \n ",T);

        fclose(fp);

}
}
