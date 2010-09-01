/* FD1D_1.3.c.  1D FDTD simulation of a dielectric slab */
/* Simulation of a pulse hitting a dielectric slab */

# include <math.h>
# include <stdlib.h>
# include <stdio.h>

#define KE  200

int main (void)
{
    float ex[KE],hy[KE];
    float cb[KE];
    int n,k,kc,ke,kstart,nsteps;
    float ddx,dt,T,epsz,epsilon,sigma,eaf;
    float t0,spread,pi,pulse;
    FILE *fp;
//     *fopen();
    float ex_low_m1,ex_low_m2,ex_high_m1,ex_high_m2;

    kc = KE/2;
    ddx = .01;                 /* Cells size */
    pi  = 3.14159;
    dt =ddx/6e8;               /* Time steps */
    epsz = 8.8e-12;

    for ( k=1; k <= KE; k++ ) {  /* Initialize to free space */
      cb[k] = .5;
    }
       printf( "Dielectric starts at --> ");
       scanf("%d", &kstart);
       printf( "Epsilon --> ");
       scanf("%f", &epsilon);
       printf("%d %6.2f   \n", kstart,epsilon);

    for ( k=kstart; k <= KE; k++ ) {
      cb[k] = .5/epsilon;
    }

       for ( k=1; k <= KE; k++ )
       { printf( "%2d   %4.2f\n",k,cb[k]); }

    /* These parameters specify the input pulse */
    t0 = 40.0;
    spread = 15.0;
    T = 0;
    nsteps = 1;

    /* Main part of the program */

while ( nsteps > 0 ) {
       printf( "nsteps --> ");
       scanf("%d", &nsteps);
       printf("%d \n", nsteps);

    for ( n=1; n <=nsteps ; n++)
    {
       T = T + 1;

       /* Calculate the Ex field */
       for ( k=0; k < KE; k++ )
       { ex[k] = ex[k] + cb[k]*( hy[k-1] - hy[k] ) ; }

       /* Put a Gaussian pulse at the low end */

       pulse =  exp(-.5*(pow((t0-T)/spread,2.0)));
       ex[5] =  ex[5] + pulse;
       printf( "%5.1f %6.2f %6.2f\n",T,pulse,ex[5]);

       /* Boundary conditions */
       ex[0]      = ex_low_m2;
       ex_low_m2  = ex_low_m1;
       ex_low_m1  = ex[1];

       ex[KE-1]      = ex_high_m2;
       ex_high_m2  = ex_high_m1;
       ex_high_m1  = ex[KE-2];

       /* Calculate the Hy field */
       for ( k=0; k < KE-1; k++ )
       { hy[k] = hy[k] + .5*( ex[k] - ex[k+1] ) ; }

     }

       for ( k=0; k < KE; k++ )
       { printf( "%2d   %6.2f \n",k,ex[k]); }

         /* Write the E field out to a file "Ex" */
       fp = fopen( "Ex","w");
       for ( k=0; k < KE; k++ )
       { fprintf( fp,"  %6.3f \n",ex[k]); }
       fclose(fp);

       printf( "%5.1f \n",T);
}

}
